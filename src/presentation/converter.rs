// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for converting between system topology and formatted results.

use crate::analysis::molecule::{AtomType, BondType, MoleculeType};
use crate::analysis::order::{add_option_order, AnalysisOrder};
use crate::analysis::ordermap::{add_option_map, Map};
use crate::analysis::timewise::{AddSum, TimeWiseAddTreatment};
use crate::analysis::topology::SystemTopology;
use crate::input::Analysis;
use crate::presentation::aaresults::{AAAtomResults, AAMoleculeResults, AAOrderResults};
use crate::presentation::cgresults::{CGMoleculeResults, CGOrderResults};
use crate::presentation::OrderMapsCollection;
use crate::presentation::{
    BondResults, GridMapF32, Order, OrderCollection, OrderResults, OrderType,
};
use crate::PANIC_MESSAGE;
use groan_rs::structures::gridmap::DataOrder;
use indexmap::IndexMap;
use std::marker::PhantomData;
use std::num::NonZeroUsize;
use std::ops::AddAssign;

/// Structure converting from system topology to formatted results.
#[derive(Debug, Clone)]
pub(crate) struct ResultsConverter<O: MolConvert> {
    analysis: Analysis,
    _type: PhantomData<O>,
}

impl<O: MolConvert> ResultsConverter<O> {
    #[inline(always)]
    pub(crate) fn new(analysis: Analysis) -> Self {
        Self {
            analysis,
            _type: PhantomData,
        }
    }

    /// Convert topology with raw calculated order parameters into a presentable structure.

    pub(crate) fn convert_topology(self, topology: SystemTopology) -> O {
        O::new(
            topology
                .molecule_types()
                .iter()
                .map(|x| x.name().clone())
                .zip(
                    topology
                        .molecule_types()
                        .iter()
                        .map(|mol| O::convert_molecule(mol, &self)),
                )
                .collect::<IndexMap<String, O::MoleculeResults>>(),
            self.analysis,
        )
    }

    /// Convert raw order from analysis into a presentable order parameter.
    #[inline]
    fn convert_order<T: TimeWiseAddTreatment>(&self, analysis_order: &AnalysisOrder<T>) -> Order {
        let min_samples = NonZeroUsize::new(self.analysis.min_samples()).expect(PANIC_MESSAGE);
        let order = analysis_order.calc_order(min_samples);
        let error = self
            .analysis
            .estimate_error()
            .as_ref()
            .map(|e| e.n_blocks())
            .map(|blocks| analysis_order.estimate_error(blocks, min_samples))
            .flatten();
        O::OrderType::convert(order, error)
    }

    /// Convert raw order calculated for a single bond type into presentable order parameters.
    fn convert_bond(&self, bond_type: &BondType, molecule: &str) -> BondResults {
        let total = self.convert_order(bond_type.total());
        let upper = bond_type.upper().as_ref().map(|x| self.convert_order(x));
        let lower = bond_type.lower().as_ref().map(|x| self.convert_order(x));

        let total_ordermap = bond_type
            .total_map()
            .as_ref()
            .map(|x| Self::convert_ordermap(x));
        let upper_ordermap = bond_type
            .upper_map()
            .as_ref()
            .map(|x| Self::convert_ordermap(x));
        let lower_ordermap = bond_type
            .lower_map()
            .as_ref()
            .map(|x| Self::convert_ordermap(x));

        BondResults::new(
            bond_type.bond_topology(),
            molecule,
            OrderCollection::new(Some(total), upper, lower),
            OrderMapsCollection::new(total_ordermap, upper_ordermap, lower_ordermap),
        )
    }

    /// Convert raw ordermap into a presentable ordermap.
    fn convert_ordermap(ordermap: &Map) -> GridMapF32 {
        let min_samples = ordermap.params().min_samples();

        let extracted: Vec<f32> = ordermap
            .values()
            .extract_convert()
            .zip(ordermap.samples().extract_convert())
            .map(|(value, samples)| {
                if samples.2 < min_samples {
                    f32::NAN
                } else {
                    O::OrderType::convert(value.2 / samples.2 as f32, None).value
                }
            })
            .collect();

        GridMapF32::from_vec(
            ordermap.values().span_x(),
            ordermap.values().span_y(),
            ordermap.values().tile_dim(),
            extracted,
            DataOrder::default(),
            f32::clone as fn(&f32) -> f32,
        )
        .unwrap_or_else(|e| {
            panic!(
                "FATAL GORDER ERROR | ResultsConverter::convert_ordermap | \
            Could not convert Map to GridMapF32 ({}). {}",
                e, PANIC_MESSAGE
            )
        })
    }

    /// Convert temporary OrderSummer structure.
    #[inline]
    fn convert_order_summer(&self, summer: &OrderSummer) -> (OrderCollection, OrderMapsCollection) {
        let total = summer.total.as_ref().map(|x| self.convert_order(x));
        let upper = summer.upper.as_ref().map(|x| self.convert_order(x));
        let lower = summer.lower.as_ref().map(|x| self.convert_order(x));

        let total_ordermap = summer
            .ordermap_total
            .as_ref()
            .map(|x| Self::convert_ordermap(x));
        let upper_ordermap = summer
            .ordermap_upper
            .as_ref()
            .map(|x| Self::convert_ordermap(x));
        let lower_ordermap = summer
            .ordermap_lower
            .as_ref()
            .map(|x| Self::convert_ordermap(x));

        (
            OrderCollection::new(total, upper, lower),
            OrderMapsCollection::new(total_ordermap, upper_ordermap, lower_ordermap),
        )
    }
}

/// Trait implemented for all `OrderResults` structures allowing to convert between raw molecules
/// and molecule results.
pub(crate) trait MolConvert: OrderResults {
    /// Convert a molecule type with raw calculated order parameters into a presentable structure.
    fn convert_molecule(
        molecule_type: &MoleculeType,
        converter: &ResultsConverter<Self>,
    ) -> Self::MoleculeResults;
}

impl MolConvert for AAOrderResults {
    fn convert_molecule(
        molecule_type: &MoleculeType,
        converter: &ResultsConverter<Self>,
    ) -> Self::MoleculeResults {
        let mut order = IndexMap::new();
        let mut summer = OrderSummer::default();

        for heavy_atom in molecule_type.order_atoms().atoms() {
            let mut relevant_bonds = Vec::new();

            for bond in molecule_type.order_bonds().bond_types() {
                if bond.contains(heavy_atom) {
                    relevant_bonds.push(bond);
                    summer += bond;
                }
            }

            let results =
                Self::convert_atom(heavy_atom, molecule_type.name(), &relevant_bonds, converter);
            // only include atoms that have some data associated with them
            if !results.is_empty() {
                order.insert(heavy_atom.clone(), results);
            }
        }

        let (average, ordermaps) = converter.convert_order_summer(&summer);
        AAMoleculeResults::new(molecule_type.name(), average, ordermaps, order)
    }
}

impl AAOrderResults {
    /// Convert results collected for individual bonds into a result for a heavy atom.
    fn convert_atom(
        heavy_atom: &AtomType,
        molecule: &str,
        bonds: &[&BondType],
        converter: &ResultsConverter<Self>,
    ) -> AAAtomResults {
        let mut results = IndexMap::new();
        let mut summer = OrderSummer::default();

        for bond in bonds {
            summer += bond;

            let bond_results = converter.convert_bond(bond, molecule);
            let hydrogen = bond.get_other_atom(heavy_atom).unwrap_or_else(|| {
                panic!(
                    "FATAL GORDER ERROR | AAOrderResults::convert_atom | Heavy atom not part of the bond. {}",
                    PANIC_MESSAGE
                )
            });

            results.insert(hydrogen.clone(), bond_results);
        }

        if results.is_empty() {
            AAAtomResults::new(
                heavy_atom.clone(),
                molecule,
                OrderCollection::default(),
                OrderMapsCollection::default(),
                results,
            )
        } else {
            let (average, ordermaps) = converter.convert_order_summer(&summer);
            AAAtomResults::new(heavy_atom.clone(), molecule, average, ordermaps, results)
        }
    }
}

impl MolConvert for CGOrderResults {
    fn convert_molecule(
        molecule_type: &MoleculeType,
        converter: &ResultsConverter<Self>,
    ) -> Self::MoleculeResults {
        let mut order = IndexMap::new();
        let mut summer = OrderSummer::default();

        for bond in molecule_type.order_bonds().bond_types() {
            summer += bond;
            let results = converter.convert_bond(bond, molecule_type.name());
            order.insert(bond.bond_topology().clone(), results);
        }

        let (average, ordermaps) = converter.convert_order_summer(&summer);
        CGMoleculeResults::new(molecule_type.name(), average, ordermaps, order)
    }
}

/// Helper struct for summing order parameters and ordermaps.
#[derive(Debug, Clone, Default)]
struct OrderSummer {
    total: Option<AnalysisOrder<AddSum>>,
    upper: Option<AnalysisOrder<AddSum>>,
    lower: Option<AnalysisOrder<AddSum>>,
    ordermap_total: Option<Map>,
    ordermap_upper: Option<Map>,
    ordermap_lower: Option<Map>,
}

impl AddAssign<&BondType> for OrderSummer {
    #[inline(always)]
    fn add_assign(&mut self, rhs: &BondType) {
        add_option_order(&mut self.total, Some(rhs.total().clone()));
        add_option_order(&mut self.upper, rhs.upper().clone());
        add_option_order(&mut self.lower, rhs.lower().clone());

        add_option_map(&mut self.ordermap_total, rhs.total_map().clone());
        add_option_map(&mut self.ordermap_upper, rhs.upper_map().clone());
        add_option_map(&mut self.ordermap_lower, rhs.lower_map().clone());
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use groan_rs::prelude::{SimBox, Vector3D};

    use crate::{
        analysis::timewise::AddExtend,
        input::{AnalysisType, EstimateError, OrderMap, Plane},
    };

    use super::*;

    fn prepare_converter_aa() -> ResultsConverter<AAOrderResults> {
        let analysis = Analysis::new()
            .structure("system.gro")
            .trajectory("md.xtc")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@hydrogens and element name carbon",
            ))
            .output("order.yaml")
            .min_samples(50)
            .estimate_error(EstimateError::default())
            .build()
            .unwrap();

        ResultsConverter::new(analysis)
    }

    fn prepare_converter_cg() -> ResultsConverter<CGOrderResults> {
        let analysis = Analysis::new()
            .structure("system.gro")
            .trajectory("md.xtc")
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .output("order.yaml")
            .min_samples(50)
            .estimate_error(EstimateError::default())
            .build()
            .unwrap();

        ResultsConverter::new(analysis)
    }

    #[test]
    fn convert_order_simple_aa() {
        let analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 56, false);
        let converter = prepare_converter_aa();
        let order = converter.convert_order(&analysis_order);
        assert_relative_eq!(order.value(), -0.8092857, epsilon = 1e-4);
        assert!(order.error().is_none());

        let analysis_order = AnalysisOrder::<AddSum>::new(45.32, 56, false);
        let converter = prepare_converter_aa();
        let order = converter.convert_order(&analysis_order);
        assert_relative_eq!(order.value(), -0.8092857, epsilon = 1e-4);
        assert!(order.error().is_none());
    }

    #[test]
    fn convert_order_simple_limit_aa() {
        let analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 48, false);
        let converter = prepare_converter_aa();
        let order = converter.convert_order(&analysis_order);
        assert!(order.value().is_nan());
        assert!(order.error().is_none());

        let analysis_order = AnalysisOrder::<AddSum>::new(45.32, 48, false);
        let converter = prepare_converter_aa();
        let order = converter.convert_order(&analysis_order);
        assert!(order.value().is_nan());
        assert!(order.error().is_none());
    }

    #[test]
    fn convert_order_with_error_aa() {
        let mut analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 60, true);
        // we add more values, but since the average is always the same, it does not matter
        analysis_order.init_new_frame();
        analysis_order += 0.812;
        analysis_order.init_new_frame();
        analysis_order += 0.684;
        analysis_order.init_new_frame();
        analysis_order += 0.692;
        analysis_order.init_new_frame();
        analysis_order += 0.766;
        analysis_order.init_new_frame();
        analysis_order += 0.802;
        // this frame will be unused in the error analysis
        analysis_order.init_new_frame();
        analysis_order += 0.776;

        let converter = prepare_converter_aa();
        let order = converter.convert_order(&analysis_order);
        assert_relative_eq!(order.value(), -0.755333, epsilon = 1e-4);
        assert_relative_eq!(order.error().unwrap(), 0.0602428, epsilon = 1e-4);
    }

    #[test]
    fn convert_order_with_error_limit_aa() {
        let mut analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 43, true);
        // we add more values, but since the average is always the same, it does not matter
        analysis_order.init_new_frame();
        analysis_order += 0.812;
        analysis_order.init_new_frame();
        analysis_order += 0.684;
        analysis_order.init_new_frame();
        analysis_order += 0.692;
        analysis_order.init_new_frame();
        analysis_order += 0.766;
        analysis_order.init_new_frame();
        analysis_order += 0.802;
        // this frame will be unused in the error analysis
        analysis_order.init_new_frame();
        analysis_order += 0.776;

        let converter = prepare_converter_aa();
        let order = converter.convert_order(&analysis_order);
        assert!(order.value().is_nan());
        assert!(order.error().unwrap().is_nan());
    }

    #[test]
    fn convert_order_simple_cg() {
        let analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 56, false);
        let converter = prepare_converter_cg();
        let order = converter.convert_order(&analysis_order);
        assert_relative_eq!(order.value(), 0.8092857, epsilon = 1e-4);
        assert!(order.error().is_none());

        let analysis_order = AnalysisOrder::<AddSum>::new(45.32, 56, false);
        let converter = prepare_converter_cg();
        let order = converter.convert_order(&analysis_order);
        assert_relative_eq!(order.value(), 0.8092857, epsilon = 1e-4);
        assert!(order.error().is_none());
    }

    #[test]
    fn convert_order_simple_limit_cg() {
        let analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 48, false);
        let converter = prepare_converter_cg();
        let order = converter.convert_order(&analysis_order);
        assert!(order.value().is_nan());
        assert!(order.error().is_none());

        let analysis_order = AnalysisOrder::<AddSum>::new(45.32, 48, false);
        let converter = prepare_converter_cg();
        let order = converter.convert_order(&analysis_order);
        assert!(order.value().is_nan());
        assert!(order.error().is_none());
    }

    #[test]
    fn convert_order_with_error_cg() {
        let mut analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 60, true);
        // we add more values, but since the average is always the same, it does not matter
        analysis_order.init_new_frame();
        analysis_order += 0.812;
        analysis_order.init_new_frame();
        analysis_order += 0.684;
        analysis_order.init_new_frame();
        analysis_order += 0.692;
        analysis_order.init_new_frame();
        analysis_order += 0.766;
        analysis_order.init_new_frame();
        analysis_order += 0.802;
        // this frame will be unused in the error analysis
        analysis_order.init_new_frame();
        analysis_order += 0.776;

        let converter = prepare_converter_cg();
        let order = converter.convert_order(&analysis_order);
        assert_relative_eq!(order.value(), 0.755333, epsilon = 1e-4);
        assert_relative_eq!(order.error().unwrap(), 0.0602428, epsilon = 1e-4);
    }

    #[test]
    fn convert_order_with_error_limit_cg() {
        let mut analysis_order = AnalysisOrder::<AddExtend>::new(45.32, 43, true);
        // we add more values, but since the average is always the same, it does not matter
        analysis_order.init_new_frame();
        analysis_order += 0.812;
        analysis_order.init_new_frame();
        analysis_order += 0.684;
        analysis_order.init_new_frame();
        analysis_order += 0.692;
        analysis_order.init_new_frame();
        analysis_order += 0.766;
        analysis_order.init_new_frame();
        analysis_order += 0.802;
        // this frame will be unused in the error analysis
        analysis_order.init_new_frame();
        analysis_order += 0.776;

        let converter = prepare_converter_cg();
        let order = converter.convert_order(&analysis_order);
        assert!(order.value().is_nan());
        assert!(order.error().unwrap().is_nan());
    }

    fn prepare_ordermap_params() -> OrderMap {
        OrderMap::new()
            .output_directory("ordermaps")
            .plane(Plane::XY)
            .bin_size([1.0, 1.0])
            .build()
            .unwrap()
    }

    #[test]
    fn convert_ordermap() {
        let mut ordermap =
            Map::new(prepare_ordermap_params(), &SimBox::from([5.0, 5.0, 5.0])).unwrap();

        ordermap.add_order(0.10, &Vector3D::new(1.0, 1.0, 0.0));
        ordermap.add_order(0.12, &Vector3D::new(1.0, 1.0, 0.0));
        ordermap.add_order(0.12, &Vector3D::new(1.0, 2.0, 0.0));
        ordermap.add_order(0.15, &Vector3D::new(1.0, 5.0, 0.0));

        ordermap.add_order(0.32, &Vector3D::new(3.0, 2.0, 0.0));
        ordermap.add_order(0.33, &Vector3D::new(3.0, 3.0, 0.0));
        ordermap.add_order(0.34, &Vector3D::new(3.0, 4.0, 0.0));

        ordermap.add_order(0.40, &Vector3D::new(4.0, 4.0, 0.0));
        ordermap.add_order(0.48, &Vector3D::new(4.0, 4.0, 0.0));

        for (map, multiplier) in [
            ResultsConverter::<AAOrderResults>::convert_ordermap(&ordermap),
            ResultsConverter::<CGOrderResults>::convert_ordermap(&ordermap),
        ]
        .into_iter()
        .zip([-1.0, 1.0].into_iter())
        {
            assert_eq!(map.extract_convert().count(), 36);

            for x in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0] {
                for y in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0] {
                    let expected = match (x, y) {
                        (1.0, 1.0) => multiplier * 0.11,
                        (1.0, 2.0) => multiplier * 0.12,
                        (1.0, 5.0) => multiplier * 0.15,
                        (3.0, 2.0) => multiplier * 0.32,
                        (3.0, 3.0) => multiplier * 0.33,
                        (3.0, 4.0) => multiplier * 0.34,
                        (4.0, 4.0) => multiplier * 0.44,
                        _ => f32::NAN,
                    };

                    if expected.is_nan() {
                        assert!(map.get_at(x, y).unwrap().is_nan())
                    } else {
                        assert_relative_eq!(*map.get_at(x, y).unwrap(), expected, epsilon = 1e-4);
                    }
                }
            }
        }
    }

    fn prepare_ordermap_params_limit() -> OrderMap {
        OrderMap::new()
            .output_directory("ordermaps")
            .plane(Plane::XY)
            .bin_size([1.0, 1.0])
            .min_samples(2)
            .build()
            .unwrap()
    }

    #[test]
    fn convert_ordermap_limit() {
        let mut ordermap = Map::new(
            prepare_ordermap_params_limit(),
            &SimBox::from([5.0, 5.0, 5.0]),
        )
        .unwrap();

        ordermap.add_order(0.10, &Vector3D::new(1.0, 1.0, 0.0));
        ordermap.add_order(0.12, &Vector3D::new(1.0, 1.0, 0.0));
        ordermap.add_order(0.12, &Vector3D::new(1.0, 2.0, 0.0));
        ordermap.add_order(0.15, &Vector3D::new(1.0, 5.0, 0.0));

        ordermap.add_order(0.32, &Vector3D::new(3.0, 2.0, 0.0));
        ordermap.add_order(0.33, &Vector3D::new(3.0, 3.0, 0.0));
        ordermap.add_order(0.34, &Vector3D::new(3.0, 4.0, 0.0));

        ordermap.add_order(0.40, &Vector3D::new(4.0, 4.0, 0.0));
        ordermap.add_order(0.48, &Vector3D::new(4.0, 4.0, 0.0));

        for (map, multiplier) in [
            ResultsConverter::<AAOrderResults>::convert_ordermap(&ordermap),
            ResultsConverter::<CGOrderResults>::convert_ordermap(&ordermap),
        ]
        .into_iter()
        .zip([-1.0, 1.0].into_iter())
        {
            assert_eq!(map.extract_convert().count(), 36);

            for x in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0] {
                for y in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0] {
                    let expected = match (x, y) {
                        (1.0, 1.0) => multiplier * 0.11,
                        (4.0, 4.0) => multiplier * 0.44,
                        _ => f32::NAN,
                    };

                    if expected.is_nan() {
                        assert!(map.get_at(x, y).unwrap().is_nan())
                    } else {
                        assert_relative_eq!(*map.get_at(x, y).unwrap(), expected, epsilon = 1e-4);
                    }
                }
            }
        }
    }
}
