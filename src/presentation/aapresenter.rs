// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for presenting the results of the analysis of all-atom order parameters.

use indexmap::IndexMap;
use serde::Serialize;
use std::io::Write;

use crate::{
    analysis::{
        molecule::{AtomType, BondType, MoleculeType},
        topology::SystemTopology,
    },
    errors::WriteError,
    PANIC_MESSAGE,
};

use super::{AAOrder, BondResults, MoleculeResults, ResultsPresenter};

#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub(crate) struct AAOrderResults {
    molecules: Vec<AAMoleculeResults>,
}

impl From<SystemTopology> for AAOrderResults {
    #[inline(always)]
    fn from(value: SystemTopology) -> Self {
        AAOrderResults {
            molecules: value
                .molecule_types()
                .iter()
                .map(|m| AAMoleculeResults::from(m))
                .collect(),
        }
    }
}

impl ResultsPresenter for AAOrderResults {
    #[inline(always)]
    fn molecules(&self) -> &[impl MoleculeResults] {
        &self.molecules
    }

    fn write_csv_header(
        writer: &mut impl Write,
        max_bonds: usize,
        assigned_leaflets: bool,
    ) -> Result<(), WriteError> {
        write_result!(writer, "molecule,residue,atom,relative index");

        if assigned_leaflets {
            write_result!(
                writer,
                ",total full membrane,total upper leaflet,total lower leaflet"
            );
        } else {
            write_result!(writer, ",total");
        }

        for i in 1..=max_bonds {
            if assigned_leaflets {
                write_result!(writer, ",hydrogen #{} full membrane,hydrogen #{} upper leaflet,hydrogen #{} lower leaflet", i, i, i);
            } else {
                write_result!(writer, ",hydrogen #{}", i);
            }
        }

        write_result!(writer, "\n");
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize)]
struct AAMoleculeResults {
    molecule: String,
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    order: IndexMap<AtomType, AAAtomResults>,
}

impl From<&MoleculeType> for AAMoleculeResults {
    fn from(value: &MoleculeType) -> Self {
        let mut order = IndexMap::new();
        for heavy_atom in value.order_atoms().atoms() {
            let mut relevant_bonds = Vec::new();

            for bond in value.order_bonds().bond_types() {
                if bond.contains(heavy_atom) {
                    relevant_bonds.push(bond);
                }
            }

            let results = AAAtomResults::new(&relevant_bonds, heavy_atom);

            // only include atoms that have some data associated with them
            if !results.is_empty() {
                order.insert(heavy_atom.to_owned(), results);
            }
        }

        AAMoleculeResults {
            molecule: value.name().to_owned(),
            order,
        }
    }
}

impl MoleculeResults for AAMoleculeResults {
    #[inline(always)]
    fn name(&self) -> &str {
        &self.molecule
    }

    fn write_tab(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        write_result!(writer, "\nMolecule type {}\n", self.molecule);

        let max_bonds = self.max_bonds();
        if max_bonds == 0 {
            log::warn!("No atoms have any order bonds associated. Unable to write output table. Aborting...");
            return Ok(());
        }

        let leaflets = self.has_assigned_leaflets();

        write_result!(writer, "         ");

        let width = if leaflets { 28usize } else { 8usize };
        write_result!(writer, " {: ^width$} |", "TOTAL", width = width);
        for i in 1..=max_bonds {
            let hydrogen = if leaflets {
                format!("HYDROGEN #{}", i)
            } else {
                format!("H #{}", i)
            };

            write_result!(writer, " {: ^width$} |", hydrogen, width = width);
        }
        write_result!(writer, "\n");

        if leaflets {
            write_result!(writer, "         ");
            for _ in 0..=max_bonds {
                write_result!(writer, "   FULL     UPPER     LOWER   |");
            }
            write_result!(writer, "\n");
        }

        for (name, results) in self.order.iter() {
            results.write_tab(writer, name, max_bonds, leaflets)?;
        }

        Ok(())
    }

    #[inline(always)]
    fn max_bonds(&self) -> usize {
        self.order
            .values()
            .map(|x| x.bonds.len())
            .max()
            .unwrap_or(0)
    }

    /// Write results for a single molecule in an xvg format. Only order parameters for atoms are written.
    /// Order parameters are written for the full membrane and for the individual leaflets (if calculated).
    fn write_xvg(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        write_result!(
            writer,
            "@    title \"Atomistic order parameters for molecule type {}\"\n",
            self.molecule
        );
        write_result!(
            writer,
            "@    xaxis label \"Atom\"\n@    yaxis label \"-Scd\"\n"
        );

        write_result!(writer, "@    s0 legend \"Full membrane\"\n");
        if self.has_assigned_leaflets() {
            write_result!(writer, "@    s1 legend \"Upper leaflet\"\n");
            write_result!(writer, "@    s2 legend \"Lower leaflet\"\n");
        }

        write_result!(writer, "@TYPE xy\n");

        for (i, (name, results)) in self.order.iter().enumerate() {
            results.write_xvg(writer, i + 1, name)?;
        }

        Ok(())
    }

    #[inline(always)]
    fn has_assigned_leaflets(&self) -> bool {
        self.order.values().map(|x| x.upper).any(|x| x.is_some())
            || self.order.values().map(|x| x.lower).any(|x| x.is_some())
    }

    #[inline(always)]
    fn write_csv(
        &self,
        writer: &mut impl Write,
        max_bonds: usize,
        assigned_leaflets: bool,
    ) -> Result<(), WriteError> {
        for (atom_type, order) in self.order.iter() {
            order.write_csv(
                writer,
                atom_type,
                &self.molecule,
                max_bonds,
                assigned_leaflets,
            )?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Serialize)]
struct AAAtomResults {
    // 'total' is set to `None`, if there are no bonds associated with the atom
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "super::round_serialize_option_f32")]
    total: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "super::round_serialize_option_f32")]
    upper: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "super::round_serialize_option_f32")]
    lower: Option<f32>,
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    bonds: IndexMap<AtomType, BondResults>,
}

impl AAAtomResults {
    fn new(bonds: &[&BondType], heavy_atom: &AtomType) -> Self {
        let mut totals = Vec::new();
        let mut uppers = Vec::new();
        let mut lowers = Vec::new();

        let mut results = IndexMap::new();
        for bond in bonds {
            let bond_results = BondResults::convert_from::<AAOrder>(*bond);
            totals.push(bond_results.total.clone());
            uppers.push(bond_results.upper.clone());
            lowers.push(bond_results.lower.clone());

            let hydrogen = bond.get_other_atom(heavy_atom).unwrap_or_else(|| {
                panic!(
                    "FATAL GORDER ERROR | AAAtomResults::new | Heavy atom not part of bond. {}",
                    PANIC_MESSAGE
                )
            });
            results.insert(hydrogen.to_owned(), bond_results);
        }

        if results.is_empty() {
            AAAtomResults {
                total: None,
                upper: None,
                lower: None,
                bonds: results,
            }
        } else {
            let total = Some(totals.iter().sum::<f32>() / totals.len() as f32);
            let upper = average_of_some(&uppers);
            let lower = average_of_some(&lowers);

            AAAtomResults {
                total,
                upper,
                lower,
                bonds: results,
            }
        }
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }

    /// Write results for an atom in human readable table format.
    fn write_tab(
        &self,
        writer: &mut impl Write,
        atom_type: &AtomType,
        max_bonds: usize,
        leaflets: bool,
    ) -> Result<(), WriteError> {
        write_result!(writer, "{:<8} ", atom_type.atom_name());

        match self.total {
            Some(val) => {
                write_result!(writer, " {: ^8.4} ", val);

                if leaflets {
                    for value in [self.upper, self.lower] {
                        match value {
                            Some(unwrapped) => write_result!(writer, " {: ^8.4} ", unwrapped),
                            None => write_result!(writer, "               "),
                        }
                    }
                }

                write_result!(writer, "|");
            }
            None => {
                if leaflets {
                    write_result!(writer, " {: ^28} |", "");
                } else {
                    write_result!(writer, " {: ^8} |", "");
                }
            }
        }

        let mut bond_results = self.bonds.values();
        for _ in 0..max_bonds {
            match bond_results.next() {
                Some(val) => val.write_tab(writer, leaflets)?,
                None => {
                    if leaflets {
                        write_result!(writer, " {: ^28} |", "");
                    } else {
                        write_result!(writer, "          |");
                    }
                }
            }
        }

        write_result!(writer, "\n");

        Ok(())
    }

    fn write_xvg(
        &self,
        writer: &mut impl Write,
        number: usize,
        atom_type: &AtomType,
    ) -> Result<(), WriteError> {
        write_result!(writer, "# Atom {}:\n", atom_type.atom_name());
        match self.total {
            Some(order) => {
                write_result!(writer, "{:<4} {: >8.4} ", number, order);
            }
            None => {
                write_result!(writer, "{:<4} NaN ", number);
            }
        }

        for order in [self.upper, self.lower].into_iter() {
            if let Some(value) = order {
                write_result!(writer, "{: >8.4} ", value);
            }
        }

        write_result!(writer, "\n");

        Ok(())
    }

    fn write_csv(
        &self,
        writer: &mut impl Write,
        atom_type: &AtomType,
        molecule_type: &str,
        max_bonds: usize,
        assigned_leaflets: bool,
    ) -> Result<(), WriteError> {
        write_result!(
            writer,
            "{},{},{},{}",
            molecule_type,
            atom_type.residue_name(),
            atom_type.atom_name(),
            atom_type.relative_index()
        );

        if let Some(val) = self.total {
            write_result!(writer, ",{:.4}", val);
        } else {
            write_result!(writer, ",");
        }

        if assigned_leaflets {
            for value in [self.upper, self.lower] {
                match value {
                    Some(unwrapped) => write_result!(writer, ",{:.4}", unwrapped),
                    None => write_result!(writer, ","),
                }
            }
        }

        let mut bond_results = self.bonds.values();
        for _ in 0..max_bonds {
            match bond_results.next() {
                Some(val) => val.write_csv(writer, assigned_leaflets)?,
                None => {
                    if assigned_leaflets {
                        write_result!(writer, ",,,");
                    } else {
                        write_result!(writer, ",");
                    }
                }
            }
        }

        write_result!(writer, "\n");

        Ok(())
    }
}

fn average_of_some(values: &[Option<f32>]) -> Option<f32> {
    let mut all_none = true;
    let mut all_nan = true;

    let (sum, count) = values
        .iter()
        .filter_map(|&x| {
            x.map(|n| {
                // at least one value is Some, not None
                all_none = false;
                // at least one value is not NaN
                if !n.is_nan() {
                    all_nan = false;
                    n
                } else {
                    f32::NAN
                }
            })
        })
        .fold((0.0, 0), |(acc, cnt), x| {
            if !x.is_nan() {
                (acc + x, cnt + 1)
            } else {
                (acc, cnt)
            }
        });

    if all_none {
        // all values were None
        None
    } else if all_nan {
        // all non-None values were NaN
        Some(f32::NAN)
    } else if count > 0 {
        // average of not-none values
        Some(sum / count as f32)
    } else {
        panic!("FATAL GORDER ERROR | presentation::aapresenter::average_of_some | Impossible branch reached. {}", PANIC_MESSAGE);
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use groan_rs::prelude::{Atom, SimBox};

    use super::*;

    fn prepare_example_atom_results_1() -> AAAtomResults {
        let mut results = AAAtomResults {
            total: Some(0.777),
            upper: None,
            lower: None,
            bonds: IndexMap::new(),
        };

        let bond1_results = BondResults {
            total: 0.721,
            upper: None,
            lower: None,
        };

        let bond2_results = BondResults {
            total: 0.833,
            upper: None,
            lower: None,
        };

        results
            .bonds
            .insert(AtomType::new_raw(1, "POPE", "HA"), bond1_results);

        results
            .bonds
            .insert(AtomType::new_raw(3, "POPE", "HB"), bond2_results);

        results
    }

    #[test]
    fn test_serialize_aa_atom_results() {
        let results = prepare_example_atom_results_1();
        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(
            serialized,
            "total: 0.777
bonds:
  POPE HA (1):
    total: 0.721
  POPE HB (3):
    total: 0.833
"
        );
    }

    #[test]
    fn test_serialize_aa_atom_results_empty() {
        let results = AAAtomResults {
            total: None,
            upper: None,
            lower: None,
            bonds: IndexMap::new(),
        };

        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(serialized, "{}\n");
    }

    fn prepare_example_atom_results_2() -> AAAtomResults {
        let mut results = AAAtomResults {
            total: Some(0.777),
            upper: Some(0.743),
            lower: Some(0.811),
            bonds: IndexMap::new(),
        };

        let bond1_results = BondResults {
            total: 0.721,
            upper: Some(0.687),
            lower: Some(0.755),
        };

        let bond2_results = BondResults {
            total: 0.833,
            upper: Some(0.854),
            lower: Some(0.812),
        };

        results
            .bonds
            .insert(AtomType::new_raw(7, "POPE", "HA"), bond1_results);

        results
            .bonds
            .insert(AtomType::new_raw(11, "POPE", "HB"), bond2_results);

        results
    }

    #[test]
    fn test_serialize_aa_atom_results_leaflets() {
        let results = prepare_example_atom_results_2();
        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(
            serialized,
            "total: 0.777
upper: 0.743
lower: 0.811
bonds:
  POPE HA (7):
    total: 0.721
    upper: 0.687
    lower: 0.755
  POPE HB (11):
    total: 0.833
    upper: 0.854
    lower: 0.812
"
        );
    }

    #[test]
    fn test_serialize_aa_molecule_results() {
        let atom_results_1 = prepare_example_atom_results_1();
        let atom_type_1 = AtomType::new_raw(0, "POPE", "CA");
        let atom_results_2 = prepare_example_atom_results_2();
        let atom_type_2 = AtomType::new_raw(5, "POPE", "CB");
        let molecule_results = AAMoleculeResults {
            molecule: "POPE".to_owned(),
            order: IndexMap::from_iter(
                [(atom_type_1, atom_results_1), (atom_type_2, atom_results_2)].into_iter(),
            ),
        };

        let serialized = serde_yaml::to_string(&molecule_results).unwrap();
        assert_eq!(
            serialized,
            "molecule: POPE
order:
  POPE CA (0):
    total: 0.777
    bonds:
      POPE HA (1):
        total: 0.721
      POPE HB (3):
        total: 0.833
  POPE CB (5):
    total: 0.777
    upper: 0.743
    lower: 0.811
    bonds:
      POPE HA (7):
        total: 0.721
        upper: 0.687
        lower: 0.755
      POPE HB (11):
        total: 0.833
        upper: 0.854
        lower: 0.812
"
        )
    }

    #[test]
    fn test_serialize_aa_results() {
        let atom_results_1 = prepare_example_atom_results_1();
        let atom_type_1 = AtomType::new_raw(0, "POPE", "CA");
        let atom_results_2 = prepare_example_atom_results_2();
        let atom_type_2 = AtomType::new_raw(5, "POPE", "CB");
        let molecule_results = AAMoleculeResults {
            molecule: "POPE".to_owned(),
            order: IndexMap::from_iter(
                [(atom_type_1, atom_results_1), (atom_type_2, atom_results_2)].into_iter(),
            ),
        };
        let results = AAOrderResults {
            molecules: vec![molecule_results],
        };

        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(
            serialized,
            "- molecule: POPE
  order:
    POPE CA (0):
      total: 0.777
      bonds:
        POPE HA (1):
          total: 0.721
        POPE HB (3):
          total: 0.833
    POPE CB (5):
      total: 0.777
      upper: 0.743
      lower: 0.811
      bonds:
        POPE HA (7):
          total: 0.721
          upper: 0.687
          lower: 0.755
        POPE HB (11):
          total: 0.833
          upper: 0.854
          lower: 0.812
"
        );
    }

    #[test]
    fn test_aaatom_results_new() {
        let heavy_atom = Atom::new(1, "POPC", 1, "C2");
        let hydrogen1 = Atom::new(1, "POPC", 3, "HA");
        let hydrogen2 = Atom::new(1, "POPC", 4, "HB");

        let mut bond1 = BondType::new(
            0,
            &heavy_atom,
            2,
            &hydrogen1,
            0,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();
        let mut bond2 = BondType::new(
            0,
            &heavy_atom,
            3,
            &hydrogen2,
            0,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        *bond1.total_mut() += 0.234;
        *bond1.total_mut() += 0.176;
        *bond1.total_mut() += 0.112;
        // average: 0.174

        *bond2.total_mut() += 0.112;
        *bond2.total_mut() += 0.245;
        *bond2.total_mut() += -0.013;
        // average: 0.114666

        let heavy_atom_type = AtomType::new_raw(0, "POPC", "C2");
        let hydrogen1_atom_type = AtomType::new_raw(2, "POPC", "HA");
        let hydrogen2_atom_type = AtomType::new_raw(3, "POPC", "HB");

        let results = AAAtomResults::new(&[&bond1, &bond2], &heavy_atom_type);

        assert_eq!(results.bonds.len(), 2);

        let bond1_results = results.bonds.get(&hydrogen1_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, -0.174);
        assert!(bond1_results.upper.is_none());
        assert!(bond1_results.lower.is_none());

        let bond1_results = results.bonds.get(&hydrogen2_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, -0.114666);
        assert!(bond1_results.upper.is_none());
        assert!(bond1_results.lower.is_none());

        assert_relative_eq!(results.total.unwrap(), -0.144333);
        assert!(results.upper.is_none());
        assert!(results.lower.is_none());
    }

    #[test]
    fn test_aaatom_results_new_leaflets() {
        let heavy_atom = Atom::new(1, "POPC", 1, "C2");
        let hydrogen1 = Atom::new(1, "POPC", 3, "HA");
        let hydrogen2 = Atom::new(1, "POPC", 4, "HB");

        let mut bond1 = BondType::new(
            0,
            &heavy_atom,
            2,
            &hydrogen1,
            0,
            true,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();
        let mut bond2 = BondType::new(
            0,
            &heavy_atom,
            3,
            &hydrogen2,
            0,
            true,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        *bond1.total_mut() += 0.234;
        *bond1.total_mut() += 0.176;
        *bond1.total_mut() += 0.112;
        // average: 0.174

        *bond1.upper_mut().as_mut().unwrap() += 0.115;
        *bond1.upper_mut().as_mut().unwrap() += 0.325;
        *bond1.upper_mut().as_mut().unwrap() += 0.008;
        // average: 0.149333

        *bond2.total_mut() += 0.112;
        *bond2.total_mut() += 0.245;
        *bond2.total_mut() += -0.013;
        // average: 0.114666

        *bond2.lower_mut().as_mut().unwrap() += 0.332;
        *bond2.lower_mut().as_mut().unwrap() += 0.087;
        *bond2.lower_mut().as_mut().unwrap() += 0.125;
        // average: 0.181333

        let heavy_atom_type = AtomType::new_raw(0, "POPC", "C2");
        let hydrogen1_atom_type = AtomType::new_raw(2, "POPC", "HA");
        let hydrogen2_atom_type = AtomType::new_raw(3, "POPC", "HB");

        let results = AAAtomResults::new(&[&bond1, &bond2], &heavy_atom_type);

        assert_eq!(results.bonds.len(), 2);

        let bond1_results = results.bonds.get(&hydrogen1_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, -0.174);
        assert_relative_eq!(bond1_results.upper.unwrap(), -0.149333);
        assert!(bond1_results.lower.unwrap().is_nan()); // no values for the lower leaflet

        let bond1_results = results.bonds.get(&hydrogen2_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, -0.114666);
        assert!(bond1_results.upper.unwrap().is_nan()); // no values for the upper leaflet
        assert_relative_eq!(bond1_results.lower.unwrap(), -0.181333);

        assert_relative_eq!(results.total.unwrap(), -0.144333);
        assert_relative_eq!(results.upper.unwrap(), -0.149333);
        assert_relative_eq!(results.lower.unwrap(), -0.181333);
    }

    #[test]
    fn test_average_of_some() {
        let values = [Some(1.0), Some(1.0), Some(3.0), Some(2.0), Some(3.0)];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [Some(1.0), Some(2.0), Some(3.0), Some(2.0), Some(f32::NAN)];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [Some(1.0), None, Some(3.0), Some(2.0), None];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [Some(1.0), Some(f32::NAN), Some(3.0), Some(2.0), None];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [None, Some(f32::NAN), None, Some(2.0), None];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [None, None, None, Some(2.0), None];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [
            Some(f32::NAN),
            Some(f32::NAN),
            Some(f32::NAN),
            Some(2.0),
            Some(f32::NAN),
        ];
        assert_relative_eq!(average_of_some(&values).unwrap(), 2.0);

        let values = [None, None, None, None, None];
        assert!(average_of_some(&values).is_none());

        let values = [
            Some(f32::NAN),
            Some(f32::NAN),
            Some(f32::NAN),
            Some(f32::NAN),
        ];
        assert!(average_of_some(&values).unwrap().is_nan());

        let values = [Some(f32::NAN), None, Some(f32::NAN), Some(f32::NAN)];
        assert!(average_of_some(&values).unwrap().is_nan());

        let values = [None, None, None, Some(f32::NAN), None];
        assert!(average_of_some(&values).unwrap().is_nan());
    }
}
