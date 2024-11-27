// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use core::fmt;
use std::{
    num::NonZeroUsize,
    ops::{Add, AddAssign},
};

use getset::{CopyGetters, Getters, MutGetters};
use groan_rs::{
    prelude::{Atom, SimBox, Vector3D},
    structures::group::Group,
    system::System,
};
use hashbrown::HashSet;

use crate::{
    errors::{AnalysisError, TopologyError},
    Leaflet, OrderMap, PANIC_MESSAGE,
};

use super::{
    leaflets::MoleculeLeafletClassification,
    ordermap::{merge_option_maps, Map},
};

#[derive(Debug, Clone, Getters, MutGetters)]
pub(crate) struct MoleculeType {
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    name: String,
    #[getset(get = "pub(super)")]
    topology: MoleculeTopology,
    #[getset(get = "pub(crate)")]
    order_bonds: OrderBonds,
    #[getset(get = "pub(crate)")]
    order_atoms: OrderAtoms,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    leaflet_classification: Option<MoleculeLeafletClassification>,
}

impl MoleculeType {
    pub(super) fn new(
        system: &System,
        name: String,
        topology: MoleculeTopology,
        order_bonds: &HashSet<(usize, usize)>,
        order_atoms: &[usize],
        min_index: usize,
        leaflet_classification: Option<MoleculeLeafletClassification>,
        ordermap_params: Option<&OrderMap>,
        min_samples: usize,
    ) -> Result<Self, TopologyError> {
        Ok(Self {
            name: name,
            topology: topology,
            order_bonds: OrderBonds::new(
                system,
                order_bonds,
                min_index,
                leaflet_classification.is_some(),
                ordermap_params,
                min_samples,
            )?,
            order_atoms: OrderAtoms::new(system, order_atoms, min_index),
            leaflet_classification,
        })
    }

    /// Insert new bond instances to the molecule.
    #[inline]
    pub(super) fn insert(
        &mut self,
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        atoms: Group,
    ) -> Result<(), TopologyError> {
        self.order_bonds.insert(
            system,
            bonds,
            atoms.get_atoms().first().expect(PANIC_MESSAGE),
        );

        if let Some(classifier) = self.leaflet_classification.as_mut() {
            classifier.insert(&atoms, system)?;
        }

        Ok(())
    }

    /// Calculate order parameters for bonds of a single molecule type from a single simulation frame.
    #[inline]
    pub(super) fn analyze_frame(
        &mut self,
        frame: &System,
        simbox: &SimBox,
        membrane_normal: &Vector3D,
    ) -> Result<(), AnalysisError> {
        self.order_bonds
            .bond_types
            .iter_mut()
            .try_for_each(|bond_type| {
                bond_type.analyze_frame(
                    frame,
                    self.leaflet_classification.as_ref(),
                    simbox,
                    membrane_normal,
                )
            })?;

        Ok(())
    }
}

impl Add<MoleculeType> for MoleculeType {
    type Output = Self;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            name: self.name,
            topology: self.topology,
            order_bonds: self.order_bonds + rhs.order_bonds,
            order_atoms: self.order_atoms,
            leaflet_classification: self.leaflet_classification,
        }
    }
}

/// Collection of all bond types in a molecule describing its topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct MoleculeTopology {
    pub(super) bonds: HashSet<BondTopology>,
}

impl MoleculeTopology {
    /// Create new molecule topology from a set of bonds (absolute indices) and the minimum index in the molecule.
    ///
    /// ## Panics
    /// - Panics if the `min_index` is higher than any index inside the `bonds` set.
    /// - Panics if there is a bond connecting the same atom (e.g. 14-14) in the `bonds` set.
    /// - Panics if an index in the `bonds` set does not correspond to an existing atom.
    pub(super) fn new(
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        min_index: usize,
    ) -> MoleculeTopology {
        bonds_sanity_check(bonds, min_index);

        let mut converted_bonds = HashSet::new();
        for &(index1, index2) in bonds {
            let (atom1, atom2) = get_atoms_from_bond(system, index1, index2);
            let bond = BondTopology::new(index1 - min_index, atom1, index2 - min_index, atom2);

            if !converted_bonds.insert(bond) {
                panic!(
                    "FATAL GORDER ERROR | MoleculeTopology::new | Bond between atoms '{}' and '{}' defined multiple times. {}", 
                    index1, index2, PANIC_MESSAGE
                );
            }
        }

        MoleculeTopology {
            bonds: converted_bonds,
        }
    }
}

/// Collection of all bonds for which the order parameters should be calculated.
#[derive(Debug, Clone, Getters, MutGetters)]
pub(crate) struct OrderBonds {
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    bond_types: Vec<BondType>,
}

impl OrderBonds {
    /// Create a new `OrderBonds` structure from a set of bonds (absolute indices) and the minimum index of the molecule.
    ///
    /// ## Panics
    /// - Panics if the `min_index` is higher than any index inside the `bonds` set.
    /// - Panics if there is a bond connecting the same atom (e.g. 14-14) in the `bonds` set.
    /// - Panics if an index in the `bonds` set does not correspond to an existing atom.
    pub(super) fn new(
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        min_index: usize,
        classify_leaflets: bool,
        ordermap: Option<&OrderMap>,
        min_samples: usize,
    ) -> Result<Self, TopologyError> {
        bonds_sanity_check(bonds, min_index);

        let simbox = system
            .get_box()
            .ok_or_else(|| TopologyError::UndefinedBox)?;

        if !simbox.is_orthogonal() {
            return Err(TopologyError::NotOrthogonalBox);
        }

        if simbox.is_zero() {
            return Err(TopologyError::ZeroBox);
        }

        let mut order_bonds = Vec::new();
        for &(index1, index2) in bonds.iter() {
            let (atom1, atom2) = get_atoms_from_bond(system, index1, index2);
            let bond = BondType::new(
                index1,
                atom1,
                index2,
                atom2,
                min_index,
                classify_leaflets,
                ordermap,
                min_samples,
                simbox,
            )?;
            order_bonds.push(bond)
        }

        // sort order bonds so that atoms with lower indices come first
        order_bonds.sort_by(|b1, b2| {
            b1.atom1_index()
                .cmp(&b2.atom1_index())
                .then_with(|| b1.atom2_index().cmp(&b2.atom2_index()))
        });

        Ok(OrderBonds {
            bond_types: order_bonds,
        })
    }

    /// Insert new real bonds to already constructed order bonds.
    pub(super) fn insert(
        &mut self,
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        min_index: usize,
    ) {
        bonds_sanity_check(bonds, min_index);

        for &(index1, index2) in bonds.iter() {
            let (atom1, atom2) = get_atoms_from_bond(system, index1, index2);
            let bond_topology =
                BondTopology::new(index1 - min_index, atom1, index2 - min_index, atom2);

            if let Some(order_bond) = self
                .bond_types
                .iter_mut()
                .find(|order_bond| order_bond.bond_topology == bond_topology)
            {
                order_bond.insert(index1, index2);
            } else {
                panic!(
                    "FATAL GORDER ERROR | OrderBonds::add | Could not find corresponding bond type for bond between atoms '{}' and '{}'. {}",
                    index1, index2, PANIC_MESSAGE
                );
            }
        }
    }
}

impl Add<OrderBonds> for OrderBonds {
    type Output = OrderBonds;

    #[inline(always)]
    fn add(self, rhs: OrderBonds) -> Self::Output {
        OrderBonds {
            bond_types: self
                .bond_types
                .into_iter()
                .zip(rhs.bond_types.into_iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<BondType>>(),
        }
    }
}

/// Checks that `min_index` is not higher than index of an atom involved in bonding.
/// Checks that there is no self-bonding.
/// Panics if any of these checks fails.
fn bonds_sanity_check(bonds: &HashSet<(usize, usize)>, min_index: usize) {
    for &(index1, index2) in bonds {
        for index in [index1, index2] {
            if index < min_index {
                panic!(
                    "FATAL GORDER ERROR | topology::bonds_sanity_check | Atom index '{}' is lower than minimum index '{}'. {}", 
                    index1, min_index, PANIC_MESSAGE
                );
            }
        }

        if index1 == index2 {
            panic!(
                "FATAL GORDER ERROR | topology::bonds_sanity_check | Bond between the same atom (index: '{}'). {}", 
                index1, PANIC_MESSAGE
            );
        }
    }
}

/// Get atoms corresponding to the provided absolute indices,
/// panicking with a suitable error message if these indices are out of range.
fn get_atoms_from_bond(system: &System, index1: usize, index2: usize) -> (&Atom, &Atom) {
    let atom1 = system
        .get_atom(index1)
        .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | topology::get_atoms_from_bond | Index '{}' does not correspond to an existing atom. {}", index1, PANIC_MESSAGE));

    let atom2 = system
        .get_atom(index2)
        .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | topology::get_atoms_from_bond | Index '{}' does not correspond to an existing atom. {}", index2, PANIC_MESSAGE));

    (atom1, atom2)
}

/// Collection of all atom types for which order parameters should be calculated.
/// In case of coarse-grained order parameters, this involves all specified atoms.
/// In case of atomistic order parameters, this only involves heavy atoms.
#[derive(Debug, Clone, PartialEq, Eq, MutGetters, Getters)]
pub(crate) struct OrderAtoms {
    /// Ordered by the increasing relative index.
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    atoms: Vec<AtomType>,
}

impl OrderAtoms {
    pub(super) fn new(system: &System, atoms: &[usize], minimum_index: usize) -> Self {
        let mut converted_atoms = atoms
            .into_iter()
            .map(|&x| AtomType::new(x - minimum_index, system.get_atom(x).expect(PANIC_MESSAGE)))
            .collect::<Vec<AtomType>>();

        converted_atoms.sort_by(|a, b| a.relative_index.cmp(&b.relative_index));

        OrderAtoms {
            atoms: converted_atoms,
        }
    }

    #[allow(unused)]
    #[inline(always)]
    pub(super) fn new_raw(atoms: Vec<AtomType>) -> Self {
        OrderAtoms { atoms }
    }
}

#[derive(Debug, Clone, Getters, MutGetters)]
pub(crate) struct BondType {
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    bond_topology: BondTopology,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    bonds: Vec<(usize, usize)>,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    total: Order,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    upper: Option<Order>,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    lower: Option<Order>,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    min_samples: usize,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    total_map: Option<Map>,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    upper_map: Option<Map>,
    #[getset(get = "pub(crate)", get_mut = "pub(crate)")]
    lower_map: Option<Map>,
}

impl BondType {
    /// Create a new bond for the calculation of order parameters.
    pub(crate) fn new(
        abs_index_1: usize,
        atom_1: &Atom,
        abs_index_2: usize,
        atom_2: &Atom,
        min_index: usize,
        classify_leaflets: bool,
        ordermap: Option<&OrderMap>,
        min_samples: usize,
        simbox: &SimBox,
    ) -> Result<Self, TopologyError> {
        let bond_topology = BondTopology::new(
            abs_index_1 - min_index,
            atom_1,
            abs_index_2 - min_index,
            atom_2,
        );
        let real_bond = if abs_index_1 < abs_index_2 {
            (abs_index_1, abs_index_2)
        } else {
            (abs_index_2, abs_index_1)
        };

        let optional_map = if let Some(map_params) = ordermap {
            Some(
                Map::new(map_params.to_owned(), simbox)
                    .map_err(|e| TopologyError::OrderMapError(e))?,
            )
        } else {
            None
        };

        let (leaflet_order, leaflet_map) = if classify_leaflets {
            (Some(Order::default()), optional_map.clone())
        } else {
            (None, None)
        };

        Ok(Self {
            bond_topology,
            bonds: vec![real_bond],
            total: Order::default(),
            upper: leaflet_order.clone(),
            lower: leaflet_order,
            min_samples: min_samples,
            total_map: optional_map,
            upper_map: leaflet_map.clone(),
            lower_map: leaflet_map,
        })
    }

    /// Insert new real bond to the current order bond.
    #[inline(always)]
    pub(super) fn insert(&mut self, abs_index_1: usize, abs_index_2: usize) {
        if abs_index_1 < abs_index_2 {
            self.bonds.push((abs_index_1, abs_index_2));
        } else {
            self.bonds.push((abs_index_2, abs_index_1));
        }
    }

    /// Get the first atom of this bond type from BondTopology.
    #[inline(always)]
    pub(crate) fn atom1(&self) -> &AtomType {
        self.bond_topology().atom1()
    }

    /// Get the second atom of this bond type from BondTopology.
    #[inline(always)]
    pub(crate) fn atom2(&self) -> &AtomType {
        self.bond_topology().atom2()
    }

    /// Get the relative index of the first atom of this bond type.
    #[inline(always)]
    pub(crate) fn atom1_index(&self) -> usize {
        self.atom1().relative_index
    }

    /// Get the relative index of the second atom of this bond type.
    #[inline(always)]
    pub(crate) fn atom2_index(&self) -> usize {
        self.atom2().relative_index
    }

    /// Does this bond involve a specific atom type?
    #[inline(always)]
    pub(crate) fn contains(&self, atom: &AtomType) -> bool {
        self.bond_topology().contains(atom)
    }

    /// Return the other atom involved in this bond.
    /// If the provided atom is not involved in the bond, return `None`.
    #[inline(always)]
    pub(crate) fn get_other_atom(&self, atom: &AtomType) -> Option<&AtomType> {
        self.bond_topology().get_other_atom(atom)
    }

    /// Calculate the current order parameter for this bond type.
    fn analyze_frame(
        &mut self,
        frame: &System,
        leaflet_classification: Option<&MoleculeLeafletClassification>,
        simbox: &SimBox,
        membrane_normal: &Vector3D,
    ) -> Result<(), AnalysisError> {
        for (molecule_index, (index1, index2)) in self.bonds.iter().enumerate() {
            let atom1 = unsafe { frame.get_atom_unchecked(*index1) };
            let atom2 = unsafe { frame.get_atom_unchecked(*index2) };

            let pos1 = atom1
                .get_position()
                .ok_or_else(|| AnalysisError::UndefinedPosition(atom1.get_index()))?;

            let pos2 = atom2
                .get_position()
                .ok_or_else(|| AnalysisError::UndefinedPosition(atom2.get_index()))?;

            let sch = super::calc_sch(pos1, pos2, simbox, membrane_normal);
            self.total += sch;

            if let Some(map) = self.total_map.as_mut() {
                map.add_order(sch, &((pos1 + pos2) / 2.0));
            }

            // get the assignment of molecule (assignment is performed earlier)
            if let Some(classifier) = leaflet_classification {
                match classifier
                    .get_assigned_leaflet(molecule_index)
                    .unwrap_or_else(|| panic!("FATAL GORDER ERROR | BondType::analyze_frame | Molecule with internal gorder index '{}' is not assigned into a leaflet.", molecule_index)) 
                {
                    Leaflet::Upper => {
                        *self.upper.as_mut().expect(PANIC_MESSAGE) += sch;
                        if let Some(map) = self.upper_map.as_mut() {
                            map.add_order(sch, &((pos1 + pos2) / 2.0));
                        }
                    }
                    Leaflet::Lower => {
                        *self.lower.as_mut().expect(PANIC_MESSAGE) += sch;
                        if let Some(map) = self.lower_map.as_mut() {
                            map.add_order(sch, &((pos1 + pos2) / 2.0));
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Calculate average order parameter for the full membrane and optionally for
    /// the upper and lower leaflet from the collected data.
    pub(crate) fn calc_order(&self) -> (f32, Option<f32>, Option<f32>) {
        let checked_min_samples = NonZeroUsize::new(self.min_samples).unwrap_or_else(|| {
            panic!(
                "FATAL GORDER ERROR | BondType::calc_order | 'min_samples' is ZERO. {}",
                PANIC_MESSAGE
            )
        });

        let total_order = self.total.calc_order(checked_min_samples);
        let upper_order = self
            .upper
            .as_ref()
            .map(|x| x.calc_order(checked_min_samples));
        let lower_order = self
            .lower
            .as_ref()
            .map(|x| x.calc_order(checked_min_samples));

        (total_order, upper_order, lower_order)
    }
}

impl Add<BondType> for BondType {
    type Output = BondType;

    #[inline]
    fn add(self, rhs: BondType) -> Self::Output {
        BondType {
            bond_topology: self.bond_topology,
            bonds: self.bonds,
            total: self.total + rhs.total,
            upper: merge_option_order(self.upper, rhs.upper),
            lower: merge_option_order(self.lower, rhs.lower),
            min_samples: self.min_samples,
            total_map: merge_option_maps(self.total_map, rhs.total_map),
            upper_map: merge_option_maps(self.upper_map, rhs.upper_map),
            lower_map: merge_option_maps(self.lower_map, rhs.lower_map),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Getters, MutGetters)]
pub(crate) struct BondTopology {
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    atom1: AtomType,
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    atom2: AtomType,
}

impl BondTopology {
    /// Construct a new BondTopology.
    /// The provided atoms can be in any order. In the constructed structure, the atom
    /// with the smaller index will always be the `atom1`.
    #[inline]
    pub(super) fn new(
        relative_index_1: usize,
        atom_1: &Atom,
        relative_index_2: usize,
        atom_2: &Atom,
    ) -> BondTopology {
        if relative_index_1 < relative_index_2 {
            BondTopology {
                atom1: AtomType::new(relative_index_1, atom_1),
                atom2: AtomType::new(relative_index_2, atom_2),
            }
        } else {
            BondTopology {
                atom1: AtomType::new(relative_index_2, atom_2),
                atom2: AtomType::new(relative_index_1, atom_1),
            }
        }
    }

    #[allow(unused)]
    #[inline(always)]
    pub(super) fn new_from_types(atom_type1: AtomType, atom_type2: AtomType) -> BondTopology {
        BondTopology {
            atom1: atom_type1,
            atom2: atom_type2,
        }
    }

    /// Does this bond type involve the provided atom type?
    #[inline(always)]
    fn contains(&self, atom: &AtomType) -> bool {
        self.atom1() == atom || self.atom2() == atom
    }

    /// Return the other atom involved in this bond.
    /// If the provided atom is not involved in the bond, return `None`.
    #[inline(always)]
    fn get_other_atom(&self, atom: &AtomType) -> Option<&AtomType> {
        if self.atom1() == atom {
            Some(self.atom2())
        } else if self.atom2() == atom {
            Some(self.atom1())
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Getters, MutGetters, CopyGetters)]
pub(crate) struct AtomType {
    #[getset(get_copy = "pub(crate)", get_mut = "pub(super)")]
    relative_index: usize,
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    residue_name: String,
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    atom_name: String,
}

impl AtomType {
    #[inline(always)]
    pub(super) fn new(relative_index: usize, atom: &Atom) -> AtomType {
        AtomType {
            relative_index,
            residue_name: atom.get_residue_name().to_owned(),
            atom_name: atom.get_atom_name().to_owned(),
        }
    }

    #[allow(unused)]
    #[inline(always)]
    pub(crate) fn new_raw(relative_index: usize, residue_name: &str, atom_name: &str) -> AtomType {
        AtomType {
            relative_index,
            residue_name: residue_name.to_owned(),
            atom_name: atom_name.to_owned(),
        }
    }
}

impl fmt::Display for AtomType {
    #[inline(always)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}-{}-{}",
            self.residue_name(),
            self.atom_name(),
            self.relative_index()
        )
    }
}

#[derive(Debug, Clone, Default, CopyGetters)]
pub(crate) struct Order {
    #[getset(get_copy = "pub(super)")]
    order: f32,
    #[getset(get_copy = "pub(super)")]
    n_samples: usize,
}

impl Order {
    #[allow(unused)]
    #[inline(always)]
    pub(super) fn new(order: f32, n_samples: usize) -> Order {
        Order { order, n_samples }
    }

    /// Calculate average order from the collected data.
    ///
    /// Return `f32::NAN` if the number of samples is lower than the required minimal number.
    #[inline(always)]
    pub(crate) fn calc_order(&self, min_samples: NonZeroUsize) -> f32 {
        if self.n_samples < min_samples.into() {
            f32::NAN
        } else {
            self.order / self.n_samples as f32
        }
    }
}

impl Add<Order> for Order {
    type Output = Order;

    #[inline(always)]
    fn add(self, rhs: Order) -> Self::Output {
        Order {
            order: self.order + rhs.order,
            n_samples: self.n_samples + rhs.n_samples,
        }
    }
}

impl AddAssign<f32> for Order {
    #[inline(always)]
    fn add_assign(&mut self, rhs: f32) {
        self.order += rhs;
        self.n_samples += 1;
    }
}

/// Helper function for merging optional Orders.
#[inline]
fn merge_option_order(lhs: Option<Order>, rhs: Option<Order>) -> Option<Order> {
    match (lhs, rhs) {
        (Some(x), Some(y)) => Some(x + y),
        (None, None) => None,
        (Some(_), None) | (None, Some(_)) => panic!(
            "FATAL GORDER ERROR | merge_option_order | Inconsistent option value. {}",
            PANIC_MESSAGE
        ),
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::input::{ordermap::Plane, GridSpan};

    use super::*;

    #[test]
    fn test_bond_type_new() {
        let atom1 = Atom::new(1, "POPE", 1, "N");
        let atom2 = Atom::new(1, "POPE", 6, "HN");

        let bond1 = BondTopology::new(0, &atom1, 5, &atom2);
        let bond2 = BondTopology::new(5, &atom2, 0, &atom1);

        assert_eq!(bond1, bond2);
    }

    #[test]
    fn test_bond_new() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let bond1 = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();
        let bond2 = BondType::new(
            460,
            &atom2,
            455,
            &atom1,
            455,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        assert_eq!(bond1.bond_topology, bond2.bond_topology);
        assert_eq!(bond1.bonds.len(), 1);
        assert_eq!(bond2.bonds.len(), 1);
        assert_eq!(bond1.bonds[0], bond2.bonds[0]);
    }

    #[test]
    fn test_bond_new_with_leaflet_classification() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let bond = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            true,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        assert!(bond.lower.is_some());
        assert!(bond.upper.is_some());
    }

    #[test]
    fn test_bond_new_with_ordermap() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let ordermap_params = OrderMap::new()
            .dim_x(GridSpan::Auto)
            .dim_y(GridSpan::Auto)
            .output_directory(".")
            .plane(Plane::XY)
            .build()
            .unwrap();
        let bond = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            false,
            Some(&ordermap_params),
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        assert!(bond.lower.is_none());
        assert!(bond.upper.is_none());
        assert!(bond.total_map.is_some());
        assert!(bond.lower_map.is_none());
        assert!(bond.upper_map.is_none());
    }

    #[test]
    fn test_bond_new_with_leaflet_classification_and_ordermap() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let ordermap_params = OrderMap::new()
            .dim_x(GridSpan::Auto)
            .dim_y(GridSpan::Auto)
            .output_directory(".")
            .plane(Plane::XY)
            .build()
            .unwrap();
        let bond = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            true,
            Some(&ordermap_params),
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        assert!(bond.lower.is_some());
        assert!(bond.upper.is_some());
        assert!(bond.total_map.is_some());
        assert!(bond.lower_map.is_some());
        assert!(bond.upper_map.is_some());
    }

    #[test]
    fn test_bond_add() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let mut bond = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        bond.insert(1354, 1359);
        bond.insert(1676, 1671);

        assert_eq!(bond.bonds.len(), 3);
        assert_eq!(bond.bonds[0], (455, 460));
        assert_eq!(bond.bonds[1], (1354, 1359));
        assert_eq!(bond.bonds[2], (1671, 1676));
    }

    #[test]
    fn test_bond_calc_order_basic() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let mut bond = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        bond.total.n_samples = 17;
        bond.total.order = 3.978;

        let (total, upper, lower) = bond.calc_order();

        assert_relative_eq!(total, 0.234);
        assert!(upper.is_none());
        assert!(lower.is_none());
    }

    #[test]
    fn test_bond_calc_order_leaflets() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let mut bond = BondType::new(
            455,
            &atom1,
            460,
            &atom2,
            455,
            true,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        )
        .unwrap();

        bond.total.n_samples = 17;
        bond.total.order = 3.978;

        bond.upper.as_mut().unwrap().n_samples = 8;
        bond.upper.as_mut().unwrap().order = 1.976;

        bond.lower.as_mut().unwrap().n_samples = 9;
        bond.lower.as_mut().unwrap().order = 1.989;

        let (total, upper, lower) = bond.calc_order();

        assert_relative_eq!(total, 0.234);
        assert_relative_eq!(upper.unwrap(), 0.247);
        assert_relative_eq!(lower.unwrap(), 0.221);
    }

    #[test]
    fn order_atoms_new() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let atoms = [133, 127, 163, 156, 145];
        let order_atoms = OrderAtoms::new(&system, &atoms, 125);

        let expected_atoms = [
            AtomType::new(2, system.get_atom(127).unwrap()),
            AtomType::new(8, system.get_atom(133).unwrap()),
            AtomType::new(20, system.get_atom(145).unwrap()),
            AtomType::new(31, system.get_atom(156).unwrap()),
            AtomType::new(38, system.get_atom(163).unwrap()),
        ];

        for (atom, expected) in order_atoms.atoms.iter().zip(expected_atoms.iter()) {
            assert_eq!(atom, expected);
        }
    }

    fn expected_bonds(system: &System) -> [BondTopology; 5] {
        [
            BondTopology::new(
                44,
                system.get_atom(169).unwrap(),
                45,
                system.get_atom(170).unwrap(),
            ),
            BondTopology::new(
                44,
                system.get_atom(169).unwrap(),
                46,
                system.get_atom(171).unwrap(),
            ),
            BondTopology::new(
                88,
                system.get_atom(213).unwrap(),
                89,
                system.get_atom(214).unwrap(),
            ),
            BondTopology::new(
                88,
                system.get_atom(213).unwrap(),
                90,
                system.get_atom(215).unwrap(),
            ),
            BondTopology::new(
                121,
                system.get_atom(246).unwrap(),
                122,
                system.get_atom(247).unwrap(),
            ),
        ]
    }

    #[test]
    fn order_bonds_new() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let bonds = [(169, 170), (169, 171), (213, 214), (213, 215), (246, 247)];
        let bonds_set = HashSet::from(bonds.clone());
        let order_bonds = OrderBonds::new(&system, &bonds_set, 125, false, None, 1).unwrap();

        let expected_bonds = expected_bonds(&system);

        // bonds can be in any order inside `order_bonds`
        for bond in order_bonds.bond_types.iter() {
            let expected = expected_bonds
                .iter()
                .enumerate()
                .find(|(_, expected)| &bond.bond_topology == *expected);

            if let Some((i, _)) = expected {
                assert_eq!(bond.bonds.len(), 1);
                assert_eq!(bond.bonds[0], bonds[i]);
            } else {
                panic!("Expected bond not found.");
            }
        }
    }

    #[test]
    fn order_bonds_add() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let bonds = [(169, 170), (169, 171), (213, 214), (213, 215), (246, 247)];
        let bonds_set = HashSet::from(bonds.clone());
        let mut order_bonds = OrderBonds::new(&system, &bonds_set, 125, false, None, 1).unwrap();

        let new_bonds = [(919, 920), (919, 921), (963, 964), (963, 965), (996, 997)];
        let new_bonds_set = HashSet::from(new_bonds.clone());
        order_bonds.insert(&system, &new_bonds_set, 875);

        let expected_bonds = expected_bonds(&system);
        for bond in order_bonds.bond_types.iter() {
            let expected = expected_bonds
                .iter()
                .enumerate()
                .find(|(_, expected)| &bond.bond_topology == *expected);

            if let Some((i, _)) = expected {
                assert_eq!(bond.bonds.len(), 2);
                assert_eq!(bond.bonds[0], bonds[i]);
                assert_eq!(bond.bonds[1], new_bonds[i]);
            } else {
                panic!("Expected bond not found.");
            }
        }
    }

    #[test]
    fn molecule_topology_new() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let bonds = [(169, 170), (169, 171), (213, 214), (213, 215), (246, 247)];
        let bonds_set = HashSet::from(bonds.clone());
        let topology = MoleculeTopology::new(&system, &bonds_set, 125);

        let expected_bonds = expected_bonds(&system);

        // bonds can be in any order inside `order_bonds`
        for bond in topology.bonds.iter() {
            if !expected_bonds.iter().any(|expected| bond == expected) {
                panic!("Expected bond not found.")
            }
        }
    }
}
