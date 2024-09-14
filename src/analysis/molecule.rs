// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::ops::{Add, AddAssign};

use getset::{CopyGetters, Getters, MutGetters};
use groan_rs::{
    prelude::{Atom, SimBox, Vector3D},
    structures::group::Group,
    system::System,
};
use hashbrown::HashSet;
use serde::{Deserialize, Serialize};

use crate::{
    errors::{AnalysisError, TopologyError},
    OrderMap, PANIC_MESSAGE,
};

use super::{
    leaflets::{Leaflet, LeafletClassifier, MoleculeLeafletClassification},
    ordermap::{merge_option_maps, Map},
};

#[derive(Debug, Clone, Getters, MutGetters)]
pub(super) struct MoleculeType {
    #[getset(get = "pub(super)")]
    name: String,
    #[getset(get = "pub(super)")]
    topology: MoleculeTopology,
    #[getset(get = "pub(super)")]
    order_bonds: OrderBonds,
    #[getset(get = "pub(super)")]
    order_atoms: OrderAtoms,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    leaflet_classification: Option<MoleculeLeafletClassification>,
}

impl MoleculeType {
    pub(super) fn new(
        system: &System,
        name: &str,
        topology: &MoleculeTopology,
        order_bonds: &HashSet<(usize, usize)>,
        order_atoms: &[usize],
        min_index: usize,
        leaflet_classification: Option<MoleculeLeafletClassification>,
        ordermap_params: Option<&OrderMap>,
    ) -> Self {
        Self {
            name: name.to_owned(),
            topology: topology.to_owned(),
            order_bonds: OrderBonds::new(
                system,
                order_bonds,
                min_index,
                leaflet_classification.is_some(),
                ordermap_params,
            ),
            order_atoms: OrderAtoms::new(system, order_atoms, min_index),
            leaflet_classification,
        }
    }

    /// Insert new bond instances to the molecule.
    pub(super) fn insert(
        &mut self,
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        atoms: Group,
    ) -> Result<(), TopologyError> {
        self.order_bonds
            .insert(system, bonds, atoms.get_atoms().first().unwrap());

        if let Some(classifier) = self.leaflet_classification.as_mut() {
            classifier.insert(&atoms, system)?;
        }

        Ok(())
    }

    /// Calculate order parameters for bonds of a single molecule type from a single simulation frame.
    pub(super) fn analyze_frame(
        &mut self,
        frame: &System,
        simbox: &SimBox,
        membrane_normal: &Vector3D,
    ) -> Result<(), AnalysisError> {
        for bond_type in self.order_bonds.bonds.iter_mut() {
            for (molecule_index, (index1, index2)) in bond_type.bonds.iter().enumerate() {
                let atom1 = unsafe { frame.get_atom_unchecked_as_ref(*index1) };
                let atom2 = unsafe { frame.get_atom_unchecked_as_ref(*index2) };

                let sch = super::calc_sch(atom1, atom2, simbox, membrane_normal)?;
                bond_type.total += sch;

                // assign molecule to leaflet
                if let Some(classifier) = &self.leaflet_classification {
                    match classifier.assign_to_leaflet(frame, molecule_index)? {
                        Leaflet::Upper => *bond_type.upper.as_mut().expect(PANIC_MESSAGE) += sch,
                        Leaflet::Lower => *bond_type.lower.as_mut().expect(PANIC_MESSAGE) += sch,
                    }
                }
            }
        }

        Ok(())
    }
}

impl Add<MoleculeType> for MoleculeType {
    type Output = Self;

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
    pub(super) bonds: HashSet<BondType>,
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
            let bond = BondType::new(index1 - min_index, atom1, index2 - min_index, atom2);

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
pub(super) struct OrderBonds {
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    bonds: Vec<Bond>,
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
    ) -> OrderBonds {
        bonds_sanity_check(bonds, min_index);

        let mut order_bonds = Vec::new();
        for &(index1, index2) in bonds.iter() {
            let (atom1, atom2) = get_atoms_from_bond(system, index1, index2);
            let bond = Bond::new(
                index1,
                atom1,
                index2,
                atom2,
                min_index,
                classify_leaflets,
                ordermap,
            );
            order_bonds.push(bond)
        }

        OrderBonds { bonds: order_bonds }
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
            let bond_type = BondType::new(index1 - min_index, atom1, index2 - min_index, atom2);

            if let Some(order_bond) = self
                .bonds
                .iter_mut()
                .find(|order_bond| order_bond.bond_type == bond_type)
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

    fn add(self, rhs: OrderBonds) -> Self::Output {
        OrderBonds {
            bonds: self
                .bonds
                .into_iter()
                .zip(rhs.bonds.into_iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<Bond>>(),
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
        .get_atom_as_ref(index1)
        .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | topology::get_atoms_from_bond | Index '{}' does not correspond to an existing atom. {}", index1, PANIC_MESSAGE));

    let atom2 = system
        .get_atom_as_ref(index2)
        .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | topology::get_atoms_from_bond | Index '{}' does not correspond to an existing atom. {}", index2, PANIC_MESSAGE));

    (atom1, atom2)
}

/// Collection of all atom types for which order parameters should be calculated.
/// In case of coarse-grained order parameters, this involves all specified atoms.
/// In case of atomistic order parameters, this only involves heavy atoms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct OrderAtoms {
    /// Ordered by the increasing relative index.
    pub(super) atoms: Vec<AtomType>,
}

impl OrderAtoms {
    pub(super) fn new(system: &System, atoms: &[usize], minimum_index: usize) -> OrderAtoms {
        let mut converted_atoms = atoms
            .into_iter()
            .map(|&x| {
                AtomType::new(
                    x - minimum_index,
                    system.get_atom_as_ref(x).expect(PANIC_MESSAGE),
                )
            })
            .collect::<Vec<AtomType>>();

        converted_atoms.sort_by(|a, b| a.relative_index.cmp(&b.relative_index));

        OrderAtoms {
            atoms: converted_atoms,
        }
    }
}

#[derive(Debug, Clone, Getters, MutGetters)]
pub(crate) struct Bond {
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    bond_type: BondType,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    bonds: Vec<(usize, usize)>,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    total: Order,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    upper: Option<Order>,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    lower: Option<Order>,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    total_map: Option<Map>,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    upper_map: Option<Map>,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    lower_map: Option<Map>,
}

impl Bond {
    /// Create a new bond for the calculation of order parameters.
    pub(super) fn new(
        abs_index_1: usize,
        atom_1: &Atom,
        abs_index_2: usize,
        atom_2: &Atom,
        min_index: usize,
        classify_leaflets: bool,
        ordermap: Option<&OrderMap>,
    ) -> Bond {
        let bond_type = BondType::new(
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
            Some(Map::new(map_params.to_owned()))
        } else {
            None
        };

        let (leaflet_order, leaflet_map) = if classify_leaflets {
            (Some(Order::default()), optional_map.clone())
        } else {
            (None, None)
        };

        Bond {
            bond_type,
            bonds: vec![real_bond],
            total: Order::default(),
            upper: leaflet_order.clone(),
            lower: leaflet_order,
            total_map: optional_map,
            upper_map: leaflet_map.clone(),
            lower_map: leaflet_map,
        }
    }

    /// Insert new real bond to the current order bond.
    pub(super) fn insert(&mut self, abs_index_1: usize, abs_index_2: usize) {
        if abs_index_1 < abs_index_2 {
            self.bonds.push((abs_index_1, abs_index_2));
        } else {
            self.bonds.push((abs_index_2, abs_index_1));
        }
    }
}

impl Add<Bond> for Bond {
    type Output = Bond;
    fn add(self, rhs: Bond) -> Self::Output {
        Bond {
            bond_type: self.bond_type,
            bonds: self.bonds,
            total: self.total + rhs.total,
            upper: merge_option_order(self.upper, rhs.upper),
            lower: merge_option_order(self.lower, rhs.lower),
            total_map: merge_option_maps(self.total_map, rhs.total_map),
            upper_map: merge_option_maps(self.upper_map, rhs.upper_map),
            lower_map: merge_option_maps(self.lower_map, rhs.lower_map),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Getters, MutGetters)]
pub(super) struct BondType {
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    atom1: AtomType,
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    atom2: AtomType,
}

impl BondType {
    /// Construct a new BondType.
    /// The provided atoms can be in any order. In the constructed structure, the atom
    /// with the smaller index will always be the `atom1`.
    pub(super) fn new(
        relative_index_1: usize,
        atom_1: &Atom,
        relative_index_2: usize,
        atom_2: &Atom,
    ) -> BondType {
        if relative_index_1 < relative_index_2 {
            BondType {
                atom1: AtomType::new(relative_index_1, atom_1),
                atom2: AtomType::new(relative_index_2, atom_2),
            }
        } else {
            BondType {
                atom1: AtomType::new(relative_index_2, atom_2),
                atom2: AtomType::new(relative_index_1, atom_1),
            }
        }
    }

    #[allow(unused)]
    pub(super) fn new_from_types(atom_type1: AtomType, atom_type2: AtomType) -> BondType {
        BondType {
            atom1: atom_type1,
            atom2: atom_type2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Getters, MutGetters)]
pub(crate) struct AtomType {
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    relative_index: usize,
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    residue_name: String,
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    atom_name: String,
}

impl AtomType {
    pub(super) fn new(relative_index: usize, atom: &Atom) -> AtomType {
        AtomType {
            relative_index,
            residue_name: atom.get_residue_name().to_owned(),
            atom_name: atom.get_atom_name().to_owned(),
        }
    }

    #[allow(unused)]
    pub(crate) fn new_raw(relative_index: usize, residue_name: &str, atom_name: &str) -> AtomType {
        AtomType {
            relative_index,
            residue_name: residue_name.to_owned(),
            atom_name: atom_name.to_owned(),
        }
    }
}

#[derive(Debug, Clone, Default, CopyGetters)]
pub(super) struct Order {
    #[getset(get_copy = "pub(super)")]
    order: f32,
    #[getset(get_copy = "pub(super)")]
    n_samples: usize,
}

impl Order {
    #[allow(unused)]
    pub(super) fn new(order: f32, n_samples: usize) -> Order {
        Order { order, n_samples }
    }
}

impl Add<Order> for Order {
    type Output = Order;
    fn add(self, rhs: Order) -> Self::Output {
        Order {
            order: self.order + rhs.order,
            n_samples: self.n_samples + rhs.n_samples,
        }
    }
}

impl AddAssign<f32> for Order {
    fn add_assign(&mut self, rhs: f32) {
        self.order += rhs;
        self.n_samples += 1;
    }
}

/// Helper function for merging optional Orders.
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
    use crate::input::GridSpan;

    use super::*;

    #[test]
    fn test_bond_type_new() {
        let atom1 = Atom::new(1, "POPE", 1, "N");
        let atom2 = Atom::new(1, "POPE", 6, "HN");

        let bond1 = BondType::new(0, &atom1, 5, &atom2);
        let bond2 = BondType::new(5, &atom2, 0, &atom1);

        assert_eq!(bond1, bond2);
    }

    #[test]
    fn test_bond_new() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let bond1 = Bond::new(455, &atom1, 460, &atom2, 455, false, None);
        let bond2 = Bond::new(460, &atom2, 455, &atom1, 455, false, None);

        assert_eq!(bond1.bond_type, bond2.bond_type);
        assert_eq!(bond1.bonds.len(), 1);
        assert_eq!(bond2.bonds.len(), 1);
        assert_eq!(bond1.bonds[0], bond2.bonds[0]);
    }

    #[test]
    fn test_bond_new_with_leaflet_classification() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let bond = Bond::new(455, &atom1, 460, &atom2, 455, true, None);

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
            .build()
            .unwrap();
        let bond = Bond::new(455, &atom1, 460, &atom2, 455, false, Some(&ordermap_params));

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
            .build()
            .unwrap();
        let bond = Bond::new(455, &atom1, 460, &atom2, 455, true, Some(&ordermap_params));

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

        let mut bond = Bond::new(455, &atom1, 460, &atom2, 455, false, None);

        bond.insert(1354, 1359);
        bond.insert(1676, 1671);

        assert_eq!(bond.bonds.len(), 3);
        assert_eq!(bond.bonds[0], (455, 460));
        assert_eq!(bond.bonds[1], (1354, 1359));
        assert_eq!(bond.bonds[2], (1671, 1676));
    }

    #[test]
    fn order_atoms_new() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let atoms = [133, 127, 163, 156, 145];
        let order_atoms = OrderAtoms::new(&system, &atoms, 125);

        let expected_atoms = [
            AtomType::new(2, system.get_atom_as_ref(127).unwrap()),
            AtomType::new(8, system.get_atom_as_ref(133).unwrap()),
            AtomType::new(20, system.get_atom_as_ref(145).unwrap()),
            AtomType::new(31, system.get_atom_as_ref(156).unwrap()),
            AtomType::new(38, system.get_atom_as_ref(163).unwrap()),
        ];

        for (atom, expected) in order_atoms.atoms.iter().zip(expected_atoms.iter()) {
            assert_eq!(atom, expected);
        }
    }

    fn expected_bonds(system: &System) -> [BondType; 5] {
        [
            BondType::new(
                44,
                system.get_atom_as_ref(169).unwrap(),
                45,
                system.get_atom_as_ref(170).unwrap(),
            ),
            BondType::new(
                44,
                system.get_atom_as_ref(169).unwrap(),
                46,
                system.get_atom_as_ref(171).unwrap(),
            ),
            BondType::new(
                88,
                system.get_atom_as_ref(213).unwrap(),
                89,
                system.get_atom_as_ref(214).unwrap(),
            ),
            BondType::new(
                88,
                system.get_atom_as_ref(213).unwrap(),
                90,
                system.get_atom_as_ref(215).unwrap(),
            ),
            BondType::new(
                121,
                system.get_atom_as_ref(246).unwrap(),
                122,
                system.get_atom_as_ref(247).unwrap(),
            ),
        ]
    }

    #[test]
    fn order_bonds_new() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let bonds = [(169, 170), (169, 171), (213, 214), (213, 215), (246, 247)];
        let bonds_set = HashSet::from(bonds.clone());
        let order_bonds = OrderBonds::new(&system, &bonds_set, 125, false, None);

        let expected_bonds = expected_bonds(&system);

        // bonds can be in any order inside `order_bonds`
        for bond in order_bonds.bonds.iter() {
            let expected = expected_bonds
                .iter()
                .enumerate()
                .find(|(_, expected)| &bond.bond_type == *expected);

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
        let mut order_bonds = OrderBonds::new(&system, &bonds_set, 125, false, None);

        let new_bonds = [(919, 920), (919, 921), (963, 964), (963, 965), (996, 997)];
        let new_bonds_set = HashSet::from(new_bonds.clone());
        order_bonds.insert(&system, &new_bonds_set, 875);

        let expected_bonds = expected_bonds(&system);
        for bond in order_bonds.bonds.iter() {
            let expected = expected_bonds
                .iter()
                .enumerate()
                .find(|(_, expected)| &bond.bond_type == *expected);

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
