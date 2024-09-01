// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use hashbrown::HashSet;

use getset::Getters;
use groan_rs::{
    prelude::{Atom, GridMap},
    system::System,
};

use crate::{auxiliary::PANIC_MESSAGE, ordermap::OrderMap};

#[derive(Getters)]
pub(crate) struct Molecule {
    #[getset(get = "pub(crate)")]
    name: String,
    #[getset(get = "pub(crate)")]
    topology: MoleculeTopology,
    #[getset(get = "pub(crate)")]
    order_bonds: OrderBonds,
    #[getset(get = "pub(crate)")]
    order_atoms: OrderAtoms,
}

impl Molecule {
    pub(crate) fn new(
        system: &System,
        name: &str,
        topology: &MoleculeTopology,
        order_bonds: &HashSet<(usize, usize)>,
        order_atoms: &[usize],
        min_index: usize,
    ) -> Molecule {
        Molecule {
            name: name.to_owned(),
            topology: topology.to_owned(),
            order_bonds: OrderBonds::new(system, order_bonds, min_index),
            order_atoms: OrderAtoms::new(system, order_atoms, min_index),
        }
    }

    /// Add new bond instances to the molecule.
    pub(crate) fn add(
        &mut self,
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        min_index: usize,
    ) {
        self.order_bonds.add(system, bonds, min_index);
    }
}

/// Collection of all bond types in a molecule describing its topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MoleculeTopology {
    pub(crate) bonds: HashSet<BondType>,
}

impl MoleculeTopology {
    /// Create new molecule topology from a set of bonds (absolute indices) and the minimum index in the molecule.
    ///
    /// ## Panics
    /// - Panics if the `min_index` is higher than any index inside the `bonds` set.
    /// - Panics if there is a bond connecting the same atom (e.g. 14-14) in the `bonds` set.
    /// - Panics if an index in the `bonds` set does not correspond to an existing atom.
    pub(crate) fn new(
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
                panic!("FATAL GORDER ERROR | MoleculeTopology::new | Bond between atoms '{}' and '{}' defined multiple times. {}", index1, index2, PANIC_MESSAGE);
            }
        }

        MoleculeTopology {
            bonds: converted_bonds,
        }
    }
}

/// Collection of all bonds for which the order parameters should be calculated.
pub(crate) struct OrderBonds {
    pub(crate) bonds: Vec<Bond>,
}

impl OrderBonds {
    /// Create a new `OrderBonds` structure from a set of bonds (absolute indices) and the minimum index of the molecule.
    ///
    /// ## Panics
    /// - Panics if the `min_index` is higher than any index inside the `bonds` set.
    /// - Panics if there is a bond connecting the same atom (e.g. 14-14) in the `bonds` set.
    /// - Panics if an index in the `bonds` set does not correspond to an existing atom.
    pub(crate) fn new(
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        min_index: usize,
    ) -> OrderBonds {
        bonds_sanity_check(bonds, min_index);

        let mut order_bonds = Vec::new();
        for &(index1, index2) in bonds.iter() {
            let (atom1, atom2) = get_atoms_from_bond(system, index1, index2);
            let bond = Bond::new(index1, atom1, index2, atom2, min_index);
            order_bonds.push(bond)
        }

        OrderBonds { bonds: order_bonds }
    }

    /// Add new real bonds to already constructed order bonds.
    pub(crate) fn add(
        &mut self,
        system: &System,
        bonds: &HashSet<(usize, usize)>,
        min_index: usize,
    ) {
        bonds_sanity_check(bonds, min_index);

        for &(index1, index2) in bonds.iter() {
            let (atom1, atom2) = get_atoms_from_bond(system, index1, index2);
            let bond_type = BondType::new(index1 - min_index, atom1, index2 - min_index, atom2);
            for order_bond in self.bonds.iter_mut() {
                if order_bond.bond_type == bond_type {
                    order_bond.add(index1, index2)
                }
            }
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
                panic!("FATAL GORDER ERROR | molecule::bonds_sanity_check | Atom index '{}' is lower than minimum index '{}'. {}", index1, min_index, PANIC_MESSAGE);
            }
        }

        if index1 == index2 {
            panic!("FATAL GORDER ERROR | molecule::bonds_sanity_check | Bond between the same atom (index: '{}'). {}", index1, PANIC_MESSAGE);
        }
    }
}

/// Get atoms corresponding to the provided absolute indices,
/// panicking with a suitable error message if these indices are out of range.
fn get_atoms_from_bond(system: &System, index1: usize, index2: usize) -> (&Atom, &Atom) {
    let atom1 = system
        .get_atom_as_ref(index1)
        .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | molecule::get_atoms_from_bond | Index '{}' does not correspond to an existing atom. {}", index1, PANIC_MESSAGE));

    let atom2 = system
        .get_atom_as_ref(index2)
        .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | molecule::get_atoms_from_bond | Index '{}' does not correspond to an existing atom. {}", index2, PANIC_MESSAGE));

    (atom1, atom2)
}

/// Collection of all atom types for which order parameters should be calculated.
/// In case of coarse-grained order parameters, this involves all specified atoms.
/// In case of atomistic order parameters, this only involves heavy atoms.
#[derive(Debug, Clone)]
pub(crate) struct OrderAtoms {
    pub(crate) atoms: Vec<AtomType>,
}

impl OrderAtoms {
    pub(crate) fn new(system: &System, atoms: &[usize], minimum_index: usize) -> OrderAtoms {
        OrderAtoms {
            atoms: atoms
                .into_iter()
                .map(|&x| {
                    AtomType::new(
                        x - minimum_index,
                        system.get_atom_as_ref(x).expect(PANIC_MESSAGE),
                    )
                })
                .collect(),
        }
    }
}

pub(crate) struct Bond {
    pub(crate) bond_type: BondType,
    pub(crate) bonds: Vec<(usize, usize)>,
    pub(crate) total: Order,
    pub(crate) upper: Option<Order>,
    pub(crate) lower: Option<Order>,
    pub(crate) total_map: Option<Map>,
    pub(crate) upper_map: Option<Map>,
    pub(crate) lower_map: Option<Map>,
}

impl Bond {
    /// Create a new bond for the calculation of order parameters.
    pub(crate) fn new(
        abs_index_1: usize,
        atom_1: &Atom,
        abs_index_2: usize,
        atom_2: &Atom,
        min_index: usize,
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

        Bond {
            bond_type,
            bonds: vec![real_bond],
            total: Order::default(),
            upper: None,
            lower: None,
            total_map: None,
            upper_map: None,
            lower_map: None,
        }
    }

    /// Add new real bond to the current order bond.
    pub(crate) fn add(&mut self, abs_index_1: usize, abs_index_2: usize) {
        if abs_index_1 < abs_index_2 {
            self.bonds.push((abs_index_1, abs_index_2));
        } else {
            self.bonds.push((abs_index_2, abs_index_1));
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct BondType {
    pub(crate) atom1: AtomType,
    pub(crate) atom2: AtomType,
}

impl BondType {
    /// Construct a new BondType.
    /// The provided atoms can be in any order. In the constructed structure, the atom
    /// with the smaller index will always be the `atom1`.
    pub(crate) fn new(
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
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct AtomType {
    pub(crate) relative_index: usize,
    pub(crate) residue_name: String,
    pub(crate) atom_name: String,
}

impl AtomType {
    pub(crate) fn new(relative_index: usize, atom: &Atom) -> AtomType {
        AtomType {
            relative_index,
            residue_name: atom.get_residue_name().to_owned(),
            atom_name: atom.get_atom_name().to_owned(),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct Order {
    order: f32,
    n_samples: usize,
}

pub(crate) struct Map {
    params: OrderMap,
    values: GridMap<f32, f32, Box<dyn Fn(&f32) -> f32>>,
    samples: GridMap<usize, usize, Box<dyn Fn(&usize) -> usize>>,
}

#[cfg(test)]
mod tests {
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

        let bond1 = Bond::new(455, &atom1, 460, &atom2, 455);
        let bond2 = Bond::new(460, &atom2, 455, &atom1, 455);

        assert_eq!(bond1.bond_type, bond2.bond_type);
        assert_eq!(bond1.bonds.len(), 1);
        assert_eq!(bond2.bonds.len(), 1);
        assert_eq!(bond1.bonds[0], bond2.bonds[0]);
    }

    #[test]
    fn test_bond_add() {
        let atom1 = Atom::new(17, "POPE", 456, "N");
        let atom2 = Atom::new(17, "POPE", 461, "HN");

        let mut bond = Bond::new(455, &atom1, 460, &atom2, 455);

        bond.add(1354, 1359);
        bond.add(1676, 1671);

        assert_eq!(bond.bonds.len(), 3);
        assert_eq!(bond.bonds[0], (455, 460));
        assert_eq!(bond.bonds[1], (1354, 1359));
        assert_eq!(bond.bonds[2], (1671, 1676));
    }
}
