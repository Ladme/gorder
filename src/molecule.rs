// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::collections::HashSet;

use groan_rs::prelude::{Atom, GridMap};

use crate::{auxiliary::PANIC_MESSAGE, ordermap::OrderMap};

pub(crate) struct Molecule {
    name: String,
    all_bonds: HashSet<BondType>,
    all_atom_types: Vec<AtomType>,
    order_bonds: Vec<Bond>,
    heavy_atom_types: Vec<AtomType>,
}

impl Molecule {
    /// Construct a new valid `Molecule` from a set of bonds (with absolute indices of atoms)
    /// and a vector of a pairs of absolute atom indices and atom types.
    ///
    /// Assumes that the `atom_types` vector contains at least one element.
    pub(crate) fn new(
        name: &str,
        bonds: HashSet<(usize, usize)>,
        mut atom_types: Vec<(usize, AtomType)>,
    ) -> Molecule {
        atom_types.sort_by_key(|&(index, _)| index);
        println!("{:?}", atom_types);
        let minimum_index = atom_types.get(0).expect(PANIC_MESSAGE).0;

        let final_atom_types: Vec<AtomType> = atom_types.iter().map(|x| x.1.clone()).collect();
        let final_bonds: HashSet<BondType> = HashSet::from_iter(
            bonds
                .into_iter()
                .map(|(i, j)| BondType::new(i - minimum_index, j - minimum_index)),
        );

        println!("{:?}", final_atom_types);
        println!("{:?}", final_bonds);

        Molecule {
            name: name.to_owned(),
            all_bonds: final_bonds,
            all_atom_types: final_atom_types,
            order_bonds: vec![],
            heavy_atom_types: vec![],
        }
    }
}

pub(crate) struct Bond {
    bond_type: BondType,
    bonds: Vec<(usize, usize)>,
    total: Order,
    upper: Option<Order>,
    lower: Option<Order>,
    total_map: Option<Map>,
    upper_map: Option<Map>,
    lower_map: Option<Map>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct BondType {
    relative_index_1: usize,
    relative_index_2: usize,
}

impl BondType {
    /// Construct a new BondType.
    /// The provided atom indices can be in any order. In the constructed structure, `relative_index_1` is
    /// always smaller than `relative_index_2`.
    ///
    /// ## Panics
    /// Panics if `index1` equals `index2`.
    fn new(index1: usize, index2: usize) -> BondType {
        if index1 == index2 {
            panic!("FATAL GORDER ERROR | BondType::new | `index1` ({}) and `index2` ({}) are the same. {}", index1, index2, PANIC_MESSAGE);
        }

        if index1 < index2 {
            BondType {
                relative_index_1: index1,
                relative_index_2: index2,
            }
        } else {
            BondType {
                relative_index_1: index2,
                relative_index_2: index1,
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct AtomType {
    residue_name: String,
    atom_name: String,
}

impl From<&Atom> for AtomType {
    fn from(value: &Atom) -> Self {
        Self {
            residue_name: value.get_residue_name().to_owned(),
            atom_name: value.get_atom_name().to_owned(),
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Order {
    order: f32,
    n_samples: usize,
}

pub(crate) struct Map {
    params: OrderMap,
    values: GridMap<f32, f32, Box<dyn Fn(&f32) -> f32>>,
    samples: GridMap<usize, usize, Box<dyn Fn(&usize) -> usize>>,
}
