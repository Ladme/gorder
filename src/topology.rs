// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::ops::{Add, AddAssign};

use hashbrown::HashSet;

use getset::Getters;
use groan_rs::{
    prelude::{Atom, GridMap, SimBox, Vector3D},
    structures::{gridmap::DataOrder, group::Group},
    system::System,
};

use crate::{
    auxiliary::{macros::group_name, PANIC_MESSAGE},
    errors::{AnalysisError, TopologyError},
    leaflets::LeafletClassification,
    ordermap::OrderMap,
};

#[derive(Debug, Clone, Getters)]
pub(crate) struct MoleculeType {
    #[getset(get = "pub(crate)")]
    name: String,
    #[getset(get = "pub(crate)")]
    topology: MoleculeTopology,
    #[getset(get = "pub(crate)")]
    order_bonds: OrderBonds,
    #[getset(get = "pub(crate)")]
    order_atoms: OrderAtoms,
    #[getset(get = "pub(crate)")]
    leaflet_classification: Option<MoleculeLeafletClassification>,
}

impl MoleculeType {
    pub(crate) fn new(
        system: &System,
        name: &str,
        topology: &MoleculeTopology,
        order_bonds: &HashSet<(usize, usize)>,
        order_atoms: &[usize],
        min_index: usize,
        leaflet_classification: Option<MoleculeLeafletClassification>,
    ) -> Self {
        Self {
            name: name.to_owned(),
            topology: topology.to_owned(),
            order_bonds: OrderBonds::new(system, order_bonds, min_index),
            order_atoms: OrderAtoms::new(system, order_atoms, min_index),
            leaflet_classification,
        }
    }

    /// Insert new bond instances to the molecule.
    pub(crate) fn insert(
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
    pub(crate) fn analyze_frame(
        &mut self,
        frame: &System,
        simbox: &SimBox,
        membrane_normal: &Vector3D,
    ) -> Result<(), AnalysisError> {
        // todo! leaflet classification

        for bond_type in self.order_bonds.bonds.iter_mut() {
            for (index1, index2) in bond_type.bonds.iter() {
                let atom1 = unsafe { frame.get_atom_unchecked_as_ref(*index1) };
                let atom2 = unsafe { frame.get_atom_unchecked_as_ref(*index2) };

                let pos1 = atom1
                    .get_position()
                    .ok_or_else(|| AnalysisError::UndefinedPosition(index1 + 1))?;

                let pos2 = atom2
                    .get_position()
                    .ok_or_else(|| AnalysisError::UndefinedPosition(index2 + 1))?;

                let vector = pos1.vector_to(pos2, simbox);
                let angle = vector.angle(membrane_normal);

                let cos = angle.cos();
                let sch = 0.5 * (1.0 - (3.0 * cos * cos));

                bond_type.total += sch;
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

#[derive(Debug, Clone)]
pub(crate) enum MoleculeLeafletClassification {
    Global(GlobalClassification),
    Local(LocalClassification),
    Individual(IndividualClassification),
}

impl MoleculeLeafletClassification {
    pub(crate) fn new(params: &LeafletClassification) -> Self {
        match params {
            LeafletClassification::Global(_) => {
                Self::Global(GlobalClassification { heads: Vec::new() })
            }
            LeafletClassification::Local(_) => Self::Local(LocalClassification {
                heads: Vec::new(),
                radius: params.get_radius().expect(PANIC_MESSAGE),
            }),
            LeafletClassification::Individual(_) => Self::Individual(IndividualClassification {
                heads: Vec::new(),
                methyls: Vec::new(),
            }),
        }
    }

    /// Insert new molecule into the classifier.
    pub(crate) fn insert(
        &mut self,
        molecule: &Group,
        system: &System,
    ) -> Result<(), TopologyError> {
        match self {
            Self::Global(x) => x.insert(molecule, system),
            Self::Local(x) => x.insert(molecule, system),
            Self::Individual(x) => x.insert(molecule, system),
        }
    }
}

fn get_reference_head(molecule: &Group, system: &System) -> Result<usize, TopologyError> {
    let group_name = group_name!("Heads");
    let mut atoms = Vec::new();
    for index in molecule.get_atoms().iter() {
        if system.group_isin(group_name, index).expect(PANIC_MESSAGE) {
            atoms.push(index);
        }
    }

    if atoms.len() == 0 {
        return Err(TopologyError::NoHead(
                molecule
                    .get_atoms()
                    .first()
                    .unwrap_or_else(|| panic!("FATAL GORDER ERROR | topology::get_reference_head | No atoms detected inside a molecule. {}", PANIC_MESSAGE))));
    }

    if atoms.len() > 1 {
        return Err(TopologyError::MultipleHeads(
            molecule.get_atoms().first().expect(PANIC_MESSAGE),
        ));
    }

    Ok(*atoms.get(0).expect(PANIC_MESSAGE))
}

fn get_reference_methyls(molecule: &Group, system: &System) -> Result<Vec<usize>, TopologyError> {
    let group_name = group_name!("Methyls");
    let mut atoms = Vec::new();

    for index in molecule.get_atoms().iter() {
        if system.group_isin(group_name, index).expect(PANIC_MESSAGE) {
            atoms.push(index);
        }
    }

    if atoms.len() == 0 {
        return Err(TopologyError::NoMethyl(
            molecule
                .get_atoms()
                .first()
                .unwrap_or_else(|| panic!("FATAL GORDER ERROR | topology::get_reference_methyls | No atoms detected inside a molecule. {}", PANIC_MESSAGE))));
    }

    Ok(atoms)
}

#[derive(Debug, Clone)]
pub(crate) struct GlobalClassification {
    /// Indices of headgroup identifiers (one per molecule).
    pub(crate) heads: Vec<usize>,
}

impl GlobalClassification {
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(molecule, system)?);
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub(crate) struct LocalClassification {
    /// Indices of headgroup identifiers (one per molecule).
    pub(crate) heads: Vec<usize>,
    /// Radius of a cylinder for the calculation of local membrane center of geometry (in nm).
    pub(crate) radius: f32,
}

impl LocalClassification {
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(molecule, system)?);
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub(crate) struct IndividualClassification {
    /// Indices of headgroup identifiers (one per molecule).
    pub(crate) heads: Vec<usize>,
    /// Indices of methyl identifiers (any number per molecule).
    pub(crate) methyls: Vec<Vec<usize>>,
}

impl IndividualClassification {
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(molecule, system)?);
        self.methyls.push(get_reference_methyls(molecule, system)?);

        // check that the number of methyls is concistent in the molecule
        if self.methyls.len() > 1 {
            let curr_methyls = self.methyls[self.methyls.len() - 1].len();
            let first_methyls = self.methyls[0].len();
            if curr_methyls != first_methyls {
                return Err(TopologyError::InconsistentNumberOfMethyls(
                    molecule.get_atoms().first().expect(PANIC_MESSAGE),
                    curr_methyls,
                    first_methyls,
                ));
            }
        }

        Ok(())
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
#[derive(Debug, Clone)]
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

    /// Insert new real bonds to already constructed order bonds.
    pub(crate) fn insert(
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
pub(crate) struct OrderAtoms {
    /// Ordered by the increasing relative index.
    pub(crate) atoms: Vec<AtomType>,
}

impl OrderAtoms {
    pub(crate) fn new(system: &System, atoms: &[usize], minimum_index: usize) -> OrderAtoms {
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

#[derive(Debug, Clone)]
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

    /// Insert new real bond to the current order bond.
    pub(crate) fn insert(&mut self, abs_index_1: usize, abs_index_2: usize) {
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

#[derive(Debug, Clone)]
pub(crate) struct Map {
    params: OrderMap,
    values: GridMap<f32, f32, fn(&f32) -> f32>,
    samples: GridMap<usize, usize, fn(&usize) -> usize>,
}

impl Add<Map> for Map {
    type Output = Map;

    fn add(self, rhs: Map) -> Self::Output {
        let joined_values_map =
            merge_grid_maps(self.values, rhs.values, f32::clone as fn(&f32) -> f32);

        let joined_samples_map = merge_grid_maps(
            self.samples,
            rhs.samples,
            usize::clone as fn(&usize) -> usize,
        );

        Map {
            params: self.params,
            values: joined_values_map,
            samples: joined_samples_map,
        }
    }
}

/// Helper function for merging optional Maps.
fn merge_option_maps(lhs: Option<Map>, rhs: Option<Map>) -> Option<Map> {
    match (lhs, rhs) {
        (Some(x), Some(y)) => Some(x + y),
        (None, None) => None,
        (Some(_), None) | (None, Some(_)) => panic!(
            "FATAL GORDER ERROR | merge_option_maps | Inconsistent option value. {}",
            PANIC_MESSAGE
        ),
    }
}

/// Helper function to merge two GridMaps and construct a new GridMap.
fn merge_grid_maps<T, U>(
    lhs: GridMap<T, U, fn(&T) -> U>,
    rhs: GridMap<T, U, fn(&T) -> U>,
    clone_fn: fn(&T) -> U,
) -> GridMap<T, U, fn(&T) -> U>
where
    T: Add<Output = T> + Copy + Default + std::fmt::Debug,
    U: std::fmt::Display,
{
    let merged = lhs
        .extract_raw()
        .zip(rhs.extract_raw())
        .map(|((_, _, &x), (_, _, &y))| x + y)
        .collect::<Vec<T>>();

    GridMap::from_vec(
        lhs.span_x(),
        lhs.span_y(),
        lhs.tile_dim(),
        merged,
        DataOrder::default(),
        clone_fn,
    )
    .unwrap_or_else(|_| {
        panic!(
            "FATAL GORDER ERROR | merge_grid_maps | Could not construct merged map. {}",
            PANIC_MESSAGE
        )
    })
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::auxiliary::create_group;

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
        let order_bonds = OrderBonds::new(&system, &bonds_set, 125);

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
        let mut order_bonds = OrderBonds::new(&system, &bonds_set, 125);

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

    #[test]
    fn merge_order() {
        let order1 = Order {
            order: 0.78,
            n_samples: 31,
        };

        let order2 = Order {
            order: 0.43,
            n_samples: 14,
        };

        let merged = order1 + order2;

        assert_relative_eq!(merged.order, 1.21);
        assert_eq!(merged.n_samples, 45);
    }

    #[test]
    fn merge_map() {
        let values = vec![1.0, 2.5, 3.0, 4.2, 5.3, 6.1, 7.3, 8.9];
        let map1_values = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::RowMajor,
            f32::to_owned as fn(&f32) -> f32,
        )
        .unwrap();

        let samples = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let map1_samples = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            samples,
            DataOrder::RowMajor,
            usize::to_owned as fn(&usize) -> usize,
        )
        .unwrap();

        let map1 = Map {
            params: OrderMap::new().output_directory(".").build().unwrap(),
            values: map1_values,
            samples: map1_samples,
        };

        let values = vec![0.7, 1.4, 2.1, 1.4, 2.3, 3.1, 3.3, 1.9];
        let map2_values = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::RowMajor,
            f32::to_owned as fn(&f32) -> f32,
        )
        .unwrap();

        let samples = vec![2, 1, 4, 3, 2, 1, 0, 2];
        let map2_samples = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            samples,
            DataOrder::RowMajor,
            usize::to_owned as fn(&usize) -> usize,
        )
        .unwrap();

        let map2 = Map {
            params: OrderMap::new().output_directory(".").build().unwrap(),
            values: map2_values,
            samples: map2_samples,
        };

        let map = map1 + map2;

        let expected_values = vec![1.7, 3.9, 5.1, 5.6, 7.6, 9.2, 10.6, 10.8];
        let expected_samples = vec![3, 3, 7, 7, 7, 7, 7, 10];

        let values = map
            .values
            .extract_raw()
            .map(|(_, _, &x)| x)
            .collect::<Vec<f32>>();
        let samples = map
            .samples
            .extract_raw()
            .map(|(_, _, &x)| x)
            .collect::<Vec<usize>>();

        assert_eq!(values.len(), expected_values.len());
        assert_eq!(samples.len(), expected_samples.len());

        for (got, exp) in values.iter().zip(expected_values.iter()) {
            assert_relative_eq!(got, exp);
        }

        for (got, exp) in samples.iter().zip(expected_samples.iter()) {
            assert_eq!(got, exp);
        }
    }

    #[test]
    fn test_global_leaflet_classification() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(&LeafletClassification::global(
            "@membrane",
            "name P",
        ));

        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();
        let group3 = Group::from_query("resid 264", &system).unwrap();

        classifier.insert(&group1, &system).unwrap();
        classifier.insert(&group2, &system).unwrap();
        classifier.insert(&group3, &system).unwrap();

        if let MoleculeLeafletClassification::Global(x) = classifier {
            assert_eq!(x.heads, vec![760, 18002, 34047]);
        } else {
            panic!("Invalid classifier type.")
        }
    }

    #[test]
    fn test_local_leaflet_classification() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(&LeafletClassification::local(
            "@membrane",
            "name P",
            3.3,
        ));

        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();
        let group3 = Group::from_query("resid 264", &system).unwrap();

        classifier.insert(&group1, &system).unwrap();
        classifier.insert(&group2, &system).unwrap();
        classifier.insert(&group3, &system).unwrap();

        if let MoleculeLeafletClassification::Local(x) = classifier {
            assert_eq!(x.heads, vec![760, 18002, 34047]);
            assert_eq!(x.radius, 3.3);
        } else {
            panic!("Invalid classifier type.")
        }
    }

    #[test]
    fn test_individual_leaflet_classification() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        create_group(&mut system, "Methyls", "name C218 C316").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316"),
        );

        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();
        let group3 = Group::from_query("resid 264", &system).unwrap();

        classifier.insert(&group1, &system).unwrap();
        classifier.insert(&group2, &system).unwrap();
        classifier.insert(&group3, &system).unwrap();

        if let MoleculeLeafletClassification::Individual(x) = classifier {
            assert_eq!(x.heads, vec![760, 18002, 34047]);
            assert_eq!(x.methyls, vec![[828, 871], [18070, 18113], [34115, 34158]]);
        } else {
            panic!("Invalid classifier type.")
        }
    }

    /// Helper function to run the classification test and check for the expected error.
    fn run_classification_test<F>(
        system: &mut System,
        classifier: LeafletClassification,
        group_query: &str,
        is_expected_error: F,
    ) where
        F: Fn(&TopologyError) -> bool,
    {
        let mut molecule_classifier = MoleculeLeafletClassification::new(&classifier);
        let group = Group::from_query(group_query, system).unwrap();
        match molecule_classifier.insert(&group, system) {
            Err(e) if is_expected_error(&e) => (),
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(e) => panic!("Incorrect error type returned: {:?}", e),
        }
    }

    #[test]
    fn test_global_leaflet_classification_fail_no_head() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P and not resid 144").unwrap();
        let classifier = LeafletClassification::global("@membrane", "name P and not resid 144");

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::NoHead(_))
        });
    }

    #[test]
    fn test_global_leaflet_classification_fail_multiple_heads() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P or (resid 144 and name P HA)").unwrap();
        let classifier =
            LeafletClassification::global("@membrane", "name P or (resid 144 and name P HA)");

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::MultipleHeads(_))
        });
    }

    #[test]
    fn test_local_leaflet_classification_fail_no_head() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P and not resid 144").unwrap();
        let classifier = LeafletClassification::local("@membrane", "name P and not resid 144", 2.5);

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::NoHead(_))
        });
    }

    #[test]
    fn test_local_leaflet_classification_fail_multiple_heads() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P or (resid 144 and name P HA)").unwrap();
        let classifier =
            LeafletClassification::local("@membrane", "name P or (resid 144 and name P HA)", 2.5);

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::MultipleHeads(_))
        });
    }

    #[test]
    fn test_individual_leaflet_classification_fail_no_head() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P and not resid 144").unwrap();
        create_group(&mut system, "Methyls", "name C218 C316").unwrap();
        let classifier =
            LeafletClassification::individual("name P and not resid 144", "name C218 C316");

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::NoHead(_))
        });
    }

    #[test]
    fn test_individual_leaflet_classification_fail_multiple_heads() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P or (resid 144 and name P HA)").unwrap();
        create_group(&mut system, "Methyls", "name C218 C316").unwrap();
        let classifier = LeafletClassification::individual(
            "name P or (resid 144 and name P HA)",
            "name C218 C316",
        );

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::MultipleHeads(_))
        });
    }

    #[test]
    fn test_individual_leaflet_classification_fail_no_methyl() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        create_group(&mut system, "Methyls", "name C218 C316 and not resid 144").unwrap();
        let classifier =
            LeafletClassification::individual("name P", "name C218 C316 and not resid 144");

        run_classification_test(&mut system, classifier, "resid 144", |e| {
            matches!(e, TopologyError::NoMethyl(_))
        });
    }

    #[test]
    fn test_individual_leaflet_classification_inconsistent_methyls() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        create_group(
            &mut system,
            "Methyls",
            "(name C218 C316 and not resid 144) or (resid 144 and name C218)",
        )
        .unwrap();
        let classifier = LeafletClassification::individual(
            "name P",
            "(name C218 C316 and not resid 144) or (resid 144 and name C218)",
        );

        let mut molecule_classifier = MoleculeLeafletClassification::new(&classifier);
        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();

        molecule_classifier.insert(&group1, &system).unwrap();

        match molecule_classifier.insert(&group2, &system) {
            Err(TopologyError::InconsistentNumberOfMethyls(_, a, b)) => {
                assert_eq!(a, 1);
                assert_eq!(b, 2);
            }
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(e) => panic!("Incorrect error type returned: {:?}", e),
        }
    }
}
