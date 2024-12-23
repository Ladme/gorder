// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use core::f32;

use super::common::{create_group, macros::group_name};
use crate::{
    errors::{AnalysisError, TopologyError},
    input::LeafletClassification,
    Leaflet, PANIC_MESSAGE,
};
use getset::{CopyGetters, Getters, MutGetters, Setters};
use groan_rs::{
    errors::{AtomError, PositionError, SimBoxError},
    prelude::{AtomIteratorWithBox, Cylinder, Dimension, Vector3D},
    structures::group::Group,
    system::System,
};

impl LeafletClassification {
    /// Create groups in the system that are required for leaflet classification.
    pub(super) fn prepare_system(&self, system: &mut System) -> Result<(), TopologyError> {
        match self {
            Self::Global(params) => {
                create_group(system, "Membrane", params.membrane())?;
                create_group(system, "Heads", params.heads())?;
            }
            Self::Local(params) => {
                create_group(system, "Membrane", params.membrane())?;
                create_group(system, "Heads", params.heads())?;
            }
            Self::Individual(params) => {
                create_group(system, "Heads", params.heads())?;
                create_group(system, "Methyls", params.methyls())?;
            }
        }

        Ok(())
    }
}

/// Type of leaflet classification method used for a molecule.
#[derive(Debug, Clone)]
pub(super) enum MoleculeLeafletClassification {
    Global(GlobalClassification, AssignedLeaflets),
    Local(LocalClassification, AssignedLeaflets),
    Individual(IndividualClassification, AssignedLeaflets),
}

impl MoleculeLeafletClassification {
    /// Convert the input `LeafletClassification` into an enum that is used in the analysis.
    pub(super) fn new(params: &LeafletClassification, membrane_normal: Dimension) -> Self {
        match params {
            LeafletClassification::Global(_) => Self::Global(
                GlobalClassification {
                    heads: Vec::new(),
                    membrane_center: Vector3D::new(0.0, 0.0, 0.0),
                    membrane_normal,
                },
                AssignedLeaflets::default(),
            ),
            LeafletClassification::Local(_) => Self::Local(
                LocalClassification {
                    heads: Vec::new(),
                    radius: params.get_radius().expect(PANIC_MESSAGE),
                    membrane_center: Vec::new(),
                    membrane_normal,
                },
                AssignedLeaflets::default(),
            ),
            LeafletClassification::Individual(_) => Self::Individual(
                IndividualClassification {
                    heads: Vec::new(),
                    methyls: Vec::new(),
                    membrane_normal,
                },
                AssignedLeaflets::default(),
            ),
        }
    }

    /// Insert new molecule into the leaflet classifier.
    #[inline(always)]
    pub(super) fn insert(
        &mut self,
        molecule: &Group,
        system: &System,
    ) -> Result<(), TopologyError> {
        match self {
            Self::Global(x, y) => {
                x.insert(molecule, system)?;
                y.assignment.push(None);
            }
            Self::Local(x, y) => {
                x.insert(molecule, system)?;
                y.assignment.push(None);
            }
            Self::Individual(x, y) => {
                x.insert(molecule, system)?;
                y.assignment.push(None);
            }
        }

        Ok(())
    }
}

/// Get index of an atom that represents the head of the given lipid molecule.
fn get_reference_head(molecule: &Group, system: &System) -> Result<usize, TopologyError> {
    let group_name = group_name!("Heads");
    let mut atoms = Vec::new();
    for index in molecule.get_atoms().iter() {
        if system.group_isin(group_name, index).expect(PANIC_MESSAGE) {
            atoms.push(index);
        }
    }

    if atoms.is_empty() {
        return Err(TopologyError::NoHead(
            molecule
                .get_atoms()
                .first()
                .unwrap_or_else(|| panic!("FATAL GORDER ERROR | leaflets::get_reference_head | No atoms detected inside a molecule. {}", PANIC_MESSAGE))));
    }

    if atoms.len() > 1 {
        return Err(TopologyError::MultipleHeads(
            molecule.get_atoms().first().expect(PANIC_MESSAGE),
        ));
    }

    Ok(*atoms.first().expect(PANIC_MESSAGE))
}

/// Get indices of atoms representing the methyls (or ends of tails) of the given lipid molecule.
fn get_reference_methyls(molecule: &Group, system: &System) -> Result<Vec<usize>, TopologyError> {
    let group_name = group_name!("Methyls");
    let mut atoms = Vec::new();

    for index in molecule.get_atoms().iter() {
        if system.group_isin(group_name, index).expect(PANIC_MESSAGE) {
            atoms.push(index);
        }
    }

    if atoms.is_empty() {
        return Err(TopologyError::NoMethyl(
            molecule
                .get_atoms()
                .first()
                .unwrap_or_else(|| panic!("FATAL GORDER ERROR | leaflets::get_reference_methyls | No atoms detected inside a molecule. {}", PANIC_MESSAGE))));
    }

    Ok(atoms)
}

impl MoleculeLeafletClassification {
    /// Identify leaflet in which the molecule with the specified index is located.
    #[inline(always)]
    #[allow(unused)]
    fn identify_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        match self {
            MoleculeLeafletClassification::Global(x, _) => {
                x.identify_leaflet(system, molecule_index)
            }
            MoleculeLeafletClassification::Local(x, _) => {
                x.identify_leaflet(system, molecule_index)
            }
            MoleculeLeafletClassification::Individual(x, _) => {
                x.identify_leaflet(system, molecule_index)
            }
        }
    }

    /// Assign all lipids into their respective leaflets.
    #[inline(always)]
    pub(super) fn assign_lipids(&mut self, system: &System) -> Result<(), AnalysisError> {
        match self {
            MoleculeLeafletClassification::Global(x, y) => y.assign_lipids(system, x),
            MoleculeLeafletClassification::Local(x, y) => y.assign_lipids(system, x),
            MoleculeLeafletClassification::Individual(x, y) => y.assign_lipids(system, x),
        }
    }

    /// Get the leaflet a molecule with target index is located in.
    /// This performs no calculation, this only returns the already calculated leaflet assignment.
    #[inline(always)]
    pub(super) fn get_assigned_leaflet(&self, molecule_index: usize) -> Option<Leaflet> {
        match self {
            MoleculeLeafletClassification::Global(_, y)
            | MoleculeLeafletClassification::Local(_, y)
            | MoleculeLeafletClassification::Individual(_, y) => {
                y.get_assigned_leaflet(molecule_index)
            }
        }
    }
}

/// Trait implemented by all leaflet classification methods.
trait LeafletClassifier {
    /// Caclulate membrane leaflet the specified molecule belongs to.
    fn identify_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError>;

    /// Get the number of molecules in the system.
    fn n_molecules(&self) -> usize;
}

/// Leaflet classification method that uses positions of lipid heads and the global membrane center of geometry.
#[derive(Debug, Clone, CopyGetters, Getters, MutGetters, Setters)]
pub(super) struct GlobalClassification {
    /// Indices of headgroup identifiers (one per molecule).
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    heads: Vec<usize>,
    /// Global membrane center of geometry.
    #[getset(get = "pub(super)", set = "pub(super)")]
    membrane_center: Vector3D,
    /// Orientation of the membrane normal.
    #[getset(get_copy = "pub(super)")]
    membrane_normal: Dimension,
}

impl GlobalClassification {
    /// Insert a new molecule into the classifier.
    #[inline(always)]
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(molecule, system)?);
        Ok(())
    }
}

impl LeafletClassifier for GlobalClassification {
    #[inline(always)]
    fn identify_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        common_identify_leaflet(
            &self.heads,
            molecule_index,
            system,
            self.membrane_center(),
            self.membrane_normal(),
        )
    }

    #[inline(always)]
    fn n_molecules(&self) -> usize {
        self.heads.len()
    }
}

/// Leaflet classification method that uses positions of lipid heads and the local membrane center of geometry.
#[derive(Debug, Clone, Getters, CopyGetters, MutGetters)]
pub(super) struct LocalClassification {
    /// Indices of headgroup identifiers (one per molecule).
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    heads: Vec<usize>,
    /// Radius of a cylinder for the calculation of local membrane center of geometry (in nm).
    #[getset(get_copy = "pub(super)", get_mut = "pub(super)")]
    radius: f32,
    /// Local membrane center of geometry for each molecule.
    #[getset(get = "pub(super)", get_mut)]
    membrane_center: Vec<Vector3D>,
    /// Orientation of the membrane normal.
    #[getset(get_copy = "pub(super)")]
    membrane_normal: Dimension,
}

impl LocalClassification {
    #[inline(always)]
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(molecule, system)?);
        self.membrane_center.push(Vector3D::new(0.0, 0.0, 0.0));
        Ok(())
    }

    /// Calculate and set local membrane center of geometry for each molecule.
    pub(super) fn set_membrane_center(
        &mut self,
        system: &System,
        membrane_normal: Dimension,
    ) -> Result<(), AnalysisError> {
        for (i, &index) in self.heads.iter().enumerate() {
            let position = unsafe { system.get_atom_unchecked(index) }
                .get_position()
                .ok_or(AnalysisError::UndefinedPosition(index))?;

            let cylinder = Cylinder::new(
                position.clone(),
                self.radius,
                f32::INFINITY,
                membrane_normal,
            );
            let center = system
                .group_iter(group_name!("Membrane"))
                .expect(PANIC_MESSAGE)
                .filter_geometry(cylinder)
                .get_center()
                .map_err(|_| AnalysisError::InvalidLocalMembraneCenter(index))?;

            self.membrane_center[i] = center;
        }

        Ok(())
    }
}

impl LeafletClassifier for LocalClassification {
    #[inline(always)]
    fn identify_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        common_identify_leaflet(
            &self.heads,
            molecule_index,
            system,
            self.membrane_center()
                .get(molecule_index)
                .expect(PANIC_MESSAGE),
            self.membrane_normal(),
        )
    }

    #[inline(always)]
    fn n_molecules(&self) -> usize {
        self.heads.len()
    }
}

/// Handles classification of lipid into a membrane leaflet for the Global and Local classification.
fn common_identify_leaflet(
    heads: &[usize],
    molecule_index: usize,
    system: &System,
    membrane_center: &Vector3D,
    membrane_normal: Dimension,
) -> Result<Leaflet, AnalysisError> {
    let head_index = *heads.get(molecule_index).expect(PANIC_MESSAGE);
    let head = unsafe { system.get_atom_unchecked(head_index) };

    let distance = head
        .distance_from_point(
            membrane_center,
            membrane_normal,
            system.get_box().ok_or(AnalysisError::UndefinedBox)?,
        )
        .map_err(|_| AnalysisError::UndefinedPosition(head_index))?;

    if distance >= 0.0 {
        Ok(Leaflet::Upper)
    } else {
        Ok(Leaflet::Lower)
    }
}

/// Leaflet classification method that uses positions of lipid heads and tail ends.
#[derive(Debug, Clone, Getters, MutGetters)]
pub(crate) struct IndividualClassification {
    /// Indices of headgroup identifiers (one per molecule).
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    heads: Vec<usize>,
    /// Indices of methyl identifiers (any number per molecule).
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    methyls: Vec<Vec<usize>>,
    /// Orientation of the membrane normal.
    #[getset(get_copy = "pub(super)")]
    membrane_normal: Dimension,
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

impl LeafletClassifier for IndividualClassification {
    fn identify_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        let head_index = self.heads.get(molecule_index).expect(PANIC_MESSAGE);

        let mut total_distance = 0.0;
        for methyl_index in self.methyls.get(molecule_index).expect(PANIC_MESSAGE) {
            total_distance += match system.atoms_distance(*head_index, *methyl_index, self.membrane_normal) {
                Ok(x) => x,
                Err(AtomError::OutOfRange(x)) => panic!("FATAL GORDER ERROR | IndividualClassification::identify_leaflet | Index '{}' out of range. {}", x, PANIC_MESSAGE),
                Err(AtomError::InvalidSimBox(SimBoxError::DoesNotExist)) => return Err(AnalysisError::UndefinedBox),
                Err(AtomError::InvalidSimBox(SimBoxError::NotOrthogonal)) => return Err(AnalysisError::NotOrthogonalBox),
                Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => return Err(AnalysisError::UndefinedPosition(x)),
                Err(e) => panic!("FATAL GORDER ERROR | IndividualClassification::identify_leaflet | Unexpected error type '{}' returned. {}", e, PANIC_MESSAGE),
            };
        }

        if total_distance >= 0.0 {
            Ok(Leaflet::Upper)
        } else {
            Ok(Leaflet::Lower)
        }
    }

    #[inline(always)]
    fn n_molecules(&self) -> usize {
        self.heads.len()
    }
}

#[derive(Debug, Clone, Default)]
pub(super) struct AssignedLeaflets {
    assignment: Vec<Option<Leaflet>>,
}

impl AssignedLeaflets {
    /// Assign all lipids into membrane leaflets.
    ///
    /// ## Panic
    /// Panics if the number of molecules in the classifier does not match the number of molecules in the assignment vector.
    fn assign_lipids(
        &mut self,
        system: &System,
        classifier: &impl LeafletClassifier,
    ) -> Result<(), AnalysisError> {
        assert_eq!(classifier.n_molecules(), self.assignment.len(), "FATAL GORDER ERROR | AssignedLeaflets::assign_lipids | Inconsistent number of molecules.");

        self.assignment
            .iter_mut()
            .enumerate()
            .try_for_each(|(index, assignment)| {
                *assignment = Some(classifier.identify_leaflet(system, index)?);
                Ok(())
            })?;

        Ok(())
    }

    /// Get the leaflet that was assigned to molecule of target index.
    fn get_assigned_leaflet(&self, molecule_index: usize) -> Option<Leaflet> {
        *self.assignment.get(molecule_index).expect(PANIC_MESSAGE)
    }
}

#[cfg(test)]
mod tests {
    use super::super::common::macros::group_name;
    use super::*;

    #[test]
    fn test_global_leaflet_classification() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            Dimension::Z,
        );

        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();
        let group3 = Group::from_query("resid 264", &system).unwrap();

        classifier.insert(&group1, &system).unwrap();
        classifier.insert(&group2, &system).unwrap();
        classifier.insert(&group3, &system).unwrap();

        if let MoleculeLeafletClassification::Global(x, _) = classifier {
            assert_eq!(x.heads, vec![760, 18002, 34047]);
        } else {
            panic!("Invalid classifier type.")
        }
    }

    #[test]
    fn test_local_leaflet_classification() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 3.3),
            Dimension::Z,
        );

        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();
        let group3 = Group::from_query("resid 264", &system).unwrap();

        classifier.insert(&group1, &system).unwrap();
        classifier.insert(&group2, &system).unwrap();
        classifier.insert(&group3, &system).unwrap();

        if let MoleculeLeafletClassification::Local(x, _) = classifier {
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
            Dimension::Z,
        );

        let group1 = Group::from_query("resid 7", &system).unwrap();
        let group2 = Group::from_query("resid 144", &system).unwrap();
        let group3 = Group::from_query("resid 264", &system).unwrap();

        classifier.insert(&group1, &system).unwrap();
        classifier.insert(&group2, &system).unwrap();
        classifier.insert(&group3, &system).unwrap();

        if let MoleculeLeafletClassification::Individual(x, _) = classifier {
            assert_eq!(x.heads, vec![760, 18002, 34047]);
            assert_eq!(x.methyls, vec![[828, 871], [18070, 18113], [34115, 34158]]);
        } else {
            panic!("Invalid classifier type.")
        }
    }

    #[test]
    fn test_prepare_system_global() {
        let mut system = System::from_file("tests/files/pepg_cg.tpr").unwrap();
        let classifier = LeafletClassification::global("@membrane", "name PO4");

        classifier.prepare_system(&mut system).unwrap();

        assert!(system.group_exists(group_name!("Membrane")));
        assert!(system.group_exists(group_name!("Heads")));
    }

    #[test]
    fn test_prepare_system_local() {
        let mut system = System::from_file("tests/files/pepg_cg.tpr").unwrap();
        let classifier = LeafletClassification::local("@membrane", "name PO4", 2.5);

        classifier.prepare_system(&mut system).unwrap();

        assert!(system.group_exists(group_name!("Membrane")));
        assert!(system.group_exists(group_name!("Heads")));
    }

    #[test]
    fn test_prepare_system_individual() {
        let mut system = System::from_file("tests/files/pepg_cg.tpr").unwrap();
        let classifier = LeafletClassification::individual("name PO4", "name C4A C4B");

        classifier.prepare_system(&mut system).unwrap();

        assert!(system.group_exists(group_name!("Heads")));
        assert!(system.group_exists(group_name!("Methyls")));
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
        let mut molecule_classifier = MoleculeLeafletClassification::new(&classifier, Dimension::Z);
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

        let mut molecule_classifier = MoleculeLeafletClassification::new(&classifier, Dimension::Z);
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

    #[test]
    fn test_global_assign_to_leaflet() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            Dimension::Z,
        );

        let membrane_center = system
            .selection_iter("@membrane")
            .unwrap()
            .get_center()
            .unwrap();

        match &mut classifier {
            MoleculeLeafletClassification::Global(x, _) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.set_membrane_center(membrane_center);
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.identify_leaflet(&system, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.identify_leaflet(&system, 1).unwrap(),
            Leaflet::Lower
        );
    }

    #[test]
    fn test_local_assign_to_leaflet() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5),
            Dimension::Z,
        );

        system
            .group_create(group_name!("Membrane"), "@membrane")
            .unwrap();

        match &mut classifier {
            MoleculeLeafletClassification::Local(x, _) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.membrane_center_mut().push(Vector3D::new(0.0, 0.0, 0.0));
                x.membrane_center_mut().push(Vector3D::new(0.0, 0.0, 0.0));
                x.set_membrane_center(&system, Dimension::Z).unwrap();
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.identify_leaflet(&system, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.identify_leaflet(&system, 1).unwrap(),
            Leaflet::Lower
        );
    }

    #[test]
    fn test_individual_assign_to_leaflet() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316"),
            Dimension::Z,
        );

        match &mut classifier {
            MoleculeLeafletClassification::Individual(x, _) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.methyls_mut().push(vec![1453, 1496]);
                x.methyls_mut().push(vec![11953, 11996]);
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.identify_leaflet(&system, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.identify_leaflet(&system, 1).unwrap(),
            Leaflet::Lower
        );
    }
}
