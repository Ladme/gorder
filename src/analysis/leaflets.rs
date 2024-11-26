// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use core::f32;

use super::auxiliary::{create_group, macros::group_name};
use crate::{
    errors::{AnalysisError, TopologyError},
    Leaflet, LeafletClassification, PANIC_MESSAGE,
};
use getset::{CopyGetters, Getters, MutGetters, Setters};
use groan_rs::{
    errors::{AtomError, PositionError, SimBoxError},
    prelude::{AtomIteratorWithBox, Cylinder, Dimension, Vector3D},
    structures::group::Group,
    system::System,
};

impl LeafletClassification {
    /// Prepare the system for leaflet classification.
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

pub(super) trait LeafletClassifier {
    /// Assign a molecule into the membrane leaflet.
    fn assign_to_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError>;
}

#[derive(Debug, Clone)]
pub(super) enum MoleculeLeafletClassification {
    Global(GlobalClassification),
    Local(LocalClassification),
    Individual(IndividualClassification),
}

impl MoleculeLeafletClassification {
    pub(super) fn new(params: &LeafletClassification, membrane_normal: Dimension) -> Self {
        match params {
            LeafletClassification::Global(_) => Self::Global(GlobalClassification {
                heads: Vec::new(),
                membrane_center: Vector3D::new(0.0, 0.0, 0.0),
                membrane_normal,
            }),
            LeafletClassification::Local(_) => Self::Local(LocalClassification {
                heads: Vec::new(),
                radius: params.get_radius().expect(PANIC_MESSAGE),
                membrane_center: Vec::new(),
                membrane_normal,
            }),
            LeafletClassification::Individual(_) => Self::Individual(IndividualClassification {
                heads: Vec::new(),
                methyls: Vec::new(),
                membrane_normal,
            }),
        }
    }

    /// Insert new molecule into the classifier.
    #[inline(always)]
    pub(super) fn insert(
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

impl LeafletClassifier for MoleculeLeafletClassification {
    #[inline(always)]
    fn assign_to_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        match self {
            MoleculeLeafletClassification::Global(x) => x.assign_to_leaflet(system, molecule_index),
            MoleculeLeafletClassification::Local(x) => x.assign_to_leaflet(system, molecule_index),
            MoleculeLeafletClassification::Individual(x) => {
                x.assign_to_leaflet(system, molecule_index)
            }
        }
    }
}

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
    #[inline(always)]
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(molecule, system)?);
        Ok(())
    }
}

impl LeafletClassifier for GlobalClassification {
    /// Assign a molecule into membrane leaflet.
    #[inline(always)]
    fn assign_to_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        common_assign_to_leaflet(
            &self.heads,
            molecule_index,
            system,
            self.membrane_center(),
            self.membrane_normal(),
        )
    }
}

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
                .ok_or_else(|| AnalysisError::UndefinedPosition(index + 1))?;

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
                .map_err(|_| AnalysisError::InvalidLocalMembraneCenter(index + 1))?;

            self.membrane_center[i] = center;
        }

        Ok(())
    }
}

impl LeafletClassifier for LocalClassification {
    /// Assign a molecule into membrane leaflet.
    #[inline]
    fn assign_to_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        common_assign_to_leaflet(
            &self.heads,
            molecule_index,
            system,
            self.membrane_center()
                .get(molecule_index)
                .expect(PANIC_MESSAGE),
            self.membrane_normal(),
        )
    }
}

/// Handles assignment of lipid into a membrane leaflet for the Global and Local classification.
fn common_assign_to_leaflet(
    heads: &Vec<usize>,
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
            system
                .get_box()
                .ok_or_else(|| AnalysisError::UndefinedBox)?,
        )
        .map_err(|_| AnalysisError::UndefinedPosition(head_index + 1))?;

    if distance >= 0.0 {
        Ok(Leaflet::Upper)
    } else {
        Ok(Leaflet::Lower)
    }
}

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
    /// Assign a molecule into membrane leaflet.
    fn assign_to_leaflet(
        &self,
        system: &System,
        molecule_index: usize,
    ) -> Result<Leaflet, AnalysisError> {
        let head_index = self.heads.get(molecule_index).expect(PANIC_MESSAGE);

        let mut total_distance = 0.0;
        for methyl_index in self.methyls.get(molecule_index).expect(PANIC_MESSAGE) {
            total_distance += match system.atoms_distance(*head_index, *methyl_index, self.membrane_normal) {
                Ok(x) => x,
                Err(AtomError::OutOfRange(x)) => panic!("FATAL GORDER ERROR | IndividualClassification::assign_to_leaflet | Index '{}' out of range. {}", x, PANIC_MESSAGE),
                Err(AtomError::InvalidSimBox(SimBoxError::DoesNotExist)) => return Err(AnalysisError::UndefinedBox),
                Err(AtomError::InvalidSimBox(SimBoxError::NotOrthogonal)) => return Err(AnalysisError::NotOrthogonalBox),
                Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => return Err(AnalysisError::UndefinedPosition(x)),
                Err(e) => panic!("FATAL GORDER ERROR | IndividualClassification::assign_to_leaflet | Unexpected error type '{}' returned. {}", e, PANIC_MESSAGE),
            };
        }

        if total_distance >= 0.0 {
            Ok(Leaflet::Upper)
        } else {
            Ok(Leaflet::Lower)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::auxiliary::macros::group_name;
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
            Dimension::Z,
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
            MoleculeLeafletClassification::Global(x) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.set_membrane_center(membrane_center);
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.assign_to_leaflet(&system, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.assign_to_leaflet(&system, 1).unwrap(),
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
            MoleculeLeafletClassification::Local(x) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.membrane_center_mut().push(Vector3D::new(0.0, 0.0, 0.0));
                x.membrane_center_mut().push(Vector3D::new(0.0, 0.0, 0.0));
                x.set_membrane_center(&system, Dimension::Z).unwrap();
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.assign_to_leaflet(&system, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.assign_to_leaflet(&system, 1).unwrap(),
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
            MoleculeLeafletClassification::Individual(x) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.methyls_mut().push(vec![1453, 1496]);
                x.methyls_mut().push(vec![11953, 11996]);
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.assign_to_leaflet(&system, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.assign_to_leaflet(&system, 1).unwrap(),
            Leaflet::Lower
        );
    }
}
