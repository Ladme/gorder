// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use core::f32;
use std::{
    num::NonZeroUsize,
    sync::Arc,
    time::{Duration, Instant},
};

use super::common::{create_group, macros::group_name};
use crate::{
    errors::{AnalysisError, TopologyError},
    input::{Frequency, LeafletClassification},
    Leaflet, PANIC_MESSAGE,
};
use getset::{CopyGetters, Getters, MutGetters, Setters};
use groan_rs::{
    errors::{AtomError, PositionError, SimBoxError},
    prelude::{AtomIteratorWithBox, Cylinder, Dimension, Vector3D},
    structures::group::Group,
    system::System,
};
use hashbrown::HashMap;
use parking_lot::Mutex;

use once_cell::{sync::Lazy, unsync::OnceCell};

/// [`TIMEOUT`] in seconds.
static TIMEOUT_SECONDS: u64 = 5;
/// Global soft timeout duration for spin-lock used when fetching data for leaflet assignment.
/// After this time a warning is logged.
static TIMEOUT: Lazy<Duration> = Lazy::new(|| Duration::from_secs(TIMEOUT_SECONDS));

/// [`HARD_TIMEOUT`] in seconds.
static HARD_TIMEOUT_SECONDS: u64 = 125;
/// Global HARD timeout duration for spin-lock used when fetching data for leaflet assignment.
/// After this time a PANIC is raised.
static HARD_TIMEOUT: Lazy<Duration> = Lazy::new(|| Duration::from_secs(HARD_TIMEOUT_SECONDS));

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

/// Type of leaflet classification method used for a molecule.
#[derive(Debug, Clone)]
pub(super) enum MoleculeLeafletClassification {
    Global(GlobalClassification, AssignedLeaflets),
    Local(LocalClassification, AssignedLeaflets),
    Individual(IndividualClassification, AssignedLeaflets),
}

impl MoleculeLeafletClassification {
    /// Convert the input `LeafletClassification` into an enum that is used in the analysis.
    pub(super) fn new(
        params: &LeafletClassification,
        membrane_normal: Dimension,
        n_threads: usize,
        step_size: usize,
    ) -> Self {
        let needs_shared_storage = match (params.get_frequency(), n_threads) {
            // shared storage is not needed if only one thread is used
            // (there is no other thread to share data with, duh)
            (Frequency::Every(_) | Frequency::Once, 1) => false,
            // shared storage is not needed if the frequency is performed for every analyzed frame
            // (each thread handles lipid assignment locally and no data need to be shared)
            (Frequency::Every(x), _) if x.get() == 1 => false,
            // shared storage is needed in all other cases
            _ => true,
        };

        match params {
            LeafletClassification::Global(_) => Self::Global(
                GlobalClassification {
                    heads: Vec::new(),
                    membrane_center: Vector3D::new(0.0, 0.0, 0.0),
                    membrane_normal,
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
            LeafletClassification::Local(_) => Self::Local(
                LocalClassification {
                    heads: Vec::new(),
                    radius: params.get_radius().expect(PANIC_MESSAGE),
                    membrane_center: Vec::new(),
                    membrane_normal,
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
            LeafletClassification::Individual(_) => Self::Individual(
                IndividualClassification {
                    heads: Vec::new(),
                    methyls: Vec::new(),
                    membrane_normal,
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
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
            Self::Global(x, _) => {
                x.insert(molecule, system)?;
            }
            Self::Local(x, _) => {
                x.insert(molecule, system)?;
            }
            Self::Individual(x, _) => {
                x.insert(molecule, system)?;
            }
        }

        Ok(())
    }

    /// Calculate the number of molecules assigned to the upper and to the lower leaflet.
    pub(super) fn statistics(&self) -> (usize, usize) {
        match self {
            Self::Global(_, y) | Self::Local(_, y) | Self::Individual(_, y) => {
                y.calc_assignment_statistics()
            }
        }
    }

    /// Get frequency at which the assignment should be performed.
    #[inline(always)]
    fn get_frequency(&self) -> Frequency {
        match self {
            Self::Global(x, _) => x.frequency,
            Self::Local(x, _) => x.frequency,
            Self::Individual(x, _) => x.frequency,
        }
    }

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

    /// Determine whether assignment should be performed for this trajectory frame.
    #[inline(always)]
    fn should_assign(&self, current_frame: usize) -> bool {
        match self.get_frequency() {
            Frequency::Once if current_frame == 0 => true,
            Frequency::Every(n) if current_frame % n == 0 => true,
            _ => false,
        }
    }

    /// Assign all lipids into their respective leaflets.
    #[inline]
    pub(super) fn assign_lipids(
        &mut self,
        system: &System,
        current_frame: usize,
        membrane_center: &OnceCell<Vector3D>, // only used for the `global` classification method
    ) -> Result<(), AnalysisError> {
        if !self.should_assign(current_frame) {
            return Ok(());
        }

        match self.get_frequency() {
            Frequency::Once if current_frame == 0 => (), // continue
            Frequency::Every(n) if current_frame % n == 0 => (), // continue
            _ => return Ok(()),                          // perform no assignment
        }

        match self {
            MoleculeLeafletClassification::Global(x, y) => {
                // calculate global membrane center of mass
                let center = membrane_center.get_or_try_init(|| {
                    system
                        .group_get_center(group_name!("Membrane"))
                        .map_err(|_| AnalysisError::InvalidGlobalMembraneCenter)
                })?;

                x.set_membrane_center(center.clone());
                y.assign_lipids(system, x, current_frame)
            }
            MoleculeLeafletClassification::Local(x, y) => {
                x.set_membrane_center(system, x.membrane_normal)?;
                y.assign_lipids(system, x, current_frame)
            }
            MoleculeLeafletClassification::Individual(x, y) => {
                y.assign_lipids(system, x, current_frame)
            }
        }
    }

    /// Get the leaflet a molecule with target index is located in.
    /// This performs no calculation, this only returns the already calculated leaflet assignment.
    #[inline(always)]
    pub(super) fn get_assigned_leaflet(
        &mut self,
        molecule_index: usize,
        current_frame: usize,
    ) -> Leaflet {
        match self {
            MoleculeLeafletClassification::Global(x, y) => {
                y.get_assigned_leaflet(molecule_index, current_frame, x.frequency)
            }
            MoleculeLeafletClassification::Local(x, y) => {
                y.get_assigned_leaflet(molecule_index, current_frame, x.frequency)
            }
            MoleculeLeafletClassification::Individual(x, y) => {
                y.get_assigned_leaflet(molecule_index, current_frame, x.frequency)
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
    /// Frequency with which the assignment should be performed.
    /// Note that this is a 'real frequency' (input frequency multiplied by the step_size).
    frequency: Frequency,
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
    /// Frequency with which the assignment should be performed.
    /// Note that this is a 'real frequency' (input frequency multiplied by the step_size).
    frequency: Frequency,
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
    /// Frequency with which the assignment should be performed.
    /// Note that this is a 'real frequency' (input frequency multiplied by the step_size).
    frequency: Frequency,
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

/// Vector of leaflet assignments for each molecule of the given type.
#[derive(Debug, Clone)]
pub(super) struct AssignedLeaflets {
    /// Locally stored leaflet assignment.
    local: Option<Vec<Leaflet>>,
    /// Leaflet assignment for each relevant frame. Shared across all threads.
    /// `None` if not needed, i.e. if the analysis is not multithreaded or if the frequency of leaflet assignment is 1.
    shared: Option<SharedAssignedLeaflets>,
    /// Index of the frame the data in `local` correspond to.
    local_frame: Option<usize>,
}

impl AssignedLeaflets {
    /// Create a new structure for storing information about positions of lipids in leaflets.
    /// `shared` specifies whether shared storage should be set-up.
    /// Use `true` if the analysis is multithreaded and the leaflet assignment is NOT performed every thread.
    /// Use `false` in all other cases.
    #[inline(always)]
    fn new(shared: bool) -> Self {
        if shared {
            AssignedLeaflets {
                local: None,
                shared: Some(SharedAssignedLeaflets::default()),
                local_frame: None,
            }
        } else {
            AssignedLeaflets {
                local: None,
                shared: None,
                local_frame: None,
            }
        }
    }

    /// Assign all lipids into membrane leaflets.
    fn assign_lipids(
        &mut self,
        system: &System,
        classifier: &impl LeafletClassifier,
        current_frame: usize,
    ) -> Result<(), AnalysisError> {
        self.local = Some(
            (0..classifier.n_molecules())
                .map(|index| classifier.identify_leaflet(system, index))
                .collect::<Result<Vec<_>, _>>()?,
        );

        self.local_frame = Some(current_frame);

        // copy the current assignment to shared assignment, so other threads can also access it
        // this should only be done if the number of threads is higher than 1 and the input frequency is NOT 1
        if let Some(shared) = self.shared.as_mut() {
            shared.copy_to(
                self.local_frame.expect(PANIC_MESSAGE),
                self.local.as_ref().expect(PANIC_MESSAGE),
            )
        }

        Ok(())
    }

    /// Get the leaflet that was assigned to molecule of target index.
    fn get_assigned_leaflet(
        &mut self,
        molecule_index: usize,
        current_frame: usize,
        assignment_frequency: Frequency,
    ) -> Leaflet {
        match assignment_frequency {
            Frequency::Every(n) => {
                // check whether a suitable local assignment is accessible
                let closest_frame = (current_frame / n.get()) * n.get();
                match self.local_frame {
                    // get leaflet assignment from the locally stored data
                    Some(x) if x == closest_frame => self.get_leaflet_from_local(molecule_index),
                    // copy assignment from shared data into locally stored data
                    Some(_) | None => {
                        self.copy_from_shared_assignment(closest_frame);
                        self.get_leaflet_from_local(molecule_index)
                    }
                }
            }
            Frequency::Once => {
                match self.local_frame {
                    Some(0) => self.get_leaflet_from_local(molecule_index),
                    Some(x) => panic!(
                        "FATAL GORDER ERROR | AssignLeaflets::get_assigned_leaflet | Local data is from frame `{}` but frequency is Once. (Why was assignment performed for this frame?) {}", 
                        x, PANIC_MESSAGE
                    ),
                    None => {
                        self.copy_from_shared_assignment(0);
                        self.get_leaflet_from_local(molecule_index)
                    }
                }
            }
        }
    }

    /// Get the leaflet that was assigned to molecule of target index from the local assignment storage.
    #[inline]
    fn get_leaflet_from_local(&self, molecule_index: usize) -> Leaflet {
        *self.local
            .as_ref()
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | AssignedLeaflets::get_leaflet_from_local | Molecule not assigned. Local assignment should exist. {}", 
                PANIC_MESSAGE)
            ).get(molecule_index)
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | AssignedLeaflets::get_leaflet_from_local | Molecule not found. Molecule with internal `gorder` index `{}` should exist. {}",
                molecule_index, PANIC_MESSAGE))
    }

    /// Copy the leaflet assignment for the specified frame into storage local for the thread.
    ///
    /// ## Panic
    /// Panics if the `shared` storage is `None`.
    #[inline]
    fn copy_from_shared_assignment(&mut self, frame: usize) {
        self.local = Some(self
            .shared
            .as_ref()
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | AssignedLeaflets::copy_from_shared_assignment | Data for frame `{}` requested from shared storage, but shared storage has not been set-up. {}", 
                frame, PANIC_MESSAGE))
            .copy_from(frame));

        self.local_frame = Some(frame);
    }

    /// Returns the number of molecules in the upper and the lower leaflet.
    ///
    /// ## Panic
    /// Panics if `local` is `None`.
    fn calc_assignment_statistics(&self) -> (usize, usize) {
        self.local
            .as_ref()
            .expect(PANIC_MESSAGE)
            .iter()
            .fold((0, 0), |(upper, lower), mol| match mol {
                Leaflet::Upper => (upper + 1, lower),
                Leaflet::Lower => (upper, lower + 1),
            })
    }
}

/// Stores leaflet assignment for each relevant frame. Shared among all threads.
#[derive(Debug, Clone, Default)]
pub(super) struct SharedAssignedLeaflets(Arc<Mutex<HashMap<usize, Vec<Leaflet>>>>);

impl SharedAssignedLeaflets {
    /// Copy the leaflet assignment for the specified frame.
    /// This requires waiting for the thread responsible for performing the leaflet assignment
    /// for the `frame_to_fetch` to finish the analysis of this frame.
    fn copy_from(&self, frame: usize) -> Vec<Leaflet> {
        let start_time = Instant::now();
        let mut warning_logged = false;

        // spin-lock: waiting for the requested frame to be available
        loop {
            let shared_data = self.0.lock();
            let assignment = shared_data.get(&frame);
            if let Some(assign) = assignment {
                return assign.clone();
            }

            // defensive check for a deadlock
            if start_time.elapsed() > *TIMEOUT {
                if !warning_logged {
                    log::warn!("DEADLOCKED? Thread has been waiting for shared leaflet assignment data (frame '{}') for more than {} seconds.
This may be due to resource contention or a bug. Ensure that your CPU is not oversubscribed and that you have not lost access to the trajectory file.
If `gorder` is causing oversubscription, reduce the number of threads used for the analysis.
If other computationally intensive software is running alongside `gorder`, consider terminating it.
If the issue persists, please report it by opening an issue at `github.com/Ladme/gorder/issues` or sending an email to `ladmeb@gmail.com`. 
(Note: If no progress is made, this thread will terminate in {} seconds to prevent resource exhaustion.)",
                    frame,
                    TIMEOUT_SECONDS,
                    HARD_TIMEOUT_SECONDS - TIMEOUT_SECONDS,
                );
                    warning_logged = true;
                }

                if start_time.elapsed() > *HARD_TIMEOUT {
                    panic!("FATAL GORDER ERROR | SharedAssignedLeaflets::copy_from | Deadlock. Could not get shared data for leaflet assignment. Spent more than `{}` seconds inside the spin-lock. {}", 
                    HARD_TIMEOUT_SECONDS, PANIC_MESSAGE)
                }
            }

            // shared data unlock here
        }
    }

    /// Copy the leaflet assignment to shared storage.
    fn copy_to(&mut self, frame: usize, assignment: &[Leaflet]) {
        let mut shared_data = self.0.lock();
        let previous = shared_data.insert(frame, assignment.to_owned());

        // defensive check; no other thread should have read the same frame
        assert!(previous.is_none(), "FATAL GORDER ERROR | SharedAssignedLeaflets::copy_to | Leaflet assignment for frame index `{}` already exists, but it should not. {}", 
            frame, PANIC_MESSAGE
        );
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
            1,
            1,
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
            1,
            1,
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
            1,
            1,
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
        let mut molecule_classifier =
            MoleculeLeafletClassification::new(&classifier, Dimension::Z, 1, 1);
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

        let mut molecule_classifier =
            MoleculeLeafletClassification::new(&classifier, Dimension::Z, 1, 1);
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
            1,
            1,
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
            1,
            1,
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
            1,
            1,
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

    macro_rules! test_classification_new {
        ($name:ident, $classification:expr, $variant:path, $frequency:pat) => {
            #[test]
            fn $name() {
                let single_threaded =
                    MoleculeLeafletClassification::new(&$classification, Dimension::Z, 1, 2);
                match single_threaded {
                    $variant(x, y) => {
                        assert!(matches!(x.frequency, $frequency));
                        assert!(y.shared.is_none());
                        assert!(y.local.is_none());
                        assert!(y.local_frame.is_none());
                    }
                    _ => panic!("Incorrect structure constructed."),
                }

                let multi_threaded =
                    MoleculeLeafletClassification::new(&$classification, Dimension::Z, 4, 2);
                match multi_threaded {
                    $variant(x, y) => {
                        assert!(matches!(x.frequency, $frequency));
                        match x.frequency {
                            Frequency::Every(n) if n.get() == 2 => assert!(y.shared.is_none()),
                            Frequency::Every(n) if n.get() == 10 => {
                                assert!(y.shared.is_some());
                            }
                            Frequency::Every(n) => panic!("Unexpected real frequency '{}': neither 2 (init 1) nor 10 (init 5) detected.", n),
                            _ => assert!(y.shared.is_some()),
                        }
                        assert!(y.local.is_none());
                        assert!(y.local_frame.is_none());
                    }
                    _ => panic!("Incorrect structure constructed."),
                }
            }
        };
    }

    test_classification_new!(
        test_global_leaflet_classification_new_once,
        LeafletClassification::global("@membrane", "name P").with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Global,
        Frequency::Once
    );

    test_classification_new!(
        test_local_leaflet_classification_new_once,
        LeafletClassification::local("@membrane", "name P", 2.5).with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Local,
        Frequency::Once
    );

    test_classification_new!(
        test_individual_leaflet_classification_new_once,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Individual,
        Frequency::Once
    );

    test_classification_new!(
        test_global_leaflet_classification_new_every_five,
        LeafletClassification::global("@membrane", "name P")
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_)
    );

    test_classification_new!(
        test_local_leaflet_classification_new_every_five,
        LeafletClassification::local("@membrane", "name P", 2.5)
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Local,
        Frequency::Every(_)
    );

    test_classification_new!(
        test_individual_leaflet_classification_new_every_five,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Individual,
        Frequency::Every(_)
    );

    test_classification_new!(
        test_global_leaflet_classification_new_every_one,
        LeafletClassification::global("@membrane", "name P")
            .with_frequency(Frequency::every(1).unwrap()),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_)
    );

    test_classification_new!(
        test_local_leaflet_classification_new_every_one,
        LeafletClassification::local("@membrane", "name P", 2.5)
            .with_frequency(Frequency::every(1).unwrap()),
        MoleculeLeafletClassification::Local,
        Frequency::Every(_)
    );

    test_classification_new!(
        test_individual_leaflet_classification_new_every_one,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_frequency(Frequency::every(1).unwrap()),
        MoleculeLeafletClassification::Individual,
        Frequency::Every(_)
    );

    #[test]
    fn test_should_assign_every() {
        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316"),
            Dimension::Z,
            1,
            1,
        );

        for i in 0..50 {
            assert!(classifier.should_assign(i));
        }

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            Dimension::Z,
            1,
            5,
        );

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(!classifier.should_assign(7));
        assert!(classifier.should_assign(10));
        assert!(classifier.should_assign(15));
        assert!(!classifier.should_assign(16));

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5),
            Dimension::Z,
            3,
            5,
        );

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(!classifier.should_assign(7));
        assert!(classifier.should_assign(10));
        assert!(classifier.should_assign(15));
        assert!(!classifier.should_assign(16));

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P")
                .with_frequency(Frequency::every(4).unwrap()),
            Dimension::Z,
            1,
            1,
        );

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(classifier.should_assign(8));
        assert!(!classifier.should_assign(9));
        assert!(!classifier.should_assign(10));
        assert!(!classifier.should_assign(11));
        assert!(classifier.should_assign(12));
        assert!(classifier.should_assign(40));

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316")
                .with_frequency(Frequency::every(5).unwrap()),
            Dimension::Z,
            1,
            7,
        );

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(!classifier.should_assign(5));
        assert!(!classifier.should_assign(7));
        assert!(classifier.should_assign(35));
        assert!(!classifier.should_assign(50));
        assert!(classifier.should_assign(70));
        assert!(!classifier.should_assign(244));
        assert!(classifier.should_assign(245));

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5)
                .with_frequency(Frequency::every(2).unwrap()),
            Dimension::Z,
            6,
            11,
        );

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(!classifier.should_assign(2));
        assert!(!classifier.should_assign(6));
        assert!(!classifier.should_assign(11));
        assert!(classifier.should_assign(22));
        assert!(!classifier.should_assign(33));
        assert!(classifier.should_assign(44));
        assert!(!classifier.should_assign(813));
        assert!(classifier.should_assign(814));
        assert!(!classifier.should_assign(815));
    }

    #[test]
    fn test_should_assign_once() {
        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P").with_frequency(Frequency::once()),
            Dimension::Z,
            1,
            1,
        );

        assert!(classifier.should_assign(0));
        for i in 1..51 {
            assert!(!classifier.should_assign(i));
        }

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316")
                .with_frequency(Frequency::once()),
            Dimension::Z,
            1,
            7,
        );

        assert!(classifier.should_assign(0));
        for i in 1..51 {
            assert!(!classifier.should_assign(i));
        }

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5)
                .with_frequency(Frequency::once()),
            Dimension::Z,
            6,
            11,
        );

        assert!(classifier.should_assign(0));
        for i in 1..51 {
            assert!(!classifier.should_assign(i));
        }
    }
}

#[cfg(test)]
mod tests_assigned_leaflets {
    use super::*;

    #[test]
    fn test_get_assigned_leaflet_local_no_shared() {
        let mut data = AssignedLeaflets {
            local: Some(vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper]),
            local_frame: Some(0),
            shared: None,
        };

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::every(4).unwrap());
        assert_eq!(assignment, Leaflet::Lower);
        assert!(data.shared.is_none());

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::once());
        assert_eq!(assignment, Leaflet::Lower);
        assert!(data.shared.is_none());
    }

    #[test]
    fn test_get_assigned_leaflet_local_and_shared() {
        // should be ignored since local data is loaded
        let shared = SharedAssignedLeaflets(Arc::new(Mutex::new(
            [(0, vec![Leaflet::Upper, Leaflet::Upper, Leaflet::Upper])].into(),
        )));

        let mut data = AssignedLeaflets {
            local: Some(vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper]),
            local_frame: Some(0),
            shared: Some(shared),
        };

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::every(4).unwrap());
        assert_eq!(assignment, Leaflet::Lower);

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::once());
        assert_eq!(assignment, Leaflet::Lower);
    }

    #[test]
    fn test_get_assigned_leaflet_shared_only() {
        let shared = SharedAssignedLeaflets(Arc::new(Mutex::new(
            [(0, vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper])].into(),
        )));

        let mut data = AssignedLeaflets {
            local: None,
            local_frame: None,
            shared: Some(shared),
        };

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::every(4).unwrap());
        assert_eq!(assignment, Leaflet::Lower);
        let local = data.local.as_ref().unwrap();
        assert_eq!(local[0], Leaflet::Upper);
        assert_eq!(local[1], Leaflet::Lower);
        assert_eq!(local[2], Leaflet::Upper);
        assert_eq!(data.local_frame.unwrap(), 0);

        // reset local data
        data.local = None;
        data.local_frame = None;

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::once());
        assert_eq!(assignment, Leaflet::Lower);
        let local = data.local.as_ref().unwrap();
        assert_eq!(local[0], Leaflet::Upper);
        assert_eq!(local[1], Leaflet::Lower);
        assert_eq!(local[2], Leaflet::Upper);
        assert_eq!(data.local_frame.unwrap(), 0);

        // if shared data change, `get_assigne_leaflet` should still provide the same assignment since it's already loaded to local
        data.shared = None;

        let assignment = data.get_assigned_leaflet(1, 0, Frequency::every(1).unwrap());
        assert_eq!(assignment, Leaflet::Lower);
        let local = data.local.as_ref().unwrap();
        assert_eq!(local[0], Leaflet::Upper);
        assert_eq!(local[1], Leaflet::Lower);
        assert_eq!(local[2], Leaflet::Upper);
        assert_eq!(data.local_frame.unwrap(), 0);
    }

    #[test]
    fn test_copy_from_shared() {
        let shared = SharedAssignedLeaflets(Arc::new(Mutex::new(
            [
                (0, vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper]),
                (10, vec![Leaflet::Lower, Leaflet::Upper, Leaflet::Lower]),
                (21, vec![Leaflet::Upper, Leaflet::Upper, Leaflet::Lower]),
            ]
            .into(),
        )));

        let mut data = AssignedLeaflets {
            local: None,
            local_frame: None,
            shared: Some(shared),
        };

        data.copy_from_shared_assignment(10);
        assert_eq!(data.get_leaflet_from_local(2), Leaflet::Lower);
        assert_eq!(data.local_frame.unwrap(), 10);

        data.copy_from_shared_assignment(0);
        assert_eq!(data.get_leaflet_from_local(2), Leaflet::Upper);
        assert_eq!(data.local_frame.unwrap(), 0);

        data.copy_from_shared_assignment(21);
        assert_eq!(data.get_leaflet_from_local(2), Leaflet::Lower);
        assert_eq!(data.local_frame.unwrap(), 21);
    }

    fn prepare_data(use_shared: bool) -> (System, GlobalClassification, AssignedLeaflets) {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            Dimension::Z,
            1,
            1,
        );

        let membrane_center = system
            .selection_iter("@membrane")
            .unwrap()
            .get_center()
            .unwrap();

        let mut extracted_classifier = match classifier {
            MoleculeLeafletClassification::Global(x, _) => x,
            _ => panic!("Unexpected classification method."),
        };

        extracted_classifier.heads_mut().push(1385);
        extracted_classifier.heads_mut().push(11885);
        extracted_classifier.set_membrane_center(membrane_center);

        let data = AssignedLeaflets {
            local: None,
            local_frame: None,
            shared: if use_shared {
                Some(SharedAssignedLeaflets::default())
            } else {
                None
            },
        };

        (system, extracted_classifier, data)
    }

    #[test]
    fn test_assign_lipids_no_shared() {
        let (system, classifier, mut leaflets) = prepare_data(false);

        leaflets.assign_lipids(&system, &classifier, 15).unwrap();

        assert_eq!(leaflets.local_frame.unwrap(), 15);
        assert_eq!(leaflets.get_leaflet_from_local(0), Leaflet::Upper);
        assert_eq!(leaflets.get_leaflet_from_local(1), Leaflet::Lower);
        assert!(leaflets.shared.is_none());
    }

    #[test]
    fn test_assign_lipids_with_shared() {
        let (system, classifier, mut leaflets) = prepare_data(true);

        leaflets.assign_lipids(&system, &classifier, 15).unwrap();

        assert_eq!(leaflets.local_frame.unwrap(), 15);
        assert_eq!(leaflets.get_leaflet_from_local(0), Leaflet::Upper);
        assert_eq!(leaflets.get_leaflet_from_local(1), Leaflet::Lower);

        let mutex = leaflets.shared.as_ref().unwrap().0.lock();
        let shared = mutex.get(&15).unwrap();
        assert_eq!(shared[0], Leaflet::Upper);
        assert_eq!(shared[1], Leaflet::Lower);
    }
}

#[cfg(test)]
mod tests_shared_assigned_leaflets {
    use std::thread;

    use super::*;

    #[test]
    fn test_copy_from() {
        let shared = SharedAssignedLeaflets(Arc::new(Mutex::new(
            [(100, vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper])].into(),
        )));

        let local = shared.copy_from(100);
        assert_eq!(local[0], Leaflet::Upper);
        assert_eq!(local[1], Leaflet::Lower);
        assert_eq!(local[2], Leaflet::Upper);
    }

    #[test]
    fn test_copy_to() {
        let local = vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper];
        let mut shared = SharedAssignedLeaflets::default();
        shared.copy_to(100, &local);

        let local2 = shared.copy_from(100);
        assert_eq!(local2[0], Leaflet::Upper);
        assert_eq!(local2[1], Leaflet::Lower);
        assert_eq!(local2[2], Leaflet::Upper);
    }

    #[test]
    fn test_producer_consumer() {
        let shared = SharedAssignedLeaflets::default();
        let mut shared_clone = shared.clone();

        // the consumer thread
        let consumer = thread::spawn(move || {
            let local = shared.copy_from(100);
            assert_eq!(local[0], Leaflet::Upper);
            assert_eq!(local[1], Leaflet::Lower);
            assert_eq!(local[2], Leaflet::Upper);
        });

        // the producer thread
        let producer = thread::spawn(move || {
            let local = vec![Leaflet::Upper, Leaflet::Lower, Leaflet::Upper];
            thread::sleep(Duration::from_millis(100)); // simulate work
            shared_clone.copy_to(100, &local);
        });

        producer.join().expect("Producer thread panicked");
        consumer.join().expect("Consumer thread panicked");
    }
}
