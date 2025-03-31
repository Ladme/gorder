// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use core::f32;
use std::{
    cmp::Ordering,
    collections::VecDeque,
    fs::read_to_string,
    num::NonZeroUsize,
    sync::Arc,
    time::{Duration, Instant},
};

use super::{
    common::{create_group, get_reference_head, macros::group_name},
    pbc::PBCHandler,
    topology::SystemTopology,
};
use crate::{
    analysis::topology::molecule::handle_moltypes,
    errors::{
        AnalysisError, ClusterError, ConfigError, ManualLeafletClassificationError,
        NdxLeafletClassificationError, TopologyError,
    },
    input::{Frequency, LeafletClassification, MembraneNormal},
    Leaflet, PANIC_MESSAGE,
};
use getset::{CopyGetters, Getters, MutGetters, Setters};
use groan_rs::{
    errors::{AtomError, PositionError},
    prelude::{Atom, Dimension, Groups, Vector3D},
    structures::group::Group,
    system::System,
};
use hashbrown::{HashMap, HashSet};
use parking_lot::Mutex;

use once_cell::{sync::Lazy, unsync::OnceCell};

/// [`TIMEOUT`] in seconds.
const TIMEOUT_SECONDS: u64 = 5;
/// Global soft timeout duration for spin-lock used when fetching data for leaflet assignment.
/// After this time a warning is logged.
static TIMEOUT: Lazy<Duration> = Lazy::new(|| Duration::from_secs(TIMEOUT_SECONDS));

/// [`HARD_TIMEOUT`] in seconds.
const HARD_TIMEOUT_SECONDS: u64 = 125;
/// Global HARD timeout duration for spin-lock used when fetching data for leaflet assignment.
/// After this time a PANIC is raised.
static HARD_TIMEOUT: Lazy<Duration> = Lazy::new(|| Duration::from_secs(HARD_TIMEOUT_SECONDS));

/// Relative number of lipids molecules that must remain in the same leaflet
/// between two consecutive trajectory frames analyzed by the same thread
/// for reliable leaflet match.
/// If this is not maintained, an error is raised.
const CLUSTER_CLASSIFICATION_LIMIT: f32 = 0.8;

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
            Self::FromNdx(params) => {
                create_group(system, "Heads", params.heads())?;
            }
            Self::FromFile(_) | Self::FromMap(_) => (),
            Self::Clustering(params) => {
                create_group(system, "ClusterHeads", params.heads())?;
            }
        }

        Ok(())
    }
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
#[allow(private_interfaces)]
pub(crate) enum MoleculeLeafletClassification {
    Global(GlobalClassification, AssignedLeaflets),
    Local(LocalClassification, AssignedLeaflets),
    Individual(IndividualClassification, AssignedLeaflets),
    Manual(ManualClassification, AssignedLeaflets),
    ManualNdx(NdxClassification, AssignedLeaflets),
    Clustering(ClusterClassification, AssignedLeaflets),
}

impl MoleculeLeafletClassification {
    /// Convert the input `LeafletClassification` into an enum that is used in the analysis.
    pub(super) fn new(
        params: &LeafletClassification,
        membrane_normal: &MembraneNormal,
        n_threads: usize,
        step_size: usize,
    ) -> Result<Self, ConfigError> {
        let needs_shared_storage = match (params.get_frequency(), n_threads) {
            // shared storage is not needed if only one thread is used
            // (there is no other thread to share data with, duh)
            (Frequency::Every(_) | Frequency::Once, 1) => false,
            // shared storage is not needed if the assignment is performed for every analyzed frame
            // (each thread handles lipid assignment locally and no data need to be shared)
            (Frequency::Every(x), _) if x.get() == 1 => false,
            // shared storage is needed in all other cases
            _ => true,
        };

        fn get_membrane_normal(
            params: &LeafletClassification,
            membrane_normal: &MembraneNormal,
        ) -> Result<Dimension, ConfigError> {
            match (params.get_membrane_normal(), membrane_normal) {
                (None, MembraneNormal::Static(y)) => Ok((*y).into()),
                (
                    None,
                    MembraneNormal::Dynamic(_)
                    | MembraneNormal::FromFile(_)
                    | MembraneNormal::FromMap(_),
                ) => Err(ConfigError::MissingMembraneNormal),
                (Some(x), _) => Ok(x.into()),
            }
        }

        let classification = match params {
            LeafletClassification::Global(_) => Self::Global(
                GlobalClassification {
                    heads: Vec::new(),
                    membrane_center: Vector3D::new(0.0, 0.0, 0.0),
                    membrane_normal: get_membrane_normal(params, membrane_normal)?,
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
                    membrane_normal: get_membrane_normal(params, membrane_normal)?,
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
            LeafletClassification::Individual(_) => Self::Individual(
                IndividualClassification {
                    heads: Vec::new(),
                    methyls: Vec::new(),
                    membrane_normal: get_membrane_normal(params, membrane_normal)?,
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
            LeafletClassification::FromFile(_) | LeafletClassification::FromMap(_) => Self::Manual(
                ManualClassification {
                    assignment: None,
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
            LeafletClassification::FromNdx(ndx_params) => Self::ManualNdx(
                NdxClassification {
                    ndx: ndx_params.ndx().clone(),
                    heads: Vec::new(),
                    groups: None,
                    last_assigned_frame: None,
                    upper_leaflet: ndx_params.upper_leaflet().clone(),
                    lower_leaflet: ndx_params.lower_leaflet().clone(),
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
            LeafletClassification::Clustering(cluster_params) => Self::Clustering(
                ClusterClassification {
                    heads: Vec::new(),
                    radius: cluster_params.radius(),
                    min_samples: cluster_params.min_samples(),
                    clusters: once_cell::sync::OnceCell::new(),
                    shared_clusters: Arc::new(Mutex::new(HashMap::new())),
                    frequency: params.get_frequency()
                        * NonZeroUsize::new(step_size).expect(PANIC_MESSAGE),
                },
                AssignedLeaflets::new(needs_shared_storage),
            ),
        };

        Ok(classification)
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
            Self::ManualNdx(x, _) => {
                x.insert(molecule, system)?;
            }
            Self::Clustering(x, _) => {
                x.insert(molecule, system)?;
            }
            // do nothing; manual "leaflet assignment file" classifier is not set-up like this
            Self::Manual(_, _) => (),
        }

        Ok(())
    }

    /// Calculate the number of molecules assigned to the upper and to the lower leaflet.
    pub(super) fn statistics(&self) -> (usize, usize) {
        match self {
            Self::Global(_, y)
            | Self::Local(_, y)
            | Self::Individual(_, y)
            | Self::Manual(_, y)
            | Self::ManualNdx(_, y)
            | Self::Clustering(_, y) => y.calc_assignment_statistics(),
        }
    }

    /// Get frequency at which the assignment should be performed.
    #[inline(always)]
    fn get_frequency(&self) -> Frequency {
        match self {
            Self::Global(x, _) => x.frequency,
            Self::Local(x, _) => x.frequency,
            Self::Individual(x, _) => x.frequency,
            Self::Manual(x, _) => x.frequency,
            Self::ManualNdx(x, _) => x.frequency,
            Self::Clustering(x, _) => x.frequency,
        }
    }

    /// Initialize the reading of a new frame.
    #[inline(always)]
    pub(super) fn init_new_frame(&mut self) {
        match self {
            Self::Global(_, _)
            | Self::Local(_, _)
            | Self::Individual(_, _)
            | Self::Manual(_, _)
            | Self::ManualNdx(_, _) => (),
            Self::Clustering(x, _) => x.clusters = once_cell::sync::OnceCell::new(),
        }
    }

    /// Identify leaflet in which the molecule with the specified index is located.
    #[inline(always)]
    #[allow(unused)]
    fn identify_leaflet<'a>(
        &mut self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        molecule_index: usize,
        current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        match self {
            MoleculeLeafletClassification::Global(x, _) => {
                x.identify_leaflet(system, pbc_handler, molecule_index, current_frame)
            }
            MoleculeLeafletClassification::Local(x, _) => {
                x.identify_leaflet(system, pbc_handler, molecule_index, current_frame)
            }
            MoleculeLeafletClassification::Individual(x, _) => {
                x.identify_leaflet(system, pbc_handler, molecule_index, current_frame)
            }
            MoleculeLeafletClassification::Manual(x, _) => {
                x.identify_leaflet(system, pbc_handler, molecule_index, current_frame)
            }
            MoleculeLeafletClassification::ManualNdx(x, _) => {
                x.identify_leaflet(system, pbc_handler, molecule_index, current_frame)
            }
            MoleculeLeafletClassification::Clustering(x, _) => {
                x.identify_leaflet(system, pbc_handler, molecule_index, current_frame)
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
    pub(super) fn assign_lipids<'a>(
        &mut self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        current_frame: usize,
        membrane_center: &OnceCell<Vector3D>, // only used for the `global` classification method
    ) -> Result<(), AnalysisError> {
        if !self.should_assign(current_frame) {
            return Ok(());
        }

        match self {
            MoleculeLeafletClassification::Global(x, y) => {
                // calculate global membrane center of mass
                let center = membrane_center.get_or_try_init(|| {
                    pbc_handler
                        .group_get_center(system, group_name!("Membrane"))
                        .map_err(|_| AnalysisError::InvalidGlobalMembraneCenter)
                })?;

                if center.x.is_nan() || center.y.is_nan() || center.z.is_nan() {
                    return Err(AnalysisError::InvalidGlobalMembraneCenter);
                }

                x.set_membrane_center(center.clone());
                y.assign_lipids(system, pbc_handler, x, current_frame)
            }
            MoleculeLeafletClassification::Local(x, y) => {
                x.set_membrane_center(system, pbc_handler, x.membrane_normal)?;
                y.assign_lipids(system, pbc_handler, x, current_frame)
            }
            MoleculeLeafletClassification::Individual(x, y) => {
                y.assign_lipids(system, pbc_handler, x, current_frame)
            }
            MoleculeLeafletClassification::Manual(x, y) => {
                y.assign_lipids(system, pbc_handler, x, current_frame)
            }
            MoleculeLeafletClassification::ManualNdx(x, y) => {
                y.assign_lipids(system, pbc_handler, x, current_frame)
            }
            MoleculeLeafletClassification::Clustering(x, y) => {
                y.assign_lipids(system, pbc_handler, x, current_frame)
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
            MoleculeLeafletClassification::Manual(x, y) => {
                y.get_assigned_leaflet(molecule_index, current_frame, x.frequency)
            }
            MoleculeLeafletClassification::ManualNdx(x, y) => {
                y.get_assigned_leaflet(molecule_index, current_frame, x.frequency)
            }
            MoleculeLeafletClassification::Clustering(x, y) => {
                y.get_assigned_leaflet(molecule_index, current_frame, x.frequency)
            }
        }
    }
}

/// Trait implemented by all leaflet classification methods.
trait LeafletClassifier {
    /// Caclulate membrane leaflet the specified molecule belongs to.
    fn identify_leaflet<'a>(
        &mut self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        molecule_index: usize,
        current_frame: usize,
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
        self.heads
            .push(get_reference_head(molecule, system, group_name!("Heads"))?);
        Ok(())
    }
}

impl LeafletClassifier for GlobalClassification {
    #[inline(always)]
    fn identify_leaflet<'a>(
        &mut self,
        system: &System,
        pbc_handler: &impl PBCHandler<'a>,
        molecule_index: usize,
        _current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        common_identify_leaflet(
            &self.heads,
            molecule_index,
            system,
            pbc_handler,
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
        self.heads
            .push(get_reference_head(molecule, system, group_name!("Heads"))?);
        self.membrane_center.push(Vector3D::new(0.0, 0.0, 0.0));
        Ok(())
    }

    /// Calculate and set local membrane center of geometry for each molecule.
    #[inline(always)]
    pub(super) fn set_membrane_center<'a>(
        &mut self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        membrane_normal: Dimension,
    ) -> Result<(), AnalysisError> {
        self.membrane_center = pbc_handler.calc_local_membrane_centers(
            system,
            &self.heads,
            self.radius,
            membrane_normal,
        )?;

        Ok(())
    }
}

impl LeafletClassifier for LocalClassification {
    #[inline(always)]
    fn identify_leaflet<'a>(
        &mut self,
        system: &System,
        pbc_handler: &impl PBCHandler<'a>,
        molecule_index: usize,
        _current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        common_identify_leaflet(
            &self.heads,
            molecule_index,
            system,
            pbc_handler,
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
fn common_identify_leaflet<'a>(
    heads: &[usize],
    molecule_index: usize,
    system: &System,
    pbc_handler: &impl PBCHandler<'a>,
    membrane_center: &Vector3D,
    membrane_normal: Dimension,
) -> Result<Leaflet, AnalysisError> {
    let head_index = *heads.get(molecule_index).expect(PANIC_MESSAGE);
    let head = unsafe { system.get_atom_unchecked(head_index) };

    let head_pos = head
        .get_position()
        .ok_or_else(|| AnalysisError::UndefinedPosition(head_index))?;
    let distance = pbc_handler.distance(head_pos, membrane_center, membrane_normal);

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
        self.heads
            .push(get_reference_head(molecule, system, group_name!("Heads"))?);
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
    fn identify_leaflet<'a>(
        &mut self,
        system: &System,
        pbc_handler: &impl PBCHandler<'a>,
        molecule_index: usize,
        _current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        let head_index = self.heads.get(molecule_index).expect(PANIC_MESSAGE);

        let mut total_distance = 0.0;
        for methyl_index in self.methyls.get(molecule_index).expect(PANIC_MESSAGE) {
            total_distance += match pbc_handler.atoms_distance(system, *head_index, *methyl_index, self.membrane_normal) {
                Ok(x) => x,
                Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => return Err(AnalysisError::UndefinedPosition(x)),
                Err(AtomError::OutOfRange(x)) => panic!("FATAL GORDER ERROR | IndividualClassification::identify_leaflet | Index '{}' out of range. {}", x, PANIC_MESSAGE),
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

/// Leaflet classification method where lipids are assigned manually.
#[derive(Debug, Clone, Getters, MutGetters)]
pub(crate) struct ManualClassification {
    /// Lipid assignment. Each inner vector corresponds to one frame.
    /// This is only set up later, since when the molecules are being constructed, we don't know their final names.
    assignment: Option<Vec<Vec<Leaflet>>>,
    /// Frequency of the classification.
    frequency: Frequency,
}

impl LeafletClassifier for ManualClassification {
    fn identify_leaflet<'a>(
        &mut self,
        _system: &System,
        _pbc_handler: &impl PBCHandler<'a>,
        molecule_index: usize,
        current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        // get the index of the frame in the assignment
        let assignment_frame = match self.frequency {
            // frame must be zero no matter the thread this is
            Frequency::Once => 0,
            Frequency::Every(n) => current_frame / n.get(),
        };

        Ok(self.assignment
            .as_ref()
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | ManualClassification::identify_leaflet | Assignment has not been set up. {}", PANIC_MESSAGE))
            .get(assignment_frame)
            .ok_or_else(||
                AnalysisError::ManualLeafletError(ManualLeafletClassificationError::FrameNotFound(
                    current_frame, assignment_frame,
                    self.assignment.as_ref().unwrap().len())
                )
            )?
            .get(molecule_index)
            .cloned()
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | ManualClassification::identify_leaflet | Molecule index `{}` not found but this should have been checked before. {}", 
                molecule_index, PANIC_MESSAGE))
            )
    }

    fn n_molecules(&self) -> usize {
        self.assignment
            .as_ref()
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | ManualClassification::n_molecules | Assignment has not been set up. {}", PANIC_MESSAGE)).first()
            .unwrap_or_else(||
                panic!("FATAL GORDER ERROR | ManualClassification::n_molecules | Assignment is empty. {}", PANIC_MESSAGE))
            .len()
    }
}

impl SystemTopology {
    /// Finalize the set up of the manual leaflet classifier.
    pub(super) fn finalize_manual_leaflet_classification(
        &mut self,
        params: &LeafletClassification,
    ) -> Result<(), ManualLeafletClassificationError> {
        let classification = match params {
            LeafletClassification::FromFile(params) => {
                &ManualClassification::map_from_file(params.file())?
            }
            LeafletClassification::FromMap(params) => params.assignment(),
            // do nothing for the other types of classification
            _ => return Ok(()),
        };

        let mut molecule_names = Vec::new();

        handle_moltypes!(self.molecule_types_mut(), x => {
            for molecule in x.iter_mut() {
                let assignment = classification.get(molecule.name()).ok_or_else(|| {
                    ManualLeafletClassificationError::MoleculeTypeNotFound(molecule.name().to_owned())
                })?;

                // perform sanity checks
                // at least one frame must be provided
                if assignment.is_empty() {
                    return Err(ManualLeafletClassificationError::EmptyAssignment(
                        molecule.name().to_owned(),
                    ));
                }

                // number of molecules must be consistent and match the system
                for (i, frame) in assignment.iter().enumerate() {
                    let n_molecules = molecule.n_molecules();
                    if frame.len() != n_molecules {
                        return Err(
                            ManualLeafletClassificationError::InconsistentNumberOfMolecules {
                                expected: n_molecules,
                                molecule: molecule.name().to_owned(),
                                got: frame.len(),
                                frame: i,
                            },
                        );
                    }
                }

                match molecule.leaflet_classification_mut() {
                    Some(MoleculeLeafletClassification::Manual(x, _)) => x.assignment = Some(assignment.clone()),
                    _ => panic!("FATAL GORDER ERROR | SystemTopology::finalize_manual_leaflet_classification | Unexpected MoleculeLeafletClassification. Expected Manual."),
                }

                molecule_names.push(molecule.name().to_owned());
            }
        });

        // check that there is no additional molecule in the leaflet classification structure
        for molecule_name in classification.keys() {
            if !molecule_names.contains(molecule_name) {
                return Err(ManualLeafletClassificationError::UnknownMoleculeType(
                    molecule_name.clone(),
                    molecule_names,
                ));
            }
        }

        Ok(())
    }

    /// Sanity check for manual leaflet classification.
    /// Checks that the number of analyzed frames matches the number of frames specified in the manual leaflet classification.
    /// Doesn't do anything if the leaflet classification is not manual.
    pub(super) fn validate_leaflet_classification(
        &self,
        step: usize,
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        handle_moltypes!(self.molecule_types(), x => {
            for molecule in x.iter() {
                match molecule.leaflet_classification().as_ref() {
                    None => continue,
                    Some(classification) => Self::validate_molecule_leaflet_classification(
                        classification,
                        step,
                        self.total_frames(),
                    )?,
                }
            }
        });

        Ok(())
    }

    fn validate_molecule_leaflet_classification(
        classification: &MoleculeLeafletClassification,
        step: usize,
        total_frames: usize,
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let (frequency, n_frames) = match classification {
            MoleculeLeafletClassification::Manual(x, _) => (
                x.frequency,
                x.assignment.as_ref().expect(PANIC_MESSAGE).len(),
            ),
            MoleculeLeafletClassification::ManualNdx(x, _) => (x.frequency, x.ndx.len()),
            _ => return Ok(()),
        };

        let frequency = match frequency {
            Frequency::Once => Frequency::Once,
            Frequency::Every(n) => {
                Frequency::Every(NonZeroUsize::new(n.get() / step).expect(PANIC_MESSAGE))
            }
        };

        let expected_n_frames = match frequency {
            Frequency::Once => 1,
            Frequency::Every(n) => (total_frames - 1) / n + 1,
        };

        if expected_n_frames != n_frames {
            if matches!(classification, &MoleculeLeafletClassification::Manual(_, _)) {
                Err(Box::from(
                    ManualLeafletClassificationError::UnexpectedNumberOfFrames {
                        assignment_frames: n_frames,
                        analyzed_frames: total_frames,
                        frequency,
                        expected_assignment_frames: expected_n_frames,
                    },
                ))
            } else {
                Err(Box::from(
                    NdxLeafletClassificationError::UnexpectedNumberOfNdxFiles {
                        ndx_files: n_frames,
                        analyzed_frames: total_frames,
                        frequency,
                        expected_ndx_files: expected_n_frames,
                    },
                ))
            }
        } else {
            Ok(())
        }
    }
}

impl ManualClassification {
    /// Read the leaflet classification structure from an input file.
    fn map_from_file(
        file: &str,
    ) -> Result<HashMap<String, Vec<Vec<Leaflet>>>, ManualLeafletClassificationError> {
        let string = read_to_string(file)
            .map_err(|_| ManualLeafletClassificationError::FileNotFound(file.to_owned()))?;
        serde_yaml::from_str(&string)
            .map_err(|e| ManualLeafletClassificationError::CouldNotParse(file.to_owned(), e))
    }
}

/// Leaflet classification method that uses positions of lipid heads and the global membrane center of geometry.
#[derive(Debug, Clone, CopyGetters, Getters, MutGetters, Setters)]
pub(super) struct NdxClassification {
    /// NDX files to read.
    #[getset(get = "pub(super)")]
    ndx: Vec<String>,
    /// Currently loaded groups from an ndx file.
    groups: Option<Groups>,
    /// Frame for which the current `groups` field has been obtained.
    last_assigned_frame: Option<usize>,
    /// Indices of headgroup identifiers (one per molecule).
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    heads: Vec<usize>,
    /// Name of the upper leaflet group.
    #[getset(get = "pub(super)")]
    upper_leaflet: String,
    /// Name of the lower leaflet group.
    #[getset(get = "pub(super)")]
    lower_leaflet: String,
    /// Frequency with which the assignment should be performed.
    /// Note that this is a 'real frequency' (input frequency multiplied by the step_size).
    frequency: Frequency,
}

impl NdxClassification {
    /// Insert a new molecule into the classifier.
    #[inline(always)]
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads
            .push(get_reference_head(molecule, system, group_name!("Heads"))?);
        Ok(())
    }

    /// Read an ndx file and store the groups from it.
    fn read_ndx_file(
        &mut self,
        frame_index: usize,
        n_atoms: usize,
    ) -> Result<(), NdxLeafletClassificationError> {
        // get the index of the frame in the assignment
        let assignment_index = match self.frequency {
            // frame must be zero no matter the thread this is
            Frequency::Once => 0,
            Frequency::Every(n) => frame_index / n.get(),
        };

        // get the ndx file to read
        let ndx_file = self.ndx.get(assignment_index).ok_or_else(|| {
            NdxLeafletClassificationError::FrameNotFound(
                frame_index,
                assignment_index,
                self.ndx.len(),
            )
        })?;

        // read the ndx file
        let (groups, invalid, duplicate) = Groups::from_ndx(ndx_file, n_atoms)
            .map_err(NdxLeafletClassificationError::CouldNotParse)?;

        // handle issues
        for invalid_name in invalid {
            if invalid_name == self.upper_leaflet || invalid_name == self.lower_leaflet {
                return Err(NdxLeafletClassificationError::InvalidName(
                    invalid_name,
                    ndx_file.clone(),
                ));
            }
            // ignore warnings and errors that are irrelevant
        }

        for duplicate_name in duplicate {
            if duplicate_name == self.upper_leaflet || duplicate_name == self.lower_leaflet {
                return Err(NdxLeafletClassificationError::DuplicateName(
                    duplicate_name,
                    ndx_file.clone(),
                ));
            }
            // ignore warnings and errors that are irrelevant
        }

        // check that the required groups exist
        if !groups.exists(self.upper_leaflet()) {
            return Err(NdxLeafletClassificationError::GroupNotFound(
                self.upper_leaflet.clone(),
                "upper-leaflet".to_owned(),
                ndx_file.clone(),
            ));
        }

        if !groups.exists(self.lower_leaflet()) {
            return Err(NdxLeafletClassificationError::GroupNotFound(
                self.lower_leaflet.clone(),
                "lower-leaflet".to_owned(),
                ndx_file.clone(),
            ));
        }

        self.groups = Some(groups);
        self.last_assigned_frame = Some(frame_index);

        Ok(())
    }

    /// Assign molecule to leaflet based on the loaded groups.
    // This could be done much more efficiently using head index -> molecule index hash map,
    // and looping through each NDX group only once
    // but the current trait system requires the assignment to be performed molecule-wise...
    fn assign_molecule(
        &self,
        molecule_index: usize,
    ) -> Result<Leaflet, NdxLeafletClassificationError> {
        let head_index = self
            .heads
            .get(molecule_index)
            .unwrap_or_else(|| panic!("FATAL GORDER ERROR | NdxClassification::assign_molecule | Could not find head identifier for molecule index `{}`.", molecule_index));

        // groups must exist; checked in Self::read_ndx_file
        let groups = self.groups.as_ref().expect(PANIC_MESSAGE);

        // check upper leaflet
        if groups
            .get(self.upper_leaflet())
            .expect(PANIC_MESSAGE)
            .get_atoms()
            .isin(*head_index)
        {
            return Ok(Leaflet::Upper);
        }

        // check lower leaflet
        if groups
            .get(self.lower_leaflet())
            .expect(PANIC_MESSAGE)
            .get_atoms()
            .isin(*head_index)
        {
            return Ok(Leaflet::Lower);
        }

        Err(NdxLeafletClassificationError::AssignmentNotFound(
            molecule_index,
            *head_index,
        ))
    }
}

impl LeafletClassifier for NdxClassification {
    #[inline(always)]
    fn n_molecules(&self) -> usize {
        self.heads.len()
    }

    fn identify_leaflet<'a>(
        &mut self,
        system: &System,
        _pbc_handler: &impl PBCHandler<'a>,
        molecule_index: usize,
        current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        // only read the NDX file if a) it is not read, or b) it is too old
        // we write it defensively, since this is sensitive
        match (&self.groups, self.last_assigned_frame) {
            (None, None) => self
                .read_ndx_file(current_frame, system.get_n_atoms())
                .map_err(AnalysisError::NdxLeafletError)?,
            (Some(_), None) | (None, Some(_)) => {
                panic!("FATAL GORDER ERROR | NdxClassification::identify_leaflet | Inconsistent state of `groups` and `last_ndx_index`.");
            }
            (Some(_), Some(index)) => match index.cmp(&current_frame) {
                Ordering::Less => self
                    .read_ndx_file(current_frame, system.get_n_atoms())
                    .map_err(AnalysisError::NdxLeafletError)?,
                Ordering::Equal => (),
                Ordering::Greater => {
                    panic!("FATAL GORDER ERROR | NdxClassification::identify_leaflet | Last read NDX file is for frame `{}`, but the current frame is `{}`. Went back in time?", index, current_frame);
                }
            },
        };

        // get the assignment for the target molecule
        self.assign_molecule(molecule_index)
            .map_err(AnalysisError::NdxLeafletError)
    }
}

#[derive(Debug, Clone)]
struct Clusters {
    upper: HashSet<usize>,
    lower: HashSet<usize>,
}

/// Leaflet classification method that uses DBSCAN clustering to identify leaflets.
#[derive(Debug, Clone)]
pub(super) struct ClusterClassification {
    /// Indices of headgroup identifiers (one per molecule).
    heads: Vec<usize>,
    /// Radius for the neighbor search.
    radius: f32,
    /// Minimal number of neighbors for a given headgroup to be classified as a core point.
    min_samples: usize,
    /// Clusters identified for the current frame.
    /// TODO: Currently performed independently for each molecule type. This can be optimized dramatically.
    clusters: once_cell::sync::OnceCell<Clusters>,
    /// Clusters identified for the other frames.
    /// Shared among all threads.
    shared_clusters: Arc<Mutex<HashMap<usize, Clusters>>>,
    /// Frequency with which the assignment should be performed.
    frequency: Frequency,
}

impl LeafletClassifier for ClusterClassification {
    fn n_molecules(&self) -> usize {
        self.heads.len()
    }

    fn identify_leaflet<'a>(
        &mut self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        molecule_index: usize,
        current_frame: usize,
    ) -> Result<Leaflet, AnalysisError> {
        let clusters = self.clusters.get_or_try_init(|| {
            let clusters = self.construct_clusters(system, pbc_handler, current_frame)?;

            let mut shared_clusters = self.shared_clusters.lock();
            let previous = shared_clusters.insert(current_frame, clusters.to_owned());

            // defensive check; no other thread should have read the same frame
            assert!(previous.is_none(), "FATAL GORDER ERROR | ClusterClassifier::identify_leaflet | Clusters for frame index `{}` already exist in shared storage, but they should not. {}", 
                current_frame, PANIC_MESSAGE
            );
            Ok(clusters)
        })?;

        if clusters
            .upper
            .contains(self.heads.get(molecule_index).expect(PANIC_MESSAGE))
        {
            Ok(Leaflet::Upper)
        } else {
            Ok(Leaflet::Lower)
        }
    }
}

impl ClusterClassification {
    /// Insert a new molecule into the classifier.
    #[inline(always)]
    fn insert(&mut self, molecule: &Group, system: &System) -> Result<(), TopologyError> {
        self.heads.push(get_reference_head(
            molecule,
            system,
            group_name!("ClusterHeads"),
        )?);
        Ok(())
    }

    /// Identify clusters and classify them.
    fn construct_clusters<'a>(
        &self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        current_frame: usize,
    ) -> Result<Clusters, AnalysisError> {
        let (cluster1, cluster2, min_index, min_index_cluster) =
            self.identify_clusters(system, pbc_handler)?;

        // check that all heads have been assigned into one of the identified clusters (and only once)
        let classified = cluster1.len() + cluster2.len();
        let total = system
            .group_get_n_atoms(group_name!("ClusterHeads"))
            .expect(PANIC_MESSAGE);
        match classified.cmp(&total) {
            Ordering::Less => return Err(AnalysisError::ClusterError(ClusterError::OutlierLipids(
                classified, total,
            ))),
            Ordering::Greater => panic!(
                "FATAL GORDER ERROR | ClusterClassification::construct_clusters | Some lipids were assigned into both leaflet clusters. {}", 
                PANIC_MESSAGE
            ),
            Ordering::Equal => (),
        }

        self.classify_clusters(
            cluster1,
            cluster2,
            min_index,
            min_index_cluster,
            current_frame,
        )
    }

    /// Identify clusters using DBSCAN.
    fn identify_clusters<'a>(
        &self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
    ) -> Result<(HashSet<usize>, HashSet<usize>, usize, i8), AnalysisError> {
        let mut assignments = HashMap::new();
        let mut queue = VecDeque::new();
        let mut curr_cluster = -1i8;

        // loop through all heads
        // NOTE: we are looping in order, not randomly, meaning that
        // the leaflet containing the first lipid will always be "dominant":
        // any lipid that is in range of both leaflets (such as a flip-flopping lipid)
        // will be classified as belonging to the "dominant" leaflet even if it is actually
        // closer to the other leaflet
        // this could be solved by randomizing the order of the heads, but `gorder` must remain
        // fully deterministic, so this is not possible
        for atom in system
            .group_iter(group_name!("ClusterHeads"))
            .expect(PANIC_MESSAGE)
        {
            let index = atom.get_index();
            // skip heads that are already assigned
            if assignments.contains_key(&index) {
                continue;
            }

            let pos = atom
                .get_position()
                .ok_or(AnalysisError::UndefinedPosition(index))?;

            // get the neighbors of this atoms
            let neighbors: Vec<&Atom> = pbc_handler
                .nearby_atoms(system, pos.clone(), self.radius)
                .collect();

            // ignore non-core atoms; these canot start a new cluster
            if neighbors.len() < self.min_samples {
                continue;
            }

            curr_cluster += 1;
            // if the number of clusters is higher than 2, raise an error (there should only be two leaflets)
            if curr_cluster > 1 {
                return Err(AnalysisError::ClusterError(ClusterError::TooManyClusters(
                    index,
                )));
            }

            self.start_new_cluster(
                index,
                &neighbors,
                curr_cluster,
                system,
                pbc_handler,
                &mut assignments,
                &mut queue,
            )?;
        }

        let (cluster1, cluster2, min_index, min_cluster) = self.process_assignments(assignments);
        Ok((cluster1, cluster2, min_index, min_cluster))
    }

    /// Assign atoms into a new cluster.
    fn start_new_cluster<'a>(
        &self,
        head_index: usize,
        neighbors: &[&Atom],
        cluster_id: i8,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        assignments: &mut HashMap<usize, i8>,
        queue: &mut VecDeque<usize>,
    ) -> Result<(), AnalysisError> {
        assignments.insert(head_index, cluster_id);

        // assign neighbors into a queue
        for neighbor in neighbors {
            let index = neighbor.get_index();

            // skip neighbors that are already assigned into a cluster
            if assignments.contains_key(&index) {
                continue;
            }
            assignments.insert(index, cluster_id);
            queue.push_back(index);
        }

        // process the neighbors
        self.process_bfs(system, pbc_handler, assignments, queue, cluster_id)
    }

    /// Assign atoms into a cluster using breadth-first traversal.
    fn process_bfs<'a>(
        &self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        assignments: &mut HashMap<usize, i8>,
        queue: &mut VecDeque<usize>,
        cluster_id: i8,
    ) -> Result<(), AnalysisError> {
        while let Some(index) = queue.pop_front() {
            // safety: index is obtained from a nearby atom
            let atom = unsafe { system.get_atom_unchecked(index) };
            let pos = atom
                .get_position()
                .ok_or(AnalysisError::UndefinedPosition(index))?;

            let neighbors: Vec<&Atom> = pbc_handler
                .nearby_atoms(system, pos.clone(), self.radius)
                .collect();

            if neighbors.len() < self.min_samples {
                continue;
            }

            for neighbor in neighbors {
                let neighbor_index = neighbor.get_index();

                // akip neighbors that are already assigned into a cluster
                if assignments.contains_key(&neighbor_index) {
                    continue;
                }

                assignments.insert(neighbor_index, cluster_id);
                queue.push_back(neighbor_index);
            }
        }
        Ok(())
    }

    /// Process assignments into two clusters.
    fn process_assignments(
        &self,
        assignments: HashMap<usize, i8>,
    ) -> (HashSet<usize>, HashSet<usize>, usize, i8) {
        let mut cluster1 = HashSet::new();
        let mut cluster2 = HashSet::new();
        let mut min_index = usize::MAX;
        let mut min_index_cluster = 0;

        for (index, cluster) in assignments {
            match cluster {
                0 => cluster1.insert(index),
                1 => cluster2.insert(index),
                x => panic!(
                    "FATAL GORDER ERROR | ClusterClassification::process_assignments | Invalid cluster index `{}`. {}",
                    x,
                    PANIC_MESSAGE
                ),
            };

            if index < min_index {
                min_index = index;
                min_index_cluster = cluster;
            }
        }

        (cluster1, cluster2, min_index, min_index_cluster)
    }

    /// Determine which cluster is `upper` and `which` is lower.
    ///
    /// In the first frame, the clusters are classified as follows:
    /// - the more populated leaflet is the `upper` leaflet,
    /// - if both leaflets are equally populated, the `upper` leaflet is the one containing
    ///   a reference atom with the lowest index (typically the first analyzed lipid).
    ///
    /// In the other frames, the clusters are classified by trying to match them with the
    /// clusters from the previous analyzed frame.
    ///
    /// In membranes with lipid flip-flop, the match is heuristic and may in extremely
    /// rare cases be incorrect:
    /// - Matching will succeed, if less than 20% of lipids have changed leaflet between two analyzed frames.
    /// - Matching will fail with an error, if 20-80% of lipids have changed leaflet between two analyzed frames.
    /// - Matching will fail silently and the results will be incorrect if over 80% of lipids have
    ///   changed leaflet between two analyzed frames! This should be basically unphysical, so it's not
    ///   a big concern.
    fn classify_clusters(
        &self,
        cluster1: HashSet<usize>,
        cluster2: HashSet<usize>,
        min_index: usize,
        min_index_cluster: i8,
        frame_index: usize,
    ) -> Result<Clusters, AnalysisError> {
        // only for the first frame
        if frame_index == 0 {
            return match cluster1.len().cmp(&cluster2.len()) {
                // more populated cluster is `upper`
                Ordering::Less => {
                    colog_info!("Clustering leaflet classification: classifying the more populated leaflet as '{}'.", "upper");
                    Ok(Clusters {
                        upper: cluster2,
                        lower: cluster1,
                    })
                }
                Ordering::Greater => {
                    colog_info!("Clustering leaflet classification: classifying the more populated leaflet as '{}'.", "upper");
                    Ok(Clusters {
                        upper: cluster1,
                        lower: cluster2,
                    })
                }
                Ordering::Equal => {
                    // if both clusters are equally populated, cluster containing `min_index` is `upper`
                    colog_info!("Clustering leaflet classification: classifying the leaflet containing lipid with reference atom index '{}' as '{}'.", min_index, "upper");
                    if min_index_cluster == 0 {
                        Ok(Clusters {
                            upper: cluster1,
                            lower: cluster2,
                        })
                    } else {
                        Ok(Clusters {
                            upper: cluster2,
                            lower: cluster1,
                        })
                    }
                }
            };
        }

        // load the clusters from the previous analyzed frame
        let previous_frame = match self.frequency {
            Frequency::Once => panic!("FATAL GORDER ERROR | ClusterClassifier::classify_clusters | Frequency is once, but the current frame is {}, not 0. {}", frame_index, PANIC_MESSAGE),
            Frequency::Every(x) if x.get() > frame_index => panic!("FATAL GORDER ERROR | ClusterClassifier::classify_clusters | Frequency is {}, but the current frame is {}. {}", x.get(), frame_index, PANIC_MESSAGE),
            Frequency::Every(x) => frame_index - x.get(),
        };

        let previous_clusters = self.get_from_shared(previous_frame);

        // match the clusters to the previous identified and classified clusters
        let overlap_cluster1_upper =
            cluster1.intersection(&previous_clusters.upper).count() as f32 / cluster1.len() as f32;
        let overlap_cluster1_lower =
            cluster1.intersection(&previous_clusters.lower).count() as f32 / cluster1.len() as f32;

        if overlap_cluster1_upper < CLUSTER_CLASSIFICATION_LIMIT
            && overlap_cluster1_lower < CLUSTER_CLASSIFICATION_LIMIT
        {
            return Err(AnalysisError::ClusterError(
                ClusterError::CouldNotMatchLeaflets(
                    ((1.0 - CLUSTER_CLASSIFICATION_LIMIT) * 100.0).round() as u8,
                ),
            ));
        }

        if overlap_cluster1_upper < overlap_cluster1_lower {
            Ok(Clusters {
                upper: cluster2,
                lower: cluster1,
            })
        } else {
            Ok(Clusters {
                upper: cluster1,
                lower: cluster2,
            })
        }
    }

    /// Get the clusters for the specified frame. The requested clusters are removed from the shared storage.
    /// This requires waiting for the thread responsible for performing the clustering to finish the analysis of this thread.
    fn get_from_shared(&self, frame: usize) -> Clusters {
        let start_time = Instant::now();
        let mut warning_logged = false;

        // spin-lock: waiting for the requested frame to become available
        loop {
            let mut shared_clusters = self.shared_clusters.lock();
            // take the cluster from the shared storage (it is no longer needed in the storage itself)
            if let Some(c) = shared_clusters.remove(&frame) {
                return c;
            }

            // defensive check for a deadlock
            if start_time.elapsed() > *TIMEOUT {
                if !warning_logged {
                    colog_warn!("DEADLOCKED? Thread has been waiting for shared clustering data (frame '{}') for more than {} seconds.
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
                    panic!("FATAL GORDER ERROR | ClusterClassifier::get_from_shared | Deadlock. Could not get shared clusters for leaflet assignment. Spent more than `{}` seconds inside the spin-lock. {}", 
                    HARD_TIMEOUT_SECONDS, PANIC_MESSAGE)
                }
            }

            // shared data unlock here
        }
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
    fn assign_lipids<'a>(
        &mut self,
        system: &'a System,
        pbc_handler: &'a impl PBCHandler<'a>,
        classifier: &mut impl LeafletClassifier,
        current_frame: usize,
    ) -> Result<(), AnalysisError> {
        self.local = Some(
            (0..classifier.n_molecules())
                .map(|index| classifier.identify_leaflet(system, pbc_handler, index, current_frame))
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
    /// to finish the analysis of this frame.
    fn copy_from(&self, frame: usize) -> Vec<Leaflet> {
        let start_time = Instant::now();
        let mut warning_logged = false;

        // spin-lock: waiting for the requested frame to become available
        loop {
            let shared_data = self.0.lock();
            let assignment = shared_data.get(&frame);
            if let Some(assign) = assignment {
                return assign.clone();
            }

            // defensive check for a deadlock
            if start_time.elapsed() > *TIMEOUT {
                if !warning_logged {
                    colog_warn!("DEADLOCKED? Thread has been waiting for shared leaflet assignment data (frame '{}') for more than {} seconds.
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
    use groan_rs::prelude::AtomIteratorWithBox;

    use crate::analysis::pbc::PBC3D;
    use crate::input::{Axis, DynamicNormal};

    use super::super::common::macros::group_name;
    use super::*;

    #[test]
    fn test_global_leaflet_classification() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        create_group(&mut system, "Heads", "name P").unwrap();
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
            MoleculeLeafletClassification::new(&classifier, &Axis::Z.into(), 1, 1).unwrap();
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
            MoleculeLeafletClassification::new(&classifier, &Axis::Z.into(), 1, 1).unwrap();
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
        let pbc = PBC3D::new(system.get_box().unwrap());
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
            classifier.identify_leaflet(&system, &pbc, 0, 1).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.identify_leaflet(&system, &pbc, 1, 1).unwrap(),
            Leaflet::Lower
        );
    }

    #[test]
    fn test_local_assign_to_leaflet() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        system
            .group_create(group_name!("Membrane"), "@membrane")
            .unwrap();

        let pbc = PBC3D::new(system.get_box().unwrap());
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5),
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

        match &mut classifier {
            MoleculeLeafletClassification::Local(x, _) => {
                x.heads_mut().push(1385);
                x.heads_mut().push(11885);
                x.membrane_center_mut().push(Vector3D::new(0.0, 0.0, 0.0));
                x.membrane_center_mut().push(Vector3D::new(0.0, 0.0, 0.0));
                x.set_membrane_center(&system, &pbc, Dimension::Z).unwrap();
            }
            _ => panic!("Unexpected classification method."),
        }

        assert_eq!(
            classifier.identify_leaflet(&system, &pbc, 0, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.identify_leaflet(&system, &pbc, 1, 2).unwrap(),
            Leaflet::Lower
        );
    }

    #[test]
    fn test_individual_assign_to_leaflet() {
        let system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        let pbc = PBC3D::new(system.get_box().unwrap());
        let mut classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316"),
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
            classifier.identify_leaflet(&system, &pbc, 0, 0).unwrap(),
            Leaflet::Upper
        );
        assert_eq!(
            classifier.identify_leaflet(&system, &pbc, 1, 4).unwrap(),
            Leaflet::Lower
        );
    }

    macro_rules! test_classification_new {
        ($name:ident, $classification:expr, $variant:path, $frequency:pat, $normal:expr, $exp_normal:pat) => {
            #[test]
            fn $name() {
                let single_threaded =
                    MoleculeLeafletClassification::new(&$classification, &$normal, 1, 2).unwrap();
                match single_threaded {
                    $variant(x, y) => {
                        assert!(matches!(x.frequency, $frequency));
                        assert!(matches!(x.membrane_normal, $exp_normal));
                        assert!(y.shared.is_none());
                        assert!(y.local.is_none());
                        assert!(y.local_frame.is_none());
                    }
                    _ => panic!("Incorrect structure constructed."),
                }

                let multi_threaded =
                    MoleculeLeafletClassification::new(&$classification, &$normal, 4, 2).unwrap();
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
        Frequency::Once,
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_local_leaflet_classification_new_once,
        LeafletClassification::local("@membrane", "name P", 2.5).with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Local,
        Frequency::Once,
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_individual_leaflet_classification_new_once,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Individual,
        Frequency::Once,
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_global_leaflet_classification_new_every_five,
        LeafletClassification::global("@membrane", "name P")
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_local_leaflet_classification_new_every_five,
        LeafletClassification::local("@membrane", "name P", 2.5)
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Local,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_individual_leaflet_classification_new_every_five,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Individual,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_global_leaflet_classification_new_every_one,
        LeafletClassification::global("@membrane", "name P")
            .with_frequency(Frequency::every(1).unwrap()),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_local_leaflet_classification_new_every_one,
        LeafletClassification::local("@membrane", "name P", 2.5)
            .with_frequency(Frequency::every(1).unwrap()),
        MoleculeLeafletClassification::Local,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_individual_leaflet_classification_new_every_one,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_frequency(Frequency::every(1).unwrap()),
        MoleculeLeafletClassification::Individual,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_global_leaflet_classification_matching_normals,
        LeafletClassification::global("@membrane", "name P").with_membrane_normal(Axis::Z),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_local_leaflet_classification_matching_normals,
        LeafletClassification::local("@membrane", "name P", 2.5).with_membrane_normal(Axis::Z),
        MoleculeLeafletClassification::Local,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_individual_leaflet_classification_matching_normals,
        LeafletClassification::individual("name P", "name C218 C316")
            .with_membrane_normal(Axis::Z)
            .with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Individual,
        Frequency::Once,
        Axis::Z.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_global_leaflet_classification_unmatching_normals,
        LeafletClassification::global("@membrane", "name P")
            .with_membrane_normal(Axis::X)
            .with_frequency(Frequency::every(5).unwrap()),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_),
        Axis::Z.into(),
        Dimension::X
    );

    test_classification_new!(
        test_local_leaflet_classification_unmatching_normals,
        LeafletClassification::local("@membrane", "name P", 2.5).with_membrane_normal(Axis::Z),
        MoleculeLeafletClassification::Local,
        Frequency::Every(_),
        Axis::X.into(),
        Dimension::Z
    );

    test_classification_new!(
        test_individual_leaflet_classification_unmatching_normals,
        LeafletClassification::individual("name P", "name C218 C316").with_membrane_normal(Axis::Y),
        MoleculeLeafletClassification::Individual,
        Frequency::Every(_),
        Axis::X.into(),
        Dimension::Y
    );

    test_classification_new!(
        test_global_leaflet_classification_dynamic_normals,
        LeafletClassification::global("@membrane", "name P").with_membrane_normal(Axis::Y),
        MoleculeLeafletClassification::Global,
        Frequency::Every(_),
        DynamicNormal::new("name P", 2.0).unwrap().into(),
        Dimension::Y
    );

    test_classification_new!(
        test_local_leaflet_classification_dynamic_normals,
        LeafletClassification::local("@membrane", "name P", 2.5)
            .with_membrane_normal(Axis::Z)
            .with_frequency(Frequency::once()),
        MoleculeLeafletClassification::Local,
        Frequency::Once,
        DynamicNormal::new("name P", 2.0).unwrap().into(),
        Dimension::Z
    );

    test_classification_new!(
        test_individual_leaflet_classification_dynamic_normals,
        LeafletClassification::individual("name P", "name C218 C316").with_membrane_normal(Axis::X),
        MoleculeLeafletClassification::Individual,
        Frequency::Every(_),
        DynamicNormal::new("name P", 2.0).unwrap().into(),
        Dimension::X
    );

    #[test]
    fn test_global_leaflet_classification_new_fail() {
        let params = LeafletClassification::global("@membrane", "name P");

        match MoleculeLeafletClassification::new(
            &params,
            &DynamicNormal::new("name P", 2.0).unwrap().into(),
            1,
            1,
        ) {
            Ok(_) => panic!("Function should have failed."),
            Err(ConfigError::MissingMembraneNormal) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_local_leaflet_classification_new_fail() {
        let params = LeafletClassification::local("@membrane", "name P", 2.5);

        match MoleculeLeafletClassification::new(
            &params,
            &DynamicNormal::new("name P", 2.0).unwrap().into(),
            1,
            1,
        ) {
            Ok(_) => panic!("Function should have failed."),
            Err(ConfigError::MissingMembraneNormal) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_individual_leaflet_classification_new_fail() {
        let params = LeafletClassification::individual("name P", "name C218 C316");

        match MoleculeLeafletClassification::new(
            &params,
            &DynamicNormal::new("name P", 2.0).unwrap().into(),
            1,
            1,
        ) {
            Ok(_) => panic!("Function should have failed."),
            Err(ConfigError::MissingMembraneNormal) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_global_leaflet_classification_new_normals_from_file_fail() {
        let params = LeafletClassification::global("@membrane", "name P");
        match MoleculeLeafletClassification::new(&params, &"orders.yaml".into(), 1, 1) {
            Ok(_) => panic!("Function should have failed."),
            Err(ConfigError::MissingMembraneNormal) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_individual_leaflet_classification_new_normals_from_map_fail() {
        let params = LeafletClassification::global("@membrane", "name P");
        let mut map = HashMap::new();
        map.insert(
            "POPE".to_owned(),
            vec![
                vec![Vector3D::new(1.0, 2.0, 3.0)],
                vec![Vector3D::new(2.0, 3.0, 4.0)],
            ],
        );

        match MoleculeLeafletClassification::new(&params, &map.into(), 1, 1) {
            Ok(_) => panic!("Function should have failed."),
            Err(ConfigError::MissingMembraneNormal) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_should_assign_every() {
        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316"),
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

        for i in 0..50 {
            assert!(classifier.should_assign(i));
        }

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P"),
            &Axis::Z.into(),
            1,
            5,
        )
        .unwrap();

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(!classifier.should_assign(7));
        assert!(classifier.should_assign(10));
        assert!(classifier.should_assign(15));
        assert!(!classifier.should_assign(16));

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5),
            &Axis::Z.into(),
            3,
            5,
        )
        .unwrap();

        assert!(classifier.should_assign(0));
        assert!(!classifier.should_assign(1));
        assert!(!classifier.should_assign(7));
        assert!(classifier.should_assign(10));
        assert!(classifier.should_assign(15));
        assert!(!classifier.should_assign(16));

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::global("@membrane", "name P")
                .with_frequency(Frequency::every(4).unwrap()),
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
            &Axis::Z.into(),
            1,
            7,
        )
        .unwrap();

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
            &Axis::Z.into(),
            6,
            11,
        )
        .unwrap();

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
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

        assert!(classifier.should_assign(0));
        for i in 1..51 {
            assert!(!classifier.should_assign(i));
        }

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::individual("name P", "name C218 C316")
                .with_frequency(Frequency::once()),
            &Axis::Z.into(),
            1,
            7,
        )
        .unwrap();

        assert!(classifier.should_assign(0));
        for i in 1..51 {
            assert!(!classifier.should_assign(i));
        }

        let classifier = MoleculeLeafletClassification::new(
            &LeafletClassification::local("@membrane", "name P", 2.5)
                .with_frequency(Frequency::once()),
            &Axis::Z.into(),
            6,
            11,
        )
        .unwrap();

        assert!(classifier.should_assign(0));
        for i in 1..51 {
            assert!(!classifier.should_assign(i));
        }
    }

    #[test]
    fn test_ndx_fail_missing_assignment() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/cg_leaflets.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![3049, 3061, 6, 25, 37, 73, 3637],
            upper_leaflet: "Upper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::once(),
        };

        classification.read_ndx_file(0, 7000).unwrap();
        match classification.assign_molecule(2) {
            Ok(_) => panic!("Function should have failed."),
            Err(NdxLeafletClassificationError::AssignmentNotFound(mol_id, head_id)) => {
                assert_eq!(mol_id, 2);
                assert_eq!(head_id, 6);
            }
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_ndx_fail_missing_frame() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/cg_leaflets.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![3049, 3061, 25, 37, 73, 3637],
            upper_leaflet: "Upper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::every(5).unwrap(),
        };

        match classification.read_ndx_file(6, 7000) {
            Ok(_) => panic!("Function should have failed."),
            Err(NdxLeafletClassificationError::FrameNotFound(
                frame_index,
                ndx_index,
                total_number,
            )) => {
                assert_eq!(frame_index, 6);
                assert_eq!(ndx_index, 1);
                assert_eq!(total_number, 1);
            }
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }
}

#[cfg(test)]
mod tests_assigned_leaflets {
    use groan_rs::prelude::AtomIteratorWithBox;

    use crate::{analysis::pbc::PBC3D, input::Axis};

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
            &Axis::Z.into(),
            1,
            1,
        )
        .unwrap();

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
        let (system, mut classifier, mut leaflets) = prepare_data(false);
        let pbc = PBC3D::new(system.get_box().unwrap());

        leaflets
            .assign_lipids(&system, &pbc, &mut classifier, 15)
            .unwrap();

        assert_eq!(leaflets.local_frame.unwrap(), 15);
        assert_eq!(leaflets.get_leaflet_from_local(0), Leaflet::Upper);
        assert_eq!(leaflets.get_leaflet_from_local(1), Leaflet::Lower);
        assert!(leaflets.shared.is_none());
    }

    #[test]
    fn test_assign_lipids_with_shared() {
        let (system, mut classifier, mut leaflets) = prepare_data(true);
        let pbc = PBC3D::new(system.get_box().unwrap());

        leaflets
            .assign_lipids(&system, &pbc, &mut classifier, 15)
            .unwrap();

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
mod tests_validation {
    use crate::input::Axis;

    use super::*;

    #[test]
    fn test_validate_not_manual() {
        let params = LeafletClassification::global("@membrane", "name PO4");
        let step = 1;
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            let classification = MoleculeLeafletClassification::new(
                &params,
                &MembraneNormal::Static(Axis::Z),
                n_threads,
                step,
            )
            .unwrap();
            assert!(SystemTopology::validate_molecule_leaflet_classification(
                &classification,
                step,
                101
            )
            .is_ok());
        }
    }

    #[test]
    fn test_validate_once() {
        let params =
            LeafletClassification::from_file("tests/files/inputs/leaflets_files/cg_once.yaml")
                .with_frequency(Frequency::once());
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            for step in [1, 3, 5, 10, 15, 20, 50, 200] {
                let mut classification = MoleculeLeafletClassification::new(
                    &params,
                    &MembraneNormal::Static(Axis::Z),
                    n_threads,
                    step,
                )
                .unwrap();

                if let MoleculeLeafletClassification::Manual(ref mut x, _) = classification {
                    x.assignment = ManualClassification::map_from_file(
                        "tests/files/inputs/leaflets_files/cg_once.yaml",
                    )
                    .unwrap()
                    .get("POPC")
                    .cloned();
                }

                assert!(SystemTopology::validate_molecule_leaflet_classification(
                    &classification,
                    step,
                    101 / step
                )
                .is_ok());
            }
        }
    }

    #[test]
    fn test_validate_once_ndx() {
        let params = LeafletClassification::from_ndx(
            &["index.ndx"],
            "name PO4",
            "UpperLeaflet",
            "LowerLeaflet",
        )
        .with_frequency(Frequency::once());

        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            for step in [1, 3, 5, 10, 15, 20, 50, 200] {
                let classification = MoleculeLeafletClassification::new(
                    &params,
                    &MembraneNormal::Static(Axis::Z),
                    n_threads,
                    step,
                )
                .unwrap();

                assert!(SystemTopology::validate_molecule_leaflet_classification(
                    &classification,
                    step,
                    101 / step
                )
                .is_ok());
            }
        }
    }

    #[test]
    fn test_validate_every() {
        let params =
            LeafletClassification::from_file("tests/files/inputs/leaflets_files/cg_every.yaml");
        let step = 1;
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            let mut classification = MoleculeLeafletClassification::new(
                &params,
                &MembraneNormal::Static(Axis::Z),
                n_threads,
                step,
            )
            .unwrap();

            if let MoleculeLeafletClassification::Manual(ref mut x, _) = classification {
                x.assignment = ManualClassification::map_from_file(
                    "tests/files/inputs/leaflets_files/cg_every.yaml",
                )
                .unwrap()
                .get("POPC")
                .cloned();
            }

            assert!(SystemTopology::validate_molecule_leaflet_classification(
                &classification,
                step,
                101
            )
            .is_ok());
        }
    }

    #[test]
    fn test_validate_every_ndx() {
        let params = LeafletClassification::from_ndx(
            &[
                "index1.ndx",
                "index2.ndx",
                "index3.ndx",
                "index4.ndx",
                "index5.ndx",
                "index6.ndx",
                "index7.ndx",
                "index8.ndx",
            ],
            "name PO4",
            "UpperLeaflet",
            "LowerLeaflet",
        );

        let step = 1;
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            let classification = MoleculeLeafletClassification::new(
                &params,
                &MembraneNormal::Static(Axis::Z),
                n_threads,
                step,
            )
            .unwrap();

            assert!(SystemTopology::validate_molecule_leaflet_classification(
                &classification,
                step,
                8
            )
            .is_ok());
        }
    }

    #[test]
    fn test_validate_every20() {
        let params =
            LeafletClassification::from_file("tests/files/inputs/leaflets_files/cg_every20.yaml")
                .with_frequency(Frequency::every(20).unwrap());
        let step = 1;
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            let mut classification = MoleculeLeafletClassification::new(
                &params,
                &MembraneNormal::Static(Axis::Z),
                n_threads,
                step,
            )
            .unwrap();

            if let MoleculeLeafletClassification::Manual(ref mut x, _) = classification {
                x.assignment = ManualClassification::map_from_file(
                    "tests/files/inputs/leaflets_files/cg_every20.yaml",
                )
                .unwrap()
                .get("POPC")
                .cloned();
            }

            assert!(SystemTopology::validate_molecule_leaflet_classification(
                &classification,
                step,
                101
            )
            .is_ok());
        }
    }

    #[test]
    fn test_validate_every20_ndx() {
        let params = LeafletClassification::from_ndx(
            &[
                "index1.ndx",
                "index2.ndx",
                "index3.ndx",
                "index4.ndx",
                "index5.ndx",
                "index6.ndx",
            ],
            "name PO4",
            "UpperLeaflet",
            "LowerLeaflet",
        )
        .with_frequency(Frequency::every(20).unwrap());

        let step = 1;
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            let classification = MoleculeLeafletClassification::new(
                &params,
                &MembraneNormal::Static(Axis::Z),
                n_threads,
                step,
            )
            .unwrap();

            assert!(SystemTopology::validate_molecule_leaflet_classification(
                &classification,
                step,
                101
            )
            .is_ok());
        }
    }

    #[test]
    fn test_validate_every_step() {
        let params =
            LeafletClassification::from_file("tests/files/inputs/leaflets_files/cg_every20.yaml");
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            for step in [1, 2, 3, 5, 10, 15, 20, 25, 30, 75] {
                let mut classification = MoleculeLeafletClassification::new(
                    &params,
                    &MembraneNormal::Static(Axis::Z),
                    n_threads,
                    step,
                )
                .unwrap();

                if let MoleculeLeafletClassification::Manual(ref mut x, _) = classification {
                    x.assignment = ManualClassification::map_from_file(
                        "tests/files/inputs/leaflets_files/cg_every20.yaml",
                    )
                    .unwrap()
                    .get("POPC")
                    .cloned();
                }

                assert!(SystemTopology::validate_molecule_leaflet_classification(
                    &classification,
                    step,
                    6
                )
                .is_ok());
            }
        }
    }

    #[test]
    fn test_validate_every_step_ndx() {
        let params = LeafletClassification::from_ndx(
            &[
                "index1.ndx",
                "index2.ndx",
                "index3.ndx",
                "index4.ndx",
                "index5.ndx",
                "index6.ndx",
            ],
            "name PO4",
            "UpperLeaflet",
            "LowerLeaflet",
        );

        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            for step in [1, 2, 3, 5, 10, 15, 20, 25, 30, 75] {
                let classification = MoleculeLeafletClassification::new(
                    &params,
                    &MembraneNormal::Static(Axis::Z),
                    n_threads,
                    step,
                )
                .unwrap();

                assert!(SystemTopology::validate_molecule_leaflet_classification(
                    &classification,
                    step,
                    6
                )
                .is_ok());
            }
        }
    }

    #[test]
    fn test_validate_every4() {
        let params =
            LeafletClassification::from_file("tests/files/inputs/leaflets_files/cg_every20.yaml")
                .with_frequency(Frequency::every(4).unwrap());
        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            for step in [1, 2, 3, 5, 10, 15, 20, 25, 30, 75] {
                let mut classification = MoleculeLeafletClassification::new(
                    &params,
                    &MembraneNormal::Static(Axis::Z),
                    n_threads,
                    step,
                )
                .unwrap();

                if let MoleculeLeafletClassification::Manual(ref mut x, _) = classification {
                    x.assignment = ManualClassification::map_from_file(
                        "tests/files/inputs/leaflets_files/cg_every20.yaml",
                    )
                    .unwrap()
                    .get("POPC")
                    .cloned();
                }

                assert!(SystemTopology::validate_molecule_leaflet_classification(
                    &classification,
                    step,
                    21
                )
                .is_ok());
            }
        }
    }

    #[test]
    fn test_validate_every4_ndx() {
        let params = LeafletClassification::from_ndx(
            &[
                "index1.ndx",
                "index2.ndx",
                "index3.ndx",
                "index4.ndx",
                "index5.ndx",
                "index6.ndx",
            ],
            "name PO4",
            "UpperLeaflet",
            "LowerLeaflet",
        )
        .with_frequency(Frequency::every(4).unwrap());

        for n_threads in [1, 2, 3, 4, 8, 64, 128] {
            for step in [1, 2, 3, 5, 10, 15, 20, 25, 30, 75] {
                let classification = MoleculeLeafletClassification::new(
                    &params,
                    &MembraneNormal::Static(Axis::Z),
                    n_threads,
                    step,
                )
                .unwrap();

                assert!(SystemTopology::validate_molecule_leaflet_classification(
                    &classification,
                    step,
                    21
                )
                .is_ok());
            }
        }
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

#[cfg(test)]
mod tests_ndx_reading {
    use super::*;

    #[test]
    fn test_ndx_fail_duplicate_group() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/leaflets_duplicate_main.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![1, 2, 3, 4],
            upper_leaflet: "Upper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::once(),
        };

        match classification.read_ndx_file(0, 10) {
            Ok(_) => panic!("Function should have failed."),
            Err(NdxLeafletClassificationError::DuplicateName(name, ndx)) => {
                assert_eq!(name, "Upper");
                assert_eq!(ndx, "tests/files/ndx/leaflets_duplicate_main.ndx");
            }
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_ndx_fail_invalid_group_name() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/leaflets_invalid_main.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![1, 2, 3, 4],
            upper_leaflet: "U!pper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::once(),
        };

        match classification.read_ndx_file(0, 10) {
            Ok(_) => panic!("Function should have failed."),
            Err(NdxLeafletClassificationError::InvalidName(name, ndx)) => {
                assert_eq!(name, "U!pper");
                assert_eq!(ndx, "tests/files/ndx/leaflets_invalid_main.ndx");
            }
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_ndx_fail_upper_not_found() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/leaflets_missing_upper.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![1, 2, 3, 4],
            upper_leaflet: "Upper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::once(),
        };

        match classification.read_ndx_file(0, 10) {
            Ok(_) => panic!("Function should have failed."),
            Err(NdxLeafletClassificationError::GroupNotFound(name, leaflet, ndx)) => {
                assert_eq!(name, "Upper");
                assert_eq!(leaflet, "upper-leaflet");
                assert_eq!(ndx, "tests/files/ndx/leaflets_missing_upper.ndx");
            }
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_ndx_fail_lower_not_found() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/leaflets_missing_lower.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![1, 2, 3, 4],
            upper_leaflet: "Upper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::once(),
        };

        match classification.read_ndx_file(0, 10) {
            Ok(_) => panic!("Function should have failed."),
            Err(NdxLeafletClassificationError::GroupNotFound(name, leaflet, ndx)) => {
                assert_eq!(name, "Lower");
                assert_eq!(leaflet, "lower-leaflet");
                assert_eq!(ndx, "tests/files/ndx/leaflets_missing_lower.ndx");
            }
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_ndx_pass_empty_leaflet() {
        let mut classification = NdxClassification {
            ndx: vec!["tests/files/ndx/leaflets_only_upper.ndx".to_string()],
            groups: None,
            last_assigned_frame: None,
            heads: vec![1, 2, 3, 4],
            upper_leaflet: "Upper".to_string(),
            lower_leaflet: "Lower".to_string(),
            frequency: Frequency::once(),
        };

        classification.read_ndx_file(0, 10).unwrap();
        let groups = classification.groups.unwrap();
        assert_eq!(groups.get_n_atoms("Upper").unwrap(), 4);
        assert_eq!(groups.get_n_atoms("Lower").unwrap(), 0);
    }
}
