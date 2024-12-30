// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementation of a structure describing the topology of the entire system.

use crate::input::{Analysis, EstimateError};
use crate::PANIC_MESSAGE;

use super::molecule::MoleculeType;
use crate::errors::ErrorEstimationError;
use crate::presentation::converter::{MolConvert, ResultsConverter};
use getset::{CopyGetters, Getters, MutGetters};
use groan_rs::prelude::Dimension;
use groan_rs::system::ParallelTrajData;
use indexmap::IndexMap;
use std::ops::Add;

/// Structure describing the topology of the system.
#[derive(Debug, Clone, CopyGetters, Getters, MutGetters)]
pub(crate) struct SystemTopology {
    /// ID of a thread working with this SystemTopology.
    #[getset(get = "pub(super)")]
    thread_id: usize,
    /// Total number of threads used in the analysis.
    n_threads: usize,
    /// Index of the current frame from the start of the iteration.
    #[getset(get_copy = "pub(super)")]
    frame: usize,
    /// Step size of the analysis.
    step_size: usize,
    /// Total number of frames analyzed by this thread.
    #[getset(get_copy = "pub(crate)")]
    total_frames: usize,
    /// All molecule types for which order parameters are calculated.
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    molecule_types: Vec<MoleculeType>,
    /// Normal of the membrane.
    #[getset(get_copy = "pub(super)")]
    membrane_normal: Dimension,
    /// Parameters of error estimation.
    #[getset(get = "pub(super)")]
    estimate_error: Option<EstimateError>,
}

impl SystemTopology {
    #[inline(always)]
    pub(crate) fn new(
        molecule_types: Vec<MoleculeType>,
        membrane_normal: Dimension,
        estimate_error: Option<EstimateError>,
        step_size: usize,
        n_threads: usize,
    ) -> SystemTopology {
        SystemTopology {
            thread_id: 0,
            frame: 0,
            n_threads,
            step_size,
            total_frames: 0,
            molecule_types,
            membrane_normal,
            estimate_error,
        }
    }

    /// Convert the topology into a results structure.
    pub(super) fn convert<O: MolConvert>(self, analysis: Analysis) -> O {
        let converter = ResultsConverter::<O>::new(analysis);
        converter.convert_topology(self)
    }

    /// Initialize reading of a new frame.
    pub(super) fn init_new_frame(&mut self) {
        self.molecule_types
            .iter_mut()
            .for_each(|mol| mol.init_new_frame());
    }

    /// Increase the frame counter by `step_size`.
    pub(super) fn increase_frame_counter(&mut self) {
        self.frame += self.step_size * self.n_threads;
        self.total_frames += 1;
    }

    /// Report basic information about the result of the error estimation.
    pub(super) fn error_info(&self) -> Result<(), ErrorEstimationError> {
        let Some(estimate_error) = self.estimate_error() else {
            return Ok(());
        };
        let n_blocks = estimate_error.n_blocks();

        for mol in &self.molecule_types {
            let Some(bond) = mol.order_bonds().bond_types().get(0) else {
                continue;
            };
            let Some(order) = bond.total().timewise() else {
                break; // if timewise data is not calculated for this molecule, it is not calculated for any
            };

            let block_size = order.block_size(n_blocks);
            let n_frames = order.n_frames();

            if block_size < 1 {
                return Err(ErrorEstimationError::NotEnoughData(n_frames, n_blocks));
            }

            if block_size < 10 {
                log::warn!("Error estimation: you probably do not have enough data for reasonable error estimation ({} frames might be too little).",
                    n_frames);
            }

            log::info!(
                    "Error estimation: collected {} blocks, each consisting of {} trajectory frames (total: {} frames).",
                    n_blocks, block_size, n_blocks * block_size
                );

            if n_frames != n_blocks * block_size {
                log::info!(
                    "Error estimation: data from {} frame(s) could not be distributed into blocks and will be excluded from error estimation.",
                    n_frames - n_blocks * block_size,
                );
            }

            break;
        }

        Ok(())
    }

    /// Log information about leaflet assignment for the first frame.
    /// This function should only be called once.
    #[cold]
    pub(super) fn log_first_frame_leaflet_assignment_info(&self) {
        let mut upper = IndexMap::new();
        let mut lower = IndexMap::new();

        for mol in self.molecule_types() {
            mol.leaflet_classification().as_ref().map(|x| {
                let stats = x.statistics();
                // we do not want `LIPID: 0` in the output
                if stats.0 > 0 {
                    upper.insert(mol.name(), stats.0);
                }

                if stats.1 > 0 {
                    lower.insert(mol.name(), stats.1);
                }
            });
        }

        fn indexmap2string(map: &IndexMap<&String, usize>) -> String {
            map.into_iter()
                .map(|(name, number)| (name, number.to_string()))
                .fold(String::new(), |mut acc, (name, number)| {
                    if !acc.is_empty() {
                        acc.push_str(", ");
                    }
                    acc.push_str(name);
                    acc.push_str(": ");
                    acc.push_str(&number);
                    acc
                })
        }

        if !upper.is_empty() {
            log::info!(
                "Upper leaflet in the first analyzed frame: {}",
                indexmap2string(&upper)
            );
        }

        if !lower.is_empty() {
            log::info!(
                "Lower leaflet in the first analyzed frame: {}",
                indexmap2string(&lower)
            );
        }
    }

    /// Log the total number of frames analyzed by this thread.
    pub(super) fn log_total_analyzed_frames(&self) {
        log::info!(
            "Trajectory reading completed. Analyzed {} trajectory frames.",
            self.total_frames,
        );
    }
}

impl Add<SystemTopology> for SystemTopology {
    type Output = Self;

    #[inline]
    fn add(self, rhs: SystemTopology) -> Self::Output {
        Self {
            thread_id: self.thread_id, // at this point, the thread_id does not matter
            n_threads: self.n_threads,
            frame: self.frame.max(rhs.frame), // get the higher frame, although it does not actually matter since the frame number is wrong at this point
            step_size: self.step_size,
            total_frames: self.total_frames + rhs.total_frames,
            molecule_types: self
                .molecule_types
                .into_iter()
                .zip(rhs.molecule_types)
                .map(|(a, b)| a + b)
                .collect::<Vec<MoleculeType>>(),
            membrane_normal: self.membrane_normal,
            estimate_error: self.estimate_error,
        }
    }
}

impl ParallelTrajData for SystemTopology {
    fn reduce(mut data: Vec<Self>) -> Self {
        // sort the data by `thread_id` - `groan_rs` does not guarantee the order
        // of data returned from the individual threads; although there is probably no reason
        // for `groan` to return a different order than a simple sequential, we should make sure
        data.sort_by(|a, b| a.thread_id.cmp(&b.thread_id));

        let mut iter = data.into_iter();
        let first = iter.next().unwrap_or_else(|| {
            panic!(
                "FATAL GORDER ERROR | SystemTopology::reduce | Vector should not be empty. {}",
                PANIC_MESSAGE
            )
        });

        iter.fold(first, |acc, top| acc + top)
    }

    fn initialize(&mut self, thread_id: usize) {
        self.thread_id = thread_id;
        self.frame = thread_id;
    }
}
