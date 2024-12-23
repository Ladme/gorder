// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementation of a structure describing the topology of the entire system.

use crate::input::EstimateError;

use super::molecule::MoleculeType;
use crate::errors::ErrorEstimationError;
use getset::{CopyGetters, Getters, MutGetters};
use groan_rs::prelude::Dimension;
use std::ops::Add;

/// Structure describing the topology of the system.
#[derive(Debug, Clone, CopyGetters, Getters, MutGetters)]
pub(crate) struct SystemTopology {
    #[getset(get = "pub(crate)", get_mut = "pub(super)")]
    molecule_types: Vec<MoleculeType>,
    #[getset(get_copy = "pub(super)")]
    membrane_normal: Dimension,
    #[getset(get = "pub(crate)")]
    estimate_error: Option<EstimateError>,
}

impl SystemTopology {
    #[inline(always)]
    pub(crate) fn new(
        molecule_types: Vec<MoleculeType>,
        membrane_normal: Dimension,
        estimate_error: Option<EstimateError>,
    ) -> SystemTopology {
        SystemTopology {
            molecule_types,
            membrane_normal,
            estimate_error,
        }
    }

    /// Report basic information about the result of the error estimation.
    pub(super) fn error_info(&self) -> Result<(), ErrorEstimationError> {
        let Some(estimate_error) = &self.estimate_error else {
            return Ok(());
        };
        let block_size = estimate_error.block_size();

        for mol in &self.molecule_types {
            let Some(bond) = mol.order_bonds().bond_types().get(0) else {
                continue;
            };
            let Some(order) = bond.total().timewise() else {
                break; // if timewise data is not calculated for this molecule, it is not calculated for any
            };

            let n_blocks = order.n_blocks(block_size);
            let n_frames = order.n_frames();

            if n_frames < 6 {
                return Err(ErrorEstimationError::NotEnoughData(n_frames));
            }
            if n_frames < 100 {
                log::warn!("Error estimation: you probably do not have enough data for reasonable error estimation ({} frames might be too little).", n_frames);
            }

            if n_blocks < 6 {
                return Err(ErrorEstimationError::NotEnoughBlocks(n_blocks, n_frames / 6));
            } else {
                log::info!(
                    "Error estimation: collected {} blocks, each consisting of {} trajectory frames (total: {} frames).",
                    n_blocks, block_size, n_blocks * block_size
                );
            }

            if n_frames != n_blocks * block_size {
                log::info!(
                    "Error estimation: data from {} frame(s) did not completely fill block #{} and will be excluded from error estimation.",
                    n_frames - n_blocks * block_size,
                    n_blocks + 1
                );
            }

            break;
        }

        Ok(())
    }
}

impl Add<SystemTopology> for SystemTopology {
    type Output = Self;

    #[inline]
    fn add(self, rhs: SystemTopology) -> Self::Output {
        Self {
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
