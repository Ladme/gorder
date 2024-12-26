// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementation of a structure describing the topology of the entire system.

use crate::input::{Analysis, EstimateError};

use super::molecule::MoleculeType;
use crate::errors::ErrorEstimationError;
use crate::presentation::converter::{MolConvert, ResultsConverter};
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
    #[getset(get = "pub(super)")]
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

    /// Convert the topology into a results structure.
    pub(super) fn convert<O: MolConvert>(self, analysis: Analysis) -> O {
        let converter = ResultsConverter::<O>::new(analysis);
        converter.convert_topology(self)
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
