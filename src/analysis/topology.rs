// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementation of a structure describing the topology of the entire system.

use super::molecule::MoleculeType;
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
}

impl SystemTopology {
    #[inline(always)]
    pub(crate) fn new(
        molecule_types: Vec<MoleculeType>,
        membrane_normal: Dimension,
    ) -> SystemTopology {
        SystemTopology {
            molecule_types,
            membrane_normal,
        }
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
                .zip(rhs.molecule_types.into_iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<MoleculeType>>(),
            membrane_normal: self.membrane_normal,
        }
    }
}
