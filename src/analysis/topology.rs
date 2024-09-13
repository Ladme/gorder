// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use super::molecule::MoleculeType;
use getset::{CopyGetters, Getters, MutGetters};
use groan_rs::prelude::Dimension;
use std::ops::Add;

/// Structure used in parallel analysis of the trajectory file.
#[derive(Debug, Clone, CopyGetters, Getters, MutGetters)]
pub(super) struct SystemTopology {
    #[getset(get = "pub(super)", get_mut = "pub(super)")]
    molecules: Vec<MoleculeType>,
    #[getset(get_copy = "pub(super)")]
    membrane_normal: Dimension,
}

impl SystemTopology {
    pub(super) fn new(molecules: Vec<MoleculeType>, membrane_normal: Dimension) -> SystemTopology {
        SystemTopology {
            molecules,
            membrane_normal,
        }
    }
}

impl Add<SystemTopology> for SystemTopology {
    type Output = Self;
    fn add(self, rhs: SystemTopology) -> Self::Output {
        Self {
            molecules: self
                .molecules
                .into_iter()
                .zip(rhs.molecules.into_iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<MoleculeType>>(),
            membrane_normal: self.membrane_normal,
        }
    }
}
