// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This module contains the implementation of the analysis logic.

use groan_rs::prelude::Vector3D;

use crate::input::{Analysis, AnalysisType, Axis};
use crate::presentation::AnalysisResults;

mod aaorder;
mod cgorder;
mod common;
pub(crate) mod geometry;
mod leaflets;
pub(crate) mod molecule;
pub(crate) mod order;
pub(crate) mod ordermap;
mod pbc;
mod structure;
pub(crate) mod timewise;
pub(crate) mod topology;

impl Analysis {
    /// Perform the analysis.
    pub fn run(mut self) -> Result<AnalysisResults, Box<dyn std::error::Error + Send + Sync>> {
        self.init_ordermap(self.membrane_normal());
        self.info();

        match self.analysis_type() {
            AnalysisType::AAOrder {
                heavy_atoms: _,
                hydrogens: _,
            } => aaorder::analyze_atomistic(self),
            AnalysisType::CGOrder { beads: _ } => cgorder::analyze_coarse_grained(self),
        }
    }

    /// Finish the ordermap plane initialization.
    fn init_ordermap(&mut self, membrane_normal: Axis) {
        if let Some(map) = self.map_mut().as_mut() {
            if map.plane().is_none() {
                map.set_plane(Some(membrane_normal.perpendicular()));
            }
        }
    }
}

/// Calculate instantenous value of order parameter of a bond defined by a vector going from atom1 to atom2.
/// Simulation box must be valid and orthogonal.
#[inline(always)]
pub(super) fn calc_sch(vector: &Vector3D, membrane_normal: &Vector3D) -> f32 {
    let angle = vector.angle(membrane_normal);

    let cos = angle.cos();

    (1.5 * cos * cos) - 0.5
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use groan_rs::prelude::{Dimension, SimBox};

    #[test]
    fn test_calc_sch() {
        let pos1 = Vector3D::new(1.7, 2.1, 9.7);
        let pos2 = Vector3D::new(1.9, 2.4, 0.8);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert_relative_eq!(
            calc_sch(&(pos1.vector_to(&pos2, &simbox)), &Dimension::Z.into()),
            0.8544775
        );
    }
}
