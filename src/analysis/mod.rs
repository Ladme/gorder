// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains the implementation of the analysis logic.

use groan_rs::prelude::{SimBox, Vector3D};

use crate::{Analysis, AnalysisType, Axis};

mod aaorder;
mod auxiliary;
mod leaflets;
pub(crate) mod molecule;
pub(crate) mod ordermap;
pub(crate) mod topology;

impl Analysis {
    /// Perform the analysis and write out the results.
    #[inline(always)]
    pub fn run(mut self) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.init_ordermap(self.membrane_normal());

        self.info();
        match self.analysis_type() {
            AnalysisType::AAOrder {
                heavy_atoms: _,
                hydrogens: _,
            } => crate::analysis::aaorder::analyze_atomistic(&self),
            AnalysisType::CGOrder { atoms: _ } => {
                todo!("CG Order calculation is not yet implemented.")
            }
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

/// Calculate instantenous value of order parameter of a bond defined by two atoms.
/// Simulation box must be valid and orthogonal.
#[inline(always)]
pub(super) fn calc_sch(
    pos1: &Vector3D,
    pos2: &Vector3D,
    simbox: &SimBox,
    membrane_normal: &Vector3D,
) -> f32 {
    let vector = pos1.vector_to(pos2, simbox);
    let angle = vector.angle(membrane_normal);

    let cos = angle.cos();
    let sch = 0.5 - (1.5 * cos * cos);

    sch
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use groan_rs::prelude::Dimension;

    #[test]
    fn test_calc_sch() {
        let pos1 = Vector3D::new(1.7, 2.1, 9.7);
        let pos2 = Vector3D::new(1.9, 2.4, 0.8);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert_relative_eq!(
            calc_sch(&pos1, &pos2, &simbox, &Dimension::Z.into()),
            -0.8544775
        );
    }
}
