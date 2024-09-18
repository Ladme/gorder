// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains the implementation of the analysis logic.

use groan_rs::prelude::{Atom, SimBox, Vector3D};

use crate::{errors::AnalysisError, Analysis, AnalysisType};

mod aaorder;
mod auxiliary;
mod leaflets;
pub(crate) mod molecule;
mod ordermap;
pub(crate) mod topology;

impl Analysis {
    /// Perform the analysis and write out the results.
    pub fn run(&self) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        match self.analysis_type() {
            AnalysisType::AAOrder => crate::analysis::aaorder::analyze_atomistic(self),
            AnalysisType::CGOrder => todo!("CG Order calculation is not yet implemented."),
        }
    }
}

/// Calculate instantenous value of order parameter of a bond defined by two atoms.
/// Simulation box must be valid and orthogonal.
pub(super) fn calc_sch(
    atom1: &Atom,
    atom2: &Atom,
    simbox: &SimBox,
    membrane_normal: &Vector3D,
) -> Result<f32, AnalysisError> {
    let pos1 = atom1
        .get_position()
        .ok_or_else(|| AnalysisError::UndefinedPosition(atom1.get_atom_number()))?;

    let pos2 = atom2
        .get_position()
        .ok_or_else(|| AnalysisError::UndefinedPosition(atom2.get_atom_number()))?;

    let vector = pos1.vector_to(pos2, simbox);
    let angle = vector.angle(membrane_normal);

    let cos = angle.cos();
    let sch = 0.5 - (1.5 * cos * cos);

    Ok(sch)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use groan_rs::prelude::Dimension;

    #[test]
    fn test_calc_sch() {
        let atom1 = Atom::new(1, "LYS", 1, "CA").with_position(Vector3D::new(1.7, 2.1, 9.7));
        let atom2 = Atom::new(1, "LYS", 2, "H").with_position(Vector3D::new(1.9, 2.4, 0.8));

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert_relative_eq!(
            calc_sch(&atom1, &atom2, &simbox, &Dimension::Z.into()).unwrap(),
            -0.8544775
        );
    }

    #[test]
    fn test_calc_sch_no_position() {
        let atom1 = Atom::new(1, "LYS", 1, "CA").with_position(Vector3D::new(1.7, 2.1, 9.7));
        let atom2 = Atom::new(1, "LYS", 2, "H");

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        match calc_sch(&atom1, &atom2, &simbox, &Dimension::Z.into()) {
            Err(AnalysisError::UndefinedPosition(x)) => assert_eq!(x, 2),
            Ok(_) => panic!("Function should have failed."),
            Err(e) => panic!("Incorrect error returned. {}", e),
        }
    }
}
