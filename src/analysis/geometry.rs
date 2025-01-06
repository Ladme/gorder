// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Handles dynamic geometry selection during the analysis run.

use groan_rs::{
    prelude::{Rectangular, Shape, SimBox, Vector3D},
    system::System,
};

use crate::{
    input::{
        geometry::{CuboidSelection, GeomReference},
        Geometry,
    },
    PANIC_MESSAGE,
};

/// Trait implemented by all structures that can be used for geometry selection.
pub(super) trait GeometrySelection {
    /// Create the structure from the provided `Geometry` properties.
    fn new(geometry: &Geometry) -> Self;

    /// Is the point inside the geometry?
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool;

    /// Calculate and set the reference position for this frame.
    fn init_reference(&mut self, system: &System);

    /// Prepare the system for geometry selection, i.e. construct the required groups.
    fn prepare_system(&self, system: &mut System);
}

/// No geometry selection requested.
#[derive(Debug, Clone)]
pub(super) struct NoSelection {}

/// Cuboid geometry selection.
#[derive(Debug, Clone)]
pub(super) struct CuboidAnalysis {
    reference: GeomReference,
    shape: Rectangular,
}

/*impl GeometrySelection for CuboidAnalysis {
    fn new(geometry: &Geometry) -> Self {
        match geometry {
            Geometry::Cuboid(cuboid) => {
                let reference_point = match cuboid.reference() {
                    GeomReference::Point(x) => x.clone(), // fixed value
                    GeomReference::Selection(_) => Vector3D::default(), // we set any value; reference will be set for each frame inside `init_reference`
                };

                CuboidAnalysis {
                    reference: cuboid.reference(),
                    shape: Rectangular::new(reference_point),
                }
            }
            _ => panic!(
                "FATAL GORDER ERROR | CuboidAnalysis::new | Unexpected geometry type `{:?}`. {}",
                geometry, PANIC_MESSAGE
            ),
        }
    }

    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool {
        self.shape.inside(point, simbox)
    }

    fn init_reference(&mut self, system: &System) {
        match self.reference {
            Point(_) => (), // nothing to do; shape has already been set
        }
    }
}
*/
