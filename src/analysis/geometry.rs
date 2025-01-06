// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Handles dynamic geometry selection during the analysis run.

use groan_rs::{
    prelude::{Rectangular, Shape, SimBox, Vector3D},
    system::System,
};

use crate::{
    errors::TopologyError,
    input::{
        geometry::{CuboidSelection, GeomReference},
        Geometry,
    },
    PANIC_MESSAGE,
};

use super::common::macros::group_name;

/// Trait implemented by all structures that can be used for geometry selection.
pub(super) trait GeometrySelection {
    /// Create the structure from the provided `Geometry` properties.
    fn new(geometry: &Geometry, simbox: &SimBox) -> Self;

    /// Prepare the system for geometry selection, i.e. construct the required groups.
    fn prepare_system(&self, system: &mut System) -> Result<(), TopologyError>;

    /// Calculate and set the reference position for this frame.
    fn init_reference(&mut self, system: &System);

    /// Is the point inside the geometry?
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool;
}

/// No geometry selection requested.
#[derive(Debug, Clone)]
pub(super) struct NoSelection {}

/// Cuboid geometry selection.
#[derive(Debug, Clone)]
pub(super) struct CuboidAnalysis {
    properties: CuboidSelection,
    shape: Rectangular,
}

impl GeometrySelection for CuboidAnalysis {
    fn new(geometry: &Geometry, simbox: &SimBox) -> Self {
        match geometry {
            Geometry::Cuboid(cuboid) => {
                let reference_point = match cuboid.reference() {
                    GeomReference::Point(x) => x.clone(), // fixed value
                    GeomReference::Selection(_) => Vector3D::default(), // we set any value; reference will be set for each frame inside `init_reference`
                };

                let shape = CuboidAnalysis::construct_shape(&cuboid, reference_point, simbox);

                CuboidAnalysis {
                    properties: cuboid.clone(),
                    shape,
                }
            }
            _ => panic!(
                "FATAL GORDER ERROR | CuboidAnalysis::new | Unexpected geometry type `{:?}`. {}",
                geometry, PANIC_MESSAGE
            ),
        }
    }

    fn prepare_system(&self, system: &mut System) -> Result<(), TopologyError> {
        match self.properties.reference() {
            GeomReference::Point(_) => Ok(()), // nothing to do
            GeomReference::Selection(query) => {
                // construct a group for geometry reference
                super::common::create_group(system, "GeomReference", query)
            }
        }
    }

    fn init_reference(&mut self, system: &System) {
        match self.properties.reference() {
            GeomReference::Point(_) => (), // nothing to do, reference position is fixed
            GeomReference::Selection(_) => {
                // calculate the new position and update the shape
                let reference_point = system
                    .group_get_center(group_name!("GeomReference"))
                    .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | CuboidAnalysis::init_reference | Group specifying geometry reference should exist. {}", PANIC_MESSAGE));

                let shape = CuboidAnalysis::construct_shape(
                    &self.properties,
                    reference_point,
                    // validity of the simulation box is checked at the start of every frame analysis (inside `common::analyze_frame`)
                    system.get_box().unwrap_or_else(|| panic!("FATAL GORDER ERROR | CuboidAnalysis::init_reference | Simulation box is undefined but this should have been handled before. {}", PANIC_MESSAGE)),
                );

                self.shape = shape;
            }
        }
    }

    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool {
        self.shape.inside(point, simbox)
    }
}

impl CuboidAnalysis {
    /// Convert the input cuboid box into `groan_rs`'s native rectangular shape.
    fn construct_shape(
        properties: &CuboidSelection,
        mut reference: Vector3D,
        simbox: &SimBox,
    ) -> Rectangular {
        let x = properties.xdim()[1] - properties.xdim()[0];
        let y = properties.ydim()[1] - properties.ydim()[0];
        let z = properties.zdim()[1] - properties.zdim()[0];

        reference.x += properties.xdim()[0];
        reference.y += properties.ydim()[0];
        reference.z += properties.zdim()[0];

        reference.wrap(simbox);

        Rectangular::new(reference, x, y, z)
    }
}
