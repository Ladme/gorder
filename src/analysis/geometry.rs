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

/// Enum encompassing all possible geometry selections.
#[derive(Debug, Clone)]
pub(crate) enum GeometrySelectionType {
    None(NoSelection),
    Cuboid(CuboidAnalysis),
}

impl Default for GeometrySelectionType {
    fn default() -> Self {
        GeometrySelectionType::None(NoSelection::default())
    }
}

impl GeometrySelectionType {
    /// Construct a geometry selection type from the input geometry.
    pub(super) fn from_geometry(geometry: Option<&Geometry>, simbox: &SimBox) -> Self {
        match geometry {
            None => GeometrySelectionType::None(NoSelection::default()),
            Some(Geometry::Cuboid(_)) => GeometrySelectionType::Cuboid(CuboidAnalysis::new(
                geometry.expect(PANIC_MESSAGE),
                simbox,
            )),
            Some(Geometry::Cylinder(_)) => todo!(),
            Some(Geometry::Sphere(_)) => todo!(),
        }
    }

    /// Log basic information about the performed geometry selection.
    pub(super) fn info(&self) {
        match self {
            GeometrySelectionType::None(_) => (),
            GeometrySelectionType::Cuboid(cuboid) => {
                log::info!(
                    "Will only consider bonds located inside a cuboid:
  x-dimension: from {} nm to {} nm
  y-dimension: from {} nm to {} nm
  z-dimension: from {} nm to {} nm
  relative to {}",
                    cuboid.properties.xdim()[0],
                    cuboid.properties.xdim()[1],
                    cuboid.properties.ydim()[0],
                    cuboid.properties.ydim()[1],
                    cuboid.properties.zdim()[0],
                    cuboid.properties.zdim()[1],
                    cuboid.properties.reference(),
                );
            }
        }
    }

    /// Initialize the reading of a new frame (calculate and set new reference position if needed).
    pub(super) fn init_new_frame(&mut self, system: &System) {
        match self {
            GeometrySelectionType::None(_) => (),
            GeometrySelectionType::Cuboid(x) => x.init_reference(system),
        }
    }
}

/// Trait implemented by all structures that can be used for geometry selection.
pub(crate) trait GeometrySelection: Send + Sync {
    /// Create the structure from the provided `Geometry` properties.
    fn new(geometry: &Geometry, simbox: &SimBox) -> Self;

    /// Prepare the system for geometry selection, i.e. construct the required groups.
    fn prepare_system(&self, system: &mut System) -> Result<(), TopologyError>;

    /// Calculate and set the reference position for this frame.
    fn init_reference(&mut self, system: &System);

    /// Is the point inside the geometry?
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool;
}

/// No geometry selection requested. Order parameters will be calculated for all bonds, no matter where they are.
#[derive(Debug, Clone, Default)]
pub(crate) struct NoSelection {}

impl GeometrySelection for NoSelection {
    #[inline(always)]
    fn new(_geometry: &Geometry, _simbox: &SimBox) -> Self {
        Self {}
    }

    #[inline(always)]
    fn prepare_system(&self, _system: &mut System) -> Result<(), TopologyError> {
        Ok(())
    }

    #[inline(always)]
    fn init_reference(&mut self, _system: &System) {
        ()
    }

    #[inline(always)]
    fn inside(&self, _point: &Vector3D, _simbox: &SimBox) -> bool {
        true
    }
}

/// Cuboid geometry selection.
#[derive(Debug, Clone)]
pub(crate) struct CuboidAnalysis {
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

    #[inline(always)]
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
        #[inline(always)]
        fn compute_dimension(dim: [f32; 2], reference_component: &mut f32) -> f32 {
            match dim {
                [f32::NEG_INFINITY, f32::INFINITY] => {
                    *reference_component = 0.0;
                    f32::INFINITY
                }
                [min, max] => {
                    *reference_component += min;
                    max - min
                }
            }
        }

        let x = compute_dimension(properties.xdim(), &mut reference.x);
        let y = compute_dimension(properties.ydim(), &mut reference.y);
        let z = compute_dimension(properties.zdim(), &mut reference.z);

        reference.wrap(simbox);

        Rectangular::new(reference, x, y, z)
    }
}

#[cfg(test)]
mod tests_cuboid {
    use approx::assert_relative_eq;
    use rand::prelude::*;

    use super::*;

    #[test]
    fn test_construct_shape_origin() {
        let cuboid =
            match Geometry::cuboid("@protein", [2.5, 3.1], [-1.5, 3.5], [-1.0, 1.0]).unwrap() {
                Geometry::Cuboid(x) => x,
                _ => panic!("Invalid geometry."),
            };

        let simbox = SimBox::from([10.0, 6.0, 8.0]);

        let shape = CuboidAnalysis::construct_shape(&cuboid, Vector3D::new(0.0, 0.0, 0.0), &simbox);

        let position = shape.get_position();
        assert_relative_eq!(position.x, 2.5);
        assert_relative_eq!(position.y, 4.5);
        assert_relative_eq!(position.z, 7.0);

        assert_relative_eq!(shape.get_x(), 0.6);
        assert_relative_eq!(shape.get_y(), 5.0);
        assert_relative_eq!(shape.get_z(), 2.0);
    }

    #[test]
    fn test_construct_shape_simple() {
        let cuboid =
            match Geometry::cuboid("@protein", [2.5, 3.1], [-1.5, 3.5], [-1.0, 1.0]).unwrap() {
                Geometry::Cuboid(x) => x,
                _ => panic!("Invalid geometry."),
            };

        let simbox = SimBox::from([10.0, 6.0, 8.0]);

        let shape = CuboidAnalysis::construct_shape(&cuboid, Vector3D::new(8.0, 5.5, 2.0), &simbox);

        let position = shape.get_position();
        assert_relative_eq!(position.x, 0.5);
        assert_relative_eq!(position.y, 4.0);
        assert_relative_eq!(position.z, 1.0);

        assert_relative_eq!(shape.get_x(), 0.6);
        assert_relative_eq!(shape.get_y(), 5.0);
        assert_relative_eq!(shape.get_z(), 2.0);
    }

    #[test]
    fn test_construct_shape_infinity() {
        let cuboid = match Geometry::cuboid(
            "@protein",
            [2.5, 3.1],
            [f32::NEG_INFINITY, f32::INFINITY],
            [-1.0, 1.0],
        )
        .unwrap()
        {
            Geometry::Cuboid(x) => x,
            _ => panic!("Invalid geometry."),
        };

        let simbox = SimBox::from([10.0, 6.0, 8.0]);

        let shape =
            CuboidAnalysis::construct_shape(&cuboid, Vector3D::new(15.0, 5.5, 1.0), &simbox);

        let position = shape.get_position();
        assert_relative_eq!(position.x, 7.5);
        assert_relative_eq!(position.y, 0.0);
        assert_relative_eq!(position.z, 0.0);

        assert_relative_eq!(shape.get_x(), 0.6);
        assert_relative_eq!(shape.get_y(), f32::INFINITY);
        assert_relative_eq!(shape.get_z(), 2.0);
    }

    #[test]
    fn test_inside_random() {
        let mut rng = StdRng::seed_from_u64(1288746347198273);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        for i in 0..100 {
            let xmin: f32 = rng.gen_range(0.0..10.0);
            let xmax = rng.gen_range(0.0..10.0);

            let ymin: f32 = rng.gen_range(0.0..10.0);
            let ymax = rng.gen_range(0.0..10.0);

            let zmin: f32 = rng.gen_range(0.0..10.0);
            let zmax = rng.gen_range(0.0..10.0);

            let mut xrange = [xmin, xmax];
            xrange.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if i % 8 == 0 {
                xrange = [f32::NEG_INFINITY, f32::INFINITY];
            }

            let mut yrange = [ymin, ymax];
            yrange.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if i % 10 == 0 {
                yrange = [f32::NEG_INFINITY, f32::INFINITY];
            }

            let mut zrange = [zmin, zmax];
            zrange.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if i % 6 == 0 {
                zrange = [f32::NEG_INFINITY, f32::INFINITY];
            }

            let cuboid = match Geometry::cuboid("@protein", xrange, yrange, zrange).unwrap() {
                Geometry::Cuboid(x) => x,
                _ => panic!("Invalid geometry."),
            };

            let shape = CuboidAnalysis::construct_shape(&cuboid, Vector3D::default(), &simbox);

            for _ in 0..1000 {
                let pos_x = rng.gen_range(0.0..10.0);
                let pos_y = rng.gen_range(0.0..10.0);
                let pos_z = rng.gen_range(0.0..10.0);

                let point = Vector3D::new(pos_x, pos_y, pos_z);

                let is_inside = point.x > cuboid.xdim()[0]
                    && point.x < cuboid.xdim()[1]
                    && point.y > cuboid.ydim()[0]
                    && point.y < cuboid.ydim()[1]
                    && point.z > cuboid.zdim()[0]
                    && point.z < cuboid.zdim()[1];

                assert_eq!(is_inside, shape.inside(&point, &simbox));
            }
        }
    }
}
