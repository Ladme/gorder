// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use gorder_core::prelude::Vector3D;
use pyo3::prelude::*;

use gorder_core::input::GeomReference as RsReference;
use gorder_core::input::Geometry as RsGeometry;

use crate::string2axis;
use crate::ConfigError;

macro_rules! try_extract {
    ($obj:expr, $( $t:ty ),*) => {
        $(
            if let Ok(geometry_type) = $obj.extract::<$t>() {
                return Ok(Self(geometry_type.0));
            }
        )*
    };
}

/// Helper structure for geometry selection.
#[derive(Clone)]
pub struct Geometry(pub(crate) RsGeometry);

impl<'source> FromPyObject<'source> for Geometry {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        try_extract!(obj, Cuboid, Cylinder, Sphere);

        Err(ConfigError::new_err(
            "expected an instance of Cuboid, Cylinder, or Sphere",
        ))
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Cuboid(RsGeometry);

#[pymethods]
impl Cuboid {
    /// Construct a cuboid.
    /// Returns an error in case any of the dimensions is invalid.
    ///
    /// ## Parameters
    /// - `xdim` (optional) - span of the cuboid along the x-axis; infinite if not specified
    /// - `ydim` (optional) - span of the cuboid along the y-axis; infinite if not specified
    /// - `zdim` (optional) - span of the cuboid along the z-axis; infinite if not specified
    /// - `reference` (optional)
    ///    - reference point relative to which the position of the cuboid is specified
    ///    - if not specified, the origin of the simulation box ([0, 0, 0]) is used
    #[new]
    #[pyo3(signature = (
        xdim = [f32::NEG_INFINITY, f32::INFINITY], 
        ydim = [f32::NEG_INFINITY, f32::INFINITY], 
        zdim = [f32::NEG_INFINITY, f32::INFINITY],
        reference = [0.0, 0.0, 0.0].into()))]
    pub fn new(
        xdim: [f32; 2],
        ydim: [f32; 2],
        zdim: [f32; 2],
        reference: GeomReference,
    ) -> PyResult<Self> {
        Ok(Self(
            RsGeometry::cuboid(reference.0, xdim, ydim, zdim)
                .map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Cylinder(RsGeometry);

#[pymethods]
impl Cylinder {
    /// Construct a cylinder.
    /// Returns an error in case the radius is negative or the span is invalid.
    ///
    /// ## Parameters
    /// - `radius` - radius of the cylinder
    /// - `orientation` - orientation of the main axis of the cylinder
    /// - `span` (optional) - span of the cylinder along its main axis; infinite if not specified
    /// - `reference` (optional)
    ///    - reference point relative to which the position and size of the cylinder is specified
    ///    - if not specified, the origin of the simulation box ([0, 0, 0]) is used
    #[new]
    #[pyo3(signature = (
        radius,
        orientation, 
        span = [f32::NEG_INFINITY, f32::INFINITY], 
        reference = [0.0, 0.0, 0.0].into()))]
    pub fn new(
        radius: f32,
        orientation: &str,
        span: [f32; 2],
        reference: GeomReference,
    ) -> PyResult<Self> {
        Ok(Self(
            RsGeometry::cylinder(reference.0, radius, span, string2axis(orientation)?)
                .map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Sphere(RsGeometry);

#[pymethods]
impl Sphere {
    /// Construct a sphere.
    /// Returns an error in case the radius is negative.
    ///
    /// ## Parameters
    /// - `radius` - radius of the sphere
    /// - `reference` 
    ///    - center of the sphere
    ///    - if not specified, the origin of the simulation box ([0, 0, 0]) is used
    #[new]
    #[pyo3(signature = (
        radius,
        reference = [0.0, 0.0, 0.0].into()))]
    pub fn new(
        radius: f32,
        reference: GeomReference,
    ) -> PyResult<Self> {
        Ok(Self(
            RsGeometry::sphere(reference.0, radius)
                .map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}

#[derive(Clone)]
pub struct GeomReference(RsReference);

impl<'source> FromPyObject<'source> for GeomReference {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        // try to extract as Vector3D
        if let Ok(pos) = obj.extract::<[f32; 3]>() {
            return Ok(Self(RsReference::Point(Vector3D::new(
                pos[0], pos[1], pos[2],
            ))));
        }

        // try to extract as a string
        if let Ok(s) = obj.extract::<String>() {
            let s_lower = s.to_lowercase();
            if &s_lower == "center" {
                return Ok(Self(RsReference::Center));
            }

            return Ok(Self(RsReference::Selection(s)));
        }

        Err(ConfigError::new_err(
            "invalid type for GeomReference constructor: expected a list or string",
        ))
    }
}

impl From<[f32; 3]> for GeomReference {
    fn from(value: [f32; 3]) -> Self {
        Self(RsReference::Point(Vector3D::new(value[0], value[1], value[2])))
    }
}