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

/// Calculate order parameters inside a cuboid.
///
/// Returns an error if any of the dimensions are invalid.
///
/// Attributes
/// ----------
/// xdim : Optional[List[float]]
///     Span of the cuboid along the x-axis. Defaults to infinite if not specified.
/// ydim : Optional[List[float]]
///     Span of the cuboid along the y-axis. Defaults to infinite if not specified.
/// zdim : Optional[List[float]]
///     Span of the cuboid along the z-axis. Defaults to infinite if not specified.
/// reference : Optional[List[float]]
///     Reference point relative to which the position of the cuboid is specified.
///     Defaults to the origin of the simulation box ([0.0, 0.0, 0.0]) if not specified.
#[pyclass]
#[derive(Clone)]
pub struct Cuboid(RsGeometry);

#[pymethods]
impl Cuboid {
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

/// Calculate order parameters inside a cylinder.
///
/// Returns an error if the radius is negative or the span is invalid.
///
/// Attributes
/// ----------
/// radius : float
///     Radius of the cylinder.
/// orientation : str
///     Orientation of the main axis of the cylinder.
/// span : Optional[List[float]]
///     Span of the cylinder along its main axis. Defaults to infinite if not specified.
/// reference : Optional[List[float]]
///     Reference point relative to which the position and size of the cylinder are specified.
///     Defaults to the origin of the simulation box ([0.0, 0.0, 0.0]) if not specified.
#[pyclass]
#[derive(Clone)]
pub struct Cylinder(RsGeometry);

#[pymethods]
impl Cylinder {
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

/// Calculate order parameters inside a sphere.
///
/// Returns an error if the radius is negative.
///
/// Attributes
/// ----------
/// radius : float
///     Radius of the sphere.
/// reference : Optional[List[float]]
///     Center of the sphere.
///     Defaults to the origin of the simulation box ([0.0, 0.0, 0.0]) if not specified.
#[pyclass]
#[derive(Clone)]
pub struct Sphere(RsGeometry);

#[pymethods]
impl Sphere {
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