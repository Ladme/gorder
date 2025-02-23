// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use gorder_core::input::{Axis, Plane};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(gorder, AnalysisError, PyException);
create_exception!(gorder, WriteError, PyException);
create_exception!(gorder, ConfigError, PyException);

#[macro_export]
macro_rules! apply_if_some {
    ($builder:expr, $( $var:ident => $method:ident ),+ $(,)?) => {
        $(
            if let Some(x) = $var { $builder.$method(x); }
        )+
    };
}

#[macro_export]
macro_rules! apply_inner_if_some {
    ($builder:expr, $( $var:ident => $method:ident ),+ $(,)?) => {
        $(
            if let Some(x) = $var { $builder.$method(x.0); }
        )+
    };
}

mod analysis;
mod estimate_error;
mod geometry;
mod leaflets;
mod normal;
mod ordermap;
mod results;

/// Convert string to axis or return a python error if the conversion fails.
fn string2axis(string: impl AsRef<str>) -> PyResult<Axis> {
    match string.as_ref().to_lowercase().as_str() {
        "x" => Ok(Axis::X),
        "y" => Ok(Axis::Y),
        "z" => Ok(Axis::Z),
        _ => Err(ConfigError::new_err(format!(
            "'{}' is not a valid axis",
            string.as_ref()
        ))),
    }
}

/// Convert string to plane or return a python error if the conversion failes.
fn string2plane(string: impl AsRef<str>) -> PyResult<Plane> {
    match string.as_ref().to_lowercase().as_str() {
        "xy" => Ok(Plane::XY),
        "xz" => Ok(Plane::XZ),
        "yz" => Ok(Plane::YZ),
        _ => Err(ConfigError::new_err(format!(
            "'{}' is not a valid plane",
            string.as_ref()
        ))),
    }
}

#[pymodule]
fn gorder(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // exceptions
    m.add("ConfigError", m.py().get_type::<ConfigError>())?;
    m.add("AnalysisError", m.py().get_type::<AnalysisError>())?;
    m.add("WriteError", m.py().get_type::<WriteError>())?;

    // global classes
    m.add_class::<analysis::Analysis>()?;
    m.add_class::<leaflets::Frequency>()?;
    m.add_class::<results::AnalysisResults>()?;

    // module: analysis_types
    let analysis_types = PyModule::new(m.py(), "analysis_types")?;
    analysis_types.add_class::<analysis::AAOrder>()?;
    analysis_types.add_class::<analysis::CGOrder>()?;
    m.add_submodule(&analysis_types)?;

    // module: membrane_normal
    let membrane_normal = PyModule::new(m.py(), "membrane_normal")?;
    membrane_normal.add_class::<normal::DynamicNormal>()?;
    m.add_submodule(&membrane_normal)?;

    // module: leaflets
    let leaflets = PyModule::new(m.py(), "leaflets")?;
    leaflets.add_class::<leaflets::GlobalClassification>()?;
    leaflets.add_class::<leaflets::LocalClassification>()?;
    leaflets.add_class::<leaflets::IndividualClassification>()?;
    leaflets.add_class::<leaflets::ManualClassification>()?;
    leaflets.add_class::<leaflets::NdxClassification>()?;
    m.add_submodule(&leaflets)?;

    // module: estimate_error
    let estimate_error = PyModule::new(m.py(), "estimate_error")?;
    estimate_error.add_class::<estimate_error::EstimateError>()?;
    m.add_submodule(&estimate_error)?;

    // module: geometry
    let geometry = PyModule::new(m.py(), "geometry")?;
    geometry.add_class::<geometry::Cuboid>()?;
    geometry.add_class::<geometry::Cylinder>()?;
    geometry.add_class::<geometry::Sphere>()?;
    m.add_submodule(&geometry)?;

    // module: ordermap
    let ordermap = PyModule::new(m.py(), "ordermap")?;
    ordermap.add_class::<ordermap::OrderMap>()?;
    m.add_submodule(&ordermap)?;

    Ok(())
}
