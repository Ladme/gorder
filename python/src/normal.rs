// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use gorder_core::input::DynamicNormal as RsDynamic;
use gorder_core::input::MembraneNormal as RsNormal;
use gorder_core::prelude::Vector3D;
use hashbrown::HashMap;
use numpy::ndarray;
use numpy::ndarray::ArrayView3;
use numpy::PyArray3;
use numpy::PyArrayMethods;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::string2axis;
use crate::ConfigError;

/// Structure describing the direction of the membrane normal
/// or properties necessary for its calculation.
#[derive(Clone)]
pub struct MembraneNormal(pub(crate) RsNormal);

impl<'source> FromPyObject<'source> for MembraneNormal {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        // try to extract as DynamicNormal
        if let Ok(dyn_norm) = obj.extract::<DynamicNormal>() {
            return Ok(MembraneNormal(RsNormal::Dynamic(dyn_norm.0)));
        }
        // try to extract as a string
        if let Ok(s) = obj.extract::<String>() {
            let s_lower = s.to_lowercase();
            if s_lower == "x" || s_lower == "y" || s_lower == "z" {
                return Ok(MembraneNormal(RsNormal::Static(string2axis(&s_lower)?)));
            } else {
                return Ok(MembraneNormal(RsNormal::FromFile(s)));
            }
        }
        // try to extract as a dictionary
        if let Ok(map) = extract_map(obj) {
            return Ok(MembraneNormal(RsNormal::FromMap(map)));
        }

        Err(ConfigError::new_err(
            "invalid type for MembraneNormal constructor: expected a str, DynamicNormal, or dict",
        ))
    }
}

/// Request a dynamic local membrane normal calculation.
///
/// Attributes
/// ----------
/// heads : str
///     Selection query specifying reference atoms representing lipid headgroups
///     (typically phosphorus atoms or phosphate beads).
///     There must be exactly one such atom/bead per lipid molecule.
/// radius : float
///     Radius of the sphere used to select nearby lipids for membrane normal estimation.
///     The recommended value is half the membrane thickness.
#[pyclass]
#[derive(Clone)]
pub struct DynamicNormal(RsDynamic);

#[pymethods]
impl DynamicNormal {
    #[new]
    pub fn new(heads: &str, radius: f32) -> PyResult<Self> {
        Ok(Self(
            RsDynamic::new(heads, radius).map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}

/// Converts a three-dimensional numpy array into a Vec<Vec<Vector3D>>.
/// The numpy array must have shape [outer, inner, 3].
fn extract_nested_vector<'py>(py_obj: &Bound<'py, PyAny>) -> PyResult<Vec<Vec<Vector3D>>> {
    // try to downcast the input to a PyArray3 of f32
    let array = py_obj.downcast::<PyArray3<f32>>().map_err(|_| {
        ConfigError::new_err("expected a 3D numpy array for dynamic membrane normals")
    })?;
    let view: ArrayView3<f32> = unsafe { array.as_array() };

    let shape = view.shape();
    if shape.len() != 3 || shape[2] != 3 {
        return Err(ConfigError::new_err(
            "expected a 3D numpy array with shape (n_frames, n_molecules, 3)",
        ));
    }
    let outer = shape[0];
    let inner = shape[1];
    let mut result = Vec::with_capacity(outer);
    for i in 0..outer {
        let mut inner_vec = Vec::with_capacity(inner);
        for j in 0..inner {
            // slice the row [i, j, :]
            let slice = view.slice(ndarray::s![i, j, ..]);
            // create a Vector3D from the slice
            let vec3 = Vector3D::new(slice[0], slice[1], slice[2]);
            inner_vec.push(vec3);
        }
        result.push(inner_vec);
    }
    Ok(result)
}

/// Converts a Python dictionary whose keys are strings and values are 3D numpy arrays
/// into a hashbrown::HashMap<String, Vec<Vec<Vector3D>>>.
fn extract_map<'py>(py_obj: &Bound<'py, PyAny>) -> PyResult<HashMap<String, Vec<Vec<Vector3D>>>> {
    let dict = py_obj.downcast::<PyDict>().map_err(|_| {
        ConfigError::new_err(
            "expected a dictionary using molecule types as keys and 3D numpy arrays with shape (n_frames, n_molecules, 3) as values",
        )
    })?;
    let mut map = HashMap::new();
    for (key, value) in dict.iter() {
        let key_str: String = key.extract()?;
        let nested_vec = extract_nested_vector(&value)?;
        map.insert(key_str, nested_vec);
    }
    Ok(map)
}
