// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use gorder_core::input::Frequency as RsFreq;
use gorder_core::input::LeafletClassification as RsLeafletClassification;
use gorder_core::Leaflet;
use hashbrown::HashMap;
use numpy::ndarray::ArrayView2;
use numpy::PyArray2;
use numpy::PyArrayMethods;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::string2axis;
use crate::ConfigError;

macro_rules! try_extract {
    ($obj:expr, $( $t:ty ),*) => {
        $(
            if let Ok(classification_type) = $obj.extract::<$t>() {
                return Ok(Self(classification_type.0));
            }
        )*
    };
}

/// Helper structure for LeafletClassification.
#[derive(Clone)]
pub struct LeafletClassification(pub(crate) RsLeafletClassification);

impl<'source> FromPyObject<'source> for LeafletClassification {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        try_extract!(
            obj,
            GlobalClassification,
            LocalClassification,
            IndividualClassification,
            ManualClassification,
            NdxClassification
        );

        Err(ConfigError::new_err(
            "expected an instance of GlobalClassification, LocalClassification, IndividualClassification, ManualClassification, or NdxClassification",
        ))
    }
}

/// Frequency of some action being performed.
#[pyclass]
#[derive(Clone)]
pub struct Frequency(RsFreq);

#[pymethods]
impl Frequency {
    /// Perform the action once.
    #[staticmethod]
    pub fn once() -> Self {
        Frequency(RsFreq::once())
    }

    /// Perform the action every N frames.
    /// Returns an error if `n_frames` is 0.
    #[staticmethod]
    pub fn every(n_frames: usize) -> PyResult<Self> {
        Ok(Frequency(
            RsFreq::every(n_frames).map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}

/// Classify lipids based on the global membrane center of geometry.
#[pyclass]
#[derive(Clone)]
pub struct GlobalClassification(RsLeafletClassification);

#[pymethods]
impl GlobalClassification {
    /// Classify lipids based on the global membrane center of geometry.
    /// Generally reliable and fast. The best option when working with disrupted membranes.
    ///
    /// ## Parameters
    /// - `membrane` - selection of all lipids forming the membrane
    /// - `heads` - reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule.
    /// - `frequency` (optional) - frequency of classification; defaults to every frame
    /// - `membrane_normal` (optional) - membrane normal for leaflet classification; defaults to
    ///    membrane normal specified at the global level
    #[new]
    #[pyo3(signature = (membrane, heads, frequency = None, membrane_normal = None))]
    pub fn new(
        membrane: &str,
        heads: &str,
        frequency: Option<Frequency>,
        membrane_normal: Option<&str>,
    ) -> PyResult<Self> {
        let classification = add_normal(
            add_freq(RsLeafletClassification::global(membrane, heads), frequency)?,
            membrane_normal,
        )?;

        Ok(Self(classification))
    }
}

/// Classify lipids based on the local membrane center of geometry.
#[pyclass]
#[derive(Clone)]
pub struct LocalClassification(RsLeafletClassification);

#[pymethods]
impl LocalClassification {
    /// Classify lipids based on the local membrane center of geometry.
    /// Useful for curved membranes, very slow.
    ///
    /// ## Parameters
    /// - `membrane` - selection of all lipids forming the membrane
    /// - `heads` - reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule.
    /// - `radius` - radius of a cylinder for the calculation of local membrane center of geometry (in nm)
    /// - `frequency` (optional) - frequency of classification; defaults to every frame
    /// - `membrane_normal` (optional) - membrane normal for leaflet classification; defaults to
    ///    membrane normal specified at the global level
    #[new]
    #[pyo3(signature = (membrane, heads, radius, frequency = None, membrane_normal = None))]
    pub fn new(
        membrane: &str,
        heads: &str,
        radius: f32,
        frequency: Option<Frequency>,
        membrane_normal: Option<&str>,
    ) -> PyResult<Self> {
        if radius <= 0.0 {
            return Err(ConfigError::new_err(format!(
                "radius must be greater than 0, not `{}`.",
                radius
            )));
        }

        let classification = add_normal(
            add_freq(
                RsLeafletClassification::local(membrane, heads, radius),
                frequency,
            )?,
            membrane_normal,
        )?;

        Ok(Self(classification))
    }
}

/// Classify lipids based on the orientation of acyl chains.
#[pyclass]
#[derive(Clone)]
pub struct IndividualClassification(RsLeafletClassification);

#[pymethods]
impl IndividualClassification {
    /// Classify lipids based on the orientation of acyl chains.
    /// Less reliable but fast.
    ///
    /// ## Parameters
    /// - `heads`: reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule.
    /// - `methyls`: reference atoms identifying methyl groups of lipid tails, i.e., the ends of lipid tails;
    ///    there should be only one such atom/bead per one acyl chain in the molecule (e.g., two for lipids with two acyl chains).
    /// - `frequency` (optional) - frequency of classification; defaults to every frame
    /// - `membrane_normal` (optional) - membrane normal for leaflet classification; defaults to
    ///    membrane normal specified at the global level
    #[new]
    #[pyo3(signature = (heads, methyls, frequency = None, membrane_normal = None))]
    pub fn new(
        heads: &str,
        methyls: &str,
        frequency: Option<Frequency>,
        membrane_normal: Option<&str>,
    ) -> PyResult<Self> {
        let classification = add_normal(
            add_freq(
                RsLeafletClassification::individual(heads, methyls),
                frequency,
            )?,
            membrane_normal,
        )?;

        Ok(Self(classification))
    }
}

#[pyclass]
#[derive(Clone)]
pub struct ManualClassification(RsLeafletClassification);

#[pymethods]
impl ManualClassification {
    /// Read leaflet assignment from an external yaml file or from a dictionary.
    ///
    /// ## Parameters
    /// - `input`: path to the input yaml file containing the leaflet assignment
    ///    or a dictionary containing the leaflet assignment.
    /// - `frequency` (optional) - frequency of classification; defaults to every frame
    #[new]
    #[pyo3(signature = (input, frequency = None))]
    pub fn new<'source>(
        input: &Bound<'source, PyAny>,
        frequency: Option<Frequency>,
    ) -> PyResult<Self> {
        let classification = if let Ok(file) = input.extract::<String>() {
            RsLeafletClassification::from_file(&file)
        } else if let Ok(map) = extract_map(input) {
            let converted_map = convert_leaflet_map(map)?;
            RsLeafletClassification::from_map(converted_map)
        } else {
            return Err(ConfigError::new_err(
                "invalid type for ManualClassification input: expected str or dict",
            ));
        };

        Ok(Self(add_freq(classification, frequency)?))
    }
}

/// Get leaflet assignment using NDX file(s).
#[pyclass]
#[derive(Clone)]
pub struct NdxClassification(RsLeafletClassification);

#[pymethods]
impl NdxClassification {
    /// Get leaflet assignment using NDX file(s).
    ///
    /// ## Parameters
    /// - `ndx`: a list of NDX files to read
    /// - `heads`: GSL query specifying the atoms to use as head identifiers of molecules
    /// - `upper_leaflet`: name of the group in the NDX file(s) specifying the atoms of the upper leaflet
    /// - `lower_leaflet`: name of the group in the NDX file(s) specifying the atoms of the lower leaflet
    /// - `frequency` (optional) - frequency of classification; defaults to every frame
    ///
    /// ## Notes
    /// - No glob expansion is performed for the NDX files.
    #[new]
    #[pyo3(signature = (ndx, heads, upper_leaflet, lower_leaflet, frequency = None))]
    pub fn new(
        ndx: Vec<String>,
        heads: &str,
        upper_leaflet: &str,
        lower_leaflet: &str,
        frequency: Option<Frequency>,
    ) -> PyResult<Self> {
        let classification = add_freq(
            RsLeafletClassification::from_ndx(
                &ndx.iter().map(String::as_str).collect::<Vec<&str>>(),
                heads,
                upper_leaflet,
                lower_leaflet,
            ),
            frequency,
        )?;

        Ok(Self(classification))
    }
}

/// Attempt to add frequency to leaflet classification.
fn add_freq(
    mut classification: RsLeafletClassification,
    frequency: Option<Frequency>,
) -> PyResult<RsLeafletClassification> {
    if let Some(freq) = frequency {
        classification = classification.with_frequency(freq.0);
    }

    Ok(classification)
}

/// Attempt to add membrane normal to leaflet classification.
fn add_normal(
    mut classification: RsLeafletClassification,
    normal: Option<&str>,
) -> PyResult<RsLeafletClassification> {
    if let Some(normal) = normal {
        classification = classification.with_membrane_normal(string2axis(normal)?);
    }

    Ok(classification)
}

/// Converts a Python dictionary whose keys are strings and values are 2D numpy arrays
/// into a hashbrown::HashMap<String, Vec<Vec<u8>>>.
fn extract_map<'py>(py_obj: &Bound<'py, PyAny>) -> PyResult<HashMap<String, Vec<Vec<u8>>>> {
    let dict = py_obj.downcast::<PyDict>().map_err(|_| {
        ConfigError::new_err(
            "expected a dictionary using molecule types as keys and 2D numpy arrays with shape (n_frames, n_molecules) as values",
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

/// Converts a two-dimensional numpy array into a Vec<Vec<u8>>.
fn extract_nested_vector<'py>(py_obj: &Bound<'py, PyAny>) -> PyResult<Vec<Vec<u8>>> {
    let array = py_obj
        .downcast::<PyArray2<u8>>()
        .map_err(|_| ConfigError::new_err("expected a 2D numpy array for leaflet assignment"))?;

    let array_view: ArrayView2<u8> = unsafe { array.as_array() };

    let vec: Vec<Vec<u8>> = array_view.outer_iter().map(|row| row.to_vec()).collect();

    Ok(vec)
}

/// Convert `HashMap<String, Vec<Vec<u8>>>` to `HashMap<String, Vec<Vec<Leaflet>>>`.
/// 1 => Leaflet::Upper.
/// 0 => Leaflet::Lower.
fn convert_leaflet_map(
    input: HashMap<String, Vec<Vec<u8>>>,
) -> PyResult<HashMap<String, Vec<Vec<Leaflet>>>> {
    input
        .into_iter()
        .map(|(key, matrix)| {
            let converted_matrix: PyResult<Vec<Vec<Leaflet>>> = matrix
                .into_iter()
                .map(|row| row.into_iter().map(|number| match number {
                    1 => Ok(Leaflet::Upper),
                    0 => Ok(Leaflet::Lower),
                    x => Err(ConfigError::new_err(
                        format!("'{}' is not a valid leaflet identifier (use 0 for lower leaflet, and 1 for upper leaflet)", x)))
                }).collect())
                .collect();

            converted_matrix.map(|m| (key, m))
        })
        .collect()
}
