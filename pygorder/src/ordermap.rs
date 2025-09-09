// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use gorder_core::input::ordermap::GridSpan as RsSpan;
use gorder_core::input::ordermap::OrderMap as RsMap;
use gorder_core::input::ordermap::OrderMapBuilder as RsMapBuilder;
use pyo3::prelude::*;

use crate::string2plane;
use crate::ConfigError;

/// Parameters for generating order maps.
///
/// Parameters
/// ----------
/// output_directory : Optional[str]
///     Directory where output files containing individual order maps will be saved.
///
/// min_samples : Optional[int], default=1
///     Minimum number of samples required in a grid tile to calculate the order parameter.
///
/// dim : Optional[List[Union[str, List[float]]]]
///     Span of the grid along the axes:
///     - First span corresponds to the x-axis (for xy or xz plane) or z-axis (for yz plane).
///     - Second span corresponds to the y-axis (for xy or yz plane) or z-axis (for xz plane).
///     If not specified, the span is derived from the simulation box size of the input structure.
///
/// bin_size : Optional[List[float]], default=[0.1, 0.1]
///     Size of the grid bin along the axes:
///     - First bin dimension corresponds to the x-axis (for xy or xz plane) or z-axis (for yz plane).
///     - Second bin dimension corresponds to the y-axis (for xy or yz plane) or z-axis (for xz plane).
///
/// plane : Optional[str]
///     Plane in which the order maps are constructed. Allowed values: `xy`, `xz`, `yz`.
///     If not specified, the plane is assumed to be perpendicular to the membrane normal.
///
/// Raises
/// ------
/// ConfigError
///     If `min_samples` <= 0, `bin_size` <= 0, any `dim` span is invalid (first value <= second),
///     or if `plane` is not one of the allowed values (`xy`, `xz`, `yz`).
#[pyclass]
#[derive(Clone)]
pub struct OrderMap(pub(crate) RsMap);

#[pymethods]
impl OrderMap {
    #[pyo3(signature = (
        output_directory = None,
        min_samples = 1,
        dim = [GridSpan::default(), GridSpan::default()],
        bin_size = [0.1, 0.1],
        plane = None)
    )]
    #[new]
    pub fn new(
        output_directory: Option<&str>,
        min_samples: usize,
        dim: [GridSpan; 2],
        bin_size: [f32; 2],
        plane: Option<&str>,
    ) -> PyResult<Self> {
        let mut builder: RsMapBuilder = RsMap::builder();

        apply_if_some!(builder, output_directory => output_directory);
        builder
            .min_samples(min_samples)
            .bin_size(bin_size)
            .dim([dim[0].0, dim[1].0]);

        if let Some(plane) = plane {
            builder.plane(string2plane(plane)?);
        }

        let inner = builder
            .build()
            .map_err(|e| ConfigError::new_err(e.to_string()))?;
        Ok(Self(inner))
    }
}

pub struct GridSpan(RsSpan);

impl<'source> FromPyObject<'source> for GridSpan {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        // try to extract as [f32; 2]
        if let Ok(manual) = obj.extract::<[f32; 2]>() {
            return Ok(Self(
                RsSpan::manual(manual[0], manual[1])
                    .map_err(|e| ConfigError::new_err(e.to_string()))?,
            ));
        }

        if let Ok(string) = obj.extract::<String>() {
            if string == "auto" {
                return Ok(Self(RsSpan::Auto));
            }
        }

        Err(ConfigError::new_err(
            "invalid type for GridSpan constructor: expected [float, float] or \"auto\"",
        ))
    }
}

impl Default for GridSpan {
    fn default() -> Self {
        Self(RsSpan::Auto)
    }
}
