// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use pyo3::prelude::*;

use gorder_core::input::EstimateError as RsEstimateError;

use crate::ConfigError;

/// Parameters for estimating the error of the analysis.
///
/// Attributes
/// ----------
/// n_blocks : int
///     Number of blocks to divide the trajectory into for error estimation.
///     Default is 5. It is recommended not to modify this value to artificially lower error estimates.
/// output_convergence : Optional[str]
///     Optional filename for writing convergence data.
///     If error estimation is enabled, convergence analysis will be performed even if this option
///     is not specified; the results just won’t be written to a file.
#[pyclass]
#[derive(Clone)]
pub struct EstimateError(pub(crate) RsEstimateError);

#[pymethods]
impl EstimateError {
    #[new]
    #[pyo3(signature = (n_blocks = 5, output_convergence = None))]
    pub fn new(n_blocks: usize, output_convergence: Option<&str>) -> PyResult<Self> {
        Ok(Self(
            RsEstimateError::new(Some(n_blocks), output_convergence)
                .map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}
