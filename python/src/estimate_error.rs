// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use pyo3::prelude::*;

use gorder_core::input::EstimateError as RsEstimateError;

use crate::ConfigError;

/// Parameters for estimating the error of the analysis.
#[pyclass]
#[derive(Clone)]
pub struct EstimateError(pub(crate) RsEstimateError);

#[pymethods]
impl EstimateError {
    /// Specify parameters for the error estimation.
    ///
    /// ## Parameters
    /// - `n_blocks` - number of blocks to divide the trajectory into for error estimation
    ///              - default value is 5 and you should not tweak it to get lower error estimates
    /// - `output_convergence` - optional name of the file where convergence data will be written
    ///                        - note that if error estimation is requested, convergence analysis
    ///                          will be performed even if this option is not specified;
    ///                          it will just not be written out
    #[new]
    #[pyo3(signature = (n_blocks = 5, output_convergence = None))]
    pub fn new(n_blocks: usize, output_convergence: Option<&str>) -> PyResult<Self> {
        Ok(Self(
            RsEstimateError::new(Some(n_blocks), output_convergence)
                .map_err(|e| ConfigError::new_err(e.to_string()))?,
        ))
    }
}
