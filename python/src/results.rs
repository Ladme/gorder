// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use gorder_core::prelude::AnalysisResults as RsResults;
use pyo3::prelude::*;

use crate::WriteError;

/// Results of the analysis.
#[pyclass]
pub struct AnalysisResults(pub(crate) RsResults);

#[pymethods]
impl AnalysisResults {
    /// Write the results of the analysis into output files.
    pub fn write(&self) -> PyResult<()> {
        if let Err(e) = self.0.write() {
            Err(WriteError::new_err(e.to_string()))
        } else {
            Ok(())
        }
    }
}
