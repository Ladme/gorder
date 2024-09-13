// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains structures and methods for the construction of maps of order parameters.

use derive_builder::Builder;
use getset::{CopyGetters, Getters};

use crate::errors::GridSpanError;

#[derive(Debug, Clone, Builder, Getters, CopyGetters)]
pub struct OrderMap {
    /// Directory where the output files containing the individual ordermaps should be saved.
    #[builder(setter(into))]
    #[getset(get = "pub")]
    output_directory: String,
    /// Minimal number of samples in a grid tile required to calculate order parameter for it.
    /// If not specified, the default value is 1.
    #[builder(default = "1")]
    #[getset(get_copy = "pub")]
    min_samples: usize,
    /// Span of the grid along the x-axis.
    #[builder(default)]
    #[getset(get_copy = "pub")]
    dim_x: GridSpan,
    /// Span of the grid along the z-axis.
    #[builder(default)]
    #[getset(get_copy = "pub")]
    dim_y: GridSpan,
    /// The size of the grid bin along the x-axis (relative to the simulation box size).
    /// If not specified, the default value is 0.01.
    #[builder(default = "0.01")]
    #[getset(get_copy = "pub")]
    bin_size_x: f32,
    /// The size of the grid bin along the y-axis (relative to the simulation box size).
    /// If not specified, the default value is 0.01.
    #[builder(default = "0.01")]
    #[getset(get_copy = "pub")]
    bin_size_y: f32,
}

impl OrderMap {
    pub fn new() -> OrderMapBuilder {
        OrderMapBuilder::default()
    }
}

/// Specifies the span of an ordermap grid.
#[derive(Debug, Clone, Copy, Default)]
pub enum GridSpan {
    /// Span should be obtained from the input structure file.
    #[default]
    Auto,
    /// Span is directly provided.
    Manual { start: f32, end: f32 },
}

impl GridSpan {
    /// Create a new valid `GridSpan` structure from the provided floats.
    ///
    /// This performs sanity checks, making sure that no value is lower than 0 and
    /// that start is not higher than end.
    pub fn manual(start: f32, end: f32) -> Result<GridSpan, GridSpanError> {
        if start < 0.0 {
            return Err(GridSpanError::Negative(start));
        }

        if end < 0.0 {
            return Err(GridSpanError::Negative(end));
        }

        if start > end {
            return Err(GridSpanError::Invalid(start, end));
        }

        Ok(GridSpan::Manual { start, end })
    }
}

#[cfg(test)]
mod grid_span {
    use super::*;

    #[test]
    fn manual_pass() {
        GridSpan::manual(0.0, 10.0).unwrap();
    }

    #[test]
    fn manual_fail_negative_start() {
        let error = GridSpan::manual(-1.0, 10.0).unwrap_err();
        assert!(matches!(error, GridSpanError::Negative(val) if val == -1.0));
    }

    #[test]
    fn manual_fail_negative_end() {
        let error = GridSpan::manual(0.0, -10.0).unwrap_err();
        assert!(matches!(error, GridSpanError::Negative(val) if val == -10.0));
    }

    #[test]
    fn manual_fail_invalid() {
        let error = GridSpan::manual(10.0, 0.0).unwrap_err();
        assert!(matches!(error, GridSpanError::Invalid(x, y) if x == 10.0 && y == 0.0));
    }
}
