// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for specifying parameters of the analysis.

pub mod analysis;
pub mod axis;
pub mod estimate_error;
pub mod leaflets;
pub mod ordermap;

pub use analysis::{Analysis, AnalysisType};
pub use axis::Axis;
pub use estimate_error::EstimateError;
pub use leaflets::LeafletClassification;
pub use ordermap::{GridSpan, OrderMap, Plane};
