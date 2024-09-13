// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for specifying settings of the analysis.

pub mod analysis;
pub mod axis;
pub mod leaflets;
pub mod ordermap;

pub use analysis::{Analysis, AnalysisType};
pub use axis::Axis;
pub use leaflets::LeafletClassification;
pub use ordermap::{GridSpan, OrderMap};
