// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This module contains structures and methods for specifying parameters of the analysis.

pub mod analysis;
pub mod axis;
pub mod estimate_error;
pub mod frequency;
pub mod geometry;
pub mod leaflets;
pub mod membrane_normal;
pub mod ordermap;

pub use analysis::{Analysis, AnalysisType};
pub use axis::Axis;
pub use estimate_error::EstimateError;
pub use frequency::Frequency;
pub use geometry::{GeomReference, Geometry};
pub use leaflets::{ClusteringMethod, LeafletClassification};
pub use membrane_normal::{DynamicNormal, MembraneNormal};
pub use ordermap::{GridSpan, OrderMap, Plane};
