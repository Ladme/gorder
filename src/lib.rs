// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

/// Version of the 'gorder' crate / program.
pub const GORDER_VERSION: &str = env!("CARGO_PKG_VERSION");

mod aaorder;
pub mod analysis;
mod auxiliary;
pub mod axis;
pub mod errors;
pub mod leaflets;
mod molecule;
pub mod ordermap;

pub use analysis::Analysis;
pub use axis::Axis;
