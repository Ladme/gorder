// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

/// Version of the 'gorder' crate / program.
pub const GORDER_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Message that should be added to every panic.
pub(crate) const PANIC_MESSAGE: &str =
    "\n\n\n            >>> THIS SHOULD NOT HAVE HAPPENED! PLEASE REPORT THIS ERROR <<<
(open an issue at 'github.com/Ladme/gorder/issues' or write an e-mail to 'ladmeb@gmail.com')\n\n";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Leaflet {
    Upper,
    Lower,
}

impl Display for Leaflet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Leaflet::Upper => write!(f, "upper"),
            Leaflet::Lower => write!(f, "lower"),
        }
    }
}

mod analysis;
pub mod errors;
pub mod input;
mod presentation;

use std::fmt::Display;

pub use input::Analysis;
pub use input::AnalysisType;
pub use input::Axis;
pub use input::LeafletClassification;
pub use input::OrderMap;

pub mod prelude {
    pub use super::Analysis;
    pub use super::AnalysisType;
    pub use super::Axis;
    pub use super::LeafletClassification;
    pub use super::OrderMap;
}
