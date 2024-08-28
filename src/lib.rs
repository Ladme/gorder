// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

/// Version of the 'gorder' crate / program.
pub const GORDER_VERSION: &str = env!("CARGO_PKG_VERSION");
/// Message added to every panic.
const PANIC_MESSAGE: &str =
    "\n\n\n            >>> THIS SHOULD NOT HAVE HAPPENED! PLEASE REPORT THIS ERROR <<<
(open an issue at 'github.com/Ladme/gorder/issues' or write an e-mail to 'ladmeb@gmail.com')\n\n";

pub mod analysis;
mod auxiliary;
pub mod axis;
pub mod errors;
pub mod leaflets;
pub mod ordermap;

pub use analysis::Analysis;
pub use axis::Axis;
