// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! # gorder: Everything you will ever need for lipid order calculations
//!
//! Crate for calculating lipid order parameters from Gromacs simulations.
//! `gorder` can analyze coarse-grained as well as atomistic lipid order parameters.
//!
//! **See the [gorder manual](https://ladme.github.io/gorder-manual/) for full information about the capabilities of the crate.**
//!
//! ## Usage
//!
//! Run:
//!
//! ```bash
//! $ cargo add gorder
//! ```
//!
//! Import the crate in your Rust code:
//!
//! ```rust
//! use gorder::prelude::*;
//! ```
//!
//! `gorder` is also available as a command line tool. You can install it using:
//! ```bash
//! $ cargo install gorder
//! ```
//!
//! ## Examples
//!
//! Basic analysis of atomistic lipid order parameters.
//! ```no_run
//! use gorder::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
//!     // construct the analysis
//!     let analysis = Analysis::new()
//!             .structure("system.tpr")                   // structure file
//!             .trajectory("md.xtc")                      // trajectory file to analyze
//!             .output("order.yaml")                      // output yaml file
//!             .analysis_type(AnalysisType::aaorder(      // analysis to perform
//!                "@membrane and element name carbon",    // selection of heavy atoms
//!                "@membrane and element name hydrogen",  // selection of hydrogens
//!             ))
//!             .build()?;                                 // constructing the analysis
//!
//!     // activate colog if you want logging (requires the `colog` crate)
//!     colog::init();
//!
//!     // run the analysis and write the output
//!     analysis.run()?;
//!
//!     Ok(())
//! }
//! ```
//!
//! `gorder` will automatically recognize molecules and relevant bonds and calculate order parameters for them.
//!
//! ***
//!
//! Basic analysis of coarse-grained lipid order parameters.
//!
//! ```no_run
//! use gorder::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
//!     // construct the analysis
//!     let analysis = Analysis::new()
//!             .structure("system.tpr")                   // structure file
//!             .trajectory("md.xtc")                      // trajectory file to analyze
//!             .output("order.yaml")                      // output yaml file
//!             .analysis_type(AnalysisType::cgorder(      // analysis to perform
//!                "@membrane",                            // selection of beads
//!             ))
//!             .build()?;                                 // constructing the analysis
//!
//!     // activate colog if you want logging (requires the `colog` crate)
//!     colog::init();
//!
//!     // run the analysis and write the output
//!     analysis.run()?;
//!
//!     Ok(())
//! }
//! ```
//!  
//! ***
//!
//! The `Analysis` structure has lot of other, optional fields.
//!
//! ```no_run
//! use gorder::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
//!     // construct the analysis
//!     let analysis = Analysis::new()
//!             .structure("system.tpr")                   // structure file
//!             .trajectory("md.xtc")                      // trajectory file to analyze
//!             .index("index.ndx")                        // input ndx file
//!             .output_yaml("order.yaml")                 // output yaml file
//!             .output_tab("order.tab")                   // output table file
//!             .output_xvg("order.xvg")                   // pattern for output xvg files
//!             .output_csv("order.csv")                   // output csv file
//!             .analysis_type(AnalysisType::cgorder(      // analysis to perform
//!                "@membrane",                            // selection of beads
//!             ))
//!             .membrane_normal(Axis::Z)                  // membrane normal
//!             .begin(100_000.0)                          // starting time of analysis
//!             .end(200_000.0)                            // ending time of analysis
//!             .step(5)                                   // analyze every Nth frame of the trajectory
//!             .min_samples(100)                          // minimal required number of samples
//!             .n_threads(4)                              // number of threads to use for the analysis
//!             .leaflets(                                 // calculate order for individual leaflets
//!                 LeafletClassification::global(         // method for classifying lipids into leaflets
//!                     "@membrane",                       // selection of lipids for the calculation of membrane center
//!                     "name PO4"                         // selection of lipid heads
//!                 )
//!             )
//!             .ordermaps(                                // construct maps of order parameters
//!                 OrderMap::new()
//!                     .output_directory("ordermaps")     // directory to write the ordermaps to
//!                     .dim([                             // dimensions of the map
//!                         GridSpan::Manual {             // span of the x-dimension of the map
//!                             start: 5.0,                // map goes from 5 nm to...
//!                             end: 10.0,                 // ...10 nm along the x-dimension
//!                         },
//!                         GridSpan::Auto,                // span of the y-dimension of the map
//!                     ])
//!                     .bin_size([0.05, 0.2])             // size of each grid bin of the map
//!                     .min_samples(30)                   // minimal required number of samples per bin
//!                     .plane(Plane::XY)                  // orientation of the map
//!                     .build()?
//!             )                                          
//!             .build()?;                                 // constructing the analysis
//!
//!     // activate colog if you want logging (requires the `colog` crate)
//!     colog::init();
//!
//!     // run the analysis and write the output
//!     analysis.run()?;
//!
//!     Ok(())
//! }
//! ```
//! ***
//!
//! **See the [gorder manual](https://ladme.github.io/gorder-manual/) for full information about the capabilities of the crate.**

/// Version of the `gorder` crate.
pub const GORDER_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Message that should be added to every panic.
pub(crate) const PANIC_MESSAGE: &str =
    "\n\n\n            >>> THIS SHOULD NOT HAVE HAPPENED! PLEASE REPORT THIS ERROR <<<
(open an issue at 'github.com/Ladme/gorder/issues' or write an e-mail to 'ladmeb@gmail.com')\n\n";

/// Specifies leaflet a lipid is part of.
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
pub mod presentation;

use std::fmt::Display;

/// This module contains re-exported public structures of the `gorder` crate.
pub mod prelude {
    pub use super::input::{
        analysis::AnalysisBuilder, Analysis, AnalysisType, Axis, EstimateError, GridSpan,
        LeafletClassification, OrderMap, Plane,
    };

    pub use super::analysis::molecule::AtomType;

    pub use super::presentation::{
        aaresults::{AAAtomResults, AAMoleculeResults, AAOrderResults},
        cgresults::{CGMoleculeResults, CGOrderResults},
        AnalysisResults, BondResults, GridMapF32, Order, OrderCollection, OrderMapsCollection,
        PublicMoleculeResults, PublicOrderResults,
    };
}
