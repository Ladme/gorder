// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains error types that can be returned by the `gorder` crate.

use std::path::Path;

use colored::{ColoredString, Colorize};
use groan_rs::errors::SelectError;
use thiserror::Error;

fn path_to_yellow(path: &Path) -> ColoredString {
    path.to_str().unwrap().yellow()
}

/// Errors that can occur when creating a `GridSpan` structure.
#[derive(Error, Debug)]
pub enum GridSpanError {
    #[error("{} the first coordinate for the grid span ('{}' nm) is higher than the second coordinate for the grid span ('{}' nm)", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow()
    )]
    Invalid(f32, f32),
}

/// Errors that can occur when analyzing system topology.
#[derive(Error, Debug)]
pub enum TopologyError {
    #[error("{}", .0)]
    InvalidQuery(SelectError),

    #[error("{} group '{}' is empty", "error:".red().bold(), .0.yellow())]
    EmptyGroup(String),

    #[error("{} {} atoms are part of both '{}' (query: '{}') and '{}' (query: '{}')", "error:".red().bold(), .n_overlapping.to_string().yellow(), .name1.yellow(), .query1.yellow(), .name2.yellow(), .query2.yellow()
    )]
    AtomsOverlap {
        n_overlapping: usize,
        name1: String,
        query1: String,
        name2: String,
        query2: String,
    },

    #[error("{} molecule starting with atom index '{}' contains multiple head group atoms", "error:".red().bold(), .0.to_string().yellow()
    )]
    MultipleHeads(usize),

    #[error("{} molecule starting with atom index '{}' contains no head group atom", "error:".red().bold(), .0.to_string().yellow()
    )]
    NoHead(usize),

    #[error("{} molecule starting with atom index '{}' contains no methyl group atom", "error:".red().bold(), .0.to_string().yellow()
    )]
    NoMethyl(usize),

    #[error("{} molecule starting with atom index '{}' contains a number of methyl group atoms ('{}') not consistent with other molecules ('{}')", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow(), .2.to_string().yellow()
    )]
    InconsistentNumberOfMethyls(usize, usize, usize),

    #[error("{} system has undefined simulation box", "error:".red().bold())]
    UndefinedBox,

    #[error("{} the simulation box is not orthogonal", "error:".red().bold())]
    NotOrthogonalBox,

    #[error("{} all dimensions of the simulation box are zero", "error:".red().bold())]
    ZeroBox,

    #[error("{}", .0)]
    OrderMapError(OrderMapConfigError),
}

/// Errors that can occur while analyzing the trajectory.
#[derive(Error, Debug)]
pub enum AnalysisError {
    #[error("{} system has undefined simulation box", "error:".red().bold())]
    UndefinedBox,

    #[error("{} the simulation box is not orthogonal", "error:".red().bold())]
    NotOrthogonalBox,

    #[error("{} all dimensions of the simulation box are zero", "error:".red().bold())]
    ZeroBox,

    #[error("{} atom with atom index '{}' has an undefined position", "error:".red().bold(), .0.to_string().yellow()
    )]
    UndefinedPosition(usize),

    #[error("{} could not calculate global membrane center", "error:".red().bold())]
    InvalidGlobalMembraneCenter,

    #[error("{} could not calculate local membrane center for molecule with a head identifier index '{}'", "error:".red().bold(), .0.to_string().yellow()
    )]
    InvalidLocalMembraneCenter(usize),
}

/// Errors that can occur while writing the results.
#[derive(Error, Debug)]
pub enum WriteError {
    #[error("{} could not create file '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreateFile(Box<Path>),

    #[error("{} could not create a backup for file '{}'", "error:".red().bold(), path_to_yellow(.0)
    )]
    CouldNotBackupFile(Box<Path>),

    #[error("{} could not write output file header into '{}'", "error:".red().bold(), path_to_yellow(.0)
    )]
    CouldNotWriteHeader(Box<Path>),

    #[error("{} could not write results in yaml format into '{}'", "error:".red().bold(), path_to_yellow(.0)
    )]
    CouldNotWriteYaml(Box<Path>),

    #[error("{} could not write results in the output file ({})", "error:".red().bold(), .0)]
    CouldNotWriteResults(std::io::Error),
}

/// Errors that can occur while writing the order maps.
#[derive(Error, Debug)]
pub enum OrderMapWriteError {
    #[error("{} could not create directory '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreateDirectory(Box<Path>),

    #[error("{} could not create a backup for directory '{}'", "error:".red().bold(), path_to_yellow(.0)
    )]
    CouldNotBackupDirectory(Box<Path>),

    #[error("{} could not remove an existing directory '{}'", "error:".red().bold(), path_to_yellow(.0)
    )]
    CouldNotRemoveDirectory(Box<Path>),

    #[error("{} could not create file '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreateFile(Box<Path>),

    #[error("{} could not write line into '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotWriteLine(Box<Path>),
}

/// Errors that can occur when constructing an `Analysis` structure from the provided configuration.
#[derive(Error, Debug)]
pub enum ConfigError {
    #[error("{} could not open the configuration file '{}'", "error:".red().bold(), .0.yellow())]
    CouldNotOpenConfig(String),

    #[error("{} could not understand the contents of the configuration file '{}' ({})", "error:".red().bold(), .0.yellow(), .1
    )]
    CouldNotParseConfig(String, serde_yaml::Error),

    #[error("{} the specified value of '{}' is invalid (must be positive)", "error:".red().bold(), "step".yellow()
    )]
    InvalidStep,

    #[error("{} the specified value of '{}' is invalid (must be positive)", "error:".red().bold(), "min_samples".yellow()
    )]
    InvalidMinSamples,

    #[error("{} the specified value of '{}' is invalid (must be positive)", "error:".red().bold(), "n_threads".yellow()
    )]
    InvalidNThreads,

    #[error("{} invalid values of '{}' and '{}' (begin is higher than end)",
            "error:".red().bold(),
            "begin".yellow(),
            "end".yellow())]
    InvalidBeginEnd,

    #[error("{}", .0)]
    InvalidOrderMap(OrderMapConfigError),
}

/// Errors that can occur when constructing an `OrderMap` structure from the provided configuration.
#[derive(Error, Debug)]
pub enum OrderMapConfigError {
    #[error("{} the specified value of '{}' inside '{}' is invalid (must be positive)", 
            "error:".red().bold(), 
            "min_samples".yellow(), 
            "ordermap".yellow())]
    InvalidMinSamples,

    #[error("{} invalid span of '{}': minimum ('{}') is higher than maximum ('{}')",
            "error:".red().bold(),
            "ordermap".yellow(),
            .0.to_string().yellow(), .1.to_string().yellow())]
    InvalidGridSpan(f32, f32),

    #[error("{} invalid bin size of '{}': value is '{}', must be positive", 
            "error:".red().bold(), 
            "ordermap".yellow(), 
            .0.to_string().yellow())]
    InvalidBinSize(f32),

    #[error("{} invalid bin size of '{}': bin size of '{}x{}' is larger than grid span of '{}x{}'",
            "error:".red().bold(),
            "ordermap".yellow(),
            .0.0.to_string().yellow(),
            .0.1.to_string().yellow(),
            .1.0.to_string().yellow(),
            .1.1.to_string().yellow())]
    BinTooLarge((f32, f32), (f32, f32)),
}

/// Errors that can occur when estimating the error of the calculation.
#[derive(Error, Debug)]
pub enum ErrorEstimationError {
    #[error("{} collected '{}' block(s) for error estimation; at least 6 blocks are needed ({} decrease the `block_size` property to {})",
    "error:".red().bold(), .0.to_string().yellow(), "hint:".purple().bold(), .1.to_string().bright_purple()
    )]
    NotEnoughBlocks(usize, usize),

    #[error("{} read {} trajectory frame(s) which is not enough to estimate an error ({} collect more data)",
    "error:".red().bold(), .0.to_string().yellow(), "hint:".purple().bold())]
    NotEnoughData(usize),
}
