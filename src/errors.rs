use std::path::Path;

use colored::{ColoredString, Colorize};
use groan_rs::errors::SelectError;
use thiserror::Error;

fn path_to_yellow(path: &Path) -> ColoredString {
    path.to_str().unwrap().yellow()
}

/// Errors that can occur inside the application itself.
#[derive(Error, Debug)]
pub enum ApplicationError {
    #[error("{} could not open the configuration file '{}'", "error:".red().bold(), .0.yellow())]
    CouldNotOpenConfig(String),
    #[error("{} could not understand the contents of the configuration file '{}' ({})", "error:".red().bold(), .0.yellow(), .1)]
    CouldNotParseConfig(String, serde_yaml::Error),
}

/// Errors that can occur when creating a `GridSpan` structure.
#[derive(Error, Debug)]
pub enum GridSpanError {
    #[error("{} the first coordinate for the grid span ('{}' nm) is higher than the second coordinate for the grid span ('{}' nm)", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    Invalid(f32, f32),

    #[error("{} grid span contains a negative value ('{}' nm)", "error:".red().bold(), .0.to_string().yellow())]
    Negative(f32),
}

/// Errors that can occur when analyzing system topology.
#[derive(Error, Debug)]
pub enum TopologyError {
    #[error("{}", .0)]
    InvalidQuery(SelectError),

    #[error("{} group '{}' is empty", "error:".red().bold(), .0.yellow())]
    EmptyGroup(String),

    #[error("{} some atoms are part of both '{}' (query: '{}') and '{}' (query: '{}')", "error:".red().bold(), .name1.yellow(), .query1.yellow(), .name2.yellow(), .query2.yellow())]
    AtomsOverlap {
        name1: String,
        query1: String,
        name2: String,
        query2: String,
    },

    #[error("{} residue number '{}' has an inconsistent name ('{}' vs '{}')", "error:".red().bold(), .0.to_string().yellow(), .1.yellow(), .2.yellow())]
    ResidueInconsistency(usize, String, String),

    #[error("{} molecule contains multiple atoms of the same name ('{}') of the same residue name ('{}')", "error:".red().bold(), .0.to_string().yellow(), .1.yellow())]
    DuplicateAtoms(usize, String),

    #[error("{} multiple bonds detected between atoms '{}' and '{}'", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    DuplicateBonds(usize, usize),

    #[error("{} \"hydrogen\" atom '{}' is bonded to two or more heavy atoms (the second heavy atom has a number '{}')", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    MultiBinding(usize, usize),

    #[error("{} molecule starting with atom number '{}' contains multiple head group atoms", "error:".red().bold(), .0.to_string().yellow())]
    MultipleHeads(usize),

    #[error("{} molecule starting with atom number '{}' contains no head group atom", "error:".red().bold(), .0.to_string().yellow())]
    NoHead(usize),

    #[error("{} molecule starting with atom number '{}' contains no methyl group atom", "error:".red().bold(), .0.to_string().yellow())]
    NoMethyl(usize),

    #[error("{} molecule starting with atom number '{}' contains a number of methyl group atoms ('{}') not consistent with other molecules ('{}')", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow(), .2.to_string().yellow())]
    InconsistentNumberOfMethyls(usize, usize, usize),

    #[error("{} system has undefined simulation box", "error:".red().bold())]
    UndefinedBox,

    #[error("{} the simulation box is not orthogonal", "error:".red().bold())]
    NotOrthogonalBox,

    #[error("{} all dimensions of the simulation box are zero", "error:".red().bold())]
    ZeroBox,
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

    #[error("{} atom with atom number '{}' has an undefined position", "error:".red().bold(), .0.to_string().yellow())]
    UndefinedPosition(usize),

    #[error("{} could not calculate global membrane center", "error:".red().bold())]
    InvalidGlobalMembraneCenter,

    #[error("{} could not calculate local membrane center for molecule with a head identifier number '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidLocalMembraneCenter(usize),
}

/// Errors that can occur while writing the results.
#[derive(Error, Debug)]
pub enum WriteError {
    #[error("{} could not create file '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreateFile(Box<Path>),

    #[error("{} could not create a backup for file '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotBackupFile(Box<Path>),

    #[error("{} could not write results in yaml format into '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotWriteYaml(Box<Path>),
}

/// Errors that can occur while writing the order maps.
#[derive(Error, Debug)]
pub enum OrderMapWriteError {
    #[error("{} could not create directory '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreateDirectory(Box<Path>),

    #[error("{} could not create a backup for directory '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotBackupDirectory(Box<Path>),

    #[error("{} could not remove an existing directory '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotRemoveDirectory(Box<Path>),

    #[error("{} could not create file '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreateFile(Box<Path>),

    #[error("{} could not write line into '{}'", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotWriteLine(Box<Path>),
}
