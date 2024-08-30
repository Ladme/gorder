use colored::Colorize;
use thiserror::Error;

/// Errors that can occur when creating a `GridSpan` structure.
#[derive(Error, Debug)]
pub enum GridSpanError {
    #[error("{} the first coordinate for the grid span ('{}' nm) is higher than the second coordinate for the grid span ('{}' nm).", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    Invalid(f32, f32),
    #[error("{} grid span contains a negative value ('{}' nm)", "error:".red().bold(), .0.to_string().yellow())]
    Negative(f32),
}

/// Errors that can occur when constructing a `Topology` structure.
#[derive(Error, Debug)]
pub enum TopologyError {
    #[error("{} groan selection language query '{}' is invalid.", "error:".red().bold(), .0.yellow())]
    InvalidQuery(String),

    #[error("{} group '{}' is empty.", "error:".red().bold(), .0.yellow())]
    EmptyGroup(String),

    #[error("{} an atom overlap detected between '{}' (query: '{}') and '{}' (query: '{}').", "error:".red().bold(), .name1.yellow(), .query1.yellow(), .name2.yellow(), .query2.yellow())]
    AtomsOverlap {
        name1: String,
        query1: String,
        name2: String,
        query2: String,
    },

    #[error("{} residue number '{}' has an inconsistent name ('{}' vs '{}').", "error:".red().bold(), .0.to_string().yellow(), .1.yellow(), .2.yellow())]
    ResidueInconsistency(usize, String, String),

    #[error("{} molecule contains multiple atoms of the same name ('{}') of the same residue name ('{}').", "error:".red().bold(), .0.to_string().yellow(), .1.yellow())]
    DuplicateAtoms(usize, String),

    #[error("{} multiple bonds detected between atoms '{}' and '{}'.", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    DuplicateBonds(usize, usize),

    #[error("{} \"hydrogen\" atom '{}' is bonded to two or more heavy atoms (the second heavy atom has a number '{}')", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    MultiBinding(usize, usize),
}
