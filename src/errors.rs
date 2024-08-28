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
