// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains the implementation of the `gorder` binary.

use clap::Parser;
use gorder::Analysis;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "Calculate order parameters for any atomistic or coarse grained system."
)]
pub struct Args {
    #[arg(
        help = "Config yaml file",
        long_help = "Configuration yaml file specifying the analysis options."
    )]
    pub config: String,
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args = Args::parse();
    // TODO: use custom error type
    let config_str = std::fs::read_to_string(&args.config)?;
    let analysis: Analysis = serde_yaml::from_str(&config_str)?;
    analysis.run()
}
