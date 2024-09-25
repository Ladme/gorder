// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains the implementation of the `gorder` binary.

use clap::Parser;
use colored::Colorize;
use gorder::{errors::ApplicationError, Analysis, GORDER_VERSION};

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
        long_help = "Configuration yaml file specifying the analysis settings."
    )]
    pub config: String,
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args = Args::parse();

    let config_str = std::fs::read_to_string(&args.config)
        .map_err(|_| ApplicationError::CouldNotReadConfig(args.config.clone()))?;
    let analysis: Analysis = serde_yaml::from_str(&config_str)
        .map_err(|_| ApplicationError::CouldNotReadConfig(args.config.clone()))?;

    if analysis.silent() {
        colog::basic_builder()
            .filter(None, log::LevelFilter::Error)
            .init();
    } else {
        colog::init();
        let header = format!(">>> GORDER v{} <<<", GORDER_VERSION).bold();
        println!("\n{}\n", header);
        log::info!("Read config file '{}'.", args.config);
    }

    let result = analysis.run();

    if !analysis.silent() {
        match &result {
            Ok(_) => {
                let prefix = format!(
                    "{}{}{}",
                    "[".to_string().blue().bold(),
                    "✔".to_string().bright_green().bold(),
                    "]".to_string().blue().bold()
                );
                let message = "ANALYSIS COMPLETED".to_string().bright_green().bold();
                println!("{} {}", prefix, message);
            }
            Err(e) => {
                log::error!("{}", e);

                let prefix = format!(
                    "{}{}{}",
                    "[".to_string().blue().bold(),
                    "✖".to_string().red().bold(),
                    "]".to_string().blue().bold()
                );
                let message = "ANALYSIS FAILED".to_string().red().bold();
                println!("{} {}", prefix, message);
            }
        }
    }

    result
}
