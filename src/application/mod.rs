// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This module contains the implementation of the `gorder` binary.

use std::{fs::File, io::BufWriter, path::Path};

use clap::Parser;
use colored::Colorize;
use gorder::{
    colog_info,
    errors::{ConfigError, WriteError},
    prelude::Analysis,
    GORDER_VERSION,
};
use std::io::Write;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "Calculate order parameters for any atomistic, united-atom, or coarse grained system."
)]
pub struct Args {
    #[arg(
        help = "Config yaml file",
        long_help = "Configuration yaml file specifying the analysis options."
    )]
    pub config: String,

    #[arg(
        long = "silent",
        action,
        help = "Suppress standard output",
        default_value_t = false,
        long_help = "Suppress all standard output generated by the 'gorder' program, except for error messages.
This option can also be specified in the configuration file."
    )]
    pub silent: bool,

    #[arg(
        long = "overwrite",
        action,
        help = "Overwrite existing files and directories with the same name",
        default_value_t = false,
        long_help = "Overwrite existing files and directories with the same name as the output. No backup copies will be created.
This option can also be specified in the configuration file."
    )]
    pub overwrite: bool,

    #[arg(
        long = "export-config",
        help = "Export a YAML file with all analysis options, including defaults",
        long_help = "Export a YAML file containing all parameters used for the analysis.

This includes the parameters specified in the input configuration file, on the command line, 
as well as the default parameters automatically set by `gorder`.
The argument should specify a path where the exported configuration is written."
    )]
    pub output_filename: Option<String>,
}

/// Get arguments, parse input file, run the analysis and write the results.
/// Returns `true` if run successfully, else returns `false`.
pub(crate) fn run() -> bool {
    let args = Args::parse();

    colog::init();

    let mut analysis = match Analysis::from_file(&args.config) {
        Ok(x) => x,
        Err(e) => {
            log::error!("{}", e);
            return false;
        }
    };

    // yaml output must always be specified in the config file for the application
    if analysis.output_yaml().is_none() {
        log::error!("{}", ConfigError::NoYamlOutput(args.config.clone()));
        return false;
    }

    // ordermaps output must always be specified in the config file for the application
    if let Some(maps) = analysis.map() {
        if maps.output_directory().is_none() {
            log::error!("{}", ConfigError::NoOrdermapsOutput(args.config.clone()));
            return false;
        }
    }

    if !analysis.silent() {
        analysis.set_silent(args.silent);
    }

    if !analysis.overwrite() {
        analysis.set_overwrite(args.overwrite);
    }

    let silent = analysis.silent();

    if silent {
        log::set_max_level(log::LevelFilter::Error);
    } else {
        let header = format!(">>> GORDER v{} <<<", GORDER_VERSION).bold();
        println!("\n{}\n", header);
        colog_info!("Read config file '{}'.", args.config);
    }

    let result = match analysis.run() {
        Err(e) => {
            log::error!("{}", e);
            display_result(false, silent);
            return false;
        }
        Ok(x) => x,
    };

    // write the results into an output file
    if let Err(e) = result.write() {
        log::error!("{}", e);
        display_result(false, silent);
        return false;
    }

    // writing output analysis, if requested
    if let Some(output) = args.output_filename {
        let analysis = result.analysis();
        if let Err(e) = export_analysis_options(analysis, output, analysis.overwrite()) {
            log::error!(
                "Analysis completed successfully, but exporting the analysis options failed!\n    {}",
                e
            );
        }
    }

    display_result(true, silent);
    true
}

fn display_result(result: bool, silent: bool) {
    if silent {
        return;
    }

    match result {
        true => {
            let prefix = format!(
                "{}{}{}",
                "[".to_string().blue().bold(),
                "✔".to_string().bright_green().bold(),
                "]".to_string().blue().bold()
            );
            let message = "ANALYSIS COMPLETED".to_string().bright_green().bold();
            println!("{} {}\n", prefix, message);
        }
        false => {
            let prefix = format!(
                "{}{}{}",
                "[".to_string().blue().bold(),
                "✖".to_string().red().bold(),
                "]".to_string().blue().bold()
            );
            let message = "ANALYSIS FAILED".to_string().red().bold();
            println!("{} {}\n", prefix, message);
        }
    }
}

fn export_analysis_options(
    analysis: &Analysis,
    filename: impl AsRef<Path>,
    overwrite: bool,
) -> Result<(), WriteError> {
    colog_info!(
        "Exporting all analysis options into file '{}'...",
        filename.as_ref().to_str().unwrap()
    );

    log::logger().flush();

    if filename.as_ref().exists() {
        if !overwrite {
            log::warn!(
                "Output file '{}' already exists. Backing it up.",
                filename.as_ref().to_str().unwrap()
            );
            backitup::backup(filename.as_ref())
                .map_err(|_| WriteError::CouldNotBackupFile(Box::from(filename.as_ref())))?;
        } else {
            log::warn!(
                "Output file '{}' already exists. It will be overwritten as requested.",
                filename.as_ref().to_str().unwrap()
            );
        }
    }

    let file = File::create(filename.as_ref())
        .map_err(|_| WriteError::CouldNotCreateFile(Box::from(filename.as_ref())))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "# Analysis options generated by 'gorder v{}'.",
        GORDER_VERSION
    )
    .map_err(|_| WriteError::CouldNotWriteLine(Box::from(filename.as_ref())))?;

    serde_yaml::to_writer(&mut writer, analysis).map_err(WriteError::CouldNotExportAnalysis)
}
