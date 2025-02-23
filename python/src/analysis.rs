// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use std::process;

use gorder_core::input::AnalysisType as RsAnalysisType;
use gorder_core::prelude::Analysis as RsAnalysis;
use gorder_core::prelude::AnalysisBuilder as RsAnalysisBuilder;
use pyo3::prelude::*;

use crate::estimate_error::EstimateError;
use crate::geometry::Geometry;
use crate::leaflets::LeafletClassification;
use crate::normal::MembraneNormal;
use crate::ordermap::OrderMap;
use crate::results::AnalysisResults;
use crate::{AnalysisError, ConfigError};

/// Structure describing the calculation of atomistic order paramters.
#[pyclass]
#[derive(Clone)]
pub struct AAOrder(RsAnalysisType);

#[pymethods]
impl AAOrder {
    /// Request calculation of atomistic order parameters.
    ///
    /// ## Arguments
    /// - `heavy_atoms` - specification of heavy atoms to be used in the analysis (usually carbon atoms of lipid tails)
    /// - `hydrogens` - specification of hydrogens to be used in the analysis (only atoms connected to heavy atoms will be used)
    ///
    /// ## Notes
    /// - To specify the atoms, use the [groan selection language](https://docs.rs/groan_rs/latest/groan_rs/#groan-selection-language).
    /// - Order parameters will be calculated for bonds between the `heavy_atoms` and `hydrogens`. These bonds are detected automatically.
    /// - Order parameters for heavy atoms are then calculated by averaging the order parameters of the relevant bonds.
    #[new]
    pub fn new(heavy_atoms: &str, hydrogens: &str) -> Self {
        Self(RsAnalysisType::aaorder(heavy_atoms, hydrogens))
    }
}

/// Structure describing the calculation of coarse-grained order paramters.
#[pyclass]
#[derive(Clone)]
pub struct CGOrder(RsAnalysisType);

#[pymethods]
impl CGOrder {
    /// Request calculation of coarse-grained order parameters.
    ///
    /// ## Arguments
    /// - `beads` - specification of coarse-grained beads
    ///
    /// ## Notes
    /// - To specify the beads, use the [groan selection language](https://docs.rs/groan_rs/latest/groan_rs/#groan-selection-language).
    /// - Order parameters will be calculated for bonds between the individual `beads`. These bonds are detected automatically.
    #[new]
    pub fn new(beads: &str) -> Self {
        Self(RsAnalysisType::cgorder(beads))
    }
}

#[derive(Clone)]
pub struct AnalysisType(RsAnalysisType);

impl<'source> FromPyObject<'source> for AnalysisType {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        if let Ok(analysis_type) = obj.extract::<AAOrder>() {
            return Ok(AnalysisType(analysis_type.0));
        }

        if let Ok(analysis_type) = obj.extract::<CGOrder>() {
            return Ok(AnalysisType(analysis_type.0));
        }

        Err(ConfigError::new_err(
            "expected an instance of AAOrder or CGOrder as AnalysisType",
        ))
    }
}

/// Structure holding all the information necessary to perform the analysis.
#[pyclass]
pub struct Analysis(RsAnalysis);

#[pymethods]
impl Analysis {
    // TODO: documentation
    #[new]
    #[pyo3(signature = (
        structure, 
        trajectory, 
        analysis_type,
        bonds=None, 
        index=None, 
        output_yaml=None, 
        output_tab=None, 
        output_xvg=None, 
        output_csv=None, 
        membrane_normal=None, 
        begin=None, 
        end=None, 
        step=None, 
        min_samples=None, 
        n_threads=None, 
        leaflets=None, 
        ordermap=None, 
        estimate_error=None, 
        geometry=None, 
        handle_pbc=None, 
        silent=None, 
        overwrite=None))]
    pub fn new(
        structure: &str,
        trajectory: &str,
        analysis_type: AnalysisType,
        bonds: Option<&str>,
        index: Option<&str>,
        output_yaml: Option<&str>,
        output_tab: Option<&str>,
        output_xvg: Option<&str>,
        output_csv: Option<&str>,
        membrane_normal: Option<MembraneNormal>,
        begin: Option<f32>,
        end: Option<f32>,
        step: Option<usize>,
        min_samples: Option<usize>,
        n_threads: Option<usize>,
        leaflets: Option<LeafletClassification>,
        ordermap: Option<OrderMap>,
        estimate_error: Option<EstimateError>,
        geometry: Option<Geometry>,
        handle_pbc: Option<bool>,
        silent: Option<bool>,
        overwrite: Option<bool>,
    ) -> PyResult<Self> {
        let mut builder: RsAnalysisBuilder = RsAnalysis::builder();
        builder
            .structure(structure)
            .trajectory(trajectory)
            .analysis_type(analysis_type.0);

        apply_if_some!(
            builder,
            bonds           => bonds,
            index           => index,
            output_yaml     => output_yaml,
            output_tab      => output_tab,
            output_xvg      => output_xvg,
            output_csv      => output_csv,
            begin           => begin,
            end             => end,
            step            => step,
            min_samples     => min_samples,
            n_threads       => n_threads,
            handle_pbc      => handle_pbc,
        );

        apply_inner_if_some!(
            builder,
            membrane_normal => membrane_normal,
            leaflets        => leaflets,
            estimate_error  => estimate_error,
            geometry        => geometry,
            ordermap        => map,
        );

        if let Some(true) = silent {
            builder.silent();
        }

        if let Some(true) = overwrite {
            builder.overwrite();
        }

        let inner = builder
            .build()
            .map_err(|e| ConfigError::new_err(e.to_string()))?;
        Ok(Analysis(inner))
    }

    /// Perform the analysis.
    pub fn run(&mut self) -> PyResult<AnalysisResults> {
        colog::init();
        if self.0.silent() {
            log::set_max_level(log::LevelFilter::Error);
        }

        // register a Ctrl-C handler
        ctrlc::set_handler(|| {
            process::exit(1);
        })
        .unwrap_or_else(|e| panic!("FATAL GORDER ERROR | python::analysis::Analysis | Could not set up the CTRL-C handler: {}.", e));

        match self.0.clone().run() {
            Err(e) => Err(AnalysisError::new_err(e.to_string())),
            Ok(x) => Ok(AnalysisResults(x)),
        }
    }
}
