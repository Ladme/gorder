// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use std::sync::Arc;

use gorder_core::input::analysis::TrajectoryInput as RsTrajectoryInput;
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

/// Request the calculation of atomistic order parameters.
///
/// Attributes
/// ----------
/// heavy_atoms : str
///     Selection query specifying the heavy atoms to be used in the analysis (typically carbon atoms in lipid tails).
/// hydrogens : str 
///     Selection query specifiying the hydrogen atoms to be used in the analysis (only those bonded to heavy atoms will be considered).
///
/// Notes
/// ------
/// - Atoms should be specified using the `groan selection language <https://ladme.github.io/gsl-guide>`_.
/// - Order parameters are calculated for bonds between `heavy_atoms` and `hydrogens`. These bonds are detected automatically.
/// - The order parameters for heavy atoms are determined by averaging the order parameters of the corresponding bonds.
#[pyclass]
#[derive(Clone)]
pub struct AAOrder(RsAnalysisType);

#[pymethods]
impl AAOrder {
    #[new]
    pub fn new(heavy_atoms: &str, hydrogens: &str) -> Self {
        Self(RsAnalysisType::aaorder(heavy_atoms, hydrogens))
    }
}

/// Request the calculation of coarse-grained order parameters.
///
/// Attributes
/// ----------
/// beads : str
///     Selection query specifying the coarse-grained beads to be used in the analysis.
///
/// Notes
/// -----
/// - Beads should be specified using the `groan selection language <https://ladme.github.io/gsl-guide>`_.
/// - Order parameters are calculated for bonds between individual `beads`. These bonds are detected automatically.
#[pyclass]
#[derive(Clone)]
pub struct CGOrder(RsAnalysisType);

#[pymethods]
impl CGOrder {
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

#[derive(Clone)]
pub struct TrajectoryInput(RsTrajectoryInput);

impl<'source> FromPyObject<'source> for TrajectoryInput {
    fn extract_bound(obj: &Bound<'source, PyAny>) -> PyResult<Self> {
        if let Ok(string) = obj.extract::<String>() {
            return Ok(TrajectoryInput(string.into()));
        }

        if let Ok(list) = obj.extract::<Vec<String>>() {
            return Ok(TrajectoryInput(list.into()));
        }

        Err(ConfigError::new_err("expected a string or a list of strings"))
    }
}

/// Class describing all the parameters of the analysis.
///
/// Attributes
/// ----------
/// structure : str
///     Path to a TPR (recommended), PDB, GRO, or PQR file containing the structure and topology of the system.
/// trajectory : Union[str, list[str]]
///     Path to an XTC (recommened), TRR, GRO, PDB, Amber NetCDF, DCD, or LAMMPSTRJ trajectory file to be analyzed.
///     You can provide multiple XTC or TRR trajectories and these will be seamlessly concatenated.
/// analysis_type : Union[AAOrder, CGOrder]
///     Type of analysis to perform (e.g., AAOrder or CGOrder).
/// bonds : Optional[str], default=None
///     Path to a file containing bonding information. If specified, this overrides bonds from the structure file.
/// index : Optional[str], default=None
///     Path to an NDX file specifying groups in the system.
/// output_yaml : Optional[str], default=None
///     Path to an output YAML file containing the full analysis results.
/// output_tab : Optional[str], default=None
///     Path to an output TABLE file with human-readable results.
/// output_xvg : Optional[str], default=None
///     Filename pattern for output XVG files storing results. Each molecule type gets a separate file.
/// output_csv : Optional[str], default=None
///     Path to an output CSV file containing analysis results.
/// membrane_normal : Optional[Union[str, dict, DynamicNormal]], default=None
///     Direction of the membrane normal. 
///     Allowed values are `x`, `y`, `z`, path to file, dictionary specifying manual membrane normals or an instance of `DynamicNormal`.
///     Defaults to the z-axis if not specified.
/// begin : Optional[float], default=None
///     Starting time of the trajectory analysis in picoseconds (ps). Defaults to the beginning of the trajectory.
/// end : Optional[float], default=None
///     Ending time of the trajectory analysis in picoseconds (ps). Defaults to the end of the trajectory.
/// step : Optional[int], default=None
///     Step size for analysis. Every Nth frame will be analyzed. Defaults to 1.
/// min_samples : Optional[int], default=None
///     Minimum number of samples required for each heavy atom or bond type to compute its order parameter. Defaults to 1.
/// n_threads : Optional[int], default=None
///     Number of threads to use for analysis. Defaults to 1.
/// leaflets : Optional[Union[GlobalClassification, LocalClassification, IndividualClassification, ManualClassification, NdxClassification]], default=None
///     Defines how lipids are assigned to membrane leaflets. If provided, order parameters are calculated per leaflet.
/// ordermap : Optional[OrderMap], default=None
///     Specifies parameters for ordermap calculations. If not provided, ordermaps are not generated.
/// estimate_error : Optional[EstimateError], default=None
///     Enables error estimation for each bond if specified.
/// geometry : Optional[Union[Cuboid, Cylinder, Sphere]], default=None
///     Defines a specific region in the simulation box for order parameter calculations. Defaults to the entire system.
/// handle_pbc : Optional[bool], default=True
///     If False, ignores periodic boundary conditions (PBC). Defaults to True.
/// silent : Optional[bool], default=False
///     If True, suppresses standard output messages during analysis.
/// overwrite : Optional[bool], default=False
///     If True, overwrites existing output files and directories without backups.
#[pyclass]
pub struct Analysis(RsAnalysis);

#[pymethods]
impl Analysis {
    #[allow(clippy::too_many_arguments)]
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
        trajectory: TrajectoryInput,
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
            .trajectory(trajectory.0)
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
        if self.0.silent() {
            log::set_max_level(log::LevelFilter::Error);
        }

        match self.0.clone().run() {
            Err(e) => Err(AnalysisError::new_err(e.to_string())),
            Ok(x) => Ok(AnalysisResults(Arc::new(x))),
        }
    }
}
