// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the main `Analysis` structure and its methods.

use std::fs::read_to_string;
use std::path::Path;

use derive_builder::Builder;
use getset::{CopyGetters, Getters};
use getset::{MutGetters, Setters};
use serde::Deserialize;

use crate::errors::ConfigError;

use super::Axis;
use super::LeafletClassification;
use super::OrderMap;

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum AnalysisType {
    AAOrder {
        heavy_atoms: String,
        hydrogens: String,
    },
    CGOrder {
        #[serde(alias = "atoms")]
        beads: String,
    },
}

impl AnalysisType {
    pub fn aaorder(heavy_atoms: &str, hydrogens: &str) -> Self {
        Self::AAOrder {
            heavy_atoms: heavy_atoms.to_owned(),
            hydrogens: hydrogens.to_owned(),
        }
    }

    pub fn cgorder(beads: &str) -> Self {
        Self::CGOrder {
            beads: beads.to_owned(),
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::AAOrder {
                heavy_atoms: _,
                hydrogens: _,
            } => "all-atom order parameters",
            Self::CGOrder { beads: _ } => "coarse-grained order parameters",
        }
    }
}

/// Structure holding all the information necessary to perform the specified analysis.
#[derive(Debug, Clone, Builder, Getters, CopyGetters, Setters, MutGetters, Deserialize)]
#[serde(deny_unknown_fields)]
#[builder(build_fn(validate = "Self::validate"))]
pub struct Analysis {
    /// Path to TPR file containing the topology of the system.
    #[builder(setter(into))]
    #[getset(get = "pub")]
    structure: String,
    /// Path to XTC trajectory file containing the trajectory to be analyzed.
    #[builder(setter(into))]
    #[getset(get = "pub")]
    trajectory: String,
    /// Path to NDX file containing the groups associated with the system.
    /// Optional parameter.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    index: Option<String>,
    /// Path to an output YAML file where the results of the analysis will be written.
    #[builder(setter(into))]
    #[getset(get = "pub")]
    #[serde(alias = "output")]
    output_yaml: String,
    /// Path to an output TABLE file where the results of the analysis will be written in a human readable format.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    output_tab: Option<String>,
    /// Filename pattern of the output XVG files where the results of the analysis will be written.
    /// One xvg file will be generated for each detected molecule type and the name of this molecule type
    /// will be appended to the provided filename pattern (between the file stem and the file extension).
    ///
    /// Example: filename pattern 'order.xvg' may be transformed into 'order_POPC.xvg'.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    output_xvg: Option<String>,
    /// Path to an output CSV file where the results of the analysis will be written.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    output_csv: Option<String>,
    /// Type of the analysis to perform (AAOrder / CGOrder).
    #[getset(get = "pub")]
    #[serde(alias = "type")]
    analysis_type: AnalysisType,
    /// Direction of the membrane normal.
    /// If not provided, the default value is 'Axis::Z'.
    #[builder(setter(into), default = "Axis::Z")]
    #[serde(default = "default_membrane_normal")]
    #[getset(get_copy = "pub")]
    membrane_normal: Axis,
    /// Starting time of the trajectory analysis (in ps).
    /// If not specified, the analysis starts at the beginning of the trajectory.
    #[builder(default = "0.0")]
    #[serde(default = "default_begin", alias = "start")]
    #[getset(get_copy = "pub")]
    begin: f32,
    /// Ending time of the trajectory analysis (in ps).
    /// If not specified, the analysis ends at the end of the trajectory.
    #[builder(default = "f32::INFINITY")]
    #[serde(default = "default_end")]
    #[getset(get_copy = "pub")]
    end: f32,
    /// Only every Nth frame of the simulation trajectory will be analyzed.
    /// If not specified, each frame of the trajectory will be analyzed.
    #[builder(default = "1")]
    #[serde(default = "default_one")]
    #[getset(get_copy = "pub")]
    step: usize,
    /// Minimal number of samples for each heavy atom required to calculate order parameter for it.
    /// If not specified, the default value is 1.
    #[builder(default = "1")]
    #[serde(default = "default_one")]
    #[getset(get_copy = "pub")]
    min_samples: usize,
    /// Number of threads to use to perform the analysis.
    /// If not specified, the default value is 1.
    #[builder(default = "1")]
    #[serde(default = "default_one")]
    #[getset(get_copy = "pub")]
    n_threads: usize,
    /// Specifies how to assign lipids into membrane leaflets.
    /// If specified, order parameters will be calcualted for the individual membrane leaflets
    /// as well as for the entire membrane.
    /// If not specified, only values for the entire membrane will be calculated.
    #[builder(setter(strip_option), default)]
    #[getset(get = "pub")]
    leaflets: Option<LeafletClassification>,
    /// If provided, specifies the properties of an ordermap that shall be calculated.
    /// If not specified, no map will be calculated.
    #[builder(setter(strip_option), default)]
    #[serde(alias = "maps", alias = "ordermap", alias = "ordermaps")]
    #[getset(get = "pub", get_mut = "pub(crate)")]
    map: Option<OrderMap>,
    /// Be silent. Print nothing to the standard output during the analysis.
    #[builder(setter(custom), default = "false")]
    #[serde(default = "default_false")]
    #[getset(get_copy = "pub", set = "pub")]
    silent: bool,
    /// Do not make backups. Overwrite all output files and directories.
    #[builder(setter(custom), default = "false")]
    #[serde(default = "default_false")]
    #[getset(get_copy = "pub", set = "pub")]
    overwrite: bool,
}

fn default_membrane_normal() -> Axis {
    Axis::Z
}

fn default_begin() -> f32 {
    0.0
}

fn default_end() -> f32 {
    f32::INFINITY
}

fn default_one() -> usize {
    1
}

fn default_false() -> bool {
    false
}

fn validate_step(step: usize) -> Result<(), ConfigError> {
    if step == 0 {
        Err(ConfigError::InvalidStep)
    } else {
        Ok(())
    }
}

fn validate_min_samples(samples: usize) -> Result<(), ConfigError> {
    if samples == 0 {
        Err(ConfigError::InvalidMinSamples)
    } else {
        Ok(())
    }
}

fn validate_n_threads(n_threads: usize) -> Result<(), ConfigError> {
    if n_threads == 0 {
        Err(ConfigError::InvalidNThreads)
    } else {
        Ok(())
    }
}

fn validate_begin_end(begin: f32, end: f32) -> Result<(), ConfigError> {
    if begin > end {
        Err(ConfigError::InvalidBeginEnd)
    } else {
        Ok(())
    }
}

impl Analysis {
    pub fn new() -> AnalysisBuilder {
        AnalysisBuilder::default()
    }

    pub fn from_file(path: impl AsRef<Path>) -> Result<Analysis, ConfigError> {
        let string = read_to_string(&path).map_err(|_| {
            ConfigError::CouldNotOpenConfig(path.as_ref().to_str().unwrap().to_owned())
        })?;
        let analysis: Analysis = serde_yaml::from_str(&string).map_err(|e| {
            ConfigError::CouldNotParseConfig(path.as_ref().to_str().unwrap().to_owned(), e)
        })?;

        analysis.validate()?;
        Ok(analysis)
    }

    /// Check that the Analysis structure is valid. Used after deserialization from config yaml file.
    fn validate(&self) -> Result<(), ConfigError> {
        validate_step(self.step)?;
        validate_min_samples(self.min_samples)?;
        validate_n_threads(self.n_threads)?;
        validate_begin_end(self.begin, self.end)?;

        // check the validity of the order map, if present
        if let Some(ref map) = self.map {
            map.validate()
                .map_err(|e| ConfigError::InvalidOrderMap(e))?;
        }

        Ok(())
    }

    pub fn heavy_atoms(&self) -> Option<&String> {
        match &self.analysis_type {
            AnalysisType::CGOrder { beads: _ } => None,
            AnalysisType::AAOrder {
                heavy_atoms,
                hydrogens: _,
            } => Some(heavy_atoms),
        }
    }

    pub fn hydrogens(&self) -> Option<&String> {
        match &self.analysis_type {
            AnalysisType::CGOrder { beads: _ } => None,
            AnalysisType::AAOrder {
                heavy_atoms: _,
                hydrogens,
            } => Some(hydrogens),
        }
    }

    pub fn beads(&self) -> Option<&String> {
        match &self.analysis_type {
            AnalysisType::CGOrder { beads } => Some(beads),
            AnalysisType::AAOrder {
                heavy_atoms: _,
                hydrogens: _,
            } => None,
        }
    }

    #[inline(always)]
    pub fn atoms(&self) -> Option<&String> {
        self.beads()
    }

    /// Alias for `output_yaml`.
    #[inline(always)]
    pub fn output(&self) -> &String {
        &self.output_yaml
    }
}

impl AnalysisBuilder {
    /// Be silent. Print nothing to the standard output during the analysis.
    #[inline(always)]
    pub fn silent(&mut self) -> &mut Self {
        self.silent = Some(true);
        self
    }

    /// Do not make backups. Overwrite all output files and directories.
    #[inline(always)]
    pub fn overwrite(&mut self) -> &mut Self {
        self.overwrite = Some(true);
        self
    }

    /// Alias for `output_yaml`.
    #[inline(always)]
    pub fn output(&mut self, value: &str) -> &mut Self {
        self.output_yaml(value)
    }

    /// Alias for `map`.
    #[inline(always)]
    pub fn maps(&mut self, value: OrderMap) -> &mut Self {
        self.map(value)
    }

    /// Alias for `map`.
    #[inline(always)]
    pub fn ordermap(&mut self, value: OrderMap) -> &mut Self {
        self.map(value)
    }

    /// Alias for `map`.
    #[inline(always)]
    pub fn ordermaps(&mut self, value: OrderMap) -> &mut Self {
        self.map(value)
    }

    /// Validate the process of analysis building.
    fn validate(&self) -> Result<(), String> {
        // check that step, min_samples and n_threads are not zero
        if let Some(step) = self.step {
            validate_step(step).map_err(|e| e.to_string())?;
        }

        if let Some(min_samples) = self.min_samples {
            validate_min_samples(min_samples).map_err(|e| e.to_string())?;
        }

        if let Some(n_threads) = self.n_threads {
            validate_n_threads(n_threads).map_err(|e| e.to_string())?;
        }

        // check that start is not larger than end
        if let (Some(begin), Some(end)) = (self.begin, self.end) {
            validate_begin_end(begin, end).map_err(|e| e.to_string())?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests_yaml {
    use approx::assert_relative_eq;

    use crate::{
        errors::OrderMapConfigError,
        input::{ordermap::Plane, GridSpan},
    };

    use super::*;

    #[test]
    fn analysis_yaml_pass_basic() {
        let analysis = Analysis::from_file("tests/files/inputs/basic.yaml").unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert!(analysis.index().is_none());
        assert_eq!(analysis.output(), "order.yaml");
        assert!(analysis.output_tab().is_none());
        assert!(analysis.output_xvg().is_none());
        assert!(analysis.output_csv().is_none());
        assert_eq!(analysis.membrane_normal(), Axis::Z);
        assert_eq!(
            analysis.heavy_atoms().unwrap(),
            "@membrane and element name carbon"
        );
        assert_eq!(
            analysis.hydrogens().unwrap(),
            "@membrane and element name hydrogen"
        );
        assert!(analysis.atoms().is_none());
        assert_eq!(analysis.begin(), 0.0);
        assert!(analysis.end().is_infinite());
        assert_eq!(analysis.step(), 1);
        assert_eq!(analysis.min_samples(), 1);
        assert_eq!(analysis.n_threads(), 1);
        assert!(analysis.leaflets().is_none());
        assert!(analysis.map().is_none());
        assert!(!analysis.silent());
        assert!(!analysis.overwrite());
    }

    #[test]
    fn analysis_yaml_pass_full() {
        let analysis = Analysis::from_file("tests/files/inputs/full.yaml").unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert_eq!(analysis.index().as_ref().unwrap(), "index.ndx");
        assert_eq!(analysis.output(), "order.yaml");
        assert_eq!(analysis.output_tab().as_ref().unwrap(), "order.dat");
        assert_eq!(analysis.output_xvg().as_ref().unwrap(), "order.xvg");
        assert_eq!(analysis.output_csv().as_ref().unwrap(), "order.csv");
        assert_eq!(analysis.membrane_normal(), Axis::X);
        assert!(analysis.heavy_atoms().is_none());
        assert!(analysis.hydrogens().is_none());
        assert_eq!(analysis.atoms().unwrap(), "@membrane");
        assert_eq!(analysis.begin(), 100.0);
        assert_eq!(analysis.end(), 10_000.0);
        assert_eq!(analysis.step(), 5);
        assert_eq!(analysis.min_samples(), 10);
        assert_eq!(analysis.n_threads(), 4);
        assert!(analysis.leaflets().is_some());

        let map = analysis.map().as_ref().unwrap();

        assert_eq!(map.output_directory(), ".");
        matches!(
            map.dim_x(),
            GridSpan::Manual {
                start: 0.5,
                end: 10.5
            }
        );
        matches!(map.dim_y(), GridSpan::Auto);
        assert_eq!(map.min_samples(), 10);
        assert_eq!(map.bin_size_x(), 0.05);
        assert_eq!(map.bin_size_y(), 0.02);
        assert_eq!(map.plane().unwrap(), Plane::XY);
    }

    #[test]
    fn analysis_yaml_pass_default_ordermap() {
        let analysis = Analysis::from_file("tests/files/inputs/default_ordermap.yaml").unwrap();

        let ordermap = analysis.map.unwrap();
        assert_eq!(ordermap.output_directory(), "ordermaps");
        assert_relative_eq!(ordermap.bin_size_x(), 0.1);
        assert_relative_eq!(ordermap.bin_size_y(), 0.1);
        assert_eq!(ordermap.min_samples(), 1);
        matches!(ordermap.dim_x(), GridSpan::Auto);
        matches!(ordermap.dim_y(), GridSpan::Auto);
    }

    #[test]
    fn analysis_yaml_fail_incomplete() {
        match Analysis::from_file("tests/files/inputs/incomplete.yaml") {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(ConfigError::CouldNotParseConfig(_, e)) => {
                assert!(e.to_string().contains("missing field"))
            }
            Err(e) => panic!("Unexpected error type returned: {}", e),
        }
    }

    fn check_analysis_error(file: &str, expected_variant: &str) {
        match Analysis::from_file(file) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(e) => match e {
                ConfigError::InvalidStep if expected_variant == "InvalidStep" => (),
                ConfigError::InvalidMinSamples if expected_variant == "InvalidMinSamples" => (),
                ConfigError::InvalidNThreads if expected_variant == "InvalidNThreads" => (),
                ConfigError::InvalidBeginEnd if expected_variant == "InvalidBeginEnd" => (),
                _ => panic!("Unexpected error type returned: {}", e),
            },
        }
    }

    #[test]
    fn analysis_yaml_fail_zero_step() {
        check_analysis_error("tests/files/inputs/zero_step.yaml", "InvalidStep");
    }

    #[test]
    fn analysis_yaml_fail_zero_min_samples() {
        check_analysis_error(
            "tests/files/inputs/zero_min_samples.yaml",
            "InvalidMinSamples",
        );
    }

    #[test]
    fn analysis_yaml_fail_zero_n_threads() {
        check_analysis_error("tests/files/inputs/zero_n_threads.yaml", "InvalidNThreads");
    }

    #[test]
    fn analysis_yaml_fail_start_higher_than_end() {
        check_analysis_error("tests/files/inputs/begin_higher.yaml", "InvalidBeginEnd");
    }

    #[test]
    fn analysis_yaml_fail_ordermap_zero_min_samples() {
        match Analysis::from_file("tests/files/inputs/ordermap_zero_min_samples.yaml") {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(ConfigError::InvalidOrderMap(OrderMapConfigError::InvalidMinSamples)) => (),
            Err(e) => panic!("Unexpected error type returned: {}", e),
        }
    }

    #[test]
    fn analysis_yaml_fail_ordermap_invalid_dim_x() {
        match Analysis::from_file("tests/files/inputs/ordermap_invalid_dim_x.yaml") {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(ConfigError::InvalidOrderMap(OrderMapConfigError::InvalidGridSpan(x, y))) => {
                assert_relative_eq!(x, 10.0);
                assert_relative_eq!(y, 7.5);
            }
            Err(e) => panic!("Unexpected error type returned: {}", e),
        }
    }

    #[test]
    fn analysis_yaml_fail_ordermap_invalid_dim_y() {
        match Analysis::from_file("tests/files/inputs/ordermap_invalid_dim_y.yaml") {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(ConfigError::InvalidOrderMap(OrderMapConfigError::InvalidGridSpan(x, y))) => {
                assert_relative_eq!(x, 10.0);
                assert_relative_eq!(y, 7.5);
            }
            Err(e) => panic!("Unexpected error type returned: {}", e),
        }
    }

    #[test]
    fn analysis_yaml_fail_ordermap_invalid_bin_size_x() {
        match Analysis::from_file("tests/files/inputs/ordermap_invalid_bin_size_x.yaml") {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(ConfigError::InvalidOrderMap(OrderMapConfigError::InvalidBinSize(x))) => {
                assert_relative_eq!(x, 0.0);
            }
            Err(e) => panic!("Unexpected error type returned: {}", e),
        }
    }

    #[test]
    fn analysis_yaml_fail_ordermap_invalid_bin_size_y() {
        match Analysis::from_file("tests/files/inputs/ordermap_invalid_bin_size_y.yaml") {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(ConfigError::InvalidOrderMap(OrderMapConfigError::InvalidBinSize(x))) => {
                assert_relative_eq!(x, -0.1);
            }
            Err(e) => panic!("Unexpected error type returned: {}", e),
        }
    }
}

#[cfg(test)]
mod tests_builder {

    use approx::assert_relative_eq;

    use crate::input::ordermap::Plane;

    use super::super::GridSpan;
    use super::*;

    #[test]
    fn analysis_builder_pass_basic() {
        let analysis = Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .build()
            .unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert!(analysis.index().is_none());
        assert_eq!(analysis.output(), "order.yaml");
        assert!(analysis.output_tab().is_none());
        assert!(analysis.output_xvg().is_none());
        assert!(analysis.output_csv().is_none());
        assert_eq!(analysis.membrane_normal(), Axis::Z);
        assert_eq!(
            analysis.heavy_atoms().unwrap(),
            "@membrane and element name carbon"
        );
        assert_eq!(
            analysis.hydrogens().unwrap(),
            "@membrane and element name hydrogen"
        );
        assert!(analysis.atoms().is_none());
        assert_eq!(analysis.begin(), 0.0);
        assert!(analysis.end().is_infinite());
        assert_eq!(analysis.step(), 1);
        assert_eq!(analysis.min_samples(), 1);
        assert_eq!(analysis.n_threads(), 1);
        assert!(analysis.leaflets().is_none());
        assert!(analysis.map().is_none());
        assert!(!analysis.silent());
        assert!(!analysis.overwrite());
    }

    #[test]
    fn analysis_builder_pass_full() {
        let analysis = Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .index("index.ndx")
            .output("order.yaml")
            .output_tab("order.dat")
            .output_xvg("order.xvg")
            .output_csv("order.csv")
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .membrane_normal(Axis::X)
            .begin(100.0)
            .end(10_000.0)
            .step(5)
            .min_samples(10)
            .n_threads(4)
            .leaflets(LeafletClassification::global("@membrane", "name P"))
            .map(
                OrderMap::new()
                    .output_directory(".")
                    .dim([
                        GridSpan::Manual {
                            start: 0.5,
                            end: 10.5,
                        },
                        GridSpan::Auto,
                    ])
                    .min_samples(10)
                    .bin_size([0.05, 0.02])
                    .plane(Plane::XY)
                    .build()
                    .unwrap(),
            )
            .build()
            .unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert_eq!(analysis.index().as_ref().unwrap(), "index.ndx");
        assert_eq!(analysis.output(), "order.yaml");
        assert_eq!(analysis.output_tab().as_ref().unwrap(), "order.dat");
        assert_eq!(analysis.output_xvg().as_ref().unwrap(), "order.xvg");
        assert_eq!(analysis.output_csv().as_ref().unwrap(), "order.csv");
        assert_eq!(analysis.membrane_normal(), Axis::X);
        assert!(analysis.heavy_atoms().is_none());
        assert!(analysis.hydrogens().is_none());
        assert_eq!(analysis.atoms().unwrap(), "@membrane");
        assert_eq!(analysis.begin(), 100.0);
        assert_eq!(analysis.end(), 10_000.0);
        assert_eq!(analysis.step(), 5);
        assert_eq!(analysis.min_samples(), 10);
        assert_eq!(analysis.n_threads(), 4);
        assert!(analysis.leaflets().is_some());

        let map = analysis.map().as_ref().unwrap();

        assert_eq!(map.output_directory(), ".");
        matches!(
            map.dim_x(),
            GridSpan::Manual {
                start: 0.5,
                end: 10.5
            }
        );
        matches!(map.dim_y(), GridSpan::Auto);
        assert_eq!(map.min_samples(), 10);
        assert_eq!(map.bin_size_x(), 0.05);
        assert_eq!(map.bin_size_y(), 0.02);
        assert_eq!(map.plane().unwrap(), Plane::XY);
    }

    #[test]
    fn analysis_builder_pass_default_ordermap() {
        let analysis = Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .map(
                OrderMap::new()
                    .output_directory("ordermaps")
                    .build()
                    .unwrap(),
            )
            .build()
            .unwrap();

        let ordermap = analysis.map.unwrap();
        assert_eq!(ordermap.output_directory(), "ordermaps");
        assert_relative_eq!(ordermap.bin_size_x(), 0.1);
        assert_relative_eq!(ordermap.bin_size_y(), 0.1);
        assert_eq!(ordermap.min_samples(), 1);
        matches!(ordermap.dim_x(), GridSpan::Auto);
        matches!(ordermap.dim_y(), GridSpan::Auto);
    }

    #[test]
    fn analysis_builder_fail_incomplete() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::UninitializedField(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_zero_step() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .step(0)
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_zero_min_samples() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .min_samples(0)
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_zero_threads() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .n_threads(0)
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_start_higher_than_end() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .begin(100_000.0)
            .end(50_000.0)
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }
}
