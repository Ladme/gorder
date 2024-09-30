// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the main `Analysis` structure and its methods.

use colored::Colorize;
use derive_builder::Builder;
use getset::Setters;
use getset::{CopyGetters, Getters};
use serde::Deserialize;
use serde_valid::Validate;

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
        atoms: String,
    },
}

impl AnalysisType {
    pub fn aaorder(heavy_atoms: &str, hydrogens: &str) -> Self {
        Self::AAOrder {
            heavy_atoms: heavy_atoms.to_owned(),
            hydrogens: hydrogens.to_owned(),
        }
    }

    pub fn cgorder(atoms: &str) -> Self {
        Self::CGOrder {
            atoms: atoms.to_owned(),
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::AAOrder {
                heavy_atoms: _,
                hydrogens: _,
            } => "all-atom order parameters",
            Self::CGOrder { atoms: _ } => "coarse-grained order parameters",
        }
    }
}

/// Structure holding all the information necessary to perform the specified analysis.
#[derive(Debug, Clone, Builder, Getters, CopyGetters, Setters, Deserialize, Validate)]
#[serde(deny_unknown_fields)]
#[validate(custom = |x| validate_begin_end(x.begin, x.end))]
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
    output: String,
    /// Type of the analysis to perform (AAOrder / CGOrder).
    #[getset(get = "pub")]
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
    #[serde(default = "default_begin")]
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
    #[validate(custom = validate_step)]
    #[getset(get_copy = "pub")]
    step: usize,
    /// Minimal number of samples for each heavy atom required to calculate order parameter for it.
    /// If not specified, the default value is 1.
    #[builder(default = "1")]
    #[serde(default = "default_one")]
    #[validate(custom = validate_min_samples)]
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
    #[getset(get = "pub")]
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

fn validate_step(step: &usize) -> Result<(), serde_valid::validation::Error> {
    if *step == 0 {
        let error = format!(
            "{} the specified value of '{}' is invalid (must be positive).",
            "error:".red().bold(),
            "step".yellow()
        );

        Err(serde_valid::validation::Error::Custom(error))
    } else {
        Ok(())
    }
}

fn validate_min_samples(samples: &usize) -> Result<(), serde_valid::validation::Error> {
    if *samples == 0 {
        let error = format!(
            "{} the specified value of '{}' is invalid (must be positive).",
            "error:".red().bold(),
            "min_samples".yellow()
        );

        Err(serde_valid::validation::Error::Custom(error))
    } else {
        Ok(())
    }
}

fn validate_n_threads(n_threads: &usize) -> Result<(), serde_valid::validation::Error> {
    if *n_threads == 0 {
        let error = format!(
            "{} the specified value of '{}' is invalid (must be positive).",
            "error:".red().bold(),
            "n_threads".yellow()
        );

        Err(serde_valid::validation::Error::Custom(error))
    } else {
        Ok(())
    }
}

fn validate_begin_end(begin: f32, end: f32) -> Result<(), serde_valid::validation::Error> {
    if begin > end {
        let error = format!(
            "{} invalid values of '{}' and '{}' (start is higher than end).",
            "error:".red().bold(),
            "start".yellow(),
            "end".yellow()
        );

        Err(serde_valid::validation::Error::Custom(error))
    } else {
        Ok(())
    }
}

impl Analysis {
    pub fn new() -> AnalysisBuilder {
        AnalysisBuilder::default()
    }

    pub fn heavy_atoms(&self) -> Option<&String> {
        match &self.analysis_type {
            AnalysisType::CGOrder { atoms: _ } => None,
            AnalysisType::AAOrder {
                heavy_atoms,
                hydrogens: _,
            } => Some(heavy_atoms),
        }
    }

    pub fn hydrogens(&self) -> Option<&String> {
        match &self.analysis_type {
            AnalysisType::CGOrder { atoms: _ } => None,
            AnalysisType::AAOrder {
                heavy_atoms: _,
                hydrogens,
            } => Some(hydrogens),
        }
    }

    pub fn atoms(&self) -> Option<&String> {
        match &self.analysis_type {
            AnalysisType::CGOrder { atoms } => Some(atoms),
            AnalysisType::AAOrder {
                heavy_atoms: _,
                hydrogens: _,
            } => None,
        }
    }
}

impl AnalysisBuilder {
    /// Be silent. Print nothing to the standard output during the analysis.
    pub fn silent(&mut self) -> &mut Self {
        self.silent = Some(true);
        self
    }

    /// Do not make backups. Overwrite all output files and directories.
    pub fn overwrite(&mut self) -> &mut Self {
        self.overwrite = Some(true);
        self
    }

    /// Validate the process of analysis building.
    fn validate(&self) -> Result<(), String> {
        // check that step, min_samples and n_threads are not zero
        if let Some(ref step) = self.step {
            validate_step(step).map_err(|x| x.to_string())?;
        }

        if let Some(ref min_samples) = self.min_samples {
            validate_min_samples(min_samples).map_err(|x| x.to_string())?;
        }

        if let Some(ref n_threads) = self.n_threads {
            validate_n_threads(n_threads).map_err(|x| x.to_string())?;
        }

        // check that start is not larger than end
        if let (Some(begin), Some(end)) = (self.begin, self.end) {
            validate_begin_end(begin, end).map_err(|x| x.to_string())?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests_builder {

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
                    .dim_x(GridSpan::Manual {
                        start: 0.5,
                        end: 10.5,
                    })
                    .dim_y(GridSpan::Auto)
                    .min_samples(10)
                    .bin_size_x(0.05)
                    .bin_size_y(0.02)
                    .build()
                    .unwrap(),
            )
            .build()
            .unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert_eq!(analysis.index().as_ref().unwrap(), "index.ndx");
        assert_eq!(analysis.output(), "order.yaml");
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
