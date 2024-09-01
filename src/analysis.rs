// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the main `Analysis` structure and its methods.

use colored::Colorize;
use derive_builder::Builder;
use getset::{CopyGetters, Getters};

use super::axis::Axis;
use super::leaflets::LeafletClassification;
use super::ordermap::OrderMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisType {
    AAOrder,
    CGOrder,
}

/// Structure holding all the information necessary to perform the specified analysis.
#[derive(Debug, Clone, Builder, Getters, CopyGetters)]
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
    #[getset(get_copy = "pub")]
    analysis_type: AnalysisType,
    /// Direction of the membrane normal.
    /// If not provided, the default value is 'Axis::Z'.
    #[builder(setter(into), default = "Axis::Z")]
    #[getset(get_copy = "pub")]
    membrane_normal: Axis,
    /// Selection of heavy atoms for which the order parameters should be calculated.
    /// Only needed if analysis_type is 'AAOrder'.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    heavy_atoms: Option<String>,
    /// Selection of hydrogens to be used for the order parameters calculation.
    /// Only hydrogens from this selection that are bonded to the `heavy_atoms` will be used.
    /// Only needed if analysis_type is 'AAOrder'.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    hydrogens: Option<String>,
    /// Atoms/beads forming bonds for which the coarse-grained order parameters should be calculated.
    /// The bonds will be automatically searched for between the selected atoms.
    /// Only needed if analysis_type is 'CGOrder'.
    #[builder(setter(into, strip_option), default)]
    #[getset(get = "pub")]
    atoms: Option<String>,
    /// Starting time of the trajectory analysis (in ps).
    /// If not specified, the analysis starts at the beginning of the trajectory.
    #[builder(default = "0.0")]
    #[getset(get_copy = "pub")]
    start: f32,
    /// Ending time of the trajectory analysis (in ps).
    /// If not specified, the analysis ends at the end of the trajectory.
    #[builder(default = "f32::INFINITY")]
    #[getset(get_copy = "pub")]
    end: f32,
    /// Only every Nth frame of the simulation trajectory will be analyzed.
    /// If not specified, each frame of the trajectory will be analyzed.
    #[builder(default = "1")]
    #[getset(get_copy = "pub")]
    step: usize,
    /// Minimal number of samples for each heavy atom required to calculate order parameter for it.
    /// If not specified, the default value is 1.
    #[builder(default = "1")]
    #[getset(get_copy = "pub")]
    min_samples: usize,
    /// Number of threads to use to perform the analysis.
    /// If not specified, the default value is 1.
    #[builder(default = "1")]
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
}

impl Analysis {
    pub fn new() -> AnalysisBuilder {
        AnalysisBuilder::default()
    }

    /// Perform the analysis and write out the results.
    pub fn analyze(&self) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        match self.analysis_type() {
            AnalysisType::AAOrder => crate::aaorder::analyze_atomistic(self),
            AnalysisType::CGOrder => todo!("Implement CG order calculation."),
        }
    }
}

impl AnalysisBuilder {
    /// Validate the process of analysis building.
    fn validate(&self) -> Result<(), String> {
        match self.analysis_type {
            // skip validation, building will fail anyway
            None => Ok(()),
            Some(AnalysisType::AAOrder) => {
                // `heavy_atoms` and `hydrogens` must be defined
                // `atoms` cannot be defined
                if self.heavy_atoms.is_none() {
                    let error = format!(
                        "{} heavy atoms were not specified and the analysis type is '{}'.",
                        "error:".red().bold(),
                        "aaorder".yellow()
                    );

                    Err(error)
                } else if self.hydrogens.is_none() {
                    let error = format!(
                        "{} hydrogens were not specified and the analysis type is '{}'.",
                        "error:".red().bold(),
                        "aaorder".yellow()
                    );

                    Err(error)
                } else if !self.atoms.is_none() {
                    let error = format!(
                        "{} atoms were specified but the analysis type is '{}' not '{}'.",
                        "error:".red().bold(),
                        "aaorder".yellow(),
                        "cgorder".yellow()
                    );

                    Err(error)
                } else {
                    Ok(())
                }
            }
            Some(AnalysisType::CGOrder) => {
                // `atoms` must be defined
                // `heavy_atoms` and `hydrogens` cannot be defined
                if self.atoms.is_none() {
                    let error = format!(
                        "{} atoms were not specified and the analysis type is '{}'.",
                        "error:".red().bold(),
                        "cgorder".yellow()
                    );

                    Err(error)
                } else if !self.heavy_atoms.is_none() {
                    let error = format!(
                        "{} heavy atoms were specified but the analysis type is '{}' not '{}'.",
                        "error:".red().bold(),
                        "cgorder".yellow(),
                        "aaorder".yellow()
                    );

                    Err(error)
                } else if !self.hydrogens.is_none() {
                    let error = format!(
                        "{} hydrogens were specified but the analysis type is '{}' not '{}'.",
                        "error:".red().bold(),
                        "cgorder".yellow(),
                        "aaorder".yellow()
                    );

                    Err(error)
                } else {
                    Ok(())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests_builder {

    use super::*;
    use crate::ordermap::GridSpan;

    #[test]
    fn analysis_builder_pass_basic() {
        let analysis = Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::AAOrder)
            .heavy_atoms("@membrane and element name carbon")
            .hydrogens("@membrane and element name hydrogen")
            .build()
            .unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert!(analysis.index().is_none());
        assert_eq!(analysis.output(), "order.yaml");
        assert_eq!(analysis.analysis_type(), AnalysisType::AAOrder);
        assert_eq!(analysis.membrane_normal(), Axis::Z);
        assert_eq!(
            analysis.heavy_atoms().as_ref().unwrap(),
            "@membrane and element name carbon"
        );
        assert_eq!(
            analysis.hydrogens().as_ref().unwrap(),
            "@membrane and element name hydrogen"
        );
        assert!(analysis.atoms().is_none());
        assert_eq!(analysis.start(), 0.0);
        assert!(analysis.end().is_infinite());
        assert_eq!(analysis.step(), 1);
        assert_eq!(analysis.min_samples(), 1);
        assert_eq!(analysis.n_threads(), 1);
        assert!(analysis.leaflets().is_none());
        assert!(analysis.map().is_none());
    }

    #[test]
    fn analysis_builder_pass_full() {
        let analysis = Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .index("index.ndx")
            .output("order.yaml")
            .analysis_type(AnalysisType::CGOrder)
            .atoms("@membrane")
            .membrane_normal(Axis::X)
            .start(100.0)
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
                    .bin_size(0.2)
                    .build()
                    .unwrap(),
            )
            .build()
            .unwrap();

        assert_eq!(analysis.structure(), "system.tpr");
        assert_eq!(analysis.trajectory(), "md.xtc");
        assert_eq!(analysis.index().as_ref().unwrap(), "index.ndx");
        assert_eq!(analysis.output(), "order.yaml");
        assert_eq!(analysis.analysis_type(), AnalysisType::CGOrder);
        assert_eq!(analysis.membrane_normal(), Axis::X);
        assert!(analysis.heavy_atoms().is_none());
        assert!(analysis.hydrogens().is_none());
        assert_eq!(analysis.atoms().as_ref().unwrap(), "@membrane");
        assert_eq!(analysis.start(), 100.0);
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
        assert_eq!(map.bin_size(), 0.2);
    }

    #[test]
    fn analysis_builder_fail_incomplete() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .atoms("@membrane")
            .output("order.yaml")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::UninitializedField(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_aaorder_missing_heavy() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::AAOrder)
            .hydrogens("@membrane and element name hydrogen")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_aaorder_missing_hydrogens() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .output("order.yaml")
            .analysis_type(AnalysisType::AAOrder)
            .heavy_atoms("@membrane and element name carbon")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_aaorder_atoms_present() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .analysis_type(AnalysisType::AAOrder)
            .heavy_atoms("@membrane and element name carbon")
            .hydrogens("@membrane and element name hydrogen")
            .atoms("@membrane")
            .output("order.yaml")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_cgorder_missing_atoms() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .analysis_type(AnalysisType::CGOrder)
            .hydrogens("@membrane and element name hydrogen")
            .output("order.yaml")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_cgorder_hydrogens_present() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .analysis_type(AnalysisType::CGOrder)
            .atoms("@membrane")
            .hydrogens("@membrane and element name hydrogen")
            .output("order.yaml")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }

    #[test]
    fn analysis_builder_fail_cgorder_heavy_atoms_present() {
        match Analysis::new()
            .structure("system.tpr")
            .trajectory("md.xtc")
            .analysis_type(AnalysisType::CGOrder)
            .atoms("@membrane")
            .heavy_atoms("@membrane and element name carbon")
            .output("order.yaml")
            .build()
        {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(AnalysisBuilderError::ValidationError(_)) => (),
            Err(_) => panic!("Incorrect error type returned."),
        }
    }
}
