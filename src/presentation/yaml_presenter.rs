// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for writing yaml files.

use crate::errors::WriteError;
use crate::presentation::aaresults::AAOrderResults;
use crate::presentation::cgresults::CGOrderResults;
use crate::presentation::{OrderResults, OutputFormat, Presenter, PresenterProperties};
use crate::PANIC_MESSAGE;
use std::io::Write;

/// Structure handling the writing of a yaml output.
#[derive(Debug, Clone)]
pub(super) struct YamlPresenter<'a, R: OrderResults> {
    /// Results of the analysis.
    results: &'a R,
    /// Parameters necessary for the writing of the yaml file.
    properties: YamlProperties,
}

/// Structure containing parameters necessary for the writing of the yaml file.
#[derive(Debug, Clone)]
pub(crate) struct YamlProperties {
    /// Name of the input structure file.
    structure: String,
    /// Name of the input trajectory file.
    trajectory: String,
}

/// Trait implemented by all structures that can be written in yaml format.
pub(crate) trait YamlWrite {
    /// Write the structure in a yaml format into an output file.
    fn write_yaml(
        &self,
        writer: &mut impl Write,
        properties: &YamlProperties,
    ) -> Result<(), WriteError>;
}

impl PresenterProperties for YamlProperties {
    #[allow(unused)]
    fn leaflets(&self) -> bool {
        panic!("FATAL GORDER ERROR | YamlProperties::leaflets | This method should never be called. {}", PANIC_MESSAGE);
    }
}

impl YamlProperties {
    /// Create new structure capturing properties of the yaml format.
    pub(super) fn new(structure: &str, trajectory: &str) -> Self {
        Self {
            structure: structure.to_owned(),
            trajectory: trajectory.to_owned(),
        }
    }
}

impl<'a, R: OrderResults> Presenter<'a, R> for YamlPresenter<'a, R> {
    type Properties = YamlProperties;

    fn new(results: &'a R, properties: YamlProperties) -> Self {
        Self {
            results,
            properties,
        }
    }

    fn file_format(&self) -> OutputFormat {
        OutputFormat::YAML
    }

    fn write_results(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        self.results.write_yaml(writer, &self.properties)
    }

    #[allow(unused)]
    fn write_empty_order(
        writer: &mut impl Write,
        properties: &YamlProperties,
    ) -> Result<(), WriteError> {
        panic!("FATAL GORDER ERROR | YamlPresenter::write_empty_order | This method should never be called. {}", PANIC_MESSAGE);
    }
}

impl YamlWrite for AAOrderResults {
    /// Write yaml data for atomistic order parameters.
    fn write_yaml(
        &self,
        writer: &mut impl Write,
        properties: &YamlProperties,
    ) -> Result<(), WriteError> {
        YamlPresenter::<AAOrderResults>::write_header(
            writer,
            &properties.structure,
            &properties.trajectory,
        )?;

        serde_yaml::to_writer(writer, self).map_err(|e| WriteError::CouldNotWriteYaml(e))
    }
}

impl YamlWrite for CGOrderResults {
    /// Write yaml data for coarse-grained order parameters.
    fn write_yaml(
        &self,
        writer: &mut impl Write,
        properties: &YamlProperties,
    ) -> Result<(), WriteError> {
        YamlPresenter::<CGOrderResults>::write_header(
            writer,
            &properties.structure,
            &properties.trajectory,
        )?;

        serde_yaml::to_writer(writer, self).map_err(|e| WriteError::CouldNotWriteYaml(e))
    }
}