// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for presenting the results of the analysis.

use crate::input::Plane;
use crate::presentation::aaresults::AAOrderResults;
use crate::presentation::cgresults::CGOrderResults;
use crate::presentation::csv_presenter::{CsvPresenter, CsvProperties, CsvWrite};
use crate::presentation::ordermaps_presenter::MapWrite;
use crate::presentation::ordermaps_presenter::{OrderMapPresenter, OrderMapProperties};
use crate::presentation::tab_presenter::{TabPresenter, TabProperties, TabWrite};
use crate::presentation::xvg_presenter::{XvgPresenter, XvgProperties, XvgWrite};
use crate::presentation::yaml_presenter::{YamlPresenter, YamlProperties, YamlWrite};
use crate::{
    analysis::{
        molecule::{AtomType, BondTopology, MoleculeType},
        topology::SystemTopology,
    },
    errors::WriteError,
    input::Analysis,
    PANIC_MESSAGE,
};
use getset::{CopyGetters, Getters};
use groan_rs::prelude::GridMap;
use indexmap::IndexMap;
use serde::ser::SerializeMap;
use serde::{Serialize, Serializer};
use std::fmt::Debug;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};
use strum_macros::Display;

macro_rules! write_result {
    ($dst:expr, $($arg:tt)*) => {
        write!($dst, $($arg)*).map_err(|e| WriteError::CouldNotWriteResults(e))?
    };
}

pub mod aaresults;
pub mod cgresults;
pub(crate) mod converter;
mod csv_presenter;
//pub(crate) mod ordermap;
pub mod ordermaps_presenter;
mod tab_presenter;
mod xvg_presenter;
mod yaml_presenter;

/// Enum representing any of the types of results that can be returned by the `gorder`.
#[derive(Debug, Clone)]
pub enum AnalysisResults {
    AA(AAOrderResults),
    CG(CGOrderResults),
}

impl AnalysisResults {
    /// Write the results of the analysis into output files.
    pub fn write(&self) -> Result<(), WriteError> {
        match self {
            AnalysisResults::AA(x) => x.write_all_results(),
            AnalysisResults::CG(x) => x.write_all_results(),
        }
    }
}

/// Type alias for a gridmap of f32 values.
pub type GridMapF32 = GridMap<f32, f32, fn(&f32) -> f32>;

/// Public trait implemented by results-containing structures.
pub trait PublicOrderResults {
    #[allow(private_bounds)]
    type MoleculeResults: MoleculeResults;

    /// Get the results for all molecules.
    fn molecules(&self) -> impl Iterator<Item = &Self::MoleculeResults>;

    /// Get the results for a molecule with the specified name.
    /// O(1) complexity.
    /// Returns `None` if such molecule does not exist.
    fn get_molecule(&self, name: &str) -> Option<&Self::MoleculeResults>;

    /// Get the parameters of the analysis.
    fn analysis(&self) -> &Analysis;
}

/// Trait implemented by all structures providing the full results of the analysis.
pub(crate) trait OrderResults:
    Debug + Clone + CsvWrite + TabWrite + YamlWrite + PublicOrderResults
{
    type OrderType: OrderType;

    /// Create an empty `OrderResults` structure (i.e., without any molecules).
    fn empty(analysis: Analysis) -> Self;

    /// Create a new `OrderResults` structure.
    fn new(molecules: IndexMap<String, Self::MoleculeResults>, analysis: Analysis) -> Self;

    /// Return the maximal number of bonds for heavy atoms in the system.
    /// Only makes sense for atomistic order. Returns 0 by default.
    #[inline(always)]
    fn max_bonds(&self) -> usize {
        0
    }

    /// Write results of the analysis into the output files.
    fn write_all_results(&self) -> Result<(), WriteError> {
        if self.molecules().count() == 0 {
            log::warn!("Nothing to write.");
            return Ok(());
        }

        let analysis = self.analysis();
        let errors = analysis.estimate_error().is_some();
        let leaflets = analysis.leaflets().is_some();
        let input_structure = analysis.structure();
        let input_trajectory = analysis.trajectory();
        let overwrite = analysis.overwrite();

        YamlPresenter::new(self, YamlProperties::new(input_structure, input_trajectory))
            .write(analysis.output_yaml(), overwrite)?;

        if let Some(tab) = analysis.output_tab() {
            TabPresenter::new(
                self,
                TabProperties::new(self, input_structure, input_trajectory, errors, leaflets),
            )
            .write(tab, overwrite)?;
        }

        if let Some(xvg) = analysis.output_xvg() {
            log::info!("Writing the order parameters into xvg file(s)...");

            XvgPresenter::new(
                self,
                XvgProperties::new(input_structure, input_trajectory, leaflets),
            )
            .write(xvg, overwrite)?;
        }

        if let Some(csv) = analysis.output_csv() {
            CsvPresenter::new(self, CsvProperties::new(self, errors, leaflets))
                .write(csv, overwrite)?;
        }

        if let Some(map) = analysis.map() {
            OrderMapPresenter::new(
                self,
                OrderMapProperties::new(map.plane().unwrap_or(Plane::XY).into()),
            )
            .write(map.output_directory(), overwrite)?;
        }

        Ok(())
    }
}

/// Trait implemented by all structures providing the results of the analysis for a single molecule type.
pub(crate) trait MoleculeResults:
    Debug + Clone + CsvWrite + TabWrite + XvgWrite + MapWrite
{
    /// Return the maximal number of bonds for heavy atoms in the molecule.
    /// Only makes sense for atomistic order. Returns 0 by default.
    #[inline(always)]
    fn max_bonds(&self) -> usize {
        0
    }

    /// Get the name of the molecule
    fn molecule(&self) -> &str;
}

/// All supported output file formats.
#[derive(Debug, Clone, Display)]
pub(crate) enum OutputFormat {
    #[strum(serialize = "yaml")]
    YAML,
    #[strum(serialize = "csv")]
    CSV,
    #[strum(serialize = "xvg")]
    XVG,
    #[strum(serialize = "tab")]
    TAB,
    #[strum(serialize = "map")]
    MAP,
}

/// Trait implemented by all structures that store properties of Presenters.
pub(crate) trait PresenterProperties: Debug + Clone {
    /// Is the data for leaflets available?
    fn leaflets(&self) -> bool;
}

/// Trait implemented by all structures presenting results of the analysis.
pub(crate) trait Presenter<'a, R: OrderResults>: Debug + Clone {
    /// Structure describing properties of the presenter.
    type Properties: PresenterProperties;

    /// Create a new presenter structure.
    fn new(results: &'a R, properties: Self::Properties) -> Self;

    /// Get the format of the output file that this Presenter creates.
    fn file_format(&self) -> OutputFormat;

    /// Write the results into an open output file.
    fn write_results(&self, writer: &mut impl Write) -> Result<(), WriteError>;

    /// Write empty (missing) order parameter into the output file.
    fn write_empty_order(
        writer: &mut impl Write,
        properties: &Self::Properties,
    ) -> Result<(), WriteError>;

    /// Write empty (missing) order parameters for an entire collection.
    fn write_empty_bond_collection(
        writer: &mut impl Write,
        properties: &Self::Properties,
    ) -> Result<(), WriteError> {
        if properties.leaflets() {
            for _ in 0..3 {
                Self::write_empty_order(writer, properties)?;
            }
        } else {
            Self::write_empty_order(writer, properties)?;
        }

        Ok(())
    }

    /// Create (and potentially back up) an output file, open it and write the results into it.
    fn write(&self, filename: impl AsRef<Path>, overwrite: bool) -> Result<(), WriteError> {
        log::info!(
            "Writing the order parameters into {} file '{}'...",
            self.file_format(),
            filename.as_ref().to_str().expect(PANIC_MESSAGE)
        );
        log::logger().flush();

        self.try_backup(&filename, overwrite)?;
        let mut writer = Self::create_and_open(&filename)?;
        self.write_results(&mut writer)
    }

    /// Create and open a file for buffered writing.
    #[inline(always)]
    fn create_and_open(filename: &impl AsRef<Path>) -> Result<BufWriter<File>, WriteError> {
        let file = File::create(filename.as_ref())
            .map_err(|_| WriteError::CouldNotCreateFile(Box::from(filename.as_ref())))?;

        Ok(BufWriter::new(file))
    }

    /// Back up an output file, if it is necessary and if it is requested.
    fn try_backup(&self, filename: &impl AsRef<Path>, overwrite: bool) -> Result<(), WriteError> {
        if filename.as_ref().exists() {
            if !overwrite {
                log::warn!(
                    "Output {} file '{}' already exists. Backing it up.",
                    self.file_format(),
                    filename.as_ref().to_str().expect(PANIC_MESSAGE)
                );
                backitup::backup(filename.as_ref())
                    .map_err(|_| WriteError::CouldNotBackupFile(Box::from(filename.as_ref())))?;
            } else {
                log::warn!(
                    "Output {} file '{}' already exists. It will be overwritten as requested.",
                    self.file_format(),
                    filename.as_ref().to_str().expect(PANIC_MESSAGE)
                );
            }
        }

        Ok(())
    }

    /// Write header into the output file specifying version of the library used and some other basic info.
    fn write_header(
        writer: &mut impl Write,
        structure: &str,
        trajectory: &str,
    ) -> Result<(), WriteError> {
        write_result!(writer, "# Order parameters calculated with 'gorder v{}' using structure file '{}' and trajectory file '{}'.\n",
        crate::GORDER_VERSION, structure, trajectory);

        Ok(())
    }
}

/// Single order parameter value, optionally with its estimated error.
#[derive(Debug, Clone, Copy, CopyGetters)]
pub struct Order {
    /// Value of the order parameter (mean from the analyzed frames).
    #[getset(get_copy = "pub")]
    value: f32,
    /// Estimated error for this order parameter (standard deviation of N blocks).
    #[getset(get_copy = "pub")]
    error: Option<f32>,
}

impl From<f32> for Order {
    #[inline(always)]
    fn from(value: f32) -> Self {
        Order { value, error: None }
    }
}

impl From<[f32; 2]> for Order {
    #[inline(always)]
    fn from(value: [f32; 2]) -> Self {
        Order {
            value: value[0],
            error: Some(value[1]),
        }
    }
}

impl Serialize for Order {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        if let Some(error) = self.error {
            // serialize as a dictionary if `error` is present
            let mut map = serializer.serialize_map(Some(2))?;
            map.serialize_entry("mean", &self.value.round_to_4())?;
            map.serialize_entry("error", &error.round_to_4())?;
            map.end()
        } else {
            // serialize as a single float rounded to 4 decimal places
            serializer.serialize_f64(self.value.round_to_4())
        }
    }
}

// Helper trait for rounding a float to 4 decimal places
trait RoundTo4 {
    fn round_to_4(self) -> f64;
}

impl RoundTo4 for f32 {
    fn round_to_4(self) -> f64 {
        (self as f64 * 10_000.0).round() / 10_000.0
    }
}

/// Collection of (up to) 3 order parameters: for the full membrane, the upper leaflet,
/// and the lower leaflet.
#[derive(Debug, Clone, Default, Serialize, Getters)]
pub struct OrderCollection {
    /// Order parameter for the full membrane.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[getset(get = "pub")]
    total: Option<Order>,
    /// Order parameter for the upper leaflet.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[getset(get = "pub")]
    upper: Option<Order>,
    /// Order parameter for the lower leaflet.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[getset(get = "pub")]
    lower: Option<Order>,
}

impl OrderCollection {
    fn new(total: Option<Order>, upper: Option<Order>, lower: Option<Order>) -> Self {
        Self {
            total,
            upper,
            lower,
        }
    }
}

/// Collection of (up to) 3 order maps: for the full membrane, the upper leaflet,
/// and the lower leaflet.
#[derive(Debug, Clone, Default, Getters)]
pub struct OrderMapsCollection {
    #[getset(get = "pub")]
    total: Option<GridMapF32>,
    #[getset(get = "pub")]
    upper: Option<GridMapF32>,
    #[getset(get = "pub")]
    lower: Option<GridMapF32>,
}

impl OrderMapsCollection {
    pub(super) fn new(
        total: Option<GridMapF32>,
        upper: Option<GridMapF32>,
        lower: Option<GridMapF32>,
    ) -> Self {
        Self {
            total,
            upper,
            lower,
        }
    }
}

/// Order parameters calculated for a single bond.
#[derive(Debug, Clone, Serialize, Getters)]
pub struct BondResults {
    /// Name of the bond.
    #[serde(skip)]
    #[getset(get = "pub(super)")] // intentionally not public
    bond: BondTopology,
    /// Name of the molecule this bond belongs to.
    #[serde(skip)]
    #[getset(get = "pub")]
    molecule: String,
    /// Order parameters calculated for this bond.
    #[serde(flatten)]
    #[getset(get = "pub")]
    order: OrderCollection,
    /// Ordermaps calculated for this bond.
    #[serde(skip)]
    #[getset(get = "pub")]
    ordermaps: OrderMapsCollection,
}

impl BondResults {
    /// Atom types involved in the bond.
    pub fn atoms(&self) -> (&AtomType, &AtomType) {
        (self.bond.atom1(), self.bond.atom2())
    }
}

/// Empty struct used as a marker.
#[derive(Default, Debug, Clone)]
pub(crate) struct AAOrder {}
/// Empty struct used as a marker.
#[derive(Default, Debug, Clone)]
pub(crate) struct CGOrder {}

/// Trait implemented only by `AAOrder` and `CGOrder` structs.
pub(crate) trait OrderType: Debug + Clone {
    /// Used to convert an order parameter to its final value depending on the analysis type.
    /// Atomistic order parameters are reported as -S_CD.
    /// Coarse grained order parameters are reported as P2.
    fn convert(order: f32, error: Option<f32>) -> Order;

    /// String to use as a label for z-axis in the ordermap.
    fn zlabel() -> &'static str;
}

impl OrderType for AAOrder {
    #[inline(always)]
    fn convert(order: f32, error: Option<f32>) -> Order {
        Order {
            value: -order,
            error,
        }
    }

    #[inline(always)]
    fn zlabel() -> &'static str {
        "order parameter ($-S_{CH}$)"
    }
}

impl OrderType for CGOrder {
    #[inline(always)]
    fn convert(order: f32, error: Option<f32>) -> Order {
        Order {
            value: order,
            error,
        }
    }

    #[inline(always)]
    fn zlabel() -> &'static str {
        "order parameter ($S$)"
    }
}

impl BondResults {
    #[inline(always)]
    fn new(
        bond: &BondTopology,
        molecule: &str,
        order: OrderCollection,
        ordermaps: OrderMapsCollection,
    ) -> Self {
        Self {
            bond: bond.clone(),
            molecule: molecule.to_owned(),
            order,
            ordermaps,
        }
    }
}

impl Serialize for AtomType {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let formatted_string = format!(
            "{} {} ({})",
            self.residue_name(),
            self.atom_name(),
            self.relative_index()
        );
        serializer.serialize_str(&formatted_string)
    }
}

impl Analysis {
    /// Print basic information about the analysis for the user.
    pub(crate) fn info(&self) {
        log::info!("Will calculate {}.", self.analysis_type().name());
        log::info!(
            "Membrane normal expected to be oriented along the {} axis.",
            self.membrane_normal()
        );
        if self.map().is_some() {
            log::info!(
                "Will calculate ordermaps in the {} plane.",
                self.map().as_ref().unwrap().plane().expect(PANIC_MESSAGE)
            );
        }
        if self.leaflets().is_some() {
            log::info!(
                "Will classify lipids into membrane leaflets using the '{}' method.",
                self.leaflets().as_ref().expect(PANIC_MESSAGE)
            )
        }
    }
}

impl MoleculeType {
    /// Print basic information about the molecule type for the user.
    #[inline(always)]
    fn info(&self) {
        log::info!(
            "Molecule type {}: {} order bonds, {} molecules.",
            self.name(),
            self.order_bonds().bond_types().len(),
            self.order_bonds()
                .bond_types()
                .first()
                .expect(PANIC_MESSAGE)
                .bonds()
                .len()
        )
    }
}

impl SystemTopology {
    /// Print basic information about the system topology for the user.
    #[inline(always)]
    pub(crate) fn info(&self) {
        log::info!(
            "Detected {} relevant molecule type(s).",
            self.molecule_types().len()
        );

        for molecule in self.molecule_types() {
            molecule.info();
        }
    }
}
