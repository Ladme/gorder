// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for presenting the results of the analysis.

use crate::{
    analysis::{
        molecule::{AtomType, BondTopology, BondType, MoleculeType},
        topology::SystemTopology,
    },
    errors::WriteError,
    input::Analysis,
    PANIC_MESSAGE,
};
use serde::ser::SerializeMap;
use serde::{Serialize, Serializer};
use std::fmt::Debug;
use std::{
    fs::File,
    io::{BufWriter, Write},
    marker::PhantomData,
    ops::AddAssign,
    path::{Path, PathBuf},
};

macro_rules! write_result {
    ($dst:expr, $($arg:tt)*) => {
        write!($dst, $($arg)*).map_err(|e| WriteError::CouldNotWriteResults(e))?
    };
}

pub(crate) mod aapresenter;
pub(crate) mod cgpresenter;
pub(crate) mod ordermap;

/// Trait for presentation of results for individual molcules.
/// General for both AAOrder and CGOrder.
pub(crate) trait MoleculeResults: Serialize {
    /// Get the name of the molecule.
    fn name(&self) -> &str;

    /// Get the maximal number of order bonds associated with any order atom in this molecule.
    /// For CG, this should be zero.
    fn max_bonds(&self) -> usize;

    /// Check whether order parameters for individual membrane
    /// leaflets have been calculated for at least one order bond of the molecule.
    fn has_assigned_leaflets(&self) -> bool;

    /// Write information about the molecule in human-readable table format.
    fn write_tab(&self, writer: &mut impl Write) -> Result<(), WriteError>;

    /// Write information about the molecule in an xvg format.
    fn write_xvg(&self, writer: &mut impl Write) -> Result<(), WriteError>;

    /// Write information about the molecule in an csv format.
    fn write_csv(
        &self,
        writer: &mut impl Write,
        max_bonds: usize,
        assigned_leaflets: bool,
    ) -> Result<(), WriteError>;
}

/// Trait implemented by structures presenting the results of the order calculation.
/// General for both AAOrder and CGOrder.
pub(crate) trait ResultsPresenter: Serialize {
    /// Get iterator over the results associated with the individual molecules.
    fn molecules(&self) -> impl Iterator<Item=&impl MoleculeResults>;

    /// Write the header for a csv file.
    fn write_csv_header(
        writer: &mut impl Write,
        max_bonds: usize,
        assigned_leaflets: bool,
    ) -> Result<(), WriteError>;

    /// Write the results of the analysis into a yaml file.
    #[inline]
    fn write_yaml(
        &self,
        filename: impl AsRef<Path>,
        input_structure: &str,
        input_trajectory: &str,
        overwrite: bool,
    ) -> Result<(), WriteError> {
        let writer = Self::prepare_file(
            &filename,
            input_structure,
            input_trajectory,
            "yaml",
            overwrite,
            true,
        )?;

        serde_yaml::to_writer(writer, self)
            .map_err(|_| WriteError::CouldNotWriteYaml(Box::from(filename.as_ref())))?;

        Ok(())
    }

    /// Write the results of the analysis into a table.
    #[inline]
    fn write_tab(
        &self,
        filename: impl AsRef<Path>,
        input_structure: &str,
        input_trajectory: &str,
        overwrite: bool,
    ) -> Result<(), WriteError> {
        let mut writer = Self::prepare_file(
            &filename,
            input_structure,
            input_trajectory,
            "tab",
            overwrite,
            true,
        )?;

        for mol in self.molecules() {
            mol.write_tab(&mut writer)?;
        }

        Ok(())
    }

    /// Write the results of the analysis for the individual molecules into xvg files.
    fn write_xvg(
        &self,
        file_pattern: impl AsRef<Path>,
        input_structure: &str,
        input_trajectory: &str,
        overwrite: bool,
    ) -> Result<(), WriteError> {
        let extension = file_pattern
            .as_ref()
            .extension()
            .and_then(|x| x.to_str())
            .map(Some)
            .unwrap_or(None);

        let path_buf = Self::strip_extension(file_pattern.as_ref());
        let file_path = path_buf.to_str().expect(PANIC_MESSAGE);

        // all molecule names must be unique
        let names: Vec<String> = self.molecules().map(|x| x.name().to_owned()).collect();

        for (i, mol) in self.molecules().enumerate() {
            let filename = match extension {
                Some(x) => format!("{}_{}.{}", file_path, names[i], x),
                None => format!("{}_{}", file_path, names[i]),
            };

            log::info!("Writing an xvg file '{}'...", filename);

            let mut writer = Self::prepare_file(
                &filename,
                input_structure,
                input_trajectory,
                "xvg",
                overwrite,
                true,
            )?;

            mol.write_xvg(&mut writer)?;
        }

        Ok(())
    }

    /// Write the results of the analysis into an csv file.
    fn write_csv(&self, filename: impl AsRef<Path>, overwrite: bool) -> Result<(), WriteError> {
        let mut writer = Self::prepare_file(&filename, "", "", "csv", overwrite, false)?;

        let max_bonds = self.molecules().map(|x| x.max_bonds()).max().unwrap_or(0);

        let leaflets = self.molecules().any(|x| x.has_assigned_leaflets());

        Self::write_csv_header(&mut writer, max_bonds, leaflets)?;

        for molecule in self.molecules() {
            molecule.write_csv(&mut writer, max_bonds, leaflets)?;
        }

        Ok(())
    }

    #[inline(always)]
    fn strip_extension(file_path: &Path) -> PathBuf {
        if let Some(stem) = file_path.file_stem() {
            if let Some(parent) = file_path.parent() {
                return parent.join(stem);
            }
        }
        file_path.to_path_buf()
    }

    /// Back up a file, create a new one and write a header into it.
    #[inline(always)]
    fn prepare_file(
        filename: &impl AsRef<Path>,
        input_structure: &str,
        input_trajectory: &str,
        file_type: &str,
        overwrite: bool,
        write_header: bool,
    ) -> Result<BufWriter<File>, WriteError> {
        Self::try_backup_file(filename, overwrite, file_type)?;
        let mut writer = Self::create_and_open_file(filename)?;
        if write_header {
            Self::write_header(&mut writer, filename, input_structure, input_trajectory)?;
        }

        Ok(writer)
    }

    /// Write header into an output file.
    #[inline(always)]
    fn write_header(
        writer: &mut BufWriter<File>,
        filename: &impl AsRef<Path>,
        input_structure: &str,
        input_trajectory: &str,
    ) -> Result<(), WriteError> {
        writeln!(
            writer,
            "# Order parameters calculated with 'gorder v{}' using structure file '{}' and trajectory file '{}'.",
            crate::GORDER_VERSION, input_structure, input_trajectory
        )
            .map_err(|_| WriteError::CouldNotWriteYaml(Box::from(filename.as_ref())))
    }

    /// Create and open file for buffered writing.
    #[inline(always)]
    fn create_and_open_file(filename: &impl AsRef<Path>) -> Result<BufWriter<File>, WriteError> {
        let file = File::create(filename.as_ref())
            .map_err(|_| WriteError::CouldNotCreateFile(Box::from(filename.as_ref())))?;

        Ok(BufWriter::new(file))
    }

    /// Back up an output file, if it is necessary and if it is requested.
    fn try_backup_file(
        filename: &impl AsRef<Path>,
        overwrite: bool,
        file_type: &str,
    ) -> Result<(), WriteError> {
        if filename.as_ref().exists() {
            if !overwrite {
                log::warn!(
                    "Output {} file '{}' already exists. Backing it up.",
                    file_type,
                    filename.as_ref().to_str().expect(PANIC_MESSAGE)
                );
                backitup::backup(filename.as_ref())
                    .map_err(|_| WriteError::CouldNotBackupFile(Box::from(filename.as_ref())))?;
            } else {
                log::warn!(
                    "Output {} file '{}' already exists. It will be overwritten as requested.",
                    file_type,
                    filename.as_ref().to_str().expect(PANIC_MESSAGE)
                );
            }
        }

        Ok(())
    }
}

/// Single order parameter value, optionally with its estimated error.
#[derive(Debug, Clone, Copy)]
pub(super) struct OrderValuePresenter {
    value: f32,
    error: Option<f32>,
}

impl Serialize for OrderValuePresenter {
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

impl From<f32> for OrderValuePresenter {
    fn from(value: f32) -> Self {
        OrderValuePresenter { value, error: None }
    }
}

impl From<[f32; 2]> for OrderValuePresenter {
    fn from(value: [f32; 2]) -> Self {
        OrderValuePresenter {
            value: value[0],
            error: Some(value[1]),
        }
    }
}

impl OrderValuePresenter {
    /// Write the order parameter (and its error) in 'table' format.
    #[inline(always)]
    fn write_tab(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        write_result!(writer, " {: ^8.4} ", self.value);
        Ok(())
    }

    /// Write empty space into 'table' with the correct length.
    #[inline(always)]
    fn write_empty_tab(writer: &mut impl Write, _errors: bool) -> Result<(), WriteError> {
        write_result!(writer, "          ");
        Ok(())
    }

    /// Write the order parameter (and its error) in 'csv' format.
    #[inline(always)]
    fn write_csv(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        write_result!(writer, ",{:.4}", self.value);
        Ok(())
    }

    /// Write empty (non-existent) order value into a 'csv' file
    #[inline(always)]
    fn write_empty_csv(writer: &mut impl Write, _errors: bool) -> Result<(), WriteError> {
        write_result!(writer, ",");
        Ok(())
    }
}

/// Write the order parameter (and its error) in 'table' format or handle missing value.
#[inline(always)]
fn write_optional_order_value_tab(
    value: Option<OrderValuePresenter>,
    writer: &mut impl Write,
    errors: bool,
) -> Result<(), WriteError> {
    match value {
        Some(val) => val.write_tab(writer),
        None => OrderValuePresenter::write_empty_tab(writer, errors),
    }
}

/// Write the order parameter (and its error) in 'csv' format or handle missing value.
#[inline(always)]
fn write_optional_order_value_csv(
    value: Option<OrderValuePresenter>,
    writer: &mut impl Write,
    errors: bool,
) -> Result<(), WriteError> {
    match value {
        Some(val) => val.write_csv(writer),
        None => OrderValuePresenter::write_empty_csv(writer, errors),
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

/// Order parameters calculated for a single bond.
#[derive(Debug, Clone, Serialize)]
struct BondResults {
    total: OrderValuePresenter,
    #[serde(skip_serializing_if = "Option::is_none")]
    upper: Option<OrderValuePresenter>,
    #[serde(skip_serializing_if = "Option::is_none")]
    lower: Option<OrderValuePresenter>,
}

/// Empty struct used as a marker.
#[derive(Default, Debug, Clone)]
pub(crate) struct AAOrder {}
/// Empty struct used as a marker.
#[derive(Default, Debug, Clone)]
pub(crate) struct CGOrder {}

/// Trait implemented only by `AAOrder` and `CGOrder` structs.
pub(crate) trait OrderType: Debug {
    /// Used to convert an order parameter to its final value depending on the analysis type.
    /// Atomistic order parameters are reported as -S_CD.
    /// Coarse grained order parameters are reported as P2.
    fn convert(order: f32, error: Option<f32>) -> OrderValuePresenter;

    /// String to use as a label for z-axis in the ordermap.
    fn zlabel() -> &'static str;
}

impl OrderType for AAOrder {
    #[inline(always)]
    fn convert(order: f32, error: Option<f32>) -> OrderValuePresenter {
        OrderValuePresenter {
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
    fn convert(order: f32, error: Option<f32>) -> OrderValuePresenter {
        OrderValuePresenter {
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
    fn convert_from<O: OrderType>(value: &BondType) -> Self {
        let (total, upper, lower) = value.calc_order();
        let (etotal, eupper, elower) = value.estimate_error();

        BondResults {
            total: O::convert(total, etotal),
            upper: upper.map(|x| O::convert(x, eupper)),
            lower: lower.map(|x| O::convert(x, elower)),
        }
    }
}

impl BondResults {
    /// Write results for a single bond in human readable table format.
    fn write_tab(&self, writer: &mut impl Write, leaflets: bool) -> Result<(), WriteError> {
        self.total.write_tab(writer)?;
        if leaflets {
            for value in [self.upper, self.lower] {
                write_optional_order_value_tab(value, writer, false)?;
            }
        }

        write_result!(writer, "|");

        Ok(())
    }

    /// Write results for a single bond in csv format.
    fn write_csv(&self, writer: &mut impl Write, leaflets: bool) -> Result<(), WriteError> {
        self.total.write_csv(writer)?;

        if leaflets {
            for value in [self.upper, self.lower] {
                write_optional_order_value_csv(value, writer, false)?;
            }
        }

        Ok(())
    }

    /// Write results for a single bond in xvg format.
    fn write_xvg(
        &self,
        writer: &mut impl Write,
        number: usize,
        bond_topology: &BondTopology,
    ) -> Result<(), WriteError> {
        write_result!(
            writer,
            "# Bond {} - {}:\n",
            bond_topology.atom1().atom_name(),
            bond_topology.atom2().atom_name()
        );

        write_result!(writer, "{:<4} {: >8.4} ", number, self.total.value);

        for order in [self.upper, self.lower].into_iter().flatten() {
            write_result!(writer, "{: >8.4} ", order.value);
        }

        Ok(())
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

/// Structure for collecting order parameters to calculate average order.
#[derive(Debug, Clone, Serialize, Default)]
struct AverageOrder<O: OrderType> {
    total: f32,
    total_error: f32,
    total_samples: usize,
    upper: f32,
    upper_error: f32,
    upper_samples: usize,
    lower: f32,
    lower_error: f32,
    lower_samples: usize,
    _type: PhantomData<O>,
}

impl<O: OrderType> AddAssign<&BondType> for AverageOrder<O> {
    #[inline(always)]
    fn add_assign(&mut self, rhs: &BondType) {
        let (total, upper, lower) = rhs.calc_order();
        let (etotal, eupper, elower) = rhs.estimate_error();

        self.total += total;
        let etotal_u = etotal.unwrap_or(f32::NAN);
        self.total_error += etotal_u * etotal_u;
        // we collect the number of bond types, not the number of samples
        // the number of bonds must be the same for each bond type across the system
        self.total_samples += 1;

        if let Some(u) = upper {
            self.upper += u;
            let eupper_u = eupper.unwrap_or(f32::NAN);
            self.upper_error += eupper_u * eupper_u;
            self.upper_samples += 1;
        }

        if let Some(l) = lower {
            self.lower += l;
            let elower_u = elower.unwrap_or(f32::NAN);
            self.lower_error += elower_u * elower_u;
            self.lower_samples += 1;
        }
    }
}

impl<O: OrderType> From<AverageOrder<O>> for BondResults {
    fn from(value: AverageOrder<O>) -> Self {
        fn error2option(val: f32, n_samples: usize) -> Option<f32> {
            if val.is_nan() {
                None
            } else {
                Some(val.sqrt() / n_samples as f32)
            }
        }

        let total_error = error2option(value.total_error, value.total_samples);
        let total = O::convert(value.total / value.total_samples as f32, total_error);

        let upper = if value.upper_samples > 0 {
            let upper_error = error2option(value.upper_error, value.upper_samples);
            Some(O::convert(
                value.upper / value.upper_samples as f32,
                upper_error,
            ))
        } else {
            None
        };

        let lower = if value.lower_samples > 0 {
            let lower_error = error2option(value.lower_error, value.lower_samples);
            Some(O::convert(
                value.lower / value.lower_samples as f32,
                lower_error,
            ))
        } else {
            None
        };

        Self {
            total,
            upper,
            lower,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serialize_atom_type() {
        let atom = AtomType::new_raw(5, "POPC", "CA");
        let serialized = serde_yaml::to_string(&atom).unwrap();
        assert_eq!(serialized, "POPC CA (5)\n");
    }

    #[test]
    fn test_serialize_bond_results() {
        let bond_results = BondResults {
            total: 0.777.into(),
            upper: None,
            lower: None,
        };
        let serialized = serde_yaml::to_string(&bond_results).unwrap();
        assert_eq!(serialized, "total: 0.777\n");
    }

    #[test]
    fn test_serialize_bond_results_leaflets() {
        let bond_results = BondResults {
            total: 0.777.into(),
            upper: Some(0.765.into()),
            lower: Some(0.789.into()),
        };
        let serialized = serde_yaml::to_string(&bond_results).unwrap();
        assert_eq!(serialized, "total: 0.777\nupper: 0.765\nlower: 0.789\n");
    }
}
