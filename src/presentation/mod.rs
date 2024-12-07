// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for presenting the results of the analysis.

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use crate::{
    analysis::{
        molecule::{AtomType, BondTopology, BondType, MoleculeType},
        topology::SystemTopology,
    },
    errors::WriteError,
    input::Analysis,
    PANIC_MESSAGE,
};
use serde::{Serialize, Serializer};

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
    /// Get reference to the results associated with the individual molecules.
    fn molecules(&self) -> &[impl MoleculeResults];

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

        for mol in self.molecules().iter() {
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
        let names: Vec<String> = self
            .molecules()
            .iter()
            .map(|x| x.name().to_owned())
            .collect();

        for (i, mol) in self.molecules().iter().enumerate() {
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

        let max_bonds = self
            .molecules()
            .iter()
            .map(|x| x.max_bonds())
            .max()
            .unwrap_or(0);

        let leaflets = self.molecules().iter().any(|x| x.has_assigned_leaflets());

        Self::write_csv_header(&mut writer, max_bonds, leaflets)?;

        for molecule in self.molecules().iter() {
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

/// Order parameters calculated for a single bond.
#[derive(Debug, Clone, Serialize)]
struct BondResults {
    #[serde(serialize_with = "round_serialize_f32")]
    total: f32,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_serialize_option_f32")]
    upper: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_serialize_option_f32")]
    lower: Option<f32>,
}

/// Empty struct used as a marker.
pub(crate) struct AAOrder {}
/// Empty struct used as a marker.
pub(crate) struct CGOrder {}

/// Trait implemented only by `AAOrder` and `CGOrder` structs.
pub(crate) trait OrderType {
    /// Used to convert an order parameter to its final value depending on the analysis type.
    /// Atomistic order parameters are reported as -S_CD.
    /// Coarse grained order parameters are reported as P2.
    fn convert(order: f32) -> f32;
}

impl OrderType for AAOrder {
    #[inline(always)]
    fn convert(order: f32) -> f32 {
        -order
    }
}

impl OrderType for CGOrder {
    #[inline(always)]
    fn convert(order: f32) -> f32 {
        order
    }
}

impl BondResults {
    #[inline(always)]
    fn convert_from<O: OrderType>(value: &BondType) -> Self {
        let (total, upper, lower) = value.calc_order();

        BondResults {
            total: O::convert(total),
            upper: upper.map(|x| O::convert(x)),
            lower: lower.map(|x| O::convert(x)),
        }
    }
}

impl BondResults {
    /// Write results for a single bond in human readable table format.
    fn write_tab(&self, writer: &mut impl Write, leaflets: bool) -> Result<(), WriteError> {
        write_result!(writer, " {: ^8.4} ", self.total);
        if leaflets {
            for value in [self.upper, self.lower] {
                match value {
                    Some(unwrapped) => write_result!(writer, " {: ^8.4} ", unwrapped),
                    None => write_result!(writer, "          "),
                }
            }
        }

        write_result!(writer, "|");

        Ok(())
    }

    /// Write results for a single bond in csv format.
    fn write_csv(&self, writer: &mut impl Write, leaflets: bool) -> Result<(), WriteError> {
        write_result!(writer, ",{:.4}", self.total);

        if leaflets {
            for value in [self.upper, self.lower] {
                match value {
                    Some(unwrapped) => write_result!(writer, ",{:.4}", unwrapped),
                    None => write_result!(writer, ","),
                }
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

        write_result!(writer, "{:<4} {: >8.4} ", number, self.total);

        for order in [self.upper, self.lower].into_iter() {
            if let Some(value) = order {
                write_result!(writer, "{: >8.4} ", value);
            }
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

/// Assumes that `None` value is ignored.
#[inline(always)]
fn round_serialize_option_f32<S>(x: &Option<f32>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_f32((x.expect(PANIC_MESSAGE) * 10000.0).round() / 10000.0)
}

#[inline(always)]
fn round_serialize_f32<S>(x: &f32, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_f32((x * 10000.0).round() / 10000.0)
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
                .get(0)
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
            total: 0.777,
            upper: None,
            lower: None,
        };
        let serialized = serde_yaml::to_string(&bond_results).unwrap();
        assert_eq!(serialized, "total: 0.777\n");
    }

    #[test]
    fn test_serialize_bond_results_leaflets() {
        let bond_results = BondResults {
            total: 0.777,
            upper: Some(0.765),
            lower: Some(0.789),
        };
        let serialized = serde_yaml::to_string(&bond_results).unwrap();
        assert_eq!(serialized, "total: 0.777\nupper: 0.765\nlower: 0.789\n");
    }
}
