// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for presenting the results of the analysis.

use std::io::Write;

use crate::{
    analysis::{
        molecule::{AtomType, Bond, MoleculeType},
        topology::SystemTopology,
    },
    errors::WriteError,
    Analysis, PANIC_MESSAGE,
};
use serde::{Serialize, Serializer};

macro_rules! write_result {
    ($dst:expr, $($arg:tt)*) => {
        write!($dst, $($arg)*).map_err(|e| WriteError::CouldNotWriteResults(e))?
    };
}

pub(crate) mod aapresenter;
pub(crate) mod ordermap;

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

impl From<&Bond> for BondResults {
    fn from(value: &Bond) -> Self {
        let (total, upper, lower) = value.calc_order();

        BondResults {
            total,
            upper,
            lower,
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
                    None => write_result!(writer, "               "),
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
fn round_serialize_option_f32<S>(x: &Option<f32>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_f32((x.expect(PANIC_MESSAGE) * 10000.0).round() / 10000.0)
}

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
        if self.map().is_some() {
            log::info!("Will calculate order maps.");
        }
        if self.leaflets().is_some() {
            log::info!(
                "Will classify lipids into membrane leaflets using the '{}' method.",
                self.leaflets().as_ref().expect(PANIC_MESSAGE)
            )
        }
        log::info!(
            "Membrane normal expected to be oriented along the {} axis.",
            self.membrane_normal()
        );
    }
}

impl MoleculeType {
    /// Print basic information about the molecule type for the user.
    fn info(&self) {
        log::info!(
            "Molecule type {}: {} order bonds, {} molecules.",
            self.name(),
            self.order_bonds().bonds().len(),
            self.order_bonds()
                .bonds()
                .get(0)
                .expect(PANIC_MESSAGE)
                .bonds()
                .len()
        )
    }
}

impl SystemTopology {
    /// Print basic information about the system topology for the user.
    pub(crate) fn info(&self) {
        log::info!(
            "Detected {} relevant molecule type(s).",
            self.molecules().len()
        );

        for molecule in self.molecules() {
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
