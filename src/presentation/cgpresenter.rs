// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for presenting the results of the analysis of coarse-grained order parameters.

use std::io::Write;

use indexmap::IndexMap;
use serde::{Serialize, Serializer};

use crate::{analysis::molecule::BondTopology, errors::WriteError};

use super::{
    BondResults, CGOrder, MoleculeResults, MoleculeType, ResultsPresenter, SystemTopology,
};

/// Results of the coarse-grained order parameters calculation.
#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub(crate) struct CGOrderResults {
    /// Results for individual molecules of the system.
    molecules: Vec<CGMoleculeResults>,
}

impl From<SystemTopology> for CGOrderResults {
    #[inline(always)]
    fn from(value: SystemTopology) -> Self {
        CGOrderResults {
            molecules: value
                .molecule_types()
                .iter()
                .map(CGMoleculeResults::from)
                .collect(),
        }
    }
}

impl ResultsPresenter for CGOrderResults {
    #[inline(always)]
    fn molecules(&self) -> &[impl MoleculeResults] {
        &self.molecules
    }

    fn write_csv_header(
        writer: &mut impl Write,
        _max_bonds: usize, // needed for compatibility with the trait
        assigned_leaflets: bool,
    ) -> Result<(), WriteError> {
        write_result!(writer, "molecule,atom 1,atom 2");

        if assigned_leaflets {
            write_result!(writer, ",full membrane,upper leaflet,lower leaflet\n");
        } else {
            write_result!(writer, ",full membrane\n");
        }
        Ok(())
    }
}

/// Coarse-grained order parameters calculated for a single molecule type.
#[derive(Debug, Clone, Serialize)]
struct CGMoleculeResults {
    /// Name of the molecule.
    molecule: String,
    /// Order parameters calculated for specific bonds.
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    order: IndexMap<BondTopology, BondResults>,
}

impl Serialize for BondTopology {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let formatted_string = format!(
            "{} {} ({}) - {} {} ({})",
            self.atom1().residue_name(),
            self.atom1().atom_name(),
            self.atom1().relative_index(),
            self.atom2().residue_name(),
            self.atom2().atom_name(),
            self.atom2().relative_index()
        );
        serializer.serialize_str(&formatted_string)
    }
}

impl From<&MoleculeType> for CGMoleculeResults {
    fn from(value: &MoleculeType) -> Self {
        let mut order = IndexMap::new();
        for bond in value.order_bonds().bond_types() {
            let results = BondResults::convert_from::<CGOrder>(bond);
            order.insert(bond.bond_topology().clone(), results);
        }

        CGMoleculeResults {
            molecule: value.name().to_owned(),
            order,
        }
    }
}

impl MoleculeResults for CGMoleculeResults {
    #[inline(always)]
    fn name(&self) -> &str {
        &self.molecule
    }

    #[inline(always)]
    fn max_bonds(&self) -> usize {
        0
    }

    #[inline(always)]
    fn has_assigned_leaflets(&self) -> bool {
        self.order.values().map(|x| x.upper).any(|x| x.is_some())
            || self.order.values().map(|x| x.lower).any(|x| x.is_some())
    }

    fn write_tab(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        write_result!(writer, "\nMolecule type {}\n", self.molecule);

        let leaflets = self.has_assigned_leaflets();

        if leaflets {
            write_result!(writer, "                   FULL     UPPER     LOWER   |\n");
        } else {
            write_result!(writer, "                   FULL   |\n");
        }

        for (bond, results) in self.order.iter() {
            let name = format!(
                "{} - {}",
                bond.atom1().atom_name(),
                bond.atom2().atom_name()
            );
            write_result!(writer, "{:<16}", name);
            results.write_tab(writer, leaflets)?;
            write_result!(writer, "\n");
        }

        Ok(())
    }

    fn write_csv(
        &self,
        writer: &mut impl Write,
        _max_bonds: usize,
        assigned_leaflets: bool,
    ) -> Result<(), WriteError> {
        for (bond, results) in self.order.iter() {
            write_result!(
                writer,
                "{},{},{}",
                self.molecule,
                bond.atom1().atom_name(),
                bond.atom2().atom_name()
            );

            results.write_csv(writer, assigned_leaflets)?;
            write_result!(writer, "\n");
        }

        Ok(())
    }

    fn write_xvg(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        write_result!(
            writer,
            "@    title \"Coarse-grained order parameters for molecule type {}\"\n",
            self.molecule
        );
        write_result!(
            writer,
            "@    xaxis label \"Bond\"\n@    yaxis label \"S\"\n"
        );

        write_result!(writer, "@    s0 legend \"Full membrane\"\n");
        if self.has_assigned_leaflets() {
            write_result!(writer, "@    s1 legend \"Upper leaflet\"\n");
            write_result!(writer, "@    s2 legend \"Lower leaflet\"\n");
        }

        write_result!(writer, "@TYPE xy\n");

        for (i, (bond, results)) in self.order.iter().enumerate() {
            results.write_xvg(writer, i, bond)?;
            write_result!(writer, "\n");
        }

        Ok(())
    }
}
