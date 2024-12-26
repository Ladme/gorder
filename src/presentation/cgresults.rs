// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for presenting the results of the analysis of coarse-grained order parameters.

use getset::Getters;
use indexmap::IndexMap;
use serde::{Serialize, Serializer};

use super::{BondResults, CGOrder, MoleculeResults, OrderCollection, OrderResults};
use crate::analysis::molecule::BondTopology;
use crate::input::Analysis;
use crate::presentation::ordermap::OrderMapsCollection;

/// Results of the coarse-grained order parameters calculation.
#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub struct CGOrderResults {
    /// Results for individual molecules of the system.
    molecules: IndexMap<String, CGMoleculeResults>,
    /// Parameters of the analysis.
    #[serde(skip)]
    analysis: Analysis,
}

impl OrderResults for CGOrderResults {
    type OrderType = CGOrder;
    type MoleculeResults = CGMoleculeResults;

    fn empty(analysis: Analysis) -> Self {
        Self {
            molecules: IndexMap::new(),
            analysis,
        }
    }

    fn new(molecules: IndexMap<String, Self::MoleculeResults>, analysis: Analysis) -> Self {
        Self {
            molecules,
            analysis,
        }
    }

    fn molecules(&self) -> impl Iterator<Item = &Self::MoleculeResults> {
        self.molecules.values()
    }

    fn analysis(&self) -> &Analysis {
        &self.analysis
    }
}

/// Coarse-grained order parameters calculated for a single molecule type.
#[derive(Debug, Clone, Serialize, Getters)]
pub struct CGMoleculeResults {
    /// Name of the molecule.
    #[serde(skip)]
    #[getset(get = "pub(super)")]
    molecule: String,
    /// Average order parameter for all bond types of the molecule.
    #[serde(rename = "average order")]
    #[getset(get = "pub(super)")]
    average_order: OrderCollection,
    /// Average order parameter maps calculated from all bond types of this molecule type.
    #[serde(skip)]
    #[getset(get = "pub(super)")]
    average_ordermaps: OrderMapsCollection,
    /// Order parameters calculated for specific bonds.
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    #[serde(rename = "order parameters")]
    #[getset(get = "pub(super)")]
    order: IndexMap<BondTopology, BondResults>,
}

impl CGMoleculeResults {
    pub(super) fn new(
        name: &str,
        average_order: OrderCollection,
        average_ordermaps: OrderMapsCollection,
        order: IndexMap<BondTopology, BondResults>,
    ) -> Self {
        CGMoleculeResults {
            molecule: name.to_owned(),
            average_order,
            average_ordermaps,
            order,
        }
    }
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

impl MoleculeResults for CGMoleculeResults {
    fn molecule(&self) -> &str {
        &self.molecule
    }
}
