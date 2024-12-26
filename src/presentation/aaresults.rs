// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for presenting the results of the analysis of atomistic order parameters.

use super::{AAOrder, BondResults, MoleculeResults, OrderCollection, OrderResults};
use crate::analysis::molecule::AtomType;
use crate::input::Analysis;
use crate::presentation::ordermap::OrderMapsCollection;
use getset::Getters;
use indexmap::IndexMap;
use serde::Serialize;

/// Results of the atomistic order parameters calculation.
#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub(crate) struct AAOrderResults {
    /// Results for individual molecules of the system.
    molecules: IndexMap<String, AAMoleculeResults>,
    /// Parameters of the analysis.
    #[serde(skip)]
    analysis: Analysis,
}

impl OrderResults for AAOrderResults {
    type OrderType = AAOrder;
    type MoleculeResults = AAMoleculeResults;

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

    /// Get the maximal number of bonds per heavy atoms in the system.
    /// Returns 0 if there are no molecules.
    fn max_bonds(&self) -> usize {
        self.molecules
            .values()
            .map(|mol| mol.max_bonds())
            .max()
            .unwrap_or(0)
    }
}

/// Atomistic order parameters calculated for a single molecule type.
#[derive(Debug, Clone, Serialize, Getters)]
pub(crate) struct AAMoleculeResults {
    /// Name of the molecule type.
    #[serde(skip)]
    molecule: String,
    /// Average order parameter calculated from all bond types of this molecule type.
    #[serde(rename = "average order")]
    #[getset(get = "pub(super)")]
    average_order: OrderCollection,
    /// Average order parameter maps calculated from all bond types of this molecule type.
    #[serde(skip)]
    #[getset(get = "pub(super)")]
    average_ordermaps: OrderMapsCollection,
    /// Order parameters calculated for individual atom types.
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    #[serde(rename = "order parameters")]
    #[getset(get = "pub(super)")]
    order: IndexMap<AtomType, AAAtomResults>,
}

impl AAMoleculeResults {
    pub(super) fn new(
        molecule: &str,
        average_order: OrderCollection,
        average_ordermaps: OrderMapsCollection,
        order: IndexMap<AtomType, AAAtomResults>,
    ) -> Self {
        Self {
            molecule: molecule.to_owned(),
            average_order,
            average_ordermaps,
            order,
        }
    }
}

impl MoleculeResults for AAMoleculeResults {
    fn molecule(&self) -> &str {
        &self.molecule
    }

    #[inline(always)]
    fn max_bonds(&self) -> usize {
        self.order
            .values()
            .map(|x| x.bonds.len())
            .max()
            .unwrap_or(0)
    }
}

/// Atomistic order parameters calculated for a single atom type and for involved bond types.
#[derive(Debug, Clone, Serialize, Getters)]
pub(super) struct AAAtomResults {
    /// Name of the atom type.
    #[serde(skip)]
    #[getset(get = "pub(super)")]
    atom: AtomType,
    /// Name of the molecule this atom is part of.
    #[serde(skip)]
    #[getset(get = "pub(super)")]
    molecule: String,
    /// Order parameters calculated for this atom.
    #[serde(flatten)]
    #[getset(get = "pub(super)")]
    order: OrderCollection,
    /// Ordermaps calculated for this atom.
    #[serde(skip)]
    #[getset(get = "pub(super)")]
    ordermaps: OrderMapsCollection,
    /// Order parameters calculated for bond types that this atom type is involved in.
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    #[getset(get = "pub(super)")]
    bonds: IndexMap<AtomType, BondResults>,
}

impl AAAtomResults {
    pub(super) fn new(
        atom: AtomType,
        molecule: &str,
        order: OrderCollection,
        ordermaps: OrderMapsCollection,
        bonds: IndexMap<AtomType, BondResults>,
    ) -> Self {
        Self {
            atom,
            molecule: molecule.to_owned(),
            order,
            ordermaps,
            bonds,
        }
    }

    /// Check whether the atom has any bonds associated with it.
    #[inline(always)]
    pub(super) fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }
}
