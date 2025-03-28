// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

use std::sync::Arc;

use gorder_core::prelude::AAAtomResults as RsAtomResults;
use gorder_core::prelude::BondResults as RsBondResults;
use gorder_core::prelude::Convergence as RsConvergence;
use gorder_core::prelude::GridMapF32;
use gorder_core::prelude::Order as RsOrder;
use gorder_core::prelude::OrderCollection as RsOrderCollection;
use gorder_core::prelude::OrderMapsCollection as RsMapsCollection;
use gorder_core::prelude::PublicMoleculeResults;
use gorder_core::prelude::UAAtomResults as RsUAAtomResults;
use gorder_core::prelude::UAMoleculeResults;
use gorder_core::prelude::{
    AAMoleculeResults, AnalysisResults as RsResults, CGMoleculeResults, PublicOrderResults,
};
use pyo3::prelude::*;

use crate::APIError;
use crate::AtomType;
use crate::WriteError;

/// Class containing all the results of the analysis.
#[pyclass]
pub struct AnalysisResults(pub(crate) Arc<RsResults>);

#[pymethods]
impl AnalysisResults {
    /// Write the results of the analysis into output files.
    pub fn write(&self) -> PyResult<()> {
        if let Err(e) = self.0.write() {
            Err(WriteError::new_err(e.to_string()))
        } else {
            Ok(())
        }
    }

    /// Get the total number of analyzed frames.
    pub fn n_analyzed_frames(&self) -> usize {
        self.0.n_analyzed_frames()
    }

    /// Get the results for each individual molecule.
    pub fn molecules(&self) -> Vec<MoleculeResults> {
        match self.0.as_ref() {
            RsResults::AA(x) => x
                .molecules()
                .map(|x| MoleculeResults {
                    results: self.0.clone(),
                    name: x.molecule().to_owned(),
                })
                .collect::<Vec<MoleculeResults>>(),
            RsResults::CG(x) => x
                .molecules()
                .map(|x| MoleculeResults {
                    results: self.0.clone(),
                    name: x.molecule().to_owned(),
                })
                .collect::<Vec<MoleculeResults>>(),
            RsResults::UA(x) => x
                .molecules()
                .map(|x| MoleculeResults {
                    results: self.0.clone(),
                    name: x.molecule().to_owned(),
                })
                .collect::<Vec<MoleculeResults>>(),
        }
    }

    /// Get the results for a molecule with the specified name.
    /// Raises an exception if such molecule does not exist.
    pub fn get_molecule(&self, name: &str) -> PyResult<MoleculeResults> {
        match self.0.as_ref() {
            RsResults::AA(x) => x
                .get_molecule(name)
                .ok_or_else(|| APIError::new_err("molecule with the given name does not exist"))
                .map(|x| MoleculeResults {
                    results: self.0.clone(),
                    name: x.molecule().to_owned(),
                }),
            RsResults::CG(x) => x
                .get_molecule(name)
                .ok_or_else(|| APIError::new_err("molecule with the given name does not exist"))
                .map(|x| MoleculeResults {
                    results: self.0.clone(),
                    name: x.molecule().to_owned(),
                }),
            RsResults::UA(x) => x
                .get_molecule(name)
                .ok_or_else(|| APIError::new_err("molecule with the given name does not exist"))
                .map(|x| MoleculeResults {
                    results: self.0.clone(),
                    name: x.molecule().to_owned(),
                }),
        }
    }

    /// Get the average order parameters calculated from all bond types of all molecule types.
    pub fn average_order(&self) -> OrderCollection {
        OrderCollection {
            results: self.0.clone(),
            molecule: None,
            identifier: OrderIdentifier::Average,
        }
    }

    /// Get the average order parameter maps calculated from all bond types of all molecule types.
    pub fn average_ordermaps(&self) -> OrderMapsCollection {
        OrderMapsCollection {
            results: self.0.clone(),
            molecule: None,
            identifier: OrderIdentifier::Average,
        }
    }
}

/// Results of the analysis for a single molecule type.
#[pyclass]
pub struct MoleculeResults {
    results: Arc<RsResults>,
    name: String,
}

#[pymethods]
impl MoleculeResults {
    /// Get the name of the molecule for which these results were calculated.
    pub fn molecule(&self) -> String {
        self.name.to_owned()
    }

    /// Get the average order parameters calculated from all bond types of this molecule type.
    pub fn average_order(&self) -> OrderCollection {
        OrderCollection {
            results: self.results.clone(),
            molecule: Some(self.name.clone()),
            identifier: OrderIdentifier::Average,
        }
    }

    /// Get the average order parameter maps calculated from all bond types of this molecule type.
    pub fn average_ordermaps(&self) -> OrderMapsCollection {
        OrderMapsCollection {
            results: self.results.clone(),
            molecule: Some(self.name.clone()),
            identifier: OrderIdentifier::Average,
        }
    }

    /// Get the results for each atom of the molecule.
    /// Not available for coarse-grained order parameters, use `bonds` instead.
    pub fn atoms(&self) -> PyResult<Vec<AtomResults>> {
        match self.results.as_ref() {
            RsResults::AA(x) => {
                Ok(x.get_molecule(&self.name).unwrap().atoms().map(|atom| AtomResults {
                    results: self.results.clone(),
                    molecule: self.name.clone(),
                    atom: atom.atom().relative_index(),
                })
                .collect())
            }
            RsResults::UA(x) => {
                Ok(x.get_molecule(&self.name).unwrap().atoms().map(|atom| AtomResults {
                    results: self.results.clone(),
                    molecule: self.name.clone(),
                    atom: atom.atom().relative_index(),
                })
                .collect())
            }
            RsResults::CG(_) => Err(APIError::new_err(
                "results for individual atoms are not available for coarse-grained order parameters; you want `bonds`",
            )),
        }
    }

    /// Get the results for each bond of the molecule.
    pub fn bonds(&self) -> Vec<BondResults> {
        match self.results.as_ref() {
            RsResults::AA(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .atoms()
                .flat_map(|atom| atom.bonds())
                .map(|bond| BondResults::new(self.results.clone(), bond, self.name.clone()))
                .collect(),
            RsResults::CG(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .bonds()
                .map(|bond| BondResults::new(self.results.clone(), bond, self.name.clone()))
                .collect(),
            RsResults::UA(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .atoms()
                .flat_map(|atom| {
                    atom.bonds().enumerate().map(|(hydrogen, _)| {
                        BondResults::new_ua(
                            self.results.clone(),
                            atom.atom().relative_index(),
                            hydrogen,
                            self.name.clone(),
                        )
                    })
                })
                .collect(),
        }
    }

    /// Get the results for a heavy atom with the specified relative index.
    /// Not available for coarse-grained order parameters, use `get_bond` instead.
    pub fn get_atom(&self, relative_index: usize) -> PyResult<AtomResults> {
        match self.results.as_ref() {
            RsResults::AA(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .get_atom(relative_index)
                .ok_or_else(|| {
                    APIError::new_err("atom with the given relative index does not exist")
                })
                .map(|atom| AtomResults {
                    results: self.results.clone(),
                    molecule: self.name.clone(),
                    atom: atom.atom().relative_index(),
                }),
            RsResults::UA(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .get_atom(relative_index)
                .ok_or_else(|| {
                    APIError::new_err("atom with the given relative index does not exist")
                })
                .map(|atom| AtomResults {
                    results: self.results.clone(),
                    molecule: self.name.clone(),
                    atom: atom.atom().relative_index(),
                }),
            RsResults::CG(_) => Err(APIError::new_err(
                "results for individual atoms are not available for coarse-grained order parameters; you want `get_bond`"
            ))
        }
    }

    /// Get the results for a bond involving atoms with the specified relative indices.
    /// The order of the atom indices does not matter.
    /// Not available for united-atom order parameters.
    pub fn get_bond(
        &self,
        relative_index_1: usize,
        relative_index_2: usize,
    ) -> PyResult<BondResults> {
        match self.results.as_ref() {
            RsResults::AA(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .get_bond(relative_index_1, relative_index_2)
                .ok_or_else(|| {
                    APIError::new_err("bond specified by the given relative indices does not exist")
                })
                .map(|bond| BondResults::new(self.results.clone(), bond, self.name.clone())),
            RsResults::CG(x) => x
                .get_molecule(&self.name)
                .unwrap()
                .get_bond(relative_index_1, relative_index_2)
                .ok_or_else(|| {
                    APIError::new_err("bond specified by the given relative indices does not exist")
                })
                .map(|bond| BondResults::new(self.results.clone(), bond, self.name.clone())),
            RsResults::UA(_) => Err(APIError::new_err(
                "united-atom results for individual bonds cannot be accesed by using relative indices because the hydrogen atoms are virtual and do not have assigned indices",
            )),
        }
    }

    /// Get data about the convergence of the calculated order parameters.
    /// Returns `None` if convergence is not available.
    pub fn convergence(&self) -> Option<Convergence> {
        match self.results.as_ref() {
            RsResults::AA(x) => {
                x.get_molecule(&self.name)
                    .unwrap()
                    .convergence()
                    .map(|_| Convergence {
                        results: self.results.clone(),
                        molecule: self.name.clone(),
                    })
            }
            RsResults::CG(x) => {
                x.get_molecule(&self.name)
                    .unwrap()
                    .convergence()
                    .map(|_| Convergence {
                        results: self.results.clone(),
                        molecule: self.name.clone(),
                    })
            }
            RsResults::UA(x) => {
                x.get_molecule(&self.name)
                    .unwrap()
                    .convergence()
                    .map(|_| Convergence {
                        results: self.results.clone(),
                        molecule: self.name.clone(),
                    })
            }
        }
    }
}

/// Results of the analysis for a single atom type.
#[pyclass]
pub struct AtomResults {
    results: Arc<RsResults>,
    molecule: String,
    atom: usize,
}

#[pymethods]
impl AtomResults {
    /// Get the type of the atom for which these results were calculated.
    pub fn atom(&self) -> AtomType {
        match self.results.as_ref() {
            RsResults::AA(_) => AtomType(self.get_atom_aa_results().atom().clone()),
            RsResults::UA(_) => AtomType(self.get_atom_ua_results().atom().clone()),
            RsResults::CG(_) => unreachable!(
                "FATAL GORDER ERROR | AtomResults::atom | AtomResults should not exist for CG."
            ),
        }
    }

    /// Get the name of the molecule type for which these results were calculated.
    pub fn molecule(&self) -> String {
        self.molecule.clone()
    }

    /// Get the results for each bond of the atom.
    pub fn bonds(&self) -> Vec<BondResults> {
        match self.results.as_ref() {
            RsResults::AA(_) => self
                .get_atom_aa_results()
                .bonds()
                .map(|bond| BondResults::new(self.results.clone(), bond, self.molecule.clone()))
                .collect(),
            RsResults::UA(_) => self
                .get_atom_ua_results()
                .bonds()
                .enumerate()
                .map(|(hydrogen, _)| {
                    BondResults::new_ua(
                        self.results.clone(),
                        self.atom,
                        hydrogen,
                        self.molecule.clone(),
                    )
                })
                .collect(),
            RsResults::CG(_) => unreachable!(
                "FATAL GORDER ERROR | AtomResults::bonds | AtomResults should not exist for CG."
            ),
        }
    }

    /// Get the results for a bond between this heavy atom and the hydrogen atom with the specified relative index (AA) or
    /// get the results for a bond between this heavy atom and the virtual hydrogen atom with the specified index (UA).
    pub fn get_bond(&self, relative_index: usize) -> PyResult<BondResults> {
        match self.results.as_ref() {
            RsResults::AA(_) => self
                .get_atom_aa_results()
                .get_bond(relative_index)
                .ok_or_else(|| {
                    APIError::new_err(
                        "heavy atom is not bonded to hydrogen with the given relative index",
                    )
                })
                .map(|bond| BondResults::new(self.results.clone(), bond, self.molecule.clone())),
            RsResults::UA(_) => self
                .get_atom_ua_results()
                .bonds()
                .nth(relative_index)
                .ok_or_else(|| {
                    APIError::new_err("carbon does not have a virtual hydrogen with this index")
                })
                .map(|_| {
                    BondResults::new_ua(
                        self.results.clone(),
                        self.atom,
                        relative_index,
                        self.molecule.clone(),
                    )
                }),
            RsResults::CG(_) => unreachable!(
                "FATAL GORDER ERROR | AtomResults::get_bond | AtomResults should not exist for CG."
            ),
        }
    }

    /// Get the order parameters calculated for this atom.
    pub fn order(&self) -> OrderCollection {
        OrderCollection {
            results: self.results.clone(),
            molecule: Some(self.molecule.clone()),
            identifier: OrderIdentifier::Atom(self.atom),
        }
    }

    /// Get the ordermaps calculated for this atom.
    pub fn ordermaps(&self) -> OrderMapsCollection {
        OrderMapsCollection {
            results: self.results.clone(),
            molecule: Some(self.molecule.clone()),
            identifier: OrderIdentifier::Atom(self.atom),
        }
    }
}

impl AtomResults {
    /// Helper method for obtaining reference to the results for this AA atom.
    fn get_atom_aa_results(&self) -> &RsAtomResults {
        match self.results.as_ref() {
            RsResults::AA(x) => x
                .get_molecule(&self.molecule)
                .unwrap()
                .get_atom(self.atom)
                .unwrap(),
            RsResults::CG(_) | RsResults::UA(_) => unreachable!(
                "FATAL GORDER ERROR | AtomResults::get_atom_results | Results should be atomistic."
            ),
        }
    }

    /// Helper method for obtaining reference to the results for this UA atom.
    fn get_atom_ua_results(&self) -> &RsUAAtomResults {
        match self.results.as_ref() {
            RsResults::AA(_) | RsResults::CG(_) => unreachable!("FATAL GORDER ERROR | AtomResults::get_atom_ua_results | Results should be united-atom."),
            RsResults::UA(x) => x.get_molecule(&self.molecule).unwrap().get_atom(self.atom).unwrap(),
        }
    }
}

/// Results of the analysis for a single bond type.
#[pyclass]
pub struct BondResults {
    results: Arc<RsResults>,
    molecule: String,
    // Relative indices of the involved atoms (CG, AA) OR
    // relative index of the involved atom and the bond index (UA).
    bond: (usize, usize),
}

#[pymethods]
impl BondResults {
    /// Get the name of the molecule type for which these results were calculated.
    pub fn molecule(&self) -> String {
        self.molecule.to_owned()
    }

    /// Get the atom types involved in this bond type.
    pub fn atoms(&self) -> Result<(AtomType, AtomType), PyErr> {
        match self.results.as_ref() {
            RsResults::AA(_) | RsResults::CG(_) => {
                let atoms = self.get_bond_results().atoms();
                Ok((AtomType(atoms.0.clone()), AtomType(atoms.1.clone())))
            }
            RsResults::UA(_) => Err(APIError::new_err(
                "cannot access information about atoms in a virtual united-atom bond; the bond only involves one real atom"))
        }
    }

    /// Get the order parameters calculated for this bond.
    pub fn order(&self) -> OrderCollection {
        let identifier = match self.results.as_ref() {
            RsResults::AA(_) | RsResults::CG(_) => OrderIdentifier::Bond(self.bond.0, self.bond.1),
            RsResults::UA(_) => OrderIdentifier::VirtualBond(self.bond.0, self.bond.1),
        };

        OrderCollection {
            results: self.results.clone(),
            molecule: Some(self.molecule.clone()),
            identifier: identifier,
        }
    }

    /// Get the maps of order parameters calculated for this bond.
    pub fn ordermaps(&self) -> OrderMapsCollection {
        let identifier = match self.results.as_ref() {
            RsResults::AA(_) | RsResults::CG(_) => OrderIdentifier::Bond(self.bond.0, self.bond.1),
            RsResults::UA(_) => OrderIdentifier::VirtualBond(self.bond.0, self.bond.1),
        };

        OrderMapsCollection {
            results: self.results.clone(),
            molecule: Some(self.molecule.clone()),
            identifier: identifier,
        }
    }
}

impl BondResults {
    /// Create a new BondResults wrapper.
    fn new(all_results: Arc<RsResults>, bond_results: &RsBondResults, molecule: String) -> Self {
        BondResults {
            results: all_results,
            molecule,
            bond: (
                bond_results.atoms().0.relative_index(),
                bond_results.atoms().1.relative_index(),
            ),
        }
    }

    /// Create a new BondResults wrapper for the united-atom bond.
    fn new_ua(
        all_results: Arc<RsResults>,
        carbon_index: usize,
        hydrogen_index: usize,
        molecule: String,
    ) -> Self {
        BondResults {
            results: all_results,
            molecule,
            bond: (carbon_index, hydrogen_index),
        }
    }

    /// Helper method for obtaining reference to the results for this bond.
    /// Panics, if used for UA BondResults.
    fn get_bond_results(&self) -> &RsBondResults {
        match self.results.as_ref() {
            RsResults::AA(x) => x
                .get_molecule(&self.molecule)
                .unwrap()
                .get_bond(self.bond.0, self.bond.1)
                .unwrap(),
            RsResults::CG(x) => x
                .get_molecule(&self.molecule)
                .unwrap()
                .get_bond(self.bond.0, self.bond.1)
                .unwrap(),
            RsResults::UA(_) => unreachable!("FATAL GORDER ERROR | BondResults::get_bond_results | This method cannot be used for united-atom results."),
        }
    }
}

/// Helper enum for identifying leaflet for which order parameters / ordermaps should be collected.
#[derive(Debug, Clone, Copy)]
enum Leaflet {
    Upper,
    Lower,
    Total,
}

impl Leaflet {
    #[inline]
    fn get_order(&self, collection: &RsOrderCollection) -> Option<Order> {
        match self {
            Self::Upper => collection.upper().map(Order),
            Self::Lower => collection.lower().map(Order),
            Self::Total => collection.total().map(Order),
        }
    }

    #[inline]
    fn get_ordermap<'a>(&self, collection: &'a RsMapsCollection) -> Option<&'a GridMapF32> {
        match self {
            Self::Upper => collection.upper().as_ref(),
            Self::Lower => collection.lower().as_ref(),
            Self::Total => collection.total().as_ref(),
        }
    }
}

/// Helper enum for identifying bond or
/// atom for which order parameters / ordermaps should be collected.
#[derive(Debug, Clone)]
enum OrderIdentifier {
    Average,
    Bond(usize, usize),
    Atom(usize),
    VirtualBond(usize, usize),
}

impl OrderIdentifier {
    #[inline]
    fn get_order_aa(&self, mol_results: &AAMoleculeResults, leaflet: Leaflet) -> Option<Order> {
        match self {
            Self::Average => leaflet.get_order(mol_results.average_order()),
            Self::Bond(x, y) => leaflet.get_order(mol_results.get_bond(*x, *y)?.order()),
            Self::Atom(x) => leaflet.get_order(mol_results.get_atom(*x)?.order()),
            Self::VirtualBond(_, _) => unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_order_aa | Virtual bond identifier cannot be used for AA."),
        }
    }

    #[inline]
    fn get_order_cg(&self, mol_results: &CGMoleculeResults, leaflet: Leaflet) -> Option<Order> {
        match self {
            Self::Average => leaflet.get_order(mol_results.average_order()),
            Self::Bond(x, y) => leaflet.get_order(mol_results.get_bond(*x, *y)?.order()),
            Self::Atom(_) => unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_order_cg | Atom identifier cannot be used for CG."),
            Self::VirtualBond(_, _) => unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_order_cg | Virtual bond identifier cannot be used for CG."),
        }
    }

    #[inline]
    fn get_order_ua(&self, mol_results: &UAMoleculeResults, leaflet: Leaflet) -> Option<Order> {
        match self {
            Self::Average => leaflet.get_order(mol_results.average_order()),
            Self::Bond(_, _) => unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_order_ua | Bond identifier cannot be used for UA."),
            Self::Atom(x) => leaflet.get_order(mol_results.get_atom(*x)?.order()),
            Self::VirtualBond(x, y) => leaflet.get_order(mol_results.get_atom(*x)?.bonds().nth(*y)?.order()),
        }
    }

    #[inline]
    fn get_ordermap_aa<'a>(
        &self,
        mol_results: &'a AAMoleculeResults,
        leaflet: Leaflet,
    ) -> Option<&'a GridMapF32> {
        match self {
            Self::Average => leaflet.get_ordermap(mol_results.average_ordermaps()),
            Self::Bond(x, y) => leaflet.get_ordermap(mol_results.get_bond(*x, *y)?.ordermaps()),
            Self::Atom(x) => leaflet.get_ordermap(mol_results.get_atom(*x)?.ordermaps()),
            Self::VirtualBond(_, _) => {
                unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_ordermap_aa | Virtual bond identifier cannot be used for AA.")
            }
        }
    }

    #[inline]
    fn get_ordermap_cg<'a>(
        &self,
        mol_results: &'a CGMoleculeResults,
        leaflet: Leaflet,
    ) -> Option<&'a GridMapF32> {
        match self {
            Self::Average => leaflet.get_ordermap(mol_results.average_ordermaps()),
            Self::Bond(x, y) => {
                leaflet.get_ordermap(mol_results.get_bond(*x, *y)?.ordermaps())
            }
            Self::Atom(_) => panic!("FATAL GORDER ERROR | OrderIdentifier::get_ordermap_cg | Atom identifier cannot be used for  CG."),
            Self::VirtualBond(_, _) => {
                unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_ordermap_cg | Virtual bond identifier cannot be used for CG.")
            }
        }
    }

    #[inline]
    fn get_ordermap_ua<'a>(
        &self,
        mol_results: &'a UAMoleculeResults,
        leaflet: Leaflet,
    ) -> Option<&'a GridMapF32> {
        match self {
            Self::Average => leaflet.get_ordermap(mol_results.average_ordermaps()),
            Self::Bond(_, _) => unreachable!("FATAL GORDER ERROR | OrderIdentifier::get_ordermap_ua | Bond identifier cannot be used for UA."),
            Self::Atom(x) => leaflet.get_ordermap(mol_results.get_atom(*x)?.ordermaps()),
            Self::VirtualBond(x, y) => leaflet.get_ordermap(mol_results.get_atom(*x)?.bonds().nth(*y)?.ordermaps()),
        }
    }
}

/// Order parameters for a single object calculated for the full membrane, the upper leaflet, and the lower leaflet.
#[pyclass]
pub struct OrderCollection {
    results: Arc<RsResults>,
    molecule: Option<String>, // `None` for average order for the entire system
    identifier: OrderIdentifier,
}

#[pymethods]
impl OrderCollection {
    /// Get order parameter calculated from the whole membrane.
    /// Returns `None` if the order parameter is not available.
    pub fn total(&self) -> Option<Order> {
        self.get_order(Leaflet::Total)
    }

    /// Get order parameter calculated from the upper leaflet.
    /// Returns `None` if the order parameter is not available.
    pub fn upper(&self) -> Option<Order> {
        self.get_order(Leaflet::Upper)
    }

    /// Get order parameter calculated from the lower leaflet.
    /// Returns `None` if the order parameter is not available.
    pub fn lower(&self) -> Option<Order> {
        self.get_order(Leaflet::Lower)
    }
}

impl OrderCollection {
    /// Helper method for getting order parameters for a leaflet.
    fn get_order(&self, leaflet: Leaflet) -> Option<Order> {
        match &self.molecule {
            None => match self.results.as_ref() {
                RsResults::AA(results) => leaflet.get_order(results.average_order()),
                RsResults::CG(results) => leaflet.get_order(results.average_order()),
                RsResults::UA(results) => leaflet.get_order(results.average_order()),
            },
            Some(mol) => match self.results.as_ref() {
                RsResults::AA(results) => self
                    .identifier
                    .get_order_aa(results.get_molecule(mol)?, leaflet),
                RsResults::CG(results) => self
                    .identifier
                    .get_order_cg(results.get_molecule(mol)?, leaflet),
                RsResults::UA(results) => self
                    .identifier
                    .get_order_ua(results.get_molecule(mol)?, leaflet),
            },
        }
    }
}

/// Single order parameter value, optionally with its estimated error.
#[pyclass]
pub struct Order(pub(crate) RsOrder);

#[pymethods]
impl Order {
    /// Get the value of the order parameter (mean from the analyzed frames).
    pub fn value(&self) -> f32 {
        self.0.value()
    }

    /// Get the estimated error for this order parameter.
    /// Returns `None` if the error is not available.
    pub fn error(&self) -> Option<f32> {
        self.0.error()
    }
}

/// Ordermaps for a single object calculated for the full membrane, the upper leaflet, and the lower leaflet.
#[pyclass]
#[derive(Clone)]
pub struct OrderMapsCollection {
    results: Arc<RsResults>,
    molecule: Option<String>, // `None` for average order for the entire system
    identifier: OrderIdentifier,
}

#[pymethods]
impl OrderMapsCollection {
    /// Get order parameter map calculated from the whole membrane.
    /// Returns `None` if the ordermap is not available.
    pub fn total(&self) -> Option<Map> {
        self.get_ordermap(Leaflet::Total).map(|_| Map {
            collection: self.clone(),
            leaflet: Leaflet::Total,
        })
    }

    /// Get order parameter map calculated from the upper leaflet.
    /// Returns `None` if the ordermap is not available.
    pub fn upper(&self) -> Option<Map> {
        self.get_ordermap(Leaflet::Upper).map(|_| Map {
            collection: self.clone(),
            leaflet: Leaflet::Upper,
        })
    }

    /// Get order parameter map calculated from the lower leaflet.
    /// Returns `None` if the ordermap is not available.
    pub fn lower(&self) -> Option<Map> {
        self.get_ordermap(Leaflet::Lower).map(|_| Map {
            collection: self.clone(),
            leaflet: Leaflet::Lower,
        })
    }
}

impl OrderMapsCollection {
    /// Helper method for getting ordermaps for a leaflet.
    fn get_ordermap(&self, leaflet: Leaflet) -> Option<&GridMapF32> {
        match &self.molecule {
            None => match self.results.as_ref() {
                RsResults::AA(results) => leaflet.get_ordermap(results.average_ordermaps()),
                RsResults::CG(results) => leaflet.get_ordermap(results.average_ordermaps()),
                RsResults::UA(results) => leaflet.get_ordermap(results.average_ordermaps()),
            },
            Some(mol) => match self.results.as_ref() {
                RsResults::AA(results) => self
                    .identifier
                    .get_ordermap_aa(results.get_molecule(mol)?, leaflet),
                RsResults::CG(results) => self
                    .identifier
                    .get_ordermap_cg(results.get_molecule(mol)?, leaflet),
                RsResults::UA(results) => self
                    .identifier
                    .get_ordermap_ua(results.get_molecule(mol)?, leaflet),
            },
        }
    }
}

/// Map of order parameters.
#[pyclass]
pub struct Map {
    collection: OrderMapsCollection,
    leaflet: Leaflet,
}

#[pymethods]
impl Map {
    /// Get the span of the map in the x-dimension.
    pub fn span_x(&self) -> (f32, f32) {
        let ordermap = self.collection.get_ordermap(self.leaflet).unwrap();
        ordermap.span_x()
    }

    /// Get the span of the map along the y-dimension.
    pub fn span_y(&self) -> (f32, f32) {
        let ordermap = self.collection.get_ordermap(self.leaflet).unwrap();
        ordermap.span_y()
    }

    /// Get the dimnesions of a single grid tile of the map.
    pub fn tile_dim(&self) -> (f32, f32) {
        let ordermap = self.collection.get_ordermap(self.leaflet).unwrap();
        ordermap.tile_dim()
    }

    /// Get the value of the order parameter at the specified coordinates.
    pub fn get_at(&self, x: f32, y: f32) -> Option<f32> {
        let ordermap = self.collection.get_ordermap(self.leaflet).unwrap();
        ordermap.get_at_convert(x, y)
    }

    /// Extract the order map into NumPy arrays.
    ///
    /// Returns
    /// -------
    /// Tuple[np.ndarray, np.ndarray, np.ndarray]
    ///     A tuple of NumPy arrays.
    ///     The first array (1D) contains positions of the grid tiles along the x-axis.
    ///     The second array (1D) contains positions of the grid tiles along the y-axis.
    ///     The third array (2D) contains the calculated order parameters.
    #[allow(clippy::type_complexity)]
    pub fn extract(
        &self,
        py: Python<'_>,
    ) -> (
        Py<numpy::PyArray1<f32>>,
        Py<numpy::PyArray1<f32>>,
        Py<numpy::PyArray2<f32>>,
    ) {
        let ordermap = self.collection.get_ordermap(self.leaflet).unwrap();
        (
            Self::x_positions(ordermap, py),
            Self::y_positions(ordermap, py),
            Self::values_array(ordermap, py),
        )
    }
}

impl Map {
    fn x_positions(ordermap: &GridMapF32, py: Python<'_>) -> Py<numpy::PyArray1<f32>> {
        let start = ordermap.span_x().0;
        let step = ordermap.tile_dim().0;
        let n = ordermap.n_tiles_x();

        let x_coords: Vec<f32> = (0..n).map(|i| start + i as f32 * step).collect();

        numpy::PyArray1::from_vec(py, x_coords).into()
    }

    fn y_positions(ordermap: &GridMapF32, py: Python<'_>) -> Py<numpy::PyArray1<f32>> {
        let start = ordermap.span_y().0;
        let step = ordermap.tile_dim().1;
        let n = ordermap.n_tiles_y();

        let y_coords: Vec<f32> = (0..n).map(|i| start + i as f32 * step).collect();

        numpy::PyArray1::from_vec(py, y_coords).into()
    }

    fn values_array(ordermap: &GridMapF32, py: Python<'_>) -> Py<numpy::PyArray2<f32>> {
        let mut converted_values: Vec<f32> =
            ordermap.extract_convert().map(|(_, _, val)| val).collect();

        let n_x = ordermap.n_tiles_x();
        let n_y = ordermap.n_tiles_y();

        let mut values_2d = Vec::with_capacity(n_x);
        for chunk in converted_values.chunks_mut(n_y) {
            values_2d.push(chunk.to_vec());
        }

        numpy::PyArray2::from_vec2(py, &values_2d).expect(
            "FATAL GORDER ERROR | OrderMap::values_array | Could not convert ndarray to numpy array.").into()
    }
}

/// Stores information about the convergence of calculations,
/// i.e. the average order parameter for this molecule collected in time.
#[pyclass]
pub struct Convergence {
    results: Arc<RsResults>,
    molecule: String,
}

#[pymethods]
impl Convergence {
    /// Extract indices of trajectory frames where order parameters were calculated.
    /// The first analyzed frame is assigned an index of 1. For instance, if the
    /// analysis begins at 200 ns, the frame at or just after 200 ns is indexed as 1.
    pub fn frames(&self) -> Vec<usize> {
        self.get_convergence().frames().clone()
    }

    /// Extract cumulative average order parameters for this molecule, calculated across the entire membrane.
    ///
    /// Each value in the vector represents the cumulative average for a molecule up to a specific frame:
    /// The first element corresponds to the value from the first frame.
    /// The second element is the average from the first two frames.
    /// The third element is the average from the first three frames, and so on.
    /// The last element is the overall average across all frames.
    pub fn total(&self) -> Option<Vec<f32>> {
        self.get_convergence().total().clone()
    }

    /// Extract umulative average order parameters for this molecule, calculated only for the upper leaflet.
    /// Follows the same format and logic as ``total``.
    pub fn upper(&self) -> Option<Vec<f32>> {
        self.get_convergence().upper().clone()
    }

    /// Extract cumulative average order parameters for this molecule, calculated only for the lower leaflet.
    /// Follows the same format and logic as ``total``.
    pub fn lower(&self) -> Option<Vec<f32>> {
        self.get_convergence().lower().clone()
    }
}

impl Convergence {
    /// Helper method for obtaining reference to the Convergence structure.
    fn get_convergence(&self) -> &RsConvergence {
        match self.results.as_ref() {
            RsResults::AA(results) => results
                .get_molecule(&self.molecule)
                .unwrap()
                .convergence()
                .unwrap(),
            RsResults::CG(results) => results
                .get_molecule(&self.molecule)
                .unwrap()
                .convergence()
                .unwrap(),
            RsResults::UA(results) => results
                .get_molecule(&self.molecule)
                .unwrap()
                .convergence()
                .unwrap(),
        }
    }
}
