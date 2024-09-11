// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use groan_rs::system::System;

use crate::{auxiliary::create_group, errors::TopologyError};

/// Parameters for the classification of lipids into membrane leaflets.
#[derive(Debug, Clone)]
#[allow(private_interfaces)]
pub enum LeafletClassification {
    Global(GlobalParams),
    Local(LocalParams),
    Individual(IndividualParams),
}

impl LeafletClassification {
    /// Classify lipids based on the global membrane center of geometry.
    /// Useful for disrupted membranes.
    ///
    /// ## Parameters
    /// - `membrane`: selection of all lipids forming the membrane
    /// - `heads`: reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///            there must only be one such atom/bead per lipid molecule.
    pub fn global(membrane: &str, heads: &str) -> LeafletClassification {
        Self::Global(GlobalParams {
            membrane: membrane.to_string(),
            heads: heads.to_string(),
        })
    }

    /// Classify lipids based on the local membrane center of geometry.
    /// Useful for curved membranes, slow.
    ///
    /// ## Parameters
    /// - `membrane`: selection of all lipids forming the membrane
    /// - `heads`: reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///            there must only be one such atom/bead per lipid molecule.
    /// - `radius`: radius of a cylinder for the calculation of local membrane center of geometry (in nm)
    pub fn local(membrane: &str, heads: &str, radius: f32) -> LeafletClassification {
        Self::Local(LocalParams {
            membrane: membrane.to_string(),
            heads: heads.to_string(),
            radius,
        })
    }

    /// Classify lipids based on the orientation of acyl chains.
    /// Less reliable, fast.
    ///
    /// ## Parameters
    /// - `heads`: reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///            there must only be one such atom/bead per lipid molecule.
    /// - `methyls`: reference atoms identifying methyl groups of lipid tails, i.e., the ends of lipid tails;
    ///              there should be only one such atom/bead per one acyl chain in the molecule (e.g., two for lipids with two acyl chains).
    pub fn individual(heads: &str, methyls: &str) -> LeafletClassification {
        Self::Individual(IndividualParams {
            heads: heads.to_string(),
            methyls: methyls.to_string(),
        })
    }

    /// Prepare the system for leaflet classification.
    pub(crate) fn prepare_system(&self, system: &mut System) -> Result<(), TopologyError> {
        match self {
            Self::Global(params) => {
                create_group(system, "Membrane", &params.membrane)?;
                create_group(system, "Heads", &params.heads)?;
            }
            Self::Local(params) => {
                create_group(system, "Membrane", &params.membrane)?;
                create_group(system, "Heads", &params.heads)?;
            }
            Self::Individual(params) => {
                create_group(system, "Heads", &params.heads)?;
                create_group(system, "Methyls", &params.methyls)?;
            }
        }

        Ok(())
    }

    /// Returns a radius of the cylinder for the calculation of local membrane center of geometry, if the method is Local.
    /// Otherwise, returns None.
    pub(crate) fn get_radius(&self) -> Option<f32> {
        match self {
            Self::Local(x) => Some(x.radius),
            _ => None,
        }
    }
}

/// Based on the global membrane center of geometry; useful for disrupted membranes; fast.
#[derive(Debug, Clone)]
pub(crate) struct GlobalParams {
    /// Selection of all lipids forming the membrane.
    membrane: String,
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    heads: String,
}

/// Parameters for classification of lipids.
/// Based on the local membrane center of geometry; useful for curved membranes; slow.
#[derive(Debug, Clone)]
pub(crate) struct LocalParams {
    /// Selection of all lipids forming the membrane.
    membrane: String,
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    heads: String,
    /// Radius of a cylinder for the calculation of local membrane center of geometry (in nm).
    radius: f32,
}

/// Parameters for classification of lipids.
/// Based on the orientation of the lipid tails; less reliable; fast.
#[derive(Debug, Clone)]
pub(crate) struct IndividualParams {
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    heads: String,
    /// Reference atoms identifying methyl groups of lipid tails, i.e., the ends of lipid tails.
    /// There should be only one such atom/bead per one acyl chain in the molecule (e.g., two for lipids with two acyl chains).
    methyls: String,
}

#[cfg(test)]
mod tests {
    use crate::auxiliary::macros::group_name;

    use super::*;

    #[test]
    fn test_prepare_system_global() {
        let mut system = System::from_file("tests/files/pepg_cg.tpr").unwrap();
        let classifier = LeafletClassification::global("@membrane", "name PO4");

        classifier.prepare_system(&mut system).unwrap();

        assert!(system.group_exists(group_name!("Membrane")));
        assert!(system.group_exists(group_name!("Heads")));
    }

    #[test]
    fn test_prepare_system_local() {
        let mut system = System::from_file("tests/files/pepg_cg.tpr").unwrap();
        let classifier = LeafletClassification::local("@membrane", "name PO4", 2.5);

        classifier.prepare_system(&mut system).unwrap();

        assert!(system.group_exists(group_name!("Membrane")));
        assert!(system.group_exists(group_name!("Heads")));
    }

    #[test]
    fn test_prepare_system_individual() {
        let mut system = System::from_file("tests/files/pepg_cg.tpr").unwrap();
        let classifier = LeafletClassification::individual("name PO4", "name C4A C4B");

        classifier.prepare_system(&mut system).unwrap();

        assert!(system.group_exists(group_name!("Heads")));
        assert!(system.group_exists(group_name!("Methyls")));
    }
}
