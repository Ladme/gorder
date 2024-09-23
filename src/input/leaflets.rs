// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use getset::{CopyGetters, Getters};
use serde::Deserialize;

/// Parameters for the classification of lipids into membrane leaflets.
#[derive(Debug, Clone, Deserialize)]
#[allow(private_interfaces)]
pub enum LeafletClassification {
    Global(GlobalParams),
    Local(LocalParams),
    Individual(IndividualParams),
}

impl LeafletClassification {
    /// Classify lipids based on the global membrane center of geometry.
    /// Generally reliable and fast. The best option when working with disrupted membranes.
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
    /// Useful for curved membranes, very slow.
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
    /// Less reliable but fast.
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
#[derive(Debug, Clone, Getters, CopyGetters, Deserialize)]
pub(crate) struct GlobalParams {
    /// Selection of all lipids forming the membrane.
    #[getset(get = "pub(crate)")]
    membrane: String,
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub(crate)")]
    heads: String,
}

/// Parameters for classification of lipids.
/// Based on the local membrane center of geometry; useful for curved membranes; slow.
#[derive(Debug, Clone, Getters, CopyGetters, Deserialize)]
pub(crate) struct LocalParams {
    /// Selection of all lipids forming the membrane.
    #[getset(get = "pub(crate)")]
    membrane: String,
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub(crate)")]
    heads: String,
    /// Radius of a cylinder for the calculation of local membrane center of geometry (in nm).
    #[getset(get_copy = "pub(crate)")]
    radius: f32,
}

/// Parameters for classification of lipids.
/// Based on the orientation of the lipid tails; less reliable; fast.
#[derive(Debug, Clone, Getters, CopyGetters, Deserialize)]
pub(crate) struct IndividualParams {
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub(crate)")]
    heads: String,
    /// Reference atoms identifying methyl groups of lipid tails, i.e., the ends of lipid tails.
    /// There should be only one such atom/bead per one acyl chain in the molecule (e.g., two for lipids with two acyl chains).
    #[getset(get = "pub(crate)")]
    methyls: String,
}
