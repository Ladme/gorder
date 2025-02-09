// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This module contains structures and methods for specifying membrane normal.

use std::fmt::{self, Display};

use getset::{CopyGetters, Getters};
use groan_rs::prelude::Dimension;
use serde::{Deserialize, Serialize};

use super::Axis;

/// Structure describing the direction of the membrane normal
/// or properties necessary for its calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum MembraneNormal {
    /// Membrane normal is oriented along a dimension of the simulation box: x, y, or z axis.
    Static(Axis),
    /// Membrane normal should be calculated dynamically for each molecule based on the shape of the membrane.
    Dynamic(DynamicNormal),
}

impl Display for MembraneNormal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Static(x) => write!(
                f,
                "Membrane normal expected to be oriented along the {} axis.",
                x
            ),
            Self::Dynamic(_) => write!(
                f,
                "Membrane normal will be dynamically calculated for each molecule."
            ),
        }
    }
}

impl From<Axis> for MembraneNormal {
    fn from(value: Axis) -> Self {
        Self::Static(value)
    }
}

impl From<DynamicNormal> for MembraneNormal {
    fn from(value: DynamicNormal) -> Self {
        Self::Dynamic(value)
    }
}

impl From<&MembraneNormal> for Dimension {
    fn from(value: &MembraneNormal) -> Self {
        match value {
            MembraneNormal::Static(axis) => axis.clone().into(),
            MembraneNormal::Dynamic(_) => todo!("Dynamic membrane normal not yet implemented."),
        }
    }
}

/// Structure describing properties of the dynamic local membrane normal calculation.
#[derive(Debug, Clone, Getters, CopyGetters, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DynamicNormal {
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub")]
    heads: String,
    /// Radius of the sphere for selecting nearby lipids for membrane normal estimation.
    /// The default value is 2 nm. The recommended value is half the membrane thickness.
    #[getset(get_copy = "pub")]
    #[serde(default = "default_dynamic_radius")]
    radius: f32,
}

impl DynamicNormal {
    /// Request a dynamic local membrane normal calculation.
    ///
    /// ## Parameters
    /// - `heads`: reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule
    /// - `radius`: radius of the sphere for selecting nearby lipids for membrane normal estimation;
    ///    the recommended value is half the membrane thickness
    pub fn new(heads: &str, radius: f32) -> DynamicNormal {
        DynamicNormal {
            heads: heads.to_owned(),
            radius,
        }
    }
}

fn default_dynamic_radius() -> f32 {
    2.0
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_deserialize() {
        match serde_yaml::from_str("!Static x").unwrap() {
            MembraneNormal::Static(axis) => assert_eq!(axis, Axis::X),
            _ => panic!("Incorrect membrane normal parsed."),
        }

        match serde_yaml::from_str("!Static y").unwrap() {
            MembraneNormal::Static(axis) => assert_eq!(axis, Axis::Y),
            _ => panic!("Incorrect membrane normal parsed."),
        }

        match serde_yaml::from_str("!Static z").unwrap() {
            MembraneNormal::Static(axis) => assert_eq!(axis, Axis::Z),
            _ => panic!("Incorrect membrane normal parsed."),
        }

        match serde_yaml::from_str("!Dynamic { heads: \"name P\" }").unwrap() {
            MembraneNormal::Dynamic(dynamic) => {
                assert_eq!(dynamic.heads(), "name P");
                assert_relative_eq!(dynamic.radius, 2.0);
            }
            _ => panic!("Incorrect membrane normal parsed."),
        }

        match serde_yaml::from_str(
            "!Dynamic
heads: \"name P\"
radius: 3.5",
        )
        .unwrap()
        {
            MembraneNormal::Dynamic(dynamic) => {
                assert_eq!(dynamic.heads(), "name P");
                assert_relative_eq!(dynamic.radius, 3.5);
            }
            _ => panic!("Incorrect membrane normal parsed."),
        }
    }

    #[test]
    fn test_serialize() {
        assert_eq!(
            serde_yaml::to_string(&MembraneNormal::Static(Axis::Z)).unwrap(),
            "!Static Z\n"
        );

        assert_eq!(
            serde_yaml::to_string(&MembraneNormal::Dynamic(DynamicNormal::new("name P", 3.0)))
                .unwrap(),
            "!Dynamic
heads: name P
radius: 3.0\n"
        );
    }
}
