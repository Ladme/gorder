// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This module contains structures and methods for presenting the results of the analysis.

use crate::{
    analysis::molecule::{AtomType, Bond},
    PANIC_MESSAGE,
};
use serde::{Serialize, Serializer};

pub(crate) mod aapresenter;

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
