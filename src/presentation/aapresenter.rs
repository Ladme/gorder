// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for presenting the results of the analysis of all-atom order parameters.

use indexmap::IndexMap;
use serde::Serialize;

use crate::analysis::molecule::AtomType;

use super::BondResults;

#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub(crate) struct AAOrderResults {
    molecules: Vec<AAMoleculeResults>,
}

#[derive(Debug, Clone, Serialize)]
struct AAMoleculeResults {
    molecule: String,
    order: IndexMap<AtomType, AAAtomResults>,
}

#[derive(Debug, Clone, Serialize)]
struct AAAtomResults {
    total: f32,
    #[serde(skip_serializing_if = "Option::is_none")]
    upper: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    lower: Option<f32>,
    bonds: IndexMap<AtomType, BondResults>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prepare_example_atom_results_1() -> AAAtomResults {
        let mut results = AAAtomResults {
            total: 0.777,
            upper: None,
            lower: None,
            bonds: IndexMap::new(),
        };

        let bond1_results = BondResults {
            total: 0.721,
            upper: None,
            lower: None,
        };

        let bond2_results = BondResults {
            total: 0.833,
            upper: None,
            lower: None,
        };

        results
            .bonds
            .insert(AtomType::new_raw(1, "POPE", "HA"), bond1_results);

        results
            .bonds
            .insert(AtomType::new_raw(3, "POPE", "HB"), bond2_results);

        results
    }

    #[test]
    fn test_serialize_aa_atom_results() {
        let results = prepare_example_atom_results_1();
        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(
            serialized,
            "total: 0.777
bonds:
  POPE HA 1:
    total: 0.721
  POPE HB 3:
    total: 0.833
"
        );
    }

    fn prepare_example_atom_results_2() -> AAAtomResults {
        let mut results = AAAtomResults {
            total: 0.777,
            upper: Some(0.743),
            lower: Some(0.811),
            bonds: IndexMap::new(),
        };

        let bond1_results = BondResults {
            total: 0.721,
            upper: Some(0.687),
            lower: Some(0.755),
        };

        let bond2_results = BondResults {
            total: 0.833,
            upper: Some(0.854),
            lower: Some(0.812),
        };

        results
            .bonds
            .insert(AtomType::new_raw(7, "POPE", "HA"), bond1_results);

        results
            .bonds
            .insert(AtomType::new_raw(11, "POPE", "HB"), bond2_results);

        results
    }

    #[test]
    fn test_serialize_aa_atom_results_leaflets() {
        let results = prepare_example_atom_results_2();
        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(
            serialized,
            "total: 0.777
upper: 0.743
lower: 0.811
bonds:
  POPE HA 7:
    total: 0.721
    upper: 0.687
    lower: 0.755
  POPE HB 11:
    total: 0.833
    upper: 0.854
    lower: 0.812
"
        );
    }

    #[test]
    fn test_serialize_aa_molecule_results() {
        let atom_results_1 = prepare_example_atom_results_1();
        let atom_type_1 = AtomType::new_raw(0, "POPE", "CA");
        let atom_results_2 = prepare_example_atom_results_2();
        let atom_type_2 = AtomType::new_raw(5, "POPE", "CB");
        let molecule_results = AAMoleculeResults {
            molecule: "POPE".to_owned(),
            order: IndexMap::from_iter(
                [(atom_type_1, atom_results_1), (atom_type_2, atom_results_2)].into_iter(),
            ),
        };

        let serialized = serde_yaml::to_string(&molecule_results).unwrap();
        assert_eq!(
            serialized,
            "molecule: POPE
order:
  POPE CA 0:
    total: 0.777
    bonds:
      POPE HA 1:
        total: 0.721
      POPE HB 3:
        total: 0.833
  POPE CB 5:
    total: 0.777
    upper: 0.743
    lower: 0.811
    bonds:
      POPE HA 7:
        total: 0.721
        upper: 0.687
        lower: 0.755
      POPE HB 11:
        total: 0.833
        upper: 0.854
        lower: 0.812
"
        )
    }

    #[test]
    fn test_serialize_aa_results() {
        let atom_results_1 = prepare_example_atom_results_1();
        let atom_type_1 = AtomType::new_raw(0, "POPE", "CA");
        let atom_results_2 = prepare_example_atom_results_2();
        let atom_type_2 = AtomType::new_raw(5, "POPE", "CB");
        let molecule_results = AAMoleculeResults {
            molecule: "POPE".to_owned(),
            order: IndexMap::from_iter(
                [(atom_type_1, atom_results_1), (atom_type_2, atom_results_2)].into_iter(),
            ),
        };
        let results = AAOrderResults {
            molecules: vec![molecule_results],
        };

        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(
            serialized,
            "- molecule: POPE
  order:
    POPE CA 0:
      total: 0.777
      bonds:
        POPE HA 1:
          total: 0.721
        POPE HB 3:
          total: 0.833
    POPE CB 5:
      total: 0.777
      upper: 0.743
      lower: 0.811
      bonds:
        POPE HA 7:
          total: 0.721
          upper: 0.687
          lower: 0.755
        POPE HB 11:
          total: 0.833
          upper: 0.854
          lower: 0.812
"
        );
    }
}
