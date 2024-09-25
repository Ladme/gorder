// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for presenting the results of the analysis of all-atom order parameters.

use std::{fs::File, io::BufWriter, path::Path};

use indexmap::IndexMap;
use serde::Serialize;
use std::io::Write;

use crate::{
    analysis::{
        molecule::{AtomType, Bond, MoleculeType},
        topology::SystemTopology,
    },
    errors::WriteError,
    PANIC_MESSAGE,
};

use super::BondResults;

#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub(crate) struct AAOrderResults {
    molecules: Vec<AAMoleculeResults>,
}

impl From<SystemTopology> for AAOrderResults {
    fn from(value: SystemTopology) -> Self {
        AAOrderResults {
            molecules: value
                .molecules()
                .iter()
                .map(|m| AAMoleculeResults::from(m))
                .collect(),
        }
    }
}

impl AAOrderResults {
    /// Write the results of the analysis into a new file with the specified name.
    pub(crate) fn write_yaml(
        &self,
        filename: &impl AsRef<Path>,
        input_structure: &str,
        input_trajectory: &str,
        overwrite: bool,
    ) -> Result<(), WriteError> {
        if filename.as_ref().exists() {
            if !overwrite {
                log::warn!(
                    "Output yaml file '{}' already exists. Backing it up.",
                    filename.as_ref().to_str().expect(PANIC_MESSAGE)
                );
                backitup::backup(filename)
                    .map_err(|_| WriteError::CouldNotBackupFile(Box::from(filename.as_ref())))?;
            } else {
                log::warn!(
                    "Output yaml file '{}' already exists. It will be overwritten as requested.",
                    filename.as_ref().to_str().expect(PANIC_MESSAGE)
                );
            }
        }

        let file = File::create(filename)
            .map_err(|_| WriteError::CouldNotCreateFile(Box::from(filename.as_ref())))?;

        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "# Order parameters calculated with 'gorder v{}' using structure file '{}' and trajectory file '{}'.",
            crate::GORDER_VERSION, input_structure, input_trajectory
        )
        .map_err(|_| WriteError::CouldNotWriteYaml(Box::from(filename.as_ref())))?;

        serde_yaml::to_writer(writer, self)
            .map_err(|_| WriteError::CouldNotWriteYaml(Box::from(filename.as_ref())))?;

        Ok(())
    }
}

#[derive(Debug, Clone, Serialize)]
struct AAMoleculeResults {
    molecule: String,
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    order: IndexMap<AtomType, AAAtomResults>,
}

impl From<&MoleculeType> for AAMoleculeResults {
    fn from(value: &MoleculeType) -> Self {
        let mut order = IndexMap::new();
        for heavy_atom in value.order_atoms().atoms() {
            let mut relevant_bonds = Vec::new();

            for bond in value.order_bonds().bonds() {
                if bond.contains(heavy_atom) {
                    relevant_bonds.push(bond);
                }
            }

            let results = AAAtomResults::new(&relevant_bonds, heavy_atom);

            // only include atoms that have some data associated with them
            if !results.is_empty() {
                order.insert(heavy_atom.to_owned(), results);
            }
        }

        AAMoleculeResults {
            molecule: value.name().to_owned(),
            order,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
struct AAAtomResults {
    // 'total' is set to `None`, if there are no bonds associated with the atom
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "super::round_serialize_option_f32")]
    total: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "super::round_serialize_option_f32")]
    upper: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "super::round_serialize_option_f32")]
    lower: Option<f32>,
    #[serde(skip_serializing_if = "IndexMap::is_empty")]
    bonds: IndexMap<AtomType, BondResults>,
}

impl AAAtomResults {
    fn new(bonds: &[&Bond], heavy_atom: &AtomType) -> Self {
        let mut totals = Vec::new();
        let mut uppers = Vec::new();
        let mut lowers = Vec::new();

        let mut results = IndexMap::new();
        for bond in bonds {
            let bond_results = BondResults::from(*bond);
            totals.push(bond_results.total.clone());
            uppers.push(bond_results.upper.clone());
            lowers.push(bond_results.lower.clone());

            let hydrogen = bond.get_other_atom(heavy_atom).unwrap_or_else(|| {
                panic!(
                    "FATAL GORDER ERROR | AAAtomResults::new | Heavy atom not part of bond. {}",
                    PANIC_MESSAGE
                )
            });
            results.insert(hydrogen.to_owned(), bond_results);
        }

        if results.is_empty() {
            AAAtomResults {
                total: None,
                upper: None,
                lower: None,
                bonds: results,
            }
        } else {
            let total = Some(totals.iter().sum::<f32>() / totals.len() as f32);
            let upper = average_of_some(&uppers);
            let lower = average_of_some(&lowers);

            AAAtomResults {
                total,
                upper,
                lower,
                bonds: results,
            }
        }
    }

    fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }
}

fn average_of_some(values: &[Option<f32>]) -> Option<f32> {
    let (sum, count) = values
        .iter()
        .filter_map(|&x| x.filter(|&n| n.is_finite())) // filter out None and NaN or infinite values
        .fold((0.0, 0), |(acc, cnt), x| (acc + x, cnt + 1));

    if count > 0 {
        Some(sum / count as f32)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use groan_rs::prelude::{Atom, SimBox};

    use super::*;

    fn prepare_example_atom_results_1() -> AAAtomResults {
        let mut results = AAAtomResults {
            total: Some(0.777),
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
  POPE HA (1):
    total: 0.721
  POPE HB (3):
    total: 0.833
"
        );
    }

    #[test]
    fn test_serialize_aa_atom_results_empty() {
        let results = AAAtomResults {
            total: None,
            upper: None,
            lower: None,
            bonds: IndexMap::new(),
        };

        let serialized = serde_yaml::to_string(&results).unwrap();
        assert_eq!(serialized, "{}\n");
    }

    fn prepare_example_atom_results_2() -> AAAtomResults {
        let mut results = AAAtomResults {
            total: Some(0.777),
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
  POPE HA (7):
    total: 0.721
    upper: 0.687
    lower: 0.755
  POPE HB (11):
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
  POPE CA (0):
    total: 0.777
    bonds:
      POPE HA (1):
        total: 0.721
      POPE HB (3):
        total: 0.833
  POPE CB (5):
    total: 0.777
    upper: 0.743
    lower: 0.811
    bonds:
      POPE HA (7):
        total: 0.721
        upper: 0.687
        lower: 0.755
      POPE HB (11):
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
    POPE CA (0):
      total: 0.777
      bonds:
        POPE HA (1):
          total: 0.721
        POPE HB (3):
          total: 0.833
    POPE CB (5):
      total: 0.777
      upper: 0.743
      lower: 0.811
      bonds:
        POPE HA (7):
          total: 0.721
          upper: 0.687
          lower: 0.755
        POPE HB (11):
          total: 0.833
          upper: 0.854
          lower: 0.812
"
        );
    }

    #[test]
    fn test_aaatom_results_new() {
        let heavy_atom = Atom::new(1, "POPC", 1, "C2");
        let hydrogen1 = Atom::new(1, "POPC", 3, "HA");
        let hydrogen2 = Atom::new(1, "POPC", 4, "HB");

        let mut bond1 = Bond::new(
            0,
            &heavy_atom,
            2,
            &hydrogen1,
            0,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        );
        let mut bond2 = Bond::new(
            0,
            &heavy_atom,
            3,
            &hydrogen2,
            0,
            false,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        );

        *bond1.total_mut() += 0.234;
        *bond1.total_mut() += 0.176;
        *bond1.total_mut() += 0.112;
        // average: 0.174

        *bond2.total_mut() += 0.112;
        *bond2.total_mut() += 0.245;
        *bond2.total_mut() += -0.013;
        // average: 0.11466666

        let heavy_atom_type = AtomType::new_raw(0, "POPC", "C2");
        let hydrogen1_atom_type = AtomType::new_raw(2, "POPC", "HA");
        let hydrogen2_atom_type = AtomType::new_raw(3, "POPC", "HB");

        let results = AAAtomResults::new(&[&bond1, &bond2], &heavy_atom_type);

        assert_eq!(results.bonds.len(), 2);

        let bond1_results = results.bonds.get(&hydrogen1_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, 0.174);
        assert!(bond1_results.upper.is_none());
        assert!(bond1_results.lower.is_none());

        let bond1_results = results.bonds.get(&hydrogen2_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, 0.11466666);
        assert!(bond1_results.upper.is_none());
        assert!(bond1_results.lower.is_none());

        assert_relative_eq!(results.total.unwrap(), 0.14433333);
        assert!(results.upper.is_none());
        assert!(results.lower.is_none());
    }

    #[test]
    fn test_aaatom_results_new_leaflets() {
        let heavy_atom = Atom::new(1, "POPC", 1, "C2");
        let hydrogen1 = Atom::new(1, "POPC", 3, "HA");
        let hydrogen2 = Atom::new(1, "POPC", 4, "HB");

        let mut bond1 = Bond::new(
            0,
            &heavy_atom,
            2,
            &hydrogen1,
            0,
            true,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        );
        let mut bond2 = Bond::new(
            0,
            &heavy_atom,
            3,
            &hydrogen2,
            0,
            true,
            None,
            1,
            &SimBox::from([10.0, 10.0, 10.0]),
        );

        *bond1.total_mut() += 0.234;
        *bond1.total_mut() += 0.176;
        *bond1.total_mut() += 0.112;
        // average: 0.174

        *bond1.upper_mut().as_mut().unwrap() += 0.115;
        *bond1.upper_mut().as_mut().unwrap() += 0.325;
        *bond1.upper_mut().as_mut().unwrap() += 0.008;
        // average: 0.14933333

        *bond2.total_mut() += 0.112;
        *bond2.total_mut() += 0.245;
        *bond2.total_mut() += -0.013;
        // average: 0.11466666

        *bond2.lower_mut().as_mut().unwrap() += 0.332;
        *bond2.lower_mut().as_mut().unwrap() += 0.087;
        *bond2.lower_mut().as_mut().unwrap() += 0.125;
        // average: 0.18133333

        let heavy_atom_type = AtomType::new_raw(0, "POPC", "C2");
        let hydrogen1_atom_type = AtomType::new_raw(2, "POPC", "HA");
        let hydrogen2_atom_type = AtomType::new_raw(3, "POPC", "HB");

        let results = AAAtomResults::new(&[&bond1, &bond2], &heavy_atom_type);

        assert_eq!(results.bonds.len(), 2);

        let bond1_results = results.bonds.get(&hydrogen1_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, 0.174);
        assert_relative_eq!(bond1_results.upper.unwrap(), 0.14933333);
        assert!(bond1_results.lower.unwrap().is_nan()); // no values for the lower leaflet

        let bond1_results = results.bonds.get(&hydrogen2_atom_type).unwrap();
        assert_relative_eq!(bond1_results.total, 0.11466666);
        assert!(bond1_results.upper.unwrap().is_nan()); // no values for the upper leaflet
        assert_relative_eq!(bond1_results.lower.unwrap(), 0.18133333);

        assert_relative_eq!(results.total.unwrap(), 0.14433333);
        assert_relative_eq!(results.upper.unwrap(), 0.14933333);
        assert_relative_eq!(results.lower.unwrap(), 0.18133333);
    }
}
