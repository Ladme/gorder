// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Contains structures and methods for the assignment of lipids into membrane leaflets.

use std::fmt;

use getset::{CopyGetters, Getters};
use hashbrown::HashMap;
use serde::{
    de::{self, MapAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};

use crate::Leaflet;

use super::frequency::Frequency;

/// Parameters for the classification of lipids into membrane leaflets.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub enum LeafletClassification {
    Global(GlobalParams),
    Local(LocalParams),
    Individual(IndividualParams),
    Manual(ManualParams),
}

impl fmt::Display for LeafletClassification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LeafletClassification::Global(_) => write!(f, "global"),
            LeafletClassification::Local(_) => write!(f, "local"),
            LeafletClassification::Individual(_) => write!(f, "individual"),
            LeafletClassification::Manual(_) => write!(f, "manual"),
        }
    }
}

impl LeafletClassification {
    /// Classify lipids based on the global membrane center of geometry.
    /// Generally reliable and fast. The best option when working with disrupted membranes.
    ///
    /// ## Parameters
    /// - `membrane` - selection of all lipids forming the membrane
    /// - `heads` - reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule.
    pub fn global(membrane: &str, heads: &str) -> LeafletClassification {
        Self::Global(GlobalParams {
            membrane: membrane.to_string(),
            heads: heads.to_string(),
            frequency: Frequency::default(),
        })
    }

    /// Classify lipids based on the local membrane center of geometry.
    /// Useful for curved membranes, very slow.
    ///
    /// ## Parameters
    /// - `membrane` - selection of all lipids forming the membrane
    /// - `heads` - reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule.
    /// - `radius` - radius of a cylinder for the calculation of local membrane center of geometry (in nm)
    ///
    /// ## Panic
    /// Panics if the radius is not positive.
    pub fn local(membrane: &str, heads: &str, radius: f32) -> LeafletClassification {
        if radius <= 0.0 {
            panic!("Radius must be greater than 0, not `{}`.", radius);
        }

        Self::Local(LocalParams {
            membrane: membrane.to_string(),
            heads: heads.to_string(),
            radius,
            frequency: Frequency::default(),
        })
    }

    /// Classify lipids based on the orientation of acyl chains.
    /// Less reliable but fast.
    ///
    /// ## Parameters
    /// - `heads`: reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead);
    ///    there must only be one such atom/bead per lipid molecule.
    /// - `methyls`: reference atoms identifying methyl groups of lipid tails, i.e., the ends of lipid tails;
    ///    there should be only one such atom/bead per one acyl chain in the molecule (e.g., two for lipids with two acyl chains).
    #[inline(always)]
    pub fn individual(heads: &str, methyls: &str) -> LeafletClassification {
        Self::Individual(IndividualParams {
            heads: heads.to_string(),
            methyls: methyls.to_string(),
            frequency: Frequency::default(),
        })
    }

    /// Read leaflet assignment from an external yaml file.
    ///
    /// ## Parameters
    /// - `file`: path to the input yaml file containing the leaflet assignment.
    #[inline(always)]
    pub fn from_file(file: &str) -> LeafletClassification {
        Self::Manual(ManualParams::FromFile {
            file: file.to_owned(),
            frequency: Frequency::default(),
        })
    }

    /// Provide leaflet assignment as a map.
    #[inline(always)]
    pub fn from_map(assignment: HashMap<String, Vec<Vec<Leaflet>>>) -> LeafletClassification {
        Self::Manual(ManualParams::FromMap {
            assignment,
            frequency: Frequency::default(),
        })
    }

    /// Assign lipids to leaflets every N analyzed trajectory frames or only once (using the first trajectory frame).
    /// (Note that this is 'analyzed trajectory frames' - if you skip some frames using `step`,
    /// they will not be counted here.)
    #[inline(always)]
    pub fn with_frequency(mut self, frequency: Frequency) -> Self {
        match &mut self {
            LeafletClassification::Global(x) => x.frequency = frequency,
            LeafletClassification::Local(x) => x.frequency = frequency,
            LeafletClassification::Individual(x) => x.frequency = frequency,
            LeafletClassification::Manual(x) => match x {
                ManualParams::FromFile {
                    file: _,
                    frequency: x,
                } => *x = frequency,
                ManualParams::FromMap {
                    assignment: _,
                    frequency: x,
                } => *x = frequency,
            },
        }

        self
    }

    /// Get the frequency of the analysis.
    #[inline(always)]
    pub fn get_frequency(&self) -> Frequency {
        match self {
            LeafletClassification::Global(x) => x.frequency(),
            LeafletClassification::Local(x) => x.frequency(),
            LeafletClassification::Individual(x) => x.frequency(),
            LeafletClassification::Manual(x) => x.frequency(),
        }
    }

    /// Returns a radius of the cylinder for the calculation of local membrane center of geometry, if the method is Local.
    /// Otherwise, returns None.
    #[inline(always)]
    pub(crate) fn get_radius(&self) -> Option<f32> {
        match self {
            Self::Local(x) => Some(x.radius),
            _ => None,
        }
    }
}

/// Based on the global membrane center of geometry; useful for disrupted membranes; fast.
#[derive(Debug, Clone, Getters, CopyGetters, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GlobalParams {
    /// Selection of all lipids forming the membrane.
    #[getset(get = "pub")]
    membrane: String,
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub")]
    heads: String,
    /// Frequency of leaflet assignment.
    #[getset(get_copy = "pub")]
    #[serde(default)]
    frequency: Frequency,
}

/// Parameters for classification of lipids.
/// Based on the local membrane center of geometry; useful for curved membranes; slow.
#[derive(Debug, Clone, Getters, CopyGetters, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LocalParams {
    /// Selection of all lipids forming the membrane.
    #[getset(get = "pub")]
    membrane: String,
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub")]
    heads: String,
    /// Radius of a cylinder for the calculation of local membrane center of geometry (in nm).
    #[getset(get_copy = "pub")]
    #[serde(deserialize_with = "validate_radius")]
    radius: f32,
    /// Frequency of leaflet assignment.
    #[getset(get_copy = "pub")]
    #[serde(default)]
    frequency: Frequency,
}

fn validate_radius<'de, D>(deserializer: D) -> Result<f32, D::Error>
where
    D: Deserializer<'de>,
{
    let radius = f32::deserialize(deserializer)?;
    if radius <= 0.0 {
        Err(de::Error::custom("radius must be greater than 0"))
    } else {
        Ok(radius)
    }
}

/// Parameters for classification of lipids.
/// Based on the orientation of the lipid tails; less reliable; fast.
#[derive(Debug, Clone, Getters, CopyGetters, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct IndividualParams {
    /// Reference atoms identifying lipid headgroups (usually a phosphorus atom or a phosphate bead).
    /// There must only be one such atom/bead per lipid molecule.
    #[getset(get = "pub")]
    heads: String,
    /// Reference atoms identifying methyl groups of lipid tails, i.e., the ends of lipid tails.
    /// There should be only one such atom/bead per one acyl chain in the molecule (e.g., two for lipids with two acyl chains).
    #[getset(get = "pub")]
    methyls: String,
    /// Frequency of leaflet assignment.
    #[getset(get_copy = "pub")]
    #[serde(default)]
    frequency: Frequency,
}

#[derive(Debug, Clone, Serialize)]
pub enum ManualParams {
    FromFile {
        file: String,
        frequency: Frequency,
    },
    FromMap {
        assignment: HashMap<String, Vec<Vec<Leaflet>>>,
        frequency: Frequency,
    },
}

impl<'de> Deserialize<'de> for ManualParams {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        enum Field {
            File,
            Frequency,
            Assignment,
        }

        impl<'de> Deserialize<'de> for Field {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where
                D: Deserializer<'de>,
            {
                struct FieldVisitor;

                impl<'de> Visitor<'de> for FieldVisitor {
                    type Value = Field;

                    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                        formatter.write_str("`file`, `frequency`, or `assignment`")
                    }

                    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
                    where
                        E: de::Error,
                    {
                        match value {
                            "file" => Ok(Field::File),
                            "frequency" => Ok(Field::Frequency),
                            "assignment" => Ok(Field::Assignment),
                            _ => Err(de::Error::unknown_field(
                                value,
                                &["file", "frequency", "assignment"],
                            )),
                        }
                    }
                }

                deserializer.deserialize_identifier(FieldVisitor)
            }
        }

        struct ManualParamsVisitor;

        impl<'de> Visitor<'de> for ManualParamsVisitor {
            type Value = ManualParams;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a string, a map with `file` or `assignment`, or a plain HashMap<String, Vec<Vec<Leaflet>>>")
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                // a simple string is interpreted as a file
                Ok(ManualParams::FromFile {
                    file: value.to_string(),
                    frequency: Frequency::default(),
                })
            }

            fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
            where
                A: MapAccess<'de>,
            {
                let mut file = None;
                let mut frequency = None;
                let mut assignment = None;
                let mut fallback_map: Option<HashMap<String, Vec<Vec<Leaflet>>>> = None;

                while let Some(key) = map.next_key::<String>()? {
                    match key.as_str() {
                        "file" => {
                            if file.is_some() {
                                return Err(de::Error::duplicate_field("file"));
                            }
                            file = Some(map.next_value()?);
                        }
                        "frequency" => {
                            if frequency.is_some() {
                                return Err(de::Error::duplicate_field("frequency"));
                            }
                            frequency = Some(map.next_value()?);
                        }
                        "assignment" => {
                            if assignment.is_some() {
                                return Err(de::Error::duplicate_field("assignment"));
                            }
                            assignment = Some(map.next_value()?);
                        }
                        _ => {
                            // collect entries for fallback as HashMap
                            if fallback_map.is_none() {
                                fallback_map = Some(HashMap::new());
                            }
                            if let Some(ref mut fallback) = fallback_map {
                                let value: Vec<Vec<Leaflet>> = map.next_value()?;
                                fallback.insert(key, value);
                            }
                        }
                    }
                }

                // handle explicit cases first
                if let Some(file) = file {
                    return Ok(ManualParams::FromFile {
                        file,
                        frequency: frequency.unwrap_or(Frequency::default()),
                    });
                }

                if let Some(assignment) = assignment {
                    return Ok(ManualParams::FromMap {
                        assignment,
                        frequency: frequency.unwrap_or(Frequency::default()),
                    });
                }

                // fallback to parsing as a plain HashMap
                if let Some(fallback_map) = fallback_map {
                    return Ok(ManualParams::FromMap {
                        assignment: fallback_map,
                        frequency: frequency.unwrap_or(Frequency::default()),
                    });
                }

                Err(de::Error::custom("Invalid structure for ManualParams"))
            }
        }

        deserializer.deserialize_any(ManualParamsVisitor)
    }
}

impl ManualParams {
    pub(crate) fn frequency(&self) -> Frequency {
        match self {
            Self::FromFile { file: _, frequency } => *frequency,
            Self::FromMap {
                assignment: _,
                frequency,
            } => *frequency,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn local_radius_negative_fail() {
        let _classification = LeafletClassification::local("@membrane", "name P", -1.4);
    }

    #[test]
    #[should_panic]
    fn local_radius_zero_fail() {
        let _classification = LeafletClassification::local("@membrane", "name P", 0.0);
    }

    #[test]
    fn test_parse_manual_file_only() {
        let string = "!Manual \"leaflets.yaml\"";
        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromFile { file, frequency }) => {
                assert_eq!(file, "leaflets.yaml");
                assert_eq!(frequency, Frequency::every(1).unwrap());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }

        let string = "!Manual leaflets.yaml";
        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromFile { file, frequency }) => {
                assert_eq!(file, "leaflets.yaml");
                assert_eq!(frequency, Frequency::every(1).unwrap());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }
    }

    #[test]
    fn test_parse_manual_file_explicit() {
        let string = "!Manual { file: leaflets.yaml }";
        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromFile { file, frequency }) => {
                assert_eq!(file, "leaflets.yaml");
                assert_eq!(frequency, Frequency::every(1).unwrap());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }
    }

    #[test]
    fn test_parse_manual_file_full() {
        let string = "!Manual { file: leaflets.yaml, frequency: !Once }";
        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromFile { file, frequency }) => {
                assert_eq!(file, "leaflets.yaml");
                assert_eq!(frequency, Frequency::once());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }
    }

    fn compare_assignment(assignment: &HashMap<String, Vec<Vec<Leaflet>>>) {
        let expected_popc = [
            [
                Leaflet::Upper,
                Leaflet::Upper,
                Leaflet::Upper,
                Leaflet::Lower,
                Leaflet::Lower,
                Leaflet::Lower,
            ],
            [
                Leaflet::Upper,
                Leaflet::Lower,
                Leaflet::Upper,
                Leaflet::Upper,
                Leaflet::Lower,
                Leaflet::Upper,
            ],
        ];

        let expected_pope = [
            [
                Leaflet::Lower,
                Leaflet::Lower,
                Leaflet::Upper,
                Leaflet::Lower,
                Leaflet::Lower,
                Leaflet::Upper,
            ],
            [
                Leaflet::Upper,
                Leaflet::Lower,
                Leaflet::Lower,
                Leaflet::Upper,
                Leaflet::Lower,
                Leaflet::Lower,
            ],
        ];

        for (mol_type, expected) in ["POPC", "POPE"]
            .into_iter()
            .zip([expected_popc, expected_pope].into_iter())
        {
            let data = assignment.get(mol_type).unwrap();
            assert_eq!(data.len(), expected.len());
            for (frame, frame_expected) in data.iter().zip(expected.iter()) {
                assert_eq!(frame.len(), frame_expected.len());
                for (mol, mol_expected) in frame.iter().zip(frame_expected.iter()) {
                    assert_eq!(mol, mol_expected);
                }
            }
        }
    }

    #[test]
    fn test_parse_manual_map_only() {
        let string = "!Manual
POPC: [[1, 1, 1, 0, 0, 0], [1, 0, 1, 1, 0, 1]]
POPE:
  - [Lower, Lower, Upper, Lower, Lower, Upper]
  - [1, 0, 0, 1, 0, 0]";

        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromMap {
                assignment,
                frequency,
            }) => {
                compare_assignment(&assignment);
                assert_eq!(frequency, Frequency::every(1).unwrap());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }
    }

    #[test]
    fn test_parse_manual_map_explicit() {
        let string = "!Manual
assignment:
  POPC: [[1, 1, 1, 0, 0, 0], [1, 0, 1, 1, 0, 1]]
  POPE:
    - [Lower, Lower, Upper, Lower, Lower, Upper]
    - [1, 0, 0, 1, 0, 0]";

        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromMap {
                assignment,
                frequency,
            }) => {
                compare_assignment(&assignment);
                assert_eq!(frequency, Frequency::every(1).unwrap());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }
    }

    #[test]
    fn test_parse_manual_map_full() {
        let string = "!Manual
assignment:
  POPC: [[1, 1, 1, 0, 0, 0], [1, 0, 1, 1, 0, 1]]
  POPE:
    - [Lower, Lower, Upper, Lower, Lower, Upper]
    - [1, 0, 0, 1, 0, 0]
frequency: !Once";

        match serde_yaml::from_str(string).unwrap() {
            LeafletClassification::Manual(ManualParams::FromMap {
                assignment,
                frequency,
            }) => {
                compare_assignment(&assignment);
                assert_eq!(frequency, Frequency::once());
            }
            _ => panic!("Invalid leaflet classification returned."),
        }
    }

    #[test]
    fn test_parse_manual_fail_1() {
        let string = "!Manual 7";
        let params: Result<LeafletClassification, _> = serde_yaml::from_str(string);
        assert!(params.is_err());
    }

    #[test]
    fn test_parse_manual_fail_2() {
        let string = "!Manual
POPC: [[1, 1, 1, 0, 0, 0], [1, 0, 1, 1, 0, 1]]
POPE:
  - [Lower, Lower, Upp3r, Lower, Lower, Upper]
  - [1, 0, 0, 1, 0, 0]";
        let params: Result<LeafletClassification, _> = serde_yaml::from_str(string);
        assert!(params.is_err());
    }

    #[test]
    fn test_parse_manual_fail_3() {
        let string = "!Manual
POPC: [[1, 1, -1, 0, 0, 0], [1, 0, 1, 1, 0, 1]]
POPE:
  - [Lower, Lower, Upper, Lower, Lower, Upper]
  - [1, 0, 0, 1, 0, 0]";
        let params: Result<LeafletClassification, _> = serde_yaml::from_str(string);
        assert!(params.is_err());
    }

    #[test]
    fn test_parse_manual_fail_4() {
        let string = "!Manual
POPC: leaflets.yaml
POPE:
  - [Lower, Lower, Upper, Lower, Lower, Upper]
  - [1, 0, 0, 1, 0, 0]";
        let params: Result<LeafletClassification, _> = serde_yaml::from_str(string);
        assert!(params.is_err());
    }

    #[test]
    fn test_parse_manual_fail_5() {
        let string = "!Manual
assignment: leaflets.yaml
frequency: !Once";

        let params: Result<LeafletClassification, _> = serde_yaml::from_str(string);
        assert!(params.is_err());
    }

    #[test]
    fn test_parse_manual_fail_6() {
        let string = "!Manual
file: leaflets.yaml
frequency: some";

        let params: Result<LeafletClassification, _> = serde_yaml::from_str(string);
        assert!(params.is_err());
    }

    #[test]
    fn test_parse_local_negative_radius_fail() {
        let string = "!Local
membrane: \"@membrane\"
heads: name P
radius: -1.5";

        let result: Result<LeafletClassification, _> = serde_yaml::from_str(string);

        match result {
            Err(e) => assert_eq!(e.to_string(), "radius must be greater than 0"),
            Ok(_) => panic!("Should have failed."),
        }
    }

    #[test]
    fn test_parse_local_zero_radius_fail() {
        let string = "!Local
membrane: \"@membrane\"
heads: name P
radius: 0.0";

        let result: Result<LeafletClassification, _> = serde_yaml::from_str(string);

        match result {
            Err(e) => assert_eq!(e.to_string(), "radius must be greater than 0"),
            Ok(_) => panic!("Should have failed."),
        }
    }
}
