// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the `EstimateError` structure and its methods.

use std::fmt;

use getset::{CopyGetters, Getters};
use serde::{
    de::{self, MapAccess, Visitor},
    Deserialize, Deserializer,
};

/// Parameters for estimating the error of the analysis.
#[derive(Debug, Clone, Getters, CopyGetters)]
pub struct EstimateError {
    /// Number of analyzed trajectory frames per one block.
    /// The default value is 500.
    #[getset(get_copy = "pub")]
    block_size: usize,

    /// Optional filename pattern for the output XVG files where convergence analysis will be written.
    /// A separate XVG file will be generated for each detected molecule type, with the molecule
    /// type name appended to the pattern.
    ///
    /// Example: 'convergence.xvg' may become 'convergence_POPC.xvg'.
    /// Default is None => no such data will be written.
    #[getset(get = "pub")]
    output_convergence: Option<String>,
}

impl Default for EstimateError {
    fn default() -> Self {
        Self {
            block_size: default_block_size(),
            output_convergence: None,
        }
    }
}

impl<'de> Deserialize<'de> for EstimateError {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct EstimateErrorVisitor;

        impl<'de> Visitor<'de> for EstimateErrorVisitor {
            type Value = EstimateError;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a boolean or a map for EstimateError")
            }

            fn visit_bool<E>(self, value: bool) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                if value {
                    Ok(EstimateError::default())
                } else {
                    Err(E::custom("false is not a valid value for estimate_error"))
                }
            }

            fn visit_map<M>(self, mut map: M) -> Result<Self::Value, M::Error>
            where
                M: MapAccess<'de>,
            {
                let mut block_size = None;
                let mut output_convergence = None;

                while let Some(key) = map.next_key::<String>()? {
                    match key.as_str() {
                        "block_size" => block_size = Some(map.next_value()?),
                        "output_convergence" | "output" => {
                            output_convergence = Some(map.next_value()?)
                        }
                        _ => {
                            return Err(de::Error::unknown_field(
                                &key,
                                &["block_size", "output_directory"],
                            ))
                        }
                    }
                }

                Ok(EstimateError {
                    block_size: block_size.unwrap_or_else(default_block_size),
                    output_convergence,
                })
            }
        }

        deserializer.deserialize_any(EstimateErrorVisitor)
    }
}

#[inline(always)]
fn default_block_size() -> usize {
    500
}

impl EstimateError {
    pub fn new(block_size: usize, output: Option<&str>) -> Self {
        Self {
            block_size,
            output_convergence: output.map(|s| s.to_string()),
        }
    }
}
