// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the `EstimateError` structure and its methods.

use std::fmt;

use crate::errors::ErrorEstimationError;
use getset::{CopyGetters, Getters};
use serde::{
    de::{self, MapAccess, Visitor},
    Deserialize, Deserializer,
};

/// Default number of blocks to use for error estimation.
const DEFAULT_N_BLOCKS: usize = 5;

/// Parameters for estimating the error of the analysis.
#[derive(Debug, Clone, Getters, CopyGetters)]
pub struct EstimateError {
    /// Number of blocks to divide the trajectory into for error estimation.
    /// Default value is 5, and it is recommended not to change this value.
    #[getset(get_copy = "pub")]
    n_blocks: usize,

    /// Optional filename pattern for the output XVG files where convergence analysis will be written.
    /// A separate XVG file will be generated for each detected molecule type, with the molecule
    /// type name appended to the pattern.
    ///
    /// Example: 'convergence.xvg' may become 'convergence_POPC.xvg'.
    /// Default is None => no such data will be written.
    output_convergence: Option<String>,
}

impl Default for EstimateError {
    /// Default parameters for the error estimation.
    fn default() -> Self {
        Self {
            n_blocks: DEFAULT_N_BLOCKS,
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
                let mut n_blocks = None;
                let mut output_convergence = None;

                while let Some(key) = map.next_key::<String>()? {
                    match key.as_str() {
                        "n_blocks" => n_blocks = Some(map.next_value()?),
                        "output_convergence" | "output" => {
                            output_convergence = Some(map.next_value()?)
                        }
                        _ => {
                            return Err(de::Error::unknown_field(
                                &key,
                                &["n_blocks", "output_convergence", "output"],
                            ))
                        }
                    }
                }

                Ok(EstimateError {
                    n_blocks: n_blocks.unwrap_or(DEFAULT_N_BLOCKS),
                    output_convergence,
                })
            }
        }

        deserializer.deserialize_any(EstimateErrorVisitor)
    }
}

impl EstimateError {
    /// Specify parameters for the error estimation.
    /// If any of them is `None`, the default value is used.
    pub fn new(
        n_blocks: Option<usize>,
        output: Option<&str>,
    ) -> Result<Self, ErrorEstimationError> {
        if let Some(n) = n_blocks {
            if n < 2 {
                return Err(ErrorEstimationError::NotEnoughBlocks(n));
            }
        }

        Ok(Self {
            n_blocks: n_blocks.unwrap_or(DEFAULT_N_BLOCKS),
            output_convergence: output.map(|s| s.to_string()),
        })
    }

    /// Log basic info about the error estimation.
    pub(crate) fn info(&self) {
        log::info!("Will estimate error using {} blocks.", self.n_blocks());
        if let Some(output) = &self.output_convergence {
            log::info!("Will write convergence data into a file '{}'.", output);
        }
    }

    /// Check that the parameters of the error estimation are valid.
    pub(super) fn validate(&self) -> Result<(), ErrorEstimationError> {
        if self.n_blocks < 2 {
            return Err(ErrorEstimationError::NotEnoughBlocks(self.n_blocks));
        }

        Ok(())
    }

    pub fn output_convergence(&self) -> Option<&str> {
        self.output_convergence.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_error_default() {
        let ee = EstimateError::default();
        assert_eq!(ee.n_blocks, DEFAULT_N_BLOCKS);
        assert_eq!(ee.output_convergence, None);
    }

    #[test]
    fn test_estimate_error_new() {
        let ee = EstimateError::new(Some(4), Some("convergence.xvg")).unwrap();
        assert_eq!(ee.n_blocks, 4);
        assert_eq!(ee.output_convergence, Some("convergence.xvg".to_string()));

        let ee = EstimateError::new(None, None).unwrap();
        assert_eq!(ee.n_blocks, DEFAULT_N_BLOCKS);
        assert_eq!(ee.output_convergence, None);

        let ee = EstimateError::new(Some(10), None).unwrap();
        assert_eq!(ee.n_blocks, 10);
        assert_eq!(ee.output_convergence, None);

        let ee = EstimateError::new(None, Some("convergence.xvg")).unwrap();
        assert_eq!(ee.n_blocks, DEFAULT_N_BLOCKS);
        assert_eq!(ee.output_convergence, Some("convergence.xvg".to_string()));

        assert!(EstimateError::new(Some(1), None).is_err());
    }

    #[test]
    fn test_estimate_error_from_yaml() {
        let ee: EstimateError =
            serde_yaml::from_str("n_blocks: 4\noutput_convergence: convergence.xvg").unwrap();
        assert_eq!(ee.n_blocks, 4);
        assert_eq!(ee.output_convergence, Some("convergence.xvg".to_string()));

        let ee: EstimateError =
            serde_yaml::from_str("n_blocks: 4\noutput: convergence.xvg").unwrap();
        assert_eq!(ee.n_blocks, 4);
        assert_eq!(ee.output_convergence, Some("convergence.xvg".to_string()));

        let ee: EstimateError = serde_yaml::from_str("true").unwrap();
        assert_eq!(ee.n_blocks, DEFAULT_N_BLOCKS);
        assert_eq!(ee.output_convergence, None);

        let ee: EstimateError = serde_yaml::from_str("n_blocks: 10").unwrap();
        assert_eq!(ee.n_blocks, 10);
        assert_eq!(ee.output_convergence, None);

        let ee: EstimateError =
            serde_yaml::from_str("output_convergence: convergence.xvg").unwrap();
        assert_eq!(ee.n_blocks, DEFAULT_N_BLOCKS);
        assert_eq!(ee.output_convergence, Some("convergence.xvg".to_string()));

        let ee: EstimateError = serde_yaml::from_str("n_blocks: 1").unwrap();
        assert!(ee.validate().is_err());
    }
}
