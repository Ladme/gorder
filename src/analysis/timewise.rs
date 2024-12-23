// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for timewise analysis.

use crate::PANIC_MESSAGE;
use getset::Getters;
use std::ops::{Add, AddAssign};

use super::orderval::OrderValue;

#[derive(Debug, Clone, Getters)]
pub(super) struct TimeWiseData {
    /// Number of threads that the data have been collected by.
    n_threads: usize,
    /// Order parameters calculated independently for each trajectory frame.
    order: Vec<OrderValue>,
    /// Number of samples collected from each trajectory frame.
    n_samples: Vec<usize>,
}

impl Default for TimeWiseData {
    #[inline(always)]
    fn default() -> Self {
        TimeWiseData {
            n_threads: 1,
            order: Vec::new(),
            n_samples: Vec::new(),
        }
    }
}

/// Note that this is not a commutative operation.
impl Add<TimeWiseData> for TimeWiseData {
    type Output = TimeWiseData;

    #[inline(always)]
    fn add(self, rhs: TimeWiseData) -> Self::Output {
        let order = interleave_vectors(&self.order, &rhs.order, self.n_threads, rhs.n_threads);
        let n_samples = interleave_vectors(
            &self.n_samples,
            &rhs.n_samples,
            self.n_threads,
            rhs.n_threads,
        );

        TimeWiseData {
            n_threads: self.n_threads + rhs.n_threads,
            order,
            n_samples,
        }
    }
}

impl TimeWiseData {
    /// Initiate reading of the next frame.
    #[inline(always)]
    pub(super) fn next_frame(&mut self) {
        self.order.push(OrderValue::from(0.0));
        self.n_samples.push(0);
    }

    /// Estimate the calculation error from the order parameters calculated for the individual blocks.
    /// Returns standard error of the mean.
    pub(super) fn estimate_error(&self, block_size: usize) -> f32 {
        let n_blocks = self.n_blocks(block_size);
        // error estimation requires at least 2 blocks
        if n_blocks < 2 {
            return 0.0;
        }

        let mut blocks_order = vec![OrderValue::from(0.0); n_blocks];
        let mut blocks_samples = vec![0; n_blocks];

        for (i, (o, s)) in self.order.iter().zip(&self.n_samples).enumerate() {
            let block_id = i / block_size;
            if block_id < n_blocks {
                blocks_order[block_id] += *o;
                blocks_samples[block_id] += s;
            }
        }

        let orders: Vec<f32> = blocks_order
            .into_iter()
            .zip(blocks_samples)
            .map(|(o, s)| (o / s).into())
            .collect();

        let std = statistical::standard_deviation(&orders, None);

        std / (n_blocks as f32).sqrt()
    }

    /// Get the number of blocks collected for the error estimation.
    #[inline(always)]
    pub(super) fn n_blocks(&self, block_size: usize) -> usize {
        self.order.len() / block_size
    }

    /// Get the number of frames read.
    #[inline(always)]
    pub(super) fn n_frames(&self) -> usize {
        self.order.len()
    }
}

impl AddAssign<f32> for TimeWiseData {
    #[inline(always)]
    fn add_assign(&mut self, rhs: f32) {
        *self.order.last_mut().expect(PANIC_MESSAGE) += rhs;
        *self.n_samples.last_mut().expect(PANIC_MESSAGE) += 1;
    }
}

/// Merge two vectors interleaving the second with the first one.
fn interleave_vectors<T: Clone>(vec1: &[T], vec2: &[T], n: usize, m: usize) -> Vec<T> {
    let mut result = Vec::with_capacity(vec1.len() + vec2.len());
    let mut iter1 = vec1.iter();
    let mut iter2 = vec2.iter();

    loop {
        for _ in 0..n {
            if let Some(value) = iter1.next() {
                result.push(value.clone());
            }
        }

        for _ in 0..m {
            if let Some(value) = iter2.next() {
                result.push(value.clone());
            }
        }

        if iter1.len() == 0 && iter2.len() == 0 {
            break;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn timewise_merge_simple() {
        let data1 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(2.2),
                OrderValue::from(3.3),
                OrderValue::from(4.4),
            ],
            n_samples: vec![1, 2, 3, 4],
        };

        let data2 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data_sum = data1 + data2;

        let expected_order: Vec<OrderValue> = [1.1, 21.1, 2.2, 22.2, 3.3, 23.3, 4.4, 24.4]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [1, 21, 2, 22, 3, 23, 4, 24];

        assert_eq!(data_sum.n_threads, 2);
        assert_eq!(data_sum.order.len(), expected_order.len());
        assert_eq!(data_sum.n_samples.len(), expected_n_samples.len());

        for (x, y) in data_sum.order.iter().zip(expected_order.iter()) {
            assert_eq!(x, y);
        }

        for (x, y) in data_sum.n_samples.iter().zip(expected_n_samples.iter()) {
            assert_eq!(x, y);
        }
    }

    #[test]
    fn timewise_merge_incomplete() {
        let data1 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(2.2),
                OrderValue::from(3.3),
                OrderValue::from(4.4),
            ],
            n_samples: vec![1, 2, 3, 4],
        };

        let data2 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
            ],
            n_samples: vec![21, 22, 23],
        };

        let data_sum = data1 + data2;

        let expected_order: Vec<OrderValue> = [1.1, 21.1, 2.2, 22.2, 3.3, 23.3, 4.4]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [1, 21, 2, 22, 3, 23, 4];

        assert_eq!(data_sum.n_threads, 2);
        assert_eq!(data_sum.order.len(), expected_order.len());
        assert_eq!(data_sum.n_samples.len(), expected_n_samples.len());

        for (x, y) in data_sum.order.iter().zip(expected_order.iter()) {
            assert_eq!(x, y);
        }

        for (x, y) in data_sum.n_samples.iter().zip(expected_n_samples.iter()) {
            assert_eq!(x, y);
        }
    }

    #[test]
    fn timewise_merge_two_and_one() {
        let data1 = TimeWiseData {
            n_threads: 2,
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(11.1),
                OrderValue::from(2.2),
                OrderValue::from(12.2),
                OrderValue::from(3.3),
                OrderValue::from(13.3),
                OrderValue::from(4.4),
                OrderValue::from(14.4),
            ],
            n_samples: vec![1, 11, 2, 12, 3, 13, 4, 14],
        };

        let data2 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data_sum = data1 + data2;

        let expected_order: Vec<OrderValue> = [
            1.1, 11.1, 21.1, 2.2, 12.2, 22.2, 3.3, 13.3, 23.3, 4.4, 14.4, 24.4,
        ]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [1, 11, 21, 2, 12, 22, 3, 13, 23, 4, 14, 24];

        assert_eq!(data_sum.n_threads, 3);
        assert_eq!(data_sum.order.len(), expected_order.len());
        assert_eq!(data_sum.n_samples.len(), expected_n_samples.len());

        for (x, y) in data_sum.order.iter().zip(expected_order.iter()) {
            assert_eq!(x, y);
        }

        for (x, y) in data_sum.n_samples.iter().zip(expected_n_samples.iter()) {
            assert_eq!(x, y);
        }
    }

    #[test]
    fn timewise_merge_two_and_one_incomplete() {
        let data1 = TimeWiseData {
            n_threads: 2,
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(11.1),
                OrderValue::from(2.2),
                OrderValue::from(12.2),
                OrderValue::from(3.3),
            ],
            n_samples: vec![1, 11, 2, 12, 3],
        };

        let data2 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data_sum = data1 + data2;

        let expected_order: Vec<OrderValue> = [
            1.1, 11.1, 21.1, 2.2, 12.2, 22.2, 3.3, 23.3, 24.4,
        ]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [1, 11, 21, 2, 12, 22, 3, 23, 24];

        assert_eq!(data_sum.n_threads, 3);
        assert_eq!(data_sum.order.len(), expected_order.len());
        assert_eq!(data_sum.n_samples.len(), expected_n_samples.len());

        for (x, y) in data_sum.order.iter().zip(expected_order.iter()) {
            assert_eq!(x, y);
        }

        for (x, y) in data_sum.n_samples.iter().zip(expected_n_samples.iter()) {
            assert_eq!(x, y);
        }
    }

    #[test]
    fn timewise_merge_one_and_two() {
        let data1 = TimeWiseData {
            n_threads: 1,
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data2 = TimeWiseData {
            n_threads: 2,
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(11.1),
                OrderValue::from(2.2),
                OrderValue::from(12.2),
                OrderValue::from(3.3),
                OrderValue::from(13.3),
                OrderValue::from(4.4),
                OrderValue::from(14.4),
            ],
            n_samples: vec![1, 11, 2, 12, 3, 13, 4, 14],
        };

        let data_sum = data1 + data2;

        let expected_order: Vec<OrderValue> = [
            21.1, 1.1, 11.1, 22.2, 2.2, 12.2, 23.3, 3.3, 13.3, 24.4, 4.4, 14.4,
        ]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [21, 1, 11, 22, 2, 12, 23, 3, 13, 24, 4, 14];

        assert_eq!(data_sum.n_threads, 3);
        assert_eq!(data_sum.order.len(), expected_order.len());
        assert_eq!(data_sum.n_samples.len(), expected_n_samples.len());

        for (x, y) in data_sum.order.iter().zip(expected_order.iter()) {
            assert_eq!(x, y);
        }

        for (x, y) in data_sum.n_samples.iter().zip(expected_n_samples.iter()) {
            assert_eq!(x, y);
        }
    }

    #[test]
    fn timewise_merge_four_and_three() {
        let data1 = TimeWiseData {
            n_threads: 4,
            order: vec![
                OrderValue::from(1.0),
                OrderValue::from(2.0),
                OrderValue::from(3.0),
                OrderValue::from(4.0),
                OrderValue::from(1.1),
                OrderValue::from(2.1),
                OrderValue::from(3.1),
                OrderValue::from(4.1),
                OrderValue::from(1.2),
                OrderValue::from(2.2),
                OrderValue::from(3.2),
                OrderValue::from(4.2),
            ],
            n_samples: vec![10, 20, 30, 40, 11, 21, 31, 41, 12, 22, 32, 42],
        };

        let data2 = TimeWiseData {
            n_threads: 3,
            order: vec![
                OrderValue::from(5.0),
                OrderValue::from(6.0),
                OrderValue::from(7.0),
                OrderValue::from(5.1),
                OrderValue::from(6.1),
                OrderValue::from(7.1),
                OrderValue::from(5.2),
                OrderValue::from(6.2),
                OrderValue::from(7.2),
            ],
            n_samples: vec![50, 60, 70, 51, 61, 71, 52, 62, 72],
        };

        let data_sum = data1 + data2;

        let expected_order: Vec<OrderValue> = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 1.2, 2.2, 3.2,
            4.2, 5.2, 6.2, 7.2,
        ]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [
            10, 20, 30, 40, 50, 60, 70, 11, 21, 31, 41, 51, 61, 71, 12, 22, 32, 42, 52, 62, 72,
        ];

        assert_eq!(data_sum.n_threads, 7);
        assert_eq!(data_sum.order.len(), expected_order.len());
        assert_eq!(data_sum.n_samples.len(), expected_n_samples.len());

        for (x, y) in data_sum.order.iter().zip(expected_order.iter()) {
            assert_eq!(x, y);
        }

        for (x, y) in data_sum.n_samples.iter().zip(expected_n_samples.iter()) {
            assert_eq!(x, y);
        }
    }
    #[test]
    fn estimate_error() {
        let order = [
            10.0, 15.0, 18.0, 12.0, 14.0, 15.0, 16.0, 20.0, 21.0, 18.0, 9.0, 11.0, 13.0, 14.0,
            19.0, 16.0, 17.0,
        ];

        let data = TimeWiseData {
            n_threads: 3,
            order: order
                .into_iter()
                .map(|x| x.into())
                .collect::<Vec<OrderValue>>(),
            n_samples: vec![
                10, 12, 15, 11, 13, 11, 11, 17, 18, 15, 8, 10, 12, 13, 17, 14, 15,
            ],
        };

        // blocks: 1.16216, 1.1714, 1.2391, 1.1515, 1.0952 (last two numbers ignored)

        let error = data.estimate_error(3);
        assert_relative_eq!(error, 0.0230077);
    }
}
