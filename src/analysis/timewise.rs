// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for timewise analysis.

use std::ops::{Add, AddAssign};

use crate::PANIC_MESSAGE;

use super::orderval::OrderValue;

#[derive(Debug, Clone)]
pub(super) struct TimeWiseData {
    /// Indices of the threads these data were collected by.
    thread_ids: Vec<usize>,
    /// Order parameters calculated independently for each trajectory frame.
    order: Vec<OrderValue>,
    /// Number of samples collected from each trajectory frame.
    n_samples: Vec<usize>,
}

impl Add<TimeWiseData> for TimeWiseData {
    type Output = TimeWiseData;
    fn add(self, rhs: TimeWiseData) -> Self::Output {
        let mut lhs_order = split_vector_by_n(&self.order, self.thread_ids.len());
        let mut lhs_samples = split_vector_by_n(&self.n_samples, self.thread_ids.len());

        let rhs_order = split_vector_by_n(&rhs.order, rhs.thread_ids.len());
        let rhs_samples = split_vector_by_n(&rhs.n_samples, rhs.thread_ids.len());

        lhs_order.extend(rhs_order.into_iter());
        lhs_samples.extend(rhs_samples.into_iter());

        let mut threads = self
            .thread_ids
            .iter()
            .chain(rhs.thread_ids.iter())
            .cloned()
            .collect::<Vec<usize>>();

        let order = order_and_merge(lhs_order, &threads);
        let n_samples = order_and_merge(lhs_samples, &threads);

        threads.sort_unstable();

        TimeWiseData {
            thread_ids: threads,
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
        let n_blocks = self.order.len() / block_size;

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
}

impl AddAssign<f32> for TimeWiseData {
    #[inline(always)]
    fn add_assign(&mut self, rhs: f32) {
        *self.order.last_mut().expect(PANIC_MESSAGE) += rhs;
        *self.n_samples.last_mut().expect(PANIC_MESSAGE) += 1;
    }
}

/// Splits vector into multiple vectors based on the provided number.
fn split_vector_by_n<T: Clone>(vec: &[T], n: usize) -> Vec<Vec<T>> {
    let mut result: Vec<Vec<T>> = vec![Vec::new(); n];

    for (index, value) in vec.iter().enumerate() {
        result[index % n].push(value.clone());
    }

    result
}

/// Joins multiple vectors together.
fn merge_vectors<T: Clone>(vectors: Vec<Vec<T>>) -> Vec<T> {
    let mut result = Vec::new();
    let max_length = vectors.iter().map(|v| v.len()).max().unwrap_or(0);

    for i in 0..max_length {
        for vec in &vectors {
            if let Some(value) = vec.get(i) {
                result.push(value.clone());
            }
        }
    }

    result
}

/// Orders vectors and merges them.
fn order_and_merge<T: Clone, U: Ord>(vectors: Vec<Vec<T>>, numbers: &[U]) -> Vec<T> {
    // sort the vectors
    let mut paired: Vec<(&U, Vec<T>)> = numbers.iter().zip(vectors).collect();
    paired.sort_by(|a, b| a.0.cmp(&b.0));

    // extract the vectors
    let sorted_vectors: Vec<Vec<T>> = paired.into_iter().map(|(_, vec)| vec).collect();

    // merge the sorted vectors
    merge_vectors(sorted_vectors)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn split_and_join() {
        let vec = vec![14, 16, 23, 39, 12, 21, 24, 214, 5, 29, 9, -43, 14, 234, 249];

        for n in 1..20 {
            let split = split_vector_by_n(&vec, n);
            let joined = merge_vectors(split);
            assert_eq!(joined.len(), vec.len());
            for (a, b) in vec.iter().zip(joined.iter()) {
                assert_eq!(a, b);
            }
        }
    }

    #[test]
    fn split_order_and_join() {
        let vec = vec![1, 3, 4, 6, 7, 9, 10, 12];

        let mut split_vec = split_vector_by_n(&vec, 2);
        split_vec.push(vec![2, 5, 8, 11]);

        let ordered = order_and_merge(split_vec, &[0, 2, 1]);
        assert_eq!(ordered, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    }

    #[test]
    fn timewise_merge_simple() {
        let data1 = TimeWiseData {
            thread_ids: vec![0],
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(2.2),
                OrderValue::from(3.3),
                OrderValue::from(4.4),
            ],
            n_samples: vec![1, 2, 3, 4],
        };

        let data2 = TimeWiseData {
            thread_ids: vec![2],
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data_sum = data2 + data1;

        let expected_order: Vec<OrderValue> = [1.1, 21.1, 2.2, 22.2, 3.3, 23.3, 4.4, 24.4]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let expected_n_samples = [1, 21, 2, 22, 3, 23, 4, 24];

        assert_eq!(data_sum.thread_ids, vec![0, 2]);
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
    fn timewise_merge_three_simple() {
        let data1 = TimeWiseData {
            thread_ids: vec![0],
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(2.2),
                OrderValue::from(3.3),
                OrderValue::from(4.4),
            ],
            n_samples: vec![1, 2, 3, 4],
        };

        let data2 = TimeWiseData {
            thread_ids: vec![2],
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data3 = TimeWiseData {
            thread_ids: vec![3],
            order: vec![
                OrderValue::from(31.1),
                OrderValue::from(32.2),
                OrderValue::from(33.3),
                OrderValue::from(34.4),
            ],
            n_samples: vec![31, 32, 33, 34],
        };

        let data_sum = data1 + data2 + data3;

        let expected_order: Vec<OrderValue> = [
            1.1, 21.1, 31.1, 2.2, 22.2, 32.2, 3.3, 23.3, 33.3, 4.4, 24.4, 34.4,
        ]
        .into_iter()
        .map(|x| x.into())
        .collect();
        let expected_n_samples = [1, 21, 31, 2, 22, 32, 3, 23, 33, 4, 24, 34];

        assert_eq!(data_sum.thread_ids, vec![0, 2, 3]);
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
    fn timewise_merge_three_complex() {
        let data1 = TimeWiseData {
            thread_ids: vec![0],
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(2.2),
                OrderValue::from(3.3),
                OrderValue::from(4.4),
            ],
            n_samples: vec![1, 2, 3, 4],
        };

        let data2 = TimeWiseData {
            thread_ids: vec![2],
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data3 = TimeWiseData {
            thread_ids: vec![1],
            order: vec![
                OrderValue::from(11.1),
                OrderValue::from(12.2),
                OrderValue::from(13.3),
                OrderValue::from(14.4),
            ],
            n_samples: vec![11, 12, 13, 14],
        };

        let data_sum = data1 + data2 + data3;

        let expected_order: Vec<OrderValue> = [
            1.1, 11.1, 21.1, 2.2, 12.2, 22.2, 3.3, 13.3, 23.3, 4.4, 14.4, 24.4,
        ]
        .into_iter()
        .map(|x| x.into())
        .collect();
        let expected_n_samples = [1, 11, 21, 2, 12, 22, 3, 13, 23, 4, 14, 24];

        assert_eq!(data_sum.thread_ids, vec![0, 1, 2]);
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
    fn timewise_merge_four_complex() {
        let data1 = TimeWiseData {
            thread_ids: vec![0],
            order: vec![
                OrderValue::from(1.1),
                OrderValue::from(2.2),
                OrderValue::from(3.3),
                OrderValue::from(4.4),
            ],
            n_samples: vec![1, 2, 3, 4],
        };

        let data2 = TimeWiseData {
            thread_ids: vec![2],
            order: vec![
                OrderValue::from(21.1),
                OrderValue::from(22.2),
                OrderValue::from(23.3),
                OrderValue::from(24.4),
            ],
            n_samples: vec![21, 22, 23, 24],
        };

        let data3 = TimeWiseData {
            thread_ids: vec![1],
            order: vec![
                OrderValue::from(11.1),
                OrderValue::from(12.2),
                OrderValue::from(13.3),
                OrderValue::from(14.4),
            ],
            n_samples: vec![11, 12, 13, 14],
        };

        let data4 = TimeWiseData {
            thread_ids: vec![3],
            order: vec![
                OrderValue::from(31.1),
                OrderValue::from(32.2),
                OrderValue::from(33.3),
                OrderValue::from(34.4),
            ],
            n_samples: vec![31, 32, 33, 34],
        };

        let data14 = data1 + data4;
        let data23 = data2 + data3;
        let data_sum = data14 + data23;

        let expected_order: Vec<OrderValue> = [
            1.1, 11.1, 21.1, 31.1, 2.2, 12.2, 22.2, 32.2, 3.3, 13.3, 23.3, 33.3, 4.4, 14.4, 24.4,
            34.4,
        ]
        .into_iter()
        .map(|x| x.into())
        .collect();
        let expected_n_samples = [1, 11, 21, 31, 2, 12, 22, 32, 3, 13, 23, 33, 4, 14, 24, 34];

        assert_eq!(data_sum.thread_ids, vec![0, 1, 2, 3]);
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
            thread_ids: vec![0, 1, 2],
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
