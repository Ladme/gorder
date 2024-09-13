// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains structures and methods for the construction of maps of order parameters.

use std::ops::Add;

use groan_rs::{prelude::GridMap, structures::gridmap::DataOrder};

use crate::{OrderMap, PANIC_MESSAGE};

#[derive(Debug, Clone)]
pub(crate) struct Map {
    params: OrderMap,
    values: GridMap<f32, f32, fn(&f32) -> f32>,
    samples: GridMap<usize, usize, fn(&usize) -> usize>,
}

impl Map {
    pub(super) fn new(params: OrderMap) -> Map {
        let binx = params.bin_size_x();
        let biny = params.bin_size_y();

        let xmin = 0.0 + (binx / 2.0);
        let xmax = 1.0 - (binx / 2.0);
        let ymin = 0.0 + (biny / 2.0);
        let ymax = 1.0 - (biny / 2.0);

        Map {
            params,
            values: GridMap::new(
                (xmin, xmax),
                (ymin, ymax),
                (binx, biny),
                f32::clone as fn(&f32) -> f32,
            )
            .unwrap_or_else(|e|
                 panic!("FATAL GORDER ERROR | Map::new | Could not create gridmap for values. Error: {}. {}", e, PANIC_MESSAGE)),
            samples: GridMap::new(
                (xmin, xmax),
                (ymin, ymax),
                (binx, biny),
                usize::clone as fn(&usize) -> usize,
            ).unwrap_or_else(|e|
                 panic!("FATAL GORDER ERROR | Map::new | Could not create gridmap for samples. Error: {}. {}", e, PANIC_MESSAGE))
        }
    }
}

impl Add<Map> for Map {
    type Output = Map;

    fn add(self, rhs: Map) -> Self::Output {
        let joined_values_map =
            merge_grid_maps(self.values, rhs.values, f32::clone as fn(&f32) -> f32);

        let joined_samples_map = merge_grid_maps(
            self.samples,
            rhs.samples,
            usize::clone as fn(&usize) -> usize,
        );

        Map {
            params: self.params,
            values: joined_values_map,
            samples: joined_samples_map,
        }
    }
}

/// Helper function for merging optional Maps.
pub(super) fn merge_option_maps(lhs: Option<Map>, rhs: Option<Map>) -> Option<Map> {
    match (lhs, rhs) {
        (Some(x), Some(y)) => Some(x + y),
        (None, None) => None,
        (Some(_), None) | (None, Some(_)) => panic!(
            "FATAL GORDER ERROR | merge_option_maps | Inconsistent option value. {}",
            PANIC_MESSAGE
        ),
    }
}

/// Helper function to merge two GridMaps and construct a new GridMap.
fn merge_grid_maps<T, U>(
    lhs: GridMap<T, U, fn(&T) -> U>,
    rhs: GridMap<T, U, fn(&T) -> U>,
    clone_fn: fn(&T) -> U,
) -> GridMap<T, U, fn(&T) -> U>
where
    T: Add<Output = T> + Copy + Default + std::fmt::Debug,
    U: std::fmt::Display,
{
    let merged = lhs
        .extract_raw()
        .zip(rhs.extract_raw())
        .map(|((_, _, &x), (_, _, &y))| x + y)
        .collect::<Vec<T>>();

    GridMap::from_vec(
        lhs.span_x(),
        lhs.span_y(),
        lhs.tile_dim(),
        merged,
        DataOrder::default(),
        clone_fn,
    )
    .unwrap_or_else(|_| {
        panic!(
            "FATAL GORDER ERROR | merge_grid_maps | Could not construct merged map. {}",
            PANIC_MESSAGE
        )
    })
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::analysis::molecule::Order;

    use super::*;

    #[test]
    fn merge_order() {
        let order1 = Order::new(0.78, 31);

        let order2 = Order::new(0.43, 14);

        let merged = order1 + order2;

        assert_relative_eq!(merged.order(), 1.21);
        assert_eq!(merged.n_samples(), 45);
    }

    #[test]
    fn merge_map() {
        let values = vec![1.0, 2.5, 3.0, 4.2, 5.3, 6.1, 7.3, 8.9];
        let map1_values = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::RowMajor,
            f32::to_owned as fn(&f32) -> f32,
        )
        .unwrap();

        let samples = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let map1_samples = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            samples,
            DataOrder::RowMajor,
            usize::to_owned as fn(&usize) -> usize,
        )
        .unwrap();

        let map1 = Map {
            params: OrderMap::new().output_directory(".").build().unwrap(),
            values: map1_values,
            samples: map1_samples,
        };

        let values = vec![0.7, 1.4, 2.1, 1.4, 2.3, 3.1, 3.3, 1.9];
        let map2_values = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::RowMajor,
            f32::to_owned as fn(&f32) -> f32,
        )
        .unwrap();

        let samples = vec![2, 1, 4, 3, 2, 1, 0, 2];
        let map2_samples = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            samples,
            DataOrder::RowMajor,
            usize::to_owned as fn(&usize) -> usize,
        )
        .unwrap();

        let map2 = Map {
            params: OrderMap::new().output_directory(".").build().unwrap(),
            values: map2_values,
            samples: map2_samples,
        };

        let map = map1 + map2;

        let expected_values = vec![1.7, 3.9, 5.1, 5.6, 7.6, 9.2, 10.6, 10.8];
        let expected_samples = vec![3, 3, 7, 7, 7, 7, 7, 10];

        let values = map
            .values
            .extract_raw()
            .map(|(_, _, &x)| x)
            .collect::<Vec<f32>>();
        let samples = map
            .samples
            .extract_raw()
            .map(|(_, _, &x)| x)
            .collect::<Vec<usize>>();

        assert_eq!(values.len(), expected_values.len());
        assert_eq!(samples.len(), expected_samples.len());

        for (got, exp) in values.iter().zip(expected_values.iter()) {
            assert_relative_eq!(got, exp);
        }

        for (got, exp) in samples.iter().zip(expected_samples.iter()) {
            assert_eq!(got, exp);
        }
    }
}
