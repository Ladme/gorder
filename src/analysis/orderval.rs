// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of a custom structure for precise calculation of order parameters.

use std::ops::{Add, AddAssign, Div};

use serde::Deserialize;

use crate::PANIC_MESSAGE;

const PRECISION: f64 = 1_000_000.0;

/// Value of the order parameter x PRECISION rounded to the nearest integer.
/// Avoids issues with floating point precision.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Deserialize, Default)]
#[serde(transparent)]
pub(super) struct OrderValue(i64);

impl From<f32> for OrderValue {
    fn from(value: f32) -> Self {
        // we do not check for overflow here as `value` should never be larger than 1.0
        OrderValue((value as f64 * PRECISION).round() as i64)
    }
}

impl From<OrderValue> for f32 {
    fn from(value: OrderValue) -> Self {
        ((value.0 as f64) / PRECISION) as f32
    }
}

impl Div<usize> for OrderValue {
    type Output = OrderValue;

    fn div(self, rhs: usize) -> Self::Output {
        OrderValue(self.0 / TryInto::<i64>::try_into(rhs)
            .unwrap_or_else(|e| panic!("FATAL GORDER ERROR | OrderValue::div | Conversion of usize to i64 failed. {}. Value of '{}' cannot be converted. {}", e, rhs, PANIC_MESSAGE))
        )
    }
}

impl Add<OrderValue> for OrderValue {
    type Output = OrderValue;

    fn add(self, rhs: OrderValue) -> Self::Output {
        OrderValue(
            self.0.checked_add(rhs.0)
                .unwrap_or_else(|| panic!("FATAL GORDER ERROR | OrderValue::Add | OrderValue overflowed. Tried adding '{}' and '{}'. {}", self.0, rhs.0, PANIC_MESSAGE))
        )
    }
}

impl AddAssign<OrderValue> for OrderValue {
    fn add_assign(&mut self, rhs: OrderValue) {
        self.0 = self.0.checked_add(rhs.0)
            .unwrap_or_else(|| panic!("FATAL GORDER ERROR | OrderValue::AddAssign | OrderValue overflowed. Tried adding '{}' and '{}'. {}", self.0, rhs.0, PANIC_MESSAGE));
    }
}

impl AddAssign<f32> for OrderValue {
    fn add_assign(&mut self, rhs: f32) {
        self.add_assign(OrderValue::from(rhs))
    }
}
