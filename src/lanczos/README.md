Eigen decomposition of Hermitian matrices using Lanczos algorithm

**This is a fork of the `lanczos` crate. Code in this directory is licensed under the MIT or Apache-2.0 License at your option. The License holder is Tomer Diament (tomer <tomer.dmnt@gmail.com>).**

## Overview

Using [Lanczos algorithm](https://en.wikipedia.org/wiki/Lanczos_algorithm) to estimate the extremal
Eigen values and Eigen vectors of an Symmetrics Hermitian matrix.

Supports both dense and sparse matrices via nalgebra_sparse.

Works well for large sparse matrices

```bash
cargo add lanczos
```

## Examples

```rust
use nalgebra::DMatrix;
use lanczos::{Hermitian, Order};

// ...
let eigen = matrix.eigsh(50, Order::Smallest);

// Sorted by eigenvalue in ascending order
eprintln!("{}", eigen.eigenvalues);
// Columns sorted according to eigenvalues
eprintln!("{}", eigen.eigenvectors);

// Second smallest eigen value
eprintln!("{}", eigen.eigenvalues[1]);
// Eigen vector corresponding to the second smallest eigen value
eprintln!("{}", eigen.eigenvectors.column(1));
```

