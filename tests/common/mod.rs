// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Functions used in various integration tests.

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use approx::assert_relative_eq;

/// Test utility. Diff the contents of two files without the first `skip` lines.
#[allow(dead_code)]
pub(super) fn diff_files_ignore_first(file1: &str, file2: &str, skip: usize) -> bool {
    let content1 = read_file_without_first_lines(file1, skip);
    let content2 = read_file_without_first_lines(file2, skip);
    content1 == content2
}

fn read_file_without_first_lines(file: &str, skip: usize) -> Vec<String> {
    let reader = BufReader::new(File::open(file).unwrap());
    reader
        .lines()
        .skip(skip) // skip the first line
        .map(|line| line.unwrap())
        .collect()
}

/// Test utility. Assert that two ordermap files match each other.
#[allow(dead_code)]
pub(super) fn assert_eq_maps(a: &str, b: &str, skip: usize) {
    let (file_a, file_b) = match (File::open(a), File::open(b)) {
        (Ok(f1), Ok(f2)) => (f1, f2),
        _ => panic!("One or both files do not exist."),
    };

    let mut lines_a = BufReader::new(file_a).lines().skip(skip);
    let mut lines_b = BufReader::new(file_b).lines().skip(skip);

    loop {
        match (lines_a.next(), lines_b.next()) {
            (Some(Ok(line_a)), Some(Ok(line_b))) => {
                let is_data = line_a
                    .split_whitespace()
                    .next()
                    .and_then(|s| s.parse::<f32>().ok())
                    .is_some();

                if is_data {
                    let p: Vec<_> = line_a.split_whitespace().collect();
                    let q: Vec<_> = line_b.split_whitespace().collect();
                    assert_eq!(p.len(), 3, "Data line must have 3 columns");
                    assert_eq!(q.len(), 3, "Data line must have 3 columns");

                    assert_eq!(p[0], q[0], "First columns differ");
                    assert_eq!(p[1], q[1], "Second columns differ");

                    match (p[2].parse::<f32>(), q[2].parse::<f32>()) {
                        (Ok(z1), Ok(z2)) if z1.is_nan() && z2.is_nan() => (),
                        (Ok(z1), Ok(z2)) => assert_relative_eq!(z1, z2, epsilon = 1e-3),
                        _ => panic!(
                            "Invalid or mismatched numeric formats: {} vs {}",
                            p[2], q[2]
                        ),
                    }
                } else {
                    assert_eq!(line_a, line_b, "Non-data lines differ");
                }
            }
            (None, None) => break,
            _ => panic!("Files have different number of lines"),
        }
    }
}
