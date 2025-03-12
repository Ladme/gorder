// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Integration tests for the calculation of united-atom order parameters.

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use gorder::input::{Analysis, AnalysisType};
use tempfile::{NamedTempFile, TempDir};

/// Test utility. Diff the contents of two files without the first `skip` lines.
pub(crate) fn diff_files_ignore_first(file1: &str, file2: &str, skip: usize) -> bool {
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

#[test]
fn test_ua_order_basic() {
    for n_threads in [1, 2, 3, 4, 8, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_yaml = output.path().to_str().unwrap();

        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let pattern = format!("{}/order.xvg", path_to_dir);

        let analysis = Analysis::builder()
        .structure("tests/files/ua.tpr")
        .trajectory("tests/files/ua.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .output_csv(path_to_csv)
        .output_xvg(&pattern)
        .analysis_type(AnalysisType::uaorder(
            Some("(resname POPC and name r'^C' and not name C15 C34 C24 C25) or (resname POPS and name r'^C' and not name C6 C18 C39 C27 C28)"),
            Some("(resname POPC and name C24 C25) or (resname POPS and name C27 C28)"),
            None
        ))
        .n_threads(n_threads)
        .silent()
        .overwrite()
        .build()
        .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_yaml,
            "tests/files/ua_order_basic.yaml",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_table,
            "tests/files/ua_order_basic.tab",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/ua_order_basic.csv",
            0
        ));

        for molecule in ["POPC", "POPS"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!("tests/files/ua_order_basic_{}.xvg", molecule);

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }
    }
}

#[test]
fn test_ua_order_basic_saturated_only() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/ua.tpr")
        .trajectory("tests/files/ua.xtc")
        .output_yaml(path_to_yaml)
        .analysis_type(AnalysisType::uaorder(
            Some("(resname POPC and name r'^C' and not name C15 C34 C24 C25) or (resname POPS and name r'^C' and not name C6 C18 C39 C27 C28)"),
            None,
            None
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/ua_order_basic_saturated.yaml",
        1
    ));
}

#[test]
fn test_ua_order_basic_unsaturated_only() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/ua.tpr")
        .trajectory("tests/files/ua.xtc")
        .output_yaml(path_to_yaml)
        .analysis_type(AnalysisType::uaorder(
            None,
            Some("(resname POPC and name C24 C25) or (resname POPS and name C27 C28)"),
            None,
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/ua_order_basic_unsaturated.yaml",
        1
    ));
}
