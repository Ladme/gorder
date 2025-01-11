// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Integration tests for the calculation of coarse-grained order parameters.

use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
};

use approx::assert_relative_eq;
use gorder::prelude::*;
use std::io::Write;
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
        .skip(skip) // Skip the first line
        .map(|line| line.unwrap())
        .collect()
}

#[test]
fn test_cg_order_basic_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_basic.yaml",
        1
    ));
}

#[test]
fn test_cg_order_basic_ndx_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .index("tests/files/cg.ndx")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("Membrane"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_basic.yaml",
        1
    ));
}

#[test]
fn test_cg_order_basic_table() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/cg_order_basic.tab",
        1
    ));
}

#[test]
fn test_cg_order_basic_xvg() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_xvg(pattern)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        let path_expected = format!("tests/files/cg_order_basic_{}.xvg", molecule);

        assert!(diff_files_ignore_first(&path, &path_expected, 1));
    }
}

#[test]
fn test_cg_order_basic_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/cg_order_basic.csv",
        0
    ));
}

#[test]
fn test_cg_order_basic_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64, 128] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_basic.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_leaflets_yaml() {
    for method in [
        LeafletClassification::global("@membrane", "name PO4"),
        LeafletClassification::local("@membrane", "name PO4", 2.5),
        LeafletClassification::individual("name PO4", "name C4A C4B"),
    ] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_leaflets.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_leaflets_yaml_from_gro() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.gro")
        .bonds("tests/files/cg.bnd")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_leaflets.yaml",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_yaml_from_pdb() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.pdb")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_leaflets.yaml",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_yaml_from_pqr() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.pqr")
        .bonds("tests/files/cg.bnd")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_leaflets.yaml",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64, 128] {
        for method in [
            LeafletClassification::global("@membrane", "name PO4"),
            LeafletClassification::local("@membrane", "name PO4", 2.5),
            LeafletClassification::individual("name PO4", "name C4A C4B"),
        ] {
            let output = NamedTempFile::new().unwrap();
            let path_to_output = output.path().to_str().unwrap();

            let analysis = Analysis::builder()
                .structure("tests/files/cg.tpr")
                .trajectory("tests/files/cg.xtc")
                .output(path_to_output)
                .analysis_type(AnalysisType::cgorder("@membrane"))
                .leaflets(method)
                .n_threads(n_threads)
                .silent()
                .overwrite()
                .build()
                .unwrap();

            analysis.run().unwrap().write().unwrap();

            assert!(diff_files_ignore_first(
                path_to_output,
                "tests/files/cg_order_leaflets.yaml",
                1
            ));
        }
    }
}

#[test]
fn test_cg_order_leaflets_yaml_multiple_threads_various_frequencies() {
    for n_threads in [1, 2, 5, 8, 128] {
        for method in [
            LeafletClassification::global("@membrane", "name PO4"),
            LeafletClassification::local("@membrane", "name PO4", 2.5),
            LeafletClassification::individual("name PO4", "name C4A C4B"),
        ] {
            for freq in [
                Frequency::every(4).unwrap(),
                Frequency::every(20).unwrap(),
                Frequency::every(200).unwrap(),
                Frequency::once(),
            ] {
                let output = NamedTempFile::new().unwrap();
                let path_to_output = output.path().to_str().unwrap();

                let analysis = Analysis::builder()
                    .structure("tests/files/cg.tpr")
                    .trajectory("tests/files/cg.xtc")
                    .output(path_to_output)
                    .analysis_type(AnalysisType::cgorder("@membrane"))
                    .leaflets(method.clone().with_frequency(freq))
                    .n_threads(n_threads)
                    .silent()
                    .overwrite()
                    .build()
                    .unwrap();

                analysis.run().unwrap().write().unwrap();

                assert!(diff_files_ignore_first(
                    path_to_output,
                    "tests/files/cg_order_leaflets.yaml",
                    1
                ));
            }
        }
    }
}

#[test]
fn test_cg_order_leaflets_table() {
    for method in [
        LeafletClassification::global("@membrane", "name PO4"),
        LeafletClassification::local("@membrane", "name PO4", 2.5),
        LeafletClassification::individual("name PO4", "name C4A C4B"),
    ] {
        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output_tab(path_to_table)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_table,
            "tests/files/cg_order_leaflets.tab",
            1
        ));
    }
}

#[test]
fn test_cg_order_leaflets_xvg() {
    for method in [
        LeafletClassification::global("@membrane", "name PO4"),
        LeafletClassification::local("@membrane", "name PO4", 2.5),
        LeafletClassification::individual("name PO4", "name C4A C4B"),
    ] {
        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let pattern = format!("{}/order.xvg", path_to_dir);

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output_xvg(pattern)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        for molecule in ["POPC", "POPE", "POPG"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!("tests/files/cg_order_leaflets_{}.xvg", molecule);

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }
    }
}

#[test]
fn test_cg_order_leaflets_csv() {
    for method in [
        LeafletClassification::global("@membrane", "name PO4"),
        LeafletClassification::local("@membrane", "name PO4", 2.5),
        LeafletClassification::individual("name PO4", "name C4A C4B"),
    ] {
        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output_csv(path_to_csv)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/cg_order_leaflets.csv",
            0
        ));
    }
}

#[test]
fn test_cg_order_limit_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .min_samples(5000)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_limit.yaml",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_limit_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .min_samples(2000)
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_leaflets_limit.yaml",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_limit_tab() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .min_samples(2000)
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/cg_order_leaflets_limit.tab",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_limit_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .min_samples(2000)
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/cg_order_leaflets_limit.csv",
        1
    ));
}

#[test]
fn test_cg_order_begin_end_step_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .begin(350_000.0)
        .end(356_000.0)
        .step(5)
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = analysis.run().unwrap();
    assert_eq!(results.n_analyzed_frames(), 13);
    results.write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_begin_end_step.yaml",
        1
    ));
}

#[test]
fn test_cg_order_begin_end_step_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64, 128] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .begin(350_000.0)
            .end(356_000.0)
            .step(5)
            .leaflets(LeafletClassification::global("@membrane", "name PO4"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        let results = analysis.run().unwrap();
        assert_eq!(results.n_analyzed_frames(), 13);
        results.write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_begin_end_step.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_begin_end_step_yaml_leaflets_multiple_threads_various_frequencies() {
    for n_threads in [1, 2, 5, 8, 128] {
        for method in [
            LeafletClassification::global("@membrane", "name PO4"),
            LeafletClassification::local("@membrane", "name PO4", 2.5),
            LeafletClassification::individual("name PO4", "name C4A C4B"),
        ] {
            for freq in [
                Frequency::every(2).unwrap(),
                Frequency::every(10).unwrap(),
                Frequency::once(),
            ] {
                let output = NamedTempFile::new().unwrap();
                let path_to_output = output.path().to_str().unwrap();

                let analysis = Analysis::builder()
                    .structure("tests/files/cg.tpr")
                    .trajectory("tests/files/cg.xtc")
                    .output(path_to_output)
                    .analysis_type(AnalysisType::cgorder("@membrane"))
                    .begin(350_000.0)
                    .end(356_000.0)
                    .step(5)
                    .leaflets(method.clone().with_frequency(freq))
                    .n_threads(n_threads)
                    .silent()
                    .overwrite()
                    .build()
                    .unwrap();

                let results = analysis.run().unwrap();
                assert_eq!(results.n_analyzed_frames(), 13);
                results.write().unwrap();

                assert!(diff_files_ignore_first(
                    path_to_output,
                    "tests/files/cg_order_begin_end_step.yaml",
                    1
                ));
            }
        }
    }
}

#[test]
fn test_cg_order_begin_end_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .begin(350_000.0)
        .end(356_000.0)
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = analysis.run().unwrap();
    assert_eq!(results.n_analyzed_frames(), 61);
    results.write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_begin_end.yaml",
        1
    ));
}

#[test]
fn test_cg_order_begin_end_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64, 128] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .begin(350_000.0)
            .end(356_000.0)
            .leaflets(LeafletClassification::global("@membrane", "name PO4"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        let results = analysis.run().unwrap();
        assert_eq!(results.n_analyzed_frames(), 61);
        results.write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_begin_end.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_no_molecules() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output("THIS_FILE_SHOULD_NOT_BE_CREATED_CG_1")
        .analysis_type(AnalysisType::cgorder("@ion"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_CG_1").exists());
}

#[test]
fn test_cg_order_empty_molecules() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output("THIS_FILE_SHOULD_NOT_BE_CREATED_CG_2")
        .analysis_type(AnalysisType::cgorder("name PO4 C1A"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_CG_2").exists());
}

macro_rules! create_file_for_backup {
    ($path:expr) => {{
        File::create($path)
            .unwrap()
            .write_all("This file will be backed up.".as_bytes())
            .unwrap()
    }};
}

fn read_and_compare_files(dir: &str, exclude_paths: &[&str], expected_content: &str) {
    let mut count = 0;
    for entry in std::fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.is_dir() || exclude_paths.contains(&path.to_str().unwrap()) {
            continue;
        }

        count += 1;

        let mut file_content = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut file_content)
            .unwrap();

        assert_eq!(file_content, expected_content);
    }

    assert_eq!(count, 6);
}

#[test]
fn test_cg_order_basic_all_formats_backup() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let file_paths = [
        format!("{}/order.yaml", path_to_dir),
        format!("{}/order.tab", path_to_dir),
        format!("{}/order.csv", path_to_dir),
        format!("{}/order_POPC.xvg", path_to_dir),
        format!("{}/order_POPE.xvg", path_to_dir),
        format!("{}/order_POPG.xvg", path_to_dir),
    ];

    for path in &file_paths {
        create_file_for_backup!(path);
    }

    let xvg_pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_yaml(&file_paths[0])
        .output_tab(&file_paths[1])
        .output_csv(&file_paths[2])
        .output_xvg(&xvg_pattern)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    let all_files = [
        ("tests/files/cg_order_basic.yaml", &file_paths[0]),
        ("tests/files/cg_order_basic.tab", &file_paths[1]),
        ("tests/files/cg_order_basic.csv", &file_paths[2]),
        ("tests/files/cg_order_basic_POPC.xvg", &file_paths[3]),
        ("tests/files/cg_order_basic_POPE.xvg", &file_paths[4]),
        ("tests/files/cg_order_basic_POPG.xvg", &file_paths[5]),
    ];

    for (expected, result) in &all_files {
        assert!(diff_files_ignore_first(result, expected, 1));
    }

    read_and_compare_files(
        path_to_dir,
        &file_paths.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
        "This file will be backed up.",
    );
}

#[test]
fn test_cg_order_maps_basic() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder(
            "resname POPC and name C1B C2B C3B C4B",
        ))
        .map(
            OrderMap::builder()
                .bin_size([1.0, 1.0])
                .output_directory(path_to_dir)
                .min_samples(10)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    let expected_file_names = [
        "ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
        "ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
        "ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
        "ordermap_average_full.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/POPC/{}", path_to_dir, file);
        let test_file = format!("tests/files/ordermaps_cg/{}", file);
        assert!(diff_files_ignore_first(&real_file, &test_file, 2));
    }

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_small.yaml",
        1
    ));
}

#[test]
fn test_cg_order_maps_leaflets() {
    for method in [
        LeafletClassification::global("@membrane", "name PO4"),
        LeafletClassification::local("@membrane", "name PO4", 2.5),
        LeafletClassification::individual("name PO4", "name C4A C4B"),
    ] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder(
                "resname POPC and name C1B C2B C3B C4B",
            ))
            .leaflets(method)
            .map(
                OrderMap::builder()
                    .bin_size([1.0, 1.0])
                    .output_directory(path_to_dir)
                    .min_samples(10)
                    .build()
                    .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
            "ordermap_POPC-C1B-8--POPC-C2B-9_upper.dat",
            "ordermap_POPC-C2B-9--POPC-C3B-10_lower.dat",
            "ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
            "ordermap_POPC-C3B-10--POPC-C4B-11_upper.dat",
            "ordermap_POPC-C1B-8--POPC-C2B-9_lower.dat",
            "ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
            "ordermap_POPC-C2B-9--POPC-C3B-10_upper.dat",
            "ordermap_POPC-C3B-10--POPC-C4B-11_lower.dat",
            "ordermap_average_full.dat",
            "ordermap_average_upper.dat",
            "ordermap_average_lower.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps_cg/{}", file);
            assert!(diff_files_ignore_first(&real_file, &test_file, 2));
        }

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_leaflets_small.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_maps_basic_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64, 128] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder(
                "resname POPC and name C1B C2B C3B C4B",
            ))
            .n_threads(n_threads)
            .map(
                OrderMap::builder()
                    .bin_size([1.0, 1.0])
                    .output_directory(path_to_dir)
                    .min_samples(10)
                    .build()
                    .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
            "ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
            "ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
            "ordermap_average_full.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps_cg/{}", file);
            assert!(diff_files_ignore_first(&real_file, &test_file, 2));
        }

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_small.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_maps_basic_backup() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let outer_directory = TempDir::new().unwrap();
    let path_to_outer_dir = outer_directory.path().to_str().unwrap();

    let directory = TempDir::new_in(path_to_outer_dir).unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let backup_file = format!("{}/to_backup.txt", path_to_dir);
    create_file_for_backup!(&backup_file);

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder(
            "resname POPC and name C1B C2B C3B C4B",
        ))
        .map(
            OrderMap::builder()
                .bin_size([1.0, 1.0])
                .output_directory(path_to_dir)
                .min_samples(10)
                .build()
                .unwrap(),
        )
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    let expected_file_names = [
        "ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
        "ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
        "ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
        "ordermap_average_full.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/POPC/{}", path_to_dir, file);
        let test_file = format!("tests/files/ordermaps_cg/{}", file);
        assert!(diff_files_ignore_first(&real_file, &test_file, 2));
    }

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_small.yaml",
        1
    ));

    // check backed up directory
    let directories = std::fs::read_dir(path_to_outer_dir)
        .unwrap()
        .map(|x| x.unwrap().path())
        .filter(|x| x.is_dir() && x.to_str().unwrap() != path_to_dir)
        .collect::<Vec<PathBuf>>();

    assert_eq!(directories.len(), 1);

    let file_path = format!("{}/to_backup.txt", directories[0].display());
    let mut file_content = String::new();
    File::open(&file_path)
        .unwrap()
        .read_to_string(&mut file_content)
        .unwrap();

    assert_eq!(file_content, "This file will be backed up.");
}

#[test]
fn test_cg_order_error_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_error.yaml",
        1
    ));
}

#[test]
fn test_cg_order_error_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .n_threads(n_threads)
            .estimate_error(EstimateError::default())
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_error.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_error_leaflets_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_error_leaflets.yaml",
        1
    ));
}

#[test]
fn test_cg_order_error_leaflets_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(LeafletClassification::global("@membrane", "name PO4"))
            .n_threads(n_threads)
            .estimate_error(EstimateError::default())
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_error_leaflets.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_error_tab() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/cg_order_error.tab",
        1
    ));
}

#[test]
fn test_cg_order_error_leaflets_tab() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/cg_order_error_leaflets.tab",
        1
    ));
}

#[test]
fn test_cg_order_error_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/cg_order_error.csv",
        1
    ));
}

#[test]
fn test_cg_order_error_leaflets_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/cg_order_error_leaflets.csv",
        1
    ));
}

#[test]
fn test_cg_order_error_xvg() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_xvg(pattern)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        // same files as when `estimate_error` is not provided - xvg files do not show error
        let path_expected = format!("tests/files/cg_order_basic_{}.xvg", molecule);

        assert!(diff_files_ignore_first(&path, &path_expected, 1));
    }
}

#[test]
fn test_cg_order_error_leaflets_xvg() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_xvg(pattern)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        // same files as when `estimate_error` is not provided - xvg files do not show error
        let path_expected = format!("tests/files/cg_order_leaflets_{}.xvg", molecule);

        assert!(diff_files_ignore_first(&path, &path_expected, 1));
    }
}

#[test]
fn test_cg_order_error_limit() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .min_samples(5000)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/cg_order_error_limit.yaml",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/cg_order_error_limit.tab",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/cg_order_error_limit.csv",
        1
    ));
}

#[test]
fn test_cg_order_error_leaflets_limit() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .min_samples(2000)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/cg_order_error_leaflets_limit.yaml",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/cg_order_error_leaflets_limit.tab",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/cg_order_error_leaflets_limit.csv",
        1
    ));
}

#[test]
fn test_cg_order_convergence() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_convergence.xvg",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_convergence() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_leaflets_convergence.xvg",
        1
    ));
}

#[test]
fn test_cg_order_convergence_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_convergence.xvg",
            1
        ));
    }
}

#[test]
fn test_cg_order_convergence_step() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
        .step(5)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/cg_order_convergence_s5.xvg",
        1
    ));
}

#[test]
fn test_cg_order_leaflets_scrambling_various_methods_and_frequencies() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    for n_threads in [1, 2, 5, 8, 128] {
        for method in [
            LeafletClassification::global("@membrane", "name PO4"),
            LeafletClassification::local("@membrane", "name PO4", 3.0),
            LeafletClassification::individual("name PO4", "name C4A C4B"),
        ] {
            for freq in [
                Frequency::every(1).unwrap(),
                Frequency::every(10).unwrap(),
                Frequency::once(),
            ] {
                let analysis = Analysis::builder()
                    .structure("tests/files/scrambling/cg_scrambling.tpr")
                    .trajectory("tests/files/scrambling/cg_scrambling.xtc")
                    .analysis_type(AnalysisType::cgorder("@membrane"))
                    .output_yaml(path_to_output)
                    .leaflets(method.clone().with_frequency(freq))
                    .n_threads(n_threads)
                    .silent()
                    .overwrite()
                    .build()
                    .unwrap();

                analysis.run().unwrap().write().unwrap();

                let test_file = match (method.clone(), freq) {
                    (LeafletClassification::Global(_), Frequency::Every(n)) if n.get() == 1 => {
                        "order_global.yaml"
                    }
                    (LeafletClassification::Global(_), Frequency::Every(n)) if n.get() == 10 => {
                        "order_global_every_10.yaml"
                    }
                    (LeafletClassification::Local(_), Frequency::Every(n)) if n.get() == 1 => {
                        "order_local.yaml"
                    }
                    (LeafletClassification::Local(_), Frequency::Every(n)) if n.get() == 10 => {
                        "order_local_every_10.yaml"
                    }
                    (LeafletClassification::Individual(_), Frequency::Every(n)) if n.get() == 1 => {
                        "order_individual.yaml"
                    }
                    (LeafletClassification::Individual(_), Frequency::Every(n))
                        if n.get() == 10 =>
                    {
                        "order_individual_every_10.yaml"
                    }
                    (_, Frequency::Once) => "order_once.yaml",
                    _ => unreachable!("Unexpected method-frequency combination."),
                };

                assert!(diff_files_ignore_first(
                    path_to_output,
                    &format!("tests/files/scrambling/{}", test_file),
                    1
                ));
            }
        }
    }
}

#[test]
fn test_cg_order_leaflets_asymmetric_multiple_threads() {
    for n_threads in [1, 2, 5, 8] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let output_tab = NamedTempFile::new().unwrap();
        let path_to_tab = output_tab.path().to_str().unwrap();

        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let xvg_pattern = format!("{}/order.xvg", path_to_dir);

        let analysis = Analysis::builder()
            .structure("tests/files/asymmetric/cg_asym.tpr")
            .trajectory("tests/files/asymmetric/cg_asym.xtc")
            .output_yaml(path_to_yaml)
            .output_tab(path_to_tab)
            .output_csv(path_to_csv)
            .output_xvg(xvg_pattern)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(LeafletClassification::global("@membrane", "name PO4"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_yaml,
            "tests/files/asymmetric/cg_order_asymmetric.yaml",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_tab,
            "tests/files/asymmetric/cg_order_asymmetric.tab",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/asymmetric/cg_order_asymmetric.csv",
            0
        ));

        for molecule in ["POPE", "POPG"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!(
                "tests/files/asymmetric/cg_order_asymmetric_{}.xvg",
                molecule
            );

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }
    }
}

#[test]
fn test_cg_order_leaflets_asymmetric_error_multiple_threads() {
    for n_threads in [1, 2, 5, 8] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let output_tab = NamedTempFile::new().unwrap();
        let path_to_tab = output_tab.path().to_str().unwrap();

        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let xvg_pattern = format!("{}/order.xvg", path_to_dir);

        let analysis = Analysis::builder()
            .structure("tests/files/asymmetric/cg_asym.tpr")
            .trajectory("tests/files/asymmetric/cg_asym.xtc")
            .output_yaml(path_to_yaml)
            .output_tab(path_to_tab)
            .output_csv(path_to_csv)
            .output_xvg(xvg_pattern)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .leaflets(LeafletClassification::global("@membrane", "name PO4"))
            .estimate_error(EstimateError::default())
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_yaml,
            "tests/files/asymmetric/cg_order_asymmetric_errors.yaml",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_tab,
            "tests/files/asymmetric/cg_order_asymmetric_errors.tab",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/asymmetric/cg_order_asymmetric_errors.csv",
            0
        ));

        for molecule in ["POPE", "POPG"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!(
                "tests/files/asymmetric/cg_order_asymmetric_{}.xvg",
                molecule
            );

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }
    }
}

#[test]
fn test_cg_order_leaflets_asymmetric_ordermaps_multiple_threads() {
    for n_threads in [1, 2, 5, 8] {
        let analysis = Analysis::builder()
            .structure("tests/files/asymmetric/cg_asym.tpr")
            .trajectory("tests/files/asymmetric/cg_asym.xtc")
            .analysis_type(AnalysisType::cgorder("@membrane and name C1B C2B C3B C4B"))
            .leaflets(LeafletClassification::global("@membrane", "name PO4"))
            .ordermaps(OrderMap::builder().bin_size([1.0, 1.0]).build().unwrap())
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        let result = match analysis.run().unwrap() {
            AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
            AnalysisResults::CG(x) => x,
        };

        for molecule in result.molecules() {
            for bond in molecule.bonds() {
                let total = bond.ordermaps().total().as_ref().unwrap();
                let upper = bond.ordermaps().upper().as_ref().unwrap();
                let lower = bond.ordermaps().lower().as_ref().unwrap();

                if molecule.molecule() == "POPE" {
                    for ((_, _, t), (_, _, u)) in
                        total.extract_convert().zip(upper.extract_convert())
                    {
                        assert_relative_eq!(t, u, epsilon = 1e-4);
                    }

                    for (_, _, l) in lower.extract_convert() {
                        assert!(l.is_nan());
                    }
                } else if molecule.molecule() == "POPG" {
                    for ((_, _, t), (_, _, l)) in
                        total.extract_convert().zip(lower.extract_convert())
                    {
                        assert_relative_eq!(t, l, epsilon = 1e-4);
                    }

                    for (_, _, u) in upper.extract_convert() {
                        assert!(u.is_nan());
                    }
                } else {
                    panic!("Unexpected molecule.")
                }
            }
        }
    }
}

#[test]
fn test_cg_order_geometry_cuboid_full() {
    for reference in [
        GeomReference::default(),
        Vector3D::new(2.0, 2.0, 3.0).into(),
        "@membrane".into(),
        GeomReference::center(),
    ] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .geometry(
                Geometry::cuboid(
                    reference,
                    [f32::NEG_INFINITY, f32::INFINITY],
                    [f32::NEG_INFINITY, f32::INFINITY],
                    [f32::NEG_INFINITY, f32::INFINITY],
                )
                .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_basic.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_geometry_cylinder_full() {
    for reference in [
        GeomReference::default(),
        Vector3D::new(2.0, 2.0, 3.0).into(),
        "@membrane".into(),
        GeomReference::center(),
    ] {
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let output = NamedTempFile::new().unwrap();
            let path_to_output = output.path().to_str().unwrap();

            let analysis = Analysis::builder()
                .structure("tests/files/cg.tpr")
                .trajectory("tests/files/cg.xtc")
                .output(path_to_output)
                .analysis_type(AnalysisType::cgorder("@membrane"))
                .geometry(
                    Geometry::cylinder(
                        reference.clone(),
                        f32::INFINITY,
                        [f32::NEG_INFINITY, f32::INFINITY],
                        axis,
                    )
                    .unwrap(),
                )
                .silent()
                .overwrite()
                .build()
                .unwrap();

            analysis.run().unwrap().write().unwrap();

            assert!(diff_files_ignore_first(
                path_to_output,
                "tests/files/cg_order_basic.yaml",
                1
            ));
        }
    }
}

#[test]
fn test_cg_order_geometry_sphere_full() {
    for reference in [
        GeomReference::default(),
        Vector3D::new(2.0, 2.0, 3.0).into(),
        "@membrane".into(),
        GeomReference::center(),
    ] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .geometry(Geometry::sphere(reference, f32::INFINITY).unwrap())
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_basic.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_geometry_cuboid_box_center_square_multiple_threads() {
    for n_threads in [1, 2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .geometry(
                Geometry::cuboid(
                    GeomReference::center(),
                    [-8.0, -2.0],
                    [2.0, 8.0],
                    [f32::NEG_INFINITY, f32::INFINITY],
                )
                .unwrap(),
            )
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_cuboid_square.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_geometry_cylinder_static_multiple_threads() {
    for n_threads in [1, 2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .geometry(
                Geometry::cylinder(
                    Vector3D::new(2.0, 1.0, 0.0),
                    3.25,
                    [f32::NEG_INFINITY, f32::INFINITY],
                    Axis::Z,
                )
                .unwrap(),
            )
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_cylinder.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_geometry_sphere_dynamic_multiple_threads() {
    for n_threads in [1, 2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/cg.tpr")
            .trajectory("tests/files/cg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::cgorder("@membrane"))
            .geometry(Geometry::sphere("resid 1", 2.5).unwrap())
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/cg_order_sphere.yaml",
            1
        ));
    }
}

#[test]
fn test_cg_order_geometry_cuboid_z() {
    let analysis_geometry = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder(
            "@membrane and name PO4 NC3 NH3 GL0 GL1 GL2 C1A C1B",
        ))
        .geometry(
            Geometry::cuboid(
                "@membrane",
                [f32::NEG_INFINITY, f32::INFINITY],
                [f32::NEG_INFINITY, f32::INFINITY],
                [0.0, 3.5],
            )
            .unwrap(),
        )
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let analysis_leaflets = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder(
            "@membrane and name PO4 NC3 NH3 GL0 GL1 GL2 C1A C1B",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results_geometry = analysis_geometry.run().unwrap();
    let results_leaflets = analysis_leaflets.run().unwrap();

    let results_geometry = match results_geometry {
        AnalysisResults::CG(x) => x,
        AnalysisResults::AA(_) => panic!("Invalid results."),
    };

    let results_leaflets = match results_leaflets {
        AnalysisResults::CG(x) => x,
        AnalysisResults::AA(_) => panic!("Invalid results."),
    };

    for (mol, mol2) in results_geometry
        .molecules()
        .zip(results_leaflets.molecules())
    {
        for (bond, bond2) in mol.bonds().zip(mol2.bonds()) {
            assert_relative_eq!(
                bond.order().total().unwrap().value(),
                bond.order().upper().unwrap().value()
            );
            assert!(bond.order().lower().unwrap().value().is_nan());

            assert_relative_eq!(
                bond.order().total().unwrap().value(),
                bond2.order().upper().unwrap().value()
            );
        }
    }
}

#[test]
fn test_cg_order_geometry_cylinder_z() {
    let analysis_geometry = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder(
            "@membrane and name PO4 NC3 NH3 GL0 GL1 GL2 C1A C1B",
        ))
        .geometry(Geometry::cylinder("@membrane", f32::INFINITY, [0.0, 3.5], Axis::Z).unwrap())
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let analysis_leaflets = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder(
            "@membrane and name PO4 NC3 NH3 GL0 GL1 GL2 C1A C1B",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results_geometry = analysis_geometry.run().unwrap();
    let results_leaflets = analysis_leaflets.run().unwrap();

    let results_geometry = match results_geometry {
        AnalysisResults::CG(x) => x,
        AnalysisResults::AA(_) => panic!("Invalid results."),
    };

    let results_leaflets = match results_leaflets {
        AnalysisResults::CG(x) => x,
        AnalysisResults::AA(_) => panic!("Invalid results."),
    };

    for (mol, mol2) in results_geometry
        .molecules()
        .zip(results_leaflets.molecules())
    {
        for (bond, bond2) in mol.bonds().zip(mol2.bonds()) {
            assert_relative_eq!(
                bond.order().total().unwrap().value(),
                bond.order().upper().unwrap().value()
            );
            assert!(bond.order().lower().unwrap().value().is_nan());

            assert_relative_eq!(
                bond.order().total().unwrap().value(),
                bond2.order().upper().unwrap().value()
            );
        }
    }
}

#[test]
fn test_cg_order_basic_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
        AnalysisResults::CG(x) => x,
    };

    assert_eq!(results.n_analyzed_frames(), 101);
    assert_eq!(results.analysis().structure(), "tests/files/cg.tpr");

    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    assert_relative_eq!(
        results.average_order().total().unwrap().value(),
        0.2962,
        epsilon = 1e-4
    );
    assert!(results.average_order().upper().is_none());
    assert!(results.average_order().lower().is_none());

    assert!(results.average_ordermaps().total().is_none());
    assert!(results.average_ordermaps().upper().is_none());
    assert!(results.average_ordermaps().lower().is_none());

    let expected_molecule_names = ["POPC", "POPE", "POPG"];
    let expected_average_orders = [0.2943, 0.2972, 0.3059];
    let expected_bond_orders = [0.3682, 0.3759, 0.3789];

    for (i, molecule) in results.molecules().enumerate() {
        assert_eq!(molecule.molecule(), expected_molecule_names[i]);

        let average_order = molecule.average_order();
        assert_relative_eq!(
            average_order.total().unwrap().value(),
            expected_average_orders[i],
            epsilon = 1e-4
        );
        assert!(average_order.total().unwrap().error().is_none());
        assert!(average_order.upper().is_none());
        assert!(average_order.lower().is_none());

        let average_maps = molecule.average_ordermaps();
        assert!(average_maps.total().is_none());
        assert!(average_maps.upper().is_none());
        assert!(average_maps.lower().is_none());

        // bonds
        assert_eq!(molecule.bonds().count(), 11);

        let bond = molecule.get_bond(4, 5).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), "C1A");
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), "D2A");
        assert_eq!(a2.relative_index(), 5);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);

        let order = bond.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_bond_orders[i],
            epsilon = 1e-4
        );
        assert!(order.total().unwrap().error().is_none());
        assert!(order.upper().is_none());
        assert!(order.lower().is_none());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // the same bond
        let bond = molecule.get_bond(5, 4).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a2.relative_index(), 5);

        // nonexistent bond
        assert!(molecule.get_bond(1, 3).is_none());
        assert!(molecule.get_bond(15, 16).is_none());
    }
}

#[test]
fn test_cg_order_error_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
        AnalysisResults::CG(x) => x,
    };

    assert_eq!(results.n_analyzed_frames(), 101);
    assert!(results
        .analysis()
        .estimate_error()
        .as_ref()
        .unwrap()
        .output_convergence()
        .is_none());

    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    assert_relative_eq!(
        results.average_order().total().unwrap().value(),
        0.2962,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().total().unwrap().error().unwrap(),
        0.0050,
        epsilon = 1e-4
    );
    assert!(results.average_order().upper().is_none());
    assert!(results.average_order().lower().is_none());

    assert!(results.average_ordermaps().total().is_none());
    assert!(results.average_ordermaps().upper().is_none());
    assert!(results.average_ordermaps().lower().is_none());

    let expected_molecule_names = ["POPC", "POPE", "POPG"];
    let expected_average_orders = [0.2943, 0.2972, 0.3059];
    let expected_average_errors = [0.0067, 0.0052, 0.0089];

    let expected_bond_orders = [0.3682, 0.3759, 0.3789];
    let expected_bond_errors = [0.0125, 0.0164, 0.0159];

    let expected_convergence_frames = (1..=101).collect::<Vec<usize>>();
    let expected_convergence_values = [
        [0.2756, 0.2902, 0.2943],
        [0.2830, 0.2995, 0.2972],
        [0.3198, 0.3066, 0.3059],
    ];

    for (i, molecule) in results.molecules().enumerate() {
        assert_eq!(molecule.molecule(), expected_molecule_names[i]);

        let average_order = molecule.average_order();
        assert_relative_eq!(
            average_order.total().unwrap().value(),
            expected_average_orders[i],
            epsilon = 1e-4
        );
        assert_relative_eq!(
            average_order.total().unwrap().error().unwrap(),
            expected_average_errors[i],
            epsilon = 1e-4
        );
        assert!(average_order.upper().is_none());
        assert!(average_order.lower().is_none());

        let average_maps = molecule.average_ordermaps();
        assert!(average_maps.total().is_none());
        assert!(average_maps.upper().is_none());
        assert!(average_maps.lower().is_none());

        // convergence (available even if `output_convergence` is not specified)
        let convergence = molecule.convergence().unwrap();
        assert_eq!(
            convergence.frames().len(),
            expected_convergence_frames.len()
        );
        for (val, exp) in convergence
            .frames()
            .iter()
            .zip(expected_convergence_frames.iter())
        {
            assert_eq!(val, exp);
        }

        for (j, frame) in [0, 50, 100].into_iter().enumerate() {
            assert_relative_eq!(
                *convergence.total().as_ref().unwrap().get(frame).unwrap(),
                expected_convergence_values.get(i).unwrap().get(j).unwrap(),
                epsilon = 1e-4
            );
        }

        assert!(convergence.upper().is_none());
        assert!(convergence.lower().is_none());

        // bonds
        assert_eq!(molecule.bonds().count(), 11);

        let bond = molecule.get_bond(4, 5).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), "C1A");
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), "D2A");
        assert_eq!(a2.relative_index(), 5);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);

        let order = bond.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_bond_orders[i],
            epsilon = 1e-4
        );
        assert_relative_eq!(
            order.total().unwrap().error().unwrap(),
            expected_bond_errors[i],
            epsilon = 1e-4
        );
        assert!(order.upper().is_none());
        assert!(order.lower().is_none());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // the same bond
        let bond = molecule.get_bond(5, 4).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a2.relative_index(), 5);

        // nonexistent bond
        assert!(molecule.get_bond(1, 3).is_none());
        assert!(molecule.get_bond(15, 16).is_none());
    }
}

#[test]
fn test_cg_order_leaflets_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
        AnalysisResults::CG(x) => x,
    };

    assert_eq!(results.n_analyzed_frames(), 101);
    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    assert_relative_eq!(
        results.average_order().total().unwrap().value(),
        0.2962,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().upper().unwrap().value(),
        0.2971,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().lower().unwrap().value(),
        0.2954,
        epsilon = 1e-4
    );

    assert!(results.average_ordermaps().total().is_none());
    assert!(results.average_ordermaps().upper().is_none());
    assert!(results.average_ordermaps().lower().is_none());

    let expected_molecule_names = ["POPC", "POPE", "POPG"];
    let expected_average_orders = [0.2943, 0.2972, 0.3059];
    let expected_average_upper = [0.2965, 0.2965, 0.3085];
    let expected_average_lower = [0.2920, 0.2980, 0.3033];

    let expected_bond_orders = [0.3682, 0.3759, 0.3789];
    let expected_bond_upper = [0.3647, 0.3713, 0.4129];
    let expected_bond_lower = [0.3717, 0.3806, 0.3449];

    for (i, molecule) in results.molecules().enumerate() {
        assert_eq!(molecule.molecule(), expected_molecule_names[i]);

        let average_order = molecule.average_order();
        assert_relative_eq!(
            average_order.total().unwrap().value(),
            expected_average_orders[i],
            epsilon = 1e-4
        );
        assert!(average_order.total().unwrap().error().is_none());

        assert_relative_eq!(
            average_order.upper().unwrap().value(),
            expected_average_upper[i],
            epsilon = 1e-4
        );
        assert!(average_order.upper().unwrap().error().is_none());

        assert_relative_eq!(
            average_order.lower().unwrap().value(),
            expected_average_lower[i],
            epsilon = 1e-4
        );
        assert!(average_order.lower().unwrap().error().is_none());

        let average_maps = molecule.average_ordermaps();
        assert!(average_maps.total().is_none());
        assert!(average_maps.upper().is_none());
        assert!(average_maps.lower().is_none());

        // bonds
        assert_eq!(molecule.bonds().count(), 11);

        let bond = molecule.get_bond(4, 5).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), "C1A");
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), "D2A");
        assert_eq!(a2.relative_index(), 5);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);

        let order = bond.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_bond_orders[i],
            epsilon = 1e-4
        );
        assert!(order.total().unwrap().error().is_none());

        assert_relative_eq!(
            order.upper().unwrap().value(),
            expected_bond_upper[i],
            epsilon = 1e-4
        );
        assert!(order.upper().unwrap().error().is_none());

        assert_relative_eq!(
            order.lower().unwrap().value(),
            expected_bond_lower[i],
            epsilon = 1e-4
        );
        assert!(order.lower().unwrap().error().is_none());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // the same bond
        let bond = molecule.get_bond(5, 4).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a2.relative_index(), 5);

        // nonexistent bond
        assert!(molecule.get_bond(1, 3).is_none());
        assert!(molecule.get_bond(15, 16).is_none());
    }
}

#[test]
fn test_cg_order_error_leaflets_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder("@membrane"))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
        AnalysisResults::CG(x) => x,
    };

    assert_eq!(results.n_analyzed_frames(), 101);
    assert!(results
        .analysis()
        .estimate_error()
        .as_ref()
        .unwrap()
        .output_convergence()
        .is_none());

    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    assert_relative_eq!(
        results.average_order().total().unwrap().value(),
        0.2962,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().total().unwrap().error().unwrap(),
        0.0050,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().upper().unwrap().value(),
        0.2971,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().upper().unwrap().error().unwrap(),
        0.0049,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().lower().unwrap().value(),
        0.2954,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        results.average_order().lower().unwrap().error().unwrap(),
        0.0056,
        epsilon = 1e-4
    );

    assert!(results.average_ordermaps().total().is_none());
    assert!(results.average_ordermaps().upper().is_none());
    assert!(results.average_ordermaps().lower().is_none());

    let expected_molecule_names = ["POPC", "POPE", "POPG"];

    for (i, molecule) in results.molecules().enumerate() {
        assert_eq!(molecule.molecule(), expected_molecule_names[i]);

        let average_order = molecule.average_order();
        assert!(average_order.total().unwrap().error().is_some());
        assert!(average_order.upper().unwrap().error().is_some());
        assert!(average_order.lower().unwrap().error().is_some());

        let average_maps = molecule.average_ordermaps();
        assert!(average_maps.total().is_none());
        assert!(average_maps.upper().is_none());
        assert!(average_maps.lower().is_none());

        // convergence (available even if `output_convergence` is not specified)
        let convergence = molecule.convergence().unwrap();
        assert_eq!(convergence.frames().len(), 101);
        assert!(convergence.total().is_some());
        assert!(convergence.upper().is_some());
        assert!(convergence.lower().is_some());

        // bonds
        assert_eq!(molecule.bonds().count(), 11);

        let bond = molecule.get_bond(4, 5).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), "C1A");
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), "D2A");
        assert_eq!(a2.relative_index(), 5);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);

        let order = bond.order();
        assert!(order.total().unwrap().error().is_some());
        assert!(order.upper().unwrap().error().is_some());
        assert!(order.lower().unwrap().error().is_some());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // the same bond
        let bond = molecule.get_bond(5, 4).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), 4);
        assert_eq!(a2.relative_index(), 5);

        // nonexistent bond
        assert!(molecule.get_bond(1, 3).is_none());
        assert!(molecule.get_bond(15, 16).is_none());
    }
}

#[test]
fn test_cg_order_ordermaps_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder(
            "resname POPC and name C1B C2B C3B C4B",
        ))
        .map(
            OrderMap::builder()
                .bin_size([1.0, 1.0])
                .min_samples(10)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
        AnalysisResults::CG(x) => x,
    };

    assert_eq!(results.n_analyzed_frames(), 101);
    assert_eq!(results.molecules().count(), 1);

    // average ordermaps for the entire molecule
    let molecule = results.get_molecule("POPC").unwrap();
    let map = molecule.average_ordermaps().total().as_ref().unwrap();
    assert!(molecule.average_ordermaps().upper().is_none());
    assert!(molecule.average_ordermaps().lower().is_none());

    let span_x = map.span_x();
    let span_y = map.span_y();
    let bin = map.tile_dim();

    assert_relative_eq!(span_x.0, 0.0);
    assert_relative_eq!(span_x.1, 12.747616);
    assert_relative_eq!(span_y.0, 0.0);
    assert_relative_eq!(span_y.1, 12.747616);
    assert_relative_eq!(bin.0, 1.0);
    assert_relative_eq!(bin.1, 1.0);

    assert_relative_eq!(
        map.get_at_convert(1.0, 8.0).unwrap(),
        0.3590,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(7.0, 0.0).unwrap(),
        0.3765,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(13.0, 11.0).unwrap(),
        0.4296,
        epsilon = 1e-4
    );

    // ordermaps for a selected bond
    let bond = molecule.get_bond(9, 10).unwrap();
    let map = bond.ordermaps().total().as_ref().unwrap();
    assert!(bond.ordermaps().upper().is_none());
    assert!(bond.ordermaps().lower().is_none());

    assert_relative_eq!(
        map.get_at_convert(1.0, 8.0).unwrap(),
        0.3967,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(7.0, 0.0).unwrap(),
        0.3213,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(13.0, 11.0).unwrap(),
        0.4104,
        epsilon = 1e-4
    );
}

#[test]
fn test_cg_order_ordermaps_leaflets_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/cg.tpr")
        .trajectory("tests/files/cg.xtc")
        .analysis_type(AnalysisType::cgorder(
            "resname POPC and name C1B C2B C3B C4B",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name PO4"))
        .map(
            OrderMap::builder()
                .bin_size([1.0, 1.0])
                .min_samples(10)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(_) => panic!("Incorrect results type returned."),
        AnalysisResults::CG(x) => x,
    };

    assert_eq!(results.n_analyzed_frames(), 101);
    assert_eq!(results.molecules().count(), 1);

    // average ordermaps for the entire molecule
    let molecule = results.get_molecule("POPC").unwrap();
    let total = molecule.average_ordermaps().total().as_ref().unwrap();

    let span_x = total.span_x();
    let span_y = total.span_y();
    let bin = total.tile_dim();

    assert_relative_eq!(span_x.0, 0.0);
    assert_relative_eq!(span_x.1, 12.747616);
    assert_relative_eq!(span_y.0, 0.0);
    assert_relative_eq!(span_y.1, 12.747616);
    assert_relative_eq!(bin.0, 1.0);
    assert_relative_eq!(bin.1, 1.0);

    assert_relative_eq!(
        total.get_at_convert(1.0, 8.0).unwrap(),
        0.3590,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        total.get_at_convert(13.0, 11.0).unwrap(),
        0.4296,
        epsilon = 1e-4
    );

    let upper = molecule.average_ordermaps().upper().as_ref().unwrap();

    assert_relative_eq!(
        upper.get_at_convert(1.0, 8.0).unwrap(),
        0.3418,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        upper.get_at_convert(13.0, 11.0).unwrap(),
        0.4051,
        epsilon = 1e-4
    );

    let lower = molecule.average_ordermaps().lower().as_ref().unwrap();

    assert_relative_eq!(
        lower.get_at_convert(1.0, 8.0).unwrap(),
        0.3662,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        lower.get_at_convert(13.0, 11.0).unwrap(),
        0.4506,
        epsilon = 1e-4
    );

    // ordermaps for a selected bond
    let bond = molecule.get_bond(9, 10).unwrap();
    let total = bond.ordermaps().total().as_ref().unwrap();

    assert_relative_eq!(
        total.get_at_convert(1.0, 8.0).unwrap(),
        0.3967,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        total.get_at_convert(13.0, 11.0).unwrap(),
        0.4104,
        epsilon = 1e-4
    );

    let upper = bond.ordermaps().upper().as_ref().unwrap();

    assert_relative_eq!(
        upper.get_at_convert(1.0, 8.0).unwrap(),
        0.3573,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        upper.get_at_convert(13.0, 11.0).unwrap(),
        0.4807,
        epsilon = 1e-4
    );

    let lower = bond.ordermaps().lower().as_ref().unwrap();

    assert_relative_eq!(
        lower.get_at_convert(1.0, 8.0).unwrap(),
        0.4118,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        lower.get_at_convert(13.0, 11.0).unwrap(),
        0.3563,
        epsilon = 1e-4
    );
}
