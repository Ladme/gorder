// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Integration tests for the `gorder` library.

use std::{
    fs::File,
    io::Read,
    path::{Path, PathBuf},
};

use approx::assert_relative_eq;
use gorder::prelude::*;
use groan_rs::prelude::GridMap;
use std::io::Write;
use tempfile::{NamedTempFile, TempDir};

#[test]
fn test_aa_order_basic_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_basic.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_basic_fail_overlap() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon or serial 876 to 1234",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    match analysis.run() {
        Ok(_) => panic!("Function should have failed."),
        Err(e) => assert!(e.to_string().contains("are part of both")),
    }
}

#[test]
fn test_aa_order_basic_table() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_table).unwrap();
    let mut expected = File::open("tests/files/aa_order_basic.tab").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_basic_xvg() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_xvg(pattern)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        let path_expected = format!("tests/files/aa_order_basic_{}.xvg", molecule);
        let mut result = File::open(path).unwrap();
        let mut expected = File::open(path_expected).unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_basic_csv() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_csv).unwrap();
    let mut expected = File::open("tests/files/aa_order_basic.csv").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_basic_xvg_weird_names() {
    for name in ["order", ".this.is.a.weird.name.xvg"] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let pattern = format!("{}/{}", path_to_dir, name);

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output_yaml(path_to_yaml)
            .output_xvg(pattern)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        for molecule in ["POPC", "POPE", "POPG"] {
            let path = if name.contains(".xvg") {
                format!("{}/.this.is.a.weird.name_{}.xvg", path_to_dir, molecule)
            } else {
                format!("{}/order_{}", path_to_dir, molecule)
            };
            println!("{:?}", path);

            let path_expected = format!("tests/files/aa_order_basic_{}.xvg", molecule);
            let mut result = File::open(path).unwrap();
            let mut expected = File::open(path_expected).unwrap();

            assert!(file_diff::diff_files(&mut result, &mut expected));
        }
    }
}

#[test]
fn test_aa_order_basic_yaml_multiple_threads() {
    for n_threads in [3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_basic.yaml").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_basic_table_multiple_threads() {
    for n_threads in [3, 5, 8, 12, 16, 64] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output_yaml(path_to_yaml)
            .output_tab(path_to_table)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_table).unwrap();
        let mut expected = File::open("tests/files/aa_order_basic.tab").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_leaflets_yaml() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.5),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_leaflets.yaml").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_leaflets_table() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.5),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output_yaml(path_to_yaml)
            .output_tab(path_to_table)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_table).unwrap();
        let mut expected = File::open("tests/files/aa_order_leaflets.tab").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_leaflets_xvg() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.5),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let pattern = format!("{}/order.xvg", path_to_dir);

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output_yaml(path_to_yaml)
            .output_xvg(pattern)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        for molecule in ["POPC", "POPE", "POPG"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!("tests/files/aa_order_leaflets_{}.xvg", molecule);
            let mut result = File::open(path).unwrap();
            let mut expected = File::open(path_expected).unwrap();

            assert!(file_diff::diff_files(&mut result, &mut expected));
        }
    }
}

#[test]
fn test_aa_order_leaflets_csv() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.5),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let output_yaml = NamedTempFile::new().unwrap();
        let path_to_yaml = output_yaml.path().to_str().unwrap();

        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output_yaml(path_to_yaml)
            .output_csv(path_to_csv)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .leaflets(method)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_csv).unwrap();
        let mut expected = File::open("tests/files/aa_order_leaflets.csv").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_leaflets_yaml_supershort() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg_selected.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_selected.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_one_different_hydrogen_numbers_table() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::aaorder(
            "(resname POPC and name C29 C210) or (resname POPE and element name carbon)",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_table).unwrap();
    let mut expected = File::open("tests/files/aa_order_different_hydrogen_numbers.tab").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_one_different_hydrogen_numbers_csv() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "(resname POPC and name C29 C210) or (resname POPE and element name carbon)",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_csv).unwrap();
    let mut expected = File::open("tests/files/aa_order_different_hydrogen_numbers.csv").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_limit_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .min_samples(2000)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_limit.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_leaflets_limit_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .min_samples(500)
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_leaflets_limit.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_leaflets_limit_tab() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .min_samples(500)
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_table).unwrap();
    let mut expected = File::open("tests/files/aa_order_leaflets_limit.tab").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_leaflets_limit_csv() {
    let output_yaml = NamedTempFile::new().unwrap();
    let path_to_yaml = output_yaml.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .min_samples(500)
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_csv).unwrap();
    let mut expected = File::open("tests/files/aa_order_leaflets_limit.csv").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_begin_end_step_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .begin(450_000.0)
        .end(450_200.0)
        .step(3)
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_begin_end_step.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_no_molecules() {
    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output("THIS_FILE_SHOULD_NOT_BE_CREATED_1")
        .analysis_type(AnalysisType::aaorder(
            "@ion",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_1").exists());
}

#[test]
fn test_aa_order_empty_molecules() {
    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output("THIS_FILE_SHOULD_NOT_BE_CREATED_2")
        .analysis_type(AnalysisType::aaorder(
            "@water and element symbol O",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_2").exists());
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
fn test_aa_order_basic_all_formats_backup() {
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

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(&file_paths[0])
        .output_tab(&file_paths[1])
        .output_csv(&file_paths[2])
        .output_xvg(&xvg_pattern)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let all_files = [
        ("tests/files/aa_order_basic.yaml", &file_paths[0]),
        ("tests/files/aa_order_basic.tab", &file_paths[1]),
        ("tests/files/aa_order_basic.csv", &file_paths[2]),
        ("tests/files/aa_order_basic_POPC.xvg", &file_paths[3]),
        ("tests/files/aa_order_basic_POPE.xvg", &file_paths[4]),
        ("tests/files/aa_order_basic_POPG.xvg", &file_paths[5]),
    ];

    for (expected, result) in &all_files {
        let mut result_file = File::open(result).unwrap();
        let mut expected_file = File::open(expected).unwrap();
        assert!(file_diff::diff_files(&mut result_file, &mut expected_file));
    }

    read_and_compare_files(
        &path_to_dir,
        &file_paths.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
        "This file will be backed up.",
    );
}

fn parse_float(string: &str) -> Option<f32> {
    string.parse::<f32>().ok()
}

fn compare_ordermaps(file1: impl AsRef<Path>, file2: impl AsRef<Path>) {
    let map2 =
        GridMap::from_file(file2, f32::clone, &[' '], parse_float, &["#", "$", "@"]).unwrap();

    let map1 =
        GridMap::from_file(file1, f32::clone, &[' '], parse_float, &["#", "$", "@"]).unwrap();

    for ((x1, y1, z1), (x2, y2, z2)) in map1.extract_convert().zip(map2.extract_convert()) {
        assert_relative_eq!(x1, x2);
        assert_relative_eq!(y1, y2);
        if !z1.is_nan() || !z2.is_nan() {
            assert_relative_eq!(z1, z2, epsilon = 0.0002);
        }
    }
}

#[test]
fn test_aa_order_maps_basic() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .map(
            OrderMap::new()
                .bin_size_y(4.0)
                .bin_size_y(4.0)
                .output_directory(path_to_dir)
                .min_samples(5)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let expected_file_names = [
        "ordermap_POPC-C218-87--POPC-H18R-88_total.dat",
        "ordermap_POPC-C218-87--POPC-H18S-89_total.dat",
        "ordermap_POPC-C218-87--POPC-H18T-90_total.dat",
        "ordermap_POPC-C218-87_total.dat",
        "ordermap_POPC-C22-32--POPC-H2R-33_total.dat",
        "ordermap_POPC-C22-32--POPC-H2S-34_total.dat",
        "ordermap_POPC-C22-32_total.dat",
        "ordermap_POPC-C24-47--POPC-H4R-48_total.dat",
        "ordermap_POPC-C24-47--POPC-H4S-49_total.dat",
        "ordermap_POPC-C24-47_total.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/POPC/{}", path_to_dir, file);
        let test_file = format!("tests/files/ordermaps/{}", file);
        compare_ordermaps(real_file, test_file);
    }

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_small.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
}

#[test]
fn test_aa_order_maps_leaflets() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.0),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "resname POPC and name C22 C24 C218",
                "@membrane and element name hydrogen",
            ))
            .leaflets(method)
            .map(
                OrderMap::new()
                    .bin_size_y(4.0)
                    .bin_size_y(4.0)
                    .output_directory(path_to_dir)
                    .min_samples(5)
                    .build()
                    .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C218-87_lower.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_upper.dat",
            "ordermap_POPC-C22-32_lower.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_upper.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_upper.dat",
            "ordermap_POPC-C218-87--POPC-H18R-88_lower.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_lower.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_lower.dat",
            "ordermap_POPC-C22-32_total.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_lower.dat",
            "ordermap_POPC-C218-87--POPC-H18R-88_total.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_total.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_total.dat",
            "ordermap_POPC-C22-32_upper.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_total.dat",
            "ordermap_POPC-C218-87--POPC-H18R-88_upper.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_upper.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_upper.dat",
            "ordermap_POPC-C24-47_lower.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_upper.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_lower.dat",
            "ordermap_POPC-C218-87_total.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_lower.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_lower.dat",
            "ordermap_POPC-C24-47_total.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_total.dat",
            "ordermap_POPC-C218-87_upper.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_total.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_total.dat",
            "ordermap_POPC-C24-47_upper.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps/{}", file);
            compare_ordermaps(real_file, test_file);
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_leaflets_small.yaml").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_maps_basic_multiple_threads() {
    for n_threads in [3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::new()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "resname POPC and name C22 C24 C218",
                "@membrane and element name hydrogen",
            ))
            .n_threads(n_threads)
            .map(
                OrderMap::new()
                    .bin_size_y(4.0)
                    .bin_size_y(4.0)
                    .output_directory(path_to_dir)
                    .min_samples(5)
                    .build()
                    .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C218-87--POPC-H18R-88_total.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_total.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_total.dat",
            "ordermap_POPC-C218-87_total.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_total.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_total.dat",
            "ordermap_POPC-C22-32_total.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_total.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_total.dat",
            "ordermap_POPC-C24-47_total.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps/{}", file);
            compare_ordermaps(real_file, test_file);
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_small.yaml").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}

#[test]
fn test_aa_order_maps_basic_weird_molecules() {
    // calculation of ordermaps for system with molecules sharing their name and being composed of multiple residues
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/multiple_resid_same_name.tpr")
        .trajectory("tests/files/multiple_resid_same_name.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "resname POPC POPE and name C1A C3A C1B C3B",
            "resname POPC POPE and name D2A C4A C2B C4B",
        ))
        .map(
            OrderMap::new()
                .bin_size_y(4.0)
                .bin_size_y(4.0)
                .output_directory(path_to_dir)
                .min_samples(1)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let expected_file_names = [
        "POPC-POPE1/ordermap_POPC-C1A-4--POPC-D2A-5_total.dat",
        "POPC-POPE1/ordermap_POPC-D2A-5--POPE-C3A-6_total.dat",
        "POPC-POPE1/ordermap_POPE-C3A-6--POPE-C4A-7_total.dat",
        "POPC-POPE1/ordermap_POPE-C1B-8--POPE-C2B-9_total.dat",
        "POPC-POPE1/ordermap_POPE-C2B-9--POPE-C3B-10_total.dat",
        "POPC-POPE1/ordermap_POPE-C3B-10--POPE-C4B-11_total.dat",
        "POPC-POPE1/ordermap_POPC-C1A-4_total.dat",
        "POPC-POPE1/ordermap_POPE-C3A-6_total.dat",
        "POPC-POPE1/ordermap_POPE-C1B-8_total.dat",
        "POPC-POPE1/ordermap_POPE-C3B-10_total.dat",
        "POPC-POPE2/ordermap_POPC-C1A-4--POPC-D2A-5_total.dat",
        "POPC-POPE2/ordermap_POPC-D2A-5--POPE-C3A-6_total.dat",
        "POPC-POPE2/ordermap_POPE-C3A-6--POPE-C4A-7_total.dat",
        "POPC-POPE2/ordermap_POPE-C3B-10--POPE-C4B-11_total.dat",
        "POPC-POPE2/ordermap_POPC-C1A-4_total.dat",
        "POPC-POPE2/ordermap_POPE-C3A-6_total.dat",
        "POPC-POPE2/ordermap_POPE-C3B-10_total.dat",
        "POPC-POPE1/ordermap_POPC-C1A-4--POPC-D2A-5_total.dat",
        "POPC/ordermap_POPC-D2A-5--POPC-C3A-6_total.dat",
        "POPC/ordermap_POPC-C3A-6--POPC-C4A-7_total.dat",
        "POPC/ordermap_POPC-C1B-8--POPC-C2B-9_total.dat",
        "POPC/ordermap_POPC-C2B-9--POPC-C3B-10_total.dat",
        "POPC/ordermap_POPC-C3B-10--POPC-C4B-11_total.dat",
        "POPC/ordermap_POPC-C1A-4_total.dat",
        "POPC/ordermap_POPC-C3A-6_total.dat",
        "POPC/ordermap_POPC-C1B-8_total.dat",
        "POPC/ordermap_POPC-C3B-10_total.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/{}", path_to_dir, file);
        assert!(Path::new(&real_file).exists());
    }
}

#[test]
fn test_aa_order_maps_basic_backup() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let outer_directory = TempDir::new().unwrap();
    let path_to_outer_dir = outer_directory.path().to_str().unwrap();

    let directory = TempDir::new_in(&path_to_outer_dir).unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let backup_file = format!("{}/to_backup.txt", path_to_dir);
    create_file_for_backup!(&backup_file);

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .map(
            OrderMap::new()
                .bin_size_y(4.0)
                .bin_size_y(4.0)
                .output_directory(path_to_dir)
                .min_samples(5)
                .build()
                .unwrap(),
        )
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let expected_file_names = [
        "ordermap_POPC-C218-87--POPC-H18R-88_total.dat",
        "ordermap_POPC-C218-87--POPC-H18S-89_total.dat",
        "ordermap_POPC-C218-87--POPC-H18T-90_total.dat",
        "ordermap_POPC-C218-87_total.dat",
        "ordermap_POPC-C22-32--POPC-H2R-33_total.dat",
        "ordermap_POPC-C22-32--POPC-H2S-34_total.dat",
        "ordermap_POPC-C22-32_total.dat",
        "ordermap_POPC-C24-47--POPC-H4R-48_total.dat",
        "ordermap_POPC-C24-47--POPC-H4S-49_total.dat",
        "ordermap_POPC-C24-47_total.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/POPC/{}", path_to_dir, file);
        let test_file = format!("tests/files/ordermaps/{}", file);
        compare_ordermaps(real_file, test_file);
    }

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_small.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));

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
