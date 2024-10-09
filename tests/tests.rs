// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Integration tests for the `gorder` library.

use std::{fs::File, path::Path};

use approx::assert_relative_eq;
use gorder::prelude::*;
use groan_rs::prelude::GridMap;
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
        let real_file = format!("{}/{}", path_to_dir, file);
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
            let real_file = format!("{}/{}", path_to_dir, file);
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
            let real_file = format!("{}/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps/{}", file);
            compare_ordermaps(real_file, test_file);
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_small.yaml").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}
