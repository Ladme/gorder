// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Integration tests for the `gorder` library.

use std::{fs::File, path::Path};

use gorder::prelude::*;
use tempfile::NamedTempFile;

#[test]
fn test_aa_order_basic_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::new()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::AAOrder)
        .heavy_atoms("@membrane and element name carbon")
        .hydrogens("@membrane and element name hydrogen")
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_basic.yaml").unwrap();

    assert!(file_diff::diff_files(&mut result, &mut expected));
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
            .analysis_type(AnalysisType::AAOrder)
            .heavy_atoms("@membrane and element name carbon")
            .hydrogens("@membrane and element name hydrogen")
            .n_threads(n_threads)
            .silent()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_basic.yaml").unwrap();

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
            .analysis_type(AnalysisType::AAOrder)
            .heavy_atoms("@membrane and element name carbon")
            .hydrogens("@membrane and element name hydrogen")
            .leaflets(method)
            .silent()
            .build()
            .unwrap();

        analysis.run().unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("tests/files/aa_order_leaflets.yaml").unwrap();

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
        .analysis_type(AnalysisType::AAOrder)
        .heavy_atoms("@membrane and element name carbon")
        .hydrogens("@membrane and element name hydrogen")
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_selected.yaml").unwrap();

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
        .analysis_type(AnalysisType::AAOrder)
        .heavy_atoms("@membrane and element name carbon")
        .hydrogens("@membrane and element name hydrogen")
        .min_samples(2000)
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap();

    let mut result = File::open(path_to_output).unwrap();
    let mut expected = File::open("tests/files/aa_order_limit.yaml").unwrap();

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
        .analysis_type(AnalysisType::AAOrder)
        .heavy_atoms("@membrane and element name carbon")
        .hydrogens("@membrane and element name hydrogen")
        .begin(450_000.0)
        .end(450_200.0)
        .step(3)
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
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
        .analysis_type(AnalysisType::AAOrder)
        .heavy_atoms("@ion")
        .hydrogens("@membrane and element name hydrogen")
        .silent()
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
        .analysis_type(AnalysisType::AAOrder)
        .heavy_atoms("@water and element symbol O")
        .hydrogens("@membrane and element name hydrogen")
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_2").exists());
}
