// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Tests of the binary application.

use std::io::BufRead;
use std::{fs::File, io::BufReader};

use assert_cmd::Command;

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
fn test_bin_aa_order_basic_yaml() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/basic_aa.yaml",
            "--silent",
            "--overwrite",
        ])
        .assert()
        .success()
        .stdout("");

    assert!(diff_files_ignore_first(
        "temp_aa_order.yaml",
        "tests/files/aa_order_basic.yaml",
        1
    ));

    std::fs::remove_file("temp_aa_order.yaml").unwrap();
}

#[test]
fn test_bin_cg_order_leaflets_yaml_tab() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args(["tests/files/inputs/leaflets_cg.yaml", "--overwrite"])
        .assert()
        .success()
        .stdout("");

    assert!(diff_files_ignore_first(
        "temp_cg_order.yaml",
        "tests/files/cg_order_leaflets.yaml",
        1
    ));

    assert!(diff_files_ignore_first(
        "temp_cg_order.tab",
        "tests/files/cg_order_leaflets.tab",
        1
    ));

    std::fs::remove_file("temp_cg_order.yaml").unwrap();
    std::fs::remove_file("temp_cg_order.tab").unwrap();
}

#[test]
fn test_bin_cg_order_maps() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args(["tests/files/inputs/maps_cg.yaml", "--overwrite"])
        .assert()
        .success()
        .stdout("");

    assert!(diff_files_ignore_first(
        "temp_cg_order_maps.yaml",
        "tests/files/cg_order_small.yaml",
        1
    ));
    std::fs::remove_file("temp_cg_order_maps.yaml").unwrap();

    let expected_file_names = [
        "ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
        "ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
        "ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
        "ordermap_average_full.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("temp_cg_ordermaps/POPC/{}", file);
        let test_file = format!("tests/files/ordermaps_cg/{}", file);
        assert!(diff_files_ignore_first(&real_file, &test_file, 2));
    }

    std::fs::remove_dir_all("temp_cg_ordermaps").unwrap();
}

#[test]
fn test_bin_estimate_error() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args(["tests/files/inputs/estimate_error_cg.yaml", "--overwrite"])
        .assert()
        .success()
        .stdout("");

    assert!(diff_files_ignore_first(
        "temp_cg_ee.yaml",
        "tests/files/cg_order_error.yaml",
        1
    ));

    assert!(diff_files_ignore_first(
        "temp_cg_ee.tab",
        "tests/files/cg_order_error.tab",
        1
    ));

    assert!(diff_files_ignore_first(
        "temp_cg_ee.csv",
        "tests/files/cg_order_error.csv",
        0
    ));

    std::fs::remove_file("temp_cg_ee.yaml").unwrap();
    std::fs::remove_file("temp_cg_ee.tab").unwrap();
    std::fs::remove_file("temp_cg_ee.csv").unwrap();
}

#[test]
fn test_bin_cg_order_maps_export_config() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/maps_cg_for_export_config.yaml",
            "--overwrite",
            "--export-config",
            "temp_analysis_out.yaml",
        ])
        .assert()
        .success()
        .stdout("");

    // remove everything and rerun the analysis using output config to check that it works
    std::fs::remove_file("temp_cg_order_maps_for_export_config.yaml").unwrap();
    std::fs::remove_dir_all("temp_cg_ordermaps_for_export_config").unwrap();

    Command::cargo_bin("gorder")
        .unwrap()
        .args(["temp_analysis_out.yaml", "--overwrite"])
        .assert()
        .success()
        .stdout("");

    assert!(diff_files_ignore_first(
        "temp_cg_order_maps_for_export_config.yaml",
        "tests/files/cg_order_small.yaml",
        1
    ));
    std::fs::remove_file("temp_cg_order_maps_for_export_config.yaml").unwrap();

    let expected_file_names = [
        "ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
        "ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
        "ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
        "ordermap_average_full.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("temp_cg_ordermaps_for_export_config/POPC/{}", file);
        let test_file = format!("tests/files/ordermaps_cg/{}", file);
        assert!(diff_files_ignore_first(&real_file, &test_file, 2));
    }

    std::fs::remove_dir_all("temp_cg_ordermaps_for_export_config").unwrap();
    std::fs::remove_file("temp_analysis_out.yaml").unwrap();
}

#[test]
fn test_bin_aa_order_writing_fail() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/writing_fail.yaml",
            "--silent",
            "--overwrite",
        ])
        .assert()
        .failure()
        .stdout("")
        .stderr(
            "[E] error: could not create file \'this_directory_does_not_exist/temp_aa_order.yaml\'\n",
        );
}

#[test]
fn test_bin_aa_order_fail() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/atom_overlap_error.yaml",
            "--silent",
            "--overwrite",
        ])
        .assert()
        .failure()
        .stdout("")
        .stderr("[E] error: 217 atoms are part of both \'HeavyAtoms\' (query: \'@membrane and element name carbon or serial 876 to 1234\') and \'Hydrogens\' (query: \'@membrane and element name hydrogen\')\n");
}

#[test]
fn test_bin_missing_output_fail() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/basic.yaml",
            "--silent",
            "--overwrite",
        ])
        .assert()
        .failure()
        .stdout("")
        .stderr("[E] error: no yaml output file specified in the configuration file \'tests/files/inputs/basic.yaml\' (hint: add \'output: output.yaml\' to your configuration file)\n");
}

#[test]
fn test_bin_missing_maps_output_fail() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/default_ordermap.yaml",
            "--silent",
            "--overwrite",
        ])
        .assert()
        .failure()
        .stdout("")
        .stderr("[E] error: no output directory for ordermaps specified in the configuration file \'tests/files/inputs/default_ordermap.yaml\'\n");
}

#[test]
fn test_bin_output_config_writing_fails() {
    Command::cargo_bin("gorder")
        .unwrap()
        .args([
            "tests/files/inputs/basic_aa_config_fails.yaml",
            "--silent",
            "--overwrite",
            "--export-config",
            "this_directory_does_not_exist/analysis_out.yaml",
        ])
        .assert()
        .success()
        .stdout("")
        .stderr(
            "[E] Analysis completed successfully, but exporting the analysis options failed!
 |      error: could not create file \'this_directory_does_not_exist/analysis_out.yaml\'\n",
        );

    assert!(diff_files_ignore_first(
        "temp_aa_order_config_fails.yaml",
        "tests/files/aa_order_basic.yaml",
        1
    ));

    std::fs::remove_file("temp_aa_order_config_fails.yaml").unwrap();
}
