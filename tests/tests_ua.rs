// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Integration tests for the calculation of united-atom order parameters.

use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use gorder::prelude::*;
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

#[test]
fn test_ua_order_from_aa_ignore() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .analysis_type(AnalysisType::uaorder(
            Some("@membrane and element name carbon and not name C29 C210 C21 C31"),
            Some("@membrane and name C29 C210"),
            Some("element name hydrogen"),
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/ua_order_from_aa.yaml",
        1
    ));
}

#[test]
fn test_ua_order_leaflets() {
    for n_threads in [1, 2, 3, 8, 64] {
        for method in [
            LeafletClassification::global("@membrane", "name r'^P'"),
            LeafletClassification::local("@membrane", "name r'^P'", 2.5),
            LeafletClassification::individual(
                "name r'^P'",
                "(resname POPC and name CA2 C50) or (resname POPS and name C36 C55)",
            ),
        ] {
            for freq in [
                Frequency::every(1).unwrap(),
                Frequency::every(5).unwrap(),
                Frequency::every(100).unwrap(),
                Frequency::once(),
            ] {
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
                .leaflets(method.clone().with_frequency(freq))
                .n_threads(n_threads)
                .silent()
                .overwrite()
                .build()
                .unwrap();

                analysis.run().unwrap().write().unwrap();

                assert!(diff_files_ignore_first(
                    path_to_yaml,
                    "tests/files/ua_order_leaflets.yaml",
                    1
                ));

                assert!(diff_files_ignore_first(
                    path_to_table,
                    "tests/files/ua_order_leaflets.tab",
                    1
                ));

                assert!(diff_files_ignore_first(
                    path_to_csv,
                    "tests/files/ua_order_leaflets.csv",
                    0
                ));

                for molecule in ["POPC", "POPS"] {
                    let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
                    let path_expected = format!("tests/files/ua_order_leaflets_{}.xvg", molecule);

                    assert!(diff_files_ignore_first(&path, &path_expected, 1));
                }
            }
        }
    }
}

#[test]
fn test_ua_order_begin_end_step() {
    for n_threads in [1, 2, 4, 8, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_yaml = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/ua.tpr")
            .trajectory("tests/files/ua.xtc")
            .output_yaml(path_to_yaml)
            .analysis_type(AnalysisType::uaorder(
                Some("(resname POPC and name r'^C' and not name C15 C34 C24 C25) or (resname POPS and name r'^C' and not name C6 C18 C39 C27 C28)"),
                Some("(resname POPC and name C24 C25) or (resname POPS and name C27 C28)"),
                None
            ))
            .n_threads(n_threads)
            .leaflets(LeafletClassification::global("@membrane", "name r'^P'"))
            .begin(199200.0)
            .end(199800.0)
            .step(3)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        let results = analysis.run().unwrap();

        assert_eq!(results.n_analyzed_frames(), 11);
        results.write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_yaml,
            "tests/files/ua_order_begin_end_step.yaml",
            1
        ));
    }
}

#[test]
fn test_ua_order_nothing_to_analyze() {
    let analysis = Analysis::builder()
        .structure("tests/files/ua.tpr")
        .trajectory("tests/files/ua.xtc")
        .output("THIS_FILE_SHOULD_NOT_BE_CREATED_UA_1")
        .analysis_type(AnalysisType::uaorder(Some("@water"), Some("name Cs"), None))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_UA_1").exists());
}

#[test]
fn test_ua_order_maps_basic() {
    for n_threads in [1, 3, 8, 32] {
        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/ua.tpr")
            .trajectory("tests/files/ua.xtc")
            .analysis_type(AnalysisType::uaorder(
                Some("resname POPC and name C50 C20 C13"),
                Some("resname POPC and name C24"),
                None,
            ))
            .map(
                OrderMap::builder()
                    .bin_size([0.5, 2.0])
                    .output_directory(path_to_dir)
                    .min_samples(5)
                    .build()
                    .unwrap(),
            )
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C13-12_full.dat",
            "ordermap_POPC-C13-12--POPC-H1-12_full.dat",
            "ordermap_POPC-C20-19_full.dat",
            "ordermap_POPC-C20-19--POPC-H1-19_full.dat",
            "ordermap_POPC-C20-19--POPC-H2-19_full.dat",
            "ordermap_POPC-C24-23_full.dat",
            "ordermap_POPC-C24-23--POPC-H1-23_full.dat",
            "ordermap_POPC-C50-49_full.dat",
            "ordermap_POPC-C50-49--POPC-H1-49_full.dat",
            "ordermap_POPC-C50-49--POPC-H2-49_full.dat",
            "ordermap_POPC-C50-49--POPC-H3-49_full.dat",
            "ordermap_average_full.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps_ua/{}", file);
            assert!(diff_files_ignore_first(&real_file, &test_file, 2));
        }

        // full map for the entire system is the same as for POPC
        let real_file = format!("{}/ordermap_average_full.dat", path_to_dir);
        let test_file = "tests/files/ordermaps_ua/ordermap_average_full.dat";
        assert!(diff_files_ignore_first(&real_file, test_file, 2));

        // check the script
        let real_script = format!("{}/plot.py", path_to_dir);
        assert!(diff_files_ignore_first(&real_script, "scripts/plot.py", 0));
    }
}

#[test]
fn test_ua_order_maps_leaflets() {
    for n_threads in [1, 3, 8, 32] {
        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/ua.tpr")
            .trajectory("tests/files/ua.xtc")
            .analysis_type(AnalysisType::uaorder(
                Some("resname POPC and name C50 C20 C13"),
                Some("resname POPC and name C24"),
                None,
            ))
            .map(
                OrderMap::builder()
                    .bin_size([0.5, 2.0])
                    .output_directory(path_to_dir)
                    .min_samples(5)
                    .build()
                    .unwrap(),
            )
            .leaflets(LeafletClassification::global("@membrane", "name r'^P'"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C13-12_upper.dat",
            "ordermap_POPC-C20-19--POPC-H2-19_upper.dat",
            "ordermap_POPC-C50-49_full.dat",
            "ordermap_POPC-C50-49--POPC-H3-49_full.dat",
            "ordermap_POPC-C20-19_full.dat",
            "ordermap_POPC-C20-19_upper.dat",
            "ordermap_POPC-C50-49_lower.dat",
            "ordermap_POPC-C50-49--POPC-H3-49_lower.dat",
            "ordermap_POPC-C20-19_lower.dat",
            "ordermap_POPC-C24-23_full.dat",
            "ordermap_POPC-C50-49--POPC-H1-49_full.dat",
            "ordermap_POPC-C50-49--POPC-H3-49_upper.dat",
            "ordermap_POPC-C13-12_full.dat",
            "ordermap_POPC-C20-19--POPC-H1-19_full.dat",
            "ordermap_POPC-C24-23_lower.dat",
            "ordermap_POPC-C50-49--POPC-H1-49_lower.dat",
            "ordermap_POPC-C50-49_upper.dat",
            "ordermap_POPC-C13-12_lower.dat",
            "ordermap_POPC-C20-19--POPC-H1-19_lower.dat",
            "ordermap_POPC-C24-23--POPC-H1-23_full.dat",
            "ordermap_POPC-C50-49--POPC-H1-49_upper.dat",
            "ordermap_POPC-C13-12--POPC-H1-12_full.dat",
            "ordermap_POPC-C20-19--POPC-H1-19_upper.dat",
            "ordermap_POPC-C24-23--POPC-H1-23_lower.dat",
            "ordermap_POPC-C50-49--POPC-H2-49_full.dat",
            "ordermap_POPC-C13-12--POPC-H1-12_lower.dat",
            "ordermap_POPC-C20-19--POPC-H2-19_full.dat",
            "ordermap_POPC-C24-23--POPC-H1-23_upper.dat",
            "ordermap_POPC-C50-49--POPC-H2-49_lower.dat",
            "ordermap_POPC-C13-12--POPC-H1-12_upper.dat",
            "ordermap_POPC-C20-19--POPC-H2-19_lower.dat",
            "ordermap_POPC-C24-23_upper.dat",
            "ordermap_POPC-C50-49--POPC-H2-49_upper.dat",
            "ordermap_average_full.dat",
            "ordermap_average_upper.dat",
            "ordermap_average_lower.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps_ua/{}", file);
            assert!(diff_files_ignore_first(&real_file, &test_file, 2));
        }

        // full map for the entire system is the same as for POPC
        let real_file = format!("{}/ordermap_average_full.dat", path_to_dir);
        let test_file = "tests/files/ordermaps_ua/ordermap_average_full.dat";
        assert!(diff_files_ignore_first(&real_file, test_file, 2));

        // check the script
        let real_script = format!("{}/plot.py", path_to_dir);
        assert!(diff_files_ignore_first(&real_script, "scripts/plot.py", 0));
    }
}

#[test]
fn test_ua_order_error_convergence() {
    for n_threads in [1, 3, 8, 32] {
        let output = NamedTempFile::new().unwrap();
        let path_to_yaml = output.path().to_str().unwrap();

        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let output_convergence = NamedTempFile::new().unwrap();
        let path_to_convergence = output_convergence.path().to_str().unwrap();

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
            .estimate_error(EstimateError::new(None, Some(path_to_convergence)).unwrap())
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_yaml,
            "tests/files/ua_order_error.yaml",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_table,
            "tests/files/ua_order_error.tab",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/ua_order_error.csv",
            0
        ));

        for molecule in ["POPC", "POPS"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!("tests/files/ua_order_basic_{}.xvg", molecule);

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }

        assert!(diff_files_ignore_first(
            path_to_convergence,
            "tests/files/ua_order_convergence.xvg",
            1
        ));
    }
}

#[test]
fn test_ua_order_error_leaflets_convergence() {
    for n_threads in [1, 3, 8, 32] {
        let output = NamedTempFile::new().unwrap();
        let path_to_yaml = output.path().to_str().unwrap();

        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let output_convergence = NamedTempFile::new().unwrap();
        let path_to_convergence = output_convergence.path().to_str().unwrap();

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
            .estimate_error(EstimateError::new(None, Some(path_to_convergence)).unwrap())
            .leaflets(LeafletClassification::global("@membrane", "name r'^P'"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_yaml,
            "tests/files/ua_order_leaflets_error.yaml",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_table,
            "tests/files/ua_order_leaflets_error.tab",
            1
        ));

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/ua_order_leaflets_error.csv",
            0
        ));

        for molecule in ["POPC", "POPS"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!("tests/files/ua_order_leaflets_{}.xvg", molecule);

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }

        assert!(diff_files_ignore_first(
            path_to_convergence,
            "tests/files/ua_order_leaflets_convergence.xvg",
            1
        ));
    }
}
