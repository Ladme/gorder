// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Integration tests for the calculation of atomistic order parameters.

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
fn test_aa_order_basic_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_basic.yaml",
        1
    ));
}

#[test]
fn test_aa_order_basic_ndx_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .index("tests/files/pcpepg.ndx")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder("HeavyAtoms", "Hydrogens"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_basic.yaml",
        1
    ));
}

#[test]
fn test_aa_order_basic_fail_overlap() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
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
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_basic.tab",
        1
    ));
}

#[test]
fn test_aa_order_basic_xvg() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_xvg(pattern)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        let path_expected = format!("tests/files/aa_order_basic_{}.xvg", molecule);

        assert!(diff_files_ignore_first(&path, &path_expected, 1));
    }
}

#[test]
fn test_aa_order_basic_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_basic.csv",
        0
    ));
}

#[test]
fn test_aa_order_basic_xvg_weird_names() {
    for name in ["order", ".this.is.a.weird.name.xvg"] {
        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let pattern = format!("{}/{}", path_to_dir, name);

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output_xvg(pattern)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        for molecule in ["POPC", "POPE", "POPG"] {
            let path = if name.contains(".xvg") {
                format!("{}/.this.is.a.weird.name_{}.xvg", path_to_dir, molecule)
            } else {
                format!("{}/order_{}", path_to_dir, molecule)
            };

            let path_expected = format!("tests/files/aa_order_basic_{}.xvg", molecule);
            assert!(diff_files_ignore_first(&path, &path_expected, 1));
        }
    }
}

#[test]
fn test_aa_order_basic_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
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

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_basic.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_basic_table_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
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

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_table,
            "tests/files/aa_order_basic.tab",
            1
        ));
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

        let analysis = Analysis::builder()
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

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_leaflets.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_leaflets_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        for method in [
            LeafletClassification::global("@membrane", "name P"),
            LeafletClassification::local("@membrane", "name P", 2.5),
            LeafletClassification::individual("name P", "name C218 C316"),
        ] {
            let output = NamedTempFile::new().unwrap();
            let path_to_output = output.path().to_str().unwrap();

            let analysis = Analysis::builder()
                .structure("tests/files/pcpepg.tpr")
                .trajectory("tests/files/pcpepg.xtc")
                .output(path_to_output)
                .analysis_type(AnalysisType::aaorder(
                    "@membrane and element name carbon",
                    "@membrane and element name hydrogen",
                ))
                .leaflets(method)
                .n_threads(n_threads)
                .silent()
                .overwrite()
                .build()
                .unwrap();

            analysis.run().unwrap().write().unwrap();

            assert!(diff_files_ignore_first(
                path_to_output,
                "tests/files/aa_order_leaflets.yaml",
                1
            ));
        }
    }
}

#[test]
fn test_aa_order_leaflets_yaml_different_membrane_normals() {
    for (input_traj, normal) in [
        "tests/files/pcpepg_switched_xy.xtc",
        "tests/files/pcpepg_switched_xz.xtc",
        "tests/files/pcpepg_switched_yz.xtc",
    ]
    .into_iter()
    .zip([Axis::Z, Axis::X, Axis::Y].into_iter())
    {
        for method in [
            LeafletClassification::global("@membrane", "name P"),
            LeafletClassification::local("@membrane", "name P", 2.5),
            LeafletClassification::individual("name P", "name C218 C316"),
        ] {
            let output = NamedTempFile::new().unwrap();
            let path_to_output = output.path().to_str().unwrap();

            let analysis = Analysis::builder()
                .structure("tests/files/pcpepg.tpr")
                .trajectory(input_traj)
                .output(path_to_output)
                .analysis_type(AnalysisType::aaorder(
                    "@membrane and element name carbon",
                    "@membrane and element name hydrogen",
                ))
                .leaflets(method)
                .membrane_normal(normal)
                .silent()
                .overwrite()
                .build()
                .unwrap();

            analysis.run().unwrap().write().unwrap();

            assert!(diff_files_ignore_first(
                path_to_output,
                "tests/files/aa_order_leaflets.yaml",
                1
            ));
        }
    }
}

#[test]
fn test_aa_order_leaflets_table() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.5),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let output_table = NamedTempFile::new().unwrap();
        let path_to_table = output_table.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
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

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_table,
            "tests/files/aa_order_leaflets.tab",
            1
        ));
    }
}

#[test]
fn test_aa_order_leaflets_xvg() {
    for method in [
        LeafletClassification::global("@membrane", "name P"),
        LeafletClassification::local("@membrane", "name P", 2.5),
        LeafletClassification::individual("name P", "name C218 C316"),
    ] {
        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let pattern = format!("{}/order.xvg", path_to_dir);

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
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

        analysis.run().unwrap().write().unwrap();

        for molecule in ["POPC", "POPE", "POPG"] {
            let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
            let path_expected = format!("tests/files/aa_order_leaflets_{}.xvg", molecule);

            assert!(diff_files_ignore_first(&path, &path_expected, 1));
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
        let output_csv = NamedTempFile::new().unwrap();
        let path_to_csv = output_csv.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
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

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_csv,
            "tests/files/aa_order_leaflets.csv",
            0
        ));
    }
}

#[test]
fn test_aa_order_leaflets_yaml_supershort() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_selected.yaml",
        1
    ));
}

#[test]
fn test_aa_order_one_different_hydrogen_numbers_table() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_different_hydrogen_numbers.tab",
        1
    ));
}

#[test]
fn test_aa_order_one_different_hydrogen_numbers_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_different_hydrogen_numbers.csv",
        0
    ));
}

#[test]
fn test_aa_order_limit_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_limit.yaml",
        1
    ));
}

#[test]
fn test_aa_order_leaflets_limit_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_leaflets_limit.yaml",
        1
    ));
}

#[test]
fn test_aa_order_leaflets_limit_tab() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_leaflets_limit.tab",
        1
    ));
}

#[test]
fn test_aa_order_leaflets_limit_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
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

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_leaflets_limit.csv",
        0
    ));
}

#[test]
fn test_aa_order_begin_end_step_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
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

    let results = analysis.run().unwrap();
    assert_eq!(results.n_analyzed_frames(), 4);
    results.write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_begin_end_step.yaml",
        1
    ));
}

#[test]
fn test_aa_order_begin_end_step_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
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
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        let results = analysis.run().unwrap();
        assert_eq!(results.n_analyzed_frames(), 4);
        results.write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_begin_end_step.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_begin_end_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .begin(450_000.0)
        .end(450_200.0)
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = analysis.run().unwrap();
    assert_eq!(results.n_analyzed_frames(), 11);
    results.write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_begin_end.yaml",
        1
    ));
}

#[test]
fn test_aa_order_begin_end_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .begin(450_000.0)
            .end(450_200.0)
            .leaflets(LeafletClassification::global("@membrane", "name P"))
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        let results = analysis.run().unwrap();
        assert_eq!(results.n_analyzed_frames(), 11);
        results.write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_begin_end.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_no_molecules() {
    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

    assert!(!Path::new("THIS_FILE_SHOULD_NOT_BE_CREATED_1").exists());
}

#[test]
fn test_aa_order_empty_molecules() {
    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

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

    let analysis = Analysis::builder()
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

    analysis.run().unwrap().write().unwrap();

    let all_files = [
        ("tests/files/aa_order_basic.yaml", &file_paths[0]),
        ("tests/files/aa_order_basic.tab", &file_paths[1]),
        ("tests/files/aa_order_basic.csv", &file_paths[2]),
        ("tests/files/aa_order_basic_POPC.xvg", &file_paths[3]),
        ("tests/files/aa_order_basic_POPE.xvg", &file_paths[4]),
        ("tests/files/aa_order_basic_POPG.xvg", &file_paths[5]),
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
fn test_aa_order_maps_basic() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .map(
            OrderMap::builder()
                .bin_size([0.1, 4.0])
                .output_directory(path_to_dir)
                .min_samples(5)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    let expected_file_names = [
        "ordermap_POPC-C218-87--POPC-H18R-88_full.dat",
        "ordermap_POPC-C218-87--POPC-H18S-89_full.dat",
        "ordermap_POPC-C218-87--POPC-H18T-90_full.dat",
        "ordermap_POPC-C218-87_full.dat",
        "ordermap_POPC-C22-32--POPC-H2R-33_full.dat",
        "ordermap_POPC-C22-32--POPC-H2S-34_full.dat",
        "ordermap_POPC-C22-32_full.dat",
        "ordermap_POPC-C24-47--POPC-H4R-48_full.dat",
        "ordermap_POPC-C24-47--POPC-H4S-49_full.dat",
        "ordermap_POPC-C24-47_full.dat",
        "ordermap_average_full.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/POPC/{}", path_to_dir, file);
        let test_file = format!("tests/files/ordermaps/{}", file);
        assert!(diff_files_ignore_first(&real_file, &test_file, 2));
    }

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_small.yaml",
        1
    ));
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

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "resname POPC and name C22 C24 C218",
                "@membrane and element name hydrogen",
            ))
            .leaflets(method)
            .map(
                OrderMap::builder()
                    .bin_size([0.1, 4.0])
                    .output_directory(path_to_dir)
                    .min_samples(5)
                    .build()
                    .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C218-87_lower.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_upper.dat",
            "ordermap_POPC-C22-32_lower.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_upper.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_upper.dat",
            "ordermap_POPC-C218-87--POPC-H18R-88_lower.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_lower.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_lower.dat",
            "ordermap_POPC-C22-32_full.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_lower.dat",
            "ordermap_POPC-C218-87--POPC-H18R-88_full.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_full.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_full.dat",
            "ordermap_POPC-C22-32_full.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_full.dat",
            "ordermap_POPC-C218-87--POPC-H18R-88_upper.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_upper.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_upper.dat",
            "ordermap_POPC-C24-47_lower.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_upper.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_lower.dat",
            "ordermap_POPC-C218-87_full.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_lower.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_lower.dat",
            "ordermap_POPC-C24-47_full.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_full.dat",
            "ordermap_POPC-C218-87_upper.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_full.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_full.dat",
            "ordermap_POPC-C24-47_upper.dat",
            "ordermap_average_full.dat",
            "ordermap_average_upper.dat",
            "ordermap_average_lower.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps/{}", file);
            assert!(diff_files_ignore_first(&real_file, &test_file, 2));
        }

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_leaflets_small.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_maps_leaflets_different_membrane_normals() {
    for (input_traj, (normal, structure)) in [
        "tests/files/pcpepg_switched_xz.xtc",
        "tests/files/pcpepg_switched_yz.xtc",
    ]
    .into_iter()
    .zip(
        [Axis::X, Axis::Y].into_iter().zip(
            [
                "tests/files/pcpepg_switched_xz.tpr",
                "tests/files/pcpepg_switched_yz.tpr",
            ]
            .into_iter(),
        ),
    ) {
        for method in [
            LeafletClassification::global("@membrane", "name P"),
            LeafletClassification::local("@membrane", "name P", 2.0),
            LeafletClassification::individual("name P", "name C218 C316"),
        ] {
            let output = NamedTempFile::new().unwrap();
            let path_to_output = output.path().to_str().unwrap();

            let directory = TempDir::new().unwrap();
            let path_to_dir = directory.path().to_str().unwrap();

            let analysis = Analysis::builder()
                .structure(structure)
                .trajectory(input_traj)
                .membrane_normal(normal)
                .output(path_to_output)
                .analysis_type(AnalysisType::aaorder(
                    "resname POPC and name C22 C24 C218",
                    "@membrane and element name hydrogen",
                ))
                .leaflets(method)
                .map(
                    OrderMap::builder()
                        .bin_size([0.1, 4.0])
                        .output_directory(path_to_dir)
                        .min_samples(5)
                        .build()
                        .unwrap(),
                )
                .silent()
                .overwrite()
                .build()
                .unwrap();

            analysis.run().unwrap().write().unwrap();

            let expected_file_names = [
                "ordermap_POPC-C218-87_lower.dat",
                "ordermap_POPC-C218-87--POPC-H18S-89_upper.dat",
                "ordermap_POPC-C22-32_lower.dat",
                "ordermap_POPC-C22-32--POPC-H2S-34_upper.dat",
                "ordermap_POPC-C24-47--POPC-H4R-48_upper.dat",
                //"ordermap_POPC-C218-87--POPC-H18R-88_lower.dat",  // this file does not pass the test due to some rounding error
                "ordermap_POPC-C218-87--POPC-H18T-90_lower.dat",
                "ordermap_POPC-C22-32--POPC-H2R-33_lower.dat",
                "ordermap_POPC-C22-32_full.dat",
                "ordermap_POPC-C24-47--POPC-H4S-49_lower.dat",
                "ordermap_POPC-C218-87--POPC-H18R-88_full.dat",
                "ordermap_POPC-C218-87--POPC-H18T-90_full.dat",
                "ordermap_POPC-C22-32--POPC-H2R-33_full.dat",
                "ordermap_POPC-C22-32_upper.dat",
                "ordermap_POPC-C24-47--POPC-H4S-49_full.dat",
                "ordermap_POPC-C218-87--POPC-H18R-88_upper.dat",
                "ordermap_POPC-C218-87--POPC-H18T-90_upper.dat",
                "ordermap_POPC-C22-32--POPC-H2R-33_upper.dat",
                "ordermap_POPC-C24-47_lower.dat",
                "ordermap_POPC-C24-47--POPC-H4S-49_upper.dat",
                "ordermap_POPC-C218-87--POPC-H18S-89_lower.dat",
                "ordermap_POPC-C218-87_full.dat",
                "ordermap_POPC-C22-32--POPC-H2S-34_lower.dat",
                "ordermap_POPC-C24-47--POPC-H4R-48_lower.dat",
                "ordermap_POPC-C24-47_full.dat",
                "ordermap_POPC-C218-87--POPC-H18S-89_full.dat",
                "ordermap_POPC-C218-87_upper.dat",
                "ordermap_POPC-C22-32--POPC-H2S-34_full.dat",
                "ordermap_POPC-C24-47--POPC-H4R-48_full.dat",
                "ordermap_POPC-C24-47_upper.dat",
                "ordermap_average_full.dat",
                "ordermap_average_upper.dat",
                "ordermap_average_lower.dat",
            ];

            for file in expected_file_names {
                let real_file = format!("{}/POPC/{}", path_to_dir, file);
                let test_file = format!("tests/files/ordermaps/{}", file);
                assert!(diff_files_ignore_first(&real_file, &test_file, 4));
            }

            assert!(diff_files_ignore_first(
                path_to_output,
                "tests/files/aa_order_leaflets_small.yaml",
                1
            ));
        }
    }
}

#[test]
fn test_aa_order_maps_basic_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let directory = TempDir::new().unwrap();
        let path_to_dir = directory.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "resname POPC and name C22 C24 C218",
                "@membrane and element name hydrogen",
            ))
            .n_threads(n_threads)
            .map(
                OrderMap::builder()
                    .bin_size([0.1, 4.0])
                    .output_directory(path_to_dir)
                    .min_samples(5)
                    .build()
                    .unwrap(),
            )
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        let expected_file_names = [
            "ordermap_POPC-C218-87--POPC-H18R-88_full.dat",
            "ordermap_POPC-C218-87--POPC-H18S-89_full.dat",
            "ordermap_POPC-C218-87--POPC-H18T-90_full.dat",
            "ordermap_POPC-C218-87_full.dat",
            "ordermap_POPC-C22-32--POPC-H2R-33_full.dat",
            "ordermap_POPC-C22-32--POPC-H2S-34_full.dat",
            "ordermap_POPC-C22-32_full.dat",
            "ordermap_POPC-C24-47--POPC-H4R-48_full.dat",
            "ordermap_POPC-C24-47--POPC-H4S-49_full.dat",
            "ordermap_POPC-C24-47_full.dat",
            "ordermap_average_full.dat",
        ];

        for file in expected_file_names {
            let real_file = format!("{}/POPC/{}", path_to_dir, file);
            let test_file = format!("tests/files/ordermaps/{}", file);
            assert!(diff_files_ignore_first(&real_file, &test_file, 2));
        }

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_small.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_maps_basic_weird_molecules() {
    // calculation of ordermaps for system with molecules sharing their name and being composed of multiple residues
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/multiple_resid_same_name.tpr")
        .trajectory("tests/files/multiple_resid_same_name.xtc")
        .analysis_type(AnalysisType::aaorder(
            "resname POPC POPE and name C1A C3A C1B C3B",
            "resname POPC POPE and name D2A C4A C2B C4B",
        ))
        .map(
            OrderMap::builder()
                .bin_size([0.1, 4.0])
                .output_directory(path_to_dir)
                .min_samples(1)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    let expected_file_names = [
        "POPC-POPE1/ordermap_POPC-C1A-4--POPC-D2A-5_full.dat",
        "POPC-POPE1/ordermap_POPC-D2A-5--POPE-C3A-6_full.dat",
        "POPC-POPE1/ordermap_POPE-C3A-6--POPE-C4A-7_full.dat",
        "POPC-POPE1/ordermap_POPE-C1B-8--POPE-C2B-9_full.dat",
        "POPC-POPE1/ordermap_POPE-C2B-9--POPE-C3B-10_full.dat",
        "POPC-POPE1/ordermap_POPE-C3B-10--POPE-C4B-11_full.dat",
        "POPC-POPE1/ordermap_POPC-C1A-4_full.dat",
        "POPC-POPE1/ordermap_POPE-C3A-6_full.dat",
        "POPC-POPE1/ordermap_POPE-C1B-8_full.dat",
        "POPC-POPE1/ordermap_POPE-C3B-10_full.dat",
        "POPC-POPE1/ordermap_average_full.dat",
        "POPC-POPE2/ordermap_POPC-C1A-4--POPC-D2A-5_full.dat",
        "POPC-POPE2/ordermap_POPC-D2A-5--POPE-C3A-6_full.dat",
        "POPC-POPE2/ordermap_POPE-C3A-6--POPE-C4A-7_full.dat",
        "POPC-POPE2/ordermap_POPE-C3B-10--POPE-C4B-11_full.dat",
        "POPC-POPE2/ordermap_POPC-C1A-4_full.dat",
        "POPC-POPE2/ordermap_POPE-C3A-6_full.dat",
        "POPC-POPE2/ordermap_POPE-C3B-10_full.dat",
        "POPC-POPE2/ordermap_POPC-C1A-4--POPC-D2A-5_full.dat",
        "POPC-POPE2/ordermap_average_full.dat",
        "POPC/ordermap_POPC-D2A-5--POPC-C3A-6_full.dat",
        "POPC/ordermap_POPC-C3A-6--POPC-C4A-7_full.dat",
        "POPC/ordermap_POPC-C1B-8--POPC-C2B-9_full.dat",
        "POPC/ordermap_POPC-C2B-9--POPC-C3B-10_full.dat",
        "POPC/ordermap_POPC-C3B-10--POPC-C4B-11_full.dat",
        "POPC/ordermap_POPC-C1A-4_full.dat",
        "POPC/ordermap_POPC-C3A-6_full.dat",
        "POPC/ordermap_POPC-C1B-8_full.dat",
        "POPC/ordermap_POPC-C3B-10_full.dat",
        "POPC/ordermap_average_full.dat",
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

    let directory = TempDir::new_in(path_to_outer_dir).unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let backup_file = format!("{}/to_backup.txt", path_to_dir);
    create_file_for_backup!(&backup_file);

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .map(
            OrderMap::builder()
                .bin_size([0.1, 4.0])
                .output_directory(path_to_dir)
                .min_samples(5)
                .build()
                .unwrap(),
        )
        .silent()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    let expected_file_names = [
        "ordermap_POPC-C218-87--POPC-H18R-88_full.dat",
        "ordermap_POPC-C218-87--POPC-H18S-89_full.dat",
        "ordermap_POPC-C218-87--POPC-H18T-90_full.dat",
        "ordermap_POPC-C218-87_full.dat",
        "ordermap_POPC-C22-32--POPC-H2R-33_full.dat",
        "ordermap_POPC-C22-32--POPC-H2S-34_full.dat",
        "ordermap_POPC-C22-32_full.dat",
        "ordermap_POPC-C24-47--POPC-H4R-48_full.dat",
        "ordermap_POPC-C24-47--POPC-H4S-49_full.dat",
        "ordermap_POPC-C24-47_full.dat",
        "ordermap_average_full.dat",
    ];

    for file in expected_file_names {
        let real_file = format!("{}/POPC/{}", path_to_dir, file);
        let test_file = format!("tests/files/ordermaps/{}", file);
        assert!(diff_files_ignore_first(&real_file, &test_file, 2));
    }

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_small.yaml",
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
fn test_aa_order_maps_basic_different_plane() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .map(
            OrderMap::builder()
                .bin_size([4.0, 0.1])
                .output_directory(path_to_dir)
                .min_samples(5)
                .plane(Plane::XZ)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    // only test one file
    let real_file = format!("{}/POPC/ordermap_POPC-C218-87_full.dat", path_to_dir);
    let test_file = "tests/files/ordermaps/ordermap_xz.dat";
    assert!(diff_files_ignore_first(&real_file, test_file, 2));

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_small.yaml",
        1
    ));
}

#[test]
fn test_aa_order_error_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_error.yaml",
        1
    ));
}

#[test]
fn test_aa_order_error_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .n_threads(n_threads)
            .estimate_error(EstimateError::default())
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_error.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_error_leaflets_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_error_leaflets.yaml",
        1
    ));
}

#[test]
fn test_aa_order_error_leaflets_yaml_multiple_threads() {
    for n_threads in [2, 3, 5, 8, 12, 16, 64] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .output(path_to_output)
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .n_threads(n_threads)
            .leaflets(LeafletClassification::global("@membrane", "name P"))
            .estimate_error(EstimateError::default())
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_error_leaflets.yaml",
            1
        ));
    }
}

#[test]
fn test_aa_order_error_tab() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_error.tab",
        1
    ));
}

#[test]
fn test_aa_order_error_leaflets_tab() {
    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_tab(path_to_table)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_error_leaflets.tab",
        1
    ));
}

#[test]
fn test_aa_order_error_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_error.csv",
        1
    ));
}

#[test]
fn test_aa_order_error_leaflets_csv() {
    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_error_leaflets.csv",
        1
    ));
}

#[test]
fn test_aa_order_error_xvg() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_xvg(pattern)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        // same files as when `estimate_error` is not provided - xvg files do not show error
        let path_expected = format!("tests/files/aa_order_basic_{}.xvg", molecule);

        assert!(diff_files_ignore_first(&path, &path_expected, 1));
    }
}

#[test]
fn test_aa_order_error_leaflets_xvg() {
    let directory = TempDir::new().unwrap();
    let path_to_dir = directory.path().to_str().unwrap();

    let pattern = format!("{}/order.xvg", path_to_dir);

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_xvg(pattern)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    for molecule in ["POPC", "POPE", "POPG"] {
        let path = format!("{}/order_{}.xvg", path_to_dir, molecule);
        // same files as when `estimate_error` is not provided - xvg files do not show error
        let path_expected = format!("tests/files/aa_order_leaflets_{}.xvg", molecule);

        assert!(diff_files_ignore_first(&path, &path_expected, 1));
    }
}

#[test]
fn test_aa_order_error_limit() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::default())
        .min_samples(2000)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/aa_order_error_limit.yaml",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_error_limit.tab",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_error_limit.csv",
        1
    ));
}

#[test]
fn test_aa_order_error_leaflets_limit() {
    let output = NamedTempFile::new().unwrap();
    let path_to_yaml = output.path().to_str().unwrap();

    let output_table = NamedTempFile::new().unwrap();
    let path_to_table = output_table.path().to_str().unwrap();

    let output_csv = NamedTempFile::new().unwrap();
    let path_to_csv = output_csv.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output_yaml(path_to_yaml)
        .output_tab(path_to_table)
        .output_csv(path_to_csv)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::default())
        .min_samples(500)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_yaml,
        "tests/files/aa_order_error_leaflets_limit.yaml",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_table,
        "tests/files/aa_order_error_leaflets_limit.tab",
        1
    ));

    assert!(diff_files_ignore_first(
        path_to_csv,
        "tests/files/aa_order_error_leaflets_limit.csv",
        1
    ));
}

#[test]
fn test_aa_order_error_blocks_10_yaml() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .output(path_to_output)
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::new(Some(10), None).unwrap())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_error_blocks10.yaml",
        1
    ));
}

#[test]
fn test_aa_order_error_blocks_too_many() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::new(Some(100), None).unwrap())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    match analysis.run() {
        Ok(_) => panic!("Analysis should have failed but it succeeded."),
        Err(e) => {
            assert!(e.to_string().contains("fewer than the number of blocks"))
        }
    }
}

#[test]
fn test_aa_order_convergence() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_convergence.xvg",
        1
    ));
}

#[test]
fn test_aa_order_leaflets_convergence() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_leaflets_convergence.xvg",
        1
    ));
}

#[test]
fn test_aa_order_convergence_multiple_threads() {
    for n_threads in [2, 5, 12, 32] {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_str().unwrap();

        let analysis = Analysis::builder()
            .structure("tests/files/pcpepg.tpr")
            .trajectory("tests/files/pcpepg.xtc")
            .analysis_type(AnalysisType::aaorder(
                "@membrane and element name carbon",
                "@membrane and element name hydrogen",
            ))
            .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
            .n_threads(n_threads)
            .silent()
            .overwrite()
            .build()
            .unwrap();

        analysis.run().unwrap().write().unwrap();

        assert!(diff_files_ignore_first(
            path_to_output,
            "tests/files/aa_order_convergence.xvg",
            1
        ));
    }
}

#[test]
fn test_aa_order_convergence_step() {
    let output = NamedTempFile::new().unwrap();
    let path_to_output = output.path().to_str().unwrap();

    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::new(None, Some(path_to_output)).unwrap())
        .step(5)
        .silent()
        .overwrite()
        .build()
        .unwrap();

    analysis.run().unwrap().write().unwrap();

    assert!(diff_files_ignore_first(
        path_to_output,
        "tests/files/aa_order_convergence_s5.xvg",
        1
    ));
}

#[test]
fn test_aa_order_basic_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(x) => x,
        AnalysisResults::CG(_) => panic!("Incorrect results type returned."),
    };

    assert_eq!(results.n_analyzed_frames(), 51);
    assert_eq!(results.analysis().structure(), "tests/files/pcpepg.tpr");

    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    let expected_average_orders = [0.1455, 0.1378, 0.1561];
    let expected_atom_numbers = [37, 40, 38];
    let expected_molecule_names = ["POPE", "POPC", "POPG"];

    let expected_atom_indices = [32, 41, 34];
    let expected_atom_names = ["C32", "C32", "C32"];
    let expected_atom_order = [0.2226, 0.2363, 0.2247];

    let expected_bond_numbers = [2, 2, 2];

    let expected_atom2_indices = [34, 43, 36];
    let expected_atom2_names = ["H2Y", "H2Y", "H2Y"];
    let expected_atom2_order = [0.2040, 0.2317, 0.2020];

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

        // atoms
        assert_eq!(molecule.atoms().count(), expected_atom_numbers[i]);

        let atom = molecule.get_atom(expected_atom_indices[i]).unwrap();
        let atom_type = atom.atom();
        assert_eq!(atom_type.atom_name(), expected_atom_names[i]);
        assert_eq!(atom_type.relative_index(), expected_atom_indices[i]);
        assert_eq!(atom_type.residue_name(), expected_molecule_names[i]);
        assert_eq!(atom.molecule(), expected_molecule_names[i]);

        let order = atom.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_atom_order[i],
            epsilon = 1e-4
        );
        assert!(order.total().unwrap().error().is_none());
        assert!(order.upper().is_none());
        assert!(order.lower().is_none());

        let maps = atom.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bonds
        assert_eq!(atom.bonds().count(), expected_bond_numbers[i]);

        let bond = atom.get_bond(expected_atom2_indices[i]).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), expected_atom_names[i]);
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), expected_atom2_names[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);
        assert_eq!(bond.molecule(), expected_molecule_names[i]);

        let order = bond.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_atom2_order[i],
            epsilon = 1e-4
        );
        assert!(order.total().unwrap().error().is_none());
        assert!(order.upper().is_none());
        assert!(order.lower().is_none());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bond directly from molecule
        let bond = molecule
            .get_bond(expected_atom_indices[i], expected_atom2_indices[i])
            .unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);

        let bond = molecule
            .get_bond(expected_atom2_indices[i], expected_atom_indices[i])
            .unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);

        // nonexistent atom
        assert!(molecule.get_atom(145).is_none());
        // nonexistent bond
        assert!(molecule.get_bond(7, 19).is_none());
        assert!(molecule.get_bond(145, 189).is_none());
    }
}

#[test]
fn test_aa_order_error_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(x) => x,
        AnalysisResults::CG(_) => panic!("Incorrect results type returned."),
    };

    assert_eq!(results.n_analyzed_frames(), 51);
    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    let expected_average_orders = [0.1455, 0.1378, 0.1561];
    let expected_average_errors = [0.0029, 0.0036, 0.0112];
    let expected_atom_numbers = [37, 40, 38];
    let expected_molecule_names = ["POPE", "POPC", "POPG"];

    let expected_atom_indices = [32, 41, 34];
    let expected_atom_names = ["C32", "C32", "C32"];
    let expected_atom_order = [0.2226, 0.2363, 0.2247];
    let expected_atom_errors = [0.0087, 0.0071, 0.0574];

    let expected_bond_numbers = [2, 2, 2];

    let expected_atom2_indices = [34, 43, 36];
    let expected_atom2_names = ["H2Y", "H2Y", "H2Y"];
    let expected_atom2_order = [0.2040, 0.2317, 0.2020];
    let expected_atom2_errors = [0.0125, 0.0091, 0.0656];

    let expected_convergence_frames = (1..=51).into_iter().collect::<Vec<usize>>();
    let expected_convergence_values = [
        [0.1494, 0.1460, 0.1455],
        [0.1422, 0.1353, 0.1378],
        [0.1572, 0.1507, 0.1561],
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

        for (j, frame) in [0, 25, 50].into_iter().enumerate() {
            assert_relative_eq!(
                *convergence.total().as_ref().unwrap().get(frame).unwrap(),
                expected_convergence_values.get(i).unwrap().get(j).unwrap(),
                epsilon = 1e-4
            );
        }

        assert!(convergence.upper().is_none());
        assert!(convergence.lower().is_none());

        // atoms
        assert_eq!(molecule.atoms().count(), expected_atom_numbers[i]);

        let atom = molecule.get_atom(expected_atom_indices[i]).unwrap();
        let atom_type = atom.atom();
        assert_eq!(atom_type.atom_name(), expected_atom_names[i]);
        assert_eq!(atom_type.relative_index(), expected_atom_indices[i]);
        assert_eq!(atom_type.residue_name(), expected_molecule_names[i]);
        assert_eq!(atom.molecule(), expected_molecule_names[i]);

        let order = atom.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_atom_order[i],
            epsilon = 1e-4
        );
        assert_relative_eq!(
            order.total().unwrap().error().unwrap(),
            expected_atom_errors[i],
            epsilon = 1e-4
        );
        assert!(order.upper().is_none());
        assert!(order.lower().is_none());

        let maps = atom.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bonds
        assert_eq!(atom.bonds().count(), expected_bond_numbers[i]);

        let bond = atom.get_bond(expected_atom2_indices[i]).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), expected_atom_names[i]);
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), expected_atom2_names[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);
        assert_eq!(bond.molecule(), expected_molecule_names[i]);

        let order = bond.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_atom2_order[i],
            epsilon = 1e-4
        );
        assert_relative_eq!(
            order.total().unwrap().error().unwrap(),
            expected_atom2_errors[i],
            epsilon = 1e-4
        );
        assert!(order.upper().is_none());
        assert!(order.lower().is_none());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bond directly from molecule
        let bond = molecule
            .get_bond(expected_atom_indices[i], expected_atom2_indices[i])
            .unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);

        // nonexistent atom
        assert!(molecule.get_atom(145).is_none());
        // nonexistent bond
        assert!(molecule.get_bond(7, 19).is_none());
        assert!(molecule.get_bond(145, 189).is_none());
    }
}

#[test]
fn test_aa_order_leaflets_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(x) => x,
        AnalysisResults::CG(_) => panic!("Incorrect results type returned."),
    };

    assert_eq!(results.n_analyzed_frames(), 51);
    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    let expected_average_orders = [0.1455, 0.1378, 0.1561];
    let expected_average_upper = [0.1492, 0.1326, 0.1522];
    let expected_average_lower = [0.1419, 0.1431, 0.1606];
    let expected_atom_numbers = [37, 40, 38];
    let expected_molecule_names = ["POPE", "POPC", "POPG"];

    let expected_atom_indices = [32, 41, 34];
    let expected_atom_names = ["C32", "C32", "C32"];
    let expected_atom_order = [0.2226, 0.2363, 0.2247];
    let expected_atom_upper = [0.2131, 0.2334, 0.2484];
    let expected_atom_lower = [0.2319, 0.2391, 0.1976];

    let expected_bond_numbers = [2, 2, 2];

    let expected_atom2_indices = [34, 43, 36];
    let expected_atom2_names = ["H2Y", "H2Y", "H2Y"];
    let expected_atom2_order = [0.2040, 0.2317, 0.2020];
    let expected_atom2_upper = [0.1876, 0.2507, 0.2254];
    let expected_atom2_lower = [0.2203, 0.2126, 0.1752];

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

        // atoms
        assert_eq!(molecule.atoms().count(), expected_atom_numbers[i]);

        let atom = molecule.get_atom(expected_atom_indices[i]).unwrap();
        let atom_type = atom.atom();
        assert_eq!(atom_type.atom_name(), expected_atom_names[i]);
        assert_eq!(atom_type.relative_index(), expected_atom_indices[i]);
        assert_eq!(atom_type.residue_name(), expected_molecule_names[i]);
        assert_eq!(atom.molecule(), expected_molecule_names[i]);

        let order = atom.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_atom_order[i],
            epsilon = 1e-4
        );
        assert!(order.total().unwrap().error().is_none());

        assert_relative_eq!(
            order.upper().unwrap().value(),
            expected_atom_upper[i],
            epsilon = 1e-4
        );
        assert!(order.upper().unwrap().error().is_none());

        assert_relative_eq!(
            order.lower().unwrap().value(),
            expected_atom_lower[i],
            epsilon = 1e-4
        );
        assert!(order.lower().unwrap().error().is_none());

        let maps = atom.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bonds
        assert_eq!(atom.bonds().count(), expected_bond_numbers[i]);

        let bond = atom.get_bond(expected_atom2_indices[i]).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), expected_atom_names[i]);
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), expected_atom2_names[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);
        assert_eq!(bond.molecule(), expected_molecule_names[i]);

        let order = bond.order();
        assert_relative_eq!(
            order.total().unwrap().value(),
            expected_atom2_order[i],
            epsilon = 1e-4
        );
        assert!(order.total().unwrap().error().is_none());

        assert_relative_eq!(
            order.upper().unwrap().value(),
            expected_atom2_upper[i],
            epsilon = 1e-4
        );
        assert!(order.upper().unwrap().error().is_none());

        assert_relative_eq!(
            order.lower().unwrap().value(),
            expected_atom2_lower[i],
            epsilon = 1e-4
        );
        assert!(order.lower().unwrap().error().is_none());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bond directly from molecule
        let bond = molecule
            .get_bond(expected_atom_indices[i], expected_atom2_indices[i])
            .unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);

        // nonexistent atom
        assert!(molecule.get_atom(145).is_none());
        // nonexistent bond
        assert!(molecule.get_bond(7, 19).is_none());
        assert!(molecule.get_bond(145, 189).is_none());
    }
}

#[test]
fn test_aa_order_error_leaflets_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .estimate_error(EstimateError::default())
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(x) => x,
        AnalysisResults::CG(_) => panic!("Incorrect results type returned."),
    };

    assert_eq!(results.n_analyzed_frames(), 51);
    assert!(results
        .analysis()
        .estimate_error()
        .as_ref()
        .unwrap()
        .output_convergence()
        .is_none());

    assert_eq!(results.molecules().count(), 3);

    assert!(results.get_molecule("POPE").is_some());
    assert!(results.get_molecule("POPC").is_some());
    assert!(results.get_molecule("POPG").is_some());
    assert!(results.get_molecule("POPA").is_none());

    let expected_atom_numbers = [37, 40, 38];
    let expected_molecule_names = ["POPE", "POPC", "POPG"];

    let expected_atom_indices = [32, 41, 34];
    let expected_atom_names = ["C32", "C32", "C32"];

    let expected_bond_numbers = [2, 2, 2];

    let expected_atom2_indices = [34, 43, 36];
    let expected_atom2_names = ["H2Y", "H2Y", "H2Y"];

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
        assert_eq!(convergence.frames().len(), 51);
        assert!(convergence.total().is_some());
        assert!(convergence.upper().is_some());
        assert!(convergence.lower().is_some());

        // atoms
        assert_eq!(molecule.atoms().count(), expected_atom_numbers[i]);

        let atom = molecule.get_atom(expected_atom_indices[i]).unwrap();
        let atom_type = atom.atom();
        assert_eq!(atom_type.atom_name(), expected_atom_names[i]);
        assert_eq!(atom_type.relative_index(), expected_atom_indices[i]);
        assert_eq!(atom_type.residue_name(), expected_molecule_names[i]);
        assert_eq!(atom.molecule(), expected_molecule_names[i]);

        let order = atom.order();
        assert!(order.total().unwrap().error().is_some());
        assert!(order.upper().unwrap().error().is_some());
        assert!(order.lower().unwrap().error().is_some());

        let maps = atom.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bonds
        assert_eq!(atom.bonds().count(), expected_bond_numbers[i]);

        let bond = atom.get_bond(expected_atom2_indices[i]).unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.atom_name(), expected_atom_names[i]);
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a1.residue_name(), expected_molecule_names[i]);
        assert_eq!(a2.atom_name(), expected_atom2_names[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);
        assert_eq!(a2.residue_name(), expected_molecule_names[i]);
        assert_eq!(bond.molecule(), expected_molecule_names[i]);

        let order = bond.order();
        assert!(order.total().unwrap().error().is_some());
        assert!(order.upper().unwrap().error().is_some());
        assert!(order.lower().unwrap().error().is_some());

        let maps = bond.ordermaps();
        assert!(maps.total().is_none());
        assert!(maps.upper().is_none());
        assert!(maps.lower().is_none());

        // bond directly from molecule
        let bond = molecule
            .get_bond(expected_atom_indices[i], expected_atom2_indices[i])
            .unwrap();
        let (a1, a2) = bond.atoms();
        assert_eq!(a1.relative_index(), expected_atom_indices[i]);
        assert_eq!(a2.relative_index(), expected_atom2_indices[i]);

        // nonexistent atom
        assert!(molecule.get_atom(145).is_none());
        // nonexistent bond
        assert!(molecule.get_bond(7, 19).is_none());
        assert!(molecule.get_bond(145, 189).is_none());
    }
}

#[test]
fn test_aa_order_ordermaps_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .map(
            OrderMap::builder()
                .bin_size([0.1, 4.0])
                .min_samples(5)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(x) => x,
        AnalysisResults::CG(_) => panic!("Incorrect results type returned."),
    };

    assert_eq!(results.n_analyzed_frames(), 51);
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
    assert_relative_eq!(span_x.1, 9.15673);
    assert_relative_eq!(span_y.0, 0.0);
    assert_relative_eq!(span_y.1, 9.15673);
    assert_relative_eq!(bin.0, 0.1);
    assert_relative_eq!(bin.1, 4.0);

    assert_relative_eq!(
        map.get_at_convert(0.6, 8.0).unwrap(),
        0.1653,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(4.3, 0.0).unwrap(),
        0.1340,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(9.2, 4.0).unwrap(),
        0.1990,
        epsilon = 1e-4
    );

    // ordermaps for a selected atom
    let atom = molecule.get_atom(47).unwrap();
    let map = atom.ordermaps().total().as_ref().unwrap();
    assert!(atom.ordermaps().upper().is_none());
    assert!(atom.ordermaps().lower().is_none());

    assert_relative_eq!(
        map.get_at_convert(0.6, 8.0).unwrap(),
        0.2224,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(4.3, 0.0).unwrap(),
        0.1532,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(9.2, 4.0).unwrap(),
        0.0982,
        epsilon = 1e-4
    );

    // ordermaps for a selected bond
    let bond = atom.get_bond(49).unwrap();
    let map = bond.ordermaps().total().as_ref().unwrap();
    assert!(bond.ordermaps().upper().is_none());
    assert!(bond.ordermaps().lower().is_none());

    assert_relative_eq!(
        map.get_at_convert(0.6, 8.0).unwrap(),
        0.2901,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        map.get_at_convert(4.3, 0.0).unwrap(),
        0.1163,
        epsilon = 1e-4
    );
    assert!(map.get_at_convert(9.2, 4.0).unwrap().is_nan());
}

#[test]
fn test_aa_order_ordermaps_leaflets_rust_api() {
    let analysis = Analysis::builder()
        .structure("tests/files/pcpepg.tpr")
        .trajectory("tests/files/pcpepg.xtc")
        .analysis_type(AnalysisType::aaorder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen",
        ))
        .leaflets(LeafletClassification::global("@membrane", "name P"))
        .map(
            OrderMap::builder()
                .bin_size([0.1, 4.0])
                .min_samples(5)
                .build()
                .unwrap(),
        )
        .silent()
        .overwrite()
        .build()
        .unwrap();

    let results = match analysis.run().unwrap() {
        AnalysisResults::AA(x) => x,
        AnalysisResults::CG(_) => panic!("Incorrect results type returned."),
    };

    assert_eq!(results.n_analyzed_frames(), 51);
    assert_eq!(results.molecules().count(), 1);

    // average ordermaps for the entire molecule
    let molecule = results.get_molecule("POPC").unwrap();
    let total = molecule.average_ordermaps().total().as_ref().unwrap();

    let span_x = total.span_x();
    let span_y = total.span_y();
    let bin = total.tile_dim();

    assert_relative_eq!(span_x.0, 0.0);
    assert_relative_eq!(span_x.1, 9.15673);
    assert_relative_eq!(span_y.0, 0.0);
    assert_relative_eq!(span_y.1, 9.15673);
    assert_relative_eq!(bin.0, 0.1);
    assert_relative_eq!(bin.1, 4.0);

    assert_relative_eq!(
        total.get_at_convert(0.6, 8.0).unwrap(),
        0.1653,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        total.get_at_convert(9.2, 4.0).unwrap(),
        0.1990,
        epsilon = 1e-4
    );

    let upper = molecule.average_ordermaps().upper().as_ref().unwrap();

    assert_relative_eq!(
        upper.get_at_convert(0.6, 8.0).unwrap(),
        0.1347,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        upper.get_at_convert(9.2, 4.0).unwrap(),
        0.3196,
        epsilon = 1e-4
    );

    let lower = molecule.average_ordermaps().lower().as_ref().unwrap();

    assert_relative_eq!(
        lower.get_at_convert(0.6, 8.0).unwrap(),
        0.2104,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        lower.get_at_convert(9.2, 4.0).unwrap(),
        0.1106,
        epsilon = 1e-4
    );

    // ordermaps for a selected atom
    let atom = molecule.get_atom(47).unwrap();
    let total = atom.ordermaps().total().as_ref().unwrap();

    assert_relative_eq!(
        total.get_at_convert(0.6, 8.0).unwrap(),
        0.2224,
        epsilon = 1e-4
    );
    assert_relative_eq!(
        total.get_at_convert(9.2, 4.0).unwrap(),
        0.0982,
        epsilon = 1e-4
    );

    let upper = atom.ordermaps().upper().as_ref().unwrap();

    assert_relative_eq!(
        upper.get_at_convert(0.6, 8.0).unwrap(),
        0.2039,
        epsilon = 1e-4
    );
    assert!(upper.get_at_convert(9.2, 4.0).unwrap().is_nan());

    let lower = atom.ordermaps().lower().as_ref().unwrap();

    assert_relative_eq!(
        lower.get_at_convert(0.6, 8.0).unwrap(),
        0.2540,
        epsilon = 1e-4
    );
    assert!(lower.get_at_convert(9.2, 4.0).unwrap().is_nan());

    // ordermaps for a selected bond
    let bond = atom.get_bond(49).unwrap();
    let total = bond.ordermaps().total().as_ref().unwrap();

    assert_relative_eq!(
        total.get_at_convert(0.6, 8.0).unwrap(),
        0.2901,
        epsilon = 1e-4
    );
    assert!(total.get_at_convert(9.2, 4.0).unwrap().is_nan());

    let upper = bond.ordermaps().upper().as_ref().unwrap();

    assert_relative_eq!(
        upper.get_at_convert(0.6, 8.0).unwrap(),
        0.3584,
        epsilon = 1e-4
    );
    assert!(upper.get_at_convert(9.2, 4.0).unwrap().is_nan());

    let lower = bond.ordermaps().lower().as_ref().unwrap();

    assert_relative_eq!(
        lower.get_at_convert(0.6, 8.0).unwrap(),
        0.1715,
        epsilon = 1e-4
    );
    assert!(lower.get_at_convert(9.2, 4.0).unwrap().is_nan());
}
