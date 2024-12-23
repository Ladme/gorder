// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the coarse-grained order parameters.

use groan_rs::{
    files::FileType,
    prelude::{ProgressPrinter, XtcReader},
    system::System,
};

use crate::{
    analysis::{
        common::{analyze_frame, macros::group_name, sanity_check_molecules, write_results},
        topology::SystemTopology,
    },
    errors::AnalysisError,
    input::Analysis,
    presentation::{cgpresenter::CGOrderResults, CGOrder},
    PANIC_MESSAGE,
};

/// Analyze the coarse-grained order parameters.
pub(super) fn analyze_coarse_grained(
    analysis: &Analysis,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut system = System::from_file_with_format(analysis.structure(), FileType::TPR)?;
    log::info!("Read molecular topology from '{}'.", analysis.structure());

    if let Some(ndx) = analysis.index() {
        system.read_ndx(ndx)?;
        log::info!(
            "Read {} group(s) from ndx file '{}'.",
            system.get_n_groups() - 2,
            ndx
        );
    }

    super::common::create_group(
        &mut system,
        "Beads",
        analysis.beads().as_ref().unwrap_or_else(||
             panic!("FATAL GORDER ERROR | cgorder::analyze_coarse_grained | Selection of order beads should be provided. {}", PANIC_MESSAGE)),
    )?;

    log::info!(
        "Detected {} beads for order calculation using a query '{}'.",
        system
            .group_get_n_atoms(group_name!("Beads"))
            .expect(PANIC_MESSAGE),
        analysis.beads().as_ref().expect(PANIC_MESSAGE)
    );

    // prepare system for leaflet classification
    if let Some(leaflet) = analysis.leaflets() {
        leaflet.prepare_system(&mut system)?;
    }

    log::info!("Detecting molecule types...");
    log::logger().flush();

    // get the relevant molecules
    let molecules = super::common::classify_molecules(
        &system,
        "Beads",
        "Beads",
        analysis.leaflets().as_ref(),
        analysis.membrane_normal().into(),
        analysis.map().as_ref(),
        analysis.min_samples(),
        analysis.estimate_error().as_ref(),
    )?;

    if !sanity_check_molecules(&molecules) {
        return Ok(());
    }

    let data = SystemTopology::new(
        molecules,
        analysis.membrane_normal().into(),
        analysis.estimate_error().clone(),
    );
    data.info();

    let progress_printer = if analysis.silent() {
        None
    } else {
        Some(ProgressPrinter::new().with_print_freq(100 / analysis.n_threads()))
    };

    log::info!(
        "Will read trajectory file '{}' (start: {} ps, end: {} ps, step: {}).",
        analysis.trajectory(),
        analysis.begin(),
        analysis.end(),
        analysis.step()
    );

    log::info!(
        "Performing the analysis using {} thread(s)...",
        analysis.n_threads()
    );

    // run the analysis in parallel
    let result = system.traj_iter_map_reduce::<XtcReader, SystemTopology, AnalysisError>(
        analysis.trajectory(),
        analysis.n_threads(),
        analyze_frame,
        data,
        Some(analysis.begin()),
        Some(analysis.end()),
        Some(analysis.step()),
        progress_printer,
    )?;

    // write out the maps
    result.handle_ordermap_directory(analysis.overwrite())?;
    result.prepare_directories()?;
    result.write_ordermaps_bonds::<CGOrder>()?;

    // write out the results
    let results = CGOrderResults::from(result);
    write_results(results, analysis)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use groan_rs::prelude::Dimension;

    use crate::{
        analysis::molecule::{BondType, MoleculeType},
        input::LeafletClassification,
    };

    use super::*;

    fn prepare_data_for_tests(
        leaflet_classification: Option<LeafletClassification>,
    ) -> (System, SystemTopology) {
        let mut system = System::from_file("tests/files/cg.tpr").unwrap();

        system
            .group_create(group_name!("Beads"), "@membrane")
            .unwrap();

        if let Some(leaflet) = &leaflet_classification {
            leaflet.prepare_system(&mut system).unwrap();
        }

        let molecules = super::super::common::classify_molecules(
            &system,
            "Beads",
            "Beads",
            leaflet_classification.as_ref(),
            Dimension::Z,
            None,
            1,
            None,
        )
        .unwrap();

        (system, SystemTopology::new(molecules, Dimension::Z, None))
    }

    fn expected_total_orders() -> [Vec<f32>; 3] {
        [
            vec![
                -39.98374, 146.05228, -43.874344, 99.498955, 121.651955, 95.177536, 87.30928,
                68.875565, 122.233536, 99.90945, 71.819855,
            ],
            vec![
                -20.425762, 150.84778, -43.26054, 93.48221, 117.66597, 89.95231, 77.81377,
                55.82394, 122.598946, 104.67025, 63.30497,
            ],
            vec![
                -1.502139, 16.002146, -6.857792, 12.772206, 12.551049, 8.470685, 10.038215,
                8.334017, 9.102596, 9.754787, 6.800507,
            ],
        ]
    }

    fn expected_upper_orders() -> [Vec<f32>; 3] {
        [
            vec![
                -21.27362, 71.48379, -19.56165, 52.70994, 63.685772, 51.68486, 38.602947,
                38.312897, 69.63991, 53.493954, 37.61219,
            ],
            vec![
                -11.361008, 72.02412, -18.039215, 46.06996, 59.89987, 41.5592, 47.412964,
                30.745993, 58.742195, 53.235092, 32.349113,
            ],
            vec![
                -0.483492, 7.674479, -3.941579, 6.569143, 7.597515, 4.912679, 5.621792, 2.704198,
                4.938721, 4.097075, 3.658832,
            ],
        ]
    }

    fn expected_lower_orders() -> [Vec<f32>; 3] {
        [
            vec![
                -18.71012, 74.56849, -24.312695, 46.789017, 57.966187, 43.492676, 48.706337,
                30.56267, 52.59362, 46.41549, 34.207664,
            ],
            vec![
                -9.064755, 78.82366, -25.221325, 47.41225, 57.766094, 48.393112, 30.400806,
                25.07795, 63.85675, 51.43516, 30.955854,
            ],
            vec![
                -1.018647, 8.327666, -2.916213, 6.203063, 4.953534, 3.558006, 4.416423, 5.629819,
                4.163875, 5.657712, 3.141675,
            ],
        ]
    }

    #[test]
    fn test_cgorder_analyze_frame_basic() {
        let (system, mut data) = prepare_data_for_tests(None);

        analyze_frame(&system, &mut data).unwrap();
        let expected_total_orders = expected_total_orders();

        for (m, molecule) in data.molecule_types().iter().enumerate() {
            let n_instances = molecule.order_bonds().bond_types()[0].bonds().len();
            let orders = molecule
                .order_bonds()
                .bond_types()
                .iter()
                .map(|b| b.total().order().into())
                .collect::<Vec<f32>>();

            let samples = molecule
                .order_bonds()
                .bond_types()
                .iter()
                .map(|b| b.total().n_samples())
                .collect::<Vec<usize>>();

            assert_eq!(orders.len(), expected_total_orders[m].len());
            for (real, expected) in orders.iter().zip(expected_total_orders[m].iter()) {
                assert_relative_eq!(real, expected);
            }

            for sample in samples {
                assert_eq!(sample, n_instances);
            }
        }
    }

    fn collect_bond_data<T, F>(molecule: &MoleculeType, func: F) -> Vec<T>
    where
        F: Fn(&BondType) -> T,
    {
        molecule
            .order_bonds()
            .bond_types()
            .iter()
            .map(func)
            .collect()
    }

    #[test]
    fn test_cgorder_analyze_frame_leaflets() {
        let classifier = LeafletClassification::global("@membrane", "name PO4");

        let (system, mut data) = prepare_data_for_tests(Some(classifier));

        analyze_frame(&system, &mut data).unwrap();
        let expected_total_orders = expected_total_orders();
        let expected_upper_orders = expected_upper_orders();
        let expected_lower_orders = expected_lower_orders();
        let expected_total_samples = [242, 242, 24];
        let expected_upper_samples = [121, 121, 12];
        let expected_lower_samples = [121, 121, 12];

        for (m, molecule) in data.molecule_types().iter().enumerate() {
            let total_orders: Vec<f32> = collect_bond_data(molecule, |b| b.total().order().into());
            let upper_orders =
                collect_bond_data(molecule, |b| b.upper().as_ref().unwrap().order().into());
            let lower_orders =
                collect_bond_data(molecule, |b| b.lower().as_ref().unwrap().order().into());
            let total_samples = collect_bond_data(molecule, |b| b.total().n_samples());
            let upper_samples =
                collect_bond_data(molecule, |b| b.upper().as_ref().unwrap().n_samples());
            let lower_samples =
                collect_bond_data(molecule, |b| b.lower().as_ref().unwrap().n_samples());

            for (order, samples, expected_order, expected_samples) in [
                (
                    &total_orders,
                    &total_samples,
                    &expected_total_orders[m],
                    expected_total_samples[m],
                ),
                (
                    &upper_orders,
                    &upper_samples,
                    &expected_upper_orders[m],
                    expected_upper_samples[m],
                ),
                (
                    &lower_orders,
                    &lower_samples,
                    &expected_lower_orders[m],
                    expected_lower_samples[m],
                ),
            ] {
                assert_eq!(order.len(), expected_order.len());
                for (real, expected) in order.iter().zip(expected_order.iter()) {
                    assert_relative_eq!(real, expected);
                }

                for &sample in samples {
                    assert_eq!(sample, expected_samples);
                }
            }
        }
    }
}
