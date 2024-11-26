// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the atomistic lipid order parameters.

use super::leaflets::MoleculeLeafletClassification;
use super::{auxiliary::macros::group_name, topology::SystemTopology};
use crate::errors::{AnalysisError, TopologyError};
use crate::presentation::aapresenter::AAOrderResults;
use crate::{Analysis, PANIC_MESSAGE};

use groan_rs::prelude::OrderedAtomIterator;
use groan_rs::{
    files::FileType,
    prelude::{ProgressPrinter, XtcReader},
    system::System,
};

/// Calculate the atomistic lipid order parameters.
pub(super) fn analyze_atomistic(
    analysis: &Analysis,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut system = System::from_file_with_format(analysis.structure(), FileType::TPR)?;
    log::info!("Read molecular topology from '{}'.", analysis.structure());

    super::auxiliary::create_group(
        &mut system,
        "HeavyAtoms",
        analysis.heavy_atoms().as_ref().unwrap_or_else(||
             panic!("FATAL GORDER ERROR | aaorder::analyze_atomistic | Selection of heavy atoms should be provided. {}", PANIC_MESSAGE)),
    )?;

    log::info!(
        "Detected {} heavy atoms using query '{}'.",
        system
            .group_get_n_atoms(group_name!("HeavyAtoms"))
            .expect(PANIC_MESSAGE),
        analysis.heavy_atoms().as_ref().expect(PANIC_MESSAGE)
    );

    super::auxiliary::create_group(
        &mut system,
        "Hydrogens", 
        analysis.hydrogens().as_ref().unwrap_or_else(||
            panic!("FATAL GORDER ERROR | aaorder::analyze_atomistic | Selection of hydrogens should be provided. {}", PANIC_MESSAGE))
    )?;

    log::info!(
        "Detected {} hydrogen atoms using query '{}'.",
        system
            .group_get_n_atoms(group_name!("Hydrogens"))
            .expect(PANIC_MESSAGE),
        analysis.hydrogens().as_ref().expect(PANIC_MESSAGE)
    );

    // check that heavy_atoms and hydrogens do not overlap
    let n_overlapping = system
        .group_iter(group_name!("HeavyAtoms"))
        .expect(PANIC_MESSAGE)
        .intersection(
            system
                .group_iter(group_name!("Hydrogens"))
                .expect(PANIC_MESSAGE),
        )
        .count();
    if n_overlapping > 0 {
        return Err(Box::from(TopologyError::AtomsOverlap {
            n_overlapping,
            name1: "HeavyAtoms".to_owned(),
            query1: analysis.heavy_atoms().expect(PANIC_MESSAGE).clone(),
            name2: "Hydrogens".to_owned(),
            query2: analysis.hydrogens().expect(PANIC_MESSAGE).clone(),
        }));
    }
    // prepare system for leaflet classification
    if let Some(leaflet) = analysis.leaflets() {
        leaflet.prepare_system(&mut system)?;
    }

    log::info!("Detecting molecule types...");
    log::logger().flush();

    // get the relevant molecules
    let molecules = super::auxiliary::classify_molecules(
        &system,
        "HeavyAtoms",
        "Hydrogens",
        analysis.leaflets().as_ref(),
        analysis.membrane_normal().into(),
        analysis.map().as_ref(),
        analysis.min_samples(),
    )?;

    // if no molecules are detected, end the analysis
    if molecules.len() == 0 {
        log::warn!("No molecules suitable for the analysis detected.");
        return Ok(());
    }

    // if only empty molecules are detected, end the analysis
    for mol in molecules.iter() {
        if mol.order_bonds().bond_types().is_empty() {
            log::warn!("No bonds suitable for the analysis detected.");
            return Ok(());
        }
    }

    let data = SystemTopology::new(molecules, analysis.membrane_normal().into());
    data.info();

    let progress_printer = if analysis.silent() {
        None
    } else {
        Some(ProgressPrinter::new())
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
    result.write_ordermaps_bonds()?;
    result.write_ordermaps_atoms()?;

    // write out the results
    log::info!(
        "Writing the order parameters into a yaml file '{}'...",
        analysis.output()
    );
    log::logger().flush();

    let results = AAOrderResults::from(result);
    results.write_yaml(
        analysis.output(),
        analysis.structure(),
        analysis.trajectory(),
        analysis.overwrite(),
    )?;

    if let Some(tab) = analysis.output_tab() {
        log::info!("Writing the order parameters into a table '{}'...", tab);

        log::logger().flush();
        results.write_tab(
            tab,
            analysis.structure(),
            analysis.trajectory(),
            analysis.overwrite(),
        )?;
    };

    if let Some(xvg) = analysis.output_xvg() {
        log::info!("Writing the order parameters into xvg file(s)...");

        log::logger().flush();
        results.write_xvg(
            xvg,
            analysis.structure(),
            analysis.trajectory(),
            analysis.overwrite(),
        )?;
    }

    if let Some(csv) = analysis.output_csv() {
        log::info!("Writing the order parameters into a csv file '{}'...", csv);
        log::logger().flush();
        results.write_csv(csv, analysis.overwrite())?;
    }

    Ok(())
}

/// Analyze order parameters in a single simulation frame.
fn analyze_frame(frame: &System, data: &mut SystemTopology) -> Result<(), AnalysisError> {
    let simbox = frame.get_box().ok_or(AnalysisError::UndefinedBox)?;

    if !simbox.is_orthogonal() {
        return Err(AnalysisError::NotOrthogonalBox);
    }

    if simbox.is_zero() {
        return Err(AnalysisError::ZeroBox);
    }

    let membrane_normal = data.membrane_normal().into();
    let membrane_center = {
        match data
            .molecule_types()
            .get(0)
            .expect(PANIC_MESSAGE)
            .leaflet_classification()
        {
            Some(MoleculeLeafletClassification::Global(_, _)) => frame
                .group_get_center(group_name!("Membrane"))
                .map(Some)
                .map_err(|_| AnalysisError::InvalidGlobalMembraneCenter),
            _ => Ok(None), // no applicable classification, return Ok with None
        }
    }?;

    // assign molecules to leaflets
    for molecule in data.molecule_types_mut().iter_mut() {
        if let Some(classifier) = molecule.leaflet_classification_mut() {
            match classifier {
                MoleculeLeafletClassification::Global(x, _) => {
                    x.set_membrane_center(membrane_center.clone().expect(PANIC_MESSAGE));
                }
                MoleculeLeafletClassification::Local(x, _) => {
                    x.set_membrane_center(frame, membrane_normal)?;
                }
                MoleculeLeafletClassification::Individual(_, _) => (),
            };

            classifier.assign_lipids(frame)?;
        }

        molecule.analyze_frame(frame, &simbox, &membrane_normal.into())?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use approx::assert_relative_eq;
    use groan_rs::prelude::Dimension;

    use crate::{
        analysis::molecule::{BondType, MoleculeType},
        LeafletClassification,
    };

    use super::*;

    fn prepare_data_for_tests(
        leaflet_classification: Option<LeafletClassification>,
    ) -> (System, SystemTopology) {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();

        system
            .group_create(
                group_name!("HeavyAtoms"),
                "@membrane and element name carbon",
            )
            .unwrap();

        system
            .group_create(
                group_name!("Hydrogens"),
                "@membrane and element name hydrogen",
            )
            .unwrap();

        if let Some(leaflet) = &leaflet_classification {
            leaflet.prepare_system(&mut system).unwrap();
        }

        let molecules = super::super::auxiliary::classify_molecules(
            &system,
            "HeavyAtoms",
            "Hydrogens",
            leaflet_classification.as_ref(),
            Dimension::Z,
            None,
            1,
        )
        .unwrap();

        (system, SystemTopology::new(molecules, Dimension::Z))
    }

    fn expected_total_orders() -> [Vec<f32>; 3] {
        [
            vec![
                11.389851, 5.8786674, -13.541704, -11.096256, 27.989347, 25.653984, 19.826765,
                12.905661, 10.441488, 22.069092, 1.5912378, 27.31636, 30.524326, 24.914057,
                31.363668, 28.719318, 30.966377, 28.59851, 30.240564, 22.511253, 31.064243,
                22.650856, 28.634373, 14.978862, 9.811056, 9.882705, 12.237981, 11.099975,
                11.654347, 19.075571, 18.815758, 23.57948, 19.751911, 18.58691, 23.117897,
                20.328485, 23.120123, 19.377031, 21.543617, 5.945527, 20.534142, 4.208125,
                6.063838, 2.4329069, 19.9338, 23.645094, 26.961494, 24.438534, 28.552168,
                25.008698, 27.270025, 36.33753, 28.786337, 36.512115, 34.37085, 31.867886,
                32.78751, 31.194407, 26.743479, 29.923029, 25.999723, 29.702309, 28.175156,
                25.97152, 24.937578, 23.082615, 25.887957, 11.512761, 19.91333, 9.057474, 5.056279,
                6.974998, 6.6677446,
            ],
            vec![
                7.6329927, 11.82006, -6.7167015, 0.9685844, -1.9114108, 1.6918764, -5.1857877,
                2.0685284, -6.9205585, -1.8942664, 5.227478, 3.6607904, -9.327742, 23.518724,
                29.831762, 19.633091, 13.166875, 16.535954, 20.407568, 4.831814, 34.152706,
                29.234451, 27.199028, 29.97675, 27.201118, 29.567423, 26.336637, 32.356857,
                26.892096, 25.80453, 25.827118, 27.205551, 16.994942, 19.094893, 11.559245,
                4.492818, 4.772769, 8.061023, 17.412926, 13.769697, 13.628067, 16.356262,
                18.402748, 15.024002, 22.206394, 14.288321, 17.652033, 13.681939, 20.933537,
                9.603974, 16.412033, 0.7050859, 2.477138, 29.031404, 21.529348, 28.327091,
                28.09423, 29.071299, 23.746437, 34.364956, 27.89767, 30.205153, 30.77279,
                32.497437, 24.186451, 31.099829, 27.633347, 32.111458, 18.132113, 30.367407,
                20.990025, 27.911137, 20.481098, 15.41196, 15.65695, 14.946924, 13.218306,
                9.822741, 14.016406, 5.663187, 0.9808493, 2.2279348,
            ],
            vec![
                3.2941704,
                -3.4283745,
                -1.0800903,
                -2.7592607,
                -0.78275836,
                3.254287,
                1.7753217,
                1.837405,
                -0.5119834,
                2.6444633,
                3.7265575,
                -3.0317092,
                3.2530794,
                2.4848883,
                -0.7066252,
                2.2119145,
                1.5368084,
                1.9389439,
                3.4290047,
                2.2273543,
                4.541318,
                -1.4675866,
                4.356757,
                0.622792,
                0.6956645,
                2.4216743,
                0.95012796,
                2.0549011,
                1.2461977,
                0.6637146,
                2.0619702,
                0.8555777,
                2.9841769,
                1.3860189,
                3.9738538,
                3.595663,
                4.880531,
                4.8376565,
                3.8616748,
                4.9901414,
                4.761696,
                2.181141,
                1.6104923,
                1.3128071,
                0.47169197,
                3.6841085,
                2.4182067,
                4.1024256,
                3.5949047,
                3.3663125,
                3.2254627,
                1.8537161,
                3.553464,
                2.3270063,
                2.9063377,
                3.51264,
                1.606179,
                2.8827128,
                1.5890352,
                3.9305735,
                1.550467,
                2.114952,
                3.5553129,
                4.157771,
                3.7806582,
                3.2784903,
                2.5110192,
                2.231767,
                1.3526614,
                -0.3523764,
                0.36293957,
                -0.38007507,
                -1.8708402,
                0.9648583,
            ],
        ]
    }

    fn expected_upper_orders() -> [Vec<f32>; 3] {
        [
            vec![
                5.291138,
                0.7893448,
                -7.4505186,
                -4.3143415,
                14.737594,
                9.5098,
                12.143466,
                10.080931,
                7.8712277,
                11.660126,
                -7.1850986,
                13.638062,
                13.475298,
                12.030589,
                14.92661,
                14.842073,
                13.062992,
                13.427263,
                9.725433,
                11.220982,
                12.190741,
                11.2161,
                13.531443,
                7.9799776,
                3.7025452,
                5.2039742,
                1.5994557,
                5.2090025,
                4.9148817,
                5.659866,
                10.39819,
                8.781141,
                10.99968,
                5.941553,
                9.561794,
                9.174992,
                12.3557,
                11.325861,
                12.095253,
                4.116472,
                15.033502,
                5.2955003,
                3.6110816,
                6.175382,
                10.74008,
                10.91692,
                15.40557,
                12.284264,
                16.144743,
                12.567815,
                14.028609,
                16.806726,
                14.612876,
                17.943256,
                16.082022,
                17.858244,
                17.361277,
                18.055008,
                16.81815,
                13.963133,
                13.338042,
                11.94195,
                10.077962,
                13.209058,
                8.454681,
                9.680062,
                10.61186,
                6.2456894,
                9.020746,
                3.8777738,
                3.0184572,
                6.515422,
                0.038672924,
            ],
            vec![
                0.09561741,
                6.5563784,
                -2.6725936,
                1.1876457,
                1.9365915,
                4.416533,
                -5.3795056,
                1.2769076,
                -6.0004873,
                -2.6329336,
                4.9884143,
                5.378984,
                -7.9673023,
                12.298699,
                16.696861,
                15.491873,
                6.2738566,
                10.584627,
                14.7118845,
                0.54840636,
                15.736791,
                13.422434,
                15.17904,
                17.256449,
                15.526822,
                14.9644785,
                14.486415,
                16.26818,
                13.833231,
                13.754543,
                12.476433,
                12.78299,
                10.746775,
                6.4811716,
                6.8892145,
                -1.4401631,
                0.35716954,
                -0.59143853,
                5.7783365,
                4.959049,
                4.514746,
                8.965897,
                5.0259604,
                7.8060284,
                8.344823,
                6.2071376,
                6.2493534,
                7.9856224,
                8.4859295,
                5.660628,
                5.7422824,
                2.854232,
                0.8136455,
                15.208016,
                8.857441,
                17.164436,
                11.207313,
                17.699806,
                11.005365,
                20.54641,
                14.422905,
                20.429094,
                13.450137,
                20.360537,
                13.401671,
                18.521997,
                14.82664,
                17.555243,
                6.645715,
                17.48721,
                6.4551787,
                14.336542,
                7.2606273,
                5.3290744,
                7.3608913,
                8.053697,
                6.7248073,
                5.1744037,
                5.2074428,
                -1.3688129,
                -2.3737435,
                5.366385,
            ],
            vec![
                1.666534,
                -2.0133758,
                0.01122576,
                -1.861103,
                0.27143353,
                1.6223149,
                1.1447518,
                1.3004985,
                -1.6272843,
                0.10428712,
                1.5828216,
                -0.23673561,
                2.0741496,
                1.1635864,
                -1.9033535,
                -0.120156795,
                1.1250267,
                0.42159596,
                2.512545,
                2.5652704,
                2.9481556,
                -0.71810645,
                2.6394992,
                0.23712501,
                1.1062751,
                1.2469394,
                0.19242239,
                1.7671094,
                1.9391097,
                -0.37342647,
                0.59781325,
                -0.17323819,
                0.74175435,
                0.94068533,
                2.1077933,
                2.1033597,
                2.9768264,
                2.7450604,
                1.8227005,
                2.99542,
                2.4076483,
                0.38178957,
                0.433581,
                0.7217841,
                -0.32690927,
                1.6457819,
                2.3393478,
                1.2538127,
                2.8450933,
                0.7929933,
                2.821195,
                1.2464004,
                2.055948,
                1.3039267,
                1.8322172,
                1.894172,
                1.94645,
                1.026144,
                1.9833946,
                1.8156502,
                2.0707102,
                1.1142877,
                2.034043,
                2.0814497,
                1.6403255,
                2.534144,
                0.73801744,
                1.3562584,
                -0.06696755,
                -0.5553895,
                -0.3246111,
                -0.8916584,
                -0.82411814,
                -0.16832197,
            ],
        ]
    }

    fn expected_lower_orders() -> [Vec<f32>; 3] {
        [
            vec![
                6.098713, 5.0893226, -6.0911913, -6.7819114, 13.251752, 16.14418, 7.6832995,
                2.824732, 2.57026, 10.408969, 8.776338, 13.6783, 17.049034, 12.883465, 16.437056,
                13.877241, 17.90339, 15.171249, 20.515127, 11.290272, 18.873503, 11.434755,
                15.10293, 6.9988832, 6.1085105, 4.67873, 10.638524, 5.890972, 6.7394624, 13.415707,
                8.41757, 14.798344, 8.752231, 12.645353, 13.556108, 11.153498, 10.764422, 8.051168,
                9.448361, 1.8290555, 5.50064, -1.0873742, 2.4527583, -3.7424753, 9.193717,
                12.728172, 11.555926, 12.154267, 12.407433, 12.440878, 13.24142, 19.530802,
                14.173456, 18.568853, 18.288828, 14.009648, 15.426229, 13.139398, 9.925323,
                15.959889, 12.66168, 17.76036, 18.097197, 12.762459, 16.482891, 13.402548,
                15.276099, 5.2670703, 10.892589, 5.1797013, 2.037821, 0.4595784, 6.629072,
            ],
            vec![
                7.5373755, 5.2636786, -4.0441113, -0.219062, -3.8480031, -2.724658, 0.1937184,
                0.79162085, -0.9200721, 0.73866636, 0.23906446, -1.7181942, -1.3604406, 11.220024,
                13.134903, 4.1412187, 6.893021, 5.9513283, 5.6956882, 4.283407, 18.415913,
                15.812028, 12.019984, 12.720299, 11.6742935, 14.602946, 11.850225, 16.08868,
                13.058862, 12.049989, 13.350681, 14.422564, 6.2481685, 12.613721, 4.6700373,
                5.932982, 4.415599, 8.652459, 11.634583, 8.810646, 9.113321, 7.39036, 13.376797,
                7.2179737, 13.861573, 8.081185, 11.402679, 5.696314, 12.447605, 3.9433446,
                10.669755, -2.1491466, 1.6634929, 13.823395, 12.671909, 11.162656, 16.886913,
                11.371497, 12.741075, 13.818545, 13.474773, 9.776059, 17.322653, 12.136895,
                10.784781, 12.577827, 12.806705, 14.556214, 11.48639, 12.880192, 14.534845,
                13.574595, 13.220472, 10.082883, 8.296055, 6.893231, 6.493497, 4.6483345, 8.808964,
                7.031999, 3.3545926, -3.138451,
            ],
            vec![
                1.6276367,
                -1.4149985,
                -1.091316,
                -0.8981575,
                -1.054192,
                1.6319721,
                0.63056993,
                0.5369064,
                1.1153008,
                2.5401762,
                2.1437356,
                -2.7949736,
                1.1789296,
                1.3213018,
                1.1967282,
                2.3320713,
                0.41178164,
                1.5173478,
                0.9164593,
                -0.33791587,
                1.5931625,
                -0.7494801,
                1.7172582,
                0.38566688,
                -0.41061056,
                1.1747347,
                0.75770557,
                0.28779185,
                -0.69291174,
                1.0371411,
                1.464157,
                1.0288159,
                2.2424226,
                0.44533324,
                1.8660605,
                1.4923036,
                1.903705,
                2.092596,
                2.038974,
                1.9947213,
                2.3540473,
                1.7993515,
                1.1769114,
                0.5910229,
                0.79860127,
                2.0383267,
                0.07885873,
                2.8486135,
                0.7498114,
                2.5733194,
                0.40426755,
                0.60731584,
                1.4975159,
                1.0230796,
                1.0741208,
                1.6184679,
                -0.3402709,
                1.8565687,
                -0.39435932,
                2.114923,
                -0.52024305,
                1.0006644,
                1.52127,
                2.0763211,
                2.1403327,
                0.7443464,
                1.7730019,
                0.87550855,
                1.419629,
                0.20301312,
                0.6875506,
                0.5115834,
                -1.046722,
                1.1331804,
            ],
        ]
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
    fn test_analyze_frame_basic() {
        let (system, mut data) = prepare_data_for_tests(None);

        analyze_frame(&system, &mut data).unwrap();
        let expected_total_orders = expected_total_orders();

        for (m, molecule) in data.molecule_types().iter().enumerate() {
            let n_instances = molecule.order_bonds().bond_types()[0].bonds().len();
            let orders = molecule
                .order_bonds()
                .bond_types()
                .iter()
                .map(|b| b.total().order())
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

    #[test]
    fn test_analyze_frame_leaflets() {
        let classifier = LeafletClassification::global("@membrane", "name P");

        let (system, mut data) = prepare_data_for_tests(Some(classifier));

        analyze_frame(&system, &mut data).unwrap();
        let expected_total_orders = expected_total_orders();
        let expected_upper_orders = expected_upper_orders();
        let expected_lower_orders = expected_lower_orders();
        let expected_total_samples = [131, 128, 15];
        let expected_upper_samples = [65, 64, 8];
        let expected_lower_samples = [66, 64, 7];

        for (m, molecule) in data.molecule_types().iter().enumerate() {
            let total_orders = collect_bond_data(&molecule, |b| b.total().order());
            let upper_orders =
                collect_bond_data(&molecule, |b| b.upper().as_ref().unwrap().order());
            let lower_orders =
                collect_bond_data(&molecule, |b| b.lower().as_ref().unwrap().order());
            let total_samples = collect_bond_data(&molecule, |b| b.total().n_samples());
            let upper_samples =
                collect_bond_data(&molecule, |b| b.upper().as_ref().unwrap().n_samples());
            let lower_samples =
                collect_bond_data(&molecule, |b| b.lower().as_ref().unwrap().n_samples());

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
