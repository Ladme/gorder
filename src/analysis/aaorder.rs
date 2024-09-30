// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the atomistic lipid order parameters.

use super::leaflets::MoleculeLeafletClassification;
use super::{auxiliary::macros::group_name, topology::SystemTopology};
use crate::errors::{AnalysisError, TopologyError};
use crate::presentation::aapresenter::AAOrderResults;
use crate::{Analysis, PANIC_MESSAGE};

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
    let _ = system.group_intersection(
        group_name!("HeavyAtoms"),
        group_name!("Hydrogens"),
        group_name!("Intersection"),
    );
    if system
        .group_get_n_atoms(group_name!("Intersection"))
        .expect(PANIC_MESSAGE)
        != 0
    {
        return Err(Box::from(TopologyError::AtomsOverlap {
            name1: "HeavyAtoms".to_owned(),
            query1: analysis
                .heavy_atoms()
                .as_ref()
                .expect(PANIC_MESSAGE)
                .clone(),
            name2: "Hydrogens".to_owned(),
            query2: analysis.hydrogens().as_ref().expect(PANIC_MESSAGE).clone(),
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
        if mol.order_bonds().bonds().is_empty() {
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

    Ok(())
}

/// Analyze order parameters in a single simulation frame.
fn analyze_frame(frame: &System, data: &mut SystemTopology) -> Result<(), AnalysisError> {
    let simbox = frame.get_box_as_ref().ok_or(AnalysisError::UndefinedBox)?;

    if !simbox.is_orthogonal() {
        return Err(AnalysisError::NotOrthogonalBox);
    }

    if simbox.is_zero() {
        return Err(AnalysisError::ZeroBox);
    }

    let membrane_normal = data.membrane_normal().into();
    let membrane_center = {
        match data
            .molecules()
            .get(0)
            .expect(PANIC_MESSAGE)
            .leaflet_classification()
        {
            Some(MoleculeLeafletClassification::Global(_)) => frame
                .group_get_center(group_name!("Membrane"))
                .map(Some)
                .map_err(|_| AnalysisError::InvalidGlobalMembraneCenter),
            _ => Ok(None), // No applicable classification, return Ok with None
        }
    }?;

    for molecule in data.molecules_mut().iter_mut() {
        match molecule.leaflet_classification_mut() {
            Some(MoleculeLeafletClassification::Global(x)) => {
                x.set_membrane_center(membrane_center.clone().expect(PANIC_MESSAGE));
            }
            Some(MoleculeLeafletClassification::Local(x)) => {
                x.set_membrane_center(frame, membrane_normal)?;
            }
            Some(MoleculeLeafletClassification::Individual(_)) | None => (),
        };

        molecule.analyze_frame(frame, &simbox, &membrane_normal.into())?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use approx::assert_relative_eq;
    use groan_rs::prelude::Dimension;

    use crate::{
        analysis::molecule::{Bond, MoleculeType},
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
                -13.541704, 6.974998, 24.937578, 28.719318, 25.008698, 26.961494, 9.811056,
                9.057474, 21.543617, 2.4329069, 11.099975, 36.33753, 32.78751, 27.270025, 18.58691,
                6.063838, 22.069092, 30.240564, 19.826765, 28.786337, 12.905661, 11.512761,
                12.237981, 31.867886, 24.438534, 25.97152, 22.511253, 1.5912378, 18.815758,
                24.914057, 28.634373, 19.075571, 29.923029, -11.096256, 28.552168, 20.534142,
                20.328485, 27.31636, 19.751911, 11.389851, 23.082615, 26.743479, 28.175156,
                29.702309, 30.524326, 23.120123, 10.441488, 19.377031, 36.512115, 28.59851,
                23.57948, 22.650856, 25.999723, 25.887957, 31.194407, 23.117897, 30.966377,
                11.654347, 25.653984, 9.882705, 5.8786674, 19.91333, 5.056279, 31.363668, 19.9338,
                34.37085, 23.645094, 6.6677446, 4.208125, 27.989347, 14.978862, 5.945527,
                31.064243,
            ],
            vec![
                29.567423, 11.559245, 18.132113, 30.205153, 13.769697, 11.82006, 18.402748,
                34.152706, 29.97675, 21.529348, 4.772769, 25.827118, 27.201118, -6.7167015,
                -1.9114108, 14.946924, 25.80453, 32.497437, 13.166875, 17.412926, 28.327091,
                14.288321, 20.481098, 9.822741, 31.099829, 16.994942, 20.933537, -9.327742,
                -5.1857877, 26.336637, 34.364956, 29.071299, 8.061023, 29.831762, 27.633347,
                2.2279348, 30.77279, 23.746437, 4.831814, 29.031404, 2.477138, 5.663187, 13.218306,
                30.367407, 15.024002, 19.633091, 15.41196, 20.990025, 28.09423, 19.094893,
                -1.8942664, 26.892096, 15.65695, 0.9685844, 23.518724, 9.603974, 2.0685284,
                27.89767, 0.9808493, 27.199028, 3.6607904, 24.186451, -6.9205585, 1.6918764,
                27.911137, 13.681939, 16.535954, 20.407568, 0.7050859, 27.205551, 17.652033,
                22.206394, 32.356857, 16.412033, 16.356262, 32.111458, 7.6329927, 4.492818,
                13.628067, 5.227478, 29.234451, 14.016406,
            ],
            vec![
                1.8537161,
                2.8827128,
                -0.3523764,
                3.8616748,
                2.0549011,
                2.0619702,
                3.7806582,
                3.2784903,
                2.2273543,
                2.4182067,
                4.356757,
                4.8376565,
                0.36293957,
                -0.38007507,
                2.3270063,
                2.114952,
                3.5553129,
                3.2530794,
                -1.4675866,
                2.9063377,
                2.2119145,
                1.2461977,
                4.541318,
                0.47169197,
                -1.0800903,
                4.9901414,
                3.9305735,
                1.3128071,
                1.7753217,
                3.553464,
                3.254287,
                -2.7592607,
                2.181141,
                1.5368084,
                -0.7066252,
                3.6841085,
                3.3663125,
                1.550467,
                0.6637146,
                -3.4283745,
                1.6104923,
                1.837405,
                2.4216743,
                2.9841769,
                3.7265575,
                1.9389439,
                4.1024256,
                3.2254627,
                1.3526614,
                0.9648583,
                0.8555777,
                -0.78275836,
                -3.0317092,
                1.3860189,
                3.51264,
                2.4848883,
                2.5110192,
                2.231767,
                3.9738538,
                0.95012796,
                3.2941704,
                3.5949047,
                1.606179,
                1.5890352,
                -1.8708402,
                4.880531,
                2.6444633,
                0.6956645,
                0.622792,
                -0.5119834,
                4.157771,
                3.4290047,
                4.761696,
                3.595663,
            ],
        ]
    }

    fn expected_upper_orders() -> [Vec<f32>; 3] {
        [
            vec![
                -7.4505186,
                6.515422,
                8.454681,
                14.842073,
                12.567815,
                15.40557,
                3.7025452,
                3.8777738,
                12.095253,
                6.175382,
                5.2090025,
                16.806726,
                17.361277,
                14.028609,
                5.941553,
                3.6110816,
                11.660126,
                9.725433,
                12.143466,
                14.612876,
                10.080931,
                6.2456894,
                1.5994557,
                17.858244,
                12.284264,
                13.209058,
                11.220982,
                -7.1850986,
                10.39819,
                12.030589,
                13.531443,
                5.659866,
                13.963133,
                -4.3143415,
                16.144743,
                15.033502,
                9.174992,
                13.638062,
                10.99968,
                5.291138,
                9.680062,
                16.81815,
                10.077962,
                11.94195,
                13.475298,
                12.3557,
                7.8712277,
                11.325861,
                17.943256,
                13.427263,
                8.781141,
                11.2161,
                13.338042,
                10.61186,
                18.055008,
                9.561794,
                13.062992,
                4.9148817,
                9.5098,
                5.2039742,
                0.7893448,
                9.020746,
                3.0184572,
                14.92661,
                10.74008,
                16.082022,
                10.91692,
                0.038672924,
                5.2955003,
                14.737594,
                7.9799776,
                4.116472,
                12.190741,
            ],
            vec![
                14.9644785,
                6.8892145,
                6.645715,
                20.429094,
                4.959049,
                6.5563784,
                5.0259604,
                15.736791,
                17.256449,
                8.857441,
                0.35716954,
                12.476433,
                15.526822,
                -2.6725936,
                1.9365915,
                8.053697,
                13.754543,
                20.360537,
                6.2738566,
                5.7783365,
                17.164436,
                6.2071376,
                7.2606273,
                5.1744037,
                18.521997,
                10.746775,
                8.4859295,
                -7.9673023,
                -5.3795056,
                14.486415,
                20.54641,
                17.699806,
                -0.59143853,
                16.696861,
                14.82664,
                5.366385,
                13.450137,
                11.005365,
                0.54840636,
                15.208016,
                0.8136455,
                -1.3688129,
                6.7248073,
                17.48721,
                7.8060284,
                15.491873,
                5.3290744,
                6.4551787,
                11.207313,
                6.4811716,
                -2.6329336,
                13.833231,
                7.3608913,
                1.1876457,
                12.298699,
                5.660628,
                1.2769076,
                14.422905,
                -2.3737435,
                15.17904,
                5.378984,
                13.401671,
                -6.0004873,
                4.416533,
                14.336542,
                7.9856224,
                10.584627,
                14.7118845,
                2.854232,
                12.78299,
                6.2493534,
                8.344823,
                16.26818,
                5.7422824,
                8.965897,
                17.555243,
                0.09561741,
                -1.4401631,
                4.514746,
                4.9884143,
                13.422434,
                5.2074428,
            ],
            vec![
                1.2464004,
                1.026144,
                -0.5553895,
                1.8227005,
                1.7671094,
                0.59781325,
                1.6403255,
                2.534144,
                2.5652704,
                2.3393478,
                2.6394992,
                2.7450604,
                -0.3246111,
                -0.8916584,
                1.3039267,
                1.1142877,
                2.034043,
                2.0741496,
                -0.71810645,
                1.8322172,
                -0.120156795,
                1.9391097,
                2.9481556,
                -0.32690927,
                0.01122576,
                2.99542,
                1.8156502,
                0.7217841,
                1.1447518,
                2.055948,
                1.6223149,
                -1.861103,
                0.38178957,
                1.1250267,
                -1.9033535,
                1.6457819,
                0.7929933,
                2.0707102,
                -0.37342647,
                -2.0133758,
                0.433581,
                1.3004985,
                1.2469394,
                0.74175435,
                1.5828216,
                0.42159596,
                1.2538127,
                2.821195,
                -0.06696755,
                -0.16832197,
                -0.17323819,
                0.27143353,
                -0.23673561,
                0.94068533,
                1.894172,
                1.1635864,
                0.73801744,
                1.3562584,
                2.1077933,
                0.19242239,
                1.666534,
                2.8450933,
                1.94645,
                1.9833946,
                -0.82411814,
                2.9768264,
                0.10428712,
                1.1062751,
                0.23712501,
                -1.6272843,
                2.0814497,
                2.512545,
                2.4076483,
                2.1033597,
            ],
        ]
    }

    fn expected_lower_orders() -> [Vec<f32>; 3] {
        [
            vec![
                -6.0911913, 0.4595784, 16.482891, 13.877241, 12.440878, 11.555926, 6.1085105,
                5.1797013, 9.448361, -3.7424753, 5.890972, 19.530802, 15.426229, 13.24142,
                12.645353, 2.4527583, 10.408969, 20.515127, 7.6832995, 14.173456, 2.824732,
                5.2670703, 10.638524, 14.009648, 12.154267, 12.762459, 11.290272, 8.776338,
                8.41757, 12.883465, 15.10293, 13.415707, 15.959889, -6.7819114, 12.407433, 5.50064,
                11.153498, 13.6783, 8.752231, 6.098713, 13.402548, 9.925323, 18.097197, 17.76036,
                17.049034, 10.764422, 2.57026, 8.051168, 18.568853, 15.171249, 14.798344,
                11.434755, 12.66168, 15.276099, 13.139398, 13.556108, 17.90339, 6.7394624,
                16.14418, 4.67873, 5.0893226, 10.892589, 2.037821, 16.437056, 9.193717, 18.288828,
                12.728172, 6.629072, -1.0873742, 13.251752, 6.9988832, 1.8290555, 18.873503,
            ],
            vec![
                14.602946, 4.6700373, 11.48639, 9.776059, 8.810646, 5.2636786, 13.376797,
                18.415913, 12.720299, 12.671909, 4.415599, 13.350681, 11.6742935, -4.0441113,
                -3.8480031, 6.893231, 12.049989, 12.136895, 6.893021, 11.634583, 11.162656,
                8.081185, 13.220472, 4.6483345, 12.577827, 6.2481685, 12.447605, -1.3604406,
                0.1937184, 11.850225, 13.818545, 11.371497, 8.652459, 13.134903, 12.806705,
                -3.138451, 17.322653, 12.741075, 4.283407, 13.823395, 1.6634929, 7.031999,
                6.493497, 12.880192, 7.2179737, 4.1412187, 10.082883, 14.534845, 16.886913,
                12.613721, 0.73866636, 13.058862, 8.296055, -0.219062, 11.220024, 3.9433446,
                0.79162085, 13.474773, 3.3545926, 12.019984, -1.7181942, 10.784781, -0.9200721,
                -2.724658, 13.574595, 5.696314, 5.9513283, 5.6956882, -2.1491466, 14.422564,
                11.402679, 13.861573, 16.08868, 10.669755, 7.39036, 14.556214, 7.5373755, 5.932982,
                9.113321, 0.23906446, 15.812028, 8.808964,
            ],
            vec![
                0.60731584,
                1.8565687,
                0.20301312,
                2.038974,
                0.28779185,
                1.464157,
                2.1403327,
                0.7443464,
                -0.33791587,
                0.07885873,
                1.7172582,
                2.092596,
                0.6875506,
                0.5115834,
                1.0230796,
                1.0006644,
                1.52127,
                1.1789296,
                -0.7494801,
                1.0741208,
                2.3320713,
                -0.69291174,
                1.5931625,
                0.79860127,
                -1.091316,
                1.9947213,
                2.114923,
                0.5910229,
                0.63056993,
                1.4975159,
                1.6319721,
                -0.8981575,
                1.7993515,
                0.41178164,
                1.1967282,
                2.0383267,
                2.5733194,
                -0.52024305,
                1.0371411,
                -1.4149985,
                1.1769114,
                0.5369064,
                1.1747347,
                2.2424226,
                2.1437356,
                1.5173478,
                2.8486135,
                0.40426755,
                1.419629,
                1.1331804,
                1.0288159,
                -1.054192,
                -2.7949736,
                0.44533324,
                1.6184679,
                1.3213018,
                1.7730019,
                0.87550855,
                1.8660605,
                0.75770557,
                1.6276367,
                0.7498114,
                -0.3402709,
                -0.39435932,
                -1.046722,
                1.903705,
                2.5401762,
                -0.41061056,
                0.38566688,
                1.1153008,
                2.0763211,
                0.9164593,
                2.3540473,
                1.4923036,
            ],
        ]
    }

    fn collect_bond_data<T, F>(molecule: &MoleculeType, func: F) -> Vec<T>
    where
        F: Fn(&Bond) -> T,
    {
        molecule.order_bonds().bonds().iter().map(func).collect()
    }

    #[test]
    fn test_analyze_frame_basic() {
        let (system, mut data) = prepare_data_for_tests(None);

        analyze_frame(&system, &mut data).unwrap();
        let expected_total_orders = expected_total_orders();

        for (m, molecule) in data.molecules().iter().enumerate() {
            let n_instances = molecule.order_bonds().bonds()[0].bonds().len();
            let orders = molecule
                .order_bonds()
                .bonds()
                .iter()
                .map(|b| b.total().order())
                .collect::<Vec<f32>>();

            let samples = molecule
                .order_bonds()
                .bonds()
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

        for (m, molecule) in data.molecules().iter().enumerate() {
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
