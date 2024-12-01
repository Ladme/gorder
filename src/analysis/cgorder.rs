// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the coarse-grained lipid order parameters.

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
    presentation::cgpresenter::CGOrderResults,
    Analysis, PANIC_MESSAGE,
};

pub(super) fn analyze_coarse_grained(
    analysis: &Analysis,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut system = System::from_file_with_format(analysis.structure(), FileType::TPR)?;
    log::info!("Read molecular topology from '{}'.", analysis.structure());

    super::common::create_group(
        &mut system,
        "Beads",
        analysis.beads().as_ref().unwrap_or_else(||
             panic!("FATAL GORDER ERROR | aaorder::analyze_coarse_grained | Selection of order beads should be provided. {}", PANIC_MESSAGE)),
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
    )?;

    if !sanity_check_molecules(&molecules) {
        return Ok(());
    }

    let data = SystemTopology::new(molecules, analysis.membrane_normal().into());
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
    result.write_ordermaps_bonds()?;

    // write out the results
    let results = CGOrderResults::from(result);
    write_results(results, &analysis)?;

    Ok(())
}
