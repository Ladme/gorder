// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the atomistic lipid order parameters.

use super::leaflets::MoleculeLeafletClassification;
use super::{auxiliary::macros::group_name, topology::SystemTopology};
use crate::errors::{AnalysisError, TopologyError};
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

    super::auxiliary::create_group(
        &mut system,
        "HeavyAtoms",
        analysis.heavy_atoms().as_ref().unwrap_or_else(||
             panic!("FATAL GORDER ERROR | aaorder::analyze_atomistic | Selection of heavy atoms should be provided. {}", PANIC_MESSAGE)),
    )?;

    super::auxiliary::create_group(
        &mut system,
        "Hydrogens", 
        analysis.hydrogens().as_ref().unwrap_or_else(||
            panic!("FATAL GORDER ERROR | aaorder::analyze_atomistic | Selection of hydrogens should be provided. {}", PANIC_MESSAGE))
    )?;

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

    // get the relevant molecules
    let molecules = super::auxiliary::classify_molecules(
        &system,
        "HeavyAtoms",
        "Hydrogens",
        analysis.leaflets().as_ref(),
        analysis.membrane_normal().into(),
        analysis.map().as_ref(),
    )?;

    // if no molecules are detected, end the analysis
    if molecules.len() == 0 {
        log::info!("No molecules detected for analysis.");
        return Ok(());
    }

    let data = SystemTopology::new(molecules, analysis.membrane_normal().into());

    // run the analysis in parallel
    let result = system.traj_iter_map_reduce::<XtcReader, SystemTopology, AnalysisError>(
        analysis.trajectory(),
        analysis.n_threads(),
        analyze_frame,
        data,
        Some(analysis.start()),
        Some(analysis.end()),
        Some(analysis.step()),
        Some(ProgressPrinter::new()),
    )?;

    // todo! average & convert to output structure
    // todo! write out

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
