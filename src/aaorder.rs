// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the atomistic lipid order parameters.

use std::{collections::HashSet, ops::Add};

use groan_rs::{
    files::FileType,
    prelude::{ProgressPrinter, SimBox, Vector3D, XtcReader},
    structures::group::Group,
    system::System,
};

use crate::{
    auxiliary::{macros::group_name, GORDER_GROUP_PREFIX, PANIC_MESSAGE},
    errors::{AnalysisError, TopologyError},
    topology::{AtomType, BondType, MoleculeType},
    Analysis,
};

#[derive(Debug, Clone)]
struct SystemTopology {
    molecules: Vec<MoleculeType>,
    membrane_normal: Vector3D,
}

impl Add<SystemTopology> for SystemTopology {
    type Output = Self;
    fn add(self, rhs: SystemTopology) -> Self::Output {
        Self {
            molecules: self
                .molecules
                .into_iter()
                .zip(rhs.molecules.into_iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<MoleculeType>>(),
            membrane_normal: self.membrane_normal,
        }
    }
}

/// Calculate the atomistic lipid order parameters.
pub(crate) fn analyze_atomistic(
    analysis: &Analysis,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut system = System::from_file_with_format(analysis.structure(), FileType::TPR)?;

    crate::auxiliary::create_group(
        &mut system,
        "HeavyAtoms",
        analysis.heavy_atoms().as_ref().unwrap_or_else(||
             panic!("FATAL GORDER ERROR | aaorder::analyze_atomistic | Selection of heavy atoms should be provided. {}", PANIC_MESSAGE)),
    )?;

    crate::auxiliary::create_group(
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
    let molecules = crate::auxiliary::classify_molecules(
        &system,
        "HeavyAtoms",
        "Hydrogens",
        analysis.leaflets().as_ref(),
    )?;

    let data = SystemTopology {
        molecules: molecules,
        membrane_normal: analysis.membrane_normal().into(),
    };

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

fn analyze_frame(frame: &System, data: &mut SystemTopology) -> Result<(), AnalysisError> {
    let simbox = frame.get_box_as_ref().ok_or(AnalysisError::UndefinedBox)?;

    if !simbox.is_orthogonal() {
        return Err(AnalysisError::NotOrthogonalBox);
    }

    if simbox.is_zero() {
        return Err(AnalysisError::ZeroBox);
    }

    for molecule in data.molecules.iter_mut() {
        molecule.analyze_frame(frame, &simbox, &data.membrane_normal)?;
    }

    Ok(())
}
