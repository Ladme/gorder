// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Contains the implementation of the calculation of the atomistic lipid order parameters.

use std::collections::HashSet;

use groan_rs::{files::FileType, structures::group::Group, system::System};

use crate::{
    auxiliary::{macros::group_name, GORDER_GROUP_PREFIX, PANIC_MESSAGE},
    errors::TopologyError,
    topology::{AtomType, BondType, MoleculeType},
    Analysis,
};

struct OrderData {
    molecules: Vec<MoleculeType>,
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
    let molecules = crate::auxiliary::classify_molecules(&system, "HeavyAtoms", "Hydrogens", analysis.leaflets().as_ref())?;

    //system.traj_iter_map_reduce(trajectory_file, n_threads, body, start_time, end_time, step, progress_printer);

    Ok(())
}

fn analyze_frame(frame: &System, data: &Vec<MoleculeType>) {
}

/*/// Modifies the HeavyAtoms group so it only contains heavy atoms that are connected to at least one hydrogen.
/// Modifies the Hydrogens group so it only contains hydrogens that are connected to exactly one heavy atom.
/// In case a hydrogen atom is connected to multiple heavy atoms, raises an error.
///
/// ## Warning
/// Only works if the `System` originates from a tpr file.
fn filter_unconnected(system: &mut System, heavy_atoms: &str, hydrogens: &str) -> Result<(), TopologyError> {
    let heavy_atoms_name = format!("{}{}", GORDER_GROUP_PREFIX, heavy_atoms);
    let hydrogens_name = format!("{}{}", GORDER_GROUP_PREFIX, hydrogens);

    let mut bonded_heavy_atoms = HashSet::new();
    let mut bonded_hydrogens = HashSet::new();
    for atom in system.group_iter(&heavy_atoms_name).expect(PANIC_MESSAGE) {
        for bonded in atom.get_bonded().iter() {
            if system.group_isin(&hydrogens_name, bonded).expect(PANIC_MESSAGE) {
                // add hydrogen to the list of bonded hydrogens
                if !bonded_hydrogens.insert(bonded) {
                    return Err(TopologyError::MultiBinding(bonded, atom.get_atom_number()))
                }

                // add heavy atom to the list of bonded heavy atoms
                bonded_heavy_atoms.insert(atom.get_atom_number());
            }
        }
    }

    // print info about the filtering
    let original_heavy_atoms_n = system.group_get_n_atoms(&heavy_atoms_name).expect(PANIC_MESSAGE);
    let original_hydrogens_n = system.group_get_n_atoms(&hydrogens_name).expect(PANIC_MESSAGE);
    log::info!("{} out of {} selected heavy atoms are bonded to a hydrogen. Ignoring the rest.", original_heavy_atoms_n, bonded_heavy_atoms.len());
    log::info!("{} out of {} selected hydrogens are bonded to a heavy atom. Ignoring the rest.", original_hydrogens_n, bonded_hydrogens.len());

    Ok(())
}*/
