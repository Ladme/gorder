// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementations of some commonly used functions.

use std::collections::HashSet;

use groan_rs::{errors::GroupError, system::System};

use crate::{
    errors::TopologyError,
    molecule::{AtomType, Molecule},
};

/// A prefix used as an identifier for Gorder groups.
pub(crate) const GORDER_GROUP_PREFIX: &str = "xxxGorderReservedxxx-";

#[macro_use]
pub(crate) mod macros {
    macro_rules! group_name {
        ($group:expr) => {
            concat!("xxxGorderReservedxxx-", $group)
        };
    }

    pub(crate) use group_name;
}

/// Message that should be added to every panic.
pub(crate) const PANIC_MESSAGE: &str =
    "\n\n\n            >>> THIS SHOULD NOT HAVE HAPPENED! PLEASE REPORT THIS ERROR <<<
(open an issue at 'github.com/Ladme/gorder/issues' or write an e-mail to 'ladmeb@gmail.com')\n\n";

/// Create group handling all potential errors. Also check that the group is not empty.
pub(crate) fn create_group(
    system: &mut System,
    group: &str,
    query: &str,
) -> Result<(), TopologyError> {
    let group_name = format!("{}{}", GORDER_GROUP_PREFIX, group);

    match system.group_create(&group_name, query) {
        Ok(_) | Err(GroupError::AlreadyExistsWarning(_)) => (),
        Err(GroupError::InvalidQuery(_)) => {
            return Err(TopologyError::InvalidQuery(query.to_owned()))
        }
        Err(e) => panic!(
            "FATAL GORDER ERROR | auxiliary::create_group | Unexpected error `{}` returned when selecting '{}' using the query '{}'. {}",
            e, group, query, PANIC_MESSAGE
        ),
    }

    if system.group_isempty(&group_name).unwrap_or_else(|_| {
        panic!(
            "FATAL GORDER ERROR | auxiliary::create_group | Group '{}' should exist. {}",
            group, PANIC_MESSAGE,
        )
    }) {
        Err(TopologyError::EmptyGroup(group.to_owned()))
    } else {
        Ok(())
    }
}

/// ## Warning
/// Only works if the `System` originates from a tpr file.
fn classify_molecules(system: &System, group: &str) -> Result<HashSet<Molecule>, TopologyError> {
    let group_name = format!("{}{}", GORDER_GROUP_PREFIX, group);

    let mut visited = HashSet::new();
    let mut molecules = HashSet::new();
    for atom in system.group_iter(&group_name).expect(PANIC_MESSAGE) {
        let index = atom.get_atom_number() - 1;
        if !visited.insert(index) {
            continue;
        }

        let mut bonds = HashSet::new();
        let mut atom_types: Vec<(usize, AtomType)> = vec![(index, AtomType::from(atom))];
        let mut residues = vec![atom.get_residue_name()];

        // iterate through the molecule
        for atom2 in system.molecule_iter(index).expect(PANIC_MESSAGE) {
            let index2 = atom2.get_atom_number() - 1;
            if !visited.insert(index2) {
                panic!("FATAL GORDER ERROR | auxiliary::classify_molecules | `atom2` (index: {}) was visited before `atom` (index: {}) which is part of the same molecule. {}", 
                index2, index, PANIC_MESSAGE);
            }

            atom_types.push((index2, AtomType::from(atom2)));
            if !residues.contains(&atom2.get_residue_name()) {
                residues.push(atom2.get_residue_name());
            }

            if !bonds.insert((index, index2)) {
                panic!("FATAL GORDER ERROR | auxiliary::classify_molecules | Bond between `atom` and `atom2` (indices: {} and {}) encountered multiple times in the molecule. {}",
                index, index2, PANIC_MESSAGE);
            }
        }

        // construct the molecule
        if atom_types.len() == 0 {
            panic!("FATAL GORDER ERROR | auxiliary::classify_molecules | `atom_types` vector must contain at least 1 element. {}", PANIC_MESSAGE);
        }

        let molecule = Molecule::new(&residues.join("-"), bonds, atom_types);
    }

    // attempt to add the molecule into the set of molecules

    Ok(molecules)
}
