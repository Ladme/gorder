// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementations of some commonly used functions.

use hashbrown::HashSet;

use groan_rs::{errors::GroupError, system::System};

use crate::{
    errors::TopologyError,
    molecule::{AtomType, BondType, Molecule, MoleculeTopology, OrderAtoms, OrderBonds},
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
fn classify_molecules(
    system: &System,
    group1: &str,
    group2: &str,
) -> Result<Vec<Molecule>, TopologyError> {
    let group1_name = format!("{}{}", GORDER_GROUP_PREFIX, group1);
    let group2_name = format!("{}{}", GORDER_GROUP_PREFIX, group2);

    let mut visited = HashSet::new();
    let mut molecules: Vec<Molecule> = Vec::new();
    for atom in system.group_iter(&group1_name).expect(PANIC_MESSAGE) {
        let index = atom.get_atom_number() - 1;
        if !visited.insert(index) {
            continue;
        }

        let mut all_bonds = HashSet::new();
        let mut residues = Vec::new();
        let mut order_atoms = Vec::new();
        let mut minimum_index = index;

        // iterate through the molecule
        for atom2 in system.molecule_iter(index).expect(PANIC_MESSAGE) {
            let index2 = atom2.get_atom_number() - 1;
            visited.insert(index2);

            if index2 < minimum_index {
                minimum_index = index2;
            }

            if !residues.contains(&atom2.get_residue_name()) {
                residues.push(atom2.get_residue_name());
            }

            if system
                .group_isin(&group1_name, index2)
                .expect(PANIC_MESSAGE)
            {
                order_atoms.push(index2);
            }

            for bonded in system.bonded_atoms_iter(index2).expect(PANIC_MESSAGE) {
                let index3 = bonded.get_atom_number() - 1;

                if visited.contains(&index3) {
                    continue;
                }

                if !all_bonds.insert((index2, index3)) {
                    panic!("FATAL GORDER ERROR | auxiliary::classify_molecules | Bond between `atom` and `atom2` (indices: {} and {}) encountered multiple times in the molecule. {}",
                    index, index2, PANIC_MESSAGE);
                }
            }
        }

        // select order bonds
        let mut order_bonds = HashSet::new();
        for &(a1, a2) in all_bonds.iter() {
            if (system.group_isin(&group1_name, a1).expect(PANIC_MESSAGE)
                && system.group_isin(&group2_name, a2).expect(PANIC_MESSAGE))
                || (system.group_isin(&group2_name, a1).expect(PANIC_MESSAGE)
                    && system.group_isin(&group1_name, a2).expect(PANIC_MESSAGE))
            {
                if !order_bonds.insert((a1, a2)) {
                    panic!("FATAL GORDER ERROR | auxiliary::classify_molecules | Order bond between '{}' and '{}' encountered multiple times in the molecule. {}", a1, a2, PANIC_MESSAGE);
                }
            }
        }

        // add molecule to vector of molecules, if it does not already exist
        let topology = MoleculeTopology::new(system, &all_bonds, minimum_index);

        let mut add_molecule = true;
        for molecule in molecules.iter_mut() {
            if *molecule.topology() == topology {
                molecule.add(system, &order_bonds, minimum_index);
                add_molecule = false;
                break;
            }
        }

        if add_molecule {
            let name = residues.join("-");
            molecules.push(Molecule::new(
                system,
                &name,
                &topology,
                &order_bonds,
                &order_atoms,
                minimum_index,
            ));
        }
    }

    Ok(molecules)
}

#[cfg(test)]
mod tests {
    use groan_rs::files::FileType;

    use super::*;

    #[test]
    fn test_classify_molecules() {
        let mut system =
            System::from_file_with_format("tests/files/pcpepg.tpr", FileType::TPR).unwrap();

        create_group(
            &mut system,
            "HeavyAtoms",
            "@membrane and element name carbon",
        )
        .unwrap();
        create_group(
            &mut system,
            "Hydrogens",
            "@membrane and element name hydrogen",
        )
        .unwrap();

        let molecules = classify_molecules(&system, "HeavyAtoms", "Hydrogens").unwrap();
        println!("{}", molecules.len());

        for molecule in molecules {
            println!("MOLECULE: {}", molecule.name());
            for bond in molecule.topology().bonds.iter() {
                println!("{}-{}", &bond.atom1.atom_name, &bond.atom2.atom_name);
            }

            for order_bond in molecule.order_bonds().bonds.iter() {
                println!(
                    "order: {}-{}",
                    &order_bond.bond_type.atom1.atom_name, &order_bond.bond_type.atom2.atom_name
                );
            }
        }
    }
}
