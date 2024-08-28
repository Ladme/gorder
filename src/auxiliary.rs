// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Implementations of some commonly used functions.

use groan_rs::{errors::GroupError, system::System};

use crate::errors::TopologyError;

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
            "Unexpected error `{}` returned when selecting '{}' using the query '{}'.",
            e, group, query,
        ),
    }

    if system
        .group_isempty(&group_name)
        .expect("Group should exist.")
    {
        Err(TopologyError::EmptyGroup(group.to_owned()))
    } else {
        Ok(())
    }
}

/// Check whether the provided two groups share any atoms. Returns `true` if at least one atom is shared.
///
/// ## Warning
/// Only works if the `System` originates from a tpr file.
pub(crate) fn any_shared(system: &System, group1: &str, group2: &str) -> bool {
    let group_name1 = format!("{}{}", GORDER_GROUP_PREFIX, group1);
    let group_name2 = format!("{}{}", GORDER_GROUP_PREFIX, group2);

    for atom in system.group_iter(&group_name1).unwrap() {
        // we get atom index from atom number
        // this only works because we are reading the system from a tpr file
        // in which the atoms are numbered sequentially
        let atom_index = atom.get_atom_number() - 1;
        if system.group_isin(&group_name2, atom_index).unwrap() {
            return true;
        }
    }

    false
}
