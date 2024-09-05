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

/// Classifies molecules in the system into distinct types based on their topology
/// and returns a list of classified molecules.
///
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

        let (all_bonds, residues, order_atoms, minimum_index) =
            process_molecule_dfs(system, &group1_name, index, &mut visited)?;

        let order_bonds = select_order_bonds(system, &group1_name, &group2_name, &all_bonds)?;

        if let Some(existing_molecule) = molecules
            .iter_mut()
            .find(|m| *m.topology() == MoleculeTopology::new(system, &all_bonds, minimum_index))
        {
            existing_molecule.add(system, &order_bonds, minimum_index);
        } else {
            let name = residues.join("-");
            molecules.push(Molecule::new(
                system,
                &name,
                &MoleculeTopology::new(system, &all_bonds, minimum_index),
                &order_bonds,
                &order_atoms,
                minimum_index,
            ));
        }
    }

    Ok(molecules)
}

/// Uses DFS to traverse a molecule, collecting bonds, residues, order atoms,
/// and determining the minimum atom index.
fn process_molecule_dfs(
    system: &System,
    group_name: &str,
    start_index: usize,
    visited: &mut HashSet<usize>,
) -> Result<(HashSet<(usize, usize)>, Vec<String>, Vec<usize>, usize), TopologyError> {
    let mut all_bonds = HashSet::new();
    let mut residues = Vec::new();
    let mut order_atoms = Vec::new();
    let mut minimum_index = start_index;

    let mut stack = vec![start_index];

    while let Some(index) = stack.pop() {
        let atom = system.get_atom_as_ref(index).expect(PANIC_MESSAGE);
        visited.insert(index);
        minimum_index = minimum_index.min(index);

        if !residues.contains(&atom.get_residue_name().to_owned()) {
            residues.push(atom.get_residue_name().to_owned());
        }

        if system.group_isin(group_name, index).expect(PANIC_MESSAGE) {
            order_atoms.push(index);
        }

        for bonded in system.bonded_atoms_iter(index).expect(PANIC_MESSAGE) {
            let bonded_index = bonded.get_atom_number() - 1;
            if visited.contains(&bonded_index) {
                continue;
            }

            all_bonds.insert((index, bonded_index));
            stack.push(bonded_index);
        }
    }

    Ok((all_bonds, residues, order_atoms, minimum_index))
}

/// Selects bonds between two groups that should be considered as order bonds,
/// ensuring they appear only once in the molecule.
fn select_order_bonds(
    system: &System,
    group1_name: &str,
    group2_name: &str,
    all_bonds: &HashSet<(usize, usize)>,
) -> Result<HashSet<(usize, usize)>, TopologyError> {
    let mut order_bonds = HashSet::new();

    for &(a1, a2) in all_bonds.iter() {
        if (system.group_isin(group1_name, a1).expect(PANIC_MESSAGE)
            && system.group_isin(group2_name, a2).expect(PANIC_MESSAGE))
            || (system.group_isin(group2_name, a1).expect(PANIC_MESSAGE)
                && system.group_isin(group1_name, a2).expect(PANIC_MESSAGE))
        {
            if !order_bonds.insert((a1, a2)) {
                panic!(
                    "FATAL GORDER ERROR | auxiliary::select_order_bonds | Order bond between '{}' and '{}' encountered multiple times in the molecule. {}",
                    a1, a2, PANIC_MESSAGE
                );
            }
        }
    }

    Ok(order_bonds)
}

#[cfg(test)]
mod tests {
    use groan_rs::files::FileType;

    use super::*;

    fn pope_topology() -> MoleculeTopology {
        let string_representation = "0 POPE N 3 POPE HN3
0 POPE N 2 POPE HN2
0 POPE N 4 POPE C12
0 POPE N 1 POPE HN1
4 POPE C12 7 POPE C11
4 POPE C12 6 POPE H12B
4 POPE C12 5 POPE H12A
7 POPE C11 14 POPE O12
7 POPE C11 9 POPE H11B
7 POPE C11 8 POPE H11A
10 POPE P 14 POPE O12
10 POPE P 11 POPE O13
10 POPE P 13 POPE O11
10 POPE P 12 POPE O14
13 POPE O11 15 POPE C1
15 POPE C1 17 POPE HB
15 POPE C1 18 POPE C2
15 POPE C1 16 POPE HA
18 POPE C2 19 POPE HS
18 POPE C2 20 POPE O21
18 POPE C2 26 POPE C3
20 POPE O21 21 POPE C21
21 POPE C21 22 POPE O22
21 POPE C21 23 POPE C22
23 POPE C22 24 POPE H2R
23 POPE C22 25 POPE H2S
23 POPE C22 35 POPE C23
26 POPE C3 29 POPE O31
26 POPE C3 27 POPE HX
26 POPE C3 28 POPE HY
29 POPE O31 30 POPE C31
30 POPE C31 31 POPE O32
30 POPE C31 32 POPE C32
32 POPE C32 33 POPE H2X
32 POPE C32 82 POPE C33
32 POPE C32 34 POPE H2Y
35 POPE C23 37 POPE H3S
35 POPE C23 38 POPE C24
35 POPE C23 36 POPE H3R
38 POPE C24 39 POPE H4R
38 POPE C24 41 POPE C25
38 POPE C24 40 POPE H4S
41 POPE C25 43 POPE H5S
41 POPE C25 42 POPE H5R
41 POPE C25 44 POPE C26
44 POPE C26 45 POPE H6R
44 POPE C26 46 POPE H6S
44 POPE C26 47 POPE C27
47 POPE C27 50 POPE C28
47 POPE C27 48 POPE H7R
47 POPE C27 49 POPE H7S
50 POPE C28 51 POPE H8R
50 POPE C28 52 POPE H8S
50 POPE C28 53 POPE C29
53 POPE C29 54 POPE H91
53 POPE C29 55 POPE C210
55 POPE C210 56 POPE H101
55 POPE C210 57 POPE C211
57 POPE C211 58 POPE H11R
57 POPE C211 59 POPE H11S
57 POPE C211 60 POPE C212
60 POPE C212 62 POPE H12S
60 POPE C212 61 POPE H12R
60 POPE C212 63 POPE C213
63 POPE C213 65 POPE H13S
63 POPE C213 64 POPE H13R
63 POPE C213 66 POPE C214
66 POPE C214 67 POPE H14R
66 POPE C214 68 POPE H14S
66 POPE C214 69 POPE C215
69 POPE C215 71 POPE H15S
69 POPE C215 70 POPE H15R
69 POPE C215 72 POPE C216
72 POPE C216 75 POPE C217
72 POPE C216 74 POPE H16S
72 POPE C216 73 POPE H16R
75 POPE C217 76 POPE H17R
75 POPE C217 78 POPE C218
75 POPE C217 77 POPE H17S
78 POPE C218 79 POPE H18R
78 POPE C218 80 POPE H18S
78 POPE C218 81 POPE H18T
82 POPE C33 84 POPE H3Y
82 POPE C33 85 POPE C34
82 POPE C33 83 POPE H3X
85 POPE C34 86 POPE H4X
85 POPE C34 87 POPE H4Y
85 POPE C34 88 POPE C35
88 POPE C35 89 POPE H5X
88 POPE C35 91 POPE C36
88 POPE C35 90 POPE H5Y
91 POPE C36 92 POPE H6X
91 POPE C36 94 POPE C37
91 POPE C36 93 POPE H6Y
94 POPE C37 95 POPE H7X
94 POPE C37 96 POPE H7Y
94 POPE C37 97 POPE C38
97 POPE C38 100 POPE C39
97 POPE C38 99 POPE H8Y
97 POPE C38 98 POPE H8X
100 POPE C39 101 POPE H9X
100 POPE C39 103 POPE C310
100 POPE C39 102 POPE H9Y
103 POPE C310 104 POPE H10X
103 POPE C310 106 POPE C311
103 POPE C310 105 POPE H10Y
106 POPE C311 109 POPE C312
106 POPE C311 108 POPE H11Y
106 POPE C311 107 POPE H11X
109 POPE C312 110 POPE H12X
109 POPE C312 111 POPE H12Y
109 POPE C312 112 POPE C313
112 POPE C313 115 POPE C314
112 POPE C313 113 POPE H13X
112 POPE C313 114 POPE H13Y
115 POPE C314 117 POPE H14Y
115 POPE C314 116 POPE H14X
115 POPE C314 118 POPE C315
118 POPE C315 121 POPE C316
118 POPE C315 119 POPE H15X
118 POPE C315 120 POPE H15Y
121 POPE C316 124 POPE H16Z
121 POPE C316 123 POPE H16Y
121 POPE C316 122 POPE H16X";

        string2topology(string_representation)
    }

    fn pope_order_bonds() -> HashSet<BondType> {
        let string_representation = "4 POPE C12 5 POPE H12A
4 POPE C12 6 POPE H12B
7 POPE C11 8 POPE H11A
7 POPE C11 9 POPE H11B
15 POPE C1 17 POPE HB
15 POPE C1 16 POPE HA
18 POPE C2 19 POPE HS
23 POPE C22 24 POPE H2R
23 POPE C22 25 POPE H2S
26 POPE C3 27 POPE HX
26 POPE C3 28 POPE HY
32 POPE C32 33 POPE H2X
32 POPE C32 34 POPE H2Y
35 POPE C23 36 POPE H3R
35 POPE C23 37 POPE H3S
38 POPE C24 39 POPE H4R
38 POPE C24 40 POPE H4S
41 POPE C25 43 POPE H5S
41 POPE C25 42 POPE H5R
44 POPE C26 45 POPE H6R
44 POPE C26 46 POPE H6S
47 POPE C27 49 POPE H7S
47 POPE C27 48 POPE H7R
50 POPE C28 52 POPE H8S
50 POPE C28 51 POPE H8R
53 POPE C29 54 POPE H91
55 POPE C210 56 POPE H101
57 POPE C211 58 POPE H11R
57 POPE C211 59 POPE H11S
60 POPE C212 62 POPE H12S
60 POPE C212 61 POPE H12R
63 POPE C213 65 POPE H13S
63 POPE C213 64 POPE H13R
66 POPE C214 67 POPE H14R
66 POPE C214 68 POPE H14S
69 POPE C215 70 POPE H15R
69 POPE C215 71 POPE H15S
72 POPE C216 74 POPE H16S
72 POPE C216 73 POPE H16R
75 POPE C217 77 POPE H17S
75 POPE C217 76 POPE H17R
78 POPE C218 81 POPE H18T
78 POPE C218 80 POPE H18S
78 POPE C218 79 POPE H18R
82 POPE C33 83 POPE H3X
82 POPE C33 84 POPE H3Y
85 POPE C34 86 POPE H4X
85 POPE C34 87 POPE H4Y
88 POPE C35 90 POPE H5Y
88 POPE C35 89 POPE H5X
91 POPE C36 93 POPE H6Y
91 POPE C36 92 POPE H6X
94 POPE C37 95 POPE H7X
94 POPE C37 96 POPE H7Y
97 POPE C38 99 POPE H8Y
97 POPE C38 98 POPE H8X
100 POPE C39 101 POPE H9X
100 POPE C39 102 POPE H9Y
103 POPE C310 105 POPE H10Y
103 POPE C310 104 POPE H10X
106 POPE C311 108 POPE H11Y
106 POPE C311 107 POPE H11X
109 POPE C312 111 POPE H12Y
109 POPE C312 110 POPE H12X
112 POPE C313 113 POPE H13X
112 POPE C313 114 POPE H13Y
115 POPE C314 117 POPE H14Y
115 POPE C314 116 POPE H14X
118 POPE C315 120 POPE H15Y
118 POPE C315 119 POPE H15X
121 POPE C316 123 POPE H16Y
121 POPE C316 122 POPE H16X
121 POPE C316 124 POPE H16Z";

        string2bonds(string_representation)
    }

    fn pope_order_atoms() -> OrderAtoms {
        let string_representation = "4 POPE C12
7 POPE C11
15 POPE C1
18 POPE C2
21 POPE C21
23 POPE C22
26 POPE C3
30 POPE C31
32 POPE C32
35 POPE C23
38 POPE C24
41 POPE C25
44 POPE C26
47 POPE C27
50 POPE C28
53 POPE C29
55 POPE C210
57 POPE C211
60 POPE C212
63 POPE C213
66 POPE C214
69 POPE C215
72 POPE C216
75 POPE C217
78 POPE C218
82 POPE C33
85 POPE C34
88 POPE C35
91 POPE C36
94 POPE C37
97 POPE C38
100 POPE C39
103 POPE C310
106 POPE C311
109 POPE C312
112 POPE C313
115 POPE C314
118 POPE C315
121 POPE C316";

        string2atoms(string_representation)
    }

    fn popc_topology() -> MoleculeTopology {
        let string_representation = "0 POPC N 4 POPC C13
0 POPC N 8 POPC C14
0 POPC N 12 POPC C15
0 POPC N 1 POPC C12
1 POPC C12 2 POPC H12A
1 POPC C12 3 POPC H12B
1 POPC C12 16 POPC C11
4 POPC C13 6 POPC H13B
4 POPC C13 7 POPC H13C
4 POPC C13 5 POPC H13A
8 POPC C14 11 POPC H14C
8 POPC C14 9 POPC H14A
8 POPC C14 10 POPC H14B
12 POPC C15 15 POPC H15C
12 POPC C15 13 POPC H15A
12 POPC C15 14 POPC H15B
16 POPC C11 18 POPC H11B
16 POPC C11 17 POPC H11A
16 POPC C11 22 POPC O12
19 POPC P 23 POPC O11
19 POPC P 22 POPC O12
19 POPC P 20 POPC O13
19 POPC P 21 POPC O14
23 POPC O11 24 POPC C1
24 POPC C1 27 POPC C2
24 POPC C1 25 POPC HA
24 POPC C1 26 POPC HB
27 POPC C2 35 POPC C3
27 POPC C2 29 POPC O21
27 POPC C2 28 POPC HS
29 POPC O21 30 POPC C21
30 POPC C21 32 POPC C22
30 POPC C21 31 POPC O22
32 POPC C22 34 POPC H2S
32 POPC C22 44 POPC C23
32 POPC C22 33 POPC H2R
35 POPC C3 36 POPC HX
35 POPC C3 38 POPC O31
35 POPC C3 37 POPC HY
38 POPC O31 39 POPC C31
39 POPC C31 40 POPC O32
39 POPC C31 41 POPC C32
41 POPC C32 91 POPC C33
41 POPC C32 43 POPC H2Y
41 POPC C32 42 POPC H2X
44 POPC C23 47 POPC C24
44 POPC C23 46 POPC H3S
44 POPC C23 45 POPC H3R
47 POPC C24 50 POPC C25
47 POPC C24 48 POPC H4R
47 POPC C24 49 POPC H4S
50 POPC C25 52 POPC H5S
50 POPC C25 51 POPC H5R
50 POPC C25 53 POPC C26
53 POPC C26 54 POPC H6R
53 POPC C26 55 POPC H6S
53 POPC C26 56 POPC C27
56 POPC C27 59 POPC C28
56 POPC C27 58 POPC H7S
56 POPC C27 57 POPC H7R
59 POPC C28 62 POPC C29
59 POPC C28 60 POPC H8R
59 POPC C28 61 POPC H8S
62 POPC C29 64 POPC C210
62 POPC C29 63 POPC H91
64 POPC C210 65 POPC H101
64 POPC C210 66 POPC C211
66 POPC C211 67 POPC H11R
66 POPC C211 68 POPC H11S
66 POPC C211 69 POPC C212
69 POPC C212 70 POPC H12R
69 POPC C212 71 POPC H12S
69 POPC C212 72 POPC C213
72 POPC C213 75 POPC C214
72 POPC C213 73 POPC H13R
72 POPC C213 74 POPC H13S
75 POPC C214 78 POPC C215
75 POPC C214 76 POPC H14R
75 POPC C214 77 POPC H14S
78 POPC C215 81 POPC C216
78 POPC C215 80 POPC H15S
78 POPC C215 79 POPC H15R
81 POPC C216 84 POPC C217
81 POPC C216 82 POPC H16R
81 POPC C216 83 POPC H16S
84 POPC C217 87 POPC C218
84 POPC C217 85 POPC H17R
84 POPC C217 86 POPC H17S
87 POPC C218 89 POPC H18S
87 POPC C218 90 POPC H18T
87 POPC C218 88 POPC H18R
91 POPC C33 94 POPC C34
91 POPC C33 92 POPC H3X
91 POPC C33 93 POPC H3Y
94 POPC C34 95 POPC H4X
94 POPC C34 97 POPC C35
94 POPC C34 96 POPC H4Y
97 POPC C35 99 POPC H5Y
97 POPC C35 100 POPC C36
97 POPC C35 98 POPC H5X
100 POPC C36 103 POPC C37
100 POPC C36 102 POPC H6Y
100 POPC C36 101 POPC H6X
103 POPC C37 106 POPC C38
103 POPC C37 105 POPC H7Y
103 POPC C37 104 POPC H7X
106 POPC C38 107 POPC H8X
106 POPC C38 109 POPC C39
106 POPC C38 108 POPC H8Y
109 POPC C39 110 POPC H9X
109 POPC C39 111 POPC H9Y
109 POPC C39 112 POPC C310
112 POPC C310 115 POPC C311
112 POPC C310 113 POPC H10X
112 POPC C310 114 POPC H10Y
115 POPC C311 117 POPC H11Y
115 POPC C311 118 POPC C312
115 POPC C311 116 POPC H11X
118 POPC C312 121 POPC C313
118 POPC C312 119 POPC H12X
118 POPC C312 120 POPC H12Y
121 POPC C313 124 POPC C314
121 POPC C313 122 POPC H13X
121 POPC C313 123 POPC H13Y
124 POPC C314 127 POPC C315
124 POPC C314 126 POPC H14Y
124 POPC C314 125 POPC H14X
127 POPC C315 128 POPC H15X
127 POPC C315 129 POPC H15Y
127 POPC C315 130 POPC C316
130 POPC C316 131 POPC H16X
130 POPC C316 133 POPC H16Z
130 POPC C316 132 POPC H16Y";

        string2topology(string_representation)
    }

    fn popc_order_bonds() -> HashSet<BondType> {
        let string_representation = "1 POPC C12 3 POPC H12B
1 POPC C12 2 POPC H12A
4 POPC C13 5 POPC H13A
4 POPC C13 7 POPC H13C
4 POPC C13 6 POPC H13B
8 POPC C14 10 POPC H14B
8 POPC C14 11 POPC H14C
8 POPC C14 9 POPC H14A
12 POPC C15 14 POPC H15B
12 POPC C15 13 POPC H15A
12 POPC C15 15 POPC H15C
16 POPC C11 18 POPC H11B
16 POPC C11 17 POPC H11A
24 POPC C1 26 POPC HB
24 POPC C1 25 POPC HA
27 POPC C2 28 POPC HS
32 POPC C22 33 POPC H2R
32 POPC C22 34 POPC H2S
35 POPC C3 37 POPC HY
35 POPC C3 36 POPC HX
41 POPC C32 42 POPC H2X
41 POPC C32 43 POPC H2Y
44 POPC C23 46 POPC H3S
44 POPC C23 45 POPC H3R
47 POPC C24 49 POPC H4S
47 POPC C24 48 POPC H4R
50 POPC C25 51 POPC H5R
50 POPC C25 52 POPC H5S
53 POPC C26 55 POPC H6S
53 POPC C26 54 POPC H6R
56 POPC C27 57 POPC H7R
56 POPC C27 58 POPC H7S
59 POPC C28 60 POPC H8R
59 POPC C28 61 POPC H8S
62 POPC C29 63 POPC H91
64 POPC C210 65 POPC H101
66 POPC C211 67 POPC H11R
66 POPC C211 68 POPC H11S
69 POPC C212 71 POPC H12S
69 POPC C212 70 POPC H12R
72 POPC C213 74 POPC H13S
72 POPC C213 73 POPC H13R
75 POPC C214 76 POPC H14R
75 POPC C214 77 POPC H14S
78 POPC C215 80 POPC H15S
78 POPC C215 79 POPC H15R
81 POPC C216 83 POPC H16S
81 POPC C216 82 POPC H16R
84 POPC C217 85 POPC H17R
84 POPC C217 86 POPC H17S
87 POPC C218 90 POPC H18T
87 POPC C218 89 POPC H18S
87 POPC C218 88 POPC H18R
91 POPC C33 93 POPC H3Y
91 POPC C33 92 POPC H3X
94 POPC C34 95 POPC H4X
94 POPC C34 96 POPC H4Y
97 POPC C35 98 POPC H5X
97 POPC C35 99 POPC H5Y
100 POPC C36 101 POPC H6X
100 POPC C36 102 POPC H6Y
103 POPC C37 104 POPC H7X
103 POPC C37 105 POPC H7Y
106 POPC C38 107 POPC H8X
106 POPC C38 108 POPC H8Y
109 POPC C39 110 POPC H9X
109 POPC C39 111 POPC H9Y
112 POPC C310 114 POPC H10Y
112 POPC C310 113 POPC H10X
115 POPC C311 116 POPC H11X
115 POPC C311 117 POPC H11Y
118 POPC C312 120 POPC H12Y
118 POPC C312 119 POPC H12X
121 POPC C313 122 POPC H13X
121 POPC C313 123 POPC H13Y
124 POPC C314 125 POPC H14X
124 POPC C314 126 POPC H14Y
127 POPC C315 128 POPC H15X
127 POPC C315 129 POPC H15Y
130 POPC C316 133 POPC H16Z
130 POPC C316 131 POPC H16X
130 POPC C316 132 POPC H16Y";

        string2bonds(string_representation)
    }

    fn popc_order_atoms() -> OrderAtoms {
        let string_representation = "1 POPC C12
4 POPC C13
8 POPC C14
12 POPC C15
16 POPC C11
24 POPC C1
27 POPC C2
30 POPC C21
32 POPC C22
35 POPC C3
39 POPC C31
41 POPC C32
44 POPC C23
47 POPC C24
50 POPC C25
53 POPC C26
56 POPC C27
59 POPC C28
62 POPC C29
64 POPC C210
66 POPC C211
69 POPC C212
72 POPC C213
75 POPC C214
78 POPC C215
81 POPC C216
84 POPC C217
87 POPC C218
91 POPC C33
94 POPC C34
97 POPC C35
100 POPC C36
103 POPC C37
106 POPC C38
109 POPC C39
112 POPC C310
115 POPC C311
118 POPC C312
121 POPC C313
124 POPC C314
127 POPC C315
130 POPC C316";

        string2atoms(string_representation)
    }

    fn popg_topology() -> MoleculeTopology {
        let string_representation = "0 POPG C13 3 POPG OC3
0 POPG C13 5 POPG C12
0 POPG C13 2 POPG H13B
0 POPG C13 1 POPG H13A
3 POPG OC3 4 POPG HO3
5 POPG C12 6 POPG H12A
5 POPG C12 7 POPG OC2
5 POPG C12 9 POPG C11
7 POPG OC2 8 POPG HO2
9 POPG C11 11 POPG H11B
9 POPG C11 10 POPG H11A
9 POPG C11 15 POPG O12
12 POPG P 14 POPG O14
12 POPG P 16 POPG O11
12 POPG P 13 POPG O13
12 POPG P 15 POPG O12
16 POPG O11 17 POPG C1
17 POPG C1 20 POPG C2
17 POPG C1 19 POPG HB
17 POPG C1 18 POPG HA
20 POPG C2 22 POPG O21
20 POPG C2 21 POPG HS
20 POPG C2 28 POPG C3
22 POPG O21 23 POPG C21
23 POPG C21 24 POPG O22
23 POPG C21 25 POPG C22
25 POPG C22 37 POPG C23
25 POPG C22 26 POPG H2R
25 POPG C22 27 POPG H2S
28 POPG C3 29 POPG HX
28 POPG C3 31 POPG O31
28 POPG C3 30 POPG HY
31 POPG O31 32 POPG C31
32 POPG C31 33 POPG O32
32 POPG C31 34 POPG C32
34 POPG C32 84 POPG C33
34 POPG C32 36 POPG H2Y
34 POPG C32 35 POPG H2X
37 POPG C23 39 POPG H3S
37 POPG C23 40 POPG C24
37 POPG C23 38 POPG H3R
40 POPG C24 43 POPG C25
40 POPG C24 41 POPG H4R
40 POPG C24 42 POPG H4S
43 POPG C25 46 POPG C26
43 POPG C25 45 POPG H5S
43 POPG C25 44 POPG H5R
46 POPG C26 49 POPG C27
46 POPG C26 47 POPG H6R
46 POPG C26 48 POPG H6S
49 POPG C27 52 POPG C28
49 POPG C27 51 POPG H7S
49 POPG C27 50 POPG H7R
52 POPG C28 54 POPG H8S
52 POPG C28 53 POPG H8R
52 POPG C28 55 POPG C29
55 POPG C29 56 POPG H91
55 POPG C29 57 POPG C210
57 POPG C210 59 POPG C211
57 POPG C210 58 POPG H101
59 POPG C211 62 POPG C212
59 POPG C211 61 POPG H11S
59 POPG C211 60 POPG H11R
62 POPG C212 63 POPG H12R
62 POPG C212 64 POPG H12S
62 POPG C212 65 POPG C213
65 POPG C213 67 POPG H13S
65 POPG C213 68 POPG C214
65 POPG C213 66 POPG H13R
68 POPG C214 69 POPG H14R
68 POPG C214 70 POPG H14S
68 POPG C214 71 POPG C215
71 POPG C215 72 POPG H15R
71 POPG C215 74 POPG C216
71 POPG C215 73 POPG H15S
74 POPG C216 77 POPG C217
74 POPG C216 75 POPG H16R
74 POPG C216 76 POPG H16S
77 POPG C217 79 POPG H17S
77 POPG C217 80 POPG C218
77 POPG C217 78 POPG H17R
80 POPG C218 82 POPG H18S
80 POPG C218 81 POPG H18R
80 POPG C218 83 POPG H18T
84 POPG C33 85 POPG H3X
84 POPG C33 86 POPG H3Y
84 POPG C33 87 POPG C34
87 POPG C34 88 POPG H4X
87 POPG C34 90 POPG C35
87 POPG C34 89 POPG H4Y
90 POPG C35 93 POPG C36
90 POPG C35 92 POPG H5Y
90 POPG C35 91 POPG H5X
93 POPG C36 96 POPG C37
93 POPG C36 95 POPG H6Y
93 POPG C36 94 POPG H6X
96 POPG C37 97 POPG H7X
96 POPG C37 98 POPG H7Y
96 POPG C37 99 POPG C38
99 POPG C38 101 POPG H8Y
99 POPG C38 102 POPG C39
99 POPG C38 100 POPG H8X
102 POPG C39 105 POPG C310
102 POPG C39 104 POPG H9Y
102 POPG C39 103 POPG H9X
105 POPG C310 108 POPG C311
105 POPG C310 107 POPG H10Y
105 POPG C310 106 POPG H10X
108 POPG C311 110 POPG H11Y
108 POPG C311 109 POPG H11X
108 POPG C311 111 POPG C312
111 POPG C312 113 POPG H12Y
111 POPG C312 112 POPG H12X
111 POPG C312 114 POPG C313
114 POPG C313 116 POPG H13Y
114 POPG C313 115 POPG H13X
114 POPG C313 117 POPG C314
117 POPG C314 118 POPG H14X
117 POPG C314 120 POPG C315
117 POPG C314 119 POPG H14Y
120 POPG C315 123 POPG C316
120 POPG C315 122 POPG H15Y
120 POPG C315 121 POPG H15X
123 POPG C316 124 POPG H16X
123 POPG C316 125 POPG H16Y
123 POPG C316 126 POPG H16Z";

        string2topology(string_representation)
    }

    fn popg_order_bonds() -> HashSet<BondType> {
        let string_representation = "0 POPG C13 2 POPG H13B
0 POPG C13 1 POPG H13A
5 POPG C12 6 POPG H12A
9 POPG C11 10 POPG H11A
9 POPG C11 11 POPG H11B
17 POPG C1 19 POPG HB
17 POPG C1 18 POPG HA
20 POPG C2 21 POPG HS
25 POPG C22 27 POPG H2S
25 POPG C22 26 POPG H2R
28 POPG C3 29 POPG HX
28 POPG C3 30 POPG HY
34 POPG C32 35 POPG H2X
34 POPG C32 36 POPG H2Y
37 POPG C23 39 POPG H3S
37 POPG C23 38 POPG H3R
40 POPG C24 41 POPG H4R
40 POPG C24 42 POPG H4S
43 POPG C25 45 POPG H5S
43 POPG C25 44 POPG H5R
46 POPG C26 48 POPG H6S
46 POPG C26 47 POPG H6R
49 POPG C27 50 POPG H7R
49 POPG C27 51 POPG H7S
52 POPG C28 54 POPG H8S
52 POPG C28 53 POPG H8R
55 POPG C29 56 POPG H91
57 POPG C210 58 POPG H101
59 POPG C211 60 POPG H11R
59 POPG C211 61 POPG H11S
62 POPG C212 63 POPG H12R
62 POPG C212 64 POPG H12S
65 POPG C213 66 POPG H13R
65 POPG C213 67 POPG H13S
68 POPG C214 69 POPG H14R
68 POPG C214 70 POPG H14S
71 POPG C215 73 POPG H15S
71 POPG C215 72 POPG H15R
74 POPG C216 75 POPG H16R
74 POPG C216 76 POPG H16S
77 POPG C217 79 POPG H17S
77 POPG C217 78 POPG H17R
80 POPG C218 83 POPG H18T
80 POPG C218 82 POPG H18S
80 POPG C218 81 POPG H18R
84 POPG C33 86 POPG H3Y
84 POPG C33 85 POPG H3X
87 POPG C34 88 POPG H4X
87 POPG C34 89 POPG H4Y
90 POPG C35 91 POPG H5X
90 POPG C35 92 POPG H5Y
93 POPG C36 94 POPG H6X
93 POPG C36 95 POPG H6Y
96 POPG C37 97 POPG H7X
96 POPG C37 98 POPG H7Y
99 POPG C38 100 POPG H8X
99 POPG C38 101 POPG H8Y
102 POPG C39 103 POPG H9X
102 POPG C39 104 POPG H9Y
105 POPG C310 106 POPG H10X
105 POPG C310 107 POPG H10Y
108 POPG C311 109 POPG H11X
108 POPG C311 110 POPG H11Y
111 POPG C312 113 POPG H12Y
111 POPG C312 112 POPG H12X
114 POPG C313 115 POPG H13X
114 POPG C313 116 POPG H13Y
117 POPG C314 119 POPG H14Y
117 POPG C314 118 POPG H14X
120 POPG C315 121 POPG H15X
120 POPG C315 122 POPG H15Y
123 POPG C316 124 POPG H16X
123 POPG C316 126 POPG H16Z
123 POPG C316 125 POPG H16Y";

        string2bonds(string_representation)
    }

    fn popg_order_atoms() -> OrderAtoms {
        let string_representation = "0 POPG C13
5 POPG C12
9 POPG C11
17 POPG C1
20 POPG C2
23 POPG C21
25 POPG C22
28 POPG C3
32 POPG C31
34 POPG C32
37 POPG C23
40 POPG C24
43 POPG C25
46 POPG C26
49 POPG C27
52 POPG C28
55 POPG C29
57 POPG C210
59 POPG C211
62 POPG C212
65 POPG C213
68 POPG C214
71 POPG C215
74 POPG C216
77 POPG C217
80 POPG C218
84 POPG C33
87 POPG C34
90 POPG C35
93 POPG C36
96 POPG C37
99 POPG C38
102 POPG C39
105 POPG C310
108 POPG C311
111 POPG C312
114 POPG C313
117 POPG C314
120 POPG C315
123 POPG C316";

        string2atoms(string_representation)
    }

    fn string2topology(string: &str) -> MoleculeTopology {
        let bonds = string2bonds(string);
        MoleculeTopology { bonds }
    }

    fn string2bonds(string: &str) -> HashSet<BondType> {
        let string_representation = string.trim().split("\n");

        let mut bonds = HashSet::new();
        for bond in string_representation.into_iter() {
            let split: Vec<&str> = bond.split_whitespace().into_iter().collect();
            let index1 = split[0].parse::<usize>().unwrap();
            let resname1 = split[1];
            let atomname1 = split[2];

            let index2 = split[3].parse::<usize>().unwrap();
            let resname2 = split[4];
            let atomname2 = split[5];

            let atom1 = AtomType {
                relative_index: index1,
                residue_name: resname1.to_owned(),
                atom_name: atomname1.to_owned(),
            };

            let atom2 = AtomType {
                relative_index: index2,
                residue_name: resname2.to_owned(),
                atom_name: atomname2.to_owned(),
            };

            let bond = BondType { atom1, atom2 };
            bonds.insert(bond);
        }

        bonds
    }

    fn string2atoms(string: &str) -> OrderAtoms {
        let string_representation = string.trim().split("\n");

        let mut atoms = Vec::new();
        for atom in string_representation.into_iter() {
            let split: Vec<&str> = atom.split_whitespace().into_iter().collect();
            let relative_index = split[0].parse::<usize>().unwrap();
            let residue_name = split[1];
            let atom_name = split[2];

            atoms.push(AtomType {
                relative_index,
                residue_name: residue_name.to_owned(),
                atom_name: atom_name.to_owned(),
            })
        }

        OrderAtoms { atoms }
    }

    #[test]
    fn test_classify_molecules() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();

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
        let expected_names = ["POPE", "POPC", "POPG"];
        let expected_topology = [pope_topology(), popc_topology(), popg_topology()];
        let expected_order_bond_types =
            [pope_order_bonds(), popc_order_bonds(), popg_order_bonds()];
        let expected_n_instances = [131, 128, 15];
        let expected_bond_instances = [
            vec![
                (7, 8),
                (132, 133),
                (257, 258),
                (382, 383),
                (507, 508),
                (632, 633),
                (757, 758),
                (882, 883),
                (1007, 1008),
                (1132, 1133),
                (1257, 1258),
                (1382, 1383),
                (1507, 1508),
                (1632, 1633),
                (1757, 1758),
                (1882, 1883),
                (2007, 2008),
                (2132, 2133),
                (2257, 2258),
                (2382, 2383),
                (2507, 2508),
                (2632, 2633),
                (2757, 2758),
                (2882, 2883),
                (3007, 3008),
                (3132, 3133),
                (3257, 3258),
                (3382, 3383),
                (3507, 3508),
                (3632, 3633),
                (3757, 3758),
                (3882, 3883),
                (4007, 4008),
                (4132, 4133),
                (4257, 4258),
                (4382, 4383),
                (4507, 4508),
                (4632, 4633),
                (4757, 4758),
                (4882, 4883),
                (5007, 5008),
                (5132, 5133),
                (5257, 5258),
                (5382, 5383),
                (5507, 5508),
                (5632, 5633),
                (5757, 5758),
                (5882, 5883),
                (6007, 6008),
                (6132, 6133),
                (6257, 6258),
                (6382, 6383),
                (6507, 6508),
                (6632, 6633),
                (6757, 6758),
                (6882, 6883),
                (7007, 7008),
                (7132, 7133),
                (7257, 7258),
                (7382, 7383),
                (7507, 7508),
                (7632, 7633),
                (7757, 7758),
                (7882, 7883),
                (8007, 8008),
                (8132, 8133),
                (8257, 8258),
                (8382, 8383),
                (8507, 8508),
                (8632, 8633),
                (8757, 8758),
                (8882, 8883),
                (9007, 9008),
                (9132, 9133),
                (9257, 9258),
                (9382, 9383),
                (9507, 9508),
                (9632, 9633),
                (9757, 9758),
                (9882, 9883),
                (10007, 10008),
                (10132, 10133),
                (10257, 10258),
                (10382, 10383),
                (10507, 10508),
                (10632, 10633),
                (10757, 10758),
                (10882, 10883),
                (11007, 11008),
                (11132, 11133),
                (11257, 11258),
                (11382, 11383),
                (11507, 11508),
                (11632, 11633),
                (11757, 11758),
                (11882, 11883),
                (12007, 12008),
                (12132, 12133),
                (12257, 12258),
                (12382, 12383),
                (12507, 12508),
                (12632, 12633),
                (12757, 12758),
                (12882, 12883),
                (13007, 13008),
                (13132, 13133),
                (13257, 13258),
                (13382, 13383),
                (13507, 13508),
                (13632, 13633),
                (13757, 13758),
                (13882, 13883),
                (14007, 14008),
                (14132, 14133),
                (14257, 14258),
                (14382, 14383),
                (14507, 14508),
                (14632, 14633),
                (14757, 14758),
                (14882, 14883),
                (15007, 15008),
                (15132, 15133),
                (15257, 15258),
                (15382, 15383),
                (15507, 15508),
                (15632, 15633),
                (15757, 15758),
                (15882, 15883),
                (16007, 16008),
                (16132, 16133),
                (16257, 16258),
            ],
            vec![
                (16422, 16424),
                (16556, 16558),
                (16690, 16692),
                (16824, 16826),
                (16958, 16960),
                (17092, 17094),
                (17226, 17228),
                (17360, 17362),
                (17494, 17496),
                (17628, 17630),
                (17762, 17764),
                (17896, 17898),
                (18030, 18032),
                (18164, 18166),
                (18298, 18300),
                (18432, 18434),
                (18566, 18568),
                (18700, 18702),
                (18834, 18836),
                (18968, 18970),
                (19102, 19104),
                (19236, 19238),
                (19370, 19372),
                (19504, 19506),
                (19638, 19640),
                (19772, 19774),
                (19906, 19908),
                (20040, 20042),
                (20174, 20176),
                (20308, 20310),
                (20442, 20444),
                (20576, 20578),
                (20710, 20712),
                (20844, 20846),
                (20978, 20980),
                (21112, 21114),
                (21246, 21248),
                (21380, 21382),
                (21514, 21516),
                (21648, 21650),
                (21782, 21784),
                (21916, 21918),
                (22050, 22052),
                (22184, 22186),
                (22318, 22320),
                (22452, 22454),
                (22586, 22588),
                (22720, 22722),
                (22854, 22856),
                (22988, 22990),
                (23122, 23124),
                (23256, 23258),
                (23390, 23392),
                (23524, 23526),
                (23658, 23660),
                (23792, 23794),
                (23926, 23928),
                (24060, 24062),
                (24194, 24196),
                (24328, 24330),
                (24462, 24464),
                (24596, 24598),
                (24730, 24732),
                (24864, 24866),
                (24998, 25000),
                (25132, 25134),
                (25266, 25268),
                (25400, 25402),
                (25534, 25536),
                (25668, 25670),
                (25802, 25804),
                (25936, 25938),
                (26070, 26072),
                (26204, 26206),
                (26338, 26340),
                (26472, 26474),
                (26606, 26608),
                (26740, 26742),
                (26874, 26876),
                (27008, 27010),
                (27142, 27144),
                (27276, 27278),
                (27410, 27412),
                (27544, 27546),
                (27678, 27680),
                (27812, 27814),
                (27946, 27948),
                (28080, 28082),
                (28214, 28216),
                (28348, 28350),
                (28482, 28484),
                (28616, 28618),
                (28750, 28752),
                (28884, 28886),
                (29018, 29020),
                (29152, 29154),
                (29286, 29288),
                (29420, 29422),
                (29554, 29556),
                (29688, 29690),
                (29822, 29824),
                (29956, 29958),
                (30090, 30092),
                (30224, 30226),
                (30358, 30360),
                (30492, 30494),
                (30626, 30628),
                (30760, 30762),
                (30894, 30896),
                (31028, 31030),
                (31162, 31164),
                (31296, 31298),
                (31430, 31432),
                (31564, 31566),
                (31698, 31700),
                (31832, 31834),
                (31966, 31968),
                (32100, 32102),
                (32234, 32236),
                (32368, 32370),
                (32502, 32504),
                (32636, 32638),
                (32770, 32772),
                (32904, 32906),
                (33038, 33040),
                (33172, 33174),
                (33306, 33308),
                (33440, 33442),
            ],
            vec![
                (33620, 33621),
                (33747, 33748),
                (33874, 33875),
                (34001, 34002),
                (34128, 34129),
                (34255, 34256),
                (34382, 34383),
                (34509, 34510),
                (34636, 34637),
                (34763, 34764),
                (34890, 34891),
                (35017, 35018),
                (35144, 35145),
                (35271, 35272),
                (35398, 35399),
            ],
        ];
        let expected_order_atoms = [pope_order_atoms(), popc_order_atoms(), popg_order_atoms()];

        for (i, molecule) in molecules.into_iter().enumerate() {
            assert_eq!(molecule.name(), expected_names[i]);
            assert_eq!(molecule.topology(), &expected_topology[i]);
            let order_bonds = molecule.order_bonds();
            let order_bond_types: HashSet<BondType> = order_bonds
                .bonds
                .iter()
                .map(|x| x.bond_type.clone())
                .collect();

            assert_eq!(order_bond_types, expected_order_bond_types[i]);

            for (b, bond_instances) in order_bonds
                .bonds
                .iter()
                .map(|x| x.bonds.clone())
                .enumerate()
            {
                assert_eq!(bond_instances.len(), expected_n_instances[i]);
                if b == 0 {
                    assert_eq!(bond_instances, expected_bond_instances[i]);
                }
            }

            assert_eq!(molecule.order_atoms(), &expected_order_atoms[i]);
        }
    }

    #[test]
    fn test_classify_molecules_cg() {
        todo!("Test classification of molecules for coarse grained order parameters calculation.")
    }

    #[test]
    fn test_classify_molecules_shared_name() {
        todo!("Test classification of molecules composed of the residue of the same name.")
    }

    #[test]
    fn test_classify_molecules_multiple_residues() {
        todo!("Test classification of molecules composed of multiple residues.")
    }
}
