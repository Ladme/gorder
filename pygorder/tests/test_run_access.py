"""
Released under MIT License.
Copyright (c) 2024-2025 Ladislav Bartos
"""

import gorder, pytest, math

def compare_orders(x: float, y: float) -> bool:
    return math.isclose(round(x, 4), round(y, 4), rel_tol = 1e-4)

def test_aa_order_basic():
    analysis = gorder.Analysis(
        structure = "../tests/files/pcpepg.tpr",
        trajectory = "../tests/files/pcpepg.xtc",
        analysis_type = gorder.analysis_types.AAOrder("@membrane and element name carbon", "@membrane and element name hydrogen"),
        silent = True,
        overwrite = True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 51
    assert len(results.molecules()) == 3
    
    assert compare_orders(results.average_order().total().value(), 0.1423)
    assert results.average_order().total().error() is None
    assert results.average_order().upper() is None
    assert results.average_order().lower() is None
    
    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    with pytest.raises(gorder.exceptions.APIError) as excinfo:
        results.get_molecule("POPA")
    assert "molecule with the given name does not exist" in str(excinfo.value)

    expected_average_orders = [0.1455, 0.1378, 0.1561]
    expected_atom_numbers = [37, 40, 38]
    expected_molecule_names = ["POPE", "POPC", "POPG"]

    expected_atom_indices = [32, 41, 34]
    expected_atom_names = ["C32", "C32", "C32"]
    expected_atom_order = [0.2226, 0.2363, 0.2247]

    expected_bond_numbers = [2, 2, 2]

    expected_atom2_indices = [34, 43, 36]
    expected_atom2_names = ["H2Y", "H2Y", "H2Y"]
    expected_atom2_order = [0.2040, 0.2317, 0.2020]

    for (i, molecule) in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        average_order = molecule.average_order()
        assert compare_orders(average_order.total().value(), expected_average_orders[i])
        assert average_order.total().error() is None
        assert average_order.upper() is None
        assert average_order.lower() is None

        average_maps = molecule.average_ordermaps()
        assert average_maps.total() is None
        assert average_maps.upper() is None
        assert average_maps.lower() is None

        # ATOMS
        assert len(molecule.atoms()) == expected_atom_numbers[i]

        atom = molecule.get_atom(expected_atom_indices[i])
        atom_type = atom.atom()
        assert atom_type.atom_name() == expected_atom_names[i]
        assert atom_type.relative_index() == expected_atom_indices[i]
        assert atom_type.residue_name() == expected_molecule_names[i]
        assert atom.molecule() == expected_molecule_names[i]

        order = atom.order()
        assert compare_orders(order.total().value(), expected_atom_order[i])
        assert order.total().error() is None
        assert order.upper() is None
        assert order.lower() is None

        maps = atom.ordermaps()
        assert maps.total() is None
        assert maps.upper() is None
        assert maps.lower() is None

        # BONDS
        assert len(atom.bonds()) == expected_bond_numbers[i]

        bond = atom.get_bond(expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.atom_name() == expected_atom_names[i]
        assert a1.relative_index() == expected_atom_indices[i]
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == expected_atom2_names[i]
        assert a2.relative_index() == expected_atom2_indices[i]
        assert a2.residue_name() == expected_molecule_names[i]
        assert bond.molecule() == expected_molecule_names[i]

        order = bond.order()
        assert compare_orders(order.total().value(), expected_atom2_order[i])
        assert order.total().error() is None
        assert order.upper() is None
        assert order.lower() is None

        maps = bond.ordermaps()
        assert maps.total() is None
        assert maps.upper() is None
        assert maps.lower() is None

        # BOND FROM MOLECULE
        bond = molecule.get_bond(expected_atom_indices[i], expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.relative_index() == expected_atom_indices[i]
        assert a2.relative_index() == expected_atom2_indices[i]

        bond = molecule.get_bond(expected_atom2_indices[i], expected_atom_indices[i])
        a1, a2 = bond.atoms()
        assert a1.relative_index() == expected_atom_indices[i]
        assert a2.relative_index() == expected_atom2_indices[i]

        # NONEXISTENT ATOM
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_atom(145)
        assert "atom with the given relative index does not exist" in str(excinfo.value)

        # NONEXISTENT BOND
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(7, 19)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(145, 189)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_cg_order_basic():
    analysis = gorder.Analysis(
        structure = "../tests/files/cg.tpr",
        trajectory = "../tests/files/cg.xtc",
        analysis_type = gorder.analysis_types.CGOrder("@membrane"),
        silent = True,
        overwrite = True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 101
    assert len(results.molecules()) == 3
    
    assert compare_orders(results.average_order().total().value(), 0.2962)
    assert results.average_order().total().error() is None
    assert results.average_order().upper() is None
    assert results.average_order().lower() is None
    
    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    with pytest.raises(gorder.exceptions.APIError) as excinfo:
        results.get_molecule("POPA")
    assert "molecule with the given name does not exist" in str(excinfo.value)

    expected_molecule_names = ["POPC", "POPE", "POPG"]
    expected_average_orders = [0.2943, 0.2972, 0.3059]
    expected_bond_orders = [0.3682, 0.3759, 0.3789]

    for (i, molecule) in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        average_order = molecule.average_order()
        assert compare_orders(average_order.total().value(), expected_average_orders[i])
        assert average_order.total().error() is None
        assert average_order.upper() is None
        assert average_order.lower() is None

        average_maps = molecule.average_ordermaps()
        assert average_maps.total() is None
        assert average_maps.upper() is None
        assert average_maps.lower() is None

        # BONDS
        assert len(molecule.bonds()) == 11

        bond = molecule.get_bond(4, 5)
        a1, a2 = bond.atoms()
        assert a1.atom_name() == "C1A"
        assert a1.relative_index() == 4
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == "D2A"
        assert a2.relative_index() == 5
        assert a2.residue_name() == expected_molecule_names[i]

        order = bond.order()
        assert compare_orders(order.total().value(), expected_bond_orders[i])
        assert order.total().error() is None
        assert order.upper() is None
        assert order.lower() is None

        maps = bond.ordermaps()
        assert maps.total() is None
        assert maps.upper() is None
        assert maps.lower() is None

        # THE SAME BOND
        bond = molecule.get_bond(5, 4)
        a1, a2 = bond.atoms()
        assert a1.relative_index() == 4
        assert a2.relative_index() == 5

        # NONEXISTENT BOND
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(1, 3)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(15, 16)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

        # ATTEMPTING TO ACCESS ATOMS
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.atoms()
        assert "results for individual atoms are not available for coarse-grained order parameters" in str(excinfo.value)

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_atom(3)
        assert "results for individual atoms are not available for coarse-grained order parameters" in str(excinfo.value)

def test_aa_order_error():
    analysis = gorder.Analysis(
        structure = "../tests/files/pcpepg.tpr",
        trajectory = "../tests/files/pcpepg.xtc",
        analysis_type = gorder.analysis_types.AAOrder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen"
        ),
        estimate_error = gorder.estimate_error.EstimateError(),
        silent = True,
        overwrite = True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 51
    assert len(results.molecules()) == 3

    avg_order = results.average_order()
    assert compare_orders(avg_order.total().value(), 0.1423)
    assert compare_orders(avg_order.total().error(), 0.0026)
    assert avg_order.upper() is None
    assert avg_order.lower() is None

    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    expected_average_orders = [0.1455, 0.1378, 0.1561]
    expected_average_errors = [0.0029, 0.0036, 0.0112]
    expected_atom_numbers = [37, 40, 38]
    expected_molecule_names = ["POPE", "POPC", "POPG"]

    expected_atom_indices = [32, 41, 34]
    expected_atom_names = ["C32", "C32", "C32"]
    expected_atom_order = [0.2226, 0.2363, 0.2247]
    expected_atom_errors = [0.0087, 0.0071, 0.0574]

    expected_bond_numbers = [2, 2, 2]

    expected_atom2_indices = [34, 43, 36]
    expected_atom2_names = ["H2Y", "H2Y", "H2Y"]
    expected_atom2_order = [0.2040, 0.2317, 0.2020]
    expected_atom2_errors = [0.0125, 0.0091, 0.0656]

    expected_convergence_frames = list(range(1, 52))
    expected_convergence_values = [
        [0.1494, 0.1460, 0.1455],
        [0.1422, 0.1353, 0.1378],
        [0.1572, 0.1507, 0.1561],
    ]

    for i, molecule in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        avg_order = molecule.average_order()
        assert compare_orders(avg_order.total().value(), expected_average_orders[i])
        assert compare_orders(avg_order.total().error(), expected_average_errors[i])
        assert avg_order.upper() is None
        assert avg_order.lower() is None

        # CONVERGENCE
        convergence = molecule.convergence()
        assert convergence.frames() == expected_convergence_frames
        
        sample_frames = [0, 25, 50]
        extracted_conv = convergence.total()
        for j, frame in enumerate(sample_frames):
            conv_value = extracted_conv[frame]
            expected = expected_convergence_values[i][j]
            assert compare_orders(conv_value, expected)

        assert convergence.upper() is None
        assert convergence.lower() is None

        # ATOMS
        assert len(molecule.atoms()) == expected_atom_numbers[i]
        atom = molecule.get_atom(expected_atom_indices[i])
        atom_type = atom.atom()
        assert atom_type.atom_name() == expected_atom_names[i]
        assert atom_type.relative_index() == expected_atom_indices[i]
        assert atom_type.residue_name() == expected_molecule_names[i]
        assert atom.molecule() == expected_molecule_names[i]

        order = atom.order()
        assert compare_orders(order.total().value(), expected_atom_order[i])
        assert compare_orders(order.total().error(), expected_atom_errors[i])
        assert order.upper() is None
        assert order.lower() is None

        maps = atom.ordermaps()
        assert maps.total() is None
        assert maps.upper() is None
        assert maps.lower() is None

        # BONDS
        assert len(atom.bonds()) == expected_bond_numbers[i]
        bond = atom.get_bond(expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.atom_name() == expected_atom_names[i]
        assert a1.relative_index() == expected_atom_indices[i]
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == expected_atom2_names[i]
        assert a2.relative_index() == expected_atom2_indices[i]
        assert a2.residue_name() == expected_molecule_names[i]
        assert bond.molecule() == expected_molecule_names[i]
        
        bond_order = bond.order()
        assert compare_orders(bond_order.total().value(), expected_atom2_order[i])
        assert compare_orders(bond_order.total().error(), expected_atom2_errors[i])
        assert bond_order.upper() is None
        assert bond_order.lower() is None

        bond_maps = bond.ordermaps()
        assert bond_maps.total() is None
        assert bond_maps.upper() is None
        assert bond_maps.lower() is None

        # BOND FROM MOLECULE
        bond = molecule.get_bond(expected_atom_indices[i], expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.relative_index() == expected_atom_indices[i]
        assert a2.relative_index() == expected_atom2_indices[i]

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_atom(145)
        assert "atom with the given relative index does not exist" in str(excinfo.value)
            
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(7, 19)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_cg_order_error_leaflets():
    analysis = gorder.Analysis(
        structure="../tests/files/cg.tpr",
        trajectory="../tests/files/cg.xtc",
        analysis_type=gorder.analysis_types.CGOrder("@membrane"),
        estimate_error=gorder.estimate_error.EstimateError(),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 101
    assert len(results.molecules()) == 3

    avg_order = results.average_order()
    assert compare_orders(avg_order.total().value(), 0.2962)
    assert compare_orders(avg_order.total().error(), 0.0050)
    assert avg_order.upper() is None
    assert avg_order.lower() is None

    # Test ordermaps absence
    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    expected_molecule_names = ["POPC", "POPE", "POPG"]
    expected_average_orders = [0.2943, 0.2972, 0.3059]
    expected_average_errors = [0.0067, 0.0052, 0.0089]
    expected_bond_orders = [0.3682, 0.3759, 0.3789]
    expected_bond_errors = [0.0125, 0.0164, 0.0159]
    expected_convergence_values = [
        [0.2756, 0.2902, 0.2943],
        [0.2830, 0.2995, 0.2972],
        [0.3198, 0.3066, 0.3059],
    ]

    for i, molecule in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        # AVERAGE
        avg_order = molecule.average_order()
        assert compare_orders(avg_order.total().value(), expected_average_orders[i])
        assert compare_orders(avg_order.total().error(), expected_average_errors[i])
        assert avg_order.upper() is None
        assert avg_order.lower() is None

        avg_maps = molecule.average_ordermaps()
        assert avg_maps.total() is None
        assert avg_maps.upper() is None
        assert avg_maps.lower() is None

        # CONVERGENCE
        convergence = molecule.convergence()
        assert len(convergence.frames()) == 101
        assert convergence.upper() is None
        assert convergence.lower() is None

        # Check specific convergence points
        sample_indices = [0, 50, 100]
        for j, idx in enumerate(sample_indices):
            conv_value = convergence.total()[idx]
            expected = expected_convergence_values[i][j]
            assert compare_orders(conv_value, expected)

        # BONDS
        assert len(molecule.bonds()) == 11
        bond = molecule.get_bond(4, 5)
        a1, a2 = bond.atoms()
        assert a1.atom_name() == "C1A"
        assert a1.relative_index() == 4
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == "D2A"
        assert a2.relative_index() == 5
        assert a2.residue_name() == expected_molecule_names[i]

        order = bond.order()
        assert compare_orders(order.total().value(), expected_bond_orders[i])
        assert compare_orders(order.total().error(), expected_bond_errors[i])
        assert order.upper() is None
        assert order.lower() is None

        reverse_bond = molecule.get_bond(5, 4)
        a1_rev, a2_rev = reverse_bond.atoms()
        assert a1_rev.relative_index() == 4
        assert a2_rev.relative_index() == 5

        # NONEXISTENT BONDS
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(1, 3)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(15, 16)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_aa_order_leaflets():
    analysis = gorder.Analysis(
        structure = "../tests/files/pcpepg.tpr",
        trajectory = "../tests/files/pcpepg.xtc",
        analysis_type = gorder.analysis_types.AAOrder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen"
        ),
        leaflets = gorder.leaflets.GlobalClassification("@membrane", "name P"),
        silent = True,
        overwrite = True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 51
    assert len(results.molecules()) == 3

    avg_order = results.average_order()
    assert compare_orders(avg_order.total().value(), 0.1423)
    assert compare_orders(avg_order.upper().value(), 0.1411)
    assert compare_orders(avg_order.lower().value(), 0.1434)
    assert avg_order.total().error() is None
    assert avg_order.upper().error() is None
    assert avg_order.lower().error() is None

    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    expected_average_orders = [0.1455, 0.1378, 0.1561]
    expected_average_upper = [0.1492, 0.1326, 0.1522]
    expected_average_lower = [0.1419, 0.1431, 0.1606]
    expected_atom_numbers = [37, 40, 38]
    expected_molecule_names = ["POPE", "POPC", "POPG"]

    expected_atom_indices = [32, 41, 34]
    expected_atom_names = ["C32", "C32", "C32"]
    expected_atom_order = [0.2226, 0.2363, 0.2247]
    expected_atom_upper = [0.2131, 0.2334, 0.2484]
    expected_atom_lower = [0.2319, 0.2391, 0.1976]

    expected_bond_numbers = [2, 2, 2]

    expected_atom2_indices = [34, 43, 36]
    expected_atom2_names = ["H2Y", "H2Y", "H2Y"]
    expected_atom2_order = [0.2040, 0.2317, 0.2020]
    expected_atom2_upper = [0.1876, 0.2507, 0.2254]
    expected_atom2_lower = [0.2203, 0.2126, 0.1752]

    for i, molecule in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        # AVERAGE
        avg_order = molecule.average_order()
        assert compare_orders(avg_order.total().value(), expected_average_orders[i])
        assert compare_orders(avg_order.upper().value(), expected_average_upper[i])
        assert compare_orders(avg_order.lower().value(), expected_average_lower[i])
        assert avg_order.total().error() is None
        assert avg_order.upper().error() is None
        assert avg_order.lower().error() is None

        avg_maps = molecule.average_ordermaps()
        assert avg_maps.total() is None
        assert avg_maps.upper() is None
        assert avg_maps.lower() is None

        # ATOMS
        assert len(molecule.atoms()) == expected_atom_numbers[i]
        atom = molecule.get_atom(expected_atom_indices[i])
        atom_type = atom.atom()
        assert atom_type.atom_name() == expected_atom_names[i]
        assert atom_type.relative_index() == expected_atom_indices[i]
        assert atom_type.residue_name() == expected_molecule_names[i]
        assert atom.molecule() == expected_molecule_names[i]

        order = atom.order()
        assert compare_orders(order.total().value(), expected_atom_order[i])
        assert compare_orders(order.upper().value(), expected_atom_upper[i])
        assert compare_orders(order.lower().value(), expected_atom_lower[i])
        assert order.total().error() is None
        assert order.upper().error() is None
        assert order.lower().error() is None

        maps = atom.ordermaps()
        assert maps.total() is None
        assert maps.upper() is None
        assert maps.lower() is None

        # BONDS
        assert len(atom.bonds()) == expected_bond_numbers[i]
        bond = atom.get_bond(expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.atom_name() == expected_atom_names[i]
        assert a1.relative_index() == expected_atom_indices[i]
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == expected_atom2_names[i]
        assert a2.relative_index() == expected_atom2_indices[i]
        assert a2.residue_name() == expected_molecule_names[i]
        assert bond.molecule() == expected_molecule_names[i]

        bond_order = bond.order()
        assert compare_orders(bond_order.total().value(), expected_atom2_order[i])
        assert compare_orders(bond_order.upper().value(), expected_atom2_upper[i])
        assert compare_orders(bond_order.lower().value(), expected_atom2_lower[i])
        assert bond_order.total().error() is None
        assert bond_order.upper().error() is None
        assert bond_order.lower().error() is None

        bond_maps = bond.ordermaps()
        assert bond_maps.total() is None
        assert bond_maps.upper() is None
        assert bond_maps.lower() is None

        # BOND FROM MOL
        bond = molecule.get_bond(expected_atom_indices[i], expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.relative_index() == expected_atom_indices[i]
        assert a2.relative_index() == expected_atom2_indices[i]

        # ERROR HANDLING
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_atom(145)
        assert "atom with the given relative index does not exist" in str(excinfo.value)
            
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(7, 19)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_cg_order_leaflets():
    analysis = gorder.Analysis(
        structure="../tests/files/cg.tpr",
        trajectory="../tests/files/cg.xtc",
        analysis_type=gorder.analysis_types.CGOrder("@membrane"),
        leaflets=gorder.leaflets.GlobalClassification("@membrane", "name PO4"),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 101
    assert len(results.molecules()) == 3

    avg_order = results.average_order()
    assert compare_orders(avg_order.total().value(), 0.2962)
    assert compare_orders(avg_order.upper().value(), 0.2971)
    assert compare_orders(avg_order.lower().value(), 0.2954)
    assert avg_order.total().error() is None
    assert avg_order.upper().error() is None
    assert avg_order.lower().error() is None

    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    expected_molecule_names = ["POPC", "POPE", "POPG"]
    expected_average_orders = [0.2943, 0.2972, 0.3059]
    expected_average_upper = [0.2965, 0.2965, 0.3085]
    expected_average_lower = [0.2920, 0.2980, 0.3033]
    expected_bond_orders = [0.3682, 0.3759, 0.3789]
    expected_bond_upper = [0.3647, 0.3713, 0.4129]
    expected_bond_lower = [0.3717, 0.3806, 0.3449]

    for i, molecule in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        # AVERAGE
        avg_order = molecule.average_order()
        assert compare_orders(avg_order.total().value(), expected_average_orders[i])
        assert compare_orders(avg_order.upper().value(), expected_average_upper[i])
        assert compare_orders(avg_order.lower().value(), expected_average_lower[i])
        assert avg_order.total().error() is None
        assert avg_order.upper().error() is None
        assert avg_order.lower().error() is None

        avg_maps = molecule.average_ordermaps()
        assert avg_maps.total() is None
        assert avg_maps.upper() is None
        assert avg_maps.lower() is None

        # BONDS
        assert len(molecule.bonds()) == 11
        bond = molecule.get_bond(4, 5)
        a1, a2 = bond.atoms()
        assert a1.atom_name() == "C1A"
        assert a1.relative_index() == 4
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == "D2A"
        assert a2.relative_index() == 5
        assert a2.residue_name() == expected_molecule_names[i]

        order = bond.order()
        assert compare_orders(order.total().value(), expected_bond_orders[i])
        assert compare_orders(order.upper().value(), expected_bond_upper[i])
        assert compare_orders(order.lower().value(), expected_bond_lower[i])
        assert order.total().error() is None
        assert order.upper().error() is None
        assert order.lower().error() is None

        reverse_bond = molecule.get_bond(5, 4)
        a1_rev, a2_rev = reverse_bond.atoms()
        assert a1_rev.relative_index() == 4
        assert a2_rev.relative_index() == 5

        # NONEXISTENT BONDS
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(1, 3)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(15, 16)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_aa_order_error_leaflets():
    analysis = gorder.Analysis(
        structure="../tests/files/pcpepg.tpr",
        trajectory="../tests/files/pcpepg.xtc",
        analysis_type=gorder.analysis_types.AAOrder(
            "@membrane and element name carbon",
            "@membrane and element name hydrogen"
        ),
        leaflets=gorder.leaflets.GlobalClassification("@membrane", "name P"),
        estimate_error=gorder.estimate_error.EstimateError(),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 51
    assert len(results.molecules()) == 3

    avg_order = results.average_order()
    assert compare_orders(avg_order.total().value(), 0.1423)
    assert compare_orders(avg_order.total().error(), 0.0026)
    assert compare_orders(avg_order.upper().value(), 0.1411)
    assert compare_orders(avg_order.upper().error(), 0.0024)
    assert compare_orders(avg_order.lower().value(), 0.1434)
    assert compare_orders(avg_order.lower().error(), 0.0031)

    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    expected_atom_numbers = [37, 40, 38]
    expected_molecule_names = ["POPE", "POPC", "POPG"]

    expected_atom_indices = [32, 41, 34]
    expected_atom_names = ["C32", "C32", "C32"]

    expected_bond_numbers = [2, 2, 2]

    expected_atom2_indices = [34, 43, 36]
    expected_atom2_names = ["H2Y", "H2Y", "H2Y"]

    for i, molecule in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        avg_order = molecule.average_order()
        assert avg_order.total().error() is not None
        assert avg_order.upper().error() is not None
        assert avg_order.lower().error() is not None

        # Test ordermaps absence
        avg_maps = molecule.average_ordermaps()
        assert avg_maps.total() is None
        assert avg_maps.upper() is None
        assert avg_maps.lower() is None

        convergence = molecule.convergence()
        assert len(convergence.frames()) == 51
        assert convergence.total() is not None
        assert convergence.upper() is not None
        assert convergence.lower() is not None

        assert len(molecule.atoms()) == expected_atom_numbers[i]
        atom = molecule.get_atom(expected_atom_indices[i])
        atom_type = atom.atom()
        assert atom_type.atom_name() == expected_atom_names[i]
        assert atom_type.relative_index() == expected_atom_indices[i]
        assert atom_type.residue_name() == expected_molecule_names[i]
        assert atom.molecule() == expected_molecule_names[i]

        order = atom.order()
        assert order.total().error() is not None
        assert order.upper().error() is not None
        assert order.lower().error() is not None

        maps = atom.ordermaps()
        assert maps.total() is None
        assert maps.upper() is None
        assert maps.lower() is None

        assert len(atom.bonds()) == expected_bond_numbers[i]
        bond = atom.get_bond(expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.atom_name() == expected_atom_names[i]
        assert a1.relative_index() == expected_atom_indices[i]
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == expected_atom2_names[i]
        assert a2.relative_index() == expected_atom2_indices[i]
        assert a2.residue_name() == expected_molecule_names[i]
        assert bond.molecule() == expected_molecule_names[i]

        bond_order = bond.order()
        assert bond_order.total().error() is not None
        assert bond_order.upper().error() is not None
        assert bond_order.lower().error() is not None

        bond_maps = bond.ordermaps()
        assert bond_maps.total() is None
        assert bond_maps.upper() is None
        assert bond_maps.lower() is None

        bond = molecule.get_bond(expected_atom_indices[i], expected_atom2_indices[i])
        a1, a2 = bond.atoms()
        assert a1.relative_index() == expected_atom_indices[i]
        assert a2.relative_index() == expected_atom2_indices[i]

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_atom(145)
        assert "atom with the given relative index does not exist" in str(excinfo.value)
            
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(7, 19)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_cg_order_error_leaflets():
    analysis = gorder.Analysis(
        structure="../tests/files/cg.tpr",
        trajectory="../tests/files/cg.xtc",
        analysis_type=gorder.analysis_types.CGOrder("@membrane"),
        leaflets=gorder.leaflets.GlobalClassification("@membrane", "name PO4"),
        estimate_error=gorder.estimate_error.EstimateError(),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 101
    assert len(results.molecules()) == 3

    avg_order = results.average_order()
    assert compare_orders(avg_order.total().value(), 0.2962)
    assert compare_orders(avg_order.total().error(), 0.0050)
    assert compare_orders(avg_order.upper().value(), 0.2971)
    assert compare_orders(avg_order.upper().error(), 0.0049)
    assert compare_orders(avg_order.lower().value(), 0.2954)
    assert compare_orders(avg_order.lower().error(), 0.0056)

    assert results.average_ordermaps().total() is None
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    expected_molecule_names = ["POPC", "POPE", "POPG"]

    for i, molecule in enumerate(results.molecules()):
        assert molecule.molecule() == expected_molecule_names[i]

        # AVERAGE
        avg_order = molecule.average_order()
        assert avg_order.total().error() is not None
        assert avg_order.upper().error() is not None
        assert avg_order.lower().error() is not None

        avg_maps = molecule.average_ordermaps()
        assert avg_maps.total() is None
        assert avg_maps.upper() is None
        assert avg_maps.lower() is None

        convergence = molecule.convergence()
        assert len(convergence.frames()) == 101
        assert convergence.total() is not None
        assert convergence.upper() is not None
        assert convergence.lower() is not None

        # BONDS
        assert len(molecule.bonds()) == 11
        bond = molecule.get_bond(4, 5)
        a1, a2 = bond.atoms()
        assert a1.atom_name() == "C1A"
        assert a1.relative_index() == 4
        assert a1.residue_name() == expected_molecule_names[i]
        assert a2.atom_name() == "D2A"
        assert a2.relative_index() == 5
        assert a2.residue_name() == expected_molecule_names[i]

        bond_order = bond.order()
        assert bond_order.total().error() is not None
        assert bond_order.upper().error() is not None
        assert bond_order.lower().error() is not None

        reverse_bond = molecule.get_bond(5, 4)
        a1_rev, a2_rev = reverse_bond.atoms()
        assert a1_rev.relative_index() == 4
        assert a2_rev.relative_index() == 5

        # NONEXISTENT BONDS
        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(1, 3)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

        with pytest.raises(gorder.exceptions.APIError) as excinfo:
            molecule.get_bond(15, 16)
        assert "bond specified by the given relative indices does not exist" in str(excinfo.value)

def test_aa_order_ordermaps():
    analysis = gorder.Analysis(
        structure = "../tests/files/pcpepg.tpr",
        trajectory = "../tests/files/pcpepg.xtc",
        analysis_type = gorder.analysis_types.AAOrder(
            "resname POPC and name C22 C24 C218", 
            "@membrane and element name hydrogen"
        ),
        ordermap = gorder.ordermap.OrderMap(bin_size = [0.1, 4.0], min_samples = 5),
        silent = True,
        overwrite = True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 51
    assert len(results.molecules()) == 1

    # SYSTEM
    assert results.average_order().total() is not None
    assert results.average_order().upper() is None
    assert results.average_order().lower() is None
    
    map = results.average_ordermaps().total()
    assert compare_orders(map.get_at(0.6, 8.0), 0.1653)
    assert compare_orders(map.get_at(4.3, 0.0), 0.1340)
    assert compare_orders(map.get_at(9.2, 4.0), 0.1990)
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    # MOLECULE
    molecule = results.get_molecule("POPC")
    map = molecule.average_ordermaps().total()
    assert molecule.average_ordermaps().upper() is None
    assert molecule.average_ordermaps().lower() is None

    span_x = map.span_x()
    span_y = map.span_y()
    bin = map.tile_dim()

    assert math.isclose(span_x[0], 0.0, rel_tol = 1e-5)
    assert math.isclose(span_x[1], 9.15673, rel_tol = 1e-5)
    assert math.isclose(span_y[0], 0.0, rel_tol = 1e-5)
    assert math.isclose(span_y[1], 9.15673, rel_tol = 1e-5)
    assert math.isclose(bin[0], 0.1, rel_tol = 1e-5)
    assert math.isclose(bin[1], 4.0, rel_tol = 1e-5)

    assert compare_orders(map.get_at(0.6, 8.0), 0.1653)
    assert compare_orders(map.get_at(4.3, 0.0), 0.1340)
    assert compare_orders(map.get_at(9.2, 4.0), 0.1990)

    # ATOM
    atom = molecule.get_atom(47)
    map = atom.ordermaps().total()
    assert atom.ordermaps().upper() is None
    assert atom.ordermaps().lower() is None

    assert compare_orders(map.get_at(0.6, 8.0), 0.2224)
    assert compare_orders(map.get_at(4.3, 0.0), 0.1532)
    assert compare_orders(map.get_at(9.2, 4.0), 0.0982)

    # BOND
    bond = atom.get_bond(49)
    map = bond.ordermaps().total()
    assert bond.ordermaps().upper() is None
    assert bond.ordermaps().lower() is None

    assert compare_orders(map.get_at(0.6, 8.0), 0.2901)
    assert compare_orders(map.get_at(4.3, 0.0), 0.1163)
    assert math.isnan(map.get_at(9.2, 4.0))

    # EXTRACT CHECK
    (extracted_x, extracted_y, extracted_values) = map.extract()
    assert len(extracted_x) == 93
    assert len(extracted_y) == 3
   
    for (real, expected) in zip(extracted_x, [x / 10 for x in range(0, 93)]):
        assert math.isclose(real, expected, rel_tol = 1e-5)
   
    for (real, expected) in zip(extracted_y, [0.0, 4.0, 8.0]):
        assert math.isclose(real, expected, rel_tol = 1e-5)
   
    for (xi, x) in enumerate(extracted_x):
        for (yi, y) in enumerate(extracted_y):
            get = map.get_at(x, y)
            ext = extracted_values[xi][yi]
            if math.isnan(get) and math.isnan(ext):
                continue
            assert compare_orders(map.get_at(x, y), extracted_values[xi][yi])

def test_cg_order_ordermaps():
    analysis = gorder.Analysis(
        structure="../tests/files/cg.tpr",
        trajectory="../tests/files/cg.xtc",
        analysis_type=gorder.analysis_types.CGOrder("resname POPC and name C1B C2B C3B C4B"),
        ordermap=gorder.ordermap.OrderMap(
            bin_size=[1.0, 1.0],
            min_samples=10
        ),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 101
    assert len(results.molecules()) == 1

    # SYSTEM
    assert results.average_order().total() is not None
    assert results.average_order().upper() is None
    assert results.average_order().lower() is None
    
    map = results.average_ordermaps().total()
    assert compare_orders(map.get_at(1.0, 8.0), 0.3590)
    assert compare_orders(map.get_at(7.0, 0.0), 0.3765)
    assert compare_orders(map.get_at(13.0, 11.0), 0.4296)
    assert results.average_ordermaps().upper() is None
    assert results.average_ordermaps().lower() is None

    # MOLECULE
    molecule = results.get_molecule("POPC")
    mol_maps = molecule.average_ordermaps()
    
    total_map = mol_maps.total()
    assert mol_maps.upper() is None
    assert mol_maps.lower() is None

    span_x = total_map.span_x()
    span_y = total_map.span_y()
    bin_size = total_map.tile_dim()
    
    assert math.isclose(span_x[0], 0.0, rel_tol=1e-5)
    assert math.isclose(span_x[1], 12.747616, rel_tol=1e-5)
    assert math.isclose(span_y[0], 0.0, rel_tol=1e-5)
    assert math.isclose(span_y[1], 12.747616, rel_tol=1e-5)
    assert math.isclose(bin_size[0], 1.0, rel_tol=1e-5)
    assert math.isclose(bin_size[1], 1.0, rel_tol=1e-5)

    assert compare_orders(total_map.get_at(1.0, 8.0), 0.3590)
    assert compare_orders(total_map.get_at(7.0, 0.0), 0.3765)
    assert compare_orders(total_map.get_at(13.0, 11.0), 0.4296)

    # BOND
    bond = molecule.get_bond(9, 10)
    bond_maps = bond.ordermaps()
    
    bond_total = bond_maps.total()
    assert bond_maps.upper() is None
    assert bond_maps.lower() is None

    assert compare_orders(bond_total.get_at(1.0, 8.0), 0.3967)
    assert compare_orders(bond_total.get_at(7.0, 0.0), 0.3213)
    assert compare_orders(bond_total.get_at(13.0, 11.0), 0.4104)

    # EXTRACT CHECK
    (extracted_x, extracted_y, extracted_values) = map.extract()
    assert len(extracted_x) == 14
    assert len(extracted_y) == 14
   
    for (real, expected) in zip(extracted_x, [x for x in range(0, 14)]):
        assert math.isclose(real, expected, rel_tol = 1e-5)
   
    for (real, expected) in zip(extracted_y, [y for y in range(0, 14)]):
        assert math.isclose(real, expected, rel_tol = 1e-5)
   
    for (xi, x) in enumerate(extracted_x):
        for (yi, y) in enumerate(extracted_y):
            get = map.get_at(x, y)
            ext = extracted_values[xi][yi]
            if math.isnan(get) and math.isnan(ext):
                continue
            assert compare_orders(map.get_at(x, y), extracted_values[xi][yi])

def test_aa_order_ordermaps_leaflets():
    analysis = gorder.Analysis(
        structure="../tests/files/pcpepg.tpr",
        trajectory="../tests/files/pcpepg.xtc",
        analysis_type=gorder.analysis_types.AAOrder(
            "resname POPC and name C22 C24 C218",
            "@membrane and element name hydrogen"
        ),
        leaflets=gorder.leaflets.GlobalClassification("@membrane", "name P"),
        ordermap=gorder.ordermap.OrderMap(
            bin_size=[0.1, 4.0],
            min_samples=5
        ),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 51
    assert len(results.molecules()) == 1

    # SYSTEM
    sys_maps = results.average_ordermaps()
    assert compare_orders(sys_maps.total().get_at(0.6, 8.0), 0.1653)
    assert compare_orders(sys_maps.total().get_at(9.2, 4.0), 0.1990)
    assert compare_orders(sys_maps.upper().get_at(0.6, 8.0), 0.1347)
    assert compare_orders(sys_maps.upper().get_at(9.2, 4.0), 0.3196)
    assert compare_orders(sys_maps.lower().get_at(0.6, 8.0), 0.2104)
    assert compare_orders(sys_maps.lower().get_at(9.2, 4.0), 0.1106)

    # MOLECULE
    molecule = results.get_molecule("POPC")
    mol_maps = molecule.average_ordermaps()
    
    total_map = mol_maps.total()
    span_x = total_map.span_x()
    span_y = total_map.span_y()
    bin_size = total_map.tile_dim()
    
    assert math.isclose(span_x[0], 0.0, rel_tol=1e-5)
    assert math.isclose(span_x[1], 9.15673, rel_tol=1e-5)
    assert math.isclose(span_y[0], 0.0, rel_tol=1e-5)
    assert math.isclose(span_y[1], 9.15673, rel_tol=1e-5)
    assert math.isclose(bin_size[0], 0.1, rel_tol=1e-5)
    assert math.isclose(bin_size[1], 4.0, rel_tol=1e-5)
    
    assert compare_orders(total_map.get_at(0.6, 8.0), 0.1653)
    assert compare_orders(total_map.get_at(9.2, 4.0), 0.1990)
    upper_map = mol_maps.upper()
    assert compare_orders(upper_map.get_at(0.6, 8.0), 0.1347)
    assert compare_orders(upper_map.get_at(9.2, 4.0), 0.3196)
    lower_map = mol_maps.lower()
    assert compare_orders(lower_map.get_at(0.6, 8.0), 0.2104)
    assert compare_orders(lower_map.get_at(9.2, 4.0), 0.1106)

    # ATOM
    atom = molecule.get_atom(47)
    atom_maps = atom.ordermaps()

    atom_total = atom_maps.total()
    assert compare_orders(atom_total.get_at(0.6, 8.0), 0.2224)
    assert compare_orders(atom_total.get_at(9.2, 4.0), 0.0982)
    atom_upper = atom_maps.upper()
    assert compare_orders(atom_upper.get_at(0.6, 8.0), 0.2039)
    assert math.isnan(atom_upper.get_at(9.2, 4.0))
    atom_lower = atom_maps.lower()
    assert compare_orders(atom_lower.get_at(0.6, 8.0), 0.2540)
    assert math.isnan(atom_lower.get_at(9.2, 4.0))

    # BOND
    bond = atom.get_bond(49)
    bond_maps = bond.ordermaps()
    
    bond_total = bond_maps.total()
    assert compare_orders(bond_total.get_at(0.6, 8.0), 0.2901)
    assert math.isnan(bond_total.get_at(9.2, 4.0))
    bond_upper = bond_maps.upper()
    assert compare_orders(bond_upper.get_at(0.6, 8.0), 0.3584)
    assert math.isnan(bond_upper.get_at(9.2, 4.0))
    bond_lower = bond_maps.lower()
    assert compare_orders(bond_lower.get_at(0.6, 8.0), 0.1715)
    assert math.isnan(bond_lower.get_at(9.2, 4.0))

def test_cg_order_ordermaps_leaflets():
    analysis = gorder.Analysis(
        structure="../tests/files/cg.tpr",
        trajectory="../tests/files/cg.xtc",
        analysis_type=gorder.analysis_types.CGOrder("resname POPC and name C1B C2B C3B C4B"),
        leaflets=gorder.leaflets.GlobalClassification("@membrane", "name PO4"),
        ordermap=gorder.ordermap.OrderMap(
            bin_size=[1.0, 1.0],
            min_samples=10
        ),
        silent=True,
        overwrite=True,
    )

    results = analysis.run()

    assert results.n_analyzed_frames() == 101
    assert len(results.molecules()) == 1

    # SYSTEM
    sys_maps = results.average_ordermaps()
    assert compare_orders(sys_maps.total().get_at(1.0, 8.0), 0.3590)
    assert compare_orders(sys_maps.total().get_at(13.0, 11.0), 0.4296)
    assert compare_orders(sys_maps.upper().get_at(1.0, 8.0), 0.3418)
    assert compare_orders(sys_maps.upper().get_at(13.0, 11.0), 0.4051)
    assert compare_orders(sys_maps.lower().get_at(1.0, 8.0), 0.3662)
    assert compare_orders(sys_maps.lower().get_at(13.0, 11.0), 0.4506)

    # MOLECULE
    molecule = results.get_molecule("POPC")
    mol_maps = molecule.average_ordermaps()

    total_map = mol_maps.total()
    span_x = total_map.span_x()
    span_y = total_map.span_y()
    bin_size = total_map.tile_dim()
    
    assert math.isclose(span_x[0], 0.0, rel_tol=1e-5)
    assert math.isclose(span_x[1], 12.747616, rel_tol=1e-5)
    assert math.isclose(span_y[0], 0.0, rel_tol=1e-5)
    assert math.isclose(span_y[1], 12.747616, rel_tol=1e-5)
    assert math.isclose(bin_size[0], 1.0, rel_tol=1e-5)
    assert math.isclose(bin_size[1], 1.0, rel_tol=1e-5)

    assert compare_orders(total_map.get_at(1.0, 8.0), 0.3590)
    assert compare_orders(total_map.get_at(13.0, 11.0), 0.4296)
    upper_map = mol_maps.upper()
    assert compare_orders(upper_map.get_at(1.0, 8.0), 0.3418)
    assert compare_orders(upper_map.get_at(13.0, 11.0), 0.4051)
    lower_map = mol_maps.lower()
    assert compare_orders(lower_map.get_at(1.0, 8.0), 0.3662)
    assert compare_orders(lower_map.get_at(13.0, 11.0), 0.4506)

    # BOND
    bond = molecule.get_bond(9, 10)
    bond_maps = bond.ordermaps()
    
    bond_total = bond_maps.total()
    assert compare_orders(bond_total.get_at(1.0, 8.0), 0.3967)
    assert compare_orders(bond_total.get_at(13.0, 11.0), 0.4104)
    bond_upper = bond_maps.upper()
    assert compare_orders(bond_upper.get_at(1.0, 8.0), 0.3573)
    assert compare_orders(bond_upper.get_at(13.0, 11.0), 0.4807)
    bond_lower = bond_maps.lower()
    assert compare_orders(bond_lower.get_at(1.0, 8.0), 0.4118)
    assert compare_orders(bond_lower.get_at(13.0, 11.0), 0.3563)