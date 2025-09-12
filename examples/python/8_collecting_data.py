# Calculation of atomistic lipid order 
# parameters while collecting leaflet
# assignment data and membrane normals.
# No output formats, access via Python only.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 8_collecting_data.py`.

# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gorder @ git+https://github.com/Ladme/gorder.git#subdirectory=pygorder",
#     "numpy>=2.3.3",
# ]
# ///

# Note that we need to have NumPy installed, but we do not need to import it.

import gorder

analysis = gorder.Analysis(
    # Input Gromacs binary run file.
    structure = "system.tpr",
    # Input Gromacs compressed trajectory file.
    trajectory = "traj.xtc",
    # Calculate atomistic order parameters.
    analysis_type = gorder.analysis_types.AAOrder(
        # Calculate order parameters for carbons of the palmitoyl and oleoyl chains.
        heavy_atoms = "@membrane and name r'C3.+|C2.+'",
        # These are the hydrogens of the lipids in the membrane. 
        # `gorder` will search for bonds between the specified heavy atoms and hydrogens.
        hydrogens = "@membrane and element name hydrogen"
    ),
    # Assign lipids into leaflets using spectral clustering.
    leaflets = gorder.leaflets.ClusteringClassification(
        # Phosphorus atoms of membrane lipids identifying lipid heads.
        # Only one such atom per molecule!
        heads = "name P",
        # Collect the leaflet assignment making it accessible (see below).
        collect = True,
        # Assignment will be performed for every analyzed trajectory frame (default).
    ),
    # Calculate membrane normal for each lipid molecule dynamically, based on the shape of the vesicle.
    membrane_normal = gorder.membrane_normal.DynamicNormal(
        # Phosphorus atoms of membrane lipids identifying lipid heads. 
        # Only one such atom per molecule!
        "name P",
        # Collect the calculated membrane normals making them accessible (see below).
        collect = True,
    ),
    # Perform the analysis using 4 threads.
    n_threads = 4,

    # We do not want any output files to be generated.
    # Everything will be handled using the Python API.
)

# Run the analysis.
results = analysis.run()

##################################
#   ACCESSING ORDER PARAMETERS   #
##################################
for molecule in results.molecules():
    print(f"Molecule type '{molecule.molecule()}':")
    for atom in molecule.atoms():
        print(f"   Atom type '{atom.atom().atom_name()}': {atom.order().total().value():.4f}") # pyright: ignore[reportOptionalMemberAccess]

##################################
#  ACCESSING LEAFLET ASSIGNMENT  #
##################################
leaflets = results.leaflets_data()
assert leaflets is not None

# Iterate through all molecule types.
for molecule in results.molecules():
    print(f"Molecule type '{molecule.molecule()}':")
    data = leaflets.get_molecule(molecule.molecule())
    assert data is not None

    # Iterate through all frames.
    for (frame_number, frame_data) in zip(leaflets.frames(), data):
        print(f"    Frame number {frame_number}:")

        # Iterate through assignment of all lipid molecules.
        for lipid in frame_data:
            print(f"      > {'upper' if lipid == 1 else 'lower'}")

##################################
#   ACCESSING MEMBRANE NORMALS   #
##################################
normals = results.normals_data()
assert normals is not None

# Iterate through all molecule types.
for molecule in results.molecules():
    print(f"Molecule type '{molecule.molecule()}':")
    data = normals.get_molecule(molecule.molecule())
    assert data is not None

    # Iterate through all frames.
    for (frame_number, frame_data) in zip(normals.frames(), data):
        print(f"    Frame number {frame_number}:")

        # Iterate through normal vectors of all lipid molecules.
        for normal in frame_data:
            print(f"      > [x: {normal[0]: 6.4f}, y: {normal[1]: 6.4f}, z: {normal[2]: 6.4f}]")