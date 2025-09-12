# Basic calculation of united-atom lipid order parameters.
# Only full membrane, all output formats.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 3_basic_united_atom.py`.

# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gorder @ git+https://github.com/Ladme/gorder.git#subdirectory=pygorder"
# ]
# ///

import gorder

analysis = gorder.Analysis(
    # Input Gromacs binary run file.
    structure = "system.tpr",
    # Input Gromacs compressed trajectory file.
    trajectory = "traj.xtc",
    # Calculate united-atom order parameters.
    analysis_type = gorder.analysis_types.UAOrder(
        # Calculate order parameters for saturated carbons of the palmitoyl and oleoyl chains.
        saturated = "resname POPC and name r'C(1[7-9]|2[0-36-9]|3[0-16-9]|4[0-9]|50|A[1-2])'",
        # Calculate order parameters for unsaturated carbons of the oleoyl chain.
        unsaturated = "resname POPC and name C24 C25"
    ),
    # Path to the output yaml file. Contains full results of the analysis.
    output_yaml = "order.yaml",
    # Path to the output "table" file. Simple to read.
    output_tab = "order.tab",
    # Path to the output xvg file. Can be plotted using `xmgrace`.
    output_xvg = "order.xvg",
    # Path to the output csv file. Open with a spreadsheet software or read with Pandas/Polars.
    output_csv = "order.csv",
    # Perform the analysis using 4 threads.
    n_threads = 4,
)

# Run the analysis.
results = analysis.run()

# Write the results of the analysis into output files.
results.write()

# Or access the results manually.
for molecule in results.molecules():
    print(f"Molecule type '{molecule.molecule()}':")
    for atom in molecule.atoms():
        print(f"   Atom type '{atom.atom().atom_name()}': {atom.order().total().value():.4f}") # pyright: ignore[reportOptionalMemberAccess]