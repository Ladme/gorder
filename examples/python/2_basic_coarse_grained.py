# Basic calculation of coarse-grained lipid order parameters.
# Only full membrane, all output formats.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 2_basic_coarse_grained.py`.

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
    # Calculate coarse-grained order parameters.
    analysis_type = gorder.analysis_types.CGOrder(
        # Analyze all bonds between particles of the membrane.
        "@membrane",
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
    for bond in molecule.bonds():
        print(f"   Bond type '{bond.atoms()[0].atom_name()}-{bond.atoms()[1].atom_name()}': {bond.order().total().value():.4f}")