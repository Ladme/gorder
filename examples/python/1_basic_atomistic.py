# Basic calculation of atomistic lipid order parameters.
# Only full membrane, all output formats.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 1_basic_atomistic.py`.

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
    trajectory = "traj.xtc".
    # Calculate atomistic order parameters.
    analysis_type = gorder.analysis_types.AAOrder(
        # Calculate order parameters for carbons of the palmitoyl and oleoyl chains.
        "@membrane and name r'C3.+|C2.+'",
        # These are the hydrogens of the lipids in the membrane. 
        # `gorder` will search for bonds between the specified heavy atoms and hydrogens.
        "@membrane and element name hydrogen"
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

# TODO
