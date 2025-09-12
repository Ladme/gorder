# Calculation of atomistic lipid order 
# parameters and construction of ordermaps.
# Only full membrane, all output formats.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 6_ordermaps.py`.

# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gorder @ git+https://github.com/Ladme/gorder.git#subdirectory=pygorder",
#     "numpy>=2.3.3",
# ]
# ///

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
    # Construct maps of order parameters.
    ordermap = gorder.ordermap.OrderMap(
        # Save the generated ordermaps into this directory.
        output_directory = "ordermaps",
        # At least 50 samples must be collected for
        # a bin of the map to calculate order parameter for it.
        # (Order parameters calculated in rarely visited areas are unreliable.)
        min_samples = 50,
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

# You can also access the the ordermaps using Python but you need to import NumPy.
import numpy as np
for molecule in results.molecules():
    for atom in molecule.atoms():
        print(f"Molecule {molecule.molecule()} / Atom {atom.atom().atom_name()}")
        ordermap = atom.ordermaps().total()
        assert ordermap is not None
        (xdim, ydim, values) = ordermap.extract()
        # `xdim` is an 1D array containing positions of the grid tiles along the x-axis.
        # `ydim` is an 1D array containing positions of the grid tiles along the y-axis.
        # `values` is a 2D array containing the calculated order parameters.
        # `values` has dimensions of (len(xdim), len(ydim))

        np.set_printoptions(precision = 4, suppress = True, linewidth = 200)
        print("x ->", xdim)
        print("y ->", ydim)
        print(values)
