# Calculation of atomistic lipid order
# parameters with error estimation.
# Only full membrane, all output formats.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 7_estimate_error.py`.

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
    # Calculate atomistic order parameters.
    analysis_type = gorder.analysis_types.AAOrder(
        # Calculate order parameters for carbons of the palmitoyl and oleoyl chains.
        heavy_atoms = "@membrane and name r'C3.+|C2.+'",
        # These are the hydrogens of the lipids in the membrane. 
        # `gorder` will search for bonds between the specified heavy atoms and hydrogens.
        hydrogens = "@membrane and element name hydrogen"
    ),
    # Estimate the error of the analysis.
    estimate_error = gorder.estimate_error.EstimateError(
        # Output the convergence of the calculated order parameters into an xvg file.
        # This file can be plotted using `xmgrace`.
        # If the order parameters stabilize, your simulation is converged.
        output_convergence = "convergence.xvg"
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
        order = atom.order().total()
        print(f"   Atom type '{atom.atom().atom_name()}': {order.value():.4f} ± {order.error():.4f}") # pyright: ignore[reportOptionalMemberAccess]

    # You can also directly access the convergence data for each molecule type.
    convergence = molecule.convergence()
    assert convergence is not None
    print("   Average order parameter in time")
    for (frame, order) in zip(convergence.frames(), convergence.total()): # pyright: ignore[reportArgumentType]
        print(f"     > Frame {frame:5}: {order:.4f}")

    
    