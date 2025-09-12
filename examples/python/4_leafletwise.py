# Leaflet-wise calculation of atomistic lipid 
# order parameters. All output formats.
# Analysis performed using 4 threads.
# Manual: ladme.github.io/gorder-manual

# Run using `uv`: `uv run 4_leafletwise.py`.

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
    # Assign lipids into leaflets based on the global center of geometry of the membrane.
    leaflets = gorder.leaflets.GlobalClassification(
        # All particles of the membrane for center of geometry calculation.
        membrane = "@membrane",
        # Phosphorus atoms of membrane lipids identifying lipid heads.
        # Only one such atom per molecule!
        heads = "name P",
        # Perform the leaflet assignment only for the first trajectory frame
        # and then use the same assignment for all the following frames.
        # DELETE THIS LINE if lipids are able to flip-flop in your membrane!
        frequency = gorder.Frequency.once(),
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
        print(f"   Atom type '{atom.atom().atom_name()}'") 
        order = atom.order()
        print(f"     > total: {order.total().value():.4f}") # pyright: ignore[reportOptionalMemberAccess]
        print(f"     > upper: {order.upper().value():.4f}") # pyright: ignore[reportOptionalMemberAccess]
        print(f"     > lower: {order.lower().value():.4f}") # pyright: ignore[reportOptionalMemberAccess]