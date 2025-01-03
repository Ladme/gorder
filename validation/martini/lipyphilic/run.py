
import numpy as np
import MDAnalysis as mda
from lipyphilic.analysis.order_parameter import SCC

def main():
    u = mda.Universe("../files/system.tpr", "../files/traj.xtc")

    analyses = [SCC(universe = u, tail_sel = "name NC3 PO4"),
    SCC(universe = u, tail_sel = "name PO4 GL1"),
    SCC(universe = u, tail_sel = "name GL1 GL2"),

    SCC(universe = u, tail_sel = "name GL1 C1A"),
    SCC(universe = u, tail_sel = "name C1A D2A"),
    SCC(universe = u, tail_sel = "name D2A C3A"),
    SCC(universe = u, tail_sel = "name C3A C4A"),

    SCC(universe = u, tail_sel = "name GL2 C1B"),
    SCC(universe = u, tail_sel = "name C1B C2B"),
    SCC(universe = u, tail_sel = "name C2B C3B"),
    SCC(universe = u, tail_sel = "name C3B C4B"),
    ]

    with open("results.dat", "w") as output:
        for (i, analysis) in enumerate(analyses):
            analysis.run(start = None, stop = None, step = None, verbose = True)

            output.write(f"{i} {np.mean(analysis.SCC):.4f}\n")
    

if __name__ == "__main__":
    main()