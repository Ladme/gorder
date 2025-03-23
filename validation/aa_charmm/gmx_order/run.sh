#!/bin/bash

# To calculate order parameters for the individual atoms, we have to read through the trajectory twice: once for the palmitoyl chain and once for the oleoyl chain.
# Note that `gmx order` is not able to calculate order parameters for unsatured carbons correctly!

hyperfine --runs 5 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'gmx_mpi order -s ../files/system.tpr -f ../files/traj.xtc -n index_palmitoyl.ndx -od order_palmitoyl.xvg && gmx_mpi order -s ../files/system.tpr -f ../files/traj.xtc -n index_oleoyl.ndx -od order_oleoyl.xvg'