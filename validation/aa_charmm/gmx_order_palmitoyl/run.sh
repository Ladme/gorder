#!/bin/bash

hyperfine --runs 5 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'gmx_mpi order -s ../files/system.tpr -f ../files/traj.xtc -n index_palmitoyl.ndx -od order_palmitoyl.xvg'