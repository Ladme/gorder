#!/bin/bash

# conda activate loos

hyperfine --runs 1 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'order_params ../files/system.psf ../files/traj.xtc "resname == \"POPC\" && name =~ \"^C2[1-9][0-9]?\"" 2 18 --1 && order_params ../files/system.psf ../files/traj.xtc "resname == \"POPC\" && name =~ \"^C3[1-9][0-9]?\"" 2 16 --1' --show-output
