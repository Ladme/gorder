#!/bin/bash

hyperfine --runs 1 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'buildH -c ../files/system.gro -t ../files/traj.xtc -l Berger_POPC -d Berger_POPC_min.def' --show-output