#!/bin/bash

hyperfine --runs 5 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'order -c ../files/system.gro -f ../files/traj.xtc -a "resname POPC" -i ../files/martini_phospholipids.itp'