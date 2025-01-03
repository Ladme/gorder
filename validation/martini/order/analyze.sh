#!/bin/bash

hyperfine --warmup 3 './order -c ../files/system.gro -f ../files/traj.xtc -i ../files/martini_phospholipids.itp -a "resname POPC"'