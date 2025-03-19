#!/bin/bash

hyperfine --runs 1 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'python do-order.py ../files/traj_whole.xtc ../files/system.tpr 0 1000000 1 0 0 1 512 POPC'