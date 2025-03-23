#!/bin/bash

# Warm-up used for installing dependencies.
hyperfine --warmup 1 --runs 1 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches ; rm -f ../files/.traj_whole.xtc_offsets.lock ../files/.traj_whole.xtc_offsets.npz' 'uv run analyze.py'