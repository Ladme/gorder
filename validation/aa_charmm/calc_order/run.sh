#!/bin/bash

# Run uv first to prepare the virtual environment.

hyperfine --runs 1 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'uv run calc_order.py' --show-output