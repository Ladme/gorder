#!/bin/bash

hyperfine --runs 1 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' 'vmd -e calc_op.tcl'