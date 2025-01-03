#!/bin/bash

module add gromacs:2021.4-plumed

hyperfine --warmup 3 'gmx order -s ../files/system.tpr -f ../files/traj.xtc' --runs 5