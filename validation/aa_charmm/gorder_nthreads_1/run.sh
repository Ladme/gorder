#!/bin/bash

hyperfine --runs 5 --prepare 'sync; echo 3 | sudo tee /proc/sys/vm/drop_caches' './gorder analyze.yaml'