#!/bin/bash

hyperfine --warmup 3 './gorder analyze.yaml' --runs 5