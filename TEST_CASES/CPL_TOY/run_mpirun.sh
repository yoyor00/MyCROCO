#!/bin/bash
# Run the mpirun command with redirection
mpirun --allow-run-as-root -np 2 toywav : -np 4 crocox > out.log 2>&1
