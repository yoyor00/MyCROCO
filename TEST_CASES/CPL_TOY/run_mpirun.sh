#!/bin/bash
# Run the mpirun command with redirection
mpirun --allow-run-as-root -np 2 toywav : -np 2 crocox > out.log 2>&1

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "Error: mpirun command failed" >&2
  exit 1
fi

echo "mpirun command completed successfully."