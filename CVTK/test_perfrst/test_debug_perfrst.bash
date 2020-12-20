#!/bin/bash

# VERIFIER QUES LES CPPDEFS.H SONT LES MEMES!!

#0
./compile_croco_PR.bash

#1
rm -Rf check_file_1_0
mkdir OUT OUT_seg12
rm -Rf OUT/* OUT_seg12/*

#2
#./croco_EXACT_RST_DBWRITE croco.in.rst_seg12_PR
mpirun -np 4 ./croco_EXACT_RST_DBWRITE croco.in.rst_seg12_PR


#3
#./croco_EXACT_RST_DBREAD croco.in
mpirun -np 4 ./croco_EXACT_RST_DBREAD croco.in
