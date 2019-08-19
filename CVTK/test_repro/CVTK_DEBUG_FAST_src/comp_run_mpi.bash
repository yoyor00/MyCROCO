#!/bin/bash
#===================================
#PBS -q omp
#PBS -l walltime=02:00:00
#PBS -l ncpus=4
#PBS -l mpiprocs=4
#PBS -l mem=20gb
#PBS -j oe
#set -x
set -e

source CONFIGURE_GLOBAL
source configure_file

test -z "$CI_CROCO_PWD" && cd $SUBMIT_DIR || cd "$CI_CROCO_PWD"
echo "$CI_CROCO_PWD"
echo $PBS_O_LOGNAME
#===================================

par1='MPI'
# Compile
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log 2>&1
/bin/mv croco croco_${par1}.exe

# Run
$MPIRUN -np $NBPROCS ./croco_${par1}.exe $CROCOIN > mpi_${NBPROCS}_${TEST_NAME}.log || $CROCO_CI_MPIRUN 

