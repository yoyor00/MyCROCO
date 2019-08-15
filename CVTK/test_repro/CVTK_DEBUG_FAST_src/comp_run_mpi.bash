#!/bin/bash
#===================================
#PBS -q omp
#PBS -l walltime=02:00:00
#PBS -l ncpus=4
#PBS -l mpiprocs=4
#PBS -l mem=20gb
#PBS -j oe
set -x
set -e
test -z "$CI_CROCO_PWD" && cd $PBS_O_WORKDIR || cd "$CI_CROCO_PWD"
echo "$CI_CROCO_PWD"
echo $PBS_O_LOGNAME
#===================================
source configure_file

par1='MPI'
# Compile
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log
/bin/mv croco croco_${par1}.exe

# Run
MPIRUN=${CROCO_CI_MPIRUN:-${MPI_LAUNCH}}

$MPIRUN -np $NBPROCS ./croco_${par1}.exe $CROCOIN > mpi_${NBPROCS}_${TEST_NAME}.log || $CROCO_CI_MPIRUN 

