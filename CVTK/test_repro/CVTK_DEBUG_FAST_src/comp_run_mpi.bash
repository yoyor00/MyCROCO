#!/bin/bash
#===================================
#PBS -q omp
#PBS -l walltime=02:00:00
#PBS -l ncpus=4
#PBS -l mpiprocs=4
#PBS -l mem=20gb
#PBS -j oe
#set -x
set -eu

source CONFIGURE_GLOBAL
source configure_file

cd $SUBMIT_DIR
echo "$SUBMIT_DIR"

#===================================

par1='MPI'
# Compile
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1 || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 1 ; }
/bin/mv croco croco_${par1}.exe

# Run
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
$MPIRUN -np $NBPROCS ./croco_${par1}.exe $CROCOIN > mpi_${NBPROCS}_${TEST_NAME}.log 2>&1  || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 2 ; }

# Additional check in case of clean stop before the end
SUCCESS1=$(tail -n 2  mpi_${NBPROCS}_${TEST_NAME}.log | head -n 1)
SUCCESS=$(echo $SUCCESS1 | sed -e " s/\ //g")
if [ "$SUCCESS" != 'MAIN:DONE' ]; then
  echo -e "   $msg2" | tee -a mylog.txt
  echo -e $msg1 
  exit 2 
fi	