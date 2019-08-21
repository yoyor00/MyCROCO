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

if [ ! -f ${TEST_NAME}_steps ]; then 
  echo 'N' > ${TEST_NAME}_steps
  echo 'N' >> ${TEST_NAME}_steps
  echo 'N' >> ${TEST_NAME}_steps
fi

#===================================

par1='MPI'
# Compile
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="$(tput setaf 1 ; tput bold)${msg1}$( tput sgr0)"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1 || { echo -e "   $msg2" > /dev/tty ; echo -e $msg1 ; exit 1 ; }
/bin/mv croco croco_${par1}.exe

# Run
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="$(tput setaf 1 ; tput bold)${msg1}$( tput sgr0)"
$MPIRUN -np $NBPROCS ./croco_${par1}.exe $CROCOIN > mpi_${NBPROCS}_${TEST_NAME}.log 2>&1  || { echo -e "   $msg2" > /dev/tty ; echo -e $msg1 ; exit 2 ; }

# Additional check in case of clean stop before the end
SUCCESS1=$(tail -n 2  mpi_${NBPROCS}_${TEST_NAME}.log | head -n 1)
SUCCESS=$(echo $SUCCESS1 | sed -e " s/\ //g")
if [ "$SUCCESS" != 'MAIN:DONE' ]; then
  echo -e "   $msg2" > /dev/tty  
  echo -e $msg1 
  exit 2 
fi	