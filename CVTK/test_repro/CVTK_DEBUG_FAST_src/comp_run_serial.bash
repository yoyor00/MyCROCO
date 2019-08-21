#!/bin/bash
#===================================
#PBS -q sequentiel
#PBS -l walltime=02:00:00
#PBS -l mem=20gb
#PBS -j oe 
#set -x
set -e
set -u

source CONFIGURE_GLOBAL
source configure_file

cd "$SUBMIT_DIR"
echo "$SUBMIT_DIR"

#===================================
par1='SERIAL'

# Compilation
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="$(tput setaf 1 ; tput bold)${msg1}$( tput sgr0)"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1  || { echo -e "   $msg2" > /dev/stdin ; echo -e $msg1 ; exit 1 ; }
mv croco croco_${par1}.exe

# Run (and produce check file)
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="$(tput setaf 1 ; tput bold)${msg1}$( tput sgr0)"
./croco_${par1}.exe $CROCOIN > serial_${TEST_NAME}.log 2>&1  && { sed -e '2c Y' ${TEST_NAME}_steps > tmp.txt ; \mv tmp.txt ${TEST_NAME}_steps ;
 } || { echo -e $msg2 > /dev/stdin ; echo -e $msg1 ; exit 2 ; }

# Additional check in case of clean stop before the end
SUCCESS1=$(tail -n 2 serial_${TEST_NAME}.log | head -n 1)
SUCCESS=$(echo $SUCCESS1 | sed -e " s/\ //g")
if [ "$SUCCESS" != 'MAIN:DONE' ]; then
  echo -e "   $msg2" > /dev/stdin  
  echo -e $msg1 
  exit 2 
fi	