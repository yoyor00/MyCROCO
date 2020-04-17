#!/bin/bash
#===================================

#set -x
set -eu

source CONFIGURE_GLOBAL
source configure_file

cd "$SUBMIT_DIR"
echo "$SUBMIT_DIR"

#===================================
par1='SERIAL'

# Compilation
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1  || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 1 ; }
mv croco croco_${par1}.exe

# Run (and produce check file)
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
./croco_${par1}.exe $CROCOIN > serial_${TEST_NAME}.log 2>&1  || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 2 ; }

# Additional check in case of clean stop before the end
SUCCESS=1
grep 'MAIN: DONE'  serial_${TEST_NAME}.log || SUCCESS=0

if [  "$SUCCESS" -eq 0 ]; then
  echo -e "   $msg2" | tee -a mylog.txt
  echo -e $msg1 
  exit 2 
fi	

