#!/bin/bash

###########################
source configure_file
###########################

# Get revision of sources
#-------------------------
numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
#echo
#echo  "Testing CROCO Rev$numrev"

# Define RVTK
#------------
rm -f Recap_*
touch Recap_${TEST_NAME}.git${numrev}
rm -f mylog.txt
touch mylog.txt

SCRIPT_RVTK=rvtk_fast_qsub_ANA.bash
echo -e "   - Run Tests" 
./$SCRIPT_RVTK > Recap_${TEST_NAME}.git${numrev}
cat mylog.txt
if [ $? -gt 0 ]; then
  exit
fi  	

