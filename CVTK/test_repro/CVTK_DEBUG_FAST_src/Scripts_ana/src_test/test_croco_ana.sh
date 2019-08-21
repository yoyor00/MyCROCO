#!/bin/bash

#echo
#echo '============================================== '
#echo
#echo '         TEST CROCO WITH DEBUG RVTK            '
#echo
#echo '============================================== '
#echo
##
## Get code source from jobcomp file
##
#sed -n -e '/SOURCE=/p' jobcomp_rvtk.bash > tmp1
#sed -n '$p' tmp1 > tmp2
#eval "SOURCE_CROCO=`sed -n -e '/SOURCE=/ s/.*\= *//p' tmp2`"
#rm -f tmp1 tmp2
#echo
#echo 'SOURCE_CROCO='$SOURCE_CROCO


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

SCRIPT_RVTK=rvtk_fast_qsub_ANA.bash
./$SCRIPT_RVTK > Recap_${TEST_NAME}.git${numrev}
if [ $? -gt 0 ]; then
  exit
fi  	

