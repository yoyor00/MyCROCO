#!/bin/bash

#==
source configure_file
numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
today=`date +%Y%m%d`

#Fill the log file
#=======================
echo '    ' >> Recap_${TEST_NAME}.git${numrev}
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" >> Recap_${TEST_NAME}.git${numrev}
echo '    ' >> Recap_${TEST_NAME}.git${numrev}

if [ ${FLAG_OPENMP} = 1 ]; then
  cat recap_omp.log >> Recap_${TEST_NAME}.git${numrev}
  cat openmp2X2_${TEST_NAME}.log >> Recap_${TEST_NAME}.git${numrev}
  echo "--------------" >> Recap_${TEST_NAME}.git${numrev}
fi


if [ ${FLAG_MPI} = 1 ]; then
  cat recap_mpi.log >> Recap_${TEST_NAME}.git${numrev}
  cat mpi2X2_${TEST_NAME}.log >> Recap_${TEST_NAME}.git${numrev}
  echo "--------------" >> Recap_${TEST_NAME}.git${numrev}
fi

cp -Rf Recap_${TEST_NAME}.git${numrev} Recap_${TEST_NAME}_${today}.git${numrev}

# Proceed the extraction
#=======================
./extract_results_croco.pl "${TEST_NAME}"
# => creation du fichier Results_${TEST_NAME}_${today}.git${numrev}

# Merge and clean the logfile
head -n 43 Recap_${TEST_NAME}_${today}.git${numrev} > header.txt
echo "===============" >> header.txt
echo "               " >> header.txt
cat Results_${TEST_NAME}_${today}.git${numrev} >> header.txt

# 
#mv header.txt Results_${TEST_NAME}_${today}.git${numrev}
#mv Results_${TEST_NAME}_${today}.git${numrev} ./Log
