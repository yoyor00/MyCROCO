#!/bin/bash

cd $1 
#==
numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
cat recap_omp.log >> Recap_REGIONAL.git${numrev}
cat recap_mpi.log >> Recap_REGIONAL.git${numrev}
today=`date +%Y%m%d`
recapfile="Recap_REGIONAL_${today}.git${numrev}"
mv Recap_REGIONAL.git${numrev} $recapfile

echo "extract_results_croco.pl"
./extract_results_croco.pl "REGIONAL"

head -n 20 $recapfile > header.txt
cat Results_REGIONAL_${today}.git${numrev} >> header.txt
mv header.txt Results_REGIONAL_${today}.git${numrev}
mv Results_REGIONAL_${today}.git${numrev} Results_REGIONAL_${1}_${today}.git${numrev}
cp Results_REGIONAL_${1}_${today}.git${numrev} ./Log
#==
cd -
