#!/bin/bash
#################################
# Gather all the log file to put them in a global one
# to have the echo on set -x
set -x

today=`date +%Y%m%d`
numrev=`echo BASIN/Recap_* | cut -d. -f2`

for testREGO in `ls Configure_Test` ; do  
    testREG=`echo $testREGO | cut -d/ -f2-`
    echo $testREG
    #echo "$testREG/Recap_${testREG}_${today}_*.git* Log_Summary"
    cp $testREG/Recap_${testREG}_${today}.git* Log_Summary/Junk
done

cd Log_Summary/Junk
for i in `ls -1 Recap_*${today}.git*` ; do 
    echo $i 
    cat $i >>  gather_recap_tmp
done
cd -

mv Log_Summary/Junk/gather_recap_tmp Log_Summary/gather_recap_${today}_${numrev}
cp Log_Summary/gather_recap_${today}_${numrev} .
rm -Rf Log_Summary/Junk/*
