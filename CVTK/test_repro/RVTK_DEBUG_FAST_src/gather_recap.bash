#!/bin/bash
#################################
# Gather all the log file to put them in a global one
set -x

today=`date +%Y%m%d`
numrev=`echo TestREF/Recap_* | cut -d. -f2`

#for testREG in TestREF TestTIDES TestBULK TestCLIM TestPSOURCE TestAGRIF1W TestAGRIF2W ; do
for testREGO in `ls Configure_Test` ; do  
    testREG=`echo $testREGO | cut -d/ -f2-`
    echo $testREG
    #echo "$testREG/Recap_${testREG}_${today}_*.git* Log_Summary"
    cp $testREG/Recap_${testREG}_${today}.git* Log_Summary
done

cd Log_Summary/
for i in `ls -1 Recap_*${today}.git*` ; do 
    echo $i 
    cat $i >>  gather_recap_tmp
done
cd -
mv Log_Summary/gather_recap_tmp ./gather_recap_${today}_${numrev}
rm -Rf Log_Summary/Recap*_.git

#cd Log_Summary
#cat Recap_${testREG}_${today}.git* >> gather_recap_${today}_${numrev}
#cp -Rf Log_Summary/gather_recap_${today}_${numrev} .
