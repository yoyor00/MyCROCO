#!/bin/bash
#################################
# Gather all the log file to put them in a global one

today=`date +%Y%m%d`
for test in TestREF TestTIDES TestBULK TestCLIM TestPSOURCE TestAGRIF1W TestAGRIF2W ; do
    echo $test
    cp $test/Recap_${test}_*.git* Log2/
    cd Log2 
    cat Recap_${test}_*.git* >> gather_recap_$today
    cd -
done

