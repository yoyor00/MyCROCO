#!/bin/bash

./git_process.bash

for testconf in `ls -1 ./Configure_Test_vort/ `;do
echo $testconf
#rm -rf $testconf
mkdir $testconf
./mk_TestDIR_vort.bash $testconf
done
#

exit


