#!/bin/bash

./git_process.bash

for testconf in `ls -1 ./Configure_Test/ `;do
echo $testconf
#rm -rf $testconf
mkdir $testconf
./mk_TestDIR_reg.bash $testconf
done
#
#gather_recap.bash 

exit


