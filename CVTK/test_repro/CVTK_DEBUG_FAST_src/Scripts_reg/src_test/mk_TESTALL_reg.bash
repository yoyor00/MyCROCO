#!/bin/bash

./git_process.bash

for testconf in `ls -1 ./Configure_Test_reg/ `;do
echo $testconf
#rm -rf $testconf
mkdir $testconf
./mk_TestDIR_reg.bash $testconf
done
#

exit


