#!/bin/bash
#echo '============================================================='
echo 'Create the link between TESTROOT . and '$PWD

dir_testroot=$DATAWORK/RVTK_DEBUG_REG_DEV/TESTROOT/KTEST

#echo 'Process input file'
ln -sf $dir_testroot/TEST_CASES .

ln -sf $dir_testroot/test_croco_ana.sh .
ln -sf $dir_testroot/rvtk_fast_qsub_ANA.bash .
ln -sf $dir_testroot/jobcomp_rvtk.bash .
ln -sf $dir_testroot/extract_results_croco.bash .

#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
