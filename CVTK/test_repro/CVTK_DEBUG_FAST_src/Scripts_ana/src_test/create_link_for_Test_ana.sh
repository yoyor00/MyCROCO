#!/bin/bash
#echo '============================================================='
echo 'Create the link between TESTROOT . and '$PWD

source CONFIGURE_ANA

#echo 'Process input file'
ln -sf $dir_home/../CONFIGURE_GLOBAL .

ln -sf $dir_home/../TEST_CASES .
ln -sf $dir_home/../test_croco_ana.sh .
ln -sf $dir_home/../rvtk_fast_qsub_ANA.bash .
ln -sf $dir_home/../jobcomp_rvtk.bash .
ln -sf $dir_home/../extract_results_croco.bash .
ln -sf $dir_home/../comp_run_*.bash .

#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
