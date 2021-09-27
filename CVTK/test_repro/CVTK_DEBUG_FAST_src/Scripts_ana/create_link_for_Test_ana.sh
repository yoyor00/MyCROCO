#!/bin/bash
#echo '============================================================='
#set -x
#set -e
echo 'Create the link between TESTROOT . and '$PWD

source ../CONFIGURE_ANA

#echo 'common scripts
ln -sf $dir_home/../CONFIGURE_GLOBAL .
ln -sf $dir_home/../CONFIGURE_ANA .

[ -d TEST_CASES ] && rm TEST_CASES
ln -sf $dir_home/../../../common/TEST_CASES_CVTK TEST_CASES
ln -sf $dir_home/../../../common/TEST_CASES_CVTK/*.nc* .
ln -sf $dir_home/../jobcomp_rvtk.bash .
ln -sf $dir_home/../extract_results_croco.bash .
ln -sf $dir_home/../comp_run_*.bash .




# specific scripts
ln -sf $dir_home/../test_croco.sh .
ln -sf $dir_home/../rvtk_fast_qsub.bash .
#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
