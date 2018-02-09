#!/bin/bash
#echo '============================================================='
echo 'Create the link between RVTK_DEBUG_REG/TESTROOT dir . and the RVTK_DEBUG_src/'
#echo '  '

dir_home=$HOME/GIT/croco/CVTK/test_repro/RVTK_DEBUG_FAST_src/Scripts_ana
dir_web=/home/datawork-croco/public/commit-check
dir_TESTROOT=$DATAWORK/RVTK_DEBUG_REG_DEV/TESTROOT/KTEST

# later rm -rf $dir_TESTROOT
mkdir $dir_TESTROOT

ln -sf $dir_home/src_test/* $dir_TESTROOT/.
#ln -sf $dir_home/src_test/*bash $dir_TESTROOT/.
#ln -sf $dir_home/src_test/test_croco_ana.sh $dir_TESTROOT/.
#ln -sf $dir_home/src_test/create_link_for_Test_ana.sh $dir_TESTROOT/.

ln -sf $dir_home/Configure_Test_ana $dir_TESTROOT/.

# common scripts and programms
ln -sf $dir_home/../gitinfo.sh $dir_TESTROOT/.
ln -sf $dir_home/../Log_Details $dir_TESTROOT/.
ln -sf $dir_home/../git_process.bash $dir_TESTROOT/.
ln -sf $dir_home/../gather_recap.bash $dir_TESTROOT/.
ln -sf $dir_home/../jobcomp_rvtk.bash $dir_TESTROOT/.
ln -sf $dir_home/../comp_run_*.bash $dir_TESTROOT/.
ln -sf $dir_home/../extract_results_croco.bash $dir_TESTROOT/.

# input files + namlist
[ ! -d $dir_TESTROOT/TEST_CASES ] && mkdir $dir_TESTROOT/TEST_CASES

ln -sf $dir_home/ANA/* $dir_TESTROOT/TEST_CASES/.
ln -sf $dir_home/IGW/* $dir_TESTROOT/TEST_CASES/.
ln -sf $dir_home/SHOREFACE/* $dir_TESTROOT/TEST_CASES/.
ln -sf $dir_home/JET/* $dir_TESTROOT/TEST_CASES/.


ln -sf $dir_web/Log_Summary $dir_TESTROOT/.
echo 'Well done: Finish linking'
