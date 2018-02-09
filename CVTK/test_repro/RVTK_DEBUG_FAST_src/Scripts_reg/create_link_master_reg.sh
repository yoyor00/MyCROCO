#!/bin/bash
#echo '============================================================='
echo 'Create the link between RVTK_DEBUG_REG/TESTROOT dir . and the RVTK_DEBUG_src/'
#echo '  '
set -x

dir_home=$HOME/GIT/croco/CVTK/test_repro/RVTK_DEBUG_FAST_src/Scripts_reg
dir_web=/home/datawork-croco/public/commit-check
dir_datafiles=$DATAWORK/CROCO_FILES_VHR_BCK
dir_TESTROOT=$DATAWORK/RVTK_DEBUG_REG_DEV/TESTROOT/REG/

mkdir $dir_TESTROOT
ln -sf $dir_home/src_test/* $dir_TESTROOT/.
ln -sf $dir_home/Configure_Test_reg $dir_TESTROOT/.

# common scripts and programms
ln -sf $dir_home/../gitinfo.sh $dir_TESTROOT/.
ln -sf $dir_home/../Log_Details $dir_TESTROOT/.
ln -sf $dir_home/../git_process.bash $dir_TESTROOT/.
ln -sf $dir_home/../gather_recap.bash $dir_TESTROOT/.
ln -sf $dir_home/../jobcomp_rvtk.bash $dir_TESTROOT/.
ln -sf $dir_home/../comp_run_*.bash $dir_TESTROOT/.
ln -sf $dir_home/../extract_results_croco.bash $dir_TESTROOT/.

# input files + namlist
ln -sf $dir_home/../VHR $dir_TESTROOT/.
#ln -sf $dir_home/m $dir_TESTROOT/.

ln -sf $dir_web/Log_Summary $dir_TESTROOT/.
echo 'Well done: Finish linking'
