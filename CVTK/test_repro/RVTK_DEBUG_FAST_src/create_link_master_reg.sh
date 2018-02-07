#!/bin/bash
#echo '============================================================='
echo 'Create the link between RVTK_DEBUG_REG/TESTROOT dir . and the RVTK_DEBUG_src/'
#echo '  '

dir_home=$HOME/GIT/croco/CVTK/test_repro/RVTK_DEBUG_FAST_src
dir_web=/home/datawork-croco/public/commit-check
dir_datafiles=$DATAWORK/CROCO_FILES_VHR_BCK
dir_TESTROOT=$DATAWORK/RVTK_DEBUG_REG_DEV/TESTROOT
# later rm -rf $dir_TESTROOT
mkdir $dir_TESTROOT

ln -sf $dir_home/mk_*.bash $dir_TESTROOT/.
ln -sf $dir_home/*bash $dir_TESTROOT/.
ln -sf $dir_home/test_croco_reg.sh $dir_TESTROOT/.
ln -sf $dir_home/create_link_for_Test_reg.sh $dir_TESTROOT/.
ln -sf $dir_home/gitinfo.sh $dir_TESTROOT/.
ln -sf $dir_home/Configure_Test $dir_TESTROOT/.
ln -sf $dir_home/Log_Details $dir_TESTROOT/.
ln -sf $dir_home/git_process.bash $dir_TESTROOT/.
ln -sf $dir_home/gather_recap.bash $dir_TESTROOT/.
ln -sf $dir_home/VHR $dir_TESTROOT/.
#ln -sf $dir_home/m $dir_TESTROOT/.

ln -sf $dir_web/Log_Summary $dir_TESTROOT/.
echo 'Well done: Finish linking'
