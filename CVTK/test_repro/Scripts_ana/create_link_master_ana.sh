#!/bin/bash

##set -x
##set -e
set -u

#echo '============================================================='
echo 'Create the link between test_repro sources and test dir'

source "$CVTKHOME/CONFIGURE_ANA"

rm -Rf $dir_test
mkdir -p $dir_test/Junk
[[ ! -d  $dir_web ]] && mkdir -p $dir_web

#
\cp -rf $CI_PROJECT_DIR/TEST_CASES/* $CVTKHOME/../common/TEST_CASES_CVTK/.

for file in $(ls $CVTKHOME/../common/TEST_CASES_CVTK/croco.in*)
do 
  line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file
done
 
# CI common scripts and programms
ln -sf $CVTKHOME/../common/CONFIGURE_GLOBAL $dir_test/
ln -sf $CVTKHOME/../common/gitinfo.sh $dir_test/
ln -sf $CVTKHOME/../common/git_process.bash $dir_test/
ln -sf $CVTKHOME/../common/mk_CLEANALL.bash $dir_test/
ln -sf $CVTKHOME/../common/mk_CHECKALL.bash $dir_test/
ln -sf $CVTKHOME/../common/jobcomp_rvtk.bash $dir_test/
ln -sf $CVTKHOME/../common/print/* $dir_test/

# test repro specific + ana specific
ln -sf $dir_web $dir_test/
ln -sf $dir_home/Configure_Test_ana $dir_test/
ln -sf $dir_home/../CONFIGURE_ANA $dir_test/
ln -sf $dir_home/../mk_TestDIR.bash $dir_test/
ln -sf $dir_home/../mk_TESTALL.bash $dir_test/
ln -sf $dir_home/../gather_recap.bash $dir_test/

# cleaning
rm -Rf $dir_test/Configure_Test; 
cd $dir_test	
ln -sf Configure_Test_ana Configure_Test
cd -
#
echo 'Well done: Finish linking'

