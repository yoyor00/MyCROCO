#!/bin/bash
#set -x
set -e
set -u

#echo '============================================================='
echo 'Create the link between CVTK_DEBUG_FAST_src and CVTK_DEBUG/test dir'

source "$CVTKHOME/CONFIGURE_GLOBAL"
source "$CVTKHOME/CONFIGURE_ANA"

rm -Rf $dir_test
mkdir -p $dir_test/Junk
[[ ! -d  $dir_web ]] && mkdir -p $dir_web

#
\cp -rf $CI_PROJECT_DIR/TEST_CASES/* $CVTKHOME/TEST_CASES_CVTK/.

for file in $(ls $CVTKHOME/TEST_CASES_CVTK/croco.in*)
do 
  line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/10/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/10/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/10/" $file > tmp.txt && \mv tmp.txt $file
done

ln -sf $dir_home/Configure_Test_ana $dir_test/

# configure files
ln -sf $dir_home/../CONFIGURE_GLOBAL $dir_test/
ln -sf $dir_home/../CONFIGURE_ANA $dir_test/

# common scripts and programms
ln -sf $dir_home/../gather_recap.bash $dir_test/
ln -sf $dir_home/../gitinfo.sh $dir_test/
ln -sf $dir_home/../git_process.bash $dir_test/
ln -sf $dir_home/../mk_CLEANALL.bash $dir_test/
ln -sf $dir_home/../mk_CHECKALL.bash $dir_test/
ln -sf $dir_web $dir_test/
ln -sf  $dir_home/../print/* $dir_test/

# ana specific and programms
ln -sf $dir_home/src_test/mk_TestDIR_ana.bash $dir_test/
ln -sf $dir_home/src_test/mk_TESTALL_ana.bash $dir_test/

# cleaning
rm -Rf $dir_test/Configure_Test; 
cd $dir_test	
ln -sf Configure_Test_ana Configure_Test
cd -
#
echo 'Well done: Finish linking'
