#!/bin/bash

#set -x
set -e
set -u

#echo '============================================================='
echo 'Create the link between test_perfrst/src and test dir'

source "$PERFRSTHOME/CONFIGURE_GLOBAL_PERFRST"
source "$PERFRSTHOME/CONFIGURE_REG_PERFRST"

rm -Rf $dir_test
mkdir -p $dir_test/Junk
[[ ! -d  $dir_web ]] && mkdir -p $dir_web

#echo $PERFRSTHOME
#ls -l $PERFRSTHOME/
#exit

#
\cp -rf $CI_PROJECT_DIR/OCEAN/croco.in* $PERFRSTHOME/TEST_CASES_CVTK/VHR/
\cp -rf $CI_PROJECT_DIR/OCEAN/AGRIF_FixedGrids.in $PERFRSTHOME/TEST_CASES_CVTK/VHR/

cp $PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in $PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in.write
cp $PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in $PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in.read

###################################################
# croco.in.write
###################################################

file=$PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in.write

line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
[ ! -z $toto ] && sed -e "${line} s/$toto/20/" $file > tmp.txt && \mv tmp.txt $file

# initial nrrec 1
line=$(($(grep -n 'initial:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s/$toto/1/" $file > tmp.txt && \mv tmp.txt $file

#restart NRST 3
line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s/$toto/3/" $file > tmp.txt && \mv tmp.txt $file

#restart NRPFRST 2
line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
[ ! -z $toto ] && sed -e "${line} s/$toto/2/" $file > tmp.txt && \mv tmp.txt $file

#restart / filename
line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +2))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s%$toto%CROCO_OUTFILES/croco_rst_perfrst.nc%" $file > tmp.txt && \mv tmp.txt $file

# time stepping NTIMES 6
line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

# time stepping dt 900
line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
[ ! -z $toto ] && sed -e "${line} s/$toto/900/" $file > tmp.txt && \mv tmp.txt $file  

###################################################
# croco.in.read
###################################################
file=$PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in.read

line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
[ ! -z $toto ] && sed -e "${line} s/$toto/20/" $file > tmp.txt && \mv tmp.txt $file

# initial nrrec 2
line=$(($(grep -n 'initial:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s/$toto/2/" $file > tmp.txt && \mv tmp.txt $file

# initial / filename CROCO_FILES/croco_rst.00000.nc
line=$(($(grep -n 'initial:' $file  |  awk -F ':' '{print $1}') +2))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s%$toto%CROCO_OUTFILES/croco_rst_perfrst.00000.nc%" $file > tmp.txt && \mv tmp.txt $file


#restart NRST 3
line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

#restart NRPFRST 2
line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
[ ! -z $toto ] && sed -e "${line} s/$toto/2/" $file > tmp.txt && \mv tmp.txt $file

# time stepping NTIMES 3
line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
[ ! -z $toto ] && sed -e "${line} s/$toto/3/" $file > tmp.txt && \mv tmp.txt $file

# time stepping dt 900
line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
[ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
[ ! -z $toto ] && sed -e "${line} s/$toto/900/" $file > tmp.txt && \mv tmp.txt $file  

####==========================================================
####==========================================================
# for file in $(ls $PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in.read)
# do 
#   line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
#   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
#   [ ! -z $toto ] && sed -e "${line} s/$toto/20/" $file > tmp.txt && \mv tmp.txt $file

#   line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
#   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
#   [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file
  
#   line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
#   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
#   [ ! -z $toto ] && sed -e "${line} s/$toto/12/" $file > tmp.txt && \mv tmp.txt $file

#   line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
#   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
#   [ ! -z $toto ] && sed -e "${line} s/$toto/900/" $file > tmp.txt && \mv tmp.txt $file  
# done


# # for file in $(ls $PERFRSTHOME/TEST_CASES_CVTK/VHR/croco.in.1)
# # do 
# #   line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
# #   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
# #   [ ! -z $toto ] && sed -e "${line} s/$toto/15/" $file > tmp.txt && \mv tmp.txt $file

# #   line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
# #   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
# #   [ ! -z $toto ] && sed -e "${line} s/$toto/15/" $file > tmp.txt && \mv tmp.txt $file

# #     line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
# #   [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
# #   [ ! -z $toto ] && sed -e "${line} s/$toto/300/" $file > tmp.txt && \mv tmp.txt $file
  
# # done

# # for file in $(ls $PERFRSTHOME/TEST_CASES_CVTK/VHR/AGRIF_FixedGrids.in )
# # do 
# #     #sed '2c 79 137 37 117 3 3 3 3' $file > tmp.txt && \mv tmp.txt $file
# #     sed "2c $nest_position_reg" $file > tmp.txt && \mv tmp.txt $file
# # done
####==========================================================

ln -sf $dir_home/Configure_Test_reg_perfrst $dir_test/

# configure files
ln -sf $dir_home/../CONFIGURE_GLOBAL_PERFRST $dir_test/
ln -sf $dir_home/../CONFIGURE_REG_PERFRST $dir_test/

# common scripts and programms
ln -sf $dir_home/../gather_recap_perfrst.bash $dir_test/
ln -sf $dir_home/../gitinfo.sh $dir_test/
ln -sf $dir_home/../git_process.bash $dir_test/
ln -sf $dir_home/../mk_CLEANALL.bash $dir_test/
ln -sf $dir_home/../mk_CHECKALL.bash $dir_test/
ln -sf $dir_web $dir_test/
ln -sf  $dir_home/../print/* $dir_test/

# ana specific and programms
ln -sf $dir_home/../mk_TestDIR.bash $dir_test/
ln -sf $dir_home/../mk_TESTALL.bash $dir_test/

# cleaning
rm -Rf $dir_test/Configure_Test; 
cd $dir_test	
ln -sf Configure_Test_reg_perfrst Configure_Test
cd -
#ls -l $dir_test
#
echo 'Well done: Finish linking'

