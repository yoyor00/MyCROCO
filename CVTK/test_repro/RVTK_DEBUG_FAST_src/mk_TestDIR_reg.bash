#!/bin/bash

# Directory creation
echo "======================================================"
[ ! -d $1 ] && mkdir $1
echo 'Create and setup dir. : ' $1
cp Configure_Test/$1 $1/configure_file 
cp gitinfos $1
cp -Rf comp_run_* clean_rvtkdir.sh create_link_for_Test_reg.sh $1
#===========================================================
cd $1
export mytest=$1
echo " Create links "
./create_link_for_Test_reg.sh
echo " test croco"
./test_croco_reg.sh
cd -
echo "======================================================"
echo "  "

