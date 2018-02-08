#!/bin/bash

# Directory creation
echo "======================================================"
[ ! -d $1 ] && mkdir $1
echo 'Create and setup dir. : ' $1
cp Configure_Test_ana/$1 $1/configure_file 
cp gitinfos $1
cp -Rf comp_run_* create_link_for_Test_ana.sh $1
#===========================================================
cd $1
export mytest=$1
echo " Create links "
./create_link_for_Test_ana.sh
echo " test croco"
./test_croco_ana.sh
cd -
echo "======================================================"
echo "  "

