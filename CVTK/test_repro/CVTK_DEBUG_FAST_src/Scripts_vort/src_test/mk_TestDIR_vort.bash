#!/bin/bash
#set -x

[ ! -d gitinfos ] && ./git_process.bash

# Directory creation
echo "======================================================"
[ ! -d $1 ] && mkdir $1
echo 'Create and setup dir. : ' $1

rm -Rf Configure_Test ; ln -sf  Configure_Test_vort Configure_Test
cp Configure_Test_vort/$1 $1/configure_file 

cp gitinfos $1
cp -Rf comp_run_* create_link_for_Test_vort.sh $1
#===========================================================
cd $1
export mytest=$1
echo " Create links "
./create_link_for_Test_vort.sh
echo " test croco"
./test_croco_vort.sh
cd -
echo "======================================================"
echo "  "

