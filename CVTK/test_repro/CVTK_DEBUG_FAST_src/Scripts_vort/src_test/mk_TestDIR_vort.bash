#!/bin/bash
#set -x
set -e
set -u

source CONFIGURE_GLOBAL
source CONFIGURE_VORT

echo "   - SOURCE_CROCO="$SOURCE_CROCO

[ ! -d gitinfos ] && ./git_process.bash

numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
echo  "   - Testing CROCO Rev$numrev"
# Directory creation
#echo "======================================================"
[ ! -d $1 ] && mkdir $1
echo '   - Create and setup dir. : ' $1

rm -Rf Configure_Test ; ln -sf  Configure_Test_vort Configure_Test

cp Configure_Test_vort/$1 $1/configure_file 
cp gitinfos $1
ln -sf $dir_home/src_test/create_link_for_Test_vort.sh $1
#===========================================================
cd $1 ;  
export mytest=$1
./create_link_for_Test_vort.sh
./test_croco_vort.sh
cd -

#echo "======================================================"
#echo "  "

