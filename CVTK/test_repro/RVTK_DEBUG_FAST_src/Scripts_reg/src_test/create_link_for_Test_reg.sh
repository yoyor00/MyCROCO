#!/bin/bash
#echo '============================================================='
echo 'Create the link between TESTROOT . and '$PWD
set -x
dir_testroot=/home7/datawork/crocodev/RVTK_DEBUG_REG_DEV/TESTROOT/REG
dir_datafiles=$DATAWORK/CROCO_FILES_VHR_BCK

#echo 'Process input file'
if [ ! -d CROCO_FILES ] ; then 
ln -sf $dir_datafiles CROCO_FILES
fi
ln -sf $dir_testroot/VHR .
ln -sf $dir_testroot/test_croco_reg.sh .
ln -sf $dir_testroot/rvtk_fast_qsub_REGIONAL.bash .
ln -sf $dir_testroot/jobcomp_rvtk.bash .
ln -sf $dir_testroot/extract_results_croco.bash .

#echo 'Process namelist files'
ln -sf ${dir_home}/VHR/AGRIF_FixedGrids.in.REGIONAL.VHR AGRIF_FixedGrids.in.REGIONAL
ln -sf VHR/croco.in.VHR.1 croco.in.1
ln -sf VHR/croco.in.VHR croco.in

cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
