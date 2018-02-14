#!/bin/bash
#echo '============================================================='
echo 'Create the link between TESTROOT . and '$PWD

dir_testroot=$DATAWORK/RVTK_DEBUG_REG_DEV/TESTROOT/VORT

#echo 'Process input file'
ln -sf $dir_testroot/VORTEX .
ln -sf $dir_testroot/test_croco_vort.sh .
ln -sf $dir_testroot/rvtk_fast_qsub_VORT.bash .
ln -sf $dir_testroot/jobcomp_rvtk.bash .
ln -sf $dir_testroot/extract_results_croco.bash .

#echo 'Process namelist files'
ln -sf ${dir_home}/VORTEX/AGRIF_FixedGrids.in.VORTEX AGRIF_FixedGrids.in
ln -sf VORTEX/croco.in.vortex.1 croco.in.1
ln -sf VORTEX/croco.in.vortex croco.in

cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
