#!/bin/bash

export TEST_PROJECT_DIR=$1
export TEST_DATA_DIR=$2

# prepare config
$TEST_PROJECT_DIR/create_config.bash -o all-prod-cpl -d $TEST_PROJECT_DIR/test/home -w $TEST_PROJECT_DIR/test/work -n BENGUELA_TOY -s $TEST_PROJECT_DIR

# copy needed files
cp $TEST_PROJECT_DIR/TEST_CASES/CPL_TOY/cppdefs.h $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/CROCO_IN/.
cp $TEST_PROJECT_DIR/TEST_CASES/CPL_TOY/param.h $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/CROCO_IN/.
cp $TEST_PROJECT_DIR/TEST_CASES/CPL_TOY/* $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/.
cp $TEST_PROJECT_DIR/TEST_CASES/CPL_TOY/Makefile.Linux $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/TOY_IN/.

# unpack data
cd $TEST_DATA_DIR/TEST_TOY/
tar -xf CROCO_FILES.tar.gz
tar -xf TOY_FILES.tar.gz
mv $TEST_DATA_DIR/TEST_TOY/CROCO_FILES/* $TEST_PROJECT_DIR/test/work/BENGUELA_TOY/CROCO_FILES/.
mv $TEST_DATA_DIR/TEST_TOY/TOY_FILES/* $TEST_PROJECT_DIR/test/work/BENGUELA_TOY/TOY_FILES/.

# compile CROCO
cd $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/CROCO_IN/
./jobcomp >& compile_coupled.log

# compile TOY
cd $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/TOY_IN/
./make_toy_compil.sh

# prepare INPUT FILES
cd $TEST_PROJECT_DIR/test/home/BENGUELA_TOY
source myenv_mypath.sh
$SCRIPTDIR/OASIS_SCRIPTS/create_oasis_toy_files.sh $TEST_PROJECT_DIR/test/work/BENGUELA_TOY/TOY_FILES/ww3_20050101_20050131.nc toy_wav.nc ww3 1,12
$SCRIPTDIR/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh grid_wav.nc wav.nc toy "TOY_T0M1 TOY___HS TOY__DIR "
$SCRIPTDIR/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh $TEST_PROJECT_DIR/test/work/BENGUELA_TOY/CROCO_FILES/croco_grd.nc oce.nc croco "CROCO_SSH CROCO_EOCE CROCO_NOCE "
ln -s $TEST_PROJECT_DIR/test/work/BENGUELA_TOY/CROCO_FILES .

# rename executables
mv $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/CROCO_IN/croco crocox
mv $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/TOY_IN/toy_wav toywav

# launch
chmod +x run_mpirun.sh
./run_mpirun.sh

# save for artifacts
mkdir -p $TEST_PROJECT_DIR/BENGUELA_TOY_log
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/CROCO_IN/compile_coupled.log $TEST_PROJECT_DIR/BENGUELA_TOY_log/.
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/TOY_IN/*.out $TEST_PROJECT_DIR/BENGUELA_TOY_log/.
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/debug.* $TEST_PROJECT_DIR/BENGUELA_TOY_log/.
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/nout.000000 $TEST_PROJECT_DIR/BENGUELA_TOY_log/.
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/OUTPUT_TOY_wav.txt $TEST_PROJECT_DIR/BENGUELA_TOY_log/.
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/*timers_0000 $TEST_PROJECT_DIR/BENGUELA_TOY_log/.
cp $TEST_PROJECT_DIR/test/home/BENGUELA_TOY/mpirun.log $TEST_PROJECT_DIR/BENGUELA_TOY_log/.