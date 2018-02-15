#!/bin/bash
#===================================
#PBS -q sequentiel
#PBS -l walltime=02:00:00
#PBS -l mem=20gb
#PBS -j oe 
#PBS -M gildas.cambon@ird.fr -m abe
cd $PBS_O_WORKDIR
echo $PBS_O_LOGNAME
#===================================
source configure_file
par1='SERIAL'

./jobcomp_rvtk.bash_v1 Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log 
mv croco croco_${par1}.exe

#Run (and produce check file)
time ./croco_${par1}.exe > serial_${TEST_NAME}.log
