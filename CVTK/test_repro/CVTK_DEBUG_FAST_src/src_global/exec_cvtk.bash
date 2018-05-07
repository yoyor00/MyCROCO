#!/bin/bash

##===================================
##PBS -q sequentiel
##PBS -l walltime=02:00:00
##PBS -j oe
#cd $PBS_O_WORKDIR
#echo $PBS_O_LOGNAME
##===================================

set -x 

#-----------------------
source CONFIGURE_GLOBAL
#-----------------------

cd ${TESTROOTDIR}/REG
./mk_TESTALL_reg.bash
cd -
#==
cd ${TESTROOTDIR}/VORT
./mk_TESTALL_vort.bash
cd -
#==
cd ${TESTROOTDIR}/KTEST
./mk_TESTALL_ana.bash
cd -
#=
