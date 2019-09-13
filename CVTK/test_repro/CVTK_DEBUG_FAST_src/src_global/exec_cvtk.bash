#!/bin/bash

#set -x 

#-----------------------
source CONFIGURE_GLOBAL
#-----------------------

cd ${TESTROOTDIR}/REG
./mk_TESTALL_reg.bash
cd -

#==
#cd ${TESTROOTDIR}/VORT
#./mk_TESTALL_vort.bash
#cd -

#==
#cd ${TESTROOTDIR}/KTEST
#./mk_TESTALL_ana.bash
#cd -
#=
