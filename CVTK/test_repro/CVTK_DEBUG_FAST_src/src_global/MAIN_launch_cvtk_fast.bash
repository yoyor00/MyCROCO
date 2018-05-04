#!/bin/bash

#=========================================================
# This scripts launch the CVTK FAST procedure
# on DATARMOR for the 3 type of tests :
# - VORT  (as VORTEX)
# - KTEST (as TEST_CASES)
# - REG   (as REGIONAL)
#
# Usage : ./Lauch_CVTK_FAST.bash > CVTK_FAST.out 2>&1 
# It will launch
# O) ./cleanup_all.bash
# 1) ./setup_env_cvtk.bash
# 2) qsub exec_cvtk.bash  (if setup_env_cvtk.bash finish OK)
# 3) qsub gather_all.bash (if exec_cvtk.bash finish OK)
# The scripts are using qsub on serial queue
#
#  Gildas Cambon, LOPS (gildas.cambon@ird.fr) February 2018
#=========================================================

set -x
 
rm -f *.o????*

#- GIT pull to have the up-to date revision
./gitpull.bash

#-O Clean-up the CVTK_FAST environement for reg, ana and vort type
./cleanup_all.bash
##qsub -h -N cleanup_all cleanup_all.bash
##myjobid_clean="`qselect -N cleanup_all -u $USER`"

#-1 Set-up the CVTK_FAST environement for reg, ana and vort type
./setup_env_cvtk.bash
##qsub -N setup_env_cvtk -W depend=afterok:$myjobid_clean setup_env_cvtk.bash
##myjobid_env="`qselect -N setup_env_cvtk -u $USER`"

#-2 Lauch the various CVTK_FAST for reg, ana and vort
echo "#-2 Lauch the various CVTK_FAST for reg, ana and vort"
./exec_cvtk.bash
#qsub -h -N exec_cvtk exec_cvtk.bash
#qsub -h -N exec_cvtk exec_cvtk.bash.pbs
#myjobid_cvtk="`qselect -N exec_cvtk -u $USER`"

#-3 Get the recap for each CVTK_FAST type 
#echo "-3 Get the recap for each CVTK_FAST type" 
#qsub -N gather_all -W depend=afterany:$myjobid_cvtk gather_all.bash.pbs
#qsub -N gather_all -W depend=afterany:`qselect -N exec_cvtk` gather_all.bash.pbs

#-End 
#qrls `qselect -N exec_cvtk`

#==========
#rm -f *.o????*
