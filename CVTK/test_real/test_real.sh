#!/bin/bash

#-
#- BEGIN USER MODIFICATION
#-
# edit the name of the test cases and the name of your croco root 
# repository (/home/croco for instance) if needed#-
LIST_EXAMPLE='BASIN CANYON_A CANYON_B EQUATOR GRAV_ADJ INNERSHELF OVERFLOW SEAMOUNT SHELFRONT SOLITON UPWELLING VORTEX JET RIP  SHOREFACE SWASH THACKER TANK'
#LIST_EXAMPLE='BASIN'
ROOTDIR=$(dirname $(dirname $PWD))
#-
#- END USER MODIFICATION
#-

[ ! -f cppdefs.h ] && \cp ${ROOTDIR}/Run/cppdefs.h .
[ ! -d TESTCASES ] && \cp -r ${ROOTDIR}/Run/TEST_CASES .
[ ! -d LOG ] && mkdir LOG

# Number of cases
NB_TEST=$(echo $LIST_EXAMPLE |wc -w )

i=0
for EXAMPLE in $LIST_EXAMPLE
  do
    ((i=$i+1))

    echo '----------------------'
    echo 
    echo "RUNNING TEST Number "${i} / ${NB_TEST}" : $EXAMPLE"
    [ -f LOG/${EXAMPLE}_run.log ] && \rm LOG/${EXAMPLE}_run.log
    ./run_real.sh $EXAMPLE $ROOTDIR &> LOG/${EXAMPLE}_run.log || { echo RUNNING TEST $EXAMPLE failed... EXITING... && exit ; } 
    echo 

    echo '----------------------'
    echo 
    echo "PLOTTING TEST $EXAMPLE"
    [ -f LOG/${EXAMPLE}_plot.log ] && \rm LOG/${EXAMPLE}_plot.log
    ./plot_real.sh $EXAMPLE $ROOTDIR $i &> LOG/${EXAMPLE}_plot.log || { echo PLOTTING TEST $EXAMPLE failed... EXITING... && exit ; } 
    echo 
  done

