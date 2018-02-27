#!/bin/bash

LIST_EXAMPLE=$1
ROOTDIR=$2

if [ ${#ROOTDIR} -eq 0 ]; then
  ROOTDIR=$(dirname $(dirname $PWD))
fi 
if [ ${#LIST_EXAMPLE} -eq 0 ]; then
  LIST_EXAMPLE='BASIN CANYON_A CANYON_B EQUATOR GRAV_ADJ INNERSHELF OVERFLOW SEAMOUNT SHELFRONT SOLITON UPWELLING VORTEX JET RIP  SHOREFACE SWASH THACKER TANK'
fi  

[ ! -f  jobcomp ] && \cp ${ROOTDIR}/Run/jobcomp .
[ ! -f cppdefs.h ] && \cp ${ROOTDIR}/Run/cppdefs.h .
[ ! -d TEST_CASES ] && \cp -r ${ROOTDIR}/Run/TEST_CASES .

# undef everything
sed -i .bak '1,/^#if defined REGIONAL/ s/define /undef /g' cppdefs.h 
# undef MPI
sed -i .bak "s/^#\(.*\)define\(.*\)MPI\(.*\)/#undef MPI/ " cppdefs.h
# undef OPENMP
sed -i .bak "s/^#\(.*\)define\(.*\)OPENMP\(.*\)/#undef OPENMP/ " cppdefs.h
# undef NBQ
sed -i .bak "s/^#\(.*\)define\(.*\)NBQ\(.*\)/#undef NBQ/ " cppdefs.h

# proper path in jobcomp
sed -i .bak "s:.*SOURCE=.*:SOURCE=$ROOTDIR/OCEAN:g" jobcomp


# main loop
for EXAMPLE in $LIST_EXAMPLE
  do 
# input if needed
#  if [ "$EXAMPLE" == "VORTEX" -o "$EXAMPLE" == "JET" ]; then
    if [ "$EXAMPLE" == "VORTEX" ]; then
      example=$(echo $EXAMPLE |tr '[:upper:]' '[:lower:]')
      if [ ! -f  ${example}_grd.nc ]; then
        mymake="make_${example}"
        matlab -nodesktop  -nosplash -nodisplay -r "addpath ./TEST_CASES; ${mymake};exit"
      fi
    fi
# define config and run
    sed -i .bak "s/^#\(.*\)undef\(.*\)$EXAMPLE\(.*\)/#define $EXAMPLE/" cppdefs.h
    ./jobcomp
    ./croco 
    sed -i .bak "s/#define $EXAMPLE /#undef  $EXAMPLE/" cppdefs.h  
  done

# qq sed utiles si besoin :

# selectionne jusque non inclus 
#sed -n '1,/^#if defined REGIONAL/ {/^#if defined REGIONAL/!p;} ' cppdefs.h

# selectionne jusque non inclus et undef
# sed -ne '1,/^#if defined REGIONAL/ s/define /undef /g ;  1,/^#if defined REGIONAL/ {/^#if defined REGIONAL/!p;} ; /\!/d' cppdefs.h
 
# extrait toutes les configs 
#sed  -e  '1,/^#if defined REGIONAL/ s/define /undef /g  ; /^#if defined REGIONAL/,$d ; /^$/d; /^!/d; /^\//d; /^\*/d; /^\ /d' cppdefs.h | awk -F' ' '{print $2}'
