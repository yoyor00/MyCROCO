#!/bin/bash

source ../myenv_mypath.sh
ROOT_DIR=${OCE}/../
RUNDIR=`pwd`

set -x
set -e

if [[ $FC == ifort ]] ; then
    CPP1="cpp -traditional -DLinux -DIfort"
elif [[ $FC == gfortran ]] ; then
    CPP1="cpp -traditional -DLinux"
fi

[[ -d ${XIOS_NAM_DIR}.back ]] && rm -Rf ${XIOS_NAM_DIR}.back
[[ -d ${XIOS_NAM_DIR} ]] && mv ${XIOS_NAM_DIR} ${XIOS_NAM_DIR}.back
mkdir $XIOS_NAM_DIR

##copy cppdefs.h, set_global_definitions.h and cppdefs_dev.h
[[ -f cppdefs.h ]] &&  \cp -f cppdefs.h $XIOS_NAM_DIR || { echo "Missing cppdefs.h files ..." ; exit ; }


if [[ -f set_global_definitions.h ]]; then
    \cp -f set_global_definitions.h $XIOS_NAM_DIR/
else
    \cp -f ${OCE}/set_global_definitions.h $XIOS_NAM_DIR/
fi
#-   
if [[ -f cppdefs_dev.h ]]; then    
    \cp -f cppdefs_dev.h  $XIOS_NAM_DIR/
else
    \cp -f ${OCE}/cppdefs_dev.h $XIOS_NAM_DIR/
fi
##

## cpp-process the xml files (except iodef.xml) for croco
cd $XIOS_NAM_DIR/
$CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/field_def_croco.xml_full_withcpp field_def_croco.xml
$CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/file_def_croco.xml_full_withcpp file_def_croco.xml
cd -

## copy the xml files (except iodef.xml) for croco and wrf that do not nedd cpp-process
\cp -f  ${ROOT_DIR}/XIOS/context_croco.xml $XIOS_NAM_DIR
\cp -f  ${ROOT_DIR}/XIOS/domain_def_croco.xml $XIOS_NAM_DIR
#
\cp -f  ${ROOT_DIR}/XIOS/context_wrf.xml $XIOS_NAM_DIR
\cp -f  ${ROOT_DIR}/XIOS/file_def_wrf.xml $XIOS_NAM_DIR

## copy the various iodef.xml*, 
\cp -f  ${ROOT_DIR}/XIOS/iodef.xml* $XIOS_NAM_DIR
