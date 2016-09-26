#!/bin/bash

# Scripts to setup your own croco configuration.
# 
# What is does : 
#   - Copy the original croco/Run with cppdefs.h, param.h and *.in files needed
#   - Copy the original romstools_param.m and start.m file from croco_tools/
#   - Copy the original run_roms*.bash file from croco_tools/Pluriannual_scripts/
#   - 
#
# G. Cambon : Sept. 2016
#-----------------------------------------------------------------------------------

SOURCES_DIR='/home/gcambon/GIT_utd/croco'
TOOLS_DIR='/home/gcambon/GIT_utd/croco_tools'

MY_CONFIG_PATH='/local/tmp/2/gcambon/CONFIGS/'
MY_CONFIG_NAME='BENGUELA_LR'

#
# END USER SECTION
#==========================================================================================
# Create the directory
if [ ! -d $MY_CONFIG_PATH ]; then 
    mkdir $MY_CONFIG_PATH
    if [ ! -d $MY_CONFIG_PATH'/'$MY_CONFIG_NAME ]; then
	mkdir $MY_CONFIG_PATH'/'$MY_CONFIG_NAME
	copy_tag=1
    else
	echo 'Already a configuration exists ...'
	echo 'Check the configuration directory ' $MY_CONFIG_PATH'/'$MY_CONFIG_NAM
	copy_tag=0
    fi
else
    if [ ! -d $MY_CONFIG_PATH'/'$MY_CONFIG_NAME ]; then
	mkdir $MY_CONFIG_PATH'/'$MY_CONFIG_NAME
	copy_tag=1
    else
	echo 'Already a configuration exists ...'
	echo 'Check the configuration directory ' $MY_CONFIG_PATH'/'$MY_CONFIG_NAM
	copy_tag=0
    fi
fi

if [[ $copy_tag == 1 ]] ; then
    
# Copy the files from croco/
    echo '=> Copy the file from '$SOURCES_DIR ' and ' $TOOLS_DIR 
    echo '   needed to setup your own simulations'
    echo '         '
    
    
    cd $MY_CONFIG_PATH'/'$MY_CONFIG_NAME
    
    
    #OCEAN
    cp -Rf $SOURCES_DIR/OCEAN/jobcomp .
    cp -Rf $SOURCES_DIR/OCEAN/param.h . 
    cp -Rf $SOURCES_DIR/OCEAN/cppdefs.h . 
    cp -Rf $SOURCES_DIR/OCEAN/roms.in . 
    cp -Rf $SOURCES_DIR/OCEAN/roms.in.1 . 
    cp -Rf $SOURCES_DIR/OCEAN/roms_inter.in .
    cp -Rf $SOURCES_DIR/OCEAN/sediment.in .
    cp -Rf $SOURCES_DIR/OCEAN/roms_forecast.in .
    cp -Rf $SOURCES_DIR/OCEAN/roms_hindcast.in .
    cp -Rf $SOURCES_DIR/OCEAN/roms_stations.in .
    
    # TEST_CASE + NAMELIST_OANALYSIS
    cp -Rf $SOURCES_DIR/OCEAN/TEST_CASES .
    cp -Rf $SOURCES_DIR/OCEAN/NAMELIST_OANALYSIS .
    
    # XIOS
    cp -Rf $SOURCES_DIR/XIOS/domain_def.xml .
    cp -Rf $SOURCES_DIR/XIOS/field_def.xml_full .
    cp -Rf $SOURCES_DIR/XIOS/xios_launch.file .
    cp -Rf $SOURCES_DIR/XIOS/README_XIOS .
 
    echo '=> Copy from '$SOURCES_DIR ' done'
    echo '         '
    # Link the files from croco_tools/
    # for romstools in matlab
    cp -Rf $TOOLS_DIR/start.m .
    cp -Rf $TOOLS_DIR/romstools_param.m .
    cp -Rf $TOOLS_DIR/Misc/town.dat .
    #
    cd -
    #
    echo '=> Copy from '$TOOLS_DIR ' done'
    
fi
