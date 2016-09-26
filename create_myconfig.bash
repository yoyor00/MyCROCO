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
    mkdir 
    
    cd $MY_CONFIG_PATH'/'$MY_CONFIG_NAME
    mkdir OCEAN XIOS Misc TEST_CASES NAMELIST_OANALYSIS
    
    #OCEAN
    DIRO='OCEAN'
    cp -Rf $SOURCES_DIR/$DIRO/jobcomp $./
    cp -Rf $SOURCES_DIR/$DIRO/param.h $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/cppdefs.h $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/roms.in $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/roms.in.1 $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/roms_inter.in $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/sediment.in $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/roms_forecast.in $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/roms_hindcast.in $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/roms_stations.in $DIRO/
    
    # TEST_CASE + NAMELIST_OANALYSIS
    cp -Rf $SOURCES_DIR/Run/TEST_CASES TEST_CASES
    cp -Rf $SOURCES_DIR/Run/NAMELIST_OANALYSIS NAMELIST_OANALYSIS/
    cp -Rf $SOURCES_DIR/Run/Plurimonths_scripts .
    
    # XIOS
    DIRO='XIOS'
    cp -Rf $SOURCES_DIR/$DIRO/domain_def.xml $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/field_def.xml_full $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/xios_launch.file $DIRO/
    cp -Rf $SOURCES_DIR/$DIRO/README_$DIRO $DIRO/
 
    echo '=> Copy from '$SOURCES_DIR ' done'
    echo '         '
    # Link the files from croco_tools/
    # for romstools in matlab
    cp -Rf $TOOLS_DIR/start.m .
    cp -Rf $TOOLS_DIR/romstools_param.m .
    cp -Rf $TOOLS_DIR/Misc/town.dat Misc/
    #
    cd -
	 
    cp create_myconfig.bash $MY_CONFIG_PATH'/'$MY_CONFIG_NAME

    #
    echo '=> Copy from '$TOOLS_DIR ' done'
    
fi
