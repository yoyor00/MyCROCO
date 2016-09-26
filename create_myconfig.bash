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
    mkdir $MY_CONFIG_PATH'/'MY_CONFIG_NAME
fi
# Link the files from croco/
echo '=> Copy the file from '$SOURCES_DIR ' and ' $TOOLS_DIR 
echo '   needed to setup your own simulations'
echo '         '

cd $MY_CONFIG_PATH'/'MY_CONFIG_NAME
#
cp -Rf $SOURCES_DIR/jobcomp
cp -Rf $SOURCES_DIR/cppdef.h .
cp -Rf $SOURCES_DIR_DIR/param.h .
cp -Rf $SOURCES_DIR/*.in .

#TEST_CASES
cp -Rf $SOURCES_DIR/TEST_CASES .
cp -Rf $SOURCES_DIR/ROMS_FILES .

#NAMELIST_OANALYSIS
cp -Rf $SOURCES_DIR/$NAMELIST_OANALYSIS .

#XIOS
cp -Rf $SOURCES_DIR/domain_def.xml .
cp -Rf $SOURCES_DIR/field_def.xml_full .
cp -Rf $SOURCES_DIR/xios_launch.file .

echo '=> Copy from '$SOURCES_DIR ' done'
echo '         '
# Link the files from croco_tools/
# for romstools in matlab
cp -Rf $TOOLS_DIR/start.m .
cp -Rf $TOOLS_DIR/romstools_param.m .
ln -sf $TOOLS_DIR/*

#
echo '=> Copy from '$TOOLS_DIR ' done'

