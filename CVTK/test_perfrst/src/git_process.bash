#!/bin/bash

#=========================
source CONFIGURE_GLOBAL_PERFRST
#=========================
#echo
#echo "SOURCE_CROCO="$SOURCE_CROCO
#echo 

#echo "PROCESS GIT UPDATE"
#cd  $SOURCE_ROMS/..
##/usr/bin/git checkout
#/usr/bin/git pull
#cd -

# Get revision of sources
#-------------------------
./gitinfo.sh $SOURCE_CROCO > gitinfos
