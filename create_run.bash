#!/bin/bash
#
# Update : Apr. 2020
# G. Cambon : Sept. 2016
#
#set -x

#==========================================================================================
# BEGIN USER MODIFICATIONS
#
# Get CROCO directory
CROCO_DIR=$(cd $(dirname "$0")/..; pwd)
#
SOURCES_DIR=$(cd $(dirname "$0"); pwd)
#
TOOLS_DIR=${CROCO_DIR}/croco_tools
#
MY_CONFIG_PATH=${SOURCES_DIR}
#
# Name of the configuration directory defined by the user
#
MY_CONFIG_NAME='Run_BENGUELA_LR'
#
#
# END USER MODIFICATIONS
#==========================================================================================

while getopts :hd:s:t:n: V
do
  case $V in
    ('h') cat <<EOF
Script to setup your own croco configuration.
    What is does :
     - Copy the original croco/Run with cppdefs.h, param.h and *.in files needed
     - Copy the original crocotools_param.m and start.m file from croco_tools/
     - Copy the original run_croco*.bash file from croco_tools/Pluriannual_scripts/
    Usage:
    Use the command line:
    ./create_run.bash -d CONFIG_DIR -n CONFIG_NAME -s SOURCES_DIR -t TOOLS_DIR
    OR
    Edit the USER SECTION of the script to define the following variables :
     - SOURCES_DIR     : location of croco directory
     - TOOLS_DIR       : location of croco_tools directory
     - MY_CONFIG_PATH  : location of the repository to store the configuration
     - MY_CONFIG_NAME  : name of the configuration
EOF
    exit 0;;
    ('d')  x_d=${OPTARG};;
    ('s')  x_s=${OPTARG};;
    ('t')  x_t=${OPTARG};;
    ('n')  x_n=${OPTARG};;
  esac
done
#shift $(($OPTIND-1));

SOURCES_DIR="${x_s-$SOURCES_DIR}"
TOOLS_DIR="${x_t-$TOOLS_DIR}"
MY_CONFIG_PATH="${x_d-$MY_CONFIG_PATH}"
MY_CONFIG_NAME="${x_n-$MY_CONFIG_NAME}"

echo ""
echo "Your choices :"
echo " - SOURCES_DIR  : ${SOURCES_DIR}"
echo " - TOOLS_DIR    : ${TOOLS_DIR}"
echo " - CONFIG_DIR   : ${MY_CONFIG_PATH}"
echo " - CONFIG_NAME  : ${MY_CONFIG_NAME}"


echo -n " Do you want to proceed ? [Y/n] "
read answer
answer=`echo $answer | sed 's/^[yY].*$/y/'`
if [  -z "$answer" -o "x$answer" = "xy" ]; then
      echo " Creating configuration ..."
      echo "  "
   else
      echo " Exiting..."
      echo "  "
      exit
fi
unset -v answer

# Check if source are there
if [ ! -d $SOURCES_DIR ]; then 
	echo 'Directory for croco not found ...'
	echo 'Check the SOURCES_DIR variable ...'
	echo 'Exiting ...'
   exit 1
fi

# Check if tools are there
if [ ! -d $TOOLS_DIR ]; then 
  echo  " WARNING : croco_tools directory not found "
  echo -n " Do you want to proceed without MATLAB tools ? [Y/n] "
  read answer
  answer=`echo $answer | sed 's/^[yY].*$/y/'`
  if [ "x$answer" = "n" ]; then
#    echo " Creating configuration ..."
#    echo "  "
#  else
    echo " Exiting..."
    echo "  "
    exit
    copy_tools=1
  else
    copy_tools=0
  fi
fi
 
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
   echo "Exiting..."
    echo "  "
    exit
    fi
fi

if [[ $copy_tag == 1 ]] ; then
    
# Copy the files from croco/
    echo '=> Copy the source files from '$SOURCES_DIR
    echo '   needed to setup your own simulations'
    echo '         '
    
    cd $MY_CONFIG_PATH'/'$MY_CONFIG_NAME
    mkdir Misc TEST_CASES NAMELIST_OANALYSIS CROCO_FILES SCRATCH DATA
    
    #OCEAN
    DIRO='OCEAN'

#    cp -Rf $SOURCES_DIR/$DIRO/jobcomp .
    PAT=$(grep ^SOURCE $SOURCES_DIR/$DIRO/jobcomp)
    sed -e "s!${PAT}!SOURCE=$SOURCES_DIR/$DIRO!g" $SOURCES_DIR/$DIRO/jobcomp > jobcomp
    chmod +x jobcomp

    cp -Rf $SOURCES_DIR/$DIRO/param.h .
    cp -Rf $SOURCES_DIR/$DIRO/cppdefs.h .
    cp -Rf $SOURCES_DIR/$DIRO/croco.in .
    cp -Rf $SOURCES_DIR/$DIRO/croco.in.1 .
    cp -Rf $SOURCES_DIR/$DIRO/croco_inter.in .
    cp -Rf $SOURCES_DIR/$DIRO/sediment.in .
    cp -Rf $SOURCES_DIR/$DIRO/croco_forecast.in Misc/
    cp -Rf $SOURCES_DIR/$DIRO/croco_hindcast.in Misc/
    cp -Rf $SOURCES_DIR/$DIRO/croco_stations.in Misc/
    
    #PISCES
    cp -Rf $SOURCES_DIR/PISCES/namelist_pisces .

    # XIOS
    DIRO='XIOS'
    cp -Rf $SOURCES_DIR/$DIRO/iodef.xml .
    cp -Rf $SOURCES_DIR/$DIRO/domain_def.xml .
    cp -Rf $SOURCES_DIR/$DIRO/field_def.xml_full .
    cp -Rf $SOURCES_DIR/$DIRO/xios_launch.file .
    cp -Rf $SOURCES_DIR/$DIRO/README_XIOS .

    # TEST_CASE + NAMELIST_OANALYSIS
    cp -Rf $SOURCES_DIR/TEST_CASES .
    cp -Rf $SOURCES_DIR/SCRIPTS/NAMELIST_OANALYSIS .
    cp -Rf $SOURCES_DIR/SCRIPTS/Plurimonths_scripts/*.bash .
    
    echo '=> Copy from '$SOURCES_DIR ' done'
    echo '         '
fi

if [[ $copy_tools == 1 ]] ; then
# Copy the files from croco_tools/
# for crocotools in matlab
    echo '=> Copy the tools from and ' $TOOLS_DIR
    echo '   needed to setup your own simulations'
    echo '         '

    cp -Rf $TOOLS_DIR/start.m .
    cp -Rf $TOOLS_DIR/crocotools_param.m .
    cp -Rf $TOOLS_DIR/Town/town.dat Misc/
    #

    #
    echo '=> Copy from '$TOOLS_DIR ' done'
fi

cp -Rf $SOURCES_DIR/create_run.bash create_run.bash.BCK
cd -
