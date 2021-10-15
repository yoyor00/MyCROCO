#!/bin/bash
#
# Update S. Jullien : Oct 2021
# Update : Apr. 2020
# G. Cambon : Sept. 2016
#
#set -x

#==========================================================================================
# BEGIN USER MODIFICATIONS

# Machine you are working on
# Known machines: Linux DATARMOR IRENE JEANZAY
# ---------------------------------------------
export MACHINE="Linux"

# CROCO parent directory 
# (where croco_tools directory and croco source directory can be found)
# -----------------
export CROCO_DIR=${HOME}/croco

# Configuration name
# ------------------
export MY_CONFIG_NAME=BENGUELA

# Home and Work configuration directories
# ---------------------------------------
export MY_CONFIG_HOME=${HOME}/CONFIGS/${MY_CONFIG_NAME}
export MY_CONFIG_WORK=${HOME}/CONFIGS/${MY_CONFIG_NAME}

# Options of your configuration
# Known options: 
# oce        : croco
# prepro     : for getting scripts for CROCO preprocessing
# inter      : for running interannual runs 
# forc       : for using forecast scripts
# test_cases : for running test cases 
# xios       : xios server
# pisces     : pisces inputs
# agrif      : inputs for nests
# sediment   : inputs for sediment 
# mustang    : mustang model
# cpl        : for coupling with OASIS
# atm wav toy: other models for coupling
# -----------------------------------------------------
models=( "oce" "inter" )

# END USER MODIFICATIONS
#==========================================================================================

x_f=0

while getopts :hfd:w:s:n:o: V
do
  case $V in
    ('h') cat <<EOF
Script to setup your own configuration.
    What is does :
     - Copy the original cppdefs.h, param.h and *.in files needed
     - Copy the original crocotools_param.m and start.m file from croco_tools/
     - Copy the original run_croco*.bash file from croco/SCRIPTS/Plurimonths_scripts/
    Usage:
    Use the command line:
    ./create_config.bash -d MY_CONFIG_HOME -w MY_CONFIG_WORK -n MY_CONFIG_NAME -s CROCO_DIR -o OPTS
    OR
    Edit the USER SECTION of the script to define the following variables :
     - CROCO_DIR     : location of croco parent directory containing croco and croco_tools
     - MY_CONFIG_HOME  : location of the repository to store the configuration
     - MY_CONFIG_WORK  : location of the repository to store the configuration large input files, and where it will be run
     - MY_CONFIG_NAME  : name of the configuration
     - OPTS            : options of your configuration, see the known keywords
EOF
    exit 0;;
    ('f')  x_f=1;;
    ('d')  x_d=${OPTARG};;
    ('w')  x_d=${OPTARG};;
    ('s')  x_s=${OPTARG};;
    ('n')  x_n=${OPTARG};;
    ('o')  x_o=${OPTARG};;
  esac
done
#shift $(($OPTIND-1));

CROCO_DIR="${x_s-$CROCO_DIR}"
MY_CONFIG_HOME="${x_d-$MY_CONFIG_HOME}"
MY_CONFIG_WORK="${x_d-$MY_CONFIG_WORK}"
MY_CONFIG_NAME="${x_n-$MY_CONFIG_NAME}"
models="${x_o-$models}"

echo ""
echo "Your choices :"
echo " - CROCO_DIR  : ${CROCO_DIR}"
echo " - CONFIG_HOME_DIR   : ${MY_CONFIG_HOME}"
echo " - CONFIG_WORK_DIR   : ${MY_CONFIG_WORK}"
echo " - CONFIG_NAME  : ${MY_CONFIG_NAME}"
echo " - OPTIONS      : ${models[@]}"

if [ $x_f -eq 0 ]; then
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
fi

SOURCES_DIR=$CROCO_DIR/croco
TOOLS_DIR=$CROCO_DIR/croco_tools

# Check if source are there
if [ ! -d $SOURCES_DIR ]; then 
	echo 'Directory for croco not found ...'
	echo 'Check the SOURCES_DIR variable ...'
	echo 'Exiting ...'
   exit 1
fi

# Check if tools are there
copy_tools=1
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
  fi
fi
 
# Create the directory
if [ ! -d $MY_CONFIG_HOME ]; then 
    mkdir -p $MY_CONFIG_HOME
else
    echo 'Already a configuration exists ...'
    echo 'You should check the configuration directory ' $MY_CONFIG_HOME
    echo -n " Do you want to proceed anyway (risk of overwriting) ? [N/y] "
    read answer
    answer=`echo $answer | sed 's/^[nN].*$/n/'`
    if [  -z "$answer" -o "x$answer" = "xn" ]; then
        echo " Exiting..."
        echo "  "
        exit
    else
        echo " Proceed..."
        echo "  "
    fi
fi

if [ "$MY_CONFIG_WORK" != "$MY_CONFIG_HOME" ]; then
    if [ ! -d $MY_CONFIG_WORK ]; then
        mkdir -p $MY_CONFIG_WORK
    else
        echo 'Already a configuration exists ...'
        echo 'You chould check the configuration directory ' $MY_CONFIG_WORK
        echo -n " Do you want to proceed anyway (risk of overwriting) ? [N/y] "
        read answer
        answer=`echo $answer | sed 's/^[nN].*$/n/'`
        if [  -z "$answer" -o "x$answer" = "xn" ]; then
            echo " Exiting..."
            echo "  "
            exit
        else
            echo " Proceed..."
            echo "  "
        fi
    fi
fi

cp create_config.bash $MY_CONFIG_HOME/create_config.bash.bck

if [[ ${models[@]} =~ "oce" ]] ; then
    echo 'Copy CROCO useful scripts and input files'
    echo '-----------------------------------------'
    # Create directories
    mkdir -p $MY_CONFIG_HOME/CROCO_IN
    mkdir -p $MY_CONFIG_WORK/CROCO_FILES
    # CROCO general
    cp -f $SOURCES_DIR/OCEAN/cppdefs.h $MY_CONFIG_HOME/CROCO_IN/.
    cp -f $SOURCES_DIR/OCEAN/param.h $MY_CONFIG_HOME/CROCO_IN/.
    cp -f $SOURCES_DIR/OCEAN/jobcomp $MY_CONFIG_HOME/CROCO_IN/.
    cp -f $SOURCES_DIR/OCEAN/croco.in $MY_CONFIG_HOME/CROCO_IN/.
    cp -f $SOURCES_DIR/OCEAN/croco_stations.in $MY_CONFIG_HOME/CROCO_IN/.
    # TEST_CASES
    if [[ ${models[@]} =~ "test_cases" ]] ; then
      cp -Rf $SOURCES_DIR/TEST_CASES $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # AGRIF
    if [[ ${models[@]} =~ "agrif" ]] ; then
      cp -f $SOURCES_DIR/OCEAN/croco.in.1 $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # INTER
    if [[ ${models[@]} =~ "inter" ]] ; then
      cp -f $SOURCES_DIR/OCEAN/croco_inter.in $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # FORECAST
    if [[ ${models[@]} =~ "forc" ]] ; then
      cp -f $SOURCES_DIR/OCEAN/croco_forecast.in $MY_CONFIG_HOME/CROCO_IN/.
      cp -f $SOURCES_DIR/OCEAN/croco_hindcast.in $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # PISCES
    if [[ ${models[@]} =~ "pisces" ]] ; then
      cp -f $SOURCES_DIR/PISCES/*namelist* $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # SEDIMENT
    if [[ ${models[@]} =~ "sediment" ]] ; then
      cp -f $SOURCES_DIR/OCEAN/sediment.in $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # MUSTANG
    if [[ ${models[@]} =~ "mustang" ]] ; then
      mkdir -p $MY_CONFIG_HOME/CROCO_IN/FIC_NAMELIST
      cp -f $SOURCES_DIR/MUSTANG/NAM_CASES/*txt $MY_CONFIG_HOME/CROCO_IN/FIC_NAMELIST/.
    fi
    # OANALYSIS
    if [[ ${models[@]} =~ "oanalysis" ]] ; then
       cp -Rf $SOURCES_DIR/SCRIPTS/NAMELIST_OANALYSIS $MY_CONFIG_HOME/CROCO_IN/.
    fi
   # XIOS
    if [[ ${models[@]} =~ "xios" ]] ; then
     mkdir -p $MY_CONFIG_HOME/XIOS_IN
     cp -Rf $SOURCES_DIR/XIOS/iodef.xml $MY_CONFIG_HOME/XIOS_IN/.
     cp -Rf $SOURCES_DIR/XIOS/domain_def.xml $MY_CONFIG_HOME/XIOS_IN/.
     cp -Rf $SOURCES_DIR/XIOS/field_def.xml_full $MY_CONFIG_HOME/XIOS_IN/.
     cp -Rf $SOURCES_DIR/XIOS/xios_launch.file $MY_CONFIG_HOME/XIOS_IN/.
     cp -Rf $SOURCES_DIR/XIOS/README_XIOS $MY_CONFIG_HOME/XIOS_IN/.
    fi
    # PREPROCESSING
    if [[ ${models[@]} =~ "prepro" ]] ; then
       cp -Rf $TOOLS_DIR/start.m $MY_CONFIG_HOME/CROCO_IN/.
       cp -Rf $TOOLS_DIR/oct_start.m $MY_CONFIG_HOME/CROCO_IN/.
       cp -Rf $TOOLS_DIR/crocotools_param.m $MY_CONFIG_HOME/CROCO_IN/.
       cp -Rf $TOOLS_DIR/Town/town.dat $MY_CONFIG_HOME/CROCO_IN/.
    fi
    # SCRIPTS FOR RUNNING
    if [[ ${models[@]} =~ "inter" ]] ; then
       cp -Rf $SOURCES_DIR/SCRIPTS/Plurimonths_scripts/*.bash $MY_CONFIG_HOME/
    fi
fi

### Coupling and other models to be coupled with ###

# OASIS
if [[ ${models[@]} =~ "cpl" ]] ; then
    echo 'Copy OASIS useful scripts and input files'
    echo '-----------------------------------------'
    mkdir -p $MY_CONFIG_HOME/OASIS_IN
    cp -r $SOURCES_DIR/SCRIPTS/SCRIPTS_COUPLING/OASIS_IN/* $MY_CONFIG_HOME/OASIS_IN/.
    if [[ ${models[@]} =~ "oce" ]] ; then
      cp -r $SOURCES_DIR/SCRIPTS/SCRIPTS_COUPLING/CROCO_IN/* $MY_CONFIG_HOME/CROCO_IN/.
    fi
    if [[ ${models[@]} =~ "prepro" ]] ; then
      mkdir -p $MY_CONFIG_HOME/PREPRO
      cp -r $TOOLS_DIR/Coupling_tools/* $MY_CONFIG_HOME/PREPRO/.
    fi
fi

# WW3
if [[ ${models[@]} =~ "wav" ]] ; then
    echo 'Copy WW3 useful scripts and input files'
    echo '-----------------------------------------'
    mkdir -p $MY_CONFIG_HOME/WW3_IN
    mkdir -p $MY_CONFIG_WORK/WW3_FILES
    cp -r $SOURCES_DIR/SCRIPTS/SCRIPTS_COUPLING/WW3_IN/* $MY_CONFIG_HOME/WW3_IN/.
fi

# WRF
if [[ ${models[@]} =~ "atm" ]] ; then
    echo 'Copy WRF useful scripts and input files'
    echo '-----------------------------------------'
    mkdir -p $MY_CONFIG_HOME/WRF_IN
    mkdir -p $MY_CONFIG_WORK/WRF_FILES
    cp -r $SOURCES_DIR/SCRIPTS/SCRIPTS_COUPLING/WRF_IN/* $MY_CONFIG_HOME/WRF_IN/.
fi

# TOY
if [[ ${models[@]} =~ "toy" ]] ; then
    echo 'Copy TOY sources, useful scripts and input files'
    echo '------------------------------------------------'
    mkdir -p $MY_CONFIG_HOME/TOY_IN
    mkdir -p $MY_CONFIG_WORK/TOY_FILES
    cp -r $SOURCES_DIR/SCRIPTS/SCRIPTS_COUPLING/TOY_IN/* $MY_CONFIG_HOME/TOY_IN/.
fi

# Coupling scripts
if [[ ${models[@]} =~ "cpl" ]] || [[ ${models[@]} =~ "wav" ]] || [[ ${models[@]} =~ "atm" ]] || [[ ${models[@]} =~ "toy" ]] ; then
    echo 'Copy scripts for coupled runs'
    echo '-----------------------------'
    [ -d $MY_CONFIG_HOME/ROUTINES ] && \rm -Rf $MY_CONFIG_HOME/ROUTINES
    cp -Rf $SOURCES_DIR/SCRIPTS/SCRIPTS_COUPLING/SCRIPTS_TOOLBOX/* $MY_CONFIG_HOME/

    # Edit myjob.sh to add CPU lines for each model
    cd $MY_CONFIG_HOME/
    [ -f myjob.tmp ] && rm -Rf myjob.tmp
    [[ ${models[@]} =~ "oce" ]] && printf "export NP_OCEX=2 \nexport NP_OCEY=2\n" >> myjob.tmp
    [[ ${models[@]} =~ "wav" ]] && printf "export NP_WAV=14 \n" >> myjob.tmp
    [[ ${models[@]} =~ "atm" ]] && printf "export NP_ATM=12 \n" >> myjob.tmp
    [[ ${models[@]} =~ "toy" ]] && printf "export NP_TOY=2 \n" >> myjob.tmp
    [[ ${models[@]} =~ "xios" ]] && printf "export NP_XIOS_ATM=1\nexport NP_XIOS_OCE=1\n" >> myjob.tmp

    if [[ ${models[@]} =~ "atm" ]] ; then
        printf "\n# additional MPI Settings for ATM (WRF)\n" >> myjob.tmp
        printf "export atm_nprocX=-1      # -1 for automatic settings\n" >> myjob.tmp
        printf "export atm_nprocY=-1      # -1 for automatic settings\n" >> myjob.tmp
        printf "export atm_niotaskpg=0    # 0 for default settings\n" >> myjob.tmp
        printf "export atm_niogp=1        # 1 for default settings\n\n" >> myjob.tmp
    fi

    sed -e "/< insert here CPU >/r myjob.tmp" \
        myjob.sh > myjob_tmp
    mv myjob_tmp myjob.sh
    chmod 755 myjob.sh
    rm -Rf myjob.tmp

    # Create the path file
    cd $MY_CONFIG_HOME/ROUTINES/PATHS
    cat ./path_base.sh >> tmppath

    # add sections for each model
    for k in `seq 0 $(( ${#models[@]} - 1))` ; do
        [[ ${models[$k]} == "cpl" ]] && printf "export CPL=\"\"\n" >> tmppath
        [[ ${models[$k]} == "oce" ]] && printf "export OCE=\"${CROCO_DIR}/croco/OCEAN\"\n" >> tmppath
        [[ ${models[$k]} == "atm" ]] && printf "export ATM=\"\"\n" >> tmppath
        [[ ${models[$k]} == "wav" ]] && printf "export WAV=\"\"\n" >> tmppath
        [[ ${models[$k]} == "toy" ]] && printf "export TOY=\"\${CHOME}/TOY_IN\"\n" >> tmppath
        [[ ${models[$k]} == "xios" ]] && printf "export XIOS=\"\"\n" >> tmppath
    done
    for k in ${models[@]} ; do
        [ -f path_${k}.sh ] && cat ./path_${k}.sh >> tmppath
    done

    # replace environment variables in path file
    sed -e "s|export MACHINE=.*|export MACHINE=\"${MACHINE}\"|g" \
        -e "s|export CONFIG=.*|export CONFIG=${MY_CONFIG_NAME}|g" \
        -e "s|export CHOME=.*|export CHOME=${MY_CONFIG_HOME}|g" \
        -e "s|export CWORK=.*|export CWORK=${MY_CONFIG_WORK}|g" \
        tmppath > tmppath1
    mv tmppath1 tmppath
    mv tmppath ${MY_CONFIG_HOME}/

    # Create the env file
    [ -d ${MY_CONFIG_HOME}/ROUTINES/MACHINE/${MACHINE} ] && cd ${MY_CONFIG_HOME}/ROUTINES/MACHINE/${MACHINE} || { echo "No environement for ${MACHINE} in ${MY_CONFIG_HOME}/SCRIPT_CPL/ROUTINES/MACHINE/${MACHINE}"; exit ;}
    cp myenv.${MACHINE} tmpenv
    for k in `seq 0 $(( ${#models[@]} - 1))` ; do
        [[ ${models[$k]} == "atm" ]] && cat ./myenv.${MACHINE}.wrf >> tmpenv 
        [[ ${models[$k]} == "wav" ]] && cat ./myenv.${MACHINE}.ww3 >> tmpenv
    done
    mv tmpenv ${MY_CONFIG_HOME}/
    cd ${MY_CONFIG_HOME}

    # concatenate path and env file
    cat tmpenv tmppath > myenv_mypath.sh
    chmod 755 myenv_mypath.sh
    rm -rf tmppath tmpenv

    # Create the namelist file
    cd ${MY_CONFIG_HOME}/ROUTINES/NAMELISTS
    cp namelist_head.sh mynamelist.sh

    if [[ ${models[@]} =~ "cpl" ]]; then
        if [[ ${models[@]} =~ "oce" ]] && [[ ${models[@]} =~ "wav" ]] && [[ ${models[@]} =~ "atm" ]] ; then
            printf "export RUNtype=owa\n#\n" >> mynamelist.sh
        elif [[ ${models[@]} =~ "oce" ]] && [[ ${models[@]} =~ "wav" ]] ; then
            printf "export RUNtype=ow\n#\n" >> mynamelist.sh
        elif [[ ${models[@]} =~ "oce" ]] && [[ ${models[@]} =~ "atm" ]]; then
            printf "export RUNtype=oa\n#\n" >> mynamelist.sh
        elif [[ ${models[@]} =~ "wav" ]] && [[ ${models[@]} =~ "atm" ]]; then
            printf "export RUNtype=aw\n#\n" >> mynamelist.sh
        elif [[ ${models[@]} =~ "toy" ]]; then
            printf "export RUNtype=Put the type here (ow/oa/aw/owa)\n#\n" >> mynamelist.sh
        else 
            printf "export RUNtype=frc\n#\n" >> mynamelist.sh
        fi
    else
        printf "export RUNtype=frc\n#\n" >> mynamelist.sh
    fi

    if [[ ${models[@]} =~ "atm" ]]; then
        printf "export USE_ATM=1\n" >> mynamelist.sh
        [[ ${models[@]} =~ "xios" ]] && printf "export USE_XIOS_ATM=0\n" >> mynamelist.sh
    fi
    if [[ ${models[@]} =~ "oce" ]]; then
        printf "export USE_OCE=1\n" >> mynamelist.sh
        [[ ${models[@]} =~ "xios" ]] && printf "export USE_XIOS_OCE=0\n" >> mynamelist.sh
    fi
    if [[ ${models[@]} =~ "wav" ]]; then
        printf "export USE_WAV=1\n" >> mynamelist.sh
    fi
    if [[ ${models[@]} =~ "toy" ]]; then
        cat ./namelist_head_toy.sh >> mynamelist.sh
    fi

    cat ./namelist_rundir.sh >> mynamelist.sh

    for k in `seq 0 $(( ${#models[@]} - 1))` ; do
        [[ ${models[$k]} == "oce" ]] && printf "export OCE_EXE_DIR=${MY_CONFIG_HOME}/CROCO_IN\n" >> mynamelist.sh
        [[ ${models[$k]} == "atm" ]] && printf "export ATM_EXE_DIR=\n" >> mynamelist.sh
        [[ ${models[$k]} == "wav" ]] && printf "export WAV_EXE_DIR=\n" >> mynamelist.sh
        [[ ${models[$k]} == "toy" ]] && printf "export TOY_EXE_DIR=${MY_CONFIG_HOME}/TOY_IN\n" >> mynamelist.sh
        [[ ${models[$k]} == "xios" ]] && printf "export XIOS_EXE_DIR=\n" >> mynamelist.sh
    done

    printf "#-------------------------------------------------------------------------------\n" >> mynamelist.sh
    printf "# Model settings\n" >> mynamelist.sh
    printf "# ------------------------------------------------------------------------------\n" >> mynamelist.sh

    for k in ${models[@]} ; do
        [ -f namelist_${k}.sh ] && cat ./namelist_${k}.sh >> mynamelist.sh
    done

    if [[ ${models[@]} =~ "toy" ]] && [[ ${models[@]} =~ "cpl" ]] ; then
        sed -e "s/export namcouplename=.*/export namcouplename=namcouple.base.\${RUNtype}\${istoy}/g" \
        mynamelist.sh > mynamelist1.sh
        mv mynamelist1.sh mynamelist.sh
        chmod 755 mynamelist.sh
    fi

    sed -e "s|export CEXPER=.*|export CEXPER=${MY_CONFIG_NAME}_exp1|g" \
        mynamelist.sh > mynamelist1.sh

    mv mynamelist1.sh mynamelist.sh
    chmod 755 mynamelist.sh
    mv mynamelist.sh ../../.
fi

