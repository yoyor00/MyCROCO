#!/bin/bash
#
# Update S. Jullien : Oct 2021
# Update : Apr. 2020
# G. Cambon : Sept. 2016
#
#set -x
#set -e
#==========================================================================================
# BEGIN USER MODIFICATIONS

# Machine you are working on
# Known machines: Linux DATARMOR IRENE JEANZAY
# ---------------------------------------------
MACHINE="Linux"

# croco source directory
# ---------------------
CROCO_DIR=$(cd $(dirname "$0"); pwd)

# croco_tools directory 
# ---------------------
TOOLS_DIR=$(cd $(dirname "$0")/../croco_tools; pwd)

# Configuration name
# ------------------
MY_CONFIG_NAME=BENGUELA

# Home and Work configuration directories
# ---------------------------------------

MY_CONFIG_HOME=${CROCO_DIR}
MY_CONFIG_WORK=${CROCO_DIR}

# Options of your configuration
# Known options: 
LIST_OPTIONS=$(cat << EOF	
         		#%%%% CROCO built-in codes and toolboxes --
			#			
                        # -- IMPORTANT : type of CONFIGS architecture 
			#
			# oce-dev     : croco :  all-in architecture, for forced croco and/or dev.  
                        #                => "classic" architecture
      			# OR 
			#					
			# oce-prod    : croco : architecture for production and/or coupled run
                        #                 => creation of a CROCO_IN dir.                      
                        #
			# --- CROCO built-in scripts and toolboxes
			# prepro      : for getting scripts for CROCO preprocessing
                        # inter       : for running interannual runs
                        # forc        : for using forecast scripts
                        # test_cases  : for running test cases
 			# xios        : xios server xml file
			#			  
                        # --- CROCO built-in codes 
			# pisces      : pisces inputs
                        # agrif       : inputs for nests
                        # sediment    : inputs for sediment
                        # mustang     : mustang model
                       	# xios        : xios server xml file
			#			
                        #%%%%% External codes and toolbox => oce-prod needed
                        # cpl         : for coupling with OASIS                                    
                        # atm wav toy : other models for coupling (atm, wav, cpl, toy)
			#
			#			 			
			#%%%% All options : 			
                        # all-dev OR all-prod for CROCO built-in codes 
			#			
			# all                 for external codes (atm, wav, cpl, toy)
                        #                                           => oce-prod needed
			#%%%%    
EOF
	    )

models_incroco=( oce-dev test_cases agrif xios inter forc pisces sediment mustang oanalysis prepro )

#models_incroco=( oce-prod test_cases agrif xios inter forc pisces sediment mustang oanalysis prepro )
#models_external=( cpl atm wav toy )

# example for coupled model
# models_incroco=( oce-prod cpl wav atm toy )
# models_external=( cpl atm wav toy )

# example for oce-dev all
# models_incroco=( all-dev )

# example for oce-prod all + all other external
# models_incroco=( all-prod )
# models_external=( all )

# END USER MODIFICATIONS
#==========================================================================================

allmodels_incroco_dev=( oce-dev xios test_cases agrif inter forc pisces sediment mustang oanalysis prepro )
allmodels_incroco_prod=( oce-prod xios test_cases agrif inter forc pisces sediment mustang oanalysis prepro )
allmodels_external=( cpl wav atm toy )

x_f=0

while getopts :hfd:w:s:t:n:o: V
do
  case $V in
    ('h') cat << EOF
Script to setup your own configuration.
    What is does :
     - Copy the original cppdefs.h, param.h and *.in files needed
     - Copy the original crocotools_param.m and start.m file from croco_tools/
     - Copy the original run_croco*.bash file from croco/SCRIPTS/Plurimonths_scripts/
    Usage:
    Use the command line:
    ./create_config.bash -d MY_CONFIG_HOME -w MY_CONFIG_WORK -n MY_CONFIG_NAME -s CROCO_DIR -t TOOLS_DIR -o OPTS
    OR
    Edit the USER SECTION of the script to define the following variables :
     - CROCO_DIR       : location of  croco sources directory
     - MY_CONFIG_HOME  : location of the repository to store the configuration
     - MY_CONFIG_WORK  : location of the repository to store the configuration large input files, and where it will be run
     - MY_CONFIG_NAME  : name of the configuration
     - TOOLS_DIR       : location of  croco_tools directory
     - OPTS            : options of your configuration, comma separated (-o OPT1,OP2,OPT3 ...), with keywords in :
$LIST_OPTIONS  

EOF
    exit 0;;
    ('f')  x_f=1;;
    ('d')  x_d=${OPTARG};;
    ('w')  x_d=${OPTARG};;
    ('s')  x_s=${OPTARG};;
    ('t')  x_t=${OPTARG};;
    ('n')  x_n=${OPTARG};;
    ('o')  x_o=${OPTARG/,/ };;
  esac
done
#shift $(($OPTIND-1));

CROCO_DIR="${x_s-$CROCO_DIR}"
TOOLS_DIR="${x_t-$TOOLS_DIR}"
MY_CONFIG_NAME=${x_n-$MY_CONFIG_NAME}
MY_CONFIG_HOME=${x_d-$MY_CONFIG_HOME}/${MY_CONFIG_NAME}
MY_CONFIG_WORK=${x_d-$MY_CONFIG_WORK}/${MY_CONFIG_NAME}
models_incroco=( ${x_o[@]-${models_incroco[@]}} )
models_external=( ${x_o[@]-${models_external[@]}} )

if [ "$models_incroco" == "all-dev" ]; then
    models_incroco=${allmodels_incroco_dev[@]}
elif [ "$models_incroco" == "all-prod" ]; then
    models_incroco=${allmodels_incroco_prod[@]}
fi
if [ "$models_external" == "all" ]; then
    models_external=${allmodels_external[@]}
fi

# some check
if [[ ${models_incroco[@]} =~ "oce-dev" ]] ; then
    echo "oce-dev is defined. all-in architecture and no external codes considered"
    models_external=( )
elif [[ ${models_incroco[@]} =~ "oce-prod" ]] ; then
    echo "oce-prod is defined. architecture for production and/or coupled run"
    echo "External codes may be considered"
fi

echo ""
echo "Your choices :"
echo " - CROCO_DIR        : ${CROCO_DIR}"
echo " - TOOLS_DIR        : ${TOOLS_DIR}"
echo " - CONFIG_HOME_DIR  : ${MY_CONFIG_HOME%$MY_CONFIG_NAME}"
echo " - CONFIG_WORK_DIR  : ${MY_CONFIG_WORK%$MY_CONFIG_NAME}"
echo " - CONFIG_NAME      : ${MY_CONFIG_NAME}"
echo " - OPTIONS_INCROCO  : ${models_incroco[@]}"
echo " - OPTIONS_EXTERNAL : ${models_external[@]}"

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

# Check if source are there
if [ ! -d ${CROCO_DIR} ]; then 
	echo 'Directory for croco not found ...'
	echo 'Check the CROCO_DIR variable ...'
	echo 'Exiting ...'
   exit 1
fi

# Check if tools are there
copy_tools=1
if [[ ! -d $TOOLS_DIR  &&  $x_f -eq 0 ]]; then 
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

if [[ ${models_incroco[@]} =~ "oce-dev" ]] || [[ ${models_incroco[@]} =~ "oce-prod" ]] ; then
    echo 'Copy CROCO useful scripts and input files'
    echo '-----------------------------------------'
    # CROCO general
    if [[ ${models_incroco[@]} =~ "oce-prod" ]] ; then
	# Create directories
	mkdir -p $MY_CONFIG_HOME/CROCO_IN
	mkdir -p $MY_CONFIG_WORK/CROCO_FILES
	MY_CROCO_DIR=$MY_CONFIG_HOME/CROCO_IN/
	MY_XIOS_DIR=$MY_CONFIG_HOME/XIOS_IN/
	
    elif [[ ${models_incroco[@]} =~ "oce-dev" ]] ; then
	# Create directories
	mkdir -p $MY_CONFIG_HOME
	mkdir -p $MY_CONFIG_WORK/CROCO_FILES
	MY_CROCO_DIR=$MY_CONFIG_HOME/
	MY_XIOS_DIR=$MY_CONFIG_HOME/
    fi
    cp -f ${CROCO_DIR}/OCEAN/cppdefs.h $MY_CROCO_DIR.
    cp -f ${CROCO_DIR}/OCEAN/cppdefs_dev.h $MY_CROCO_DIR.
    cp -f ${CROCO_DIR}/OCEAN/param.h $MY_CROCO_DIR.
    
     PAT=$(grep ^SOURCE ${CROCO_DIR}/OCEAN/jobcomp)
     sed -e "s!${PAT}!SOURCE=${CROCO_DIR}/OCEAN!g" $CROCO_DIR/OCEAN/jobcomp > $MY_CROCO_DIR/jobcomp
     chmod +x $MY_CROCO_DIR/jobcomp

    cp -f ${CROCO_DIR}/OCEAN/croco.in $MY_CROCO_DIR.
    cp -f ${CROCO_DIR}/OCEAN/croco_stations.in $MY_CROCO_DIR.
    # TEST_CASES
    if [[ ${models_incroco[@]} =~ "test_cases" ]] ; then
      cp -Rf ${CROCO_DIR}/TEST_CASES $MY_CROCO_DIR.
    fi
    # AGRIF
    if [[ ${models_incroco[@]} =~ "agrif" ]] ; then
      cp -f ${CROCO_DIR}/OCEAN/croco.in.1 $MY_CROCO_DIR.
    fi
    # INTER
    if [[ ${models_incroco[@]} =~ "inter" ]] ; then
      cp -f ${CROCO_DIR}/OCEAN/croco_inter.in $MY_CROCO_DIR.
    fi
    # FORECAST
    if [[ ${models_incroco[@]} =~ "forc" ]] ; then
      cp -f ${CROCO_DIR}/OCEAN/croco_forecast.in $MY_CROCO_DIR.
      cp -f ${CROCO_DIR}/OCEAN/croco_hindcast.in $MY_CROCO_DIR.
    fi
    # PISCES
    if [[ ${models_incroco[@]} =~ "pisces" ]] ; then
      cp -f ${CROCO_DIR}/PISCES/*namelist* $MY_CROCO_DIR.
    fi
    # SEDIMENT
    if [[ ${models_incroco[@]} =~ "sediment" ]] ; then
      cp -f ${CROCO_DIR}/OCEAN/sediment.in $MY_CROCO_DIR.
    fi
    # MUSTANG
    if [[ ${models_incroco[@]} =~ "mustang" ]] ; then
      mkdir -p $MY_CROCO_DIR/MUSTANG_NAMELIST
      cp -f ${CROCO_DIR}/MUSTANG/NAM_CASES/*txt $MY_CROCO_DIR/MUSTANG_NAMELIST/.
    fi
    # OANALYSIS
    if [[ ${models_incroco[@]} =~ "oanalysis" ]] ; then
       cp -Rf ${CROCO_DIR}/SCRIPTS/NAMELIST_OANALYSIS $MY_CROCO_DIR.
    fi
   # XIOS
    if [[ ${models_incroco[@]} =~ "xios" ]] ; then
     mkdir -p $MY_XIOS_DIR
     cp -Rf ${CROCO_DIR}/XIOS/*.xml* $MY_XIOS_DIR.
     cp -Rf ${CROCO_DIR}/XIOS/xios_launch.file $MY_XIOS_DIR.
     cp -Rf ${CROCO_DIR}/XIOS/README_XIOS $MY_XIOS_DIR.
    fi
    # PREPROCESSING
    if [[ ${models_incroco[@]} =~ "prepro" ]] ; then
       cp -Rf $TOOLS_DIR/start.m $MY_CROCO_DIR.
       cp -Rf $TOOLS_DIR/oct_start.m $MY_CROCO_DIR.
       cp -Rf $TOOLS_DIR/crocotools_param.m $MY_CROCO_DIR.
       cp -Rf $TOOLS_DIR/Town/town.dat $MY_CROCO_DIR.
    fi
    # SCRIPTS FOR RUNNING
    if [[ ${models_incroco[@]} =~ "inter" ]] ; then
       cp -Rf ${CROCO_DIR}/SCRIPTS/Plurimonths_scripts/*.bash $MY_CONFIG_HOME/
    fi
fi

### Coupling and other models to be coupled with ###

# OASIS
if [[ ${models_external[@]} =~ "cpl" ]] ; then
    echo 'Copy OASIS useful scripts and input files'
    echo '-----------------------------------------'
    mkdir -p $MY_CONFIG_HOME/OASIS_IN
    cp -r ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/OASIS_IN/* $MY_CONFIG_HOME/OASIS_IN/.
    if [[ ${models_incroco[@]} =~ "oce-prod" ]] ; then
      cp -r ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/CROCO_IN/* $MY_CONFIG_HOME/CROCO_IN/.
    fi
    if [[ ${models_incroco[@]} =~ "prepro" ]] ; then
      mkdir -p $MY_CONFIG_HOME/PREPRO
      cp -r $TOOLS_DIR/Coupling_tools/* $MY_CONFIG_HOME/PREPRO/.
    fi
fi

# WW3
if [[ ${models_external[@]} =~ "wav" ]] ; then
    echo 'Copy WW3 useful scripts and input files'
    echo '-----------------------------------------'
    mkdir -p $MY_CONFIG_HOME/WW3_IN
    mkdir -p $MY_CONFIG_WORK/WW3_FILES
    cp -r ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/WW3_IN/* $MY_CONFIG_HOME/WW3_IN/.
fi

# WRF
if [[ ${models_external[@]} =~ "atm" ]] ; then
    echo 'Copy WRF useful scripts and input files'
    echo '-----------------------------------------'
    mkdir -p $MY_CONFIG_HOME/WRF_IN
    mkdir -p $MY_CONFIG_WORK/WRF_FILES
    cp -r ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/WRF_IN/* $MY_CONFIG_HOME/WRF_IN/.
fi

# TOY
if [[ ${models_external[@]} =~ "toy" ]] ; then
    echo 'Copy TOY sources, useful scripts and input files'
    echo '------------------------------------------------'
    mkdir -p $MY_CONFIG_HOME/TOY_IN
    mkdir -p $MY_CONFIG_WORK/TOY_FILES
    cp -r ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/TOY_IN/* $MY_CONFIG_HOME/TOY_IN/.
fi

# Coupling scripts
if [[ ${models_external[@]} =~ "cpl" ]] || [[ ${models_external[@]} =~ "wav" ]] || [[ ${models_external[@]} =~ "atm" ]] || [[ ${models_external[@]} =~ "toy" ]] ; then
    echo 'Copy scripts for coupled runs'
    echo '-----------------------------'
    [ -d $MY_CONFIG_HOME/ROUTINES ] && \rm -Rf $MY_CONFIG_HOME/ROUTINES
    cp -Rf ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/*.sh $MY_CONFIG_HOME/
    cp -Rf ${CROCO_DIR}/SCRIPTS/SCRIPTS_COUPLING/SCRIPTS_TOOLBOX/ $MY_CONFIG_HOME/ROUTINES

    # Edit myjob.sh to add CPU lines for each model
    cd $MY_CONFIG_HOME/
    [ -f myjob.tmp ] && rm -Rf myjob.tmp
    [[ ${models_incroco[@]} =~ "oce-prod" ]] && printf "export NP_OCEX=2 \nexport NP_OCEY=2\n" >> myjob.tmp
    [[ ${models_external[@]} =~ "wav" ]] && printf "export NP_WAV=14 \n" >> myjob.tmp
    [[ ${models_external[@]} =~ "atm" ]] && printf "export NP_ATM=12 \n" >> myjob.tmp
    [[ ${models_external[@]} =~ "toy" ]] && printf "export NP_TOY=2 \n" >> myjob.tmp
    [[ ${models_incroco[@]} =~ "xios" ]] && printf "export NP_XIOS_ATM=1\nexport NP_XIOS_OCE=1\n" >> myjob.tmp

    if [[ ${models_external[@]} =~ "atm" ]] ; then
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
    [[ ${models_external[@]} =~ "cpl" ]] && printf "export CPL=\"\"\n" >> tmppath
    [[ ${models_incroco[@]} =~ "oce-prod" ]] && printf "export OCE=\"${CROCO_DIR}/OCEAN\"\n" >> tmppath
    [[ ${models_external[@]} =~ "atm" ]] && printf "export ATM=\"\"\n" >> tmppath
    [[ ${models_external[@]} =~ "wav" ]] && printf "export WAV=\"\"\n" >> tmppath
    [[ ${models_external[@]} =~ "toy" ]] && printf "export TOY=\"\${CHOME}/TOY_IN\"\n" >> tmppath
    [[ ${models_incroco[@]} =~ "xios" ]] && printf "export XIOS=\"\"\n" >> tmppath

    [[ ${models_external[@]} =~ "cpl" ]] && cat ./path_cpl.sh >> tmppath
    [[ ${models_incroco[@]} =~ "oce-prod" ]] && cat ./path_oce.sh >> tmppath
    [[ ${models_external[@]} =~ "atm" ]] && cat ./path_atm.sh >> tmppath
    [[ ${models_external[@]} =~ "wav" ]] && cat ./path_wav.sh >> tmppath
    [[ ${models_external[@]} =~ "toy" ]] && cat ./path_toy.sh >> tmppath
    [[ ${models_incroco[@]} =~ "xios" ]] && cat ./path_xios.sh>> tmppath

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

    [[ ${models_external[@]} =~ "atm" ]] && cat ./myenv.${MACHINE}.wrf >> tmpenv 
    [[ ${models_external[@]} =~ "wav" ]] && cat ./myenv.${MACHINE}.ww3 >> tmpenv

    mv tmpenv ${MY_CONFIG_HOME}/
    cd ${MY_CONFIG_HOME}

    # concatenate path and env file
    cat tmpenv tmppath > myenv_mypath.sh
    chmod 755 myenv_mypath.sh
    rm -rf tmppath tmpenv

    # Create the namelist file
    cd ${MY_CONFIG_HOME}/ROUTINES/NAMELISTS
    cp namelist_head.sh mynamelist.sh

    if [[ ${models_external[@]} =~ "cpl" ]]; then
        if [[ ${models_incroco[@]} =~ "oce-prod" ]] && [[ ${models_external[@]} =~ "wav" ]] && [[ ${models_external[@]} =~ "atm" ]] ; then
            printf "export RUNtype=owa\n#\n" >> mynamelist.sh
        elif [[ ${models_incroco[@]} =~ "oce-prod" ]] && [[ ${models_external[@]} =~ "wav" ]] ; then
            printf "export RUNtype=ow\n#\n" >> mynamelist.sh
        elif [[ ${models_incroco[@]} =~ "oce-prod" ]] && [[ ${models_external[@]} =~ "atm" ]]; then
            printf "export RUNtype=oa\n#\n" >> mynamelist.sh
        elif [[ ${models_external[@]} =~ "wav" ]] && [[ ${models_external[@]} =~ "atm" ]]; then
            printf "export RUNtype=aw\n#\n" >> mynamelist.sh
        elif [[ ${models_external[@]} =~ "toy" ]]; then
            printf "export RUNtype=Put the type here (ow/oa/aw/owa)\n#\n" >> mynamelist.sh
        else 
            printf "export RUNtype=frc\n#\n" >> mynamelist.sh
        fi
    else
        printf "export RUNtype=frc\n#\n" >> mynamelist.sh
    fi

    if [[ ${models_external[@]} =~ "atm" ]]; then
        printf "export USE_ATM=1\n" >> mynamelist.sh
        [[ ${models_incroco[@]} =~ "xios" ]] && printf "export USE_XIOS_ATM=0\n" >> mynamelist.sh
    fi
    if [[ ${models_incroco[@]} =~ "oce-prod" ]]; then
        printf "export USE_OCE=1\n" >> mynamelist.sh
        [[ ${models_incroco[@]} =~ "xios" ]] && printf "export USE_XIOS_OCE=0\n" >> mynamelist.sh
    fi
    if [[ ${models_external[@]} =~ "wav" ]]; then
        printf "export USE_WAV=1\n" >> mynamelist.sh
    fi
    if [[ ${models_external[@]} =~ "toy" ]]; then
        cat ./namelist_head_toy.sh >> mynamelist.sh
    fi

    cat ./namelist_rundir.sh >> mynamelist.sh

    [[ ${models_incroco[@]} =~ "oce-prod" ]] && printf "export OCE_EXE_DIR=${MY_CONFIG_HOME}/CROCO_IN\n" >> mynamelist.sh
    [[ ${models_external[@]} =~ "atm" ]] && printf "export ATM_EXE_DIR=\n" >> mynamelist.sh
    [[ ${models_external[@]} =~ "wav" ]] && printf "export WAV_EXE_DIR=\n" >> mynamelist.sh
    [[ ${models_external[@]} =~ "toy" ]] && printf "export TOY_EXE_DIR=${MY_CONFIG_HOME}/TOY_IN\n" >> mynamelist.sh
    [[ ${models_incroco[@]} =~ "xios" ]] && printf "export XIOS_EXE_DIR=\n" >> mynamelist.sh

    printf "#-------------------------------------------------------------------------------\n" >> mynamelist.sh
    printf "# Model settings\n" >> mynamelist.sh
    printf "# ------------------------------------------------------------------------------\n" >> mynamelist.sh


    [[ ${models_incroco[@]} =~ "oce-prod" ]] && cat ./namelist_oce.sh >> mynamelist.sh
    for k in ${models_external[@]} ; do
        [ -f namelist_${k}.sh ] && cat ./namelist_${k}.sh >> mynamelist.sh
    done
    [[ ${models_incroco[@]} =~ "xios" ]] && cat ./namelist_xios.sh >> mynamelist.sh

    if [[ ${models_external[@]} =~ "toy" ]] && [[ ${models_external[@]} =~ "cpl" ]] ; then
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

