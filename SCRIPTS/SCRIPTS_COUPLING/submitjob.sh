#!/bin/bash  

source ./myenv_mypath.sh
set -u
umask 022

#-------------------------------------------------------------------------------
#  namelist of the experiment
#-------------------------------------------------------------------------------
#cat mypath.sh >> mynamelist.tmp
cat mynamelist.sh > mynamelist.tmp
cat ./SCRIPTS_TOOLBOX/NAMELISTS/namelist_tail.sh >> mynamelist.tmp
cat myjob.sh >> mynamelist.tmp
cat ./SCRIPTS_TOOLBOX/common_definitions.sh >> mynamelist.tmp

. ./mynamelist.tmp

[ ! -d ${JOBDIR_ROOT} ] && mkdir -p ${JOBDIR_ROOT}  # for the first submitjob.sh call

echo "  "
echo "jobs and logs directory is $JOBDIR_ROOT"
echo "  "

# copy the base scripts into the jobdir
cpfile myenv_mypath.sh ${JOBDIR_ROOT}
cpfile mynamelist.sh ${JOBDIR_ROOT}
cpfile myjob.sh ${JOBDIR_ROOT}

cd ${JOBDIR_ROOT} 
ls ${jobname}  > /dev/null  2>&1 
if [ "$?" -eq "0" ] ; then
   if [ ${CHAINED_JOB} == "FALSE" ]; then 
       echo "  "
       echo "!!!!!!!! WARNING !!!!!!!!!"
       echo "  "
       echo "A ${jobname} job already exists in  ${JOBDIR_ROOT}"
       echo -n "  Do you want to remove it and launch the new job? [y/n]"
       read answer
       if [  "x$answer" = "xy" ]; then
          echo " " 
          echo "Creating and launching new job"
          echo "   "
       else
          echo "  " 
          echo "Exiting..."
          echo "   "
          exit
       fi
       unset -v answer
   elif [ ${CHAINED_JOB} == "TRUE" ] && [ ${DATE_BEGIN_JOB} -eq ${DATE_BEGIN_EXP} ]; then
       echo "  "
       echo "!!!!!!!! WARNING !!!!!!!!!"
       echo "  "
       echo "A ${jobname} job already exists in  ${JOBDIR_ROOT}"
       echo -n "  Do you want to remove it and launch the new job? [y/n]"
       read answer
       if [  "x$answer" = "xy" ]; then
          echo " " 
          echo "Creating and launching new job"
          echo "   "
       else
          echo "  " 
          echo "Exiting..."
          echo "   "
          exit
       fi
       unset -v answer
   fi
fi
cd -

#-------------------------------------------------------------------------------
# Checking experiment options
#-------------------------------------------------------------------------------
echo "  "
echo " Checking experiment options... "

if [[ $RUNtype =~ .*a.* ]] ; then
  if [ $USE_ATM == 0 ] && [ $USE_TOYATM == 0 ]; then
    echo "  "
    echo "ERROR. RUNtype=$RUNtype and no atmosphere use defined (USE_ATM or USE_TOYATM)"
    echo " Exit"
    exit
  fi
fi
if [[ $RUNtype =~ .*o.* ]] ; then 
  if [ $USE_OCE == 0 ] && [ $USE_TOYOCE == 0 ]; then
    echo "  "
    echo "ERROR. RUNtype=$RUNtype and no ocean use defined (USE_OCE or USE_TOYOCE)"
    echo " Exit"
    exit
  fi
fi
if [[ $RUNtype =~ .*w.* ]] ; then 
  if [ $USE_WAV == 0 ] && [ $USE_TOYWAV == 0 ]; then
    echo "  "
    echo "ERROR. RUNtype=$RUNtype and no wave use defined (USE_WAV or USE_TOYWAV)"
    echo " Exit"
    exit
  fi
fi

if [ ${USE_CPL} -ge 1 ]; then
  if [ $(( ${CPL_FREQ} % ${DT_ATM} )) -ne 0 ] || \
     [ $(( ${CPL_FREQ} % ${DT_OCE} )) -ne 0 ] || \
     [ $(( ${CPL_FREQ} % ${DT_WAV} )) -ne 0 ] ; then
     printf "\n\n Problem of consistency between Coupling Frequency and Time Step with ATM, OCE or WAV model, we stop...\n\n" && exit 1
  fi
  if [ ${USE_TOY} -eq 1 ]; then 
      for k in `seq 0 $(( ${nbtoy} - 1))` ; do
          if [ $(( ${CPL_FREQ} % ${DT_TOY[$k]} )) -ne 0 ] ; then
              printf "\n\n Problem of consistency between Coupling Frequency and Time Step for TOY model, we stop...\n\n" && exit 1
          fi
      done
  fi
fi

#-------------------------------------------------------------------------------
# calendar computations (to check dates consistency)
#-------------------------------------------------------------------------------

. ${SCRIPTDIR}/caltools.sh

#-------------------------------------------------------------------------------
# create job and submit it
#-------------------------------------------------------------------------------

if [ ${USE_OCE}  -eq 1 ]; then
    if [[ ${MPI_NOLAND} == "TRUE" ]]; then
        TOTOCE=${NP_OCE}
    else
        TOTOCE=$(( $NP_OCEX * $NP_OCEY )) 
   fi
else
    TOTOCE=0
fi
[ ${USE_ATM}  -eq 1 ] && TOTATM=$NP_ATM  || TOTATM=0
[ ${USE_WAV}  -eq 1 ] && TOTWAV=$NP_WAV  || TOTWAV=0
[ ${USE_TOY}  -ge 1 ] && { TOTTOY=0 ; for k in `seq 0 $(( ${nbtoy} - 1))`; do TOTTOY=$(( $TOTTOY + $NP_TOY)) ; done;}  || TOTTOY=0
[ ${USE_XIOS_ATM} -eq 1 ] && TOTXIO=$NP_XIOS_ATM || TOTXIO=0
[ ${USE_XIOS_OCE} -eq 1 ] && TOTXIO=$(( ${TOTXIO} + ${NP_XIOS_OCE} ))
totalcore=$(( $TOTOCE + $TOTATM + $TOTWAV + $TOTTOY + $TOTXIO ))

if [[ ${MACHINE} == "DATARMOR" ]]; then
    nbnode=$(( $totalcore /29 +1))
    [[ ${totalcore} -ge 28 ]] && totalcore=28
else
     nbnode=0
fi


if [ ${MACHINE} == "IRENE" ]; then
    timedur=${TIMEJOB}
else
    timedur=$( sec2hour ${TIMEJOB} )
fi

sed -e "/< insert here variables definitions >/r mynamelist.tmp" \
    -e "s/<exp>/${ROOT_NAME_1}/g" \
    -e "s/<nbnode>/${nbnode}/g" \
    -e "s/<nmpi>/${totalcore}/g" \
    -e "s/<projectid>/${projectid}/g" \
    -e "s/<timedur>/${timedur}/g" \
    ./SCRIPTS_TOOLBOX/MACHINE/${MACHINE}/header.${MACHINE} > HEADER_tmp
    cat HEADER_tmp ./SCRIPTS_TOOLBOX/job.base.sh >  ${JOBDIR_ROOT}/${jobname}
    \rm HEADER_tmp
    \rm ./mynamelist.tmp


cd ${JOBDIR_ROOT}
chmod 755 ${jobname}

printf "\n  HOSTNAME: `hostname`\n     =>    COMPUTER: ${MACHINE}\n"  
echo
printf "  CONFIG: ${CONFIG}\n"  
printf "  CEXPER: ${CEXPER}\n"  
echo
printf "  jobname: ${jobname}\n"  
echo
#printf "  ROOT_NAME_1: ${ROOT_NAME_1}\n"  
#printf "  ROOT_NAME_2: ${ROOT_NAME_2}\n"  
#printf "  ROOT_NAME_3: ${ROOT_NAME_3}\n"  
printf "  EXEDIR: ${EXEDIR_ROOT}\n"  
printf "  OUTPUTDIR: ${OUTPUTDIR_ROOT}\n"  
printf "  RESTDIR_OUT: ${RESTDIR_ROOT}\n"  
printf "  JOBDIR: ${JOBDIR_ROOT}\n"  

if [ "${SCRIPT_DEBUG}" == "TRUE" ] ; then
   printf "\n\n\n\n  SCRIPT_DEBUG=${SCRIPT_DEBUG}  Script debug mode => No submission in the queue\n\n\n\n"
else 
    if [ ${CHAINED_JOB} == "TRUE" ]; then
#        [[ ${RESTART_FLAG} == "FALSE" ]] && . ${SCRIPTDIR}/chained_job.sh
        . ${SCRIPTDIR}/chained_job.sh
    else
       ${QSUB}${jobname}
    fi 
#
    if [ "${MODE_TEST}" != "" ] ; then
        printf "\n\n\n\n  MODE_TEST=${MODE_TEST}  Test mode and non production => No job chaining.\n\n\n\n"
    fi
fi

