#!/bin/bash

if [ ${MACHINE} == "DATARMOR" ] ; then
    ${QSUB} -h -N ${ROOT_NAME_1} ${jobname}
    sub_ext=".pbs"
    launchcmd="-W depend=afterany:$( echo $( qselect -N ${ROOT_NAME_1}) | cut -d "." -f 1 )"
elif [ ${MACHINE} == "IRENE" ]; then
    ${QSUB} ${jobname}
    sub_ext=".sh"
    prejobname=${ROOT_NAME_1}
    echo "$( ccc_mstat -f -u $USER )" >tmp.text
    cnt=$(( $(wc -l tmp.text | awk '{ print $1 }') + 1 )) # check if there are already jobs running
    \rm -f tmp.text  
elif [ ${MACHINE} == "JEANZAY" ] || [ ${MACHINE} == "NEA" ]; then
    ${QSUB} -H ${jobname}
    sub_ext=".sh"
    prejobname=${ROOT_NAME_1}
    echo "$( squeue --format="%.18i %.100j" -u $USER )" >tmp.text
    cnt=$(( $(wc -l tmp.text | awk '{ print $1 }') + 1 ))
    for linenb in `seq 1 $cnt`; do 
        line=$(sed -n "${linenb}p" tmp.text  )
	#[[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -c 1-7 ); break; }      
        # previous line supposes that jobid has 7 characters. If not, the following line could replace it
        [[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -d "." -f 1 | cut -d " " -f 1 ); break; } 
    done
    firstjobid=${prejobid}
else
    printf "\n\n Chained job for your machine is not set up yet, we stop...\n\n" && exit 1

fi


cd "${SCRIPTDIR}/"

# Start from the end of the current job (exclusive).
months=${MONTH_BEGIN_JOBp1}
days=${DAY_BEGIN_JOBp1}
years=${YEAR_BEGIN_JOBp1}
hours=${HOUR_BEGIN_JOBp1}

# NBJOB includes the current job: only NBJOB-1 jobs remain to submit.
for ((jobleft = NBJOB - 1; jobleft > 0; jobleft--)); do
    newsdate=$(makedate "$months" "$days" "$years" "$hours")
    read -r monthe daye yeare houre <<< "$(next_job_date "$months" "$days" "$years" "$hours")"
    newedate=$(makedate "$monthe" "$daye" "$yeare" "$houre")

    # Preserve historic job names for midnight/day-or-month runs.
    # Include hours for hourly runs or non-midnight boundaries to avoid collisions.
    if (( JOB_DUR_HOUR > 0 || hours != 0 || houre != 0 )); then
        printf -v newsdate '%s%02d' "$newsdate" "$hours"
        printf -v newedate '%s%02d' "$newedate" "$houre"
    fi

    future_date="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
    future_job="job_${future_date}${sub_ext}"
    cd "${JOBDIR_ROOT}" || exit 1

    # Preserve the original sed replacements, adding the starting hour.
    sed -e "s/YEAR_BEGIN_JOB=${YEAR_BEGIN_JOB}/YEAR_BEGIN_JOB=${years}/" \
        -e "s/MONTH_BEGIN_JOB=${MONTH_BEGIN_JOB}/MONTH_BEGIN_JOB=${months}/" \
        -e "s/DAY_BEGIN_JOB=${DAY_BEGIN_JOB}/DAY_BEGIN_JOB=${days}/" \
        -e "s/HOUR_BEGIN_JOB=${HOUR_BEGIN_JOB}/HOUR_BEGIN_JOB=${hours}/" \
        -e "s/export NBJOB=.*/export NBJOB=${jobleft}/" \
        -e "s/export RESTART_FLAG=.*/export RESTART_FLAG=\"TRUE\"/" \
        -e "s/export CHAINED_JOB=.*/export CHAINED_JOB=\"FALSE\"/" \
        ${jobname} > ${future_job}

    chmod 755 ${future_job}
    #
    if [ ${MACHINE} == "DATARMOR" ] ; then 
        newjobname="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
        ${QSUB} -N ${newjobname} ${launchcmd} ${future_job} 
        launchcmd="-W depend=afterany:$( echo $( qselect -N ${newjobname}) | cut -d "." -f 1 )"
        cd ${SCRIPTDIR}/
#
    elif [ ${MACHINE} == "IRENE" ] ; then
        # 
        newjobname="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
        sed -e "s/#MSUB -r .*/#MSUB -r ${newjobname}/" \
            -e "s/#MSUB -o .*/#MSUB -o ${newjobname}.jobid_%I.o/" \
            -e "s/#MSUB -e .*/#MSUB -e ${newjobname}.jobid_%I.e/" \
        ${future_job} > ${future_job}.tmp
        mv ${future_job}.tmp ${future_job}
        chmod 755 ${future_job}
        #
        echo "$( ccc_mstat -f -u $USER )" >tmp.text
        for linenb in `seq 3 $cnt`; do
            line=$(sed -n "${linenb}p" tmp.text  )
            [[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -c 1-7 ); break; }
        done
        \rm -f tmp.text
        #
        cnt=$(( ${cnt} + 1 ))
        ${QSUB} -a ${prejobid} ${future_job} 
        prejobname="${newjobname}"
        cd ${SCRIPTDIR}/
#
    elif [ ${MACHINE} == "JEANZAY" ]  || [ ${MACHINE} == "NEA" ] ; then
	newjobname="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
	sed -e "s/#SBATCH --job-name=.*/#SBATCH --job-name=${newjobname}/" \
	    -e "s/#SBATCH --output=.*/#SBATCH --output=${newjobname}.out/" \
	    -e "s/#SBATCH --error=.*/#SBATCH --error=${newjobname}.out/" \
	${future_job} > ${future_job}.tmp
	mv ${future_job}.tmp ${future_job}
	chmod 755 ${future_job}
	echo "$( squeue --format="%.18i %.100j" -u $USER )" >tmp.text
        for linenb in `seq 1 $cnt`; do
            line=$(sed -n "${linenb}p" tmp.text  )
            #[[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -c 1-7 ); break; }
            # previous line supposes that jobid has 7 characters. If not, the following line could replace it
            [[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -d "." -f 1 | cut -d " " -f 1 ); break; } 
        done
        \rm -f tmp.text
        #
        cnt=$(( ${cnt} + 1 ))
	launchcmd="--dependency=afterok:${prejobid}"
        ${QSUB} ${launchcmd} ${future_job}
	prejobname="${newjobname}"
        cd ${SCRIPTDIR}

    fi
#    
    # Next job begins exactly when the previous job ends.
    months=$monthe
    days=$daye
    years=$yeare
    hours=$houre

done
#
cd ${JOBDIR_ROOT}

if [ "${MACHINE}" == "DATARMOR" ]; then
    jobid=$( echo $( qselect -N ${ROOT_NAME_1}) | cut -d "." -f 1 )
    qrls ${jobid}
elif [ "${MACHINE}" == "JEANZAY" ] || [ "${MACHINE}" == "NEA" ]; then
    scontrol release ${firstjobid}
else
    printf "\n\n Chained job for your machine ${MACHINE} is not set up yet, we stop...\n\n" && exit 1
fi

