#!/bin/bash

#
if [ ${interponline} -eq 1 ]; then 

    if [ ${frc_ext} == "ECMWF" ]; then
        vnames='T2M SSTK U10M V10M Q STR STRD SSR TP EWSS NSSS'
    elif [ ${frc_ext} == "AROME" ]; then
        vnames='AROME'
    else
        vnames='Temperature_height_above_ground Specific_humidity Precipitation_rate Downward_Short-Wave_Rad_Flux_surface Upward_Short-Wave_Rad_Flux_surface Downward_Long-Wave_Rad_Flux Upward_Long-Wave_Rad_Flux_surface U-component_of_wind V-component_of_wind'
    fi
#    
    printf "Creating link to data for the job duration\n"
#          
    echo "Checking if Previous month is needed"
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
    if [[ ${RESTART_FLAG} == "FALSE" ]]; then
        filefrom="${OCE_FILES_DIR}/croco_${ini_ext}_Y${cur_Y}M${cur_M}.nc"
    else
        filefrom="${RESTDIR_IN}/croco_rst_${DATE_END_JOBm1}.nc"
    fi
    # scrum_time of ini file
    tstartinsec=$( echo $( ncdump -v scrum_time ${filefrom} | grep 'scrum_time =' | cut -d '=' -f 2| cut -d ' ' -f 2 ))
    tstartinsec=$(( ${tstartinsec} + 0.5 * ${DT_OCE})) # =0.5*dt like in croco 
    # Find first time value in forcing file
    fieldname=( echo "$vnames" | awk '{print $1}' )
    ncdump -v time "${OCE_FILES_ONLINEDIR}/$fieldname_Y${cur_Y}M${cur_M}.nc" | grep -n 'time =' > tmp$$
    tstartfrc=$(( $( sed -n -e "3 p" text | cut -d '=' -f 2 | cut -d ',' -f 1 ) * 86400 ))
    rm -rf tmp$$
    [[${tstartinsec} -le ${tstartfrc} ]] && { echo "Previous month is needed!"; loopstrt=-1 ;} || { loopstrt=0 ;}      
#
    for i in `seq ${loopstrt} $(( ${JOB_DUR_MTH} ))`; do
        [ ${i} -eq -1 ] && printf "Adding link to the previous month (for temporal interpolation)\n"
        [ ${i} -eq ${JOB_DUR_MTH} ] && printf "Adding link to the following month (for temporal interpolation)\n"

        mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
        cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
        cur_M=$( printf "%01d\n"  $( echo $mdy | cut -d " " -f 1) )
  
        for varname in ${vnames} ; do
            [[ ! -f "${OCE_FILES_ONLINEDIR}/${varname}_Y${cur_Y}M${cur_M}.nc" ]] && { echo "File ${varname}_Y${cur_Y}M${cur_M}.nc is missing for online interpolation, we stop..." ; exit ;}
            ${io_getfile} ${OCE_FILES_ONLINEDIR}/${varname}_Y${cur_Y}M${cur_M}.nc ./
        done
    done
    if [ ${JOB_DUR_MTH} -eq 0 ] ; then
        mdy=$( valid_date $(( $MONTH_BEGIN_JOB + 1 )) $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
        cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
        cur_M=$( printf "%01d\n"  $( echo $mdy | cut -d " " -f 1) )
        for varname in ${vnames} ; do
            ${io_getfile} ${OCE_FILES_ONLINEDIR}/${varname}_Y${cur_Y}M${cur_M}.nc ./
        done
    fi
#
else
    for nn in $( seq 0 ${AGRIFZ} ); do
        if [ ${nn} -gt 0 ]; then
            agrif_ext=".${nn}"
        else
            agrif_ext=""
        fi
        for i in `seq 0 $(( ${JOB_DUR_MTH}-1 ))`; do
            if [ ${JOB_DUR_MTH} -eq 1 ]; then
                cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
                cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
                extend=$( echo $frc_ext | cut -c 1-3 )
                ln -sf ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext} croco_${extend}.nc${agrif_ext}
            else
                printf "\n\nSurface forcing is not yet implemented for job duration longer than a month (need to solve overlapping problem)=>  we stop...\n\n"
                exit 1
            fi
        done

        if [ ${JOB_DUR_MTH} -eq 0 ] ; then
            printf "Job duration is less than a month ---> Using netcdf of the current month\n"
            cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
            cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
            extend=$( echo $frc_ext | cut -c 1-3 )
            ln -sf ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext} croco_${extend}.nc${agrif_ext}
        fi
    done
fi
