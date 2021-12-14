#!/bin/bash
#set -ue
#set -vx
#
##======================================================================
##======================================================================
## This script automatically modify the CROCO "namelist"
##======================================================================
##======================================================================
#
#
##
##======================================================================
##----------------------------------------------------------------------
##    I. Calendar computations
##----------------------------------------------------------------------
##======================================================================
##
#
# some calendar tools...
#
#                           *** WARNING ***
#   To get back the functions defined in caltools.sh, we have to call 
#   it with ". caltools.sh" instruction. If we directly call "caltools.sh",
#   we will not have the functions back
#
. ${SCRIPTDIR}/caltools.sh
#
##
##======================================================================
##----------------------------------------------------------------------
##    II. modify namelist
##----------------------------------------------------------------------
##======================================================================
##
#
for nn in $( seq 0 ${AGRIFZ} ) 
do
    if [ ${nn} -gt 0 ];    then
	namfile=croco.in.${nn}
        agrif_ext=".${nn}"
	cpfile ${OCE_NAM_DIR}/croco.in.base.${nn} ${namfile}
        cpfile ${OCE_NAM_DIR}/AGRIF_FixedGrids.in ./
	SUBTIME=$( sed -n -e "$(( 2 * ${nn} )) p" AGRIF_FixedGrids.in | awk '{print $7 }' )
    else
	namfile=croco.in
	agrif_ext=""
	cp ${OCE_NAM_DIR}/croco.in.base ${namfile}
	SUBTIME=1
    fi
    DT_OCE_2=$(( ${DT_OCE} / ${SUBTIME} ))
#-------
## Number of time step per day
#-------
    OCE_NDT_DAY=$(( 86400 / ${DT_OCE_2} ))
    OCE_NTIMES=$(( ( ${JDAY_END_JOB} - ${JDAY_BEGIN_JOB} + 1 ) * ${OCE_NDT_DAY}     ))
#
#-------
# change some namelist values
#-------
# Change in endding date for online interpolation

    mdy=$( valid_date $(( $MONTH_END_JOB + 1 )) $DAY_END_JOB $YEAR_END_JOB )
    end_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
    end_M=$( printf "%01d\n"  $( echo $mdy | cut -d " " -f 1) )
#

    printf "Computing the origin_date from start_date and scrum_time\n"
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 ) 
    cur_D=$( echo $DATE_BEGIN_JOB | cut -c 7-8 )
    scrumt=$( ncdump -v scrum_time croco_ini.nc${agrif_ext}| grep "scrum_time = " | cut -d '=' -f 2 | cut -d ' ' -f 2)

    scrumtindays=$(( $scrumt/86400))

    mdy=$( valid_date $MONTH_BEGIN_JOB $(( $DAY_BEGIN_JOB - $scrumtindays )) $YEAR_BEGIN_JOB )
    or_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
    or_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )
    or_D=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 2) )

# find vertical streching values
if [[ ${RESTART_FLAG} == "FALSE" ]]; then
    ts=$(ncdump -v theta_s croco_ini.nc${agrif_ext}| grep "theta_s = " | cut -d '=' -f 2 | cut -d ' ' -f 2)
    tb=$(ncdump -v theta_b croco_ini.nc${agrif_ext}| grep "theta_b = " | cut -d '=' -f 2 | cut -d ' ' -f 2)
    hc=$(ncdump -v hc croco_ini.nc${agrif_ext}| grep "hc = " | cut -d '=' -f 2 | cut -d ' ' -f 2)
else
    ts=$(ncdump -h croco_ini.nc${agrif_ext}| grep "theta_s = " | cut -d '=' -f 2 | cut -d ' ' -f 2)
    tb=$(ncdump -h croco_ini.nc${agrif_ext}| grep "theta_b = " | cut -d '=' -f 2 | cut -d ' ' -f 2)
    hc=$(ncdump -h croco_ini.nc${agrif_ext}| grep "hc = " | cut -d '=' -f 2 | cut -d ' ' -f 2)
fi
# find recordperdays in online bulk
if [[ ${interponline} -eq 1 ]]; then
    localmth=${cur_M}
    while [ `echo $localmth | cut -b 1` -eq 0 ]; do
        localmth=`echo $localmth | cut -b 2-`
    done
    lmonth=( 1 3 5 7 8 10 12 )
    [[ $((${cur_Y} % 4)) -eq 0  && ( $((${cur_Y} % 100)) -ne 0  ||  $((${cur_Y} % 400)) -eq 0 )]] && { leapyear=1 ;} || { leapyear=0 ;}
    [[ ${lmonth[@]} =~ ${localmth} ]] && { dayinmth=31 ;} || { dayinmth=30 ;}
    [[ ${localmth} == 2 && $leapyear == 1 ]] && { dayinmth=29 ;}
    [[ ${localmth} == 2 && $leapyear == 0 ]] && { dayinmth=28 ;}
    if [[ ${frc_ext} == "ECMWF" ]]; then
        fieldname='T2M'
    elif [[ ${frc_ext} == "AROME" ]];then
        fieldname='AROME'
    else
        fieldname='Temperature_height_above_ground'
    fi
    recpmth=$( ncdump -h "${OCE_FILES_ONLINEDIR}/$fieldname_Y${cur_Y}M${cur_M}.nc" | grep -m 1'time =' | cut -d '=' -f 2 | cut -d ";" -f 1 )
    rpd=$(( ${recpmth} / ${dayinmth}  ))
else
    rpd=4
#
sed -e "s/<ocentimes>/${OCE_NTIMES}/g" -e "s/<ocedt>/${DT_OCE_2}/g"   -e "s/<ocendtfast>/${NDTFAST}/g" \
    -e "s/<theta_s>/${ts}/g" -e "s/<theta_b>/${tb}/g" -e "s/<hc>/${hc}/g" \
    -e "s/<oce_nrst>/${OCE_NTIMES}/g" \
    -e "s|<oce_nhis>|$(( ${oce_his_sec}/ ${DT_OCE_2} ))|g" -e "s|<oce_navg>|$(( ${oce_avg_sec}/${DT_OCE_2} ))|g" \
    -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${MONTH_BEGIN_JOB}/g" -e "s/<rpd>/${rpd}/g" \
    -e "s/<dstart>/${cur_D}/g"  -e "s/<mstart>/${cur_M}/g" -e "s/<ystart>/${cur_Y}/g" \
    -e "s/<dorig>/${or_D}/g"  -e "s/<morig>/${or_M}/g" -e "s/<yorig>/${or_Y}/g" \
    -e "s/<yr2>/${end_Y}/g"             -e "s/<mo2>/${end_M}/g"           \
    ${namfile} > namelist.tmp

if [ ${USE_XIOS_OCE} -eq 1 ]; then
    sed -e "s/<title>/${OCE_OUTPUT_PREFIX}/g" \
    namelist.tmp > tmp$$
else
    sed -e "s/<title>/${CEXPER}/g" \
    namelist.tmp > tmp$$
fi
    mv tmp$$ namelist.tmp
    mv namelist.tmp ${namfile}
#
#cat namelist
#
done

