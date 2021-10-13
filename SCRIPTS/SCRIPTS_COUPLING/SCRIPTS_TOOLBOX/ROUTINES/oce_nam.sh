#!/bin/bash
set -ue
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

    printf "Computing the origin_date from start_date and scrum_time"
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 ) 
    scrumt=$( ncdump -v scrum_time ${OCE_FILES_DIR}/croco_${ini_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext}| grep "scrum_time = " | cut -d '=' -f 2 | cut -d ' ' -f 2)

    scrumtindays=$(( $scrumt/86400))

    mdy=$( valid_date $MONTH_BEGIN_JOB $(( $DAY_BEGIN_JOB - $scrumtindays )) $YEAR_BEGIN_JOB )
    or_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
    or_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )
    or_D=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 2) )



sed -e "s/<ocentimes>/${OCE_NTIMES}/g" -e "s/<ocedt>/${DT_OCE_2}/g"   -e "s/<ocendtfast>/${NDTFAST}/g" \
    -e "s/<oce_nrst>/${OCE_NTIMES}/g"   -e "s/<oce_nhis>/${oce_nhis}/g" -e "s/<oce_navg>/${oce_navg}/g"     \
    -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"             -e "s/<mo1>/${MONTH_BEGIN_JOB}/g"           \
    -e "s/<dstart>/${DAY_BEGIN_JOB}/g"  -e "s/<mstart>/${MONTH_BEGIN_JOB}/g" -e "s/<ystart>/${YEAR_BEGIN_JOB}/g" \
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


exit






