#!/bin/bash
#set -ue
#set -vx
#
. ${SCRIPTDIR}/caltools.sh
##
##======================================================================
##----------------------------------------------------------------------
##    II. modify namelist
##----------------------------------------------------------------------
##======================================================================
##
#

if [[ ${RESTART_FLAG} == "FALSE" ]]; then
  rst="false"
else
  rst="true"
fi

sed -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"   -e "s/<yr2>/${YEAR_END_JOB}/g"  \
    -e "s/<mo1>/${MONTH_BEGIN_JOB}/g"   -e "s/<mo2>/${MONTH_END_JOB}/g"  \
    -e "s/<dy1>/${DAY_BEGIN_JOB}/g"   -e "s/<dy2>/${DAY_END_JOB}/g"  \
    -e "s/<hr1>/0/g"   -e "s/<hr2>/24/g"  \
    -e "s/<rst>/$rst/g"              -e "s/<rst_int_h>/$(( ${TOTAL_JOB_DUR} * 24 ))/g"            \
    -e "s/<his_int_h>/${atm_his_h}/g"         -e "s/<his_nb_out>/${atm_his_frames}/g"    \
    -e "s/<xtrm_int_m>/${atm_diag_int_m}/g"   -e "s/<xtrm_nb_out>/${atm_diag_frames}/g"  \
    -e "s/<nproc_x>/${atm_nprocX}/g"            -e "s/<nproc_y>/${atm_nprocY}/g"             \
    -e "s/<niotaskpg>/${atm_niotaskpg}/g"       -e "s/<niogp>/${atm_niogp}/g"                \
    -e "s/time_step                           =.*/time_step                           =${DT_ATM} ,/g" \
    -e "s/<max_domains>/${NB_dom}/g" \
    -e "s/<interval_s>/${interval_seconds}/g"  -e "s/<sst_int_m>/${auxinput4_interval}/g" \
    -e "s/<nbmetsoil>/4/g" \
    $ATM_NAM_DIR/${atmnamelist} > ./namelist.input

for domnb in `seq 1 $NB_dom` ; do
    dom="d$( printf "%02d" ${domnb})"
    if [[ ${RESTART_FLAG} == "FALSE" ]]; then
        file="wrfinput_${dom}"
    else
        file="wrfrst_${dom}*"
    fi
#
    searchf=("<xdim_${dom}>" "<ydim_${dom}>" "<nbvertlev>" "<dx_${dom}>" "<dy_${dom}>" "<i_str_${dom}>" "<j_str_${dom}>" "<coef_${dom}>")
#
    dimx=$( ncdump -h  $file  | grep "west_east_stag =" | cut -d ' ' -f 3)
    dimy=$( ncdump -h  $file  | grep "south_north_stag =" | cut -d ' ' -f 3)
    dimz=$( ncdump -h  $file  | grep "bottom_top_stag =" | cut -d ' ' -f 3)
    dx=$( ncdump -h  $file | grep "DX =" | cut -d ' ' -f 3 | cut -d '.' -f 1)
    dy=$( ncdump -h  $file  | grep "DY =" | cut -d ' ' -f 3 | cut -d '.' -f 1)
    istrt=$( ncdump -h  $file  | grep "I_PARENT_START =" | cut -d ' ' -f 3)
    jstrt=$( ncdump -h  $file  | grep "J_PARENT_START =" | cut -d ' ' -f 3)
    coef=$( ncdump -h  $file  | grep "PARENT_GRID_RATIO =" | cut -d ' ' -f 3)

    sed -e "s/${searchf[0]}/${dimx}/g" \
        -e "s/${searchf[1]}/${dimy}/g" \
        -e "s/${searchf[2]}/${dimz}/g" \
        -e "s/${searchf[3]}/${dx}/g" \
        -e "s/${searchf[4]}/${dy}/g" \
        -e "s/${searchf[5]}/${istrt}/g" \
        -e "s/${searchf[6]}/${jstrt}/g" \
        -e "s/${searchf[7]}/${coef}/g" \
        ./namelist.input > ./namelist.input.tmp
        mv namelist.input.tmp namelist.input
    chmod 755 namelist.input
done

if [[ ${nestfeedback} == "FALSE" ]]; then
    sed -e "s/feedback                            = 1,/feedback                            = 0,/g" \
    ./namelist.input > ./namelist.input.tmp
    mv namelist.input.tmp namelist.input
fi
   
maxatmdom=$( echo "$wrfcpldom" | wc -w )
if [[ $maxatmdom == 1 && $AGRIFZ > 1 ]];then 
    numextmod=$(( $AGRIFZ + 1 ))
else
    [[ $maxatmdom > $(( $AGRIFZ +1 )) ]] && { numextmod=$maxatmdom ;} || { numextmod=$(( $AGRIFZ + 1 )) ;}
fi
sed -e "s/num_ext_model_couple_dom            = 1,/num_ext_model_couple_dom            = $numextmod,/g" \
./namelist.input > ./namelist.input.tmp
mv namelist.input.tmp namelist.input

if [ $USE_WAV -eq 1 ] || [ $USE_TOYWAV -eq 1 ]; then
    cplwavdom=""
    for dom in `seq 1 $NB_dom`; do
        cplwavdom="$cplwavdom 5,"
    done
    sed -e "s/isftcflx                            = 0,/isftcflx                            = $cplwavdom/g" \
    ./namelist.input > ./namelist.input.tmp
    mv namelist.input.tmp namelist.input
fi

sed -e "s/<xdim_d.*>/1/g"  -e "s/<ydim_d.*>/1/g"\
    -e "s/<dx_d.*>/1/g"  -e "s/<dy_d.*>/1/g" \
    -e "s/<i_str_d.*>/1/g"  -e "s/<j_str_d.*>/1/g" \
    -e "s/<coef_d.*>/1/g"  -e "s/<coef_d.*>/1/g" \
./namelist.input > ./namelist.input.tmp
mv namelist.input.tmp namelist.input
chmod 755 namelist.input

