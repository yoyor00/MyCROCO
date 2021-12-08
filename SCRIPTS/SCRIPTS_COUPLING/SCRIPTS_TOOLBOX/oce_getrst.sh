#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------

if [[ ${RESTART_FLAG} == "FALSE" ]]; then
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
    for nn in $( seq 0 ${AGRIFZ} ); do
       if [ ${nn} -gt 0 ]; then
           agrif_ext=".${nn}"
       else
           agrif_ext=""
       fi
           ${io_getfile} ${OCE_FILES_DIR}/croco_${ini_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext} croco_ini.nc${agrif_ext}
    done
else
    for nn in $( seq 0 ${AGRIFZ} ); do
       if [ ${nn} -gt 0 ]; then
           agrif_ext=".${nn}"
       else
           agrif_ext=""
       fi
       ${io_getfile} ${RESTDIR_IN}/croco_rst_${DATE_END_JOBm1}.nc${agrif_ext} croco_ini.nc${agrif_ext}
    done
fi
#
if [[ ${USE_ATM} == 1 ]]; then
    if [[ -f "${OCE_FILES_DIR}/coupling_masks*.nc" ]]; then
        cp ${OCE_FILES_DIR}/coupling_masks*.nc .
    elif [[ $( echo "$wrfcpldom" | wc -w ) >1 ]]; then
        module load ${ncomod}
        for nn in `seq 0 $AGRIFZ`; do
            if [ ${nn} -gt 0 ]; then
                agrif_ext=".${nn}"
            else
                agrif_ext=""
            fi
            [ -f coupling_masks${nn}.nc ] && rm -rf coupling_masks${nn}.nc
            echo "Creating coupling mask for CROCO"
            ncks -v lon_rho,lat_rho croco_grd.nc${agrif_ext} coupling_masks${nn}.nc
            for atmdom in $wrfcpldom; do
                if [[ ${RESTART_FLAG} == "FALSE" ]]; then
                    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
                    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 ) 
                    filefrom="${ATM_FILES_DIR}/wrfinput_d0"
                else
                    cur_Y=$( echo DATE_END_JOBm1 | cut -c 1-4 )
                    cur_M=$( echo DATE_END_JOBm1 | cut -c 5-6 )
                    filefrom="${RESTDIR_IN}/wrfrst_d0"
                fi
                maxdom=$( echo "$wrfcpldom" | wc -w )
                dom=$( echo "$atmdom" | cut -c 3 )
                ncap2 -A -v -s "cplmask${dom}=mask_rho*0+1" croco_grd.nc${agrif_ext} coupling_masks${nn}.nc
                for sdloop in `seq 1 $maxdom`; do
                    ncap2 -O -v -s 'latmin=XLAT.min();latmax=XLAT.max();lonmin=XLONG.min();lonmax=XLONG.max()' ${filefrom}${sdloop}_${cur_Y}_${cur_M}* tmp.nc
                    lonmin=$( ncdump -v lonmin tmp.nc  | grep "lonmin =" | cut -d ' ' -f 4)
                    latmin=$( ncdump -v latmin tmp.nc  | grep "latmin =" | cut -d ' ' -f 4)
                    lonmax=$( ncdump -v lonmax tmp.nc  | grep "lonmax =" | cut -d ' ' -f 4)
                    latmax=$( ncdump -v latmax tmp.nc  | grep "latmax =" | cut -d ' ' -f 4)
                    \rm tmp.nc
                    if [[ $sdloop > $dom ]]; then
                        ncap2 -O -s "var_tmp=cplmask${dom}; where( lat_rho > $latmin && lon_rho > $lonmin && lat_rho < $latmax && lon_rho < $lonmax ) var_tmp=0; cplmask${dom}=var_tmp" coupling_masks${nn}.nc coupling_masks${nn}.nc
                    elif [[ $sdloop -le $dom ]]; then
                        ncap2 -O -s "var_tmp=cplmask${dom}; where( lat_rho < $latmin || lon_rho < $lonmin || lat_rho > $latmax || lon_rho > $lonmax ) var_tmp=0; cplmask${dom}=var_tmp" coupling_masks${nn}.nc coupling_masks${nn}.nc
                    fi

                    ncks -O -v var_tmp -x coupling_masks${nn}.nc coupling_masks${nn}.nc
                done
            done
        done
        module unload ${ncomod}
    fi
fi

