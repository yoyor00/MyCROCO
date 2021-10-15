#set -ue
#set -vx
#
. ${SCRIPTDIR}/caltools.sh
##
##======================================================================
##----------------------------------------------------------------------
##    II. modify namcouple
##----------------------------------------------------------------------
##======================================================================
##
#

echo ' '
echo '-- OASIS inputs --------------'
echo 'fill oasis namcouple'

sed -e "s/<runtime>/$(( ${TOTAL_JOB_DUR} * 86400 ))/g" \
    -e "s/<cpldt>/${CPL_FREQ}/g" \
    ${CPL_NAM_DIR}/${namcouplename} > ./namcouple

### For ATM ###
if [ ${USE_ATM} == 1 ]; then
    for dom in $wrfcpldom ; do
#
        if [ ${DATE_BEGIN_JOB} -eq ${DATE_BEGIN_EXP} ]; then
            file="wrfinput_${dom}"
        else
            file="wrfrst_${dom}*"
        fi 
#
        if [ $dom == "d01" ]; then
            searchf=("<atmdt>" "<atmnx>" "<atmny>")
        else 
            searchf=("<atmdt${dom}>" "<atmnx${dom}>" "<atmny${dom}>")
        fi
#
        dimx=$( ncdump -h  $file  | grep "west_east_stag =" | cut -d ' ' -f 3)
        dimy=$( ncdump -h  $file  | grep "south_north_stag =" | cut -d ' ' -f 3)
        coef=$( ncdump -h  $file  | grep "PARENT_GRID_RATIO =" | cut -d ' ' -f 3)        

        sed -e "s|${searchf[0]}|$(( ${DT_ATM} / ${coef} ))|g"   -e "s/${searchf[1]}/${dimx}/g"   -e "s/${searchf[2]}/${dimy}/g" \
        ./namcouple>tmp$$
        mv tmp$$ namcouple
    done
fi

### For WAV ###
if [ ${USE_WAV} == 1 ]; then
    sed -e "s/<wavdt>/${DT_WAV}/g"   -e "s/<wavnx>/${wavnx}/g"   -e "s/<wavny>/${wavny}/g"  \
    ./namcouple>tmp$$
    mv tmp$$ namcouple
fi

### For OCE ###
if [ ${USE_OCE} == 1 ]; then
    for nn in $( seq 0 ${AGRIFZ} ); do
        if [ ${nn} -gt 0 ];    then
            agrif_ext=".${nn}"
            SUBTIME=$( sed -n -e "$(( 2 * ${nn} )) p" AGRIF_FixedGrids.in | awk '{print $7 }' )
            searchf=("<ocedt${nn}>" "<ocenx${nn}>" "<oceny${nn}>" )
            tsp=$(( ${DT_OCE} / ${SUBTIME} ))
        else
            agrif_ext=""
            searchf=("<ocedt>" "<ocenx>" "<oceny>" )
            tsp=${DT_OCE}
        fi

        dimx=$( ncdump -h  croco_grd.nc${agrif_ext}  | grep "xi_rho =" | cut -d ' ' -f 3)
        dimy=$( ncdump -h  croco_grd.nc${agrif_ext}  | grep "eta_rho =" | cut -d ' ' -f 3)

        sed -e "s/${searchf[0]}/${tsp}/g"   -e "s/${searchf[1]}/$(( ${dimx} - 2 ))/g"   -e "s/${searchf[2]}/$(( ${dimy} - 2 ))/g" \
        ./namcouple>tmp$$
        mv tmp$$ namcouple    
    done
fi

### For TOY ###
if [ ${USE_TOY} -ge 1 ]; then
    for k in `seq 0 $(( ${nbtoy} - 1))` ; do
# 
        if [ ${toytype[$k]} == "oce" ]; then
            dimx=$( ncdump -h  ${toyfile[$k]}  | grep "xi_rho =" | cut -d ' ' -f 3)
            dimy=$( ncdump -h  ${toyfile[$k]}  | grep "eta_rho =" | cut -d ' ' -f 3)
            sed -e "s/<ocedt>/${DT_TOY[$k]}/g"   -e "s/<ocenx>/$(( ${dimx} - 2 ))/g"   -e "s/<oceny>/$(( ${dimy} - 2 ))/g" \
            ./namcouple>tmp$$
            mv tmp$$ namcouple
        fi
#
        if [ ${toytype[$k]} == "atm" ]; then
            dimx=$( ncdump -h  ${toyfile[$k]}  | grep "west_east_stag =" | cut -d ' ' -f 3)
            dimy=$( ncdump -h  ${toyfile[$k]}  | grep "south_north_stag =" | cut -d ' ' -f 3)
            sed -e "s/<atmdt>/${DT_TOY[$k]}/g"   -e "s/<atmnx>/${dimx}/g"   -e "s/<atmny>/${dimy}/g" \
            ./namcouple>tmp$$
            mv tmp$$ namcouple
        fi
#
        if [ ${toytype[$k]} == "wav" ]; then
            dimx=$( ncdump -h  ${toyfile[$k]}  | grep "longitude =" | head -1 | cut -d ' ' -f 3)
            dimy=$( ncdump -h  ${toyfile[$k]}  | grep "latitude =" | head -1 |cut -d ' ' -f 3)
            sed -e "s/<wavdt>/${DT_TOY[$k]}/g"   -e "s/<wavnx>/${dimx}/g"   -e "s/<wavny>/${dimy}/g" \
            ./namcouple>tmp$$
            mv tmp$$ namcouple
        fi
    done
fi

