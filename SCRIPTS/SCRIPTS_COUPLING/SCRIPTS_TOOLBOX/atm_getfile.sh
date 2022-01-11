#-------------------------------------------------------------------------------
#                                                          Configuration files
#-------------------------------------------------------------------------------
# link data files necessary for running wrf in a dedicated directory $wrf/data
if [ ! -d ${ATM_EXE_DIR}/../data ] ; then
 mkdir ${ATM_EXE_DIR}/../data
 ln -s ${ATM_EXE_DIR}/../run/* ${ATM_EXE_DIR}/../data/.
 # remove executables that could exist and namelist file
 rm -f ${ATM_EXE_DIR}/../data/*.exe
 rm -f ${ATM_EXE_DIR}/../data/namelist.input*
fi
echo 'link wrf data files'
echo "ln -sf ${ATM_EXE_DIR}/../data/* ."
ln -sf ${ATM_EXE_DIR}/../data/* .

#-------------------------------------------------------------------------------
#                                                          BDY
#-------------------------------------------------------------------------------
if [[ ${JOB_DUR_MTH} -le 1 ]]; then
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
    ${io_getfile} ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} wrfbdy_d01
else
    module load ${ncomod}
    for i in `seq 0 $(( ${JOB_DUR_MTH} ))`; do
        mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
        cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
        cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )
        if [[ $i == 0 ]]; then
            ncrcat ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} wrfbdy_d01
        else
            [[ ! -f ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} to build atm bdy file."; exit ;} 
            ncrcat -O wrfbdy_d01 ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} wrfbdy_d01
        fi
    done
    module unload ${ncomod}
fi
    
#-------------------------------------------------------------------------------
#                                            Forcing fields (interannual case)
#-------------------------------------------------------------------------------
filelist='wrflowinp_d01'
 if [ $NB_dom -ge 2 ] ; then
  filelist="$filelist wrflowinp_d02"
  if [ $NB_dom -eq 3 ] ; then
   filelist="$filelist wrflowinp_d03"
  fi
 fi

for file in ${filelist}; do
    if [[ ${JOB_DUR_MTH} -le 1 ]]; then
        cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
        cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
        ${io_getfile} ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file} 
    else
        module load ${ncomod}
        for i in `seq 0 $(( ${JOB_DUR_MTH}-1 ))`; do
            mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
            cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
            cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )
            if [[ $i == 0 ]]; then
                ncrcat -O -F -d Time,1,-2 ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
            elif [[ $i == ${JOB_DUR_MTH} ]]; then
                ncrcat  -O ${file} ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
            else
                ncrcat -O -F -d Time,1,-2 ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}.tmp
                ncrcat -O ${file} ${file}.tmp ${file}
                rm -rf ${file}.tmp
            fi
        done
        module unload ${ncomod}
    fi
done
