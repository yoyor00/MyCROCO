#!/bin/bash

module load $ncomod

if [[ ${bdy_ext} == *'clm'* ]]; then
    bryfile="croco_clm.nc"
    timevar="tclm_time"
else
    bryfile="croco_bry.nc"
    timevar="bry_time"
fi

for i in `seq 0 $(( ${JOB_DUR_MTH}-1 ))`; do
  if [ ${JOB_DUR_MTH} -eq 1 ]; then
      cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
      cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
      ln -sf ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc ${bryfile}
  elif [ ${i} -eq 0 ]; then
      cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
      cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )

      varlist='spherical Vtransform Vstretching tstart theta_s theta_b Tcline hc sc_r sc_w Cs_r Cs_w' # One dimension stuffs that don't change 
      for varn in ${varlist} ; do
          ncks -A -v ${varn} ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc ${bryfile}
      done
      string=$( ncdump -h ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc | grep double )
      ns=$( ncdump -h ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc | grep -c double )

      for j in `seq ${vartopass} $ns`; do
          if [ ${j} -eq ${vartopass} ] ; then
              var=${timevar}
              dimt=${timevar}
          else
              end=$( echo $string | cut -d';' -f ${j})
              var1=$( echo $end | cut -d'(' -f 1)
              var=$(echo $var1 | cut -d' ' -f 2)
              dimt1=$( echo $end | cut -d'(' -f 2)
              dimt=$(echo $dimt1 | cut -d',' -f 1)
          fi

          ncks -O -F -v ${var} -d $dimt,1,2 ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc out_${var}_${i}.nc
          ncks -O --mk_rec_dmn $dimt out_${var}_${i}.nc out_${var}.nc
          \rm -f out_${var}_${i}.nc
      done

  else
      mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
      cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
      cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )

      string=$( ncdump -h ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc | grep double )
      ns=$( ncdump -h ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc | grep -c double )
# Find how many var to pass before arriving on 2-3d var
      echo "$( ncdump -h ${OCE_FILES_DIR}/croco_${bry_ext}_Y${cur_Y}M${cur_M}.nc  | grep double )"  > text.tmp
      cnt=0
      while [[ $cnt -le $ns ]]; do
          cnt=$(( $cnt +1 ))
          line=$( sed -n "${cnt}p" text.tmp )
          tmpval=$( echo "$line" | grep ',' | wc -l )
          if [[ ${tmpval} -gt 0 ]]; then
              break
          fi
      done
      rm -rf text.tmp
      vartopass=$(( $cnt - 1 ))
#
      for j in `seq ${vartopass} $ns`; do
          if [ ${j} -eq ${vartopass} ] ; then
              var=${timevar}
              dimt=${timevar}
          else
              end=$( echo $string | cut -d';' -f ${j})
              var1=$( echo $end | cut -d'(' -f 1)
              var=$(echo $var1 | cut -d' ' -f 2)
              dimt1=$( echo $end | cut -d'(' -f 2)
              dimt=$(echo $dimt1 | cut -d',' -f 1)
          fi

          if [ $i -eq $(( ${JOB_DUR_MTH}-1 )) ]; then
              ncks -O -F -v ${var} -d $dimt,2,3 ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc out_${var}_${i}.nc
              ncrcat -O out_${var}.nc out_${var}_${i}.nc out_${var}.nc
              \rm -f out_${var}_${i}.nc
              ncks -O --fix_rec_dmn $dimt out_${var}.nc out_${var}.nc
              ncks -A out_${var}.nc ${bryfile} ; \rm -f out_${var}.nc
              [ ${j} -eq 21 ] && ncks -A -v tend ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc ${bryfile}

          else
              ncks -O -F -v ${var} -d $dimt,2 ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc out_${var}_${i}.nc
              ncrcat -O out_${var}.nc out_${var}_${i}.nc out_${var}.nc
              \rm -f out_${var}_${i}.nc
          fi
      done
  fi
done

if [ ${JOB_DUR_MTH} -eq 0 ] ; then
  printf "Job duration is less than a month ---> Using netcdf of the current month\n"
  cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
  cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )    
  ln -sf ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc ${bryfile}
fi

module unload $ncomod
