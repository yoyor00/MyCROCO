#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------
if [[ ${RESTART_FLAG} == "FALSE" ]]
then

 filelist='wrfinput_d01' 
 if [ $NB_dom -ge 2 ] ; then
  filelist="$filelist wrfinput_d02"
  if [ $NB_dom -eq 3 ] ; then
   filelist="$filelist wrfinput_d03"
  fi
 fi
 for file in $filelist
  do
   echo "ln -sf ${ATM_FILES_DIR}${file} ./$file"
   ln -sf ${ATM_FILES_DIR}/${file} ./$file
  done

  for dom in $wrfcpldom ; do
     if [[ ${dom} == "d01" ]]; then
         echo 'set CPLMASK to 1 in coupled domain '$dom
         echo "ncap2 -O -s "CPLMASK(:,0,:,:)= LU_INDEX == 17" ./wrfinput_$dom ./wrfinput_$dom"
         module load $ncomod
         ncap2 -O -s "CPLMASK(:,0,:,:)= LU_INDEX == 17" ./wrfinput_$dom ./wrfinput_$dom
         module unload $ncomod
     else
         module load $ncomod
         echo 'set CPLMASK to 1 in coupled domain '$dom
         num_ext_mod=$( ncdump -h wrfinput_d01 | grep "num_ext_model_couple_dom_stag = " | cut -d ' ' -f 3)
         echo "Increase size of mum_ext_model by one for ${dom} (in case some domains already exist)"
         ncpdq -O -v CPLMASK -a num_ext_model_couple_dom_stag,Time wrfinput_d01 tmp.nc 
         cp tmp.nc tmp2.nc
         ncrcat -O tmp.nc tmp2.nc tmp2.nc 
         ncpdq -O -a Time,num_ext_model_couple_dom_stag tmp2.nc tmp2.nc 
         ncks -A wrfinput_d01 tmp2.nc
         mv tmp2.nc wrfinput_d01
         rm -rf tmp.nc
         num_ext_mod=$( ncdump -h wrfinput_d01 | grep "num_ext_model_couple_dom_stag = " | cut -d ' ' -f 3)
         
         echo "Find limits for domain $dom"
         ncap2 -O -v -s 'latmin=XLAT.min();latmax=XLAT.max();lonmin=XLONG.min();lonmax=XLONG.max()' wrfinput_${dom} tmp.nc
	 lonmin=$( ncdump -v lonmin tmp.nc  | grep "lonmin =" | cut -d ' ' -f 4)
         latmin=$( ncdump -v latmin tmp.nc  | grep "latmin =" | cut -d ' ' -f 4)
	 lonmax=$( ncdump -v lonmax tmp.nc  | grep "lonmax =" | cut -d ' ' -f 4)
	 latmax=$( ncdump -v latmax tmp.nc  | grep "latmax =" | cut -d ' ' -f 4)
         rm -rf tmp.nc
	 printf "Limits for domain ${dom} are:\n Lon min:$lonmin \n Lat min:$latmin \n Lon max:$lonmax \n Lat max:$latmax \n"
         ncap2 -F -O -s "var_tmp=CPLMASK(:,0,:,:); where( XLAT < $latmin || XLONG < $lonmin || XLAT > $latmax || XLONG > $lonmax ) var_tmp=0; CPLMASK(:,${num_ext_mod},:,:)=var_tmp" wrfinput_d01 wrfinput_d01.tmp
	 ncks -O -v var_tmp -x wrfinput_d01.tmp wrfinput_d01         
         module unload $ncomod
     fi

   done

else
    touch ls_l/getfile_atm_restarts.txt
    for file in `${MACHINE_STOCKAGE} ls ${RESTDIR_IN}/wrfrst_d0?_*`
# for i in ${RESTDIR_IN}/wrfrst_d01_*
    do
	${io_getfile} ${file} . >> ls_l/getfile_atm_restarts.txt
    done
fi
