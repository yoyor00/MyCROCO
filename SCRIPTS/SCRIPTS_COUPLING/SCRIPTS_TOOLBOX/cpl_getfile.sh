#-------------------------------------------------------------------------------
#                                                                   mozaic files
#-------------------------------------------------------------------------------
#${io_getfile} ${INPUTDIRC}/mozaic_o2a_pech12c_pech025w.nc    mozaic_o_ad01.nc 
#${io_getfile} ${INPUTDIRC}/mozaic_a2o_pech025w_pech12c.nc    mozaic_ad01_o.nc
if [[ ${WEIGHT_FLAG} == 1 ]]; then
    weight_files="${weight_o2a} ${weight_a2o}"
    for file in ${weight_files}; do
        ${io_getfile}  ${SCRIPTDIR}/../OASIS_FILES/${file} .
    done
fi
