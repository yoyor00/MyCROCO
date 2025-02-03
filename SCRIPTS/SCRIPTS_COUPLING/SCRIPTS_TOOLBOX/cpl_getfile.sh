#-------------------------------------------------------------------------------
#                                                                   mozaic files
#-------------------------------------------------------------------------------
if [[ ${WEIGHT_FLAG} == TRUE ]]; then
    printf "WEIGHT_FLAG is TRUE. Getting weight files for OASIS"
    weight_files="${weight_oce} ${weight_atm} ${weight_wav}"
    for file in ${weight_files}; do
        ${io_getfile}  ${CPL_FILES_DIR}/${file} .
    done
fi

if [[ ${GRIDS_FLAG} == 1 ]]; then
    printf "GRIDS_FLAG is TRUE. Getting grid files for OASIS"
    ${io_getfile}  ${CPL_FILES_DIR}/${grids_file} grids.nc
    ${io_getfile}  ${CPL_FILES_DIR}/${masks_file} masks.nc
    ${io_getfile}  ${CPL_FILES_DIR}/${areas_file} areas.nc
fi

