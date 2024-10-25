#-------------------------------------------------------------------------------
#                                                                   mozaic files
#-------------------------------------------------------------------------------
if [[ ${WEIGHT_FLAG} == TRUE ]]; then
    printf "WEIGHT_FLAG is TRUE. Getting weight files for OASIS"
    weight_files="${weight_o2a} ${weight_a2o}"
    for file in ${weight_files}; do
        ${io_getfile}  ${CPL_FILES_DIR}/${file} .
    done
fi

if [[ ${GRIDS_FLAG} == 1 ]]; then
    printf "GRIDS_FLAG is TRUE. Getting grid files for OASIS"
    ${io_getfile}  ${CPL_FILES_DIR}/${grids_file} grids.nc
fi

