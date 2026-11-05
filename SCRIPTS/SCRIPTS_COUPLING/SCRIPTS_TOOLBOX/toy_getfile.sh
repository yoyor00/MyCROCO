#!/bin/bash
#-------------------------------------------------------------------------------
#                                                          Configuration files
#-------------------------------------------------------------------------------

# Get grid file

echo 'create input files for TOY'
for k in `seq 0 $(( ${nbtoy} - 1))` ; do
    echo ". ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_toy_files.sh ${toyfile[$k]} toy_${toytype[$k]}.nc ${model_to_toy[$k]} ${toytimerange[$k]}"
    [[ -n ${ncomod} ]] && module load $ncomod
    ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_toy_files.sh ${toyfile[$k]} toy_${toytype[$k]}.nc ${model_to_toy[$k]} ${toytimerange[$k]}
    [[ -n ${ncomod} ]] && module unload $ncomod
done



