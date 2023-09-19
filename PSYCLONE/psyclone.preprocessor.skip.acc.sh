#!/bin/bash
##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# This script simply aims (for now) on removing the ACC
# from some files
# TODO : this can be removed when work will be finished

##########################################################
# enable strict mode
set -u
set -e

##########################################################
# error
if [[ $# != 2 ]]; then
    echo "Usage: $0 {IN_FILE} {OUT_FILE}" 1>&2
    exit 1
fi

##########################################################
# set vars
PREP_PATH=$(dirname $0)
PREP_IN_FILE=$1
PREP_IN_FILENAME=$(basename ${PREP_IN_FILE})
PREP_IN_DIRNAME=$(dirname ${PREP_IN_FILE})
PREP_OUT_FILE=$2

##########################################################
# switch to cppdef which has not ACC enabled
# TODO: to remove in final version
echo ${PREP_IN_FILENAME} 1>&2
if egrep "^${PREP_IN_FILENAME}\$" ${PREP_PATH}/skip.openacc.rules.lst; then
    cat ${PREP_IN_DIRNAME}/cppdefs.h | sed -e "s/# define OPENACC/# undef OPENACC/g" -e "s/#include \"config_post.h\"/# undef OPENACC/g" > ${PREP_OUT_FILE}.cppdefs-no-acc.h
    sed -e "s#cppdefs.h#${PREP_OUT_FILE}.cppdefs-no-acc.h#g" ${PREP_IN_FILE} > ${PREP_OUT_FILE}
else
    cp ${PREP_IN_FILE} ${PREP_OUT_FILE}
fi
