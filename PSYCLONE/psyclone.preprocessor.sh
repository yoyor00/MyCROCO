#!/bin/bash

#!/bin/bash

###########################################################
# This script simply aims at calling psyclone on the
# files which need to be applied depending on the
# tuning configuration

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
PREP_OUT_FILE=$2

##########################################################
# switch to cppdef which has not ACC enabled
# TODO: to remove in final version
if egrep "^${PREP_IN_FILENAME}" ${PREP_PATH}/psyclone.rules.lst; then
    psyclone -api nemo -l output -s ${PREP_PATH}/scripts/$(egrep "^${PREP_IN_FILENAME}" ${PREP_PATH}/psyclone.rules.lst | cut -f 2) -opsy ${PREP_OUT_FILE} ${PREP_IN_FILE}
else
    cp ${PREP_IN_FILE} ${PREP_OUT_FILE}
fi
