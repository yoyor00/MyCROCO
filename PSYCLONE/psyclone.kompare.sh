#!/bin/bash

###########################################################
# This scripts aimed at helping to compare the various
# openacc version (manual versus psyclone) to ease debugging.
#
# Usage :
#   1) Compare no-acc with psyclone-acc 
#      ----------------------------------
#      cd build-psyclone
#      ../PSYCLONE/psyclone.kompare.sh .
#      ----------------------------------
#  2) Compare orig-acc with psyclone-acc
#      ----------------------------------
#      cd build-psyclone
#      ../PSYCLONE/psyclone.kompare.sh . ../build-native-acc
#      ----------------------------------

###########################################################
# set strict
set -u
set -e

###########################################################
# handle options
if [[ $# < 1 ]]; then
	echo "Usage : $0 {BUILD_DIR} [NATIVE_BUILD_DIR]"
	exit 1
fi

###########################################################
# get path
BUILD_DIR=$1
NATIVE_BUILD_DIR='OFF'
if [[ $# == 2 ]]; then
	NATIVE_BUILD_DIR=${2}
fi

# calc paths
PREP_DIR=${BUILD_DIR}/OCEAN/prepared_sources
NATIVE_PREP_DIR=${NATIVE_BUILD_DIR}/OCEAN/prepared_sources

# make abos
PREP_DIR=$(cd "${PREP_DIR}" && pwd -P)
if [[ ${NATIVE_BUILD_DIR} != 'OFF' ]]; then
	NATIVE_PREP_DIR=$(cd "${NATIVE_PREP_DIR}" && pwd -P)
fi

# check
if [[ ! -d ${PREP_DIR} ]]; then
	echo "Non existing directory : ${PREP_DIR} !"
	exit 1
fi

# check
if [[ ${NATIVE_BUILD_DIR} != 'OFF' && ! -d ${NATIVE_PREP_DIR} ]]; then
	echo "Non existing directory : ${NATIVE_PREP_DIR} !"
	exit 1
fi

# rm
rm -f ${PREP_DIR}/compare/no-acc/*
rm -f ${PREP_DIR}/compare/acc/*
rm -f ${PREP_DIR}/compare/orig-acc/*

# create dests
mkdir -p ${PREP_DIR}/compare
mkdir -p ${PREP_DIR}/compare/no-acc
mkdir -p ${PREP_DIR}/compare/acc
mkdir -p ${PREP_DIR}/compare/orig-acc

# function
function run()
{
	echo "-- $@"
	"$@"
}

# copy files (acc)
for file in $(find ${PREP_DIR} -iname "*.psyclone.F90"); do
	orig_name=$(basename $file | sed -e 's/\.no-acc\.cpp\.mpc\.loops\.mpc\.fix\.psyclone//g')
	fname=$(basename $file)
	run ln -sf ${PREP_DIR}/${fname} ${PREP_DIR}/compare/acc/${orig_name}
done

# copy files (no acc)
for file in $(find ${PREP_DIR} -iname "*.psyclone.dummy.F90"); do
	orig_name=$(basename $file | sed -e 's/\.no-acc\.cpp\.mpc\.loops\.mpc\.fix\.psyclone\.dummy//g')
	fname=$(basename $file)
	run ln -sf ${PREP_DIR}/${fname} ${PREP_DIR}/compare/no-acc/${orig_name}
done

# copy files (acc)
if [[ ${NATIVE_BUILD_DIR} != 'OFF' ]]; then
	for file in $(find ${NATIVE_PREP_DIR} -iname "*.fix.F"); do
		orig_name=$(basename $file | sed -e 's/\.cpp\.mpc\.loops\.mpc\.fix//g')
		fname=$(basename $file)
		if [[ -f ${PREP_DIR}/compare/acc/${orig_name}90 ]]; then
			#ln -sf ${NATIVE_PREP_DIR}/${fname} ${PREP_DIR}/compare/orig-acc/${orig_name}90
			run psyclone -api nemo -l output -opsy ${PREP_DIR}/compare/orig-acc/${orig_name}90 ${NATIVE_PREP_DIR}/${fname} &
		fi
	done
	wait
fi

# kaunch
if [[ ${NATIVE_BUILD_DIR} != 'OFF' ]]; then
	set -x
	kompare ${PREP_DIR}/compare/orig-acc ${PREP_DIR}/compare/acc
else
	set -x
	kompare ${PREP_DIR}/compare/no-acc ${PREP_DIR}/compare/acc
fi
