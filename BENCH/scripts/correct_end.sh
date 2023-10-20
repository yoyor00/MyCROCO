#!/bin/bash

# To check that croco has correct end
# => because with nvfortran return ieee_inexact even on success

# enforce strict
set -u

# vars
logfile=$(mktemp)

# run
"$@" 2>&1 | tee ${logfile}

# check
grep -q 'MAIN: DONE' ${logfile}
status=$?

# remove log
rm -f ${logfile}

# display error
if [[ ${status} != 0 ]]; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 1>&2
    echo "CRCO didn't reached end of execution !" 1>&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 1>&2
    exit 1
fi
