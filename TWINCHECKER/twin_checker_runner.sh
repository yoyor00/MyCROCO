#!/bin/bash

if [[ $# != 2 ]]; then
    echo "Usage: $0 {exe_1} {exe_2}" 1>&2
    exit 1
fi

exe_1=${1}
exe_2=${2}

export TWIN_CHECKER_NAME=croco
export TWIN_CHECKER_MEETPOINT=/tmp/twin_checker_${TWIN_CHECKER_NAME}

rm -f ${TWIN_CHECKER_MEETPOINT}
mkfifo ${TWIN_CHECKER_MEETPOINT}

echo "============ TWIN CHECK | ${exe_1} | ${exe_2} =================="

# launch master
#TWIN_CHECKER_MODE=master xterm -e gdb ${exe_1} &
# launch slave
#TWIN_CHECKER_MODE=slave xterm -e gdb ${exe_2} &

TWIN_CHECKER_MODE=master xterm -e gdb -ex run ${exe_1} &
pid=$!

TWIN_CHECKER_MODE=slave gdb -ex run ${exe_2}

kill $pid

# wait both
wait
exit $?
