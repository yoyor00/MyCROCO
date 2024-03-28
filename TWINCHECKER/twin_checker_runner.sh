#!/bin/bash

set -e
set -u


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

dir_1=$(dirname ${exe_1})
dir_2=$(dirname ${exe_2})

function run_gdb()
{
    cd ${dir_1}
    TWIN_CHECKER_MODE=master xterm -e gdb -ex run ${exe_1} &
    pid=$!

    cd ${dir_2}
    TWIN_CHECKER_MODE=slave gdb -ex run ${exe_2}
}

function run_valgrind()
{
    cd ${dir_1}
    TWIN_CHECKER_MODE=master xterm -e valgrind ${exe_1} &
    pid=$!

    cd ${dir_2}
    TWIN_CHECKER_MODE=slave valgrind ${exe_2}
}

function run_std()
{
    cd ${dir_1}
    TWIN_CHECKER_MODE=master xterm -e ${exe_1} &
    pid=$!

    cd ${dir_2}
    TWIN_CHECKER_MODE=slave ${exe_2}
}

run_gdb
#run_valgrind
#run_std

kill $pid

# wait both
wait
exit $?
