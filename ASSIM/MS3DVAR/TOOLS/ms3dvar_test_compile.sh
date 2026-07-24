#!/bin/bash
# ms3dvar_test_compile.sh -- end-to-end compilation check for MS3DVAR
#
# Drives the same interactive procedure documented in docs/ms3dvar/INSTALL.md:
#   1. create_config.bash -o ms3dvar   (sets up a work config with pre-patched jobcomp)
#   2. jobcomp <scale...>            (builds the requested scale(s))
#
# Does NOT use ms3dvar_compile_scale.sh (that one is for WOOM only).
#
# Usage:
#   ./ms3dvar_test_compile.sh [SCALE...] [--keep] [--work-dir DIR] [-- JOBCOMP_OPTS...]
#
#   SCALE...    lr | mr | filter | ms | all   (default: all)
#   --keep        keep the generated work config instead of deleting it
#   --work-dir    parent directory for the generated config
#                 (default: a fresh subdir of ${HOME}/Work/CROCO)
#   -- ...        extra options forwarded to jobcomp (e.g. --fc ifort --jobs 4)

set -e

DEFAULT_WORK_ROOT="${HOME}/Work/CROCO"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CROCO_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

KEEP=0
WORK_PARENT=""
SCALES=()
JOBCOMP_OPTS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)     KEEP=1; shift ;;
        --work-dir) WORK_PARENT="$2"; shift 2 ;;
        --)         shift; JOBCOMP_OPTS+=("$@"); break ;;
        -h|--help)
            sed -n '2,20p' "$0"
            exit 0
            ;;
        *) SCALES+=("$1"); shift ;;
    esac
done
[ ${#SCALES[@]} -eq 0 ] && SCALES=(all)

if [ -z "$WORK_PARENT" ]; then
    mkdir -p "${DEFAULT_WORK_ROOT}"
    WORK_PARENT="$(mktemp -d "${DEFAULT_WORK_ROOT}/ms3dvar_test_compile.XXXXXX")"
fi
CONFIG_NAME="ms3dvar_test"
CONFIG_DIR="${WORK_PARENT}/${CONFIG_NAME}"

cleanup() {
    if [ $KEEP -eq 0 ]; then
        rm -rf "${WORK_PARENT}"
    else
        echo "Work config kept at: ${CONFIG_DIR}"
    fi
}
trap cleanup EXIT

echo "========================================================================"
echo "Step 1: create_config.bash -o ms3dvar"
echo "========================================================================"
"${CROCO_DIR}/create_config.bash" -f -s "${CROCO_DIR}" -d "${WORK_PARENT}" -w "${WORK_PARENT}" -n "${CONFIG_NAME}" -o ms3dvar

echo ""
echo "========================================================================"
echo "Step 2: jobcomp ${SCALES[*]} ${JOBCOMP_OPTS[*]}"
echo "========================================================================"
cd "${CONFIG_DIR}/ASSIM/MS3DVAR"
./jobcomp "${SCALES[@]}" "${JOBCOMP_OPTS[@]}"
