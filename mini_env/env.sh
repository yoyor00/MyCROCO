#!/bin/bash
# Source this file to use the experimental XIOS/OASIS installation.
#
#   source mini_env/env.sh
#
# This intentionally does not install anything system-wide.

if [[ -n "${BASH_SOURCE[0]:-}" ]]; then
    MINI_ENV_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
    MINI_ENV_ROOT="$(pwd)"
fi
export MINI_ENV_ROOT
export MINI_ENV_PREFIX="${MINI_ENV_PREFIX:-$MINI_ENV_ROOT/install}"

export XIOS_ROOT_DIR="${XIOS_ROOT_DIR:-$MINI_ENV_PREFIX/xios}"
export CROCO_PRISM_ROOT_DIR="${CROCO_PRISM_ROOT_DIR:-$MINI_ENV_PREFIX/oasis}"

# Names used by older CROCO coupling scripts.
export XIOS="${XIOS_ROOT_DIR}"
export CPL="${CROCO_PRISM_ROOT_DIR}"
export PRISM_ROOT_DIR="${CROCO_PRISM_ROOT_DIR}"
export CROCO_XIOS_ROOT_DIR="${XIOS_ROOT_DIR}"
export CPL_NAM_DIR="${CPL_NAM_DIR:-$MINI_ENV_ROOT/../SCRIPTS/SCRIPTS_COUPLING/OASIS_IN}"
export CPL_FILES_DIR="${CPL_FILES_DIR:-$MINI_ENV_ROOT/work/OASIS_FILES}"

if [[ -x "$MINI_ENV_ROOT/.venv/bin/python" ]]; then
    export PATH="$MINI_ENV_ROOT/.venv/bin:$PATH"
fi
export PATH="$XIOS_ROOT_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$XIOS_ROOT_DIR/lib:$CROCO_PRISM_ROOT_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

mkdir -p "$CPL_FILES_DIR"
