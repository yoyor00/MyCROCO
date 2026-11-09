#!/bin/bash
# Experimental OASIS toy-model build.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
# shellcheck source=env.sh
source "$HERE/env.sh"

SRC="$REPO_ROOT/SCRIPTS/SCRIPTS_COUPLING/TOY_IN"
BUILD="$HERE/build/toy"
PREFIX="$HERE/install/toy"
JOBS="${JOBS:-$(nproc)}"

[[ -d "$SRC" ]] || { echo "missing TOY_IN sources: $SRC" >&2; exit 1; }
[[ -d "$CROCO_PRISM_ROOT_DIR" ]] || { echo "missing OASIS install: $CROCO_PRISM_ROOT_DIR" >&2; exit 1; }

rm -rf "$BUILD" "$PREFIX"
mkdir -p "$BUILD" "$PREFIX" "$HERE/logs"

for toytype in atm oce wav; do
    echo "Building toy_${toytype}"
    rm -rf "$BUILD/$toytype"
    mkdir -p "$BUILD/$toytype"
    cp -a "$SRC/." "$BUILD/$toytype/"
    (
        cd "$BUILD/$toytype"
        ln -sf "read_namelist.F90_${toytype}" read_namelist.F90
        ln -sf "toy_model.F90_${toytype}" toy_model.F90
        make clean >/dev/null 2>&1 || true
        make -j "$JOBS" \
            F90=mpif90 \
            NETCDF_CONFIG=nf-config \
            CPL="$CROCO_PRISM_ROOT_DIR" \
            F90COMP="-O2 -ffree-line-length-none -fallow-argument-mismatch" \
            F90LINK=""
    ) > "$HERE/logs/toy-${toytype}-build.log" 2>&1
    cp "$BUILD/$toytype/toy_model" "$PREFIX/toy_${toytype}"
done

echo "TOY models installed in $PREFIX"
