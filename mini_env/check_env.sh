#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=env.sh
source "$HERE/env.sh"

echo "MINI_ENV_ROOT=$MINI_ENV_ROOT"
echo "XIOS_ROOT_DIR=$XIOS_ROOT_DIR"
echo "CROCO_PRISM_ROOT_DIR=$CROCO_PRISM_ROOT_DIR"
echo
echo "Compilers:"
for tool in gfortran mpif90 mpicc mpicxx cpp make; do
    printf '  %-10s %s\n' "$tool" "$(command -v "$tool" || echo MISSING)"
done
echo
echo "NetCDF:"
nf-config --version 2>/dev/null || echo "  nf-config MISSING"
nc-config --version 2>/dev/null || echo "  nc-config MISSING"
echo
echo "XIOS:"
[[ -x "$XIOS_ROOT_DIR/bin/xios_server.exe" ]] && echo "  OK: $XIOS_ROOT_DIR/bin/xios_server.exe" || echo "  MISSING xios_server.exe"
[[ -f "$XIOS_ROOT_DIR/lib/libxios.a" ]] && echo "  OK: $XIOS_ROOT_DIR/lib/libxios.a" || echo "  MISSING libxios.a"
echo
echo "OASIS:"
for lib in libpsmile.MPI1.a libmct.a libmpeu.a libscrip.a; do
    [[ -f "$CROCO_PRISM_ROOT_DIR/lib/$lib" ]] && echo "  OK: $lib" || echo "  MISSING: $lib"
done
for dir in psmile.MPI1 mct; do
    [[ -d "$CROCO_PRISM_ROOT_DIR/build/lib/$dir" ]] && echo "  OK: build/lib/$dir" || echo "  MISSING: build/lib/$dir"
done
echo
echo "Python/BENCH:"
"$HERE/.venv/bin/python" -c 'import json5, matplotlib, netCDF4, xarray, SCons; print("  OK: only local .venv")' 2>/dev/null || echo "  missing; run ./setup_python.sh"
