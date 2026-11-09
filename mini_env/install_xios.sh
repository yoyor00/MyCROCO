#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=versions.env
source "$HERE/versions.env"

SRC="${XIOS_SRC_DIR:-$HERE/src/xios}"
BUILD_ROOT="${XIOS_BUILD_ROOT:-$HERE/build}"
BUILD="${XIOS_BUILD_DIR:-$BUILD_ROOT/xios}"
PREFIX="${XIOS_PREFIX:-$HERE/install/xios}"
JOBS="${JOBS:-$(nproc)}"
ARCH=mini_env
NETCDF_MODE="${XIOS_NETCDF_LIB:-netcdf4_seq}"

for tool in mpicc mpif90 make nf-config nc-config; do
    command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
[[ -d "$SRC" ]] || { echo "XIOS source missing; run ./download_xios.sh first" >&2; exit 1; }

NC_INC="$(nf-config --includedir) $(nc-config --includedir)"
NC_INCDIR_FLAGS="-I$(nf-config --includedir) -I$(nc-config --includedir)"
NC_LIBDIR="$(nc-config --libdir)"
NC_LIB="$(nf-config --flibs) -lnetcdf"

mkdir -p "$BUILD_ROOT" "$PREFIX" "$HERE/logs"
mkdir -p "$SRC/arch"

cat > "$SRC/arch/arch-${ARCH}.env" <<ENV
export NETCDF_INC_DIR="$NC_INC"
export NETCDF_LIB_DIR="$NC_LIBDIR"
export HDF5_INC_DIR=""
export HDF5_LIB_DIR=""
export BOOST_INC_DIR=""
export BOOST_LIB_DIR=""
ENV

cat > "$SRC/arch/arch-${ARCH}.fcm" <<'FCM'
%CCOMPILER      mpicc
%FCOMPILER      mpif90
%LINKER         mpif90

%BASE_CFLAGS    -ansi -w
%PROD_CFLAGS    -O2
%DEV_CFLAGS     -g -O2
%DEBUG_CFLAGS   -g

%BASE_FFLAGS    -D__NONE__
%PROD_FFLAGS    -O2
%DEV_FFLAGS     -g -O2
%DEBUG_FFLAGS   -g

%BASE_INC       -D__NONE__
%BASE_LD        -lstdc++
%CPP            cpp
%FPP            cpp -P
%MAKE           make
FCM

cat > "$SRC/arch/arch-${ARCH}.path" <<PATH
NETCDF_INCDIR="$NC_INCDIR_FLAGS"
NETCDF_LIBDIR="-L$NC_LIBDIR"
NETCDF_LIB="$NC_LIB"
MPI_INCDIR=""
MPI_LIBDIR=""
MPI_LIB=""
HDF5_INCDIR=""
HDF5_LIBDIR=""
HDF5_LIB=""
BOOST_INCDIR=""
BOOST_LIBDIR=""
BOOST_LIB=""
OASIS_INCDIR=""
OASIS_LIBDIR=""
OASIS_LIB=""
PATH

rm -rf "$BUILD"
echo "Building XIOS in $BUILD"
(
  cd "$SRC"
  ./make_xios --full --prod --arch "$ARCH" --netcdf_lib "$NETCDF_MODE" \
      --build_path "$BUILD_ROOT" --build_dir "$(basename "$BUILD")" --job "$JOBS"
) > "$HERE/logs/xios-build.log" 2>&1

for item in bin lib inc; do
    [[ -e "$BUILD/$item" ]] || { echo "XIOS build did not produce $BUILD/$item" >&2; tail -60 "$HERE/logs/xios-build.log" >&2; exit 1; }
    cp -a "$BUILD/$item" "$PREFIX/"
done
echo "XIOS installed in $PREFIX"
