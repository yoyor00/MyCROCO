#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=versions.env
source "$HERE/versions.env"

SRC="$HERE/src/oasis3-mct"
BUILD="$HERE/build/oasis"
PREFIX="$HERE/install/oasis"
JOBS="${JOBS:-$(nproc)}"

for tool in mpicc mpif90 make nf-config nc-config; do
    command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
[[ -d "$SRC" ]] || { echo "OASIS source missing; run ./download_oasis.sh first" >&2; exit 1; }

NC_INC="$(nf-config --includedir)"
NC_LIBDIR="$(nc-config --libdir)"
NC_LIBS="$(nf-config --flibs) $(nc-config --libs)"

rm -rf "$BUILD"
mkdir -p "$BUILD" "$PREFIX/lib" "$PREFIX/build/lib" "$HERE/logs"

cat > "$SRC/util/make_dir/make.inc" <<MAKEINC
CHAN            = MPI1
COUPLE          = $SRC
ARCHDIR         = $BUILD

MPIDIR          = /usr
MPIBIN          = /usr/bin
MPI_INCLUDE     =
MPILIB          =

NETCDF_INCLUDE1 = $NC_INC
NETCDF_INCLUDE2 =
NETCDF_LIBRARY  = $NC_LIBS

MAKE            = make
F90             = mpif90
F               = \$(F90)
f90             = \$(F90)
f               = \$(F90)
CC              = mpicc
LD              = mpif90
AR              = ar
ARFLAGS         = -ruv

CPPDEF          = -Duse_comm_\$(CHAN) -D__VERBOSE -DTREAT_OVERLAY
F90FLAGS_1      = -ffree-line-length-none -fallow-argument-mismatch
f90FLAGS_1      = \$(F90FLAGS_1)
FFLAGS_1        = \$(F90FLAGS_1)
fFLAGS_1        = \$(F90FLAGS_1)
CCFLAGS_1       =
LDFLAGS         =

FLIBS           = \$(NETCDF_LIBRARY)
BINDIR          = \$(ARCHDIR)/bin
LIBBUILD        = \$(ARCHDIR)/build/lib
INCPSMILE       = -I\$(LIBBUILD)/psmile.\$(CHAN) -I\$(LIBBUILD)/scrip -I\$(LIBBUILD)/mct

F90FLAGS        = \$(F90FLAGS_1) \$(INCPSMILE) \$(CPPDEF) -I\$(NETCDF_INCLUDE1)
f90FLAGS        = \$(f90FLAGS_1) \$(INCPSMILE) \$(CPPDEF) -I\$(NETCDF_INCLUDE1)
FFLAGS          = \$(FFLAGS_1) \$(INCPSMILE) \$(CPPDEF) -I\$(NETCDF_INCLUDE1)
fFLAGS          = \$(fFLAGS_1) \$(INCPSMILE) \$(CPPDEF) -I\$(NETCDF_INCLUDE1)
CCFLAGS         = \$(CCFLAGS_1) \$(INCPSMILE) \$(CPPDEF) -I\$(NETCDF_INCLUDE1)
MAKEINC

echo "Building OASIS3-MCT in $BUILD"
(
  cd "$SRC/util/make_dir"
  export FCFLAGS="-g -O2 -fallow-argument-mismatch"
  make -f TopMakefileOasis3
) > "$HERE/logs/oasis-build.log" 2>&1

for lib in "$BUILD"/lib/*.a; do
    [[ -e "$lib" ]] || { echo "OASIS build produced no libraries" >&2; tail -80 "$HERE/logs/oasis-build.log" >&2; exit 1; }
    cp -a "$lib" "$PREFIX/lib/"
done
for subdir in psmile.MPI1 mct scrip; do
    [[ -d "$BUILD/build/lib/$subdir" ]] && cp -a "$BUILD/build/lib/$subdir" "$PREFIX/build/lib/"
done
echo "OASIS installed in $PREFIX"
