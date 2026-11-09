#!/bin/bash
# Experimental parallel NetCDF/HDF5 stack for XIOS.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOADS="$HERE/downloads"
SRC="$HERE/src"
BUILD="$HERE/build"
PREFIX="$HERE/install/parallel-netcdf"
JOBS="${JOBS:-$(nproc)}"

HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
NETCDF_C_VERSION="${NETCDF_C_VERSION:-4.9.3}"
NETCDF_FORTRAN_VERSION="${NETCDF_FORTRAN_VERSION:-4.6.2}"

mkdir -p "$DOWNLOADS" "$SRC" "$BUILD" "$PREFIX" "$HERE/logs"

fetch() {
    local name="$1" url="$2" archive="$DOWNLOADS/$1" srcdir="$SRC/$3"
    if [[ -d "$srcdir" ]]; then
        echo "$3 already extracted"
        return
    fi
    if [[ ! -f "$archive" ]]; then
        echo "Downloading $url"
        curl -L --fail --retry 3 -o "$archive" "$url"
    fi
    tar -xzf "$archive" -C "$SRC"
}

fetch "hdf5-${HDF5_VERSION}.tar.gz" \
      "https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz" \
      "hdf5-${HDF5_VERSION}"
fetch "netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
      "https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NETCDF_C_VERSION}.tar.gz" \
      "netcdf-c-${NETCDF_C_VERSION}"
fetch "netcdf-fortran-${NETCDF_FORTRAN_VERSION}.tar.gz" \
      "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NETCDF_FORTRAN_VERSION}.tar.gz" \
      "netcdf-fortran-${NETCDF_FORTRAN_VERSION}"

echo "Building parallel HDF5"
rm -rf "$BUILD/hdf5"
mkdir -p "$BUILD/hdf5"
(
  cd "$BUILD/hdf5"
  CC=mpicc FC=mpif90 "$SRC/hdf5-${HDF5_VERSION}/configure" \
      --prefix="$PREFIX" \
      --enable-parallel \
      --enable-shared \
      --disable-static \
      --disable-fortran \
      --disable-tools \
      --disable-tests
  make -j "$JOBS"
  make install
) > "$HERE/logs/parallel-hdf5-build.log" 2>&1

echo "Building parallel NetCDF-C"
rm -rf "$BUILD/netcdf-c"
mkdir -p "$BUILD/netcdf-c"
(
  cd "$BUILD/netcdf-c"
  CC=mpicc CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib" \
      LIBS="-lhdf5" "$SRC/netcdf-c-${NETCDF_C_VERSION}/configure" \
      --prefix="$PREFIX" \
      --enable-netcdf-4 \
      --enable-shared \
      --disable-static \
      --disable-dap \
      --disable-byterange \
      --disable-utilities \
      --disable-tests
  make -j "$JOBS"
  make install
) > "$HERE/logs/parallel-netcdf-c-build.log" 2>&1

echo "Building parallel NetCDF-Fortran"
rm -rf "$BUILD/netcdf-fortran"
mkdir -p "$BUILD/netcdf-fortran"
(
  cd "$BUILD/netcdf-fortran"
  CC=mpicc FC=mpif90 CPPFLAGS="-I$PREFIX/include" \
      LDFLAGS="-L$PREFIX/lib" LIBS="-lnetcdf" \
      "$SRC/netcdf-fortran-${NETCDF_FORTRAN_VERSION}/configure" \
      --prefix="$PREFIX" \
      --enable-shared \
      --disable-static
  make -j "$JOBS"
  make install
) > "$HERE/logs/parallel-netcdf-fortran-build.log" 2>&1

echo "Parallel NetCDF installed in $PREFIX"
