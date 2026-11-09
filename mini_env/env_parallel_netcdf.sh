#!/bin/bash
# Source after mini_env/env.sh to prefer the parallel NetCDF stack.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=env.sh
source "$HERE/env.sh"
PREFIX="$HERE/install/parallel-netcdf"
export MINI_ENV_PARALLEL_NETCDF="$PREFIX"
export PATH="$PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CPPFLAGS="-I$PREFIX/include${CPPFLAGS:+ $CPPFLAGS}"
export LDFLAGS="-L$PREFIX/lib${LDFLAGS:+ $LDFLAGS}"
