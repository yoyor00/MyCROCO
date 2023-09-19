#!/bin/bash

###########################################################
# set strict
set -u
set -e
set -x

############################################################
# default
BENCH_PREFIX=${PWD}/venv
BENCH_USE_NVHPC=no

############################################################
# help
for arg in "$@"
do
    case "${arg}" in
        -h|--help)
            echo "Usage : $0 [--nvhpc] [PREFIX_DIR]"
            exit 0
            ;;
        --nvhpc)
            BENCH_USE_NVHPC=yes
            ;;
        *)
            BENCH_PREFIX="${arg}"
            ;;
    esac
esac

###########################################################
# install NetCDF
BENCH_SOURCES=${BENCH_PREFIX}/sources

###########################################################
# create a prefix
if [[ ! -e ./venv/bin/activate ]]; then
    # create a prefix
    python3 -m venv ./venv
    source ./venv/bin/activate

    # install python deps
    pip install --upgrade pip
else
    # load prefix
    source ./venv/bin/activate
fi

###########################################################
# selection versions
NETCDF_C_VERSION=4.9.2
NETCDF_FORT_VERSION=4.5.3
NETCDF_CXX_VERSION=4.3.1
HDF5_VERSION=1_14_2

###########################################################
# download sources
mkdir -p ${BENCH_SOURCES}
cd ${BENCH_SOURCES}

###########################################################
# Clone psyclone
if [[ ! -f ${BENCH_PREFIX}/bin/psyclone ]]; then
    git clone https://github.com/svalat/PSyclone.git
    pushd PSyclone
    git checkout async-and-merge-master-ok
    pip install .
    popd
fi

###########################################################
# Clone poseidon
if [[ ! -f ${BENCH_PREFIX}/bin/poseidon ]]; then
    git clone git@gitlab.inria.fr:svalat/croco-psyclone.git
    pushd croco-psyclone/psyclone-poseidon
    pip install .
    popd
fi

###########################################################
# build HDF-5
if [[ ! -f ${BENCH_PREFIX}/include/hdf5.h ]]; then
    wget --continue https://github.com/HDFGroup/hdf5/releases/download/hdf5-${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz
    tar -xvf hdf5-${HDF5_VERSION}.tar.gz
    pushd hdfsrc
    ./configure --prefix=$BENCH_PREFIX
    make -j8 install
    popd
fi

###########################################################
# build netcdf-c
if [[ ! -f ${BENCH_PREFIX}/include/netcdf.h ]]; then
    wget --continue https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF_C_VERSION}/netcdf-c-${NETCDF_C_VERSION}.tar.gz
    tar -xvf netcdf-c-${NETCDF_C_VERSION}.tar.gz
    pushd  netcdf-c-${NETCDF_C_VERSION}
    ./configure --prefix=$BENCH_PREFIX --disable-byterange CFLAGS=-I${BENCH_PREFIX}/include LDFLAGS=-L${BENCH_PREFIX}/lib
    make -j8 install
    popd
fi

###########################################################
# build netcdf-fortran
if [[ ! -f ${BENCH_PREFIX}/include/netcdf.inc ]]; then
    wget --continue https://downloads.unidata.ucar.edu/netcdf-fortran/${NETCDF_FORT_VERSION}/netcdf-fortran-${NETCDF_FORT_VERSION}.tar.gz
    tar -xvf netcdf-fortran-${NETCDF_FORT_VERSION}.tar.gz
    pushd netcdf-fortran-${NETCDF_FORT_VERSION}
    ./configure --prefix=$BENCH_PREFIX CFLAGS=-I$(nc-config --includedir) LDFLAGS=-L$(nc-config --libdir)
    make -j8 install
    popd
fi

###########################################################
# build netcdf-cxx
if [[ ! -f ${BENCH_PREFIX}/include/ncFile.h ]]; then
    wget --continue https://downloads.unidata.ucar.edu/netcdf-cxx/${NETCDF_CXX_VERSION}/netcdf-cxx4-${NETCDF_CXX_VERSION}.tar.gz
    tar -xvf netcdf-cxx4-${NETCDF_CXX_VERSION}.tar.gz
    pushd netcdf-cxx4-${NETCDF_CXX_VERSION}
    ./configure --prefix=$BENCH_PREFIX CFLAGS=-I$(nc-config --includedir) CXXFLAGS=-I$(nc-config --includedir) LDFLAGS=-L$(nc-config --libdir)
    make -j8 install
    popd
fi

###########################################################
# build load script
cat > ${BENCH_PREFIX}/activate.sh <<EOF
#!/bin/bash

BENCH_PREFIX=${BENCH_PREFIX}

# venv
source \${BENCH_PREFIX}/bin/activate

EOF

###########################################################
# import nvidia SDK
if [[ ${BENCH_USE_NVHPC} == 'yes' ]]; then
    # install
    wget --continue https://developer.download.nvidia.com/hpc-sdk/23.7/nvhpc_2023_237_Linux_x86_64_cuda_12.2.tar.gz
    tar xpzf nvhpc_2023_237_Linux_x86_64_cuda_12.2.tar.gz
    nvhpc_2023_237_Linux_x86_64_cuda_12.2/install <<EOF

1
${BENCH_PREFIX}/nvhpc_sdk
EOF

    # add to load
    cat >> ${BENCH_PREFIX}/activate.sh <<EOF
# cuda
export NVARCH=`uname -s`_`uname -m`; export NVARCH
export NVCOMPILERS=${BENCH_PREFIX}/nvhpc_sdk; export NVCOMPILERS
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.7/compilers/man; export MANPATH
export PATH=$NVCOMPILERS/$NVARCH/23.7/compilers/bin:$PATH; export PATH
#export PATH=$NVCOMPILERS/$NVARCH/23.7/comm_libs/mpi/bin:$PATH
EOF
fi

###########################################################
# Finished
echo "================= DONE ==================="
echo "You can now source the env by using :"
echo ""
echo "------------------------------------------"
echo "source ${BENCH_PREFIX}/activate.sh"
echo "------------------------------------------"
