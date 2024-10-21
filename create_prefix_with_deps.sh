#!/bin/bash

###########################################################
# set strict
set -u
set -e
set -x

############################################################
# default
INST_CROCO_SOURCE_DIR=$(realpath $(dirname $0))
INST_PREFIX=${PWD}/venv
INST_USE_NVHPC=no
INST_USE_NETCDF=no
INST_PSYCLONE=no
INST_USE_CMAKE=no

###########################################################
# selection versions
INST_NETCDF_C_VERSION=4.9.2
INST_NETCDF_FORT_VERSION=4.5.3
INST_NETCDF_CXX_VERSION=4.3.1
INST_HDF5_VERSION=1_14_2
INST_CMAKE_VERSION=3.27.7
# we use git version for psyclone
INST_PSYCLONE_GITHUB_USER=svalat
INST_PSYCLONE_BRANCH=async-and-merge-master-ok
# When patches will be submitted we might fallback to something close to :
#INST_PSYCLONE_GITHUB_USER=stfc
#INST_PSYCLONE_BRANCH=2.4.0

############################################################
# help
for arg in "$@"
do
	case "${arg}" in
		-h|--help)
			echo "Usage : $0 [--nvhpc] [--netcdf] [--psyclone] [--cmake] [--all] [PREFIX_DIR]" 1>&2
			exit 0
			;;
		--nvhpc)
			INST_USE_NVHPC=yes
			;;
		--netcdf)
			INST_USE_NETCDF=yes
			;;
		--psyclone)
			INST_PSYCLONE=yes
			;;
		--cmake)
			INST_USE_CMAKE=yes
			;;
		--all)
			INST_USE_NVHPC=yes
			INST_USE_NETCDF=yes
			INST_PSYCLONE=yes
			;;
		*)
			INST_PREFIX=$(realpath ${arg})
			;;
	esac
done

###########################################################
# install NetCDF
INST_SOURCES=${INST_PREFIX}/sources

###########################################################
# create a prefix
if [[ ! -e ${INST_PREFIX}/bin/activate ]]; then
	# create a prefix
	python3 -m venv ${INST_PREFIX}
	source ${INST_PREFIX}/bin/activate

	# install python deps
	pip install --upgrade pip
else
	# load prefix
	source ${INST_PREFIX}/bin/activate
fi

###########################################################
# download sources
mkdir -p ${INST_SOURCES}
cd ${INST_SOURCES}

###########################################################
# Clone psyclone
if [[ ${INST_PSYCLONE} == 'yes' && ! -f ${INST_PREFIX}/bin/psyclone ]]; then
	git clone https://github.com/${INST_PSYCLONE_GITHUB_USER}/PSyclone.git
	pushd PSyclone
	git checkout ${INST_PSYCLONE_BRANCH}
	pip install .
	pip install -r ${INST_CROCO_SOURCE_DIR}/PSYCLONE/requirements.txt
	popd
fi

###########################################################
# build HDF-5
if [[ ${INST_USE_NETCDF} == 'yes' && ! -f ${INST_PREFIX}/include/hdf5.h ]]; then
	wget --continue https://github.com/HDFGroup/hdf5/releases/download/hdf5-${INST_HDF5_VERSION}/hdf5-${INST_HDF5_VERSION}.tar.gz
	tar -xvf hdf5-${INST_HDF5_VERSION}.tar.gz
	pushd hdfsrc
	./configure --prefix=$INST_PREFIX
	make -j8 install
	popd
fi

###########################################################
# build netcdf-c
if [[ ${INST_USE_NETCDF} == 'yes' && ! -f ${INST_PREFIX}/include/netcdf.h ]]; then
	wget --continue https://downloads.unidata.ucar.edu/netcdf-c/${INST_NETCDF_C_VERSION}/netcdf-c-${INST_NETCDF_C_VERSION}.tar.gz
	tar -xvf netcdf-c-${INST_NETCDF_C_VERSION}.tar.gz
	pushd  netcdf-c-${INST_NETCDF_C_VERSION}
	./configure --prefix=$INST_PREFIX --disable-byterange CFLAGS=-I${INST_PREFIX}/include LDFLAGS=-L${INST_PREFIX}/lib
	make -j8 install
	popd
fi

###########################################################
# build netcdf-fortran
if [[ ${INST_USE_NETCDF} == 'yes' && ! -f ${INST_PREFIX}/include/netcdf.inc ]]; then
	wget --continue https://downloads.unidata.ucar.edu/netcdf-fortran/${INST_NETCDF_FORT_VERSION}/netcdf-fortran-${INST_NETCDF_FORT_VERSION}.tar.gz
	tar -xvf netcdf-fortran-${INST_NETCDF_FORT_VERSION}.tar.gz
	pushd netcdf-fortran-${INST_NETCDF_FORT_VERSION}
	./configure --prefix=$INST_PREFIX CFLAGS=-I$(nc-config --includedir) LDFLAGS=-L$(nc-config --libdir)
	make -j8 install
	popd
fi

###########################################################
# build netcdf-cxx
if [[ ${INST_USE_NETCDF} == 'yes' && ! -f ${INST_PREFIX}/include/ncFile.h ]]; then
	wget --continue https://downloads.unidata.ucar.edu/netcdf-cxx/${INST_NETCDF_CXX_VERSION}/netcdf-cxx4-${INST_NETCDF_CXX_VERSION}.tar.gz
	tar -xvf netcdf-cxx4-${INST_NETCDF_CXX_VERSION}.tar.gz
	pushd netcdf-cxx4-${INST_NETCDF_CXX_VERSION}
	./configure --prefix=$INST_PREFIX CFLAGS=-I$(nc-config --includedir) CXXFLAGS=-I$(nc-config --includedir) LDFLAGS=-L$(nc-config --libdir)
	make -j8 install
	popd
fi

###########################################################
# build cmake
if [[ ${INST_USE_CMAKE} == 'yes' && ! -f ${INST_PREFIX}/bin/cmake ]]; then
	wget --continue https://github.com/Kitware/CMake/releases/download/v${INST_CMAKE_VERSION}/cmake-${INST_CMAKE_VERSION}.tar.gz
	tar -xvf cmake-${INST_CMAKE_VERSION}.tar.gz
	pushd cmake-${INST_CMAKE_VERSION}
	./configure --prefix=$INST_PREFIX
	make -j8
	make install
	popd
fi

###########################################################
# build load script
cat > ${INST_PREFIX}/activate.sh <<EOF
#!/bin/bash

INST_PREFIX=${INST_PREFIX}

# venv
source \${INST_PREFIX}/bin/activate

EOF

###########################################################
# import nvidia SDK
if [[ ${INST_USE_NVHPC} == 'yes' ]]; then
	if [[ ! -d ${INST_PREFIX}/nvhpc_sdk ]]; then
		# install
		wget --continue https://developer.download.nvidia.com/hpc-sdk/23.7/nvhpc_2023_237_Linux_x86_64_cuda_12.2.tar.gz
		tar xpzf nvhpc_2023_237_Linux_x86_64_cuda_12.2.tar.gz
		nvhpc_2023_237_Linux_x86_64_cuda_12.2/install <<EOF

1
${INST_PREFIX}/nvhpc_sdk
EOF
	fi

	# add to load
	cat >> ${INST_PREFIX}/activate.sh <<EOF
# cuda
export NVARCH=\`uname -s\`_\`uname -m\`; export NVARCH
export NVCOMPILERS=\${INST_PREFIX}/nvhpc_sdk; export NVCOMPILERS
export MANPATH=\$MANPATH:\$NVCOMPILERS/\$NVARCH/23.7/compilers/man; export MANPATH
export PATH=\$NVCOMPILERS/\$NVARCH/23.7/compilers/bin:\$PATH; export PATH
#export PATH=\$NVCOMPILERS/\$NVARCH/23.7/comm_libs/mpi/bin:\$PATH
EOF
fi

###########################################################
# Finished
set +x
echo
echo
echo "================= DONE ==================="
echo "You can now source the env by using :"
echo ""
echo "------------------------------------------"
echo "source ${INST_PREFIX}/activate.sh"
echo "------------------------------------------"
