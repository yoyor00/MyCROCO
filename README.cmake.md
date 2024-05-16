Building CROCO with CMake
=========================

This page descrived how to build CROCO with CMake in place of the default
manual script.

**CAUTION :** It currently still lack some support (eg. AGRIF & XIOS).

Dependencies
------------

The CROCO **required** dependecies are :

 * CMake (Tested : 3.22.1) : https://cmake.org/
 * NetCDF-Fortran (Tested: 4.5.3) : https://downloads.unidata.ucar.edu/netcdf/

The CROCO **optional** dependecies are :

 * NVidia HPC SDK (Tested: 2023 / 23.7 ) : https://developer.nvidia.com/hpc-sdk-downloads
 * PSyClone (Tested : 2.3.1, branch: async-and-merge-master-ok of fork https://github.com/svalat/PSyclone/) : https://psyclone.readthedocs.io/

Installing the dependencies
---------------------------

If you want to quickly build and install the dependencies localy in a subdirectory
you can simply use the provided script :

```sh
# by default it install in ./venv
#
# Note: to have cmake installed you really need to force --cmake,
#       but you might oftenly prefer to use the system one in practice.
./create_prefix_with_deps.sh [--netcdf] [--nvhpc] [--psyclone] [--cmake] [--all] [PREFIX_DIR]
```

Then afterward while you will need in each terminal session to run:

```sh
# load the venv paths in your shell env
source venv/activate.sh
```

Installing most dependecies via the system
------------------------------------------

If you prefer to install the dependencies via your system, under `debian` flavor
you will want to do:

```sh
# minimum requirement
sudo apt install cmake libnetcdff-dev python3 python3-virtualenv
# if you follow the nvidia setup to get the repo for package
sudo apt install nvhpc-23-7
# You will still need the venv to get the patched psyclone
./create_prefix_with_deps.sh --psyclone
# source
source venv/bin/activate
```

Building
--------

In order to build CROCO you need to create a `build` directory at the place
you want (can be a subdir of the sources and whatever name you want).

```sh
mkdir build
cd build
```

Remark that each build directory you create is autonomous and can use
different `case` and `parallelism` mode.

Then you need to call the `configure` script which is just a simple wrapper
arround the `cmake` command to ease its use and provide an interface similar
to `autotools`.

To build the sequential version you can :

```sh
# simplest way
../configure

# if needs to say where to find NetCDF :
../configure --with-netcdf=$HOME/usr-netcdf/

# If all the required dependencies are in the same prefix
# you can simply (not to use --with-XXX for each one)
../configure --prefix=$HOME/usr-all/

# Default CASE is BASIN, but you can alter with :
../configure --with-case=BASIN

# Pass any variable to CMake directly bypassing the configure script.
../configure --with-netcdf=$HOME/usr-netcdf/ -DWITH_CASE=CANYON

# Enable some cppkey features via command line
../configure --with-case=RIP --with-keys=+NBQ,-BISCA
```

Finally you can simply build :

```sh
# build
make -j8

# point the cases config files so CROCO find them
ln -s ../TEST_CASES ./

# run
./croco
```

Building the OpenMP version
---------------------------

Simply play with the `configure` options :

```sh
../configure --with-optim=openmp --with-threads=8
```

Building the GPU/OpenACC version
--------------------------------

Simply play with the `configure` options by :

 * Select the OpenACC mode : `--with-optim={MODE}` which can be either :
    * `openacc-psycone` to use the OpenACC auto-generated version with PSyClone.
    * `openacc-native` to use the hand made version.
 * Change the compiler for the NVHPC one.

```sh
# manual version
../configure --with-optim=openacc-native FC=nvfortran

# psyclone version
../configure --with-optim=openacc-psyclone FC=nvfortran
```

Building the MPI version
------------------------

Simply play with `--with-optim` :

```sh
# MPI with 8 tasks
../configure --with-optim=mpi --with-splitting=2x4 FC=mpifort
```

Custom compile flags
--------------------

Depending on where you are running, you might want to use some specific
compile flags.

You can simply pass them to the script :

```sh
# force some optim
../configure FFLAGS=-march=native
# build
make
```

The default flags are computed in `cmake/croco_compiler_options.cmake`

More advanced configure options
-------------------------------

You can get all the `configure` option via

```sh
../configure --help
```

Handling CMake
--------------

The `configure` script is just a simple wrapper arround CMake, if you look
closely you will see the first line printed by the script showing the exact
`cmake` command used.

You can of course directly call `cmake` as you want to play with more advanced
options.

Seeing the actual CMake configuration
-------------------------------------

You can get a basic summary on the setup by calling `cmake` without options :

```sh
cmake ..
```

If you can use the terminal cmake browser `ccmake` to get access to all the
defined variables and see their values:

```sh
# the curses version
ccmake ..
# With the QT GUi
cmake-qui ..
```

If you used the `configure` wrapper script, you can also see the history of
calls to the script in :

```sh
# if you used ./configure
cat configure.log
```

Playing directly with cmake
---------------------------

Of course you can bypass the `configure` script if you already know well
`cmake`.

Here an example :

```sh
mkdir build
cd build
cmake .. -DNETCDF_PREFIX=$HOME/usr-netcdf -DWITH_OPTIM=openmp -DBUILD_TYPE=Debug
make
```

You can look in file `options` to get the mapping between `./configure` options
and `cmake` variables.

Also notice that the `configure` script display the exact `cmake` command call
it generates so you can also starts from there is you want to extend.
