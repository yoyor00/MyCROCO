Building CROCO with CMake
=========================

This page descrived how to build CROCO with CMake in place of the default
manual script.

**CAUTION :** It currently still lack some support.

Dependencies
------------

The CROCO **required** dependecies are :

 * CMake (Tested : 3.22.1) : https://cmake.org/
 * NetCDF-Fortran (Tested: 4.5.3) : https://downloads.unidata.ucar.edu/netcdf/

The CROCO **optional** dependecies are :

 * NVidia HPC SDK (Tested: 2023 / 23.7 ) : https://developer.nvidia.com/hpc-sdk-downloads
 * PSyClone (Tested : 2.3.1, branch: async-and-merge-master-ok) : https://psyclone.readthedocs.io/
 * psyclone-poseidon (Tested : master) : https://gitlab.inria.fr/svalat/croco-psyclone

Installing the dependencies
---------------------------

If you want to quickly build and install the dependencies localy in a subdirectory
you can simply use the provided script :

```sh
# by default it install in ./venv
./create_prefix_with_deps.sh [--nvhpc] [PREFIX_DIR]
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
../configure --with-case=CANYON
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
../configure --with-parallel=openmp
```

Building the GPU/OpenACC version
--------------------------------

Simply play with the `configure` options by :

 * Select the OpenACC mode : `--with-openacc={MODE}` which can be either :
    * `psyclone` to use the OpenACC auto-generated version with PSyClone.
    * `native` to use the hand made version.
 * Change the compiler for the NVHPC one.

```sh
# manual version
../configure --with-parallel=openacc-native FC=nvfortran

# psyclone version
../configure --with-parallel=openacc-psyclone FC=nvfortran
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

The `configure` script is just a simple wrapper arround CMake, is you look
closely you will see the first line printed by the script showing the exact
`cmake` command used.

You can of course directly call `cmake` as you want to play with more advanced
options.

Seeing tu actual CMake configuration
------------------------------------

You can get a basic summary on the setup by calling `cmake` without options :

```sh
cmake ..
```

If you can use the terminal cmake browser `ccmake` to get access to all the
defined variables and see their values:

```sh
# the curses version
ccmake ..
```
