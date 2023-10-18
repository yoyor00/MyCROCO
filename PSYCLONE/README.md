PSyClone CROCO GPU / OpenACC
============================

This directory contains the PSyClone (https://github.com/stfc/PSyclone/) sripts to apply on CROCO to produce an OpenACC version of the code for
GPUs.

Installing psylcone
-------------------

To use this content, you will need to install psyclone plus some extra
python libs.

You can use the script provided in one shot :

```sh
# the first time
./create_prefix_with_deps.sh --psyclone ./venv
# in every new terminal used to build
source ./venv/bin/activate
```

**Otherwise** the requirement is to install :

- https://github.com/svalat/PSyclone/ branch `async-and-merge-master-ok`
- Python package `deepdiff` and `pytest` if you want to run the unit tests.

Installing nvfortran
--------------------

If you do not yet have `nvfortran` command on your system, you might also
use (caution, it takes a large space in the `venv` dir):

```sh
./create_prefix_with_deps.sh --nvhpc ./venv
source ./venv/bin/activate
```

**Otherwise** use the manual procedure :

- https://developer.nvidia.com/hpc-sdk

Building & running GPU version - jobcomp
----------------------------------------

If you want to build the psyclone GPU version with jobcomp you will need
to :

1. Create the `Run` directory as usual.
1. Edit the cppkeys in `cppdefs.h` to activate the key `OPENACC`
1. By default it triggers the `OPENACC_PSYCLONE` key at the end of `cppdefs.h`. If you want the manual version, disable this line.
1. Edit `jobcomp` the change :

```sh
#FC=gfortran
FC=nvfortran
```

Run `jobcomp` and enjoy.

Building & running GPU version - cmake
--------------------------------------

If you want to use the CMake build you can symply just :

```sh
mkdir build
cd build
../configure --with-parallel=openacc-psyclone FC=nvfortran
make
```

Running the unit tests
----------------------

You can simply run the unit test checking the internal behavor of the scripts
with :

```sh
# either (with cmake)
ctest
# or full manual
pytest ./PSYCLONE
```

More details
------------

If you want more details on how the port has been made, look into [./doc](./doc).
