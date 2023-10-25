Benchmarking scripts
====================

The given directory aimed at providing some basic benchmarking scripts to
compare various parallelization mode on CROCO.

Concepts
--------

First you can give a look to the `config.jsonc` file. You will see the concepts
of :

- `cases` : The cases to run (BASIN, SOLITON...) and some case size bu changing
  the grid size (eg. BASIN-LARGE, BASIN-HUGE).
- `variants`: The parallelization mode (openmp, mpi, gpu) with specification on
  the ressources to use because CROCO needs to be rebuilt for every case du to
  some currently static definitions in the CROCO code on the use of those ressources.
- `hosts`: You will find some variables to adapt the running to the host you
  are running on.

Dependencies
------------

**Note**: It currently uses the `cmake` build and not the `jobcomp` one.

As it will build CROCO for you, you will need to have the basic dependencies
available (see main README.cmake.md to know what is required).

You might need some extended dependences so if you built the deps prefix for
CROCO with `create_prefix_with_deps.sh`, then just run:

```sh
# if you already have an activated venv
pip install -r requirement.txt
```

Otherwise you first needs to create a `virtualenv` :

```sh
# only the first time
python3 -m venv venv

# for every terminal used to run the bench script
source venv/bin/activate

# only once
pip install -r requirement.txt
```

Running
-------

The simplest way to run is :

```sh
# simplest way
./bench-croco.py

# use your own config file
./bench-croco.py -c config-mine.jsonc
```

It will run the `cases` and `variants` combination defined in `default_selection`
of the config file.

You can also choose by hand some of those :

```sh
# just a single
./bench-croco.py --case BASIN --variants sequential

# or many
./bench-croco.py --case BASIN,BASIN-LARGE --variants sequential,openmp-8,mpi-8,openacc-psyclone
```

Build & results
---------------

The script will build CROCO for each `case`/`variant` in the workdir (`./rundir`
by default).

The results are stored in `./results` (by default) and a graph will be produced
for each case.

By default the benchmarking script will run CROCO 4 times for each case to plot
the error bars. You can change this easily by using :

```sh
./bench-croco.py --runs 1
```

The script will keep the results of each run by `{hostname}-{date}` but you can
also ask the script to add a subdirectory with a given name if you are making different tries.

```sh
# while you are developping your optimizations (trash dir)
./bench-croco.py --title dev
# Produce clean final perf ready for paper
./bench-croco.py --title final-perf
```

Partial re-reun. If you make a new run on the same machine not running all the
cases the script will search the missing results in the previous result
directories in order to aggregate no the charts (with same `title`). This is
usefull in dev process to compare without re-running everything.

Configuration
-------------

The configuration file can be provided under `json` or `yaml` form as you
prefer. It can also be splitted in subdirectory mode by merging sevarl subfiles.

In the main `config.yaml` file, the inclusion is described by the `imports`
section.

Meta variants
-------------

You will notice in the `config.jsonc` file the definition of `meta_variants`
which allow you to define sets of variants which you might use often for your
run.

Then you can just use them easily:

```sh
./bench-croco.py --case BASIN --variants @cpu
```

Meta case
---------

If you want to simply run all available cases and all variants, you can just:

```sh
./bench-croco.py --case @all --variants @all --auto-skip
```

Variables
---------

You might notice in the configuration file the use of some variables to be
replaced at parsing time on the form of : `{tuning.gnu}`.

It avoids duplicating too many entries in the file.

You will mostly find them at two places : 

- On the cases to unpack king of template not to fully duplicate all cases for
  OpenMP and MPI.
- To inject some values from `host` into the `configure` cmake call.
- To inject the mpi extra options in `command_prefix`.

**Note**:

- For the template base system you can recurse the template var (`{thread}`, `tasks`)
  to build the final var name (as it is done in the default config).
- For the rest you cannot currently recurse even it can be trivial to patch the code to enable it
  (just that I didn't had time to validate it really works by doing it).
  For this case, you have acces to two base keys : `tuning` & `case`
  corresponding to the corresponding par of the config file you are running
  (host & current case).

Imports
-------

You will also notice at the head of the configuration file the `imports` section
which allow to split the config file in multiple sub files and merge them
at loading time.

This is usefull for exemple for the `hosts` definitions as it permits a user
to define its own hosts without requirement to patch the files tracked by
git and to keep its definition for himself until he finally decide the commit
it.

Author
------

- Sébastien Valat - INRIA / LJK - 2023
