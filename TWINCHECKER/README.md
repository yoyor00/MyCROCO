Croco - PsyClone - Twin Checker
===============================

What it is
----------

The twin-checker is a draft psyclone tool to compare on the fly the value
computed between two run of the same app made in parallel.

WARNING - POC
-------------

It is currently a POC so code need to be cleaned up.

Parts
-----

The tool comes in two parts :

1. An instrumentation script (`twin_psyclone_script.py`) which inject calls
   to the checking library at the point we want to check variables.
1. A library which is responsible at runtime to link the two twins and track
   the value compare along the run via a shared memory segment.
1. A script to launch the two twins linked by the library (`twin_checker_runner.sh`).

Usage in CROCO
--------------

The use in CROCO is currently only available via the CMake build system, you need
to :

```sh
./configure --enable-twin-checker
```

Make the both build you want to compare in there respective directories.

Then run the two programs via :

```sh
# jump in master
cd build-sequential-gcc

# launch both
../TWINCHECKER/twin_checker_runner.sh ./croco ./../build-sequential-nvfortran/croco
```

Configuration
-------------

You can configure it via `TWINCHECKER/config.jsonc`.

It list the files to instrument and define the options to enable / disable.

The basic case is :

```jsonc
{
    "defaults": {
        /* Profile all input variables */
        "profile_in_values": true,
        /* Profile all output variables (assigns) */
        "profile_out_values": true,
        /* Also profile integers, not just reals */
        "profile_integers": false,
        /* Profile the full arrays at the exit of a kernel (after it). Usefull for GPU where we cannot instrument inside */
        "profile_after_kernel": false,
        /* Profile the full arrays before a kernel. Usefull for GPU where we cannot instrument inside */
        "profile_before_kernel": false,
        /* Before and after kernels, instrument only the arrays, not scalars (usefull for GPU case) */
        "before_after_only_arrays": true,
        /* If want to profile the full arrays or only scalars (considering x(i,j,k) as scalar) */
        "profile_arrays": false,
        /* If there is a diff fix the salve value by the master one so we continue as it would be on master. */
        "fixable": true
    },
    "files": [
        "step2d.F",
        "zetabc.F",
        "u2dbc.F",
    ]
}
```

You can also use the `specific` section to define the parameters to apply specificaly
to some files by using groups.

```jsonc
{
    "specific" : {
        /* Define a group with a name we want */
        "steps_2d": {
            /* Files to override conf */
            "files": [
                "step2d.F",
                "rho_eos.F"
            ],
            /* Config values to override over default */
            "config": {
                "profile_integers": false
            }
        },
    }
}
```

Principle
---------

The twin checker is based on a shared memory segment which is used to sync the
two twins on every check. During a check, each one exchange its value with the
other and compare.

If there is a difference, then run is failure. In that case we can either :

 - Fail
 - Report the error and fix it taking again the master value and continue so see other errors.

Limitations
-----------

It is currently valid only when the operations are done in the same order on both
twins. it currently forbid comparing a sequential & parallel version.

Also caution, as it is called as a PSyClone pass it is currently relativeley
incompatible with using it for GPU. The only working case is applying twin-checker
on file which are in a second step GPU-ised by psyclone as it will inject again
the ACC directives.

Exemple of output
-----------------

```plain
---------------------- Twin Checker Error -------------------------
Not maching value in /home/svalat/Projects/minicroco/BENCH/rundir/BASIN-SMALL/sequential-nvfortran/croco :
 - master = -1.76676247346918605187783818e-06
 - slave  = -1.76676247346918584011960136e-06
 - err    =                  ^^^^^^^^^^^^    
 - diff   = 2.11758236813575084767080625e-22
 - local  = -1.76676247346918584011960136e-06
-------------------------------------------------------------------
-cff1 * COS(cff2 * yr(i,j))
-------------------------------------------------------------------
At .../minicroco/BENCH/rundir/BASIN-SMALL/sequential-nvfortran/OCEAN/prepared_sources/OCEAN/analytical.cpp.mpc.twin_check.F90:837
-------------------------------------------------------------------
```
