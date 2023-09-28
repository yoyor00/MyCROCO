PSyClone scripts for CROCO
==========================

This directory contains the PSyClone (https://github.com/stfc/PSyclone/) sripts to apply
on CROCO to produce an OpenACC version of the code for GPUs.

Calling psyclone
----------------

Calling psyclone is just done via the command line as :

```sh
psyclone -api nemo -l output -s {SCRIPT} -opsy {OUTPUT_FORT_FILE} {INPUT_FORT_FILE}
```

Wrapper
-------

In order to ease the injection of the filtering in the compiling instructions
this directory provide the `psyclone.compiler.wrapper.py` script which does
three things :

 * Read the `psyclone.rules.json` to know which script to apply on the given source file.
 * Call `psyclone` with the given script.
 * For debugging and ease looking on what the `psyclone` script has injected, produce
   a file with a dummy script to get just the source reformated in the `psyclone` way.
   So we can just `diff` with the produced file.

The way to use it is just prepending the script in fron of the compiling command :

```sh
./psyclone.compiler.wrapper.py gfortran {ALL_FLAGS_OPTIONS} {SOURCE_FILE}
```

Note
----

It helps on the path to handle the fact that the F77 files which goes through
psyclone are renamed as F90 files and not the others which is hidden from the
top compiling system.
