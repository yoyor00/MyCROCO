Compiler wrapper
================

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

Configuration
-------------

The compiler wrapper script will automatically search and load the
`psyclone.rules.json` file to get the list of transformation scripts to apply
on each source file.
