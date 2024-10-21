Testing
=======

I wrote some unit tests to validate the PSyClone transformation rules applies.
You can simply run them by calling :

```sh
# jump in PSYCLONE directory
cd PSYCLONE
# call pytest
pytest
```

Do more
-------

The tests currently contains very basic checkings (at least one per transformation)
and probably need to be re-enfoced with more checkings.

Snippet extraction
------------------

The tests about `i,j,k` loop tranformations have a bit more advanced tests
base on snippets extracted from CROCO (currently from BASIN case only).

There is an automated way to extract the snippets while running the PSyClone
scripts.

You need to patch the code (need to make a cleaner way, eg. env var) to activate
the dump :

```python
# into scripts/lib/extensions/loops.py
ENABLE_SNIPPET_DUMPS=True
```

Then just make a build for GPU and you will get the snippets extracted to
easily extend the tests with more known loops if the code changed.

TODO : Testing snippets on GPU
------------------------------

In `tests/scripts/lib/extensions/test_running.py` you will find the necessary
to test the transformed snippets by running them on CPU.

There is really little work to extend those tests to also valide to GPU.

In practice the only blocker with the given code is the make a replacement
of `random_seed` and `random_numbers` which are GNU and not supported by
nvfortran. Then just stay to generate a third function with the
psychlone generated version of the loop.

Dumping all loops & kernels
---------------------------

You can enable the loop snippet dumping easily at any time by setting
environnement variable before calling `make`.

```sh
# will dump all loops inside ./OCEAN/loops-snippets/:
CROCO_PSYCLONE_SNIPPET_DUMPS=true make
```

Todo : Extend
-------------

We might want to extend this snipped based strategy to the full CROCO so we
can easily valide that all transformations are valid. Need to patch the other
transformations to make same.

Refactoring
-----------

The unit tests are currently written in a bit verbose way, maybe we can avoid
some code duplication by making to refactoring of the tests.
