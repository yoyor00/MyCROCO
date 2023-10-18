Using PSyClone
==============

This page recap shortly how to use [PSyClone](https://github.com/stfc/PSyclone/)
if someone whant to use it on CROCO for other stuff than what is done.

API
---

PSyClone can be used in two ways:

- either as a lib and calling the internal API functions to build our own program.
- either by using the `psyclone` command reusing some internals done for NEMO.

Calling psyclone
----------------

Calling psyclone with the command line just work like this (for what is required
for CROCO):

```sh
psyclone -api nemo -l output -s {SCRIPT} -opsy {OUTPUT_FORT_FILE} {INPUT_FORT_FILE}
```

The `{SCRIPT}` should just be a python file providing a `trans()` function
as an entry point. This function will reshape the recived intermediate
representation to `psyclone` can regenerate the code as output.

Dymmy trans
-----------

A dummy transformation doing nothing is symply :

```python
def tran(psy):
    pass
```

More
----

To look on how to transform the code, go in the documentation :
https://psyclone.readthedocs.io/
