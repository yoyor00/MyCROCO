
# Contributing


This guide describes the conventions and workflow expected from all
contributors to the project.

Most rules follow common **modern Fortran** community practices that
CROCO is gradually moving toward, even though some parts of the codebase
still follow older conventions.

Because CROCO historically evolved from **Fortran 77 code**, parts of the
codebase may still contain legacy patterns. Contributors are encouraged
to follow the modern conventions described in this guide when adding new
code or refactoring existing routines, while preserving compatibility
with the current structure of the model.

The long-term goal is to progressively improve readability, modularity,
and maintainability without introducing disruptive changes to the
existing workflow.


------------------------------------------------------------------------

# Syntax and Style

Follow common **modern Fortran style conventions** and enforce them
using tools such as **fprettify** or **findent** when formatting code.

-   Indent with **2 to 4 spaces** (do not use tabs).
-   Use **free-form source format** (`.f90`, `.f95`, `.f03`, `.f08`).
-   Line length should ideally stay below **132 characters**.

------------------------------------------------------------------------

# Code Structure

Follow the spirit of **clarity and simplicity**.

-   Always use `implicit none`.
-   Prefer **modules** over global state.
-   Avoid duplicated logic; factor reusable routines into helper
    procedures.
-   Keep procedures **short and composable** rather than creating overly
    complex abstractions.
-   Check what is **already implemented** before writing new code.
-   Explicitly specify `intent(in)`, `intent(out)`, or `intent(inout)`
    for arguments.
-   Avoid hidden side effects.
-   Avoid unnecessary generalisation or configuration unless needed.

------------------------------------------------------------------------

# Project Structure and Components

Understanding where things live helps you put new code in the right
place.

    OCEAN/
      Core Fortran source code, initialisation, IO, atmospheric 
      forcing (BULK, ABL), sediment model using USGS...

    BENCH/
      Functionnal tests script in python.

    TEST_CASES/
      Input test_cases configurations.

    AGRIF/
      Nesting capabilities using 
      AGRIF (Adaptive Grid Refinement in Fortran)

    PISCES/
      BGC model PISCES

    MUSTANG/
      Sediment model MUSTANG

    OBSTRUCTION/
      Flow in presence of various obstructions

    XIOS/
      XIOS IO part

    SCRIPTS/
      Shell scripts to help run the code

    MPI_NOLAND/
      Tools to use MPI_NOLAND capabilities


Guidelines:

-   Place new code in **existing modules** if it fits logically.
-   Create new modules only when functionality clearly belongs in a
    separate component.
-   Avoid circular dependencies between modules.


------------------------------------------------------------------------

# Git Workflow

## Branches

-   **release-vx.x.x** --- tagged stable releases
-   **master** --- integration branch
-   **feature branches** --- one branch per issue or feature

Use descriptive names such as:

    feature_new_adv_scheme
    123_fix_mpi_buffer
    fix_bulkonline_computation
    dev_2026_stochastic
    

------------------------------------------------------------------------

# Development Workflow

1.  **Open an issue** describing the bug or feature.
2.  **Create a branch** from `master`.
3.  **Develop on your branch**, committing regularly.
4.  Periodically **rebase or merge `master`** to keep the branch up to
    date.
5.  Before submitting ensure the code is running trhough test suite and 
    update `CHANGELOG.md` if needed
6.  **Open a Merge Request / Pull Request** targeting `master`.
7.  **Squash commits** if necessary into clear logical commits.
8.  After review, the MR/PR is merged.
9.  Delete the branch once merged.

------------------------------------------------------------------------

# Merge Requests

-   Target branch is usually **master**.

-   Request at least **one reviewer**.

-   Ensure that:

    -   tests pass
    -   documentation is updated if necessary in dedicated repository

A good MR description should include:

-   purpose of the change
-   related issue
-   summary of implementation
-   potential side effects

------------------------------------------------------------------------

# Testing

Testing is essential for maintaining numerical reliability.

They are performed in BENCH directory see 
https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.bench.html
for more details.

A minimal test suite is perform by :

```
./bench-croco.py \
      --case @base,@nbq,@sediment_mustang,@other \
      --variants @debug --debug --report
```

------------------------------------------------------------------------

# General Philosophy

-   Prefer **clarity over cleverness**.
-   Keep routines **simple and well documented**.
-   Avoid premature optimization unless profiling demonstrates the need.
-   Write code that future developers (including yourself) can easily
    understand.