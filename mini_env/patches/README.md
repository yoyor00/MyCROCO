# Local mini_env patches

This directory contains experimental overlays used by `mini_env` scripts.
They are not applied to the CROCO source tree.

## `XIOS/file_def_croco.xml_full_withcpp`

A copy of `XIOS/file_def_croco.xml_full_withcpp` for the BASIN setup.

BASIN uses `dt = 9600 s`.  The upstream XML activates output groups with
`output_freq="3h"` (10800 s), which is not an integer multiple of `dt`.
Unpatched XIOS can segfault in that configuration.

The patched copy sets both `file_group id="3h"` definitions to
`output_freq="9600s"`.  `mini_env/run_xios_runtime.sh` preprocesses this copy
with the run's `cppdefs.h` and writes the result to the run directory,
overwriting the XML generated from the upstream template.

Keep this file in sync with `XIOS/file_def_croco.xml_full_withcpp` when the
upstream template changes.

## `OCEAN/`

Experimental CROCO source overlays used by `run_checks.sh`.  The script
creates a symlink farm under `work/patched_source`, replaces the files below
with these copies, and points `jobcomp --src` at it.

### `OCEAN/grid.h`

Under `SPHERICAL`, upstream `grid.h` declares `xr`, `yr`, `xp`, `yp` as
allocatable but never allocates them, while `ana_grid.F` still writes them.
This causes a segfault in `ana_grid_lengths_`.  The patched copy uses static
`GLOBAL_2D_ARRAY` arrays instead.

### `OCEAN/croco_namelist_read.F90`

With `PSOURCE_NSRC=0`, the per-source arrays are zero-sized.  Writing the
`croco_psource_data` / `croco_psource_tracer` namelists then crashes some
gfortran versions.  The patched copy skips those writes when
`psource_Nsrc == 0`.

`run_checks.sh` also enables `ANA_SSFLUX` for the XIOS BASIN build so the
analytic surface salinity flux is used instead of a missing
`CROCO_FILES/croco_frc.nc`.
