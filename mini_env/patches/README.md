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
