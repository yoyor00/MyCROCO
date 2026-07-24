# WOOM user files

Per-scale source file overrides for WOOM-driven MS3DVAR compilation
(`TOOLS/ms3dvar_compile_scale.sh`, invoked by the `compile_filter`/`compile_lr`/
`compile_mr`/`compile_ms` tasks in `WOOM/tasks.cfg`).

## Usage

Drop a file in the subdirectory matching the scale (lowercase):

```
WOOM/user_files/
├── lr/
│   └── ms3dvar_lr_param.h   # overrides the normal LR param.h
├── mr/
├── ms/
└── filter/
```

Anything placed in `<scale>/` is copied into the build directory **last**,
after every other source (`OCEAN/`, `COMMON/`, and the scale's own files in
`ASSIM/MS3DVAR/<SCALE>/`), so it takes precedence over all of them.

An empty or missing subdirectory is a no-op — no override happens for that
scale.

## Why here, not somewhere else

This directory lives inside the WOOM workflow directory (alongside
`tasks.cfg`/`hosts.cfg`/`workflow.cfg`) so overrides ship and version
together with the workflow that needs them, rather than depending on a
machine-local path a user has to remember to set up separately. `tasks.cfg`
points `ms3dvar_compile_scale.sh` at it via `--user-files {{ workflow_dir }}/user_files`,
appended to each `compile_<scale>` task's commandline (WOOM's `{{ workflow_dir }}`
template variable resolves to this directory) — see the comment block at the
top of `tasks.cfg`.

This is independent from `jobcomp`'s own "local files" convenience (for
interactive/manual builds) and unrelated to `MS3DVAR_SRC` — these files are
not part of the versioned MS3DVAR sources.
