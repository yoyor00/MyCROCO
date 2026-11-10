# Experimental mini environment for XIOS / OASIS

This directory is only for local experimentation with the optional XIOS and
OASIS3-MCT parts of the CROCO build.  It is not a production environment and
its results are not scientifically validated.

Downloaded sources, build trees, installed libraries and Python virtualenvs
live under this directory but are ignored by git.  No source tarballs are
committed to the repository.

## Quick start

```bash
cd mini_env
./install_all.sh --all
source env.sh
./check_env.sh
```

The default Python/BENCH environment can be created separately:

```bash
./setup_python.sh
source env.sh
```

## What is installed

- XIOS under `mini_env/install/xios`
- OASIS3-MCT under `mini_env/install/oasis`
- downloaded sources under `mini_env/src` (ignored)
- build trees under `mini_env/build` (ignored)
- logs under `mini_env/logs` (ignored)

The default XIOS source is the official IPSL repository:

```text
https://gitlab.in2p3.fr/ipsl/projets/xios-projects/xios.git
```

Override `XIOS_REPO` / `XIOS_REF` to test another checkout or mirror:

```bash
XIOS_REPO=https://github.com/barriern/xios2-nico.git XIOS_REF=xios-2.5.x \
    ./download_xios.sh
```

The OASIS default is OASIS3-MCT 4.0 from `gitlab.com`.

## Parallel NetCDF experiment

For the XIOS runtime investigation, a local parallel HDF5/NetCDF stack can be
built under `mini_env/install/parallel-netcdf`:

```bash
./install_parallel_netcdf.sh
source env_parallel_netcdf.sh
```

It downloads HDF5, NetCDF-C and NetCDF-Fortran archives into the ignored
`mini_env/downloads/` directory and builds everything locally.

An XIOS 3 parallel build can then be produced without replacing the default
XIOS installation:

```bash
git clone --depth 1 --branch xios-3.0.x \
    https://gitlab.in2p3.fr/ipsl/projets/xios-projects/xios.git src/xios3
source env_parallel_netcdf.sh
XIOS_SRC_DIR="$PWD/src/xios3" \
XIOS_BUILD_DIR="$PWD/build/xios3-par" \
XIOS_PREFIX="$PWD/install/xios3-par" \
XIOS_NETCDF_LIB=netcdf4_par ./install_xios.sh
```

The current local XIOS runtime still fails inside the XIOS write path even
with the parallel NetCDF stack; the investigation is documented as
experimental.

## OASIS toy model

The OASIS toy models can be built with:

```bash
source env.sh
./install_toy.sh
```

They are installed under `mini_env/install/toy` as `toy_atm`, `toy_oce` and
`toy_wav`.  This only builds the executables; a full coupled run still needs
the matching `namcouple`, grid/forcing files and a launch configuration.

## Build checks

```bash
source env.sh
./run_checks.sh
```

Optional XIOS runtime attempt:

```bash
./run_checks.sh
./run_xios_runtime.sh work/xios
```

The XIOS runtime helper generates `iodef.xml` and tries to launch CROCO plus
`xios_server.exe`.  It is experimental and may fail with the current
XIOS/OpenMPI combination.

### Current local status

- XIOS build/link: PASS with make and SCons
- OASIS build/link: PASS with make and SCons
- XIOS runtime: PASS with a parallel-NetCDF XIOS build.  BASIN completes
  10 model steps and writes `Basin_Example_9600s_inst_*.nc` and
  `Basin_Example_5d_aver_*.nc` with data.
- The previous NaN/`MPI_ERR_TRUNCATE` failure was caused by the combination
  `SPONGE + NO_SPONGE_GRID`.  The XIOS interface requires `x_sponge`, so the
  build now uses `SPONGE` without `NO_SPONGE_GRID`; the local
  `croco_namelist.F90` overlay declares `x_sponge` for that case.
- The XML overlay uses the non-spherical BASIN coordinates `x_rho`/`y_rho`.
- XIOS 2.5 with serial NetCDF still has the `multiple_file` vertical-axis
  split.  Use a parallel-NetCDF XIOS build for runtime.
- XIOS server mode (`using_server=true`) still crashes in the local OpenMPI
  launch setup, so client-only mode is the only usable configuration here.
- The earlier `CScattererConnector::transfer` segfault is gone once the output
  frequency is compatible with the BASIN timestep (`9600s` instead of `3h`).
- `run_checks.sh` builds from a symlink source overlay in
  `work/patched_source` with the experimental `patches/OCEAN/` files.  The
  XIOS BASIN build also enables `ANA_SSFLUX`, avoiding a missing
  `CROCO_FILES/croco_frc.nc`.
- The runtime helper applies a local patched copy of
  `XIOS/file_def_croco.xml_full_withcpp` from `patches/XIOS/`.  It sets the
  active instant output groups to `9600s`, which is an integer multiple of
  the BASIN `dt`.  The original `XIOS/` and `TEST_CASES/` files are left
  untouched.  `XIOS_OUTPUT_FREQ` can still override the generated XML for
  experiments.
- The runtime helper defaults to `using_server=false` for the local checks.
- OASIS runtime: not implemented; build/link check only

### Working XIOS runtime

Build a parallel-NetCDF XIOS, then build CROCO against it:

```bash
source env_parallel_netcdf.sh
XIOS_SRC_DIR="$PWD/src/xios" \
XIOS_BUILD_DIR="$PWD/build/xios-par" \
XIOS_PREFIX="$PWD/install/xios-par" \
XIOS_NETCDF_LIB=netcdf4_par ./install_xios.sh

XIOS_ROOT_DIR="$PWD/install/xios-par" \
MINI_ENV_BACKENDS=make ./run_checks.sh
```

Run a short smoke test:

```bash
source env_parallel_netcdf.sh
XIOS_ROOT_DIR="$PWD/install/xios-par" XIOS_NTIMES=10 \
    ./run_xios_runtime.sh work/xios
```

`XIOS_NTIMES` limits the run to the requested number of steps; omit it for
the full BASIN configuration.

This builds a BASIN configuration with:

- XIOS enabled, linking against `install/xios`
- MPI + OASIS coupling enabled, linking against `install/oasis`

and records logs in `mini_env/logs`.  It is a build/link smoke test, not a
full coupled simulation.  A full XIOS/OASIS run needs external XML,
`namcouple`, restart/forcing files and a multi-executable launch setup.

## Environment

`env.sh` sets:

```text
XIOS_ROOT_DIR
CROCO_PRISM_ROOT_DIR
XIOS / CPL
PATH
LD_LIBRARY_PATH
CPL_NAM_DIR / CPL_FILES_DIR
```

It optionally prepends `mini_env/.venv/bin` if `setup_python.sh` has been run.

## Removal

Remove temporary work/logs/downloads but keep the installed libraries:

```bash
./cleanup.sh
# preview first
./cleanup.sh --dry-run
```

Remove everything generated by the experimental setup:

```bash
./cleanup.sh --all
```

Individual pieces can be removed with `--work`, `--logs`, `--env`, `--build`
or `--source`.
