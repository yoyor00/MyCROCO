# Bit-reproducibility of CROCO vs number of MPI ranks

Investigation of the `BENCH` test environment (BENCH's own `bench-croco.py`
build/run/check pipeline) on the `SOLITON` and `BASIN` test cases with
**gfortran**.

Environment
-----------
| item | value |
|---|---|
| compiler | GNU Fortran 15.2.0 (Ubuntu 15.2.0-16ubuntu1) |
| MPI | Open MPI 5.0.10 |
| netCDF-Fortran | 4.6.2 (`nf-config`) |
| CPU | Intel Xeon W-2295 (36 HW threads), `-march=native` |
| CROCO | `master` @ ff87b114 |
| bench flags (host `@default`/`gnu`) | `-O2 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy -march=native` |
| check semantics | `libs/check.py` compares every variable of `check_outputs` with the `sequential` run **bit for bit** (`numpy.any(ref != actual)`, `allclose` only reported in addition) |

How to reproduce
----------------
```sh
cd BENCH
python3 -m venv venv && ./venv/bin/pip install -r <(deps: numpy matplotlib tqdm json5 pyyaml jsonschema netcdf4 termcolor xarray pandas f90nml)
./venv/bin/python ./bench-croco.py \
    --case SOLITON,BASIN --variants sequential,mpi-2,mpi-4,mpi-8 \
    -r --report -N --auto-skip --host local-16cpu
```
(`local-16cpu` is a local, untracked host entry in `config.d/hosts/` that only
raises the advertised core count so that `mpi-16`/`mpi-32` are not filtered
out. `@default` (8 cores) does not change any compiler flag.)

Verdict
-------

### BENCH `--report` summary (stock case definitions, gfortran)

```
CASE-VARIANT        Build    Run      Check
SOLITON-sequential   True    True      -            <- reference
SOLITON-mpi-2        True    True      True         <- bit-identical
SOLITON-mpi-4        True    True      True         <- bit-identical
BASIN-sequential     True    True      -            <- reference
BASIN-mpi-2          True    True      True         <- bit-identical
BASIN-mpi-4          True    True      False        <- 1-ULP differences
BASIN-mpi-8          True    True      False        <- 1-ULP differences
```

* **SOLITON is bit-reproducible** for every decomposition tested
  (`2x1`=2 ranks, `2x2`=4, `2x4`=8, `4x4`=16; also `1x2`).
* **BASIN is *not* bit-reproducible** as soon as the domain is split in **XI**:
  * `1x2` (2 ranks, ETA only)  -> bit-identical to sequential ✔
  * `1x4` (4 ranks, ETA only)  -> bit-identical to sequential ✔
  * `2x1` (2 ranks, XI)       -> differs ✘
  * `2x2` (4), `2x4` (8), `4x4` (16) -> differ ✘

Notes / limits of the stock BENCH configuration:
* `SOLITON` only declares splittings for 1, 2 and 4 ranks
  (`config.d/cases/soliton.jsonc`), so `mpi-8`/`mpi-16` are silently dropped.
  The 8/16-rank results above come from a local copy of the case with the same
  physics and an extended splitting table.
* `BASIN` declares `16: "1x16"`, but that decomposition **cannot even run**:
  CROCO aborts on the first history write with `nf_put_vara` error `-57`
  (invalid coordinates) — a case/decomposition-table problem, unrelated to
  reproducibility.

Magnitude of the BASIN differences (`basin_his.nc`, 3240 steps, 37 records)
--------------------------------------------------------------------------
`sequential` vs `mpi-4` (identical numbers for `mpi-8`, `2x1`, `4x4` are
within the same 1-ULP class):

| variable | # values differing | max abs diff | max rel diff | first record |
|---|---|---|---|---|
| `zeta` | 47 / 119288 (0.04 %) | 1.86e-9 | 9.7e-7 | 1 |
| `ubar` | 10 / 117364 (0.009 %) | 7.3e-12 | 2.0e-7 | 12 |
| `vbar` | 15 / 116994 (0.013 %) | 5.8e-11 | 2.0e-7 | 4 |
| `u` | 5315 / 1173640 (0.45 %) | 3.7e-9 | 7.1e-5 | 1 |
| `v` | 7937 / 1169940 (0.68 %) | 3.7e-9 | 1.3e-4 | 1 |
| `temp` | 3 / 1192880 | 2.4e-7 | 8.5e-8 | 6 |

All the differences are **1 ULP of the stored float32** (e.g.
`zeta = 1.2010913e-3` vs `1.2010915e-3`). The grid (`h`, `pm`, `pn`, `f`,
`x_rho`, `y_rho`, `Cs_*`) and `time` are bit-identical. The runs do not
diverge, the error only grows very slowly (round-off level).

Two important properties established separately:

1. **Runs are deterministic**: re-running the same binary with the same number
   of ranks reproduces `basin_his.nc` byte-for-byte (sha256 unchanged).
2. **The result depends on NP_XI only, not on the number of ranks**: `2x1`,
   `2x2` and `2x4` produce *bit-identical* output to each other; `4x4` (NP_XI=4)
   gives another (equally close) result. So it is a property of the compile-time
   sub-domain width, not of the parallel scheduling.

Where does the difference come from?
------------------------------------
Localised with CROCO's own CVTK instrumentation (BENCH `--cvtk`, which dumps
every `check_tab2d/3d` checkpoint in double precision and stops at the first
mismatch), 8-step runs, `sequential` as writer / `mpi-4` as reader:

* `h initialisation #1` : identical
* `u`, `v` at the first `set_vbc` : identical
* `temp` at the first `set_vbc` (i.e. the **analytical initial temperature,
  before any dynamics**) : differs by 1 ULP, only in the columns sitting on
  the internal XI boundary (and their ghost cells), all `j`, `k = 1..8`.
* Then everything downstream inherits the 1-ULP seed.

Since the grid arrays (`h`, `pm`, `Cs_r`) are identical and there is **no
`MPI_ALLREDUCE`/`MPI_REDUCE` anywhere in the time-stepping path** (only
`MPI_ALLGATHER` for diagnostics such as `zonavg`/`ana_jet`), the difference is
not a reduction-order effect and not an error in the halo exchange values: it
is the same expression evaluated with a different instruction sequence.

Compiler bisection (same case, same decomposition, only the flag changes)
------------------------------------------------------------------------
| compile flags | sequential vs mpi-4 (2x2) |
|---|---|
| `-O2 ... -march=native` (bench default) | **differs** |
| `-O2 ... -march=native -ffp-contract=off` | **differs** (FMA contraction is not the cause) |
| `-O2 ... -fno-tree-vectorize` | **bit-identical** |
| `-O2 ... -fno-tree-vectorize -ffp-contract=off` | **bit-identical** |
| `-O0` (no optimisation at all) | **bit-identical** |

=> The non-reproducibility is produced by **GCC auto-vectorisation at -O2**.
`param.h` contains compile-time sub-domain sizes (`NP_XI`, `NP_ETA`,
`size_XI`, `padd_X`, `N2d`, ...), so the number of sub-domains in XI changes
the loop trip counts of the generated code; the vectoriser then emits a
different sequence of (contracted/vector) floating-point operations for
`NP_XI=1` than for `NP_XI>1`, which shows up as 1-ULP differences. This
explains all the observations: same NP_XI -> same code shape -> identical
results; `-O0` -> same arithmetic -> identical results; SOLITON unaffected
(its kernels happen to compile to the same sequence for all tested splits).

### Exact mechanism: vectorised `exp`/`tanh` -> glibc **libmvec**, which is
### not bit-identical to the scalar libm

CROCO is full of `exp()`/`tanh()` in vectorisable loops. With
`-O2 -march=native` gfortran does not just vectorise the arithmetic, it
replaces the scalar libm calls by the **vector math** entry points of glibc
libmvec (`_ZGVdN4v_exp`, `_ZGVeN8v_tanh`, ...), which are different
implementations and are *not* bit-identical to scalar `exp`/`tanh`:

```
$ nm -u Compile/ana_initial.o | grep ZGV
   U _ZGVdN4v_cos   U _ZGVdN4v_exp   U _ZGVdN4v_tanh   U _ZGVeN8v_exp ...
```

Minimal standalone proof (20-line program, `y(i)=2*exp(x/0.8)*(0.6-0.4*tanh(x/0.8))`):

| build | result |
|---|---|
| `-O2 -march=native` (vectorised, libmvec) vs `-O2 -march=native -fno-tree-vectorize` (scalar libm) | **295 / 600 values differ, 1-2 ULP, max 7.1e-15** |

And this is exactly what happens in BASIN. The analytic initial temperature
loop (`ana_initial.F`:
`t(i,j,k,1,itemp)=cff2*exp(z_r(i,j,k)/800)*(0.6-0.4*tanh(z_r(i,j,k)/800))`)
is vectorised in one build and not in the other, simply because the loop
length changes with the decomposition:

```
# sequential build (IstrR..IendR = 62 columns)
ana_initial_.f:944:19: optimized: loop vectorized using 32 byte vectors

# mpi-4 build (2x2 -> 30 columns per rank): no vectorisation report for that loop
```

Full causal chain:

1. `NP_XI` is a compile-time constant (`param.h`) -> sub-domain width 60
   (sequential) vs 30 (`NP_XI=2`) vs 15 (`NP_XI=4`);
2. GCC's vectoriser therefore takes a different decision on the
   analytic-initialisation loops;
3. vectorised loops call libmvec `exp`/`tanh`, ~1-2 ULP away from scalar libm;
4. a 1-ULP seed is injected in the initial state at the columns processed by
   the vector path;
5. it propagates without amplifying through the 3240 steps and shows up as
   1 float32 ULP in `basin_his.nc`;
6. BENCH's strict `check` (any `!=` fails) therefore reports `False`.

No reduction-order problem, no halo-exchange bug, no error in the numerics:
the same binary is bit-reproducible run after run, and only the *number of
sub-domains in XI* - i.e. the *shape of the compiled code* - changes the last
bits.

Does the flag "fix" it? (verification + cost)
--------------------------------------------
`-fno-tree-vectorize` (nothing else) restores bit-reproducibility for **all**
BASIN decompositions that used to fail (BENCH `--report`, gfortran):

| case-variant | split | default `-O2 -march=native` | `-O2 -fno-tree-vectorize` |
|---|---|---|---|
| BASIN-mpi-2 | 1x2 | True | True |
| BASIN-mpi-4 | 2x2 | **False** | True |
| BASIN-mpi-8 | 2x4 | **False** | True |
| BASIN-mpi-2 (XL) | 2x1 | **False** | True |
| BASIN-mpi-16 (XL) | 4x4 | **False** | True |

* `-ffp-contract=off` is neither necessary nor sufficient: alone it changes
  nothing, together with `-fno-tree-vectorize` it changes nothing either.
* `-O0` also works but is obviously unusable for a benchmark.

Cost of the work-around (same machine, same run, BASIN 3240 steps):

| variant | `-O2` (vectorised) | `-fno-tree-vectorize` | cost |
|---|---|---|---|
| BASIN-sequential | 40.2 s | 58.9 s | **+46 %** |
| BASIN-mpi-2 | 30.3 s | 57.4 s | **+89 %** |
| BASIN-mpi-4 | 21.0 s | 30.5 s | **+45 %** |
| BASIN-mpi-8 | 16.8 s | 27.8 s | **+66 %** |

So it is a valid flag set for an *exactness* build (reference generation,
debug, CI validation) but **not** a good default for a performance benchmark.
Cheaper alternatives, in order of preference:

1. relax BENCH's `check` so that "same up to a few ULP" passes (the code
   already computes `numpy.allclose`; only the strict `!=` test fails), and
   keep the strict test as an opt-in;
2. apply `-fno-tree-vectorize` (or `!GCC$ novector`) only to the
   initialisation kernels that inject the seed (`OCEAN/ana_initial.F`) — cheap
   to validate, was not tested here;
3. keep the strict test but compare against a reference compiled with the
   *same* `NP_XI` — only proves "same decomposition reproduces itself".

Chosen setup (no source-code change): BENCH host `bitrepro`
----------------------------------------------------------
The two flags are bundled in a host entry (`config.d/hosts/local_16cpu.jsonc`)
so they can be used without touching anything else:

```sh
./bench-croco.py --case SOLITON,BASIN --variants sequential,mpi-2,mpi-4 \
                 --host bitrepro -r --report
# or permanently: add the two flags to hosts.@default.tuning.gnu.standard
# in BENCH/config.jsonc
```

```json
"standard": "-O2 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy -fno-tree-vectorize -ffp-contract=off",
"extra":    "-march=native"
```

What each flag does:

* `-fno-tree-vectorize` — do not turn loops into SIMD (AVX) code; **this is the
  one that fixes it**, because it stops GCC from replacing the scalar libm
  `exp`/`tanh` by glibc libmvec's vector versions (1-2 ULP apart) in some
  sub-domain widths only.
* `-ffp-contract=off` — never fuse `a*b+c` into one FMA instruction (one
  rounding instead of two); a safety belt so no multiply-add can be contracted
  differently depending on the code shape. No effect on its own here, harmless.

Verified with `--host bitrepro` (BENCH `--report`):

```
SOLITON-sequential  True  |  SOLITON-mpi-2 (2x1) True  |  SOLITON-mpi-4 (2x2) True
BASIN-sequential    True  |  BASIN-mpi-2  (1x2) True  |  BASIN-mpi-4 (2x2)  True  |  BASIN-mpi-8 (2x4) True
BASIN-XL-mpi-2 (2x1) True |  BASIN-XL-mpi-16 (4x4) True
```

Status of the 2D compile guards (are they "code modifications"?)
----------------------------------------------------------------
For a 3D build (`SOLVE3D`, e.g. BASIN) the guards preprocess to **the same
executable code as HEAD** — verified by diffing the `cpp -P` output of
`HEAD:OCEAN/<file>` against the patched file with the BASIN `cppdefs.h`:

| file | difference in the preprocessed 3D source |
|---|---|
| `croco_namelist_init.F90` | none |
| `croco_namelist_read.F90` | `use param, ONLY: stdout, NT` split into two `use` lines (same symbols) |
| `ana_initial.F` | only the added comment lines |

Still, they *are* 3 lines of diff in `OCEAN/`, and there is **no compiler flag
that can replace them**: the 2D build fails on `implicit none` semantic errors
(`z_r` undeclared, `NT` undeclared, namelist group not declared). So the options
are: keep the (numerically inert) guards to be able to validate SOLITON, or
revert them and restrict the bit-reproducibility validation to 3D cases.

Recommendations
---------------
* For CI "is my parallel branch correct?" runs, the strict `check` step of
  BENCH is too strict as soon as `NP_XI > 1` in a baroclinic case: either
  accept a tolerance (`check.py` already computes `allclose`, but still fails
  on the strict test), or
* build the reference and the compared variants with
  `-fno-tree-vectorize` (cheap: this is enough to get bit-identical results
  across decompositions), or compare against a reference produced with the
  *same* NP_XI, or
* keep vectorisation but avoid the libmvec/scalar-libm mismatch, e.g.
  `-fno-tree-vectorize` only for the affected files, or
  `-ffp-contract=off -fno-builtin-exp -fno-builtin-tanh` (to be validated),
  or use `!GCC$ novector` on the analytic-initialisation loops;
* for real bit-reproducibility (typically what you want for an exact-restart /
  debugging build), add `-fno-tree-vectorize` (or `-ffp-contract=off
  -oss`/vendor equivalents) to the "debug" tuning in the BENCH host files.
* Independently: `BASIN`'s `16 -> 1x16` splitting entry should be fixed or
  removed (it aborts with a netCDF error).

Work-arounds that had to be applied to run SOLITON at all
--------------------------------------------------------
On this HEAD, every **2D** configuration (`SOLVE3D` undefined) that defines
`ANA_INITIAL` fails to compile; SOLITON could not be built at all. Three
minimal, behaviour-preserving guards were added (`git diff`, uncommitted):

* `OCEAN/ana_initial.F` – the `ANA_JET` initialisation (`ana_jet_initial_tile`
  / `ana_jet_tile`) is baroclinic (`z_r`, `Hz`, `u`, `v`, `t`) but was compiled
  unconditionally: wrapped in `#if defined ANA_INITIAL && defined SOLVE3D`
  (both the subroutines and the call in `ana_initial_tile`).
* `OCEAN/croco_namelist_init.F90` – `use param, ONLY: NT` guarded with
  `#ifdef TRACERS` (3 places); `NT` only exists for 3D builds.
* `OCEAN/croco_namelist_read.F90` – `use param, ONLY: ... NT` split/guarded,
  and the `croco_auxiliary_history_fields` namelist block guarded with
  `#if defined SOLVE3D || defined RIP` (matching its declaration).

No change was made to any numerical kernel, so the reproducibility findings
are those of the unmodified code.

Helper scripts written for this investigation (untracked)
---------------------------------------------------------
* `BENCH/bitcheck.py` – per-variable sha256 / diff report between netCDF files.
* `BENCH/cvtkdiff.py` – reads CROCO `CVTK_DEBUG_WRITE` dumps and compares a
  sequential dump with the dumps of an MPI run (interior vs ghost cells).
* `BENCH/config.d/cases/local_*.jsonc` – extended splitting tables and short
  CVTK cases.
* `BENCH/config.d/hosts/local_16cpu.jsonc` – local hosts (`local-16cpu`,
  `local-nofma`, `local-o0`, `local-novec`, `local-novec-fma`).
