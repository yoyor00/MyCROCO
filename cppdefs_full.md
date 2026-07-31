# CROCO CPP Keys Reference

Companion to `croco_nml_full.md`. 
Advanced options (rarely changed) are in `cppdefs_dev.h`.
Keys automatically derived from others should not be set here.

Indented keys are only meaningful when their parent key is defined.

---


## Parallelization and I/O Server

Only one of `OPENMP` / `MPI` / `OPENACC` may be defined at a time. 

- `OPENMP` — Open-MP shared-memory parallelism
- `MPI` — MPI distributed-memory parallelism
  - `PARALLEL_FILES` — Parallel I/O: each MPI task writes its own file
  - `NC4PAR` — NetCDF-4 parallel I/O (requires HDF5/NC4)
  - `MPI_NOLAND` — Skip land-only tiles (load balancing)
  - `MPI_TIME` — Time MPI communications
  - `XIOS` — Use the external XIOS I/O server (needs MPI)
- `OPENACC` — OpenACC GPU 

Output : 

- `FILLVAL` — Write fill value in masked land points

---

## Calendar

- `USE_CALENDAR` — Absolute date/time mode

---

## Logging

- `LOGFILE` — Redirect STDOUT to a dedicated `croco.log` file

---

## Non-Hydrostatic Option

- `NBQ` — Non-Boussinesq (non-hydrostatic) solver
  - `NBQ_PRECISE` — More precise NBQ time-stepping (slower)
  - `NOT_NBQ_AM4` — Disable AM4 time filter in the NBQ solver
  - `ACOUSTIC_FORCING` — Activate acoustic source forcing in NBQ mode (ACOUSTIC test case)
- `CROCO_QH` — Quasi-hydrostatic option (advanced)

---

## Nesting (AGRIF)

- `AGRIF` — Adaptive Mesh Refinement / nesting (1-way by default)
- `AGRIF_2WAY` — 2-way nesting: update parent solution from child

---

## Coupling

`OA_COUPLING` and `OW_COUPLING` both require OASIS-MCT and MPI.

- `OA_COUPLING` — Ocean-Atmosphere coupling via OASIS-MCT _(auto-activates `READ_PATM` and `OBC_PATM`)_
  - `OA_GRID_UV` — Interpolate coupling fields on U/V grids
- `OW_COUPLING` — Ocean-Wave coupling via OASIS-MCT (WAVEWATCH III)
  - `OW_COUPLING_FULL` — Full wave-current coupling (otherwise one-way)
  - `WAVE_SMFLUX` — Apply wave-induced surface momentum flux

---

## Open Boundary Conditions — Direction Flags

Define the sides of the domain that have open boundaries.

- `OBC_EAST` — Open eastern boundary
- `OBC_WEST` — Open western boundary
- `OBC_NORTH` — Open northern boundary
- `OBC_SOUTH` — Open southern boundary

---

## Open Boundary Conditions — Scheme Choices

Choose one scheme per variable class (M2, M3, T).

- `OBC_M2CHARACT` — Characteristic methods for barotropic velocities
- `OBC_M2ORLANSKI` — Radiative OBC for barotropic velocities
- `OBC_M2SPECIFIED` — Specified (Dirichlet) OBC for barotropic velocities
- `OBC_VOLCONS` — Enforce volume conservation at open boundaries (use with `OBC_M2ORLANSKI`)
- `OBC_SPECIFIED_WEST` — Specified (Dirichlet) OBC on western boundary
- `OBC_M3ORLANSKI` — Radiative OBC for baroclinic velocities
- `OBC_M3SPECIFIED` — Specified OBC for baroclinic velocities
- `OBC_TORLANSKI` — Radiative OBC for tracers
- `OBC_TSPECIFIED` — Specified OBC for tracers
- `OBC_TUPWIND` — Upwind OBC for tracers

---

## Tides

`TIDES` activates the tidal forcing machinery; sub-keys select which tidal fields are processed.

- `TIDES` — Activate tidal forcing at open boundaries
  - `SSH_TIDES` — Read and apply tidal sea-level data; if neither `ZCLIMATOLOGY` nor
    `Z_FRC_BRY` is already set, auto-activates `ZCLIMATOLOGY` + `ANA_SSH` as the default
  - `UV_TIDES` — Read and apply tidal current data; if neither `M2CLIMATOLOGY` nor
    `M2_FRC_BRY` is already set, auto-activates `M2CLIMATOLOGY` + `ANA_M2CLIMA` as the default
  - `POT_TIDES` — Include tidal potential forcing
  - `TIDES_MAS` — Harmonic composition to build tide elevation from SHOM (Simon method);
    requires `USE_CALENDAR`; auto-activates `MASKING`
  - `TIDERAMP` — Ramp tidal forcing over 1 day at start
  - `OBC_REDUCED_PHYSICS` — Compute tidal velcocity from tidal elevation in case of 
    tidal current is not available (undef `UV_TIDES`)

---

## Applications

- `BIOLOGY` — Biogeochemical modeling (choose model in the Biology section below)
- `STATIONS` — High-frequency output at fixed station points
- `PASSIVE_TRACER` — Add one passive (non-reactive) tracer
- `SUBSTANCE`- Other way to add tracers (required for MUSTANG)
- `MUSTANG` — MUSTANG sediment dynamics model
- `SEDIMENT` — USGS sediment transport model
- `BBL` — Bottom Boundary Layer parameterization
- `OBSTRUCTION` — Sub-grid flow obstruction (e.g. seagrass)

---

## Stochastic / Ensemble

- `STOGEN` — Stochastic forcing generator
- `ENSEMBLE` — Ensemble mode

---

## Grid Configuration

- `CURVGRID` — Curvilinear coordinate transformation
- `SPHERICAL` — Longitude/latitude grid positioning
- `MASKING` — Land masking
- `WET_DRY` — Wetting-and-drying scheme
- `NEW_S_COORD` — New (Song & Haidvogel 1994) vertical S-coord
- `EW_PERIODIC` — East-West periodic boundary conditions
- `NS_PERIODIC` — North-South periodic boundary conditions

---

## Model Dynamics

- `SOLVE3D` — Solve 3D primitive equations (baroclinic)
- `UV_COR` — Coriolis terms in momentum equations
- `UV_ADV` — Advection terms in momentum equations
- `NO_FRCFILE` — Run without any external forcing files
- `NO_TEMPERATURE` — Suppress temperature tracer
- `NO_TRACER` — Suppress all tracers (barotropic only)
- `CONST_TRACERS` — Freeze tracers (hold T/S constant; no tracer advection or diffusion)
- `ZONAL_NUDGING` — Zonal mean flow nudging (used with JET)

Barotropic time filter (choose one; `M2FILTER_POWER` is the default):

- `M2FILTER_POWER` — Power-law barotropic time filter _(default)_
- `M2FILTER_COSINE` — Cosine-shaped barotropic time filter
- `M2FILTER_FLAT` — Flat (top-hat) barotropic time filter
- `M2FILTER_NONE` — Disable the smoothed time-averaging weights and switch the barotropic
  mode to a generalized, tunable AB3-AM4 forward-backward time-stepping scheme instead
  (weighted by the `m2filter_alpha` namelist parameter)

Pressure gradient formulation (advanced; default is the density-Jacobian formulation with cubic
polynomial fit, Shchepetkin & McWilliams 2003):

- `PGF_BASIC_JACOBIAN` — Cheaper standard Jacobian formulation (for smooth topography)
  - `WJ_GRADP` — Weighted-Jacobian formulation weight (e.g. `# define WJ_GRADP 0.125`); 
    requires `PGF_BASIC_JACOBIAN` (has no effect otherwise)
- `PGF_FLAT_BOTTOM` — Cheaper pressure gradient for flat-bottom cases (drops z-grid x/y 
  gradient terms; only used without `PGF_BASIC_JACOBIAN`)

---

## Equation of State

- `SALINITY` — Salinity as an active tracer
- `NONLIN_EOS` — Nonlinear equation of state
  - `SPLIT_EOS` — Split EOS into adiabatic + compressible parts to reduce pressure-gradient errors; requires `NONLIN_EOS` 
- `RESET_RHO0` — Reset background density `rho0` to the initial volume-averaged sigma-t (+1000)
  each run, to minimize Boussinesq errors; auto-enabled by default with linear EOS (see
  `NO_RESET_RHO0`)
- `NO_RESET_RHO0` — Disable the automatic `RESET_RHO0` behavior when using linear EOS
- `GRAVITY` — Override default gravity constant (9.81 m/s²); set as `# define GRAVITY <value>` in param.h file

---

## Lateral Forcing (Climatology & Open Boundary)

`CLIMATOLOGY` activates processing of 3D climatological files used for interior nudging.
`FRC_BRY` activates open-boundary forcing files.

- `CLIMATOLOGY` — Process 2D/3D climatological data for nudging and open boundary forcing
  - `ZCLIMATOLOGY` — Process sea-level climatology
  - `M2CLIMATOLOGY` — Process barotropic velocity climatology
  - `M3CLIMATOLOGY` — Process baroclinic velocity climatology
  - `TCLIMATOLOGY` — Process tracer (T, S, ...) climatology
  - `ZNUDGING` — Nudging layer for sea level
  - `M2NUDGING` — Nudging layer for barotropic velocities
  - `M3NUDGING` — Nudging layer for baroclinic velocities
  - `TNUDGING` — Nudging layer for tracers
  - `ROBUST_DIAG` — Strong tracer nudging over the whole domain (diagnostic / relaxation simulations)
- `FRC_BRY` — Process 1D/2D boundary forcing files
  - `Z_FRC_BRY` — Open-boundary forcing for sea level
  - `M2_FRC_BRY` — Open-boundary forcing for barotropic velocities
  - `M3_FRC_BRY` — Open-boundary forcing for baroclinic velocities
  - `T_FRC_BRY` — Open-boundary forcing for tracers

---

## Surface Forcing

`BULK_FLUX`: compute air-sea fluxes from atmospheric variables (default: COARE3p0 with gustiness).
Without `BULK_FLUX`: use prescribed flux fields + `QCORRECTION`.

- `BULK_FLUX` — Bulk air-sea flux formulation
  - `BULK_ECUMEV0` — ECUME_v0 parameterization _(mutually exclusive with other BULK_ schemes)_
  - `BULK_ECUMEV6` — ECUME_v6 parameterization _(mutually exclusive with other BULK_ schemes)_
  - `BULK_WASP` — WASP parameterization _(mutually exclusive with other BULK_ schemes)_
  - `BULK_GUSTINESS` — Add gustiness correction to wind speed
  - `BULK_LW` — Long-wave radiation feedback from model SST
  - `SST_SKIN` — Prognostic skin SST (Zeng & Beljaars 2005)
  - `ANA_DIURNAL_SW` — Analytical diurnal cycle of shortwave radiation (use only if forcing has no diurnal cycle)
  - `ONLINE` — Read native atm files and interpolate online to the CROCO grid (default: cubic)
    - `AROME` — AROME high-resolution atmospheric model input
    - `ERA_ECMWF` — ERA / ECMWF reanalysis input
    - `BULK_MONTH_1DIGIT` — Use a single-digit month in online bulk forcing file names
  - `READ_PATM` — Read surface atmospheric pressure
    - `OBC_PATM` — Apply atm pressure gradient at open boundaries
  - `ABL1D` — 1D ABL model (column model coupled to ocean); auto-activates `BULK_FLUX`
    (mandatory whenever `ABL1D` is set — setting `ABL1D` alone is sufficient)
    - `ANA_ABL_LSDATA` — Analytical large-scale ABL data
    - `ANA_ABL_VGRID` — Analytical ABL vertical grid
    - `STRESS_AT_RHO_POINTS` — Compute surface stress at rho points _(auto-derived)_
    - `ABL_NUDGING` — Nudge ABL towards large-scale data _(auto-derived)_
    - `ABL_NUDGING_DYN` — Nudge dynamics (wind) _(auto-derived)_
    - `ABL_NUDGING_TRA` — Nudge thermodynamics (temperature) _(auto-derived)_
    - `ABL_DYN_RESTORE_EQ` — Restore ABL dynamics towards equilibrium
    - `ABL_NO_OCEAN_FEEDBACK` — Suppress ocean SST feedback to the ABL (one-way ABL > ocean coupling, used in test case KILPATRICK)
  - `QCORRECTION` — Linear bulk SST correction of heat flux _(without BULK\_FLUX only)_
  - `SFLX_CORR` — Freshwater flux correction towards model SSS _(without BULK\_FLUX only)_
  - `SFLX_CORR_COEF` — Prescribed correction coefficient _(without BULK\_FLUX only)_
- `SFLUX_CFB` — Current feedback on wind stress; 
  if `BULK_FLUX` additionally activates `CFB_STRESS`/`CFB_WIND_TRA`
- `SEA_ICE_NOFLUX` — Zero flux under sea ice
- `RAIN_FLUX` — Add rain-induced E-P (evaporation minus precipitation) flux into the mass budget (requires `SALINITY`)

---

## Lateral Momentum Advection

Choose exactly one horizontal advection scheme for momentum. Default (none defined): 3rd-order upstream (UP3).

- `UV_HADV_UP3` — 3rd-order upstream biased advection (default)
- `UV_HADV_UP5` — 5th-order upstream biased advection
- `UV_HADV_C2` — 2nd-order centred (use with explicit mixing)
- `UV_HADV_C4` — 4th-order centred (use with explicit mixing)
- `UV_HADV_C6` — 6th-order centred (use with explicit mixing)
- `UV_HADV_WENO5` — 5th-order WENO quasi-monotone advection
- `NO_M2_HADV_UP3` — Disable automatic UP3 advection for 2D (barotropic-only) runs (used in SOLITON test case)

---

## Lateral Momentum Mixing

- `UV_VIS2` — Laplacian horizontal viscosity
  - `UV_VIS_SMAGO` — Smagorinsky parameterization of viscosity (2D)
- `UV_VIS4` — Bilaplacian horizontal viscosity
- `UV_VIS_SMAGO_3D` — 3D Smagorinsky parameterization
- `UV_MIX_GEO` — Mix momentum on geopotential (z) surfaces
- `UV_MIX_S` — Mix momentum on iso-sigma surfaces

---

## Vertical Momentum Advection

- `UV_VADV_SPLINES` — Splines vertical advection for momentum _(default)_
- `UV_VADV_WENO5` — WENO5 vertical advection for momentum
- `UV_VADV_C2` — 2nd-order centered vertical advection for momentum

---

## Vertical Velocity Advection (W)

Only relevant with NBQ or wave-resolving simulations. 
Set automatically by `cppdefs_dev.h` when `UV_HADV_WENO5` is defined.

- `W_HADV_WENO5` — 5th-order WENOZ horizontal advection for W _(default)_
- `W_HADV_UP3` — 3rd-order upwind horizontal advection for W
- `W_HADV_UP5` — 5th-order upwind horizontal advection for W
- `W_HADV_C2` — 2nd-order centered horizontal advection for W
- `W_HADV_C4` — 4th-order centered horizontal advection for W
- `W_HADV_C6` — 6th-order centered horizontal advection for W
- `W_VADV_WENO5` — 5th-order WENOZ vertical advection for W _(default)_
- `W_VADV_SPLINES` — Splines vertical advection for W
- `W_VADV_C2` — 2nd-order centered vertical advection for W

---

## Lateral Tracer Advection

Choose exactly one horizontal advection scheme for tracers.

- `TS_HADV_UP3` — 3rd-order upstream biased advection
- `TS_HADV_RSUP3` — Rotated split 3rd-order upstream biased (recommended)
- `TS_HADV_RSUP5` — Rotated split 5th-order upstream biased, reduced dispersion/diffusion
- `TS_HADV_UP5` — 5th-order upstream biased advection
- `TS_HADV_C4` — 4th-order centred advection
- `TS_HADV_C6` — 6th-order centred advection
- `TS_HADV_WENO5` — 5th-order WENOZ quasi-monotone for all tracers
- `BIO_HADV_WENO5` — 5th-order WENOZ for passive/biology/sediment tracers only _(auto-set when PASSIVE\_TRACER, BIOLOGY, SEDIMENT or MUSTANG is defined)_

---

## Lateral Tracer Mixing

- `TS_DIF2` — Laplacian horizontal diffusion of tracers
- `TS_DIF4` — Bilaplacian horizontal diffusion of tracers
- `TS_MIX_ISO` — Mix along isopycnal (isoneutral) surfaces
- `TS_MIX_GEO` — Mix along geopotential (z) surfaces
- `TS_MIX_S` — Mix along iso-sigma surfaces
- `TS_MIX_IMP` — Stabilizing correction for rotated diffusion (use with `TS_MIX_ISO` or `TS_MIX_GEO`)

---

## Vertical Tracer Advection

- `TS_VADV_SPLINES` — Splines vertical advection for tracers _(default)_
- `TS_VADV_WENO5` — WENO5 vertical advection for tracers
- `TS_VADV_C2` — 2nd-order centered vertical advection for tracers

---

## Semi-Implicit Vertical Advection

- `VADV_ADAPT_IMP` — Adaptive semi-implicit vertical advection
  - `VADV_ADAPT_PRED` — Apply the adaptive semi-implicit scheme to both the predictor and 
    corrector steps (default: corrector step only)

---

## Sponge Layers

Enhance viscosity and diffusivity near open boundaries to damp spurious reflections.

- `SPONGE` — Sponge layers near lateral open boundaries;
  - `NO_SPONGE_GRID` — Disable reading sponge coefficients from the grid file (use analytical, used for INNERSHELF)
  - `SPONGE_GRID` — Read sponge width/value from the grid file _(auto-derived: on by default unless `NO_SPONGE_GRID` is set)_
  - `SPONGE_DIF2` — Enhanced tracer diffusivity in the sponge layer _(auto-derived, on by default)_
  - `SPONGE_VIS2` — Enhanced momentum viscosity in the sponge layer _(auto-derived, on by default)_
  - `SPONGE_SED` — Enhanced sediment diffusivity in the sponge layer _(auto-derived when SEDIMENT is defined)_

---

## Bottom Stress

- `LIMIT_BSTRESS` — Limit bottom stress to prevent negative transport _(auto-derived: on by
  default, unless `NO_LIMIT_BSTRESS` or `BSTRESS_FAST` is set — mutually exclusive with
  `BSTRESS_FAST`)_
  - `NO_LIMIT_BSTRESS` — Disable the bottom stress limiter (overrides `LIMIT_BSTRESS` default, used for INNERSHELF)
- `BSTRESS_FAST` — Compute bottom stress at 3D fast (baroclinic) time steps; also auto-activates
  `M3FAST` (and thus `SOLVE3D`, `M2FILTER_NONE`) and disables `LIMIT_BSTRESS`

---

## Vertical Mixing

Choose exactly one of `LMD_MIXING` or `GLS_MIXING` (or neither for analytical).

- `BODYFORCE` — Apply surface and bottom stresses as body forces instead of boundary conditions
- `ANA_VMIX` — Analytical viscosity/diffusivity (no TKE scheme)
- `LMD_MIXING` — KPP: Large-McWilliams-Doney turbulence closure
  - `LMD_SKPP` — KPP surface boundary layer mixing
  - `LMD_BKPP` — KPP bottom boundary layer mixing
  - `LMD_RIMIX` — Shear-instability interior mixing
  - `LMD_CONVEC` — Convective adjustment interior mixing
  - `LMD_NONLOCAL` — Nonlocal (counter-gradient) transport in SKPP
  - `LMD_DDMIX` — Double-diffusion interior mixing
  - `LMD_LANGMUIR` — Langmuir-cell parameterization
- `GLS_MIXING` — Generic Length Scale turbulence scheme
  - `GLS_MIXING_3D` — Full 3D GLS (vs. columnar 1D version); 
    auto-activates `GLS_MIXING` (and
    `UV_VIS2`/`VIS_COEF_3D`, disables `TS_DIF2`/`DIF_COEF_3D`) 
  - `GLS_KOMEGA` — K-omega closure _(choose one of the two GLS closures below; K-epsilon is the default if neither is defined)_
  - `GLS_KEPSILON` — K-epsilon closure _(default GLS closure)_
  - `GLS_COASTAL_ROUGHNESS` — Enhanced bottom roughness near coastal boundaries

  Stability function for the GLS closure (choose one; `CANUTO_A` is the default):

  - `CANUTO_A` — Canuto et al. (2001) stability function, version A _(default)_
  - `CANUTO_B` — Canuto et al. (2001) stability function, version B
  - `GibLau_78` — Gibson & Launder (1978) stability function
  - `MelYam_82` — Mellor & Yamada (1982) stability function
  - `KanCla_94` — Kantha & Clayson (1994) stability function
  - `Luyten_96` — Luyten et al. (1996) stability function
  - `Cheng_02` — Cheng et al. (2002) stability function

---

## Wave-Current Interactions

Wave forcing source (choose one; also used with `BBL` or `MUSTANG` independently of `MRL_WCI`):

- `WAVE_OFFLINE` — Read pre-computed wave fields from a file _(default when `MRL_WCI` is active
  without another source; also force-undefines `WAVE_ROLLER` in that auto-default case)_
- `ANA_WWAVE` — Analytical (constant) wave forcing (Hs, Tp, Dir)
- `WKB_WWAVE` — Internal WKB spectral wave model
  - `WKB_OBC_NORTH` — WKB open boundary – north side
  - `WKB_OBC_SOUTH` — WKB open boundary – south side
  - `WKB_OBC_WEST` — WKB open boundary – west side
  - `WKB_OBC_EAST` — WKB open boundary – east side
  - `WAVE_ROLLER` — Wave roller model
  - `WAVE_FRICTION` — Bottom friction
  - `ANA_BRY_WKB` — Analytical wave boundary forcing
  - `WKB_SHORELIKE_BRY` — Shorelike (zero-energy) inner boundary condition for WKB
  - `WKB_ADD_DIFF` — Extra WKB diffusion term _(auto-derived: always defined when `WKB_WWAVE`
    is active)_
  - `WKB_ADD_DIFFRACTION` — Wave diffraction term _(auto-derived: always defined when
    `WKB_WWAVE` is active)_
  - `WKB_KZ_FILTER` — Vertical-diffusivity filter for WKB _(force-undefined when `WKB_WWAVE`
    and `MRL_CEW` are both active)_
  - `WKB_TIME_FILTER` — Time filter for WKB _(force-undefined when `WKB_WWAVE` and `MRL_CEW`
    are both active)_
- `WAVE_IO` — Enable wave-field NetCDF I/O _(auto-derived: defined when any of `WKB_WWAVE`,
  `OW_COUPLING`, `ANA_WWAVE`, or (`WAVE_OFFLINE` and `MRL_WCI`) is active)_

Wave breaking parameterization (used with `MRL_WCI` or `WKB_WWAVE`; choose one):

- `WAVE_BREAK_CT93` — Thornton & Guza (1983) wave breaking _(default)_
- `WAVE_BREAK_TG86` — Church & Thornton (1993) wave breaking
- `WAVE_BREAK_TG86A` — Alternate Church-Thornton formulation
- `WAVE_BREAK_R93` — Roelvink (1993) wave breaking formulation
- `WAVE_BREAK_BJ78` — Battjes & Janssen (1978) wave breaking formulation

Wave-current interaction model:

- `MRL_WCI` — Wave-Current Interaction model (Mellor-type)
  - `STOKES_DRIFT` — Stokes drift _(auto-derived: always defined when `MRL_WCI` is active)_
  - `MRL_CEW` — Wave-driven longshore current energetics
  - `WAVE_STREAMING` — Bottom wave streaming (wave-induced current)
  - `WAVE_RAMP` — Ramp wave forcing at start
  - `WAVE_DRY` — Suppress wave forcing in cells shallower than `D_wavedry` (requires `WET_DRY`)

---

## Wave Maker

Wave-resolving internal wave generation (ideal cases, e.g. SANDBAR, SWASH, RIP).

- `WAVE_MAKER` — Wave-maker boundary source
  - `WAVE_MAKER_JONSWAP` — JONSWAP spectrum input _(auto-defines `WAVE_MAKER_SPECTRUM`; also the
    default chosen if `WAVE_MAKER_SPECTRUM` is set directly without picking JONSWAP or GAUSSIAN)_
  - `WAVE_MAKER_GAUSSIAN` — Gaussian spectrum input _(auto-defines `WAVE_MAKER_SPECTRUM`)_
  - `WAVE_MAKER_SPECTRUM` — Multi-frequency spectrum input _(auto-derived from
    `WAVE_MAKER_JONSWAP` or `WAVE_MAKER_GAUSSIAN`; can also be set directly, in which case
    `WAVE_MAKER_JONSWAP` is force-defaulted in; if none of these three is set, `WAVE_MAKER`
    falls back to `STOKES_WAVES`)_
  - `WAVE_MAKER_BICHROMATIC` — Bichromatic (two-frequency) wave input (e.g. SWASH GLOBEX B2/B3)
  - `WAVE_MAKER_DSPREAD` — Directional spreading
    - `WAVE_MAKER_DSPREAD_PER` — Corrects wave direction for N/S periodicity _(auto-derived:
      defined when `WAVE_MAKER_DSPREAD` and `NS_PERIODIC` are both active)_
  - `STOKES_WAVES` — Monochromatic Stokes-wave maker _(auto-derived: defined when `WAVE_MAKER`
    is active but none of `WAVE_MAKER_JONSWAP`, `WAVE_MAKER_GAUSSIAN` or `WAVE_MAKER_SPECTRUM`
    is chosen)_
- `WAVE_MAKER_INTERNAL` — Internal (immersed) wave maker

---

## Bottom Forcing

- `ANA_BSFLUX` — Analytical bottom salinity flux (typically 0)
- `ANA_BTFLUX` — Analytical bottom temperature flux (typically 0)

---

## Analytical Forcing

Use these keys when running test cases without external NetCDF files. Each `ANA_` key replaces a file read with an inline Fortran function.

- `ANA_GRID` — Analytical grid (no grid NetCDF file)
- `ANA_INITIAL` — Analytical initial conditions
- `ANA_BRY` — Analytical open boundary conditions
- `ANA_SSH` — Analytical sea level (replaces ZCLIMATOLOGY)
- `ANA_M2CLIMA` — Analytical barotropic velocity climatology
- `ANA_M3CLIMA` — Analytical baroclinic velocity climatology
- `ANA_TCLIMA` — Analytical tracer climatology
- `ANA_SMFLUX` — Analytical surface momentum (wind) flux
- `ANA_SRFLUX` — Analytical shortwave radiation flux
- `ANA_STFLUX` — Analytical surface heat flux
- `ANA_SSFLUX` — Analytical surface salinity flux
- `ANA_SST` — Analytical SST (for QCORRECTION)
- `ANA_TIDES` — Analytical tidal forcing (RIP / TIDAL_FLAT)
- `ANA_JET` — Analytical jet initial condition (JET case)
- `ANA_MORPHODYN` — Analytical morphodynamic forcing
- `ANA_PSOURCE` — Analytical vertical profiles for point sources

---

## Point Sources — Rivers

- `PSOURCE` — River / point-source mass and tracer input
- `PSOURCE_NCFILE` — Read time-varying river transports from NetCDF
  - `PSOURCE_NCFILE_TS` — Read time-varying river concentrations from NetCDF
- `PSOURCE_MASS` — Inject river as vertical volume flux at rho points (instead of flux through
  u,v faces); mutually exclusive with `PSOURCE` (auto-undefines it)

---

## Input / Output

- `AVERAGES` — Time-averaged output fields
- `AVERAGES_K` — Time-averaged vertical mixing fields
- `OUTPUTS_SURFACE` — Surface-only output snapshot fields
- `HOURLY_VELOCITIES` — Additional hourly velocity output

---

## Restart

- `EXACT_RESTART` — Write extra fields for exact bit-reproducible restart

---

## Diagnostics

Budget diagnostics for tracers, momentum, vorticity, energy, etc.

- `DO_NOT_OVERWRITE` — Keep existing diagnostic files (append mode)
- `RESTART_DIAGS` — Restart diagnostic accumulation from file
- `DIAGNOSTICS_TS` — Tracer equation budget terms
  - `DIAGNOSTICS_TS_ADV` — Use advective (flux) formulation
  - `DIAGNOSTICS_TS_MLD` — Integrate budgets over the mixed layer
    - `DIAGNOSTICS_TS_MLD_DENS` — Use density criterion for MLD
- `DIAGNOSTICS_TSVAR` — Tracer variance budget _(auto-activates DIAGNOSTICS\_TS and DIAGNOSTICS\_TS\_ADV)_
- `DIAGNOSTICS_UV` — Momentum equation budget terms
- `DIAGNOSTICS_VRT` — Depth-mean vorticity budget
- `DIAGNOSTICS_KE` — Kinetic energy budget
- `DIAGNOSTICS_BARO` — Barotropic kinetic energy
- `DIAGNOSTICS_PV` — Potential vorticity budget
- `DIAGNOSTICS_DISS` — Dissipation _(auto-activates DIAGNOSTICS\_PV)_
- `DIAGNOSTICS_EDDY` — Reynolds-stress eddy terms

---

## Biology

Exactly one biogeochemical model must be chosen when `BIOLOGY` is defined.

- `BIOLOGY` — Biogeochemical modeling
  - `PISCES` — 24-component PISCES biogeochemical model
    - `DIURNAL_INPUT_SRFLX` — Diurnal cycle modulation of shortwave input
    - `key_pisces` — Activates PISCES code _(auto-derived)_
    - `key_ligand` — Include ligand cycling _(auto-derived)_
    - `key_pisces_quota` — Quota-based nutrient model
    - `key_pisces_npzd` — Simplified NPZD-PISCES variant
    - `key_sediment` — PISCES sediment module
  - `BIO_NChlPZD` — 5-component NPZD model
    - `OXYGEN` — Oxygen tracer
  - `BIO_N2ChlPZD2` — 7-component NPZD model
  - `BIO_BioEBUS` — 12-component NPZD model (EBUS-oriented)
    - `NITROUS_OXIDE` — N2O tracer
  - `DIAGNOSTICS_BIO` — Store biological flux diagnostics; force-undefined under `PISCES` when
    `XIOS` is also defined
    - `key_trc_diaadd` — Additional PISCES tracer diagnostics _(auto-derived from
      `DIAGNOSTICS_BIO`, requires PISCES; also force-undefined under `PISCES`+`XIOS`)_

---

## Stations Output

- `STATIONS` — High-frequency output at fixed station points
  - `ALL_SIGMA` — Output all sigma levels at stations

---

## USGS Sediment Model

- `SEDIMENT` — USGS sediment transport model
  - `SUSPLOAD` — Suspended load transport
  - `BEDLOAD` — Bedload transport

    Bedload transport formula (choose one; default is `BEDLOAD_VANDERA` with wave forcing,
    `BEDLOAD_WULIN` otherwise):

    - `BEDLOAD_VANDERA` — Van der A formulation
    - `BEDLOAD_MPM` — Meyer-Peter-Muller formulation
    - `BEDLOAD_WULIN` — Wu & Lin formulation
    - `BEDLOAD_MARIEU` — Marieu formulation

    Bedload interpolation scheme (choose one; `BEDLOAD_UP1` is the default):

    - `BEDLOAD_UP1` — 1st-order upwind interpolation _(default)_
    - `BEDLOAD_UP5` — 5th-order upwind interpolation
    - `BEDLOAD_WENO5` — 5th-order WENO interpolation

    Avalanching (slope) scheme (choose one; `SLOPE_LESSER` is the default):

    - `SLOPE_LESSER` — Avalanching formulation, Lesser et al. (2004) _(default)_
    - `SLOPE_NEMETH` — Avalanching: Nemeth et al. (2006)
    - `SLOPE_KIRWAN` — Avalanching: Kirwan formulation

  - `MORPHODYN` — Morphodynamics (bed evolution)
  - `ANA_SEDIMENT` — Analytical sediment ripple and bed parameters
  - `ANA_BPFLUX` — Analytical kinematic bottom flux of sediment
  - `COHESIVE_BED` — Purely cohesive (mud) bed
  - `MIXED_BED` — Mixed cohesive/non-cohesive bed
  - `SED_FLOCS` — Flocculation of cohesive particles _(requires COHESIVE\_BED or MIXED\_BED)_
  - `SED_DEFLOC` — De-flocculation process
  - `FLOC_TURB_DISS` — Turbulent dissipation of flocs
  - `FLOC_BBL_DISS` — BBL dissipation of flocs
  - `SED_TAU_CD_CONST` — Constant drag coefficient for bed stress
  - `TAU_CRIT_WULIN` — Wu & Lin critical shear stress formulation
  - `SED_DENS` — Activate the effect of suspended sediment on the density (Warner et al. 2008)

Other:

- `NO_BSTRESS_UPWIND` — Disable the upwind bottom-stress scheme in sediment transport; only has
  an effect with `SEDIMENT` && `BEDLOAD` defined and `BEDLOAD_VANDERA` undefined

---

## MUSTANG Sediment Model

- `MUSTANG` — MUSTANG sediment dynamics model _(auto-activates SUBSTANCE, USE\_CALENDAR,
  TEMPERATURE, SALINITY, key\_noTSdiss\_insed and key\_nofluxwat\_IWS; 
  force-undefines SEDIMENT —  mutually exclusive with the USGS SEDIMENT model)
  - `key_noTSdiss_insed` — Temperature, salinity and other dissolved variables are not
    computed in the sediment: they keep constant values and no dissolved-variable fluxes
    are exchanged between water and sediment _(auto-derived)_
  - `key_nofluxwat_IWS` — No water flux exchange between water and sediment
    (recommended together with `key_noTSdiss_insed`) _(auto-derived)_
  - `key_MUSTANG_V2` — MUSTANG version 2 (without this key, version 1 is used)
  - `key_MUSTANG_bedload` — MUSTANG bedload transport (requires `key_MUSTANG_V2`)
  - `key_MUSTANG_flocmod` — MUSTANG flocculation model
  - `key_MUSTANG_slipdeposit` — Sliding (avalanching) fluxes for deposited sediment
  - `key_MUSTANG_splitlayersurf` — Split surface sediment layers for a regular, precise
    discretization at the surface
  - `key_tauskin_c_upwind` — Upwind scheme for current-induced bottom shear stress

  - `key_tauskin_c_center` — Compute bottom shear stress from the friction velocity
    directly at the rho-point (cell centre)
  - `key_tauskin_c_ubar` — Use depth-averaged velocity in the bottom shear stress
    computation (combinable with either scheme above)
  - `MORPHODYN` — Morphodynamics (bed evolution)
  - `SED_DENS` — Activate the effect of suspended sediment on the density

---

## Bottom Boundary Layer (BBL) Options

- `BBL` — Bottom Boundary Layer parameterization
  - `BBL_BREAKING_STIR` — Wave breaking stirring in BBL (requires `MRL_WCI` — guarded by
    `#if defined MRL_WCI && defined BBL_BREAKING_STIR` in bbl.F; has no effect without `MRL_WCI`)
  - `COASTAL_BSTRESS_TKE` — Inject wave-breaking TKE into the BBL bottom-stress computation
    (requires `SEDIMENT` and `GLS_MIXING`; only active inside the `#ifdef SEDIMENT` skin-stress
    block of get_vbc.F)
  - `ANA_BSEDIM` — Analytical bed parameters (when SEDIMENT not defined)
  - `ANA_WWAVE` — see Wave-Current Interactions above _(auto-defined as the BBL wave-forcing
    default when none of `OW_COUPLING`, `WAVE_OFFLINE` or `WKB_WWAVE` is defined)_
  - `Z0_BL` — Bedload roughness for ripple predictor _(auto-derived: on if `SEDIMENT` is also
    defined — not independently settable)_
  - `Z0_RIP` — Bedform (ripple) roughness for sandy beds _(auto-derived from `Z0_BL`)_
  