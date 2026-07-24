# CROCO CPP Keys Reference

Companion to `croco_nml_full.md`. 
Advanced options (rarely changed) are in `cppdefs_dev.h`.
Keys automatically derived from others should not be set here.

Indented keys are only meaningful when their parent key is defined.

---


## Parallelization and I/O Server

Only one of `OPENMP` / `MPI` may be defined at a time. 

- `OPENMP` — Open-MP shared-memory parallelism
- `MPI` — MPI distributed-memory parallelism
  - `PARALLEL_FILES` — Parallel I/O: each MPI task writes its own file
  - `NC4PAR` — NetCDF-4 parallel I/O (requires HDF5/NC4)
  - `MPI_NOLAND` — Skip land-only tiles (load balancing)
  - `MPI_TIME` — Time MPI communications
  - `XIOS` — Use the external XIOS I/O server (needs MPI)
- `OPENACC` — OpenACC GPU 

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
  - `SSH_TIDES` — Read and apply tidal sea-level data
  - `UV_TIDES` — Read and apply tidal current data
  - `POT_TIDES` — Include tidal potential forcing
  - `TIDES_MAS` — Harmonic composition to build tide elevation from SHOM (Simon method). 
    (requires `USE_CALENDAR`)
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
- `FILLVAL` — Write fill value in masked land points

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
- `M2FILTER_NONE` — Disable barotropic time filter; 
- `ZONAL_NUDGING` — Zonal mean flow nudging (used with JET)

---

## Equation of State

- `SALINITY` — Salinity as an active tracer
- `NONLIN_EOS` — Nonlinear equation of state
- `SPLIT_EOS` — Split EOS into adiabatic + compressible parts to reduce pressure-gradient errors
- `NO_RESET_RHO0` — Disable automatic reset of rho0 when using linear EOS
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
  - `READ_PATM` — Read surface atmospheric pressure
    - `OBC_PATM` — Apply atm pressure gradient at open boundaries
  - `ABL1D` — 1D ABL model (column model coupled to ocean)
    - `ANA_ABL_LSDATA` — Analytical large-scale ABL data
    - `ANA_ABL_VGRID` — Analytical ABL vertical grid
    - `STRESS_AT_RHO_POINTS` — Compute surface stress at rho points _(auto-derived)_
    - `ABL_NUDGING` — Nudge ABL towards large-scale data _(auto-derived)_
    - `ABL_NUDGING_DYN` — Nudge dynamics (wind) _(auto-derived)_
    - `ABL_NUDGING_TRA` — Nudge thermodynamics (temperature) _(auto-derived)_
    - `ABL_DYN_RESTORE_EQ` — Restore ABL dynamics towards equilibrium
    - `ABL_NO_OCEAN_FEEDBACK` — Suppress ocean SST feedback to the ABL (one-way ABL > ocean coupling, used in test case KILPATRICK)
    - `SFLUX_CFB` — Current feedback on wind stress (modifies surface stress based on ocean surface velocity)
  - `QCORRECTION` — Linear bulk SST correction of heat flux _(without BULK\_FLUX only)_
  - `SFLX_CORR` — Freshwater flux correction towards model SSS _(without BULK\_FLUX only)_
  - `SFLX_CORR_COEF` — Prescribed correction coefficient _(without BULK\_FLUX only)_
- `SFLUX_CFB` — Current feedback on wind stress (standalone, without `BULK_FLUX`)
- `SEA_ICE_NOFLUX` — Zero flux under sea ice

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

- `UV_VADV_SPLINES` — Splines vertical advection for momentum
- `UV_VADV_WENO5` — WENO5 vertical advection for momentum

---

## Vertical Velocity Advection (W)

Only relevant with NBQ or wave-resolving simulations. 
Set automatically by `cppdefs_dev.h` when `UV_HADV_WENO5` is defined.

- `W_HADV_WENO5` — 5th-order WENOZ horizontal advection for W
- `W_VADV_WENO5` — 5th-order WENOZ vertical advection for W

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

- `TS_VADV_SPLINES` — Splines vertical advection for tracers
- `TS_VADV_WENO5` — WENO5 vertical advection for tracers

---

## Semi-Implicit Vertical Advection

- `VADV_ADAPT_IMP` — Adaptive semi-implicit vertical advection

---

## Sponge Layers

Enhance viscosity and diffusivity near open boundaries to damp spurious reflections.

- `SPONGE` — Sponge layers near lateral open boundaries;
  - `NO_SPONGE_GRID` — Disable reading sponge coefficients from the grid file (use analytical, used for INNERSHELF)

---

## Bottom Stress

- `LIMIT_BSTRESS` — Limit bottom stress to prevent negative transport
  - `NO_LIMIT_BSTRESS` — Disable the bottom stress limiter (overrides `LIMIT_BSTRESS` default, used for INNERSHELF)
- `NO_BSTRESS_UPWIND` — Disable upwind bottom-stress scheme in sediment transport
- `BSTRESS_FAST` — Compute bottom stress at 3D fast (baroclinic) time steps

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
  - `GLS_MIXING_3D` — Full 3D GLS (vs. columnar 1D version)
  - `GLS_KOMEGA` — K-omega model (K-epsilon if not defined)
  - `GLS_COASTAL_ROUGHNESS` — Enhanced bottom roughness near coastal boundaries

---

## Wave-Current Interactions

Wave forcing source (choose one; also used with `BBL` or `MUSTANG` independently of `MRL_WCI`):

- `WAVE_OFFLINE` — Read pre-computed wave fields from a file _(default when `MRL_WCI` is active without other source)_
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

Wave breaking parameterization (used with `MRL_WCI` or `WKB_WWAVE`; choose one):

- `WAVE_BREAK_CT93` — Thornton & Guza (1983) wave breaking _(default)_
- `WAVE_BREAK_TG86` — Church & Thornton (1993) wave breaking
- `WAVE_BREAK_TG86A` — Alternate Church-Thornton formulation

Wave-current interaction model:

- `MRL_WCI` — Wave-Current Interaction model (Mellor-type)
  - `MRL_CEW` — Wave-driven longshore current energetics
  - `WAVE_STREAMING` — Bottom wave streaming (wave-induced current)
  - `WAVE_RAMP` — Ramp wave forcing at start
  - `WAVE_DRY` — Suppress wave forcing in cells shallower than `D_wavedry` (requires `WET_DRY`)

---

## Wave Maker

Wave-resolving internal wave generation (ideal cases, e.g. SANDBAR, SWASH, RIP).

- `WAVE_MAKER` — Wave-maker boundary source
  - `WAVE_MAKER_SPECTRUM` — Multi-frequency spectrum input
  - `WAVE_MAKER_BICHROMATIC` — Bichromatic (two-frequency) wave input (e.g. SWASH GLOBEX B2/B3)
  - `WAVE_MAKER_JONSWAP` — JONSWAP spectrum input
  - `WAVE_MAKER_DSPREAD` — Directional spreading
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
- `PSOURCE_MASS` — Inject river as vertical volume flux at rho points (instead of flux through u,v faces)

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
  - `DIAGNOSTICS_BIO` — Store biological flux diagnostics _(auto-derived)_
    - `key_trc_diaadd` — Additional PISCES tracer diagnostics _(requires PISCES)_

---

## Stations Output

- `STATIONS` — High-frequency output at fixed station points
  - `ALL_SIGMA` — Output all sigma levels at stations

---

## USGS Sediment Model

- `SEDIMENT` — USGS sediment transport model
  - `SUSPLOAD` — Suspended load transport
  - `BEDLOAD` — Bedload transport
  - `MORPHODYN` — Morphodynamics (bed evolution)
  - `ANA_SEDIMENT` — Analytical sediment ripple and bed parameters
  - `ANA_BPFLUX` — Analytical kinematic bottom flux of sediment
  - `BEDLOAD_MPM` — Bedload: Meyer-Peter-Muller formulation _(mutually exclusive)_
  - `BEDLOAD_WULIN` — Bedload: Wu & Lin formulation _(mutually exclusive)_
  - `BEDLOAD_MARIEU` — Bedload: Marieu formulation _(mutually exclusive)_
  - `BEDLOAD_WENO5` — Bedload: WENO 5th-order advection _(mutually exclusive)_
  - `SLOPE_NEMETH` — Avalanching: Nemeth et al. (2006)
  - `SLOPE_LESSER` — Avalanching: Lesser et al. (2004)
  - `COHESIVE_BED` — Purely cohesive (mud) bed
  - `MIXED_BED` — Mixed cohesive/non-cohesive bed
  - `SED_FLOCS` — Flocculation of cohesive particles _(requires COHESIVE\_BED or MIXED\_BED)_
  - `SED_DEFLOC` — De-flocculation process
  - `FLOC_TURB_DISS` — Turbulent dissipation of flocs
  - `FLOC_BBL_DISS` — BBL dissipation of flocs
  - `SED_TAU_CD_CONST` — Constant drag coefficient for bed stress
  - `TAU_CRIT_WULIN` — Wu & Lin critical shear stress formulation

---

## MUSTANG Sediment Model

- `MUSTANG` — MUSTANG sediment dynamics model
  - `key_MUSTANG_V2` — MUSTANG version 2
  - `key_MUSTANG_bedload` — MUSTANG bedload transport
  - `key_MUSTANG_flocmod` — MUSTANG flocculation model
  - `key_tauskin_c_upwind` — Upwind skin friction for MUSTANG
  - `MORPHODYN` — Morphodynamics (usually off)

---

## Bottom Boundary Layer (BBL) Options

- `BBL` — Bottom Boundary Layer parameterization
  - `BBL_BREAKING_STIR` — Wave breaking stirring in BBL
  - `COASTAL_BSTRESS_TKE` — Inject wave-breaking TKE into the BBL (requires `GLS_MIXING`)
  - `ANA_BSEDIM` — Analytical bed parameters (when SEDIMENT not defined)
  - `ANA_WWAVE` — Analytical wave forcing for BBL
  - `Z0_BL` — Compute bedload roughness for ripple predictor
  - `Z0_RIP` — Bedform (ripple) roughness for sandy beds
  - `Z0_BIO` — Biogenic bedform roughness for silty beds
