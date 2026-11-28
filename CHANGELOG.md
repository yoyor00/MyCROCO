# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added

- LICENSE : Clarify license (#7)

- STOGEN : add stochastic parametrizations (Issue #301)

- BENCH : Add performance tracking (Issue #378 and #423)
- BENCH : Add missing plot scripts for test cases (#379)

- STATION : Add TEMPERATURE cppkey for stations (#445)

- TIDES_MAS : Add online tidal prediction (SHOM Simon method) for currents and
  potential, in addition to elevation; works now with both
  ZCLIMATOLOGY/M2CLIMATOLOGY and Z_FRC_BRY/M2_FRC_BRY (#471)

- SCRIPTS_COUPLING : deal with hours (in addition to year, month, day). (#593)

- OUTPUTS : add z\_r and z\_w to history and averages files, activated with
  out\_his\_z\_r, out\_his\_z\_w, out\_avg\_z\_r, out\_avg\_z\_w in the
  auxiliary fields namelist blocks (Issue #543)

- SUBSTANCE : substance boundary conditions (`obc_cv_name`) now support
  `'nullGradient'` and `'none'`; `'none'` falls back to the uniform
  `cobc_wat` concentration. Left empty, `obc_cv_name` now defaults to the
  substance's own name (`name_var`), reading it from the bry file if
  present, same as `init_cv_name` for initial conditions below (#590)

- SUBSTANCE : substance initial conditions (`init_cv_name`) now support
  `'none'`, which skips reading the initial file for that substance and
  falls back on the uniform `cini_wat` concentration, mirroring the
  `obc_cv_name`/`cobc_wat` mechanism above. Left empty, `init_cv_name`
  keeps defaulting to the substance's own name, as before. Whenever a
  substance's initial condition is not read from file (`'none'`, or a
  named variable missing from the file), it now always falls back on
  `cini_wat` directly; the previous test-case/depth-dependent gap-filling
  for substances is removed (#590)

- SUBSTANCE : fixed (non-advected) substances now have a working initial
  condition too: `init_cv_name_fix` and `cini_wat_fix` were declared and
  read from the namelist but never actually used (their array, `cvfix_wat`,
  was unconditionally hard-set to 0) (#590)

- SUBSTANCE : `cv_rain` (substance concentration in rainwater) is now
  applied in `bulk_flux.F`, adding `prate * cv_rain` to each substance's
  surface flux, the same way rain already contributes to the temperature
  flux there; only active if CPP keys SALINITY and RAIN_FLUX are both
  defined (#590)

- SUBSTANCE : `flx_atm` (uniform atmospheric dry deposition) is now
  applied in `bulk_flux.F` too, as a constant term added to each
  substance's surface flux (CPP key SUBSTANCE only; combines with the
  `cv_rain` wet-deposition term above when RAIN_FLUX is also defined).
  (#590)

- SUBSTANCE : new `t90` parameter (mortality/decay time, in hours) for
  dissolved, non-constitutive particulate and sorbed particulate
  substances, porting MARS3D's `t90`/`f90` formulation: the
  concentration is multiplied every time step by a decay factor
  `f90 = exp(-dt*ln(10)/(t90*3600))` `t90` = 0. (default)
  disables decay, as before this was added. (#590)
  
- MUSTANG : add Van der A et al. (2013) wave skewness/asymmetry bedload
  formulation (zero-current case), activated with l\_vandera in
  namsedim\_bedload (key\_MUSTANG\_V2 only); add VILAINE\_VANDERA test
  case in BENCH (#589)

- TIMINGS : add wall-clock timing regions in OCEAN/timings.F90 (nested regions,
  per-rank JSON output), enabled at build time with `jobcomp --timings`

### Fixed
- RAIN_FLUX : fix history and average outputs of freshwater flux when RAIN_FLUX is activated.
  Before the fix, the outputs were zero everywhere (Issue #549) 

- RAIN_FLUX : fix MPI reproducibility when BULK_FLUX and RAIN_FLUX were activated (Issue #549)

- BULK_FLUX : fix compilation error when SALINITY was not defined (Issue #549)

- SUBSTANCE/MUSTANG/RAIN_FLUX : fix restart not bit-exact (Issue #610)

- OUTPUTS : fix output time axis inconsistency when appending to an existing
  restart/history/averages/diagnostics/surf/abl/station file with a different
  `start_date` than the one used when the file was created (#611)

- PSOURCE_NCFILE : fix PSOURCE being defined together with PSOURCE_MASS,
  which caused the river discharge Qsrc to be applied twice to u/v in
  step3d_uv2.F (once through PSOURCE_MASS's mass-flux treatment and once
  through PSOURCE's momentum-based treatment) (Issue #599)

- NC4PAR : fix wrong record index on MPI ranks > 0 when reopening an
  existing history/restart/diagnostics file to append records (e.g. on
  restart); non-root ranks computed total_rec one off from rank 0,
  desyncing collective scalar writes (time, time_step) across ranks
  (Issue #596)

- USE_CALENDAR : fix off-by-one in ntimes computed from start_date/end_date
  (missing -1, and int() instead of NINT() truncating float rounding
  noise) (Issue #596)

- WET_DRY : fix restart not bit-exact with WET_DRY; umask_wet/vmask_wet
  are now saved to and restored from the restart file (Issue #262)

- ANA_INITIAL : fix build failure with AGRIF + ANA_INITIAL (#609)

- TIDES_MAS : stop with a clear error instead of silently using an
  uninitialized index when a tide period has no match in the MAS
  catalogue; POT_TIDES : add missing tides.h include in prsgrd.F
  (Issue #595)

- TIDES_MAS : fix AGRIF compilation failure caused by `mod_tides_mas`
  re-exporting `param.h` constants and a module array named `omega`
  colliding with the `omega` subroutine (#605)

- AGRIF `conv` failure, 72-column fixed-form rewrap in ana_grid.F and ana_initial.F (#613 #612)

- DIAGNOSTICS : tracer budget fixes (Issue #320)
	- Taver\_mld averages were not accumulated, and Taver\_mld was not computed without t3dmix routines
	- DIAGNOSTICS\_TSVAR terms were multiplied twice by the cell volume
	- DIAGNOSTICS\_PV right-hand side (Trhs) kept in (unit of tracers) s-1
	- MLD budget closure: same MLD index at both time levels with GLS\_MIXING, THmix from split upstream advection, TVmix with implicit sponge
	- sediment settling/erosion diagnostics are now cell-volume integrals like the other terms
	- TForc nudging diagnostic uses the actual nudging increment (and tclmzon with ZONAL\_NUDGING)
	- no empty-named Tcrit\_mld variables written for passive tracers
	- MLD index of the previous time step (kbl\_nstp) saved in and read from the restart file, so that the first MLD budget after a restart is consistent

- LMD_SKPP2005 / LMD_BKPP2005 : hbl and hbbl were always read from the first record of the restart file instead of the restart record (#320 and #591)

- MUSTANG : lateral erosion feature fluxes in "dry cell" were counting twice in 
  water concentration and last index of current was wrong (Issue #349)
- MUSTANG : fix vertical axis in sediment bed mismatch when using choice_nivsed_out 
  non equal to 1 and initialisation from file and restart (Issue #469)
- MUSTANG : removed the redefinition of Hm in initMUSTANG to prevent silent 
  restart inconsistencies with MORPHODYN, update testcase plot script 
  accordingly (#470)
- MUSTANG : Fix MUSTANG restart bit-exactness (#562)

- Cleaning : typo in ncscrum.h SALINTY instead of SALINITY (#397)
- Cleaning : remove module_qsort.F90 never used            (#394)
- Cleaning : useless sponge option in croco.in.1 (#436)
  
- SCRIPTS: fix EXACT_RESTART handling in Plurimonths_scripts (#475)
- SCRIPTS: create_config.bash -f no longer asks to overwrite an existing work
  directory, so non-interactive runs do not block on stdin

- BENCH : Fix report check status in case of several files (#498)
- BENCH : Fix label in plot_realist.py (#494)
- BENCH : Fix AGRIF handling with a per-test-case AGRIF_FixedGrids.in,
  detected from the resolved cppdefs.h instead of sharing OCEAN's file
  across all cases (#489 and #576)
- BENCH : Fix --build-ref crashing when a check-output or namelist path has subdirectories (#606)
- BENCH : Fix --build-ref crashing when results/last is a dangling symlink, 
  replace it with helpers.replace_symlink() instead (#607)

- NBQ : Fix index when computing total depth cff2 while enforcing consistency between 
  2d and 3d U-momentum for northern open boundary conditions when 
  QDM_OBC_TANG_CORRECT is activated (#508).

- XIOS : fix wrong name for mask_rho in field_def_croco.xml_full_withcpp (#513)

- WAVEMAKER : fix use of wavemaker spectrum from data (bulk wave parameters not 
  initialized in this case). This was done through key ROGUE_WAVES, now changed 
  to WAVE_MAKER_DATA (#518)

- SCRIPTS_COUPLING : Fix NCO module load/unload handling when module is not available (#529)
- COUPLING : Corrects typo in WW3 modified function w3ogcmmd.ftn for BBL related variables (#537)
- COUPLING : Fix missing restart file update for OASIS (#501). 
- COUPLING : Fix PSMILE_INCDIR definition in jobcomp to support all OASIS3-MCT versions (#480)
- ABL1D : Fix ABL1D exact restart and dates in ABL his/avg files (#553)

- OBC_PATM : fix wrong i-index in get_bry.F, zetabry_west atmospheric pressure correction used Istr instead of IstrR (#530)

- Fix time in surf average output file (#388)
- Fix averages over a single time step (nwrt\*\_avg = 1) in diags, diags\_pv, diags\_eddy,
  surf, ABL and bio-diags average files: fields kept accumulating instead of being reset (#320)
- Fix grid variables writing in average and diagnostic file (#522)
- Fix bvf fields units in history and average files (#567)

- Fix ABORT MPI when problem in reading netcdf (#167)

- PISCES : 
  - Switch from relative to potential density for surface pH proxy, inconsistent with what done for calcite dissolution (#531)
  - Fix in the calculation of NEW primary productivity with XIOS (#535)
  - Fix the calculation of nitrogen fixation rate with XIOS (#541)
  - Set XIOS write timestep for PISCES identical to ocean (#548)
  - Correct degree to radian conversion (#559)
  - Add the read of pisces namelist for child grid with AGRIF (#561)
  - Enables use of high-frequency qsr (<24 hrs) with PISCES in the case of offline reading (#563) - online calculation of daily mean of qsr needed by PISCES must be implemented
  - Fix no update of POC/GOC prod by diazotropy - minor effect (#583)
  - Fix micro and meso zooplankton o2 consumption with and without XIOS (#619)
  - Fix bug on O2 production diagnostics (PPNEWo2 and PPRego2) introduced since v2.1.0 (#632) 

- jobcomp : Apply CROCO_CFT1 before compiler-branch selection. 
  Update BENCH accordingly.

- VADV_ADAPT_IMP : add a missing endif in pre_step (#554)

- Restart : Fix restart with OBC_PATM (#496)

- KPP : bug fix for the non-local term in surface KPP parameterisation
        wrong sign and + correction to recover the empirical free convection scaling of 0.2

- TIDES : Adds an error message if POT_TIDES is defined but no potentials are present in the frc file (#431)

- MPI_NOLAND :adding rmask for volume, ke, pe global calculations to avoid tiling dependancy (#476)

- NH : Fix rhs typo (#538)

### Changed
- Default namelist values changed following EXACT_RESTART removal: `nrst`
  now defaults to -1 (restart file written once, at the end of the run,
  instead of every 720 time steps) and `nrpfrst` now defaults to 0 (#553)

- OUTPUTS : output variable indices in ncscrum.h no longer differ between
  2D and 3D configurations (#500)

- OUTPUTS : output files are no longer overwritten on restart. Record
  counters of the history, averages and diagnostics files are saved in
  the restart file and restored, so a restart appends to the existing
  file instead of starting a new one. This used to require
  `DO_NOT_OVERWRITE` and `RESTART_DIAGS`, which are now the default
  behavior (and no longer exist). A new namelist option `lresetrec`
  (`&croco_initial`, default false) can be used to reset the record
  counters to zero on restart, e.g. to change an output frequency or
  start a new branch run (Issue #108)

- Vertical coordinate : only the new s-coordinate (Vtransform = 2,
  Shchepetkin & McWilliams 2009) is kept (#396).
  Input files (initial, boundary, climatology) are checked when opened:
  CROCO stops with an explicit error if VertCoordType is OLD/SH94 or if
  Vtransform /= 2. Files without VertCoordType nor Vtransform are assumed
  to use Vtransform = 2. Old input files must be regenerated with the
  current croco_tools or croco_pytools, or used with croco <= v2.1.3.
  See also croco_tools #69 and croco_pytools #118.
  Test cases that previously ran with the old s-coordinate (ANA_JET,
  BASIN, CANYON, EQUATOR, INNERSHELF, INTERNAL, KILPATRICK, OVERFLOW,
  RIVER, SEAMOUNT, SHELFRONT, SOLITON, UPWELLING, VORTEX) now use the new
  one: their results differ from previous versions. VORTEX input files
  have been regenerated with Vtransform = 2.

- CALENDAR : Remove `USE_CALENDAR` CPP key.
  Calendar support is now always active and controlled at runtime via `calendar_type` in the namelist (#452).
  Supported values (case-insensitive): `gregorian` (default) / `standard` /
  `proleptic_gregorian`, `360_day`, `noleap`/`365_day`/`no_leap`,
  `all_leap`/`366_day`, `julian`.
  The five standalone calendar utility files (`toolsectodat.F90`, `tooldatosec.F90`,
  `tooldecompdat.F90`, `tooldatetosec.F90`, `toolorigindate.F90`) are replaced by a
  single Fortran module `tools_calendar.F90`, which:
    - provides explicit interfaces (eliminates scattered `character*19 tool_sectodat`
      external declarations in all `get_*.F` forcing readers);
    - extends `tool_datetosec` to support all calendar types (was Gregorian-only);
    - renames all French variable names to English throughout in this module.
  All `get_*.F` forcing readers now log dates as `yyyy-mm-dd hh:mm:ss` strings
  (was raw floating-point days since origin).
  TOOL_ORIGINDATE error/warning messages now report the affected netCDF
  file and variable (#600).
  TOOL_DECOMPDATE stops with a clear error instead of crashing with a raw Fortran
  "Bad value during integer read" in case of non CF-compliant date (#620)

- Input file croco.in replace by a standard namelist (#497)
  Replace the fixed-format croco.in reader with a Fortran namelist system
  (croco.nml). All configuration parameters are now in structured &croco_*
  namelists.
  A convert_in_to_nml.py utility is provided to easily convert 
  previous croco.in into croco.nml.
  A detailed error message (raw content, iostat/iomsg, expected array size)
  is now written in case of a malformed namelist (#555)
  CROCO now stops with an explicit error if the namelist file is missing,
  instead of silently continuing with default values for every parameter
  (Issue #580)
  The per-source namelist groups are not printed when there is no source (#608)
  The station nc_sta.h file is updated (#622)
  (#552, #568, #570, #572, #573, #578)

- Test cases : Reorganize TEST_CASES/ into per-case subdirectories (TEST_CASES/<CASE>/)
  with standardized uppercase filenames. Update all BENCH jsonc configs
  and plot scripts accordingly, as well as production run scripts.

- AGRIF : update conv version (#510)

- SUBSTANCE : submassbalance feature is now activated only by namelist
  (Issue #347)

- Compilation : update on jobcomp (support for ifx and different version of gfortran, 
  cleaning exit status, see !172 and Issue#176), 
  update NETCDF paths default from nf-config and nc-config (Issue #473)

- MUSTANG, SUBSTANCE : separate reading of substance and mustang
  namelist (Issue #354)
- MUSTANG : review lateral erosion feature (Issue #349)
- MUSTANG : change activation of horizontal fluxes correction for sand (Issue #352)
- MUSTANG : add MRL_WCI and OW_COUPLING handling (#348 and #464)
- MUSTANG : change init and roughness length computation (#562)
- MUSTANG : dispatch per-test-case namelists into TEST_CASES/ and declare defaults in code (#353, #594)

- LOGFILE : Change LOGFILE cppkey behavior by enabling to choose filename in
  croco.in (Issue #330)

- DIAGNOSTICS : 
	- cleaning, simplifications and updates of momentum-based diagnostics (DIAGNOSTICS\_KE, DIAGNOSTICS\_VRT, DIAGNOSTICS\_M) (Issue #388)
	- kinetic energy budget is now 3d
	- momentum, energy and tracer diagnostics (DIAGNOSTICS\_TS and DIAGNOSTICS\_TS\_MLD) terms are now volume-integrated (Issue #320)
	- DIAGNOSTICS\_TS\_MLD\_DENS : MLD depths (density and temperature criteria) are both computed at w-levels, and the density criterion is used whatever the vertical mixing scheme (Issue #320)
	- MLD depths (density and temperature criteria) are both computed starting from a reference depth (10 m by default) instead of the surface. The reference depth can be set with the namelist parameter mld_depth_ref (Issue #545)

- BIOLOGY : PISCES is now the default biogeochemical model (Issue #461)

- WAVES : wepb output now includes roller delay when roller is activated (#579)
          the old wepb output can be diagnose from wepb and wepr

- BULK_FLUX : Update wasp bulk flux parametrization, 
  cppkey BULK_WASP (Issue #453)

- RIVER test case updated to pass PSOURCE_MASS with an EXP_SHAPE vertical 
  distribution of flow, enabling a transition from the AKIMA scheme to 
  SPLINES (#478)

- OMEGA : Add a condition on the NBQ_MASS key for some terms of the first 
  part of the computation of omega (#447)

- BIOLOGY : Improvements and bug fix (sedmat+sedinorg) in the PISCES sediment module (#468)
- XIOS : Align density anomaly computation with native netCDF writer,
  now with respect to 1000 kg m-3 (Issue #484)

- XIOS : Unmask output grid/geometry variables so that they keep valid
  values on land. Useful for offline diagnostics mimicking CROCO's way. 

- Compilation : replace the legacy `include 'mpif.h'` by the standard
  `use mpi` Fortran module in OCEAN, MUSTANG, XIOS, PISCES, AGRIF and
  STOGEN (#372)

- Build : replace cross_matrix.F, the Fortran tool generating Make.depend,
  by make_depend.py, a stdlib-only Python replacement; it fixes missing
  dependencies that made incremental (re)builds unreliable and comes with
  unit tests run in CI (Issue #566)

### Deprecated


### Removed
- EXACT_RESTART cpp key has been removed, exact bit-reproducible restart
  is now always active by default (#553)

- AVERAGES_K cpp key has been removed (Issue #374). 
   If deactivated, averaged vertical mixing coefficients AKv, AKt and AKs were not saved.
   Removed because redundant with booleans out_avg_akv, out_avg_akt and out_avg_aks in namelist. 

- DO_NOT_OVERWRITE and RESTART_DIAGS cpp keys have been removed, their
  behavior (appending to existing output files on restart) is now the
  default (Issue #108)

- NEW_S_COORD, cpp key has been removed (if deactivated, old s-coordinate (Vtransform = 1, 
    Song & Haidvogel 1994) was chosen , Now Vtransform=2 is used by default (Shchepetkin & McWilliams 2009))
    the key has also been removed from all TEST_CASES cppdefs (#396)

- Remove duplicate `scrum_time` netCDF variable from history, restart, 
  averages, diagnostics, surf and station output files; only `time` is 
  written now (matches the `time` dimension name). get_initial.F still 
  falls back to reading `scrum_time` when restarting from older files 
  (Issue #189)

- SUBSTANCE_SUBMASSBALANCE cpp key has been removed, feature is activated 
  by boolean in namelist (Issue #347)
- MUSTANG : 
  - remove key_MUSTANG_lateralerosion replace by a boolean in 
    namelist (Issue #349)
  - remove key_sand2D, activation only by a boolean in 
    namelist (Issue #351 and #525)
  - remove MUSTANG_CORFLUX replace by a boolean in 
    namelist (Issue #352)
  - remove key_MUSTANG_debug cppkey (Issue #346)
  - remove file scalars_F90.h, not used (Issue #382)
  - remove key_tauskin_c_ubar key_tauskin_c_center key_tauskin_c_upwind
    replace by booleans in namelist (Issue #348)
  - remove key_MUSTANG_slipdeposit replace by boolean l_slipdeposit in
    namelist (Issue #350)
  - remove key_MUSTANG_bedload, bedload transport (l_bedload_n()) now only
    conditioned by key_MUSTANG_V2 (Issue #385)
  - remove key_MUSTANG_splitlayersurf replace by boolean l_splitlayersurf
    in namelist (Issue #350)
  - remove key_MUSTANG_flocmod replace by boolean l_flocmod in namelist
    (Issue #350)

- Test cases CPP keys replace by namelist parameter (#497)
  Remove all test-case CPP guards from the solver. Test-case selection is
  now done at runtime via testcase_name in the namelist, not at compile
  time. cppdefs.h and cppdefs_dev.h are cleaned of all per-case guards.
  cppdefs.h and param.h are now provided by case in TEST_CASES. The regional 
  Benguela example is also copied to OCEAN.
  ana_initial.F : Add comment explaining why testcase_name == '' (realistic
  config, no idealized dispatch) is handled separately from case default
  (fallback for an unrecognized testcase_name) (#560)

- Obsolete, unused or undocumented CPP keys : 
  - FLOATS, deprecated (#296)
  - TS_VADV_FCT was always undef, never used (#390)
  - WET_DRY0 (#393) never used 
  - UV_HADV_TVD, UV_VADV_TVD, W_HADV_TVD, W_VADV_TVD (#391)
  - BVF_MIXING (#398)
  - LMD_NUW_GARGETT, obsolete (#402)
  - ROBUST_DIURNAL_SRFLUX (#405)
  - DUKO_2001 was always def (#407)
  - PRED_COUPLED_MODE was always def (#408) 
  - START_DATE (#417)
  - ICE (#416)
  - DECALPHA (#414)
  - CRAY, VAX, SGI, AIX (#413)
  - AUTOTILING (#411)
  - DEBUG_ARMOR, DEBUG, DIAGNOSTICS_DEBUG, NBQ_HZCORR_DEBUG (#415)
  - PP_MIXING, MY2_MIXING, MY25_MIXING (#418)
  - XCOMM_FORMAT (#419)
  - TR (#395)
  - LMD_SKPP_MONOB never define (#400)
  - LIMIT_UNSTABLE_ONLY is always define (#401)
  - MLCONVEC (#399)
  - TS_VADV_AKIMA and TS_HADV_AKIMA (#392)
  - TENDENCY, DIAGNOSTICS_EK_FULL, DIAGNOSTICS_EK_MLD (Issue #388)

### Other

- Cleaning :
  - remove TEST_CASES/IGW_OA directory with PDFs, namelist, XIOS XML files, and README  (#337)
  - remove files dynparam_f77.h, agrif_ext.h, diag_vars.h, not used (Issue #386)
  - remove files parameter.passivetrc.pisces.h, not used (Issue #387)
  - comments refering to BASIN in step2D.F (#409)
  - remove routine set_HUV1, not used (#410)
  - remove ZETA_DRY_IO cpp key and avoid modifying zeta with bathymetry in output (#406 and #384)
  - typo in diag.F CALENDAR instead of USE_CALENDAR (#412)
  - avoid hard coded define of RI_[H/V]SMOOTH in code moved 
    in cppdefs_dev.h (#403)
  - remove hard coded keys in mpc.F (#404)
  - typo and file mode (#499)

- Support :
  - upgrade ci env (ubuntu, hdf5, netcdf versions, ifx compilers) (#463)
  - use matrix capabilities in gitlab-ci (#519)
  - avoid misnaming ifort/ifx in gitlab-ci (#492)

### Contributors on this release

- Contributors already on board : 
  R. Benshila, M. Caillaud, G. Cambon, N. Ducousso, F. Dufois, S. Jullien, 
  S. Le Gac, P. Marchesiello, C. Nguyen, R. Person, J. Pianezze, S. Treillou, 
  J. Gula, C. Mazoyer, G. Morvan, M. Le Corre, P. Penven, M. Bremond

- New contributors : 
  M. Belharet, J.-M. Brankart, D. Gourves, Q. Jamet, L. Weiss,
  M. Plus, M. Schreiber, A. Zribi, B. Lemieux-Dudon, C. Menu, E. Le Bouedec
  S. Theetten, A. Mittermair, L. Maillard, N. Herry
