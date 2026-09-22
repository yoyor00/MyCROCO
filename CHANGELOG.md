# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

  ####  ### Added                                                                                                  
                                                                                                                   
    - KNBQ3 : Ajout d'une documentation détaillée (`DOC_KNBQ3.md`) sur le solveur non-hydrostatique compressible   
  KNBQ3.
    - BENCH : Ajout des configurations et cas tests pour AgAc (`AgAc_KNBQ`, `AgAc_KAGRIF`, `AgAc_KABGRIF`) et      
  Canon2D (`CANON2D_HWI`, `CANON2D_KHI`, `CANON2D_TCI`).
    - BENCH : Ajout d'un nouveau guide pratique (`BENCH/Guide_Ajout_Cas_Tests.md`) et d'une variante de compilation
  `mpi-nvfortran`.
  
  ####  ### Fixed 
  
    - KNBQ3 : Correction de l'échange périodique/MPI pour la densité `rho_nbq` dans `k3fast_mass_update.h`.        
    - OCEAN : Correction d'une division par zéro potentielle dans `set_nudgcof.F` (relaxation et éponge) lorsque le
  nombre de points `isp <= 0`.
    - jobcomp : Lancement automatique du script `prepro_KXX.py` lors de la compilation avec AGRIF et les noyaux KXX.
    - BENCH : Prise en charge des patches booléens (ex: activation d'AGRIF) et transmission de l'option AGRIF à    
  `create_config.bash`.
  

### Added

- MERGE : Fusion de la branche `dev2026_NBQ3_merge` comprenant l'intégration du nouveau noyau non-hydrostatique KNBQ3 (K3FAST), l'archivage du solveur historique KNBQ2, l'ajout du mélange vertical TKE3D, et de nouvelles configurations physiques (LES, CANON2D, ISOLITON_DJL).
- NBQ : Integrate new non-hydrostatic kernel (KNBQ3) and archive old KNBQ baseline to KNBQ2.
- LICENSE : Clarify license (#7)

- STOGEN : add stochastic parametrizations (Issue #301)

- BENCH : Add performance tracking (Issue #378 and #423)
- BENCH : Add missing plot scripts for test cases (#379)

- STATION : Add TEMPERATURE cppkey for stations (#445)

### Fixed

- NBQ : Fix compilation errors of the Adams-Bashforth 3 scheme (correction of `K3FAST_AB3` macros in `nbq.h`, `k3fast_qdmuv_update.h` and `k3fast_qdmw_update.h`).
- OMEGA : Resolve vertical velocity differences in `omega.F` by making grid respiration calculation unconditional (not restricted to `#ifdef NBQ_MASS`).
- XIOS : Fix runtime XIOS crashes by adding missing XML field definitions (`Cs_rho`, `Cs_rho_total`, `rho_surf`, `wN`) and correcting grid reference for vertical velocity `w` (from `w_3D` to `rho_3D`) in XML templates.
- MUSTANG : lateral erosion feature fluxes in "dry cell" were counting twice in 
  water concentration and last index of current was wrong (Issue #349)
- MUSTANG : fix vertical axis in sediment bed mismatch when using choice_nivsed_out 
  non equal to 1 and initialisation from file and restart (Issue #469)
- MUSTANG : removed the redefinition of Hm in initMUSTANG to prevent silent 
  restart inconsistencies with MORPHODYN, update testcase plot script 
  accordingly (#470)

- Cleaning : typo in ncscrum.h SALINTY instead of SALINITY (#397)
- Cleaning : remove module_qsort.F90 never used            (#394)
- Cleaning : useless sponge option in croco.in.1 (#436)

- BENCH : Fix report check status in case of several files (#498)
- BENCH : Fix label in plot_realist.py (#494)

- NBQ : Fix index when computing total depth cff2 while enforcing consistency between 
  2d and 3d U-momentum for northern open boundary conditions when 
  QDM_OBC_TANG_CORRECT is activated (#508).

- XIOS : fix wrong name for mask_rho in field_def_croco.xml_full_withcpp (#513)

- WAVEMAKER : fix use of wavemaker spectrum from data (bulk wave parameters not 
  initialized in this case). This was done through key ROGUE_WAVES, now changed 
  to WAVE_MAKER_DATA (#518)

- Fix time in surf average output file (#388)

### Changed

- SUBSTANCE : submassbalance feature is now activated only by namelist
  (Issue #347)

- Compilation : update on jobcomp (support for ifx and different version of gfortran, 
  cleaning exit status, see !172 and Issue#176)

- MUSTANG, SUBSTANCE : separate reading of substance and mustang
  namelist (Issue #354)

- MUSTANG : review lateral erosion feature (Issue #349)

- MUSTANG : change activation of horizontal fluxes correction for sand (Issue #352)

- LOGFILE : Change LOGFILE cppkey behavior by enabling to choose filename in
  croco.in (Issue #330)

- DIAGNOSTICS : 
	- cleaning, simplifications and updates of momentum-based diagnostics (DIAGNOSTICS\_KE, DIAGNOSTICS\_VRT, DIAGNOSTICS\_M) (Issue #388)
	- kinetic energy budget is now 3d
	- momentum and energy diagnostics are saved as cell-volume integrals

- BIOLOGY : PISCES is now the default biogeochemical model (Issue #461)

- WKB_WWAVE : variable name wepb0 or wepb directly manage in wrt_his 
  and not in cppdefs_dev.h (#465)

- BULK_FLUX : Update wasp bulk flux parametrization, 
  cppkey BULK_WASP (Issue #453)

- RIVER test case updated to pass PSOURCE_MASS with an EXP_SHAPE vertical 
  distribution of flow, enabling a transition from the AKIMA scheme to 
  SPLINES (#478)

- OMEGA : Add a condition on the NBQ_MASS key for some terms of the first 
  part of the computation of omega (#447)

- BIOLOGY : Improvements and bug fix (sedmat+sedinorg) in the PISCES sediment module (#468)

### Deprecated


### Removed

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

### Contributors on this release

- Contributors already on board : 
  R. Benshila, M. Caillaud, G. Cambon, N. Ducousso, F. Dufois, S. Jullien, 
  S. Le Gac, P. Marchesiello, C. Nguyen, R. Person, J. Pianezze, S. Treillou, 
  J. Gula

- New contributors : 
  J.-M. Brankart, D. Gourves, Q. Jamet, L. Weiss,
  M. Plus, M. Schreiber, A. Zribi, B. Lemieux-Dudon, C. Menu, E. Le Bouedec
  S. Theetten
