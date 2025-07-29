# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added
- Issue #361 : Integration of slight reformulation of zooplankton grazing 
  according to prey size done in PISCES standard version (#340) into quota version

- BENCH : Add performance tracking (Issue #378 and #423)

### Fixed

- COUPLING : fixes to prevent runtime crash when compiled in full debug mode (Issue #376)
- COUPLING : patm2D was declared twice in case of OW_COUPLING and READ_PATM (Issue #383)
- BENCH : do not put report status to True for reference variant to avoid
  to mark test passed even if not (Issue #342)
- BENCH : put jobcomp.log in results directory even if build fail (Issue #341)
- BENCH : remove openmp reproducibility check on SHOREFACE case (Incident #358)
- BENCH : fix typo AGRIF_2W to AGRIF_2WAY in realist and vortex json files (Issue #368)

- PISCES : Fix dummy line (kt variable) at end in p4zche.F90  (Issue #420)
- PISCES : Fixed error on diagnostic ligands and add Fe2+ oxydation rate (Issue #371)

- DIAGNOSTICS_EDDY & not XIOS : fix double comma in ncscrum.h (Issue #362)

- MUSTANG : lateral erosion feature fluxes in "dry cell" were counting twice in 
  water concentration and last index of current was wrong (Issue #349)

- Cleaning : typo in ncscrum.h SALINTY instead of SALINITY (#397)
- Cleaning : remove module_qsort.F90 never used            (#394)

- Compilation : fix cat "croco_ascii.txt" command in case of relative path


### Changed

- SUBSTANCE : submassbalance feature is now activated only by namelist
  (Issue #347)
- Compilation : update on jobcomp (support for ifx and different version of gfortran, 
  cleaning exit status, see !172 and Issue#176)
- MUSTANG, SUBSTANCE : separate reading of substance and mustang
  namelist (Issue #354)
- MUSTANG : review lateral erosion feature (Issue #349)
- LOGFILE : Change LOGFILE cppkey behavior by enabling to choose filename in
  croco.in (Issue #330)

### Deprecated


### Removed


- SUBSTANCE_SUBMASSBALANCE cpp key has been removed, feature is activated 
  by boolean in namelist (Issue #347)
- MUSTANG : 
  - remove key_MUSTANG_lateralerosion replace by a boolean in 
    namelist (Issue #349)
  - remove key_MUSTANG_debug cppkey (Issue #346)
  - remove file scalars_F90.h, not used (Issue #382)

- Obsolete, unused or undocumented CPP keys : 
  - FLOATS, deprecated (#296)
  - TS_VADV_FCT was always undef, never used (#390)
  - UV_HADV_TVD, UV_VADV_TVD, W_HADV_TVD, W_VADV_TVD (#391)
  - BVF_MIXING (#398)
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

### Other

- Cleaning :
  - remove files dynparam_f77.h, agrif_ext.h, diag_vars.h, not used (Issue #386)
  - remove files parameter.passivetrc.pisces.h, not used (Issue #387)
  - comments refering to BASIN in step2D.F (#409)


### Contributors on this release

- Contributors already on board : 
  R. Benshila, M. Caillaud, G.Cambon, S. Jullien, S. Le Gac, 
  P. Marchesiello, C. Nguyen, R. Person

- New contributors : 
  M. Plus, M. Schreiber 
