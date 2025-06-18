# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added
- Issue #361 : Integration of slight reformulation of zooplankton grazing 
  according to prey size done in PISCES standard version (#340) into quota version
- BENCH : Add performance tracking

### Fixed

- BENCH : do not put report status to True for reference variant to avoid
  to mark test passed even if not (Issue #342)
- BENCH : put jobcomp.log in results directory even if build fail (Issue #341)
- BENCH : remove openmp reproducibility check on SHOREFACE case (Incident #358)
- BENCH : fix typo AGRIF_2W to AGRIF_2WAY in realist and vortex json files (Issue #368)

- PISCES : Fixed error on diagnostic ligands and add Fe2+ oxydation rate (Issue #371)

- DIAGNOSTICS_EDDY & not XIOS : fix double comma in ncscrum.h (Issue #362)


### Changed

- Compilation : update on jobcomp (support for ifx and different version of gfortran, 
  cleaning exit status, see !172 and Issue#176)
- MUSTANG, SUBSTANCE : separate reading of substance and mustang
  namelist (Issue #354)
- LOGFILE : Change LOGFILE cppkey behavior by enabling to choose filename in
  croco.in (Issue #330)

### Deprecated


### Removed

- MUSTANG : remove key_MUSTANG_debug cppkey (Issue #346)

### Other


### Contributors on this release

- Contributors already on board : 
  R. Benshila, G. Cambon, M. Caillaud, S. Le Gac, P. Marchesiello 

- New contributors : 
  M. Plus, M. Schreiber 