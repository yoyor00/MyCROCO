# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added


### Fixed

- BENCH : do not put report status to True for reference variant to avoid
  to mark test passed even if not (Issue #342)
- BENCH : put jobcomp.log in results directory even if build fail (Issue #341)
- BENCH : remove openmp reproducibility check on SHOREFACE case (Incident #358)

- DIAGNOSTICS_EDDY & not XIOS : fix double comma in ncscrum.h (Issue #362)


### Changed

- SUBSTANCE : submassbalance feature is now activated only by namelist
  (Issue #347)
- Issue #330 : Change LOGFILE cppkey behavior by enabling to choose filename in
  croco.in

### Deprecated


### Removed

- SUBSTANCE_SUBMASSBALANCE cpp key has been removed, feature is activated 
  by boolean in namelist (Issue #347)

### Other


### Contributors on this release

- Contributors already on board : 
  R. Benshila, M.Caillaud, S. Le Gac, P. Marchesiello 

- New contributors : 
