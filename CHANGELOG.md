# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.X] - xxxx-xx-xx

### Added
 
- SCRIPTS : update run_croco_inter.bash to handle "USE_CALENDAR" outputs logic (Issue #455)

### Fixed

- BOTTOM STRESS : Incorrect definition of loop indices for calculating
  the bottom stress components in the case rdrg2>0 (Issue #441)

- COUPLING : exchange of some v-grid variables in cpl_prism_get when defined OA_GRID_UV (#450) 

- MUSTANG : error in case of SAND only type and key_MUSTANG_V2 
  not defined (Issue #451)

- SUBSTANCE: submassbalance error if no closed border (Issue #449)

- Cleaning : remove printing variable time_mars (#374)

- AGRIF: sponge keyword missing: put it back (even if theoretically useless, 
  it creates a read_inp error), path for online corrected in croco_inter.in,
  AGRIF_Fixed.in not copied from the right directory: corrected, copy croco_frc for all domains for tides  (solve #438)
- AGRIF : make it usable with USE_CALENDAR (Issue #339)

- BENCH : fix hwloc-ls not mandatory (solve #442)
- BENCH: fix result_pattern used to retrive performance data for performance plot, fix
  results directory for copying jobcomp.log in case of failed compilation (solve #446)

- OUTPUT : fix output with XIOS and using ABL1D (Issue #425)


### Changed

- BULK : change input reading to search for uppercase variable 
  before exiting in error (#454)

### Deprecated


### Removed


### Other

- CI :
  - change container managment for CI test (#437)

### Contributors on this release

- Contributors already on board : 
  N. Ducousso, S. Jullien, S. Le Gac, P. Marchesiello, J. Pianezze

- New contributors : 
