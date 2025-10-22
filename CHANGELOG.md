# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.X] - xxxx-xx-xx

### Added

### Fixed

- BOTTOM STRESS : Incorrect definition of loop indices for calculating
  the bottom stress components in the case rdrg2>0 (Issue #441)

- AGRIF: sponge keyword missing: put it back (even if theoretically useless, it creates a read_inp error), path for online corrected in croco_inter.in,
  AGRIF_Fixed.in not copied from the right directory: corrected, copy croco_frc for all domains for tides  (solve #438)

- BENCH : fix hwloc-ls not mandatory (solve #442)


### Changed


### Deprecated


### Removed


### Other


### Contributors on this release

- Contributors already on board : 
  N. Ducousso, S. Jullien, S. Le Gac, P. Marchesiello

- New contributors : 
