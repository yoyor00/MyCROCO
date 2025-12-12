# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.x] - 2026-xx-xx

### Added
 

### Fixed


- OCEAN : Fix unclosed parenthesis when TS_DIF4 is defined without DIF_COEF_3D (#482)

- COUPLING : missing mpi_cpl.h in get_grid.F in case of variable Z0 (Z0B_VAR) (#466)


### Changed

- BIOLOGY : Bug fix + add of PISCES diagnostics without XIOS (#474)

- OMEGA : Add a condition on the NBQ_MASS key for some terms of the first 
  part of the computation of omega (#447)

- CPP keys : restore the default definition for LIMIT_BSTRESS (the key is activated 
  in cppdefs_dev.h unless BSTRESS_FAST is previously defined) (#456) 

### Deprecated


### Removed


### Other



### Contributors on this release

- Contributors already on board : 
  G. Cambon, N. Ducousso, S. Le Gac, P. Marchesiello, R. Person

- New contributors : E. Le Bouedec
