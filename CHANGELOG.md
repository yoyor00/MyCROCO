# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.3] - 2026-04-07

### Added
 

### Fixed


OCEAN :
- Fix unclosed parenthesis when TS_DIF4 is defined without DIF_COEF_3D (#482)

AGRIF : 
- Fix allocation of message passing arrays (ibuf...) when 3 ghost points
  needed (UP5, WENO) (#310 and #458)

COUPLING : 
- Fix missing mpi_cpl.h in get_grid.F in case of variable Z0 (Z0B_VAR) (#466)

PSOURCE_NCFILE : 
- Make it usable with NO_TRACER (#459)

DIAGNOSTICS : 
- Fix integration error in diagnostics on MLD in (density, temp, n2) (#429) 
- Keep only density criterion (#429)

WET_DRY : 
- Add the correct masking of grid stiffness ratios rx0 and rx1 (#373)

### Changed

BIOLOGY : 
- Bug fix (#474)
- Add of PISCES diagnostics without XIOS (#474)

CPP keys :
- Restore default definition of `LIMIT_BSTRESS`
  (enabled in `cppdefs_dev.h` unless `BSTRESS_FAST` is defined) (#456) 

### Deprecated


### Removed


### Other


### Contributors on this release

- Contributors already on board : 
  G. Cambon, N. Ducousso, S. Jullien, M. Le Corre, S. Le Gac, P. Marchesiello, R. Person, 
  A.-L. Schaefer

- New contributors : E. Le Bouedec
