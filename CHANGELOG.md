# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.x.x] - 2024-xx-xx
### Added


### Fixed

- Fixes issue with OBC_COM_M2CHARACT_XXX (modified by mistake into OBC_COM_M2ORLANSKI in u2bc and v2bc) : fixes [#186](https://gitlab.inria.fr/croco-ocean/croco/-/issues/186)

- Fix the mismatch between PISCES code and input name variables for BSi, DSi and GSi : [#188](https://gitlab.inria.fr/croco-ocean/croco/-/issues/188)
  
- The diffusion of wave numbers in coastal areas was resulting in instability 
  and strong currents (several m/s). Resolved by modifications in WKB_ADD_DIFF :
  fixes [#197] 

### Changed

### Deprecated

### Removed

### Other
