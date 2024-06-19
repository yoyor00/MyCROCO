# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx
### Added

### Fixed

- Fixes issue with OBC_COM_M2CHARACT_XXX (modified by mistake into OBC_COM_M2ORLANSKI in u2bc and v2bc) : fixes [#186]((https://gitlab.inria.fr/croco-ocean/croco/-/issues/186)

### Changed

- Issue #184 : change roller contribution on Stokes Drift in very shallow water to avoid instabilities
- Issue #163 : for MUSTANG output, 
  - avoid possibility of overlapping in vname by 
    using a separate array vname_must
  - use l_out_subs from substance namelist to allow the output of only wanted 
    substance
  - adding boolean in namelist for choosing which variables to output and 
    allocate only the needed arrays

### Deprecated

### Removed

- Issue #163 : remove cppkeys key_MUSTANG_specif_outputs and 
  key_MUSTANG_add_consol_outputs, MUSTANG outputs are now all 
  specified by namelist

### Other
