# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.x.x] - 2024-xx-xx
### Added


### Fixed

- Fixes issue with OBC_COM_M2CHARACT_XXX (modified by mistake into OBC_COM_M2ORLANSKI in u2bc and v2bc) : fixes [#186]

- Fix the mismatch between PISCES code and input name variables for BSi, DSi and GSi : [#188]
  
- Fix test case ESTUARY grid with MPI decomposition [#193]

- The diffusion of wave numbers in coastal areas was resulting in instability 
  and strong currents (several m/s). Resolved by modifications in WKB_ADD_DIFF :
  fixes [#197] 

- Fix conversion error in ripple_dim subroutine [#200]
  
- Fix typo in analytical.F for MRL_WCI (Wave-current interactions) [#208]

- Fix the reading of initial state in case of MPINOLAND and different MPI 
  decomposition [#217]

- Fix typo in step.F on RVTK cppkey and fix reproductibility test 
  for MUSTANG [#222]

- Fix CPL scripts: 
     - when running in chained_job: remove test on deprecated variable, 
       select jobid on DATARMOR and JEAN-ZAY [#209]
     - when using online OCE_COMPIL with new python-generated ini files [#205]
     - for computing time steps when using nesting [#226]

### Changed

### Deprecated

### Removed

### Other
