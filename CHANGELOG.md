# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.0.1] - 2024-10-04

### Added

### Fixed

- Fixes issue with OBC_COM_M2CHARACT_XXX (modified by mistake into 
  OBC_COM_M2ORLANSKI in u2bc and v2bc) : fixes [#186]

- Fix the mismatch between PISCES code and input name variables 
  for BSi, DSi and GSi : [#188] [#213]
  
- Fix test case ESTUARY grid with MPI decomposition [#193]

- The diffusion of wave numbers in coastal areas was resulting in instability 
  and strong currents (several m/s). Resolved by modifications in WKB_ADD_DIFF :
  fixes [#197] 

- Fix conversion error in ripple_dim subroutine [#200]
  
- Fix typo in analytical.F for MRL_WCI (Wave-current interactions) [#208]

- Fix the reading of initial state of MUSTANG in case of MPINOLAND 
  and different MPI decomposition [#217]

- Fix typo in step.F on RVTK cppkey and fix reproductibility test 
  for MUSTANG [#222]

- Fix CPL scripts: 
     - create_config.bash when using without MATLAB tools [#204]   
     - when running in chained_job: remove test on deprecated variable, 
       select jobid on DATARMOR and JEAN-ZAY [#209]
     - when using online OCE_COMPIL with new python-generated ini files [#205]
     - when using blk (not online) and nesting [#232]
     - for computing time steps when using nesting [#226]

- Fix various typo and f77 fixed-form file with line length > 72 
  [#169] [#206] [#216] [#225] [#228]

- Fix exact restart not tested for SFLUX_CFB and exact restart not achieve 
  with PSOURCE [#219]

- Fix DIAGNOSTICS_EK_FULL compatibility with XIOS (in send_xios_diags.F) [#220]

### Changed

- Change roller contribution on Stokes Drift in very shallow water 
  to avoid instabilities [#184]

- Avoid code repetition in CI for gfortran/ifort [#201]

### Deprecated

### Removed

### Other
