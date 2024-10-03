# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx
### Added

### Fixed

- Fixes issue with OBC_COM_M2CHARACT_XXX (modified by mistake into 
  OBC_COM_M2ORLANSKI in u2bc and v2bc) : fixes [#186]

- Fix the mismatch between PISCES code and input name variables for BSi, 
  DSi and GSi : [#188]

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

- Issue #184 : change roller contribution on Stokes Drift in very shallow water 
  to avoid instabilities

- Issue #163 : for MUSTANG output, 
  - avoid possibility of overlapping in vname by 
    using a separate array vname_must
  - use l_out_subs from substance namelist to allow the output of only wanted 
    substance
  - adding boolean in namelist for choosing which variables to output and 
    allocate only the needed arrays
  - update paraMUSTANG_defaults.txt with new booleans availables
  - update XIOS output to have the save available variables in all MUSTANG
    output options

### Deprecated

### Removed

- Issue #163 : remove cppkeys key_MUSTANG_specif_outputs and 
  key_MUSTANG_add_consol_outputs, MUSTANG outputs are now all 
  specified by namelist
- remove cppkeys key_CROCO and MORPHODYN_MUSTANG_byHYDRO

### Other
