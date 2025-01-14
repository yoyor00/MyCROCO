# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added

### Fixed


- Issue #252 : fix PSOURCE_MASS capabilities broken by previous change
- Issue #258 : fix misuse of temporary WFe,WFx arrays for horizontal w 
  advection in NBQ
- Issue #264 : fix wrong hbl test on XIOS field activation with GLS_MIXING 
- Issue #267 : fix several bugs with MUSTANG output and initialisation

### Changed

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

- Issue #141 : PISCES improvements and changes
  - Phasing of the PISCES version with that used for the CMIP7 exercise (NEMO 4.2.*)
  - Update on the PISCES interfacing module between NEMO and CROCO
  - Diagenetic module improvements: performance and diagenetic processes (e.g. increased number of POC classes, ...)
  - Added creation of an independent pisces restart file (managed in namelist_pisces_ref) to improve restartability
  - Rename the simplified version of PISCES, cpp key pisces_npzd
  - Correction of some bug fixes
  
### Deprecated

### Removed

- Issue #163 : remove cppkeys key_MUSTANG_specif_outputs and 
  key_MUSTANG_add_consol_outputs, MUSTANG outputs are now all 
  specified by namelist
- remove cppkeys key_CROCO and MORPHODYN_MUSTANG_byHYDRO

### Other
