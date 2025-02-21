# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added

- Issue #298 : Benchmarking python tools to perform tests runs, see documentation or 
  BENCH/README.md for more details
- Issue #291 : add new test case FLASH_RIP

### Fixed

- Issue #239 : back to previous default option in create_config
- Issue #252 : fix PSOURCE_MASS capabilities broken by previous change
- Issue #258 : fix misuse of temporary WFe,WFx arrays for horizontal w 
  advection in NBQ
- Issue #259 : remove unused file OCEAN/spkitlocal_nh.F90 
- Issue #260 : fix wrong array name for XIOS and WKB_WAVE
- Issue #263 : cleaning unused variables in set_diags_ek.F and set_diags_pv.F
- Issue #264 : fix wrong hbl test on XIOS field activation with GLS_MIXING
- Issue #267 : fix several bugs with MUSTANG output, initialisation, add a 
  minimim porosity in deposit in V1 and fix transition with ero_option=3
- Issue #272 : fix CPP key incompatible for ABL1D
- Issue #279 : fix for :
  - sea surface currents in ABL1D
  - wrong mixing length computation
  - ABL1D perfect restart
  - ECUME6/ABL1D
- Issue #285 : fix NBQ+XIOS compilation
- Issue #299 : cleaning, comment in sta.h
- Issue #305 : fix MPI repro with PSOURCE when a source is 
  in one MPI domain rejecting towards another MPI domain
- Issue #309 : fix parallel compilation (except for AGRIF)
- Issue #318 : fix typo in order 5 scheme
- Issue #319 : fix MPI reproducibility using MUSTANG_CORFLUX

### Changed

- Issue #298 : Jobcomp script can now be executed with terminal options, see 
  ```./jobcomp -h``` for more details

- Issue #306 : in MRL_WCI change limits of wave height in surfzone

- Issue #291 : change wavemaker (periodic boundaries + single-sum)

- Issue #281 : optimization of the PISCES code on
  - representation of the lability of the particle pool
  - several optimizations to the calculation of certain variables (performance).

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
  - Correction of some bugs
  
### Deprecated

### Removed

- Issue #163 : remove cppkeys key_MUSTANG_specif_outputs and 
  key_MUSTANG_add_consol_outputs, MUSTANG outputs are now all 
  specified by namelist
- remove cppkeys key_CROCO and MORPHODYN_MUSTANG_byHYDRO

### Other
