# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [x.x.x] - xxxx-xx-xx

### Added

- Issue #250 : add dredging feature in MUSTANG
- Issues #291, #323 : add new test case FLASH_RIP

### Fixed

- Issue #172 : fix error in wrt_his.F when nrechis=0 and nrpfhis=1
- Issue #191 : fix compilation error with TVD and NBQ
- Issue #239 : back to previous default option in create_config
- Issue #252 : fix PSOURCE_MASS capabilities broken by previous change
- Issue #258 : fix misuse of temporary WFe,WFx arrays for horizontal w 
  advection in NBQ
- Issue #259 : remove unused file OCEAN/spkitlocal_nh.F90 
- Issue #260 : fix wrong array name for XIOS and WKB_WAVE
- Issue #261 : fix initialisation of non-declared arrays
- Issue #263 : cleaning unused variables in set_diags_ek.F and set_diags_pv.F
- Issue #264 : fix wrong hbl test on XIOS field activation with GLS_MIXING
- Issue #267 : fix several bugs with MUSTANG output, initialisation, add a 
  minimim porosity in deposit in V1 and fix transition with ero_option=3
- Issues #272, #279 and #321 : fix for ABL1D :
  - sea surface currents in ABL1D
  - wrong mixing length computation
  - ABL1D perfect restart
  - ECUME6/ABL1D
  - fix cppdefs organisation and incompatibility with update of ABL1D
- Issue #285 : fix NBQ+XIOS compilation
- Issue #299 : cleaning, comment in sta.h
- Issue #305 : fix MPI repro with PSOURCE when a source is 
  in one MPI domain rejecting towards another MPI domain
- Issue #308 : fix and cleaning with AGRIF and TIDE_MAS
- Issue #309 : fix parallel compilation (except for AGRIF)
- Issue #311 and #328 : fix data instruction replace by a declaration 
  of array and initialisation of it
- Issue #316 : fix MPI reproducibility with NBQ
- Issue #318 : fix typo in order 5 scheme
- Issue #319 : fix MPI reproducibility using MUSTANG_CORFLUX
- Issue #322 : fix SUBSTANCE_SUBMASSBALANCE feature on river fluxes
- Issue #324 : fix FILLVAL, mask was misdone on variable value 
  instead of mask value
- Issue #329 : fix writing of sediment layers (SEDIMENT)
- Issue #331 : fix rstTime netcdf index mismatch between PISCES and CROCO


### Changed

- Issue #330 : Change LOGFILE cppkey behavior by enabling to choose filename in
  croco.in

- Issue #306 : in MRL_WCI change limits of wave height in surfzone

- Issues #291, #323: change wavemaker (periodic boundaries + single-sum), 
  check of mean angle after the procedure and periodization is 
  abandoned if the mean angle change is too high

- Issue #281 : optimization of the PISCES code on
  - representation of the lability of the particle pool
  - several optimizations to the calculation of certain variables (performance).

- Issue #250 : for MUSTANG initialisation, review of namelist parameters

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

- Issue #271 : add information about the branch status (stable or not)
  in the README file
