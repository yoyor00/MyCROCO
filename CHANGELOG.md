# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.0] - 2025-03-20

### Added

- NEW TEST CASE : FLASH_RIP (Issues #291 and #323)

- NEW DIAGNOSTIC : add mixed layer diagnostics computation (criterion : 
  density, temperature, bvf) (Issue #218)

- NEW CPPKEY : add RAIN_FLUX cppkey to take into account water and temp flux 
  from rain (Issue #195)

- MUSTANG : add dredging feature (Issue #250)

- COUPLED RUN SCRIPTS
  - adding possibility to give a pre-built grid file to OASIS
  - new to take into account the ideal WRF configurations
  - possibility to have different WRF output frequencies between domains

### Fixed

- IO
  - fix error in wrt_his.F when nrechis=0 and nrpfhis=1 (Issue #172)
  - fix wrong array name for XIOS and WKB_WAVE (Issue #260)
  - fix wrong hbl test on XIOS field activation with GLS_MIXING (Issue #264)
  - fix FILLVAL, mask was misdone on variable value 
    instead of mask value (Issue #324)
  - fix writing of sediment layers (SEDIMENT) (Issue #329)
  - fix rstTime netcdf index mismatch between PISCES and CROCO (Issue #331)

- COMPILATION
  - fix compilation error with TVD and NBQ (Issue #191)
  - fix NBQ+XIOS compilation (Issue #285)
  - fix parallel compilation (except for AGRIF) (Issue #309)

- COUPLED RUN SCRIPTS
  - back to previous default option in create_config interactive q/a (Issue #239)

- RIVER
  - fix PSOURCE_MASS capabilities broken by previous change (Issue #252)
  - fix MPI reproducibility with PSOURCE when a source is 
    in one MPI domain rejecting towards another MPI domain (Issue #305)

- NUMERIC
  - fix misuse of temporary WFe,WFx arrays for horizontal w 
    advection in NBQ (Issue #258)
  - fix typo in order 5 scheme (Issue #318)

- NBQ
  - fix MPI reproducibility with NBQ (Issue #316)

- ABL1D (Issues #272, #279 and #321)
  - sea surface currents in ABL1D
  - wrong mixing length computation
  - ABL1D perfect restart
  - ECUME6/ABL1D
  - fix cppdefs organisation and incompatibility with update of ABL1D

- MUSTANG
  - fix several bugs with MUSTANG output, initialisation, add a 
    minimim porosity in deposit in V1 and fix transition with 
    ero_option=3 (Issue #267)
  - fix MPI reproducibility using MUSTANG_CORFLUX (Issue #319)
  - fix SUBSTANCE_SUBMASSBALANCE feature on river fluxes (Issue #322)

- OTHERS
  - fix initialisation of non-declared arrays (Issue #261)
  - cleaning unused variables in set_diags_ek.F and set_diags_pv.F (Issue #263)
  - cleaning, comment in sta.h (Issue #299)
  - fix and cleaning with AGRIF and TIDE_MAS (Issue #308)
  - fix data instruction replace by a declaration 
    of array and initialisation of it (Issue #311 and #328)
  - fix step3d_fast with wave but no wet and dry (Issue #332)
  - fix wrongly place mask in u2bc (Issue #338)

### Changed

- WAVE
  - in MRL_WCI change limits of wave height in surfzone (Issue #306)
  - change wavemaker (periodic boundaries + single-sum), 
    check of mean angle after the procedure and periodization is 
    abandoned if the mean angle change is too high (Issues #291 and #323)

- MUSTANG 
  - for MUSTANG output (Issue #163)
    - avoid possibility of overlapping in vname by 
      using a separate array vname_must
    - use l_out_subs from substance namelist to allow the output of only wanted 
      substance
    - adding boolean in namelist for choosing which variables to output and 
      allocate only the needed arrays
    - update paraMUSTANG_defaults.txt with new booleans availables
    - update XIOS output to have the save available variables in all MUSTANG
      output options
  - for MUSTANG initialisation, review of namelist parameters (Issue #250)

- PISCES improvements and changes (Issue #141 and Issue #281 )
  - Phasing of the PISCES version with that used for the CMIP7 exercise (NEMO 4.2.*)
  - Update on the PISCES interfacing module between NEMO and CROCO
  - Diagenetic module improvements: performance and diagenetic processes (e.g. increased number of POC classes, ...)
  - Added creation of an independent pisces restart file (managed in namelist_pisces_ref) to improve restartability
  - Rename the simplified version of PISCES, cpp key pisces_npzd
  - Correction of some bugs
  - Optimization of the PISCES code on
    - representation of the lability of the particle pool
    - several optimizations to the calculation of certain variables (performance).

- COUPLED RUN SCRIPTS
  - improving possibility to add weight file for OASIS
  - removing WRF rsl files because too heavy 
  - modify atm_his_h in hours as atm_his in minutes
  - update for suporting both WRFV4.2.1 and WRFV4.6 OASIS names
  - update create_config for including pytools

### Deprecated

### Removed

- Removed CPPKEYS :
  - key_MUSTANG_specif_outputs, key_MUSTANG_add_consol_outputs 
    (see Issue #163 MUSTANG outputs are now all specified by namelist)
  - key_CROCO, MORPHODYN_MUSTANG_byHYDRO

- Remove unused file OCEAN/spkitlocal_nh.F90 (Issue #259)

### Other

- Add information about the branch status (stable or not)
  in the README file (Issue #271)
