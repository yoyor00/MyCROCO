# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.1] - 2025-09-26

### Added

- PISCES  : Integration of slight reformulation of zooplankton grazing 
  according to prey size done in PISCES standard version (#340) into quota version 
  (Issue #361)

- SCRIPTS: in run_croco_inter add TIDE_FILES and ONLINE (Freq, path) missing information 

### Fixed

- PISCES : Fix dummy line (kt variable) at end in p4zche.F90  (Issue #420)
- PISCES : Fixed error on diagnostic ligands and add Fe2+ oxydation rate (Issue #371)
- COUPLING : fixes to prevent runtime crash when compiled in full debug mode (Issue #376)
- COUPLING : patm2D was declared twice in case of OW_COUPLING and READ_PATM (Issue #383)
- COUPLING : undef TIDES_MAS in oce_compile.sh, as it does not work without USE_CALENDAR 
             (link with #233)
- COUPLING : update croco.in.base accordingly to original croco.in (solve #366)
- COUPLING : path to local myenv_mypath in job dir (solve #421)
- COUPLING : add "if" conditions for cpl restart file path definitions to avoid issues when
             one of the model was not requested in create_config (solve issue #428)
- COUPLING and XIOS : correct path and  requested copy of cppdefs_dev (solve #364)
- COUPLING and EXACT RESTART : manage EXACT_RESTART option in coupling scripts (solve #92)
- BENCH : do not put report status to True for reference variant to avoid
  to mark test passed even if not (Issue #342)
- BENCH : put jobcomp.log in results directory even if build fail (Issue #341)
- BENCH : remove openmp reproducibility check on SHOREFACE case (Incident #358)
- BENCH : fix typo AGRIF_2W to AGRIF_2WAY in realist and vortex json files (Issue #368)

- DIAGNOSTICS_EDDY & not XIOS : fix double comma in ncscrum.h (Issue #362)

- Compilation : fix cat "croco_ascii.txt" command in case of relative path (Issue #424)

- Wavemaker : Fix boundary forcing in case of eastern boundary (Issue #432)

- AGRIF : incompatibility of AGRIF with cppkeys
  BULK_ECUMEV0 or BULK_ECUMEV6 (Issue #422)
- AGRIF and psource_ncfile : in croco.in.1 file should be croco_runoff.nc.1 (solve #436)
- AGRIF: sponge keyword missing: put it back, path for online corrected in croco_inter.in,
  AGRIF_Fixed.in not copied from the right directory: corrected, copy croco_frc for all
  domains for tides (even if theoretically useless, it creates a read_inp error) (solve #438)

- OUTPUT : incompatibility when activating 
  PISCES + PSOURCE_NCFILE_TS + DIURNAL_INPUT_SRFLX (Issue #435)

### Changed

- SCRIPTS : create_config now allows to deal with 2 new options for preprocessing scripts: 
            pytools and mattools (instead of prepro) - solve issue #295 
            and 2 distinct options for run scripts (either inter for Plurimonth_scripts type: 
            run_croco_inter, or runcpl for SCRIPTS_COUPLING type: submitjob, mynamelist...)

### Deprecated


### Removed


### Other


### Contributors on this release

- Contributors already on board : 
  R. Benshila, M. Caillaud, G.Cambon, S. Jullien, S. Le Gac, 
  P. Marchesiello, C. Nguyen, R. Person, J. Pianezze, S. Treillou

- New contributors : 
  M. Schreiber, A. Zribi  
