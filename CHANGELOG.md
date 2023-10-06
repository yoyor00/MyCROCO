# Changelog

## [Unreleased] - xxxx-xx-xx

### Added

- ABL1D : 
  - TODO add description
  - Add KILPATRICK test-case

- SUBMASSBALANCE : Add submassbalance feature as in MARS, see issue 
  [#65](https://gitlab.inria.fr/croco-ocean/croco/-/issues/65) ; 
  check merge request if interested

- PISCES :
  - TODO add description
  - Update PISCES (XIOS on-the-fly, ...)

- SED_DENS : Add suspended sediment contribution to density whith cpp key 
  #SED_DENS. #SED_DENS applies to #SEDIMENT and #MUSTANG cases. For MUSTANG 
  case, variable using #key_sand2D are exclude of the contribution. See 
  issue [#124](https://gitlab.inria.fr/croco-ocean/croco/-/issues/124)

### Fixed

- TKE : Changing the coefficient 'invG' used in the tridiagonal solving of 
  the tke and gls variables in gls_mixing.F. Using the Patankar trick in the 
  tridiagonal solving consists of moving the Bprod term from the forcing 
  vector S to the diagonal D when the total forcing Sprod + Bprod is 
  negative. Doing this change needs to add to the term Bprod a 'minus' sign 
  and a division by the variable (tke for the tke loop and psi for the gls 
  loop). This was done for the gls loop (invG = 1/psi in this case) but not 
  for the tke loop (ingG = 1). Hence, invG is set to 1/tke in the tke loop.

- EDDY_DIAGS : fixed see issues: 
  [#107](https://gitlab.inria.fr/croco-ocean/croco/-/issues/107), 
  [#117](https://gitlab.inria.fr/croco-ocean/croco/-/issues/117)

- KPP when ANA_DIURNAL_SW is activated : see issue 
  [#115](https://gitlab.inria.fr/croco-ocean/croco/-/issues/115)

- MRL_WCI : correct references to Qiao et al. non-breaking wave diffusion
  (thanks to ChangShui Xia )

- Run issue of RIP and SANDBAR test cases : see issue
  [#130](https://gitlab.inria.fr/croco-ocean/croco/-/issues/130)

### Changed

- CI runs through docker container from croco-ocean/croco_docker project

### Deprecated

### Removed

### Other

## [v.1.3.1](https://gitlab.inria.fr/croco-ocean/croco/-/releases/v1.3.1) - 2023-09-21

### Added

- COUPLING
  - Add EXACT_RESTART + OA/OW_COUPLING capability : the coupler initialisation 
    cpl_prism_define.F was wrongly done twice in this case (line 1147 in 
    get_initial.F). It was leading to a coupler crash during the initilization 
    phase.
  - Add the capability to cutoff the heat fluxes in case of water under the 
    freezing point (-1.8°C) => Simplest sea ice correction to T,S fluxes.

### Fixed

- OW_COUPLING: missing if defined OW_COUPLING in forces.h for smstr (that can 
  be updated in case of OW COUPLING). Fix issue #127

- USE_CALENDAR : 
  - Compatibility with XIOS
  - Change format of date format : 
    - in croco.in : from dd/mm/yyyy to yyyy-mm-dd
    - read origin_date in time netcdf attribute for netCDF input files.

- DIAGNOSTICS_BIO : 
  - fix netcdf parralel (NC4PAR) netcdf files writing with DIAGNOSTICS_BIO 
    (bgc fluxes)
  - fix time and scrum_time in avergae files of croco_bio_diags_avg.nc

- WET_DRY : solve blowup problem with dry cells at boundaries

- AGRIF Compilation : fix jobcomp.bash with gcc > 10

- PISCES : fix missing local variables initialization

### Changed

### Deprecated

### Removed

### Other

## [v.1.3] - 2022-11-28
