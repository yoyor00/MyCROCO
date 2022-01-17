CROCO v1.2
-----------
Released date : 18 January 2022

Previous release : CROCO v1.1 (October 2019)


## Environment

The whole structure of the repository has been updated, please read carefully. It comes with a new set of scripts to help to launch production simulations. See <https://croco-ocean.gitlabpages.inria.fr/croco_doc/tutos/tutos.02.contents.html> for details.

- create_config.bash: script to create a configuration environment. 2 typical modes are proposed: all-dev for the usual ("all-in") architecture for forced CROCO runs and/or dev, and all-prod/all-prod-cpl for production run architecture and coupling with external models
- SCRIPTS/ directory: where scripts for running croco in production mode are provided
- croco input and output files for interannual simulations can now be named systematically with 2 digits for month: croco_..._Y????M??.nc instead of croco_..._Y????M?.nc (the old option is however still available)
- croco/XIOS/process_xios_xml.sh : stand-alone cpp-process of the XIOS xml files. It is now separated from the jobcomp script. It is deployed in _Run/_ in all-dev setup, and in _CONFIGS/PREPRO/XIOS/_ in all-prod setup

## Air-sea interactions

_full description :_ [_https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.forcing.html_](https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.forcing.html)

- 1D Atmospheric Boundary Layer model (ABL)
- Updated Bulks
- Updated current feedback

## Coupling

_full description :_ [_https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.coupling.html_](https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.coupling.html) and [_https://croco-ocean.gitlabpages.inria.fr/croco_doc/tutos/tutos.16.coupling.html_](https://croco-ocean.gitlabpages.inria.fr/croco_doc/tutos/tutos.16.coupling.html)

- SCRIPTS/SCRIPTS_COUPLING: new coupling toolbox for running CROCO interannual simulations coupled with atmosphere and wave models. For more details see <https://croco-ocean.gitlabpages.inria.fr/croco_doc/tutos/tutos.16.coupling.workflow.html>
- Now manage nesting in the ocean and in the atmosphere
- Additional wave parameters can now be exchanged

## Sediment

_full description :_ [_https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.modules.sediment.html_](https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.modules.sediment.html)

- MUSTANG sediment model : MUSTANG is a sediment model developed at Ifremer and proposed in addition to the existing one from USGS
- Bedload : improvement of robustness and new parametrization available
- Cohesive bed : introduction of bed stratigraphy with cohesive or mixed cohesive/non-cohesive case

## Numerics

_full description :_ [_https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.numerics.html_](https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.numerics.html)

- quasi-hydro CROCO_QH : a quasi-hydrostatic option (non-traditional Coriolis terms)
- MORPHODYN correction : improvement of robustness for morphodynamics
- NO_TEMPERATURE and NO_TRACERS option : use only salinity as active tracer, run in fully homogeneous mode

## New Configurations

_full description :_ [_https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.test_cases.html_](https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.test_cases.html)

- 3 new test cases for sediment processes (can be used with USGS or MUSTANG sediment model)
  * TIDAL_FLAT : a 2DV tidal flat with sediment mixture
  * SED_TOY : a single column case
  * DUNE : migration of a dune with different sand classes
- 1 new realistic configuration COASTAL : a coastal configuration with a finer resolution than BENGUELA and corresponding/different settings (rivers, GLS mixing ...)

## Miscellaneous

- Multiple passive tracers : the number of passive tracers is now a parameter up to 100
- XIOS : former XIOS2 option is now the only option
- Restart
- NC4 for all diagnostics : all diagnostics are compatible with NetCDF4 parallel writing
- NO_LAND and PARALLEL_FILES compatibility
- EXACT_RESTART : bit to bit restartability

## Detail list of CPP keys

Suppressed keys : ''ANA_BPFLUX', 'ANA_SPFLUX', 'BEDLOAD_SOULSBY', 'BULK_EP', 'BULK_FAIRALL', 'BULK_SMFLUX', 'BULK_SM_UPDATE', 'CFB_STRESS2', 'CFB_WIND', 'MLCONVEC', 'PLUME', 'SMFLUX_CFB', 'STFLUX_CFB', 'WAVE_BODY_STREAMING', 'XIOS2'

New keys (extented version...) : 'ANA_DUNE', 'ANA_MORPHODYN', 'BEDLOAD_MARIEU', 'BEDLOAD_UP1', 'BEDLOAD_UP5', 'BEDLOAD_VANDERA', 'BEDLOAD_WENO5', 'BEDLOAD_WULIN', 'BED_ARMOR', 'BED_HIDEXP', **~~'BHFLUX'~~,** 'BULK_ECUMEV0', 'BULK_ECUMEV6', 'BULK_GUSTINESS', 'BULK_WASP','CFB_WIND_TRA', 'COASTAL', 'COHESIVE_BED', 'CROCO_QH', 'DIAGNOSTICS_BARO', 'DIAGNOSTICS_TSVAR', 'DO_NOT_OVERWRITE', 'DUNE', 'DUNE3D', 'ECUMEv0', 'ECUMEv6', 'EXACT_RESTART', 'FLOC_BBL_DISS', 'FLOC_TURB_DISS', 'GUSTINESS', 'HOURLY_VELOCITIES', 'ISOLITON', 'LMD_LANGMUIR', 'M3FAST_REINIT', 'MIXED_BED', 'MORPHODYN_MUSTANG_byHYDRO', 'MPI_TIME', 'MUSTANG', 'MUSTANG_CORFLUX', 'NO_TEMPERATURE', 'NO_TRACER', 'RVTK_DEBUG', 'RVTK_DEBUG_ADVANCED', 'RVTK_DEBUG_PERFRST', 'SANDBAR_OFFSHORE', 'SANDBAR_ONSHORE', 'SED_DEFLOC', 'SED_FLOCS', 'SED_TAU_CD_CONST', 'SED_TOY', 'SED_TOY_CONSOLID', 'SED_TOY_FLOC', 'SED_TOY_RESUSP', 'SED_TOY_ROUSE', 'SFLUX_CFB', 'SLOPE_KIRWAN', 'START_DATE', 'STOKES_DRIFT', 'SUBSTANCE', 'TAU_CRIT_WULIN', 'TEMPERATURE', 'TIDAL_FLAT', 'TRACERS', 'VILAINE', 'WASP', 'WAVE_BREAK_BJ78', 'ZETA_DRY_IO', 'key_CROCO', 'key_MUSTANG_V2', 'key_MUSTANG_bedload', 'key_MUSTANG_specif_outputs', 'key_noTSdiss_insed', 'key_nofluxwat_IWS', 'key_sand2D', 'key_tenfon_upwind'


 _   _                      __                 _
| | | | __ ___   _____     / _|_   _ _ __     | |
| |_| |/ _` \ \ / / _ \   | |_| | | | '_ \    | |
|  _  | (_| |\ V /  __/   |  _| |_| | | | |   |_|
|_| |_|\__,_| \_/ \___|   |_|  \__,_|_| |_|   (_)
