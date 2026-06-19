!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

!=======================================================================
!  MODULE croco_namelist
!
!  Purpose : Declaration and storage of all namelist parameters.
!            Default values are set here.
!=======================================================================

MODULE croco_namelist
   implicit none
   save
   public

   ! namelist filename (set via read_nml_fname from command-line arg 2)
   character(len=200) :: fname_nml = 'croco.nml'

   ! &croco_title
   character(len=80) :: title = "CROCO simulation"
   !! Configuration name

   ! &croco_logfile
   character(len=180) :: logname = "croco.log"
   !! File path for STDOUT (used with the `LOGFILE` CPP option)

   ! &croco_time_stepping
   real    :: dt = 0.0
   !! Baroclinic time step [s]
   integer :: ntimes = 0
   !! Number of time-steps required for the simulation
   integer :: ndtfast = 20
   !! Number of barotropic time-steps between each baroclinic time step.
   !! For 2D configurations, ndtfast should be unity.
   integer :: ninfo = 1
   !! Number of time-steps between printing of information to standard output

   ! &croco_history
   logical :: ldefhis = .true.
   !! Flag (T/F) for writing history files.
   !! If `.true.`, a new history file is created.
   !! If `.false.`, data is appended to an existing history file.
   integer :: nwrt = 72
   !! Number of time-steps between the writing of history fields
   integer :: nrpfhis = 0
   !! History file cycling switch:
   !!
   !!  - 0: write several records every NWRT time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFHIS records per file
   !!  - -1: overwrite record every NWRT time steps
   character(len=180) :: hisname = "CROCO_FILES/croco_his.nc"
   !! Name of history file

   ! &croco_initial
   integer :: nrrec = 1
   !! Switch to indicate start or restart from a previous solution.
   !! `nrrec` is the time index of the initial or restart NetCDF file
   !! assigned for initialization. If `nrrec` is negative (e.g. `nrrec = -1`),
   !! the model will start from the most recent time record available
   !! in the file (the initialization record is assigned internally).
   character(len=180) :: ininame = "CROCO_FILES/croco_ini.nc"
   !! Name of file containing the initial state

   ! &croco_restart
   integer :: nrst = 720
   !! Number of time-steps between writing of restart fields
   integer :: nrpfrst = -1
   !! Restart file cycling switch:
   !!
   !!  - 0: write several records every NRST time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFRST records per file
   !!  - -1: overwrite record every NRST time steps
   character(len=180) :: rstname = "CROCO_FILES/croco_rst.nc"
   !! Name of restart file

#ifdef NBQ
   ! &croco_time_stepping_nbq
   integer :: ndtnbq = 1
   !! Number of non-hydrostatic (NBQ) sub-steps per baroclinic time step
   real    :: csound_nbq = 1000.0
   !! Speed of sound used in the NBQ solver [m/s]
   real    :: visc2_nbq = 0.01
   !! Laplacian viscosity coefficient for the NBQ pressure solver [m2/s]
#endif

#ifdef SOLVE3D
   ! &croco_s_coord
   real :: theta_s = 7.0d0
   !! S-coordinate surface control parameter (stretching toward surface when > 0)
   real :: theta_b = 2.0d0
   !! S-coordinate bottom control parameter (stretching toward bottom when > 0)
   real :: Tcline = 200.0d0
   !! Width of the surface or bottom boundary layer [m] in which higher vertical
   !! resolution is required during S-coordinate stretching
#endif

#ifdef USE_CALENDAR
   ! &croco_use_calendar
   character(len=19) :: start_date = '2000-01-01 00:00:00'
   !! Run start date, format `YYYY-MM-DD HH:MM:SS` (used with `USE_CALENDAR`)
   character(len=19) :: end_date = '2000-02-01 00:00:00'
   !! Run end date, format `YYYY-MM-DD HH:MM:SS` (used with `USE_CALENDAR`)
   real :: dt_his = 1.0
   !! Time interval between history output records [hours] (used with `USE_CALENDAR`)
   real :: dt_avg = 6.0
   !! Time interval between averages output records [hours] (used with `USE_CALENDAR`)
   real :: dt_rst = 12.0
   !! Time interval between restart output records [hours] (used with `USE_CALENDAR`)
#endif

#ifndef ANA_GRID
   ! &croco_grid
   character(len=180) :: grdname = "CROCO_FILES/croco_grd.nc"
   !! Grid filename
#endif

   ! &croco_forcing
   character(len=180) :: frcname = "CROCO_FILES/croco_frc.nc"
   !! Forcing filename
   ! TODO : clean this mess of cpp keys
#if !defined(NO_FRCFILE) && ( \
   defined(TIDES) ||\
   (defined(MRL_WCI) &&\
   !defined(ANA_WWAVE) && \
   !defined(WKB_WWAVE) && \
   !defined(OW_COUPLING) ) || \
   (defined(SFLX_CORR) &&\
   !defined(ANA_SSS) && \
   !defined(BULK_FLUX) && \
   !defined(OA_COUPLING) ) || \
   (defined(QCORRECTION) &&\
   !defined(ANA_SST) && \
   !defined(OA_COUPLING) ) || \
   (defined(BBL) &&\
   !defined(ANA_BSEDIM) && \
   !defined(SEDIMENT) ) || \
   ( !defined(ANA_STFLUX) && \
   !defined(BULK_FLUX) && \
   !defined(OA_COUPLING) && \
   defined(SOLVE3D)) ||\
   (defined(SALINITY) &&\
   !defined(ANA_SSFLUX) && \
   !defined(BULK_FLUX) && \
   !defined(OA_COUPLING) && \
   defined(SOLVE3D)) ||\
   ( !defined(ANA_SRFLUX) && \
   !defined(BULK_FLUX) && \
   !defined(OA_COUPLING) && \
   defined(SOLVE3D)) ||\
   ( !defined(ANA_SMFLUX) && \
   !defined(BULK_FLUX) && \
   !defined(OA_COUPLING) ) || \
   (defined(BULK_FLUX) &&\
   !defined(SALINITY) ) || \
   (defined(WAVE_OFFLINE) &&\
   !defined(MUSTANG) ) )
   logical :: use_frcname = .true.
#else
   logical :: use_frcname = .false.
#endif

   ! &croco_bulk_forcing
   character(len=180) :: bulkname = "CROCO_FILES/croco_blk.nc"
   !! Bulk forcing filename (used with `BULK_FLUX`)

#if (  defined TCLIMATOLOGY  && \
   !defined ANA_TCLIMA) || \
   (defined ZCLIMATOLOGY&&\
   !defined ANA_SSH) || \
   (defined M2CLIMATOLOGY&&\
   !defined ANA_M2CLIMA) || \
   (defined M3CLIMATOLOGY&&\
   !defined ANA_M3CLIMA)
   ! &croco_climatology
   character(len=180) :: clmname = "CROCO_FILES/croco_clm.nc"
   !! Climatology filename for boundary and nudging (used with `CLIMATOLOGY`)
#endif

#if !defined ANA_BRY && defined FRC_BRY
   ! &croco_boundary
   character(len=180) :: bry_file = "CROCO_FILES/croco_bry.nc"
   !! Boundary conditions filename (used with `FRC_BRY`)
#endif

#if defined WKB_WWAVE && !defined ANA_BRY_WKB
   ! &croco_wkb_boundary
   character(len=180) :: brywkb_file = "CROCO_FILES/croco_wkb.nc"
   !! WKB wave boundary conditions filename (used with `WKB_WWAVE`)
#endif

#ifdef WKB_WWAVE
   ! &croco_wkb_wwave  – primary waves and empirical breaking model parameters
   real :: wkb_amp = 0.25
   !! Offshore wave amplitude [m]
   real :: wkb_ang = 190.0
   !! Offshore wave angle [deg]
   real :: wkb_prd = 8.0
   !! Offshore wave period [s]
   real :: wkb_tide = -2.0
   !! Constant offshore water level [m]
   real :: wkb_btg = 1.3
   !! B parameter for the empirical wave breaking model
   real :: wkb_gam = 0.38
   !! Gamma parameter: Hrms/h ratio threshold for wave breaking
#  ifdef WAVE_ROLLER
   ! &croco_wkb_roller  – Svendsen (1984) surface roller model parameters
   real :: wkb_rsb = 0.1
   !! sin(beta): roller dissipation slope parameter
   real :: wkb_roller = 0.5
   !! Fraction of breaking energy transferred to the surface roller [0, 1]
#  endif
#endif

#ifdef WAVE_MAKER
   ! &croco_wave_maker  – offshore wavemaker parameters for wave-resolving simulations
   real :: wmaker_amp = 0.0
   !! RMS wave amplitude [m]
   real :: wmaker_prd = 8.0
   !! Peak wave period [s]
   real :: wmaker_dir = 0.0
   !! Mean wave angle [deg]
   real :: wmaker_dsp = 0.0
   !! Directional spread [deg]
   real :: wmaker_fsp = 3.3
   !! Frequency spread: gamma parameter in JONSWAP spectrum
#endif

#ifdef AVERAGES
   ! &croco_averages
   integer :: ntsavg = 1
   !! Starting time-step for the accumulation of output time-averaged data.
   !! For instance, set to start averaging over the last day of a 30-day run.
   integer :: navg = 48
   !! Number of time-steps between writing of time-averaged fields
   integer :: nrpfavg = 0
   !! Averages file cycling switch:
   !!
   !!  - 0: write several records every NAVG time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFAVG records per file
   !!  - -1: overwrite record every NAVG time steps
   character(len=180) :: avgname = "CROCO_FILES/croco_avg.nc"
   !! Name of averages file
#endif

#if defined OUTPUTS_SURFACE && !defined XIOS
   ! &croco_surf
   logical :: ldefsurf = .true.
   !! Flag (T/F) for writing surface output files.
   !! If `.true.`, a new file is created; if `.false.`, data is appended.
   integer :: nwrtsurf = 0
   !! Number of time-steps between writing of surface fields.
   !! 0: use the same frequency as history (nwrt).
   integer :: nrpfsurf = 0
   !! Surface output file cycling switch:
   !!
   !!  - 0: write several records every NWRTSURF time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFSURF records per file
   !!  - -1: overwrite record every NWRTSURF time steps
   character(len=180) :: surfname = "CROCO_FILES/croco_surf.nc"
   !! Name of surface output file
   ! &croco_surf_history_fields
   logical :: out_his_surf = .true.
   !! Flag (T/F) to write the surface output fields to the surface history file
#  ifdef AVERAGES
   ! &croco_surf_avg
   logical :: ldefsurf_avg = .true.
   !! Flag (T/F) for writing surface averaged output files.
   integer :: ntssurf_avg = 1
   !! Starting time-step for the accumulation of surface time-averaged data
   integer :: nwrtsurf_avg = 0
   !! Number of time-steps between writing of averaged surface fields.
   !! 0: use the same frequency as averages (navg).
   integer :: nrpfsurf_avg = 0
   !! Surface averages file cycling switch:
   !!
   !!  - 0: write several records every NWRTSURF_AVG time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFSURF_AVG records per file
   !!  - -1: overwrite record every NWRTSURF_AVG time steps
   character(len=180) :: surfname_avg = "CROCO_FILES/croco_surf_avg.nc"
   !! Name of surface averages file
   ! &croco_surf_average_fields
   logical :: out_avg_surf = .true.
   !! Flag (T/F) to write the surface output fields to the surface averages file
#  endif
#endif

#if defined ABL1D && !defined XIOS
   ! &croco_abl
   logical :: ldefablhis = .true.
   !! Flag (T/F) for creating the ABL history file.
   !! If `.true.`, a new file is created; if `.false.`, data is appended.
   integer :: nwrtablhis = 36
   !! Number of time-steps between writing of ABL history fields
   integer :: nrpfablhis = 0
   !! ABL history file cycling switch:
   !!
   !!  - 0: write several records every NWRTABLHIS time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFABLHIS records per file
   !!  - -1: overwrite record every NWRTABLHIS time steps
   character(len=180) :: ablname = "croco_abl_his.nc"
   !! Name of ABL history file
#  ifdef AVERAGES
   ! &croco_abl_averages
   logical :: ldefablavg = .false.
   !! Flag (T/F) for writing time-averaged ABL fields
   integer :: ntsablavg = 1
   !! Starting time-step for the accumulation of ABL time-averaged data
   integer :: nwrtablavg = 0
   !! Number of time-steps between writing of ABL averaged fields
   integer :: nrpfablavg = 0
   !! ABL averages file cycling switch:
   !!
   !!  - 0: write several records every NWRTABLAVG time steps
   !!  - >0: create more than one file (with sequential numbers) and write NRPFABLAVG records per file
   !!  - -1: overwrite record every NWRTABLAVG time steps
   character(len=180) :: ablname_avg = "croco_abl_avg.nc"
   !! Name of ABL averages file
#  endif
#endif

#if defined WAVE_OFFLINE && defined MUSTANG
   ! &croco_wave_offline
   character(len=180) :: wave_file
   !! Offline wave forcing file for MUSTANG sediment model
#endif

#if defined BIOLOGY && defined PISCES
   ! &croco_biology
   character(len=180) :: bioname
   !! Name of file containing the iron dust forcing for the PISCES biogeochemical model
#endif

#ifdef BODYFORCE
   ! &croco_bodyforce
   integer :: levsfrc = 1
   !! Deepest sigma level to apply surface stress as a body force
   integer :: levbfrc = 1
   !! Shallowest sigma level to apply bottom stress as a body force
#endif

#if !defined NONLIN_EOS
   ! &croco_lin_eos
   real :: R0 = 1027.0
   !! Background (reference) density [kg/m3] used in the linear equation of state
   real :: T0 = 14.0
   !! Background (reference) potential temperature [Celsius]
   real :: S0 = 35.0
   !! Background (reference) salinity [PSU]
   real :: Tcoef = 1.7e-4
   !! Thermal expansion coefficient [kg/m3/Celsius]
   real :: Scoef = 7.6e-4
   !! Saline contraction coefficient [kg/m3/PSU]
#endif

#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
   ! &croco_abl_nudg_tra
   real :: ltra_min = 1.0e5
   !! ABL tracer nudging timescale at the bottom of the ABL [s]
   real :: ltra_max = 1.0e4
   !! ABL tracer nudging timescale at the top of the ABL [s]
#endif

#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
   ! &croco_abl_nudg_dyn
   real :: ldyn_min = 1.0e5
   !! ABL dynamics nudging timescale at the bottom of the ABL [s]
   real :: ldyn_max = 1.0e4
   !! ABL dynamics nudging timescale at the top of the ABL [s]
#endif

#ifdef SEDIMENT
   ! &croco_sediments
   character(len=180) :: sedname = "sediment.in"
   !! Sediment model parameters input file
#endif

#ifdef MUSTANG
   ! &croco_sediments_mustang
   character(len=180) :: sedname_must = "paraMUSTANG.txt"
   !! MUSTANG sediment model parameters input file
#endif

#ifdef SUBSTANCE
   ! &croco_substance
   character(len=180) :: subsfilename = "parasubstance.txt"
   !! Substance module parameters input file
#endif

#ifdef OBSTRUCTION
   ! &croco_obstruction
   character(len=180) :: obstname = "obstruction.in"
   !! Sub-grid obstruction module parameters input file
#endif

#ifdef XIOS
   ! &croco_xios_origin_date
   character(len=80) :: xios_origin_date = "1900-01-01 00:00:00"
   !! XIOS time origin date (format: `YYYY-MM-DD HH:MM:SS`)
#endif

#ifdef ASSIMILATION
   ! &croco_assimilation
   character(len=180) :: aparnam = "assimilation.par"
   !! Assimilation parameters file
   character(len=180) :: assname = "assimilation.nc"
   !! Assimilation data file
#endif

   ! &croco_rho0
   real :: rho0 = 1025.0
   !! Mean density [kg/m3] used in the Boussinesq approximation

   ! &croco_bottom_drag
   real :: rdrg = 3.e-4
   !! Drag coefficient for the linear bottom stress formulation [m/s]
   real :: rdrg2 = 3.e-3
   !! Drag coefficient for the constant quadratic bottom stress formulation
   real :: Zobt = 0.02
   !! Bottom roughness length for the Von-Karman (logarithmic) bottom stress formulation [m]
   real :: Cdb_min = 1.e-4
   !! Minimum drag coefficient value for the Von-Karman quadratic bottom stress formulation
   real :: Cdb_max = 1.e-1
   !! Maximum drag coefficient value for the Von-Karman quadratic bottom stress formulation

   ! &croco_gamma2
   real :: gamma2 = -1.0
   !! Lateral boundary slipperiness parameter.
   !! `+1` = free-slip, `-1` = no-slip, intermediate values give partial-slip.

   ! &croco_lateral_visc
   real :: visc2 = 0.0
   !! Laplacian background horizontal momentum viscosity [m2/s] (used with `UV_VIS2`)
   real :: visc4 = 0.0
   !! Bilaplacian background horizontal momentum viscosity [m4/s] (used with `UV_VIS4`)

#ifdef SOLVE3D
#  ifdef TRACERS
   ! &croco_tracer_diff2  – Horizontal Laplacian mixing coefficients for tracers
   real, allocatable :: tnu2(:)
   !! Laplacian background horizontal diffusivity for each tracer [m2/s] (used with `TS_DIF2`).
   !! Array of size NT (number of tracers).

   ! &croco_tracer_diff4  – Horizontal biharmonic mixing coefficients for tracers
   real, allocatable :: tnu4(:)
   !! Bilaplacian background horizontal diffusivity for each tracer [m4/s] (used with `TS_DIF4`).
   !! Array of size NT (number of tracers).
#  endif

#  if !defined LMD_MIXING
   ! &croco_vertical_mixing  – Background vertical mixing coefficients
   real :: Akv_bak = 1.e-5
   !! Background vertical viscosity coefficient [m2/s] (used with `ANA_VMIX`)
#    ifdef TRACERS
   real, allocatable :: Akt_bak(:)
   !! Background vertical diffusivity coefficient for each tracer [m2/s] (used with `ANA_VMIX`).
   !! Array of size NT (number of tracers).
#    endif
#  endif
#endif

#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
   ! &croco_psource
   integer :: psource_Nsrc = 0
   !! Number of point sources/sinks

   ! &croco_psource_ncfile_data / &croco_psource_data
   integer, allocatable :: psource_Isrc(:)
   !! I-grid indices of point sources (size psource_Nsrc)
   integer, allocatable :: psource_Jsrc(:)
   !! J-grid indices of point sources (size psource_Nsrc)
   integer, allocatable :: psource_Dsrc(:)
   !! Direction flag per source: 0 = along XI, >0 = along ETA (size psource_Nsrc)
   character(len=180) :: psource_qbarname = ""
   !! Path to the river runoff NetCDF file (used with PSOURCE_NCFILE)
   real, allocatable :: psource_qbardir(:)
   !! Sign of runoff direction at each source (size psource_Nsrc)
   real, allocatable :: psource_Qbar(:)
   !! Vertically integrated transport [m3/s] at each source (size psource_Nsrc)
# if defined TRACERS
   ! &croco_psource_tracer
   logical, allocatable :: psource_Lsrc(:, :)
   !! Tracer switch per source/tracer (psource_Nsrc,NT): .true. = apply source
   real, allocatable :: psource_Tsrc0(:, :)
   !! Initial tracer value at each source (psource_Nsrc,NT)
# endif
#endif
#if defined DIAGNOSTICS_TS
   ! &croco_diagnostics_ts
   logical :: ldefdia = .true.
   !! Flag (T/F) for writing the tracer diagnostics file
   integer :: nwrtdia = 72
   !! Number of time-steps between writing of tracer diagnostic fields (0 = use nwrt)
   integer :: nrpfdia = 0
   !! Tracer diagnostics file cycling switch (0: multi-record, >0: multi-file, -1: overwrite)
   character(len=180) :: dianame = "CROCO_FILES/croco_dia.nc"
   !! Name of tracer diagnostics file
#  ifdef AVERAGES
   ! &croco_diag_avg
   logical :: ldefdia_avg = .true.
   !! Flag (T/F) for writing averaged tracer diagnostics
   integer :: ntsdia_avg = 1
   !! Starting time-step for accumulation of averaged tracer diagnostics
   integer :: nwrtdia_avg = 72
   !! Number of time-steps between writing of averaged tracer diagnostics (0 = use navg)
   integer :: nrpfdia_avg = 0
   !! Averaged tracer diagnostics file cycling switch
   character(len=180) :: dianame_avg = "CROCO_FILES/croco_dia_avg.nc"
   !! Name of averaged tracer diagnostics file
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
   ! &croco_diag_mld_dens
   real :: mld_crit_D = 0.03
   !! Density criterion to define the mixed layer depth [kg/m3]
   real :: mld_crit_T = 0.2
   !! Temperature criterion to define the mixed layer depth [Celsius]
#  endif
#endif

#if defined DIAGNOSTICS_UV
   ! &croco_diagnosticsM
   logical :: ldefdiaM = .true.
   !! Flag (T/F) for writing momentum diagnostics file
   integer :: nwrtdiaM = 72
   !! Number of time-steps between writing of momentum diagnostics (0 = use nwrt)
   integer :: nrpfdiaM = 0
   !! Momentum diagnostics file cycling switch
   character(len=180) :: dianameM = "CROCO_FILES/croco_diaM.nc"
   !! Name of momentum diagnostics file
#  ifdef AVERAGES
   ! &croco_diagM_avg
   logical :: ldefdiaM_avg = .true.
   !! Flag (T/F) for writing averaged momentum diagnostics
   integer :: ntsdiaM_avg = 1
   !! Starting time-step for accumulation of averaged momentum diagnostics
   integer :: nwrtdiaM_avg = 72
   !! Number of time-steps between writing of averaged momentum diagnostics (0 = use navg)
   integer :: nrpfdiaM_avg = 0
   !! Averaged momentum diagnostics file cycling switch
   character(len=180) :: dianameM_avg = "CROCO_FILES/croco_diaM_avg.nc"
   !! Name of averaged momentum diagnostics file
#  endif
#endif

#ifdef DIAGNOSTICS_VRT
   ! &croco_diags_vrt
   logical :: ldefdiags_vrt = .true.
   !! Flag (T/F) for writing vorticity diagnostics file
   integer :: nwrtdiags_vrt = 72
   !! Number of time-steps between writing of vorticity diagnostics (0 = use nwrt)
   integer :: nrpfdiags_vrt = 0
   !! Vorticity diagnostics file cycling switch
   character(len=180) :: diags_vrtname = "CROCO_FILES/croco_diags_vrt.nc"
   !! Name of vorticity diagnostics file
#  ifdef AVERAGES
   ! &croco_diags_vrt_avg
   logical :: ldefdiags_vrt_avg = .true.
   !! Flag (T/F) for writing averaged vorticity diagnostics
   integer :: ntsdiags_vrt_avg = 1
   !! Starting time-step for accumulation of averaged vorticity diagnostics
   integer :: nwrtdiags_vrt_avg = 72
   !! Number of time-steps between writing of averaged vorticity diagnostics (0 = use navg)
   integer :: nrpfdiags_vrt_avg = 0
   !! Averaged vorticity diagnostics file cycling switch
   character(len=180) :: diags_vrtname_avg = "CROCO_FILES/croco_diags_vrt_avg.nc"
   !! Name of averaged vorticity diagnostics file
#  endif
#endif

#ifdef DIAGNOSTICS_KE
   ! &croco_diags_ek
   logical :: ldefdiags_ek = .true.
   !! Flag (T/F) for writing energy diagnostics file
   integer :: nwrtdiags_ek = 72
   !! Number of time-steps between writing of energy diagnostics (0 = use nwrt)
   integer :: nrpfdiags_ek = 0
   !! Energy diagnostics file cycling switch
   character(len=180) :: diags_ekname = "CROCO_FILES/croco_diags_ek.nc"
   !! Name of energy diagnostics file
#  ifdef AVERAGES
   ! &croco_diags_ek_avg
   logical :: ldefdiags_ek_avg = .true.
   !! Flag (T/F) for writing averaged energy diagnostics
   integer :: ntsdiags_ek_avg = 1
   !! Starting time-step for accumulation of averaged energy diagnostics
   integer :: nwrtdiags_ek_avg = 72
   !! Number of time-steps between writing of averaged energy diagnostics (0 = use navg)
   integer :: nrpfdiags_ek_avg = 0
   !! Averaged energy diagnostics file cycling switch
   character(len=180) :: diags_ekname_avg = "CROCO_FILES/croco_diags_ek_avg.nc"
   !! Name of averaged energy diagnostics file
#  endif
#endif

#ifdef DIAGNOSTICS_PV
   ! &croco_diags_pv
   logical :: ldefdiags_pv = .true.
   !! Flag (T/F) for writing potential vorticity diagnostics file
   integer :: nwrtdiags_pv = 72
   !! Number of time-steps between writing of PV diagnostics (0 = use nwrt)
   integer :: nrpfdiags_pv = 0
   !! PV diagnostics file cycling switch
   character(len=180) :: diags_pvname = "CROCO_FILES/croco_diags_pv.nc"
   !! Name of potential vorticity diagnostics file
#  ifdef AVERAGES
   ! &croco_diags_pv_avg
   logical :: ldefdiags_pv_avg = .true.
   !! Flag (T/F) for writing averaged PV diagnostics
   integer :: ntsdiags_pv_avg = 1
   !! Starting time-step for accumulation of averaged PV diagnostics
   integer :: nwrtdiags_pv_avg = 72
   !! Number of time-steps between writing of averaged PV diagnostics (0 = use navg)
   integer :: nrpfdiags_pv_avg = 0
   !! Averaged PV diagnostics file cycling switch
   character(len=180) :: diags_pvname_avg = "CROCO_FILES/croco_diags_pv_avg.nc"
   !! Name of averaged PV diagnostics file
#  endif
#endif

#if defined DIAGNOSTICS_EDDY && !defined XIOS
#  ifdef AVERAGES
   ! &croco_diags_eddy_avg
   logical :: ldefdiags_eddy_avg = .true.
   !! Flag (T/F) for writing averaged eddy diagnostics
   integer :: ntsdiags_eddy_avg = 1
   !! Starting time-step for accumulation of averaged eddy diagnostics
   integer :: nwrtdiags_eddy_avg = 72
   !! Number of time-steps between writing of averaged eddy diagnostics (0 = use navg)
   integer :: nrpfdiags_eddy_avg = 0
   !! Averaged eddy diagnostics file cycling switch
   character(len=180) :: diags_eddyname_avg = "CROCO_FILES/croco_diags_eddy_avg.nc"
   !! Name of averaged eddy diagnostics file
#  endif
#endif

#ifdef DIAGNOSTICS_BIO
   ! &croco_diagnostics_bio
   logical :: ldefdiabio = .true.
   !! Flag (T/F) for writing biogeochemical diagnostics file
   integer :: nwrtdiabio = 72
   !! Number of time-steps between writing of bio diagnostics (0 = use nwrt)
   integer :: nrpfdiabio = 0
   !! Bio diagnostics file cycling switch
   character(len=180) :: dianamebio = "CROCO_FILES/croco_diabio.nc"
   !! Name of biogeochemical diagnostics file
#  ifdef AVERAGES
   ! &croco_diagbio_avg
   logical :: ldefdiabio_avg = .true.
   !! Flag (T/F) for writing averaged bio diagnostics
   integer :: ntsdiabio_avg = 1
   !! Starting time-step for accumulation of averaged bio diagnostics
   integer :: nwrtdiabio_avg = 72
   !! Number of time-steps between writing of averaged bio diagnostics (0 = use navg)
   integer :: nrpfdiabio_avg = 0
   !! Averaged bio diagnostics file cycling switch
   character(len=180) :: dianamebio_avg = "CROCO_FILES/croco_diabio_avg.nc"
   !! Name of averaged bio diagnostics file
#  endif
#endif

#ifdef STATIONS
   ! &croco_stations
   logical :: ldefsta = .true.
   !! Flag (T/F) for creating a new stations output file
   integer :: nsta = 0
   !! Number of time-steps between storage of station data
   integer :: nrpfsta = 0
   !! Stations file cycling switch
   character(len=180) :: staposname = "stations.in"
   !! Station positions input file
   character(len=180) :: staname = "CROCO_FILES/stations.nc"
   !! Stations output file
   ! &croco_station_fields
   logical :: sta_grd = .true.
   !! Write grid metrics (h, f, pm, pn, angle) to stations file
   logical :: sta_temp = .true.
   !! Write temperature to stations file
   logical :: sta_salt = .true.
   !! Write salinity to stations file
   logical :: sta_rho = .true.
   !! Write density to stations file
   logical :: sta_vel = .true.
   !! Write velocity (u, v, w, ubar, vbar) to stations file
#endif

#if (defined SPONGE && !defined SPONGE_GRID)
   ! &croco_sponge
   real :: x_sponge = 0.0
   !! Thickness of sponge/nudging layer [m]
#endif

#ifdef WET_DRY
   ! &croco_wetdry
   real :: D_wetdry = 0.2
   !! Critical depth for wetting/drying cells [m]
#endif

#if defined MRL_WCI && defined WAVE_DRY
   ! &croco_wavedry
   real :: D_wavedry = 1.0
   !! Minimum water depth above which wave forcing is applied [m]
   !! (should be > D_wetdry when WET_DRY is also active)
#endif

#if defined T_FRC_BRY     || defined M2_FRC_BRY    || \
   defined M3_FRC_BRY||defined Z_FRC_BRY||\
   defined W_FRC_BRY||defined NBQ_FRC_BRY||\
   defined TCLIMATOLOGY||defined M2CLIMATOLOGY||\
   defined M3CLIMATOLOGY||defined ZCLIMATOLOGY||\
   defined WCLIMATOLOGY||defined NBQCLIMATOLOGY||\
   defined TNUDGING||defined ZNUDGING||\
   defined M2NUDGING||defined M3NUDGING||\
   defined NBQ_NUDGING||defined SPONGE
   ! &croco_nudging
   real :: tauT_in = 1.0
   !! Inflow nudging time-scale for tracers [days]
   real :: tauT_out = 360.0
   !! Outflow nudging time-scale for tracers [days]
   real :: tauM_in = 3.0
   !! Inflow nudging time-scale for momentum [days]
   real :: tauM_out = 360.0
   !! Outflow nudging time-scale for momentum [days]
#endif

#if defined BHFLUX || (defined BWFLUX && defined SALINITY)
   ! &croco_bottom_forcing
   character(len=180) :: btfname = "CROCO_FILES/croco_btf.nc"
   !! Name of bottom hydrothermal forcing file
#endif

   ! &croco_primary_history_fields / &croco_primary_history_3d_fields
   logical :: out_his_zeta = .true.
   logical :: out_his_ubar = .true.
   logical :: out_his_vbar = .true.
#ifdef SOLVE3D
   logical :: out_his_u = .true.
   logical :: out_his_v = .true.
# ifdef TRACERS
   logical, allocatable :: out_his_tracer(:)
# endif
#endif

#ifdef AVERAGES
   ! &croco_primary_average_fields / &croco_primary_3d_average_fields
   logical :: out_avg_zeta = .true.
   logical :: out_avg_ubar = .true.
   logical :: out_avg_vbar = .true.
# ifdef SOLVE3D
   logical :: out_avg_u = .true.
   logical :: out_avg_v = .true.
#  ifdef TRACERS
   logical, allocatable :: out_avg_tracer(:)
#  endif
# endif
#endif

#ifdef DIAGNOSTICS_TS
# ifdef TRACERS
   ! &croco_diag3D_history_fields
   logical, allocatable :: out_his_dia3D_tracer(:)
#  ifdef DIAGNOSTICS_TS_MLD
   ! &croco_diag2D_history_fields
   logical, allocatable :: out_his_dia2D_tracer(:)
#  endif
# endif
# if defined AVERAGES && defined TRACERS
   ! &croco_diag3D_average_fields
   logical, allocatable :: out_avg_dia3D_tracer(:)
#  ifdef DIAGNOSTICS_TS_MLD
   ! &croco_diag2D_average_fields
   logical, allocatable :: out_avg_dia2D_tracer(:)
#  endif
# endif
#endif

#ifdef DIAGNOSTICS_UV
   ! &croco_diagM_history_fields
   logical :: out_his_diagM_u = .true.
   logical :: out_his_diagM_v = .true.
# ifdef AVERAGES
   ! &croco_diagM_average_fields
   logical :: out_avg_diagM_u = .true.
   logical :: out_avg_diagM_v = .true.
# endif
#endif

#ifdef DIAGNOSTICS_VRT
   ! &croco_diags_vrt_history_fields
   logical :: out_his_diags_vrt = .true.
# ifdef AVERAGES
   ! &croco_diags_vrt_average_fields
   logical :: out_avg_diags_vrt = .true.
# endif
#endif

#ifdef DIAGNOSTICS_KE
   ! &croco_diags_ek_history_fields
   logical :: out_his_diags_ek = .true.
# ifdef AVERAGES
   ! &croco_diags_ek_average_fields
   logical :: out_avg_diags_ek = .true.
# endif
#endif

#if defined DIAGNOSTICS_PV && defined TRACERS
   ! &croco_diags_pv_history_fields
   logical, allocatable :: out_his_diags_pv_tracer(:)
# ifdef AVERAGES
   ! &croco_diags_pv_average_fields
   logical, allocatable :: out_avg_diags_pv_tracer(:)
# endif
#endif

#if defined DIAGNOSTICS_EDDY && !defined XIOS
# ifdef AVERAGES
   ! &croco_diags_eddy_average_fields
   logical :: out_avg_diags_eddy = .true.
# endif
#endif

#ifdef DIAGNOSTICS_BIO
   ! &croco_diagbioFlux/VSink/GasExc_history_fields
   logical, allocatable :: out_his_diagbioFlux(:)
   logical, allocatable :: out_his_diagbioVSink(:)
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
   logical, allocatable :: out_his_diagbioGasExc(:)
# endif
# ifdef AVERAGES
   ! &croco_diagbioFlux/VSink/GasExc_average_fields
   logical, allocatable :: out_avg_diagbioFlux(:)
   logical, allocatable :: out_avg_diagbioVSink(:)
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
   logical, allocatable :: out_avg_diagbioGasExc(:)
#  endif
# endif
#endif

#ifdef ONLINE
   ! &croco_online
   integer :: yearnum = 0
   !! Start year for online bulk forcing
   integer :: monthnum = 1
   !! Start month for online bulk forcing
   integer :: recordsperday = 1
   !! Number of forcing records per day
   integer :: yearend = 0
   !! End year for online bulk forcing
   integer :: monthend = 12
   !! End month for online bulk forcing
   character(len=250) :: pathbulk = "./"
   !! Path to bulk forcing data files
#endif

#ifdef STOGEN
   ! &croco_stochastic_history_fields
   logical :: out_his_xi2d = .true.
   logical :: out_his_xi3d = .true.
#endif

#if defined SOLVE3D && defined GLS_MIXING
   ! &croco_gls_history_fields
   logical :: out_his_tke = .true.
   logical :: out_his_gls = .true.
   logical :: out_his_lscale = .true.
# ifdef AVERAGES
   ! &croco_gls_averages_fields
   logical :: out_avg_tke = .true.
   logical :: out_avg_gls = .true.
   logical :: out_avg_lscale = .true.
# endif
#endif

#if defined SOLVE3D && defined SEDIMENT
   ! &croco_sediment_history_fields
   logical :: out_his_sed_athk = .true.
   logical :: out_his_sed_bthk = .true.
   logical :: out_his_sed_bpor = .true.
   ! &croco_sediment_bfra_history_fields
   logical, allocatable :: out_his_sed_bfra(:)
# ifdef SUSPLOAD
   ! &croco_sediment_suspload_history_fields
   logical, allocatable :: out_his_sed_dflx(:)
   logical, allocatable :: out_his_sed_eflx(:)
# endif
# ifdef BEDLOAD
   ! &croco_sediment_bedload_history_fields
   logical, allocatable :: out_his_sed_bdlu(:)
   logical, allocatable :: out_his_sed_bdlv(:)
# endif
# if defined MIXED_BED || defined COHESIVE_BED
   ! &croco_sediment_cohesive_history_fields
   logical :: out_his_sed_btcr = .true.
# endif
#endif

#ifdef BBL
   ! &croco_bbl_history_fields
   logical :: out_his_abed = .true.
   logical :: out_his_hripple = .true.
   logical :: out_his_lripple = .true.
   logical :: out_his_zbnot = .true.
   logical :: out_his_zbapp = .true.
   logical :: out_his_bostrw = .true.
#endif

#ifdef MRL_WCI
   ! &croco_wci_history_fields
   logical :: out_his_sup = .true.
   logical :: out_his_ust2d = .true.
   logical :: out_his_vst2d = .true.
# ifdef SOLVE3D
   ! &croco_wci_history_3d_fields
   logical :: out_his_ust = .true.
   logical :: out_his_vst = .true.
   logical :: out_his_wst = .true.
   logical :: out_his_akb = .true.
   logical :: out_his_akw = .true.
   logical :: out_his_kvf = .true.
   logical :: out_his_calp = .true.
   logical :: out_his_kaps = .true.
# endif
# ifdef AVERAGES
   ! &croco_wci_average_fields
   logical :: out_avg_sup = .true.
   logical :: out_avg_ust2d = .true.
   logical :: out_avg_vst2d = .true.
#  ifdef SOLVE3D
   ! &croco_wci_average_3d_fields
   logical :: out_avg_ust = .true.
   logical :: out_avg_vst = .true.
   logical :: out_avg_wst = .true.
   logical :: out_avg_akb = .true.
   logical :: out_avg_akw = .true.
   logical :: out_avg_kvf = .true.
   logical :: out_avg_calp = .true.
   logical :: out_avg_kaps = .true.
#  endif
# endif
#endif

#if defined MRL_WCI || defined OW_COUPLING
   ! &croco_wave_history_fields
   logical :: out_his_hrm = .true.
   logical :: out_his_frq = .true.
   logical :: out_his_action = .true.
   logical :: out_his_k_xi = .true.
   logical :: out_his_k_eta = .true.
   logical :: out_his_eps_b = .true.
   logical :: out_his_eps_d = .true.
   logical :: out_his_erol = .true.
   logical :: out_his_eps_r = .true.
# ifdef AVERAGES
   ! &croco_wave_average_fields
   logical :: out_avg_hrm = .true.
   logical :: out_avg_frq = .true.
   logical :: out_avg_action = .true.
   logical :: out_avg_k_xi = .true.
   logical :: out_avg_k_eta = .true.
   logical :: out_avg_eps_b = .true.
   logical :: out_avg_eps_d = .true.
   logical :: out_avg_erol = .true.
   logical :: out_avg_eps_r = .true.
# endif
#endif

#if defined ABL1D && !defined XIOS
   ! &croco_abl_history_fields
   logical :: out_his_abl_pu_dta = .true.
   logical :: out_his_abl_pv_dta = .true.
   logical :: out_his_abl_pt_dta = .true.
   logical :: out_his_abl_pq_dta = .true.
   logical :: out_his_abl_pgu_dta = .true.
   logical :: out_his_abl_pgv_dta = .true.
   logical :: out_his_abl_u_abl = .true.
   logical :: out_his_abl_v_abl = .true.
   logical :: out_his_abl_t_abl = .true.
   logical :: out_his_abl_q_abl = .true.
   logical :: out_his_abl_tke_abl = .true.
   logical :: out_his_abl_mxlm_abl = .true.
   logical :: out_his_abl_mxld_abl = .true.
   logical :: out_his_abl_avm_abl = .true.
   logical :: out_his_abl_avt_abl = .true.
   logical :: out_his_abl_ablh_abl = .true.
   logical :: out_his_abl_zr_abl = .true.
   logical :: out_his_abl_zw_abl = .true.
   logical :: out_his_abl_Hzr_abl = .true.
   logical :: out_his_abl_Hzw_abl = .true.
# ifdef AVERAGES
   ! &croco_abl_averages_fields
   logical :: out_avg_abl_pu_dta = .true.
   logical :: out_avg_abl_pv_dta = .true.
   logical :: out_avg_abl_pt_dta = .true.
   logical :: out_avg_abl_pq_dta = .true.
   logical :: out_avg_abl_pgu_dta = .true.
   logical :: out_avg_abl_pgv_dta = .true.
   logical :: out_avg_abl_u_abl = .true.
   logical :: out_avg_abl_v_abl = .true.
   logical :: out_avg_abl_t_abl = .true.
   logical :: out_avg_abl_q_abl = .true.
   logical :: out_avg_abl_tke_abl = .true.
   logical :: out_avg_abl_mxlm_abl = .true.
   logical :: out_avg_abl_mxld_abl = .true.
   logical :: out_avg_abl_avm_abl = .true.
   logical :: out_avg_abl_avt_abl = .true.
   logical :: out_avg_abl_ablh_abl = .true.
   logical :: out_avg_abl_zr_abl = .true.
   logical :: out_avg_abl_zw_abl = .true.
   logical :: out_avg_abl_Hzr_abl = .true.
   logical :: out_avg_abl_Hzw_abl = .true.
# endif
#endif

# if defined SOLVE3D || defined RIP
   ! &croco_auxiliary_history_fields  (SOLVE3D basics + surface stress)
   logical :: out_his_rho = .false.
   logical :: out_his_omega = .false.
   logical :: out_his_w = .false.
   logical :: out_his_akv = .false.
   logical :: out_his_bostr = .false.
   logical :: out_his_bustr = .false.
   logical :: out_his_bvstr = .false.
   logical :: out_his_wstr = .false.
   logical :: out_his_ustr = .false.
   logical :: out_his_vstr = .false.
# endif
# if defined SOLVE3D && defined TEMPERATURE
   ! &croco_temperature_history_fields
   logical :: out_his_akt = .false.
   logical :: out_his_shflx = .false.
   logical :: out_his_shflx_rsw = .false.
# endif
# if defined SOLVE3D && defined SALINITY
   ! &croco_salinity_history_fields
   logical :: out_his_aks = .false.
   logical :: out_his_swflx = .false.
# endif
# if defined SOLVE3D && defined BULK_FLUX
   ! &croco_bulk_flux_history_fields
   logical :: out_his_shflx_rlw = .false.
   logical :: out_his_shflx_lat = .false.
   logical :: out_his_shflx_sen = .false.
# endif
# if defined SOLVE3D && defined BHFLUX
   ! &croco_bhflux_history_fields
   logical :: out_his_bhflx = .false.
# endif
# if defined SOLVE3D && defined BWFLUX && defined SALINITY
   ! &croco_bwflux_history_fields
   logical :: out_his_bwflx = .false.
# endif
# if defined SOLVE3D && \
   (defined ANA_VMIX||defined LMD_MIXING||\
   defined LMD_SKPP||defined LMD_BKPP||\
   defined GLS_MIXING)
   ! &croco_bvf_history_fields
   logical :: out_his_bvf = .false.
# endif
# if defined SOLVE3D && (defined LMD_SKPP || defined GLS_MIXING)
   ! &croco_hbl_history_fields
   logical :: out_his_hbl = .false.
# endif
# if defined SOLVE3D && defined LMD_BKPP
   ! &croco_lmd_bkpp_history_fields
   logical :: out_his_hbbl = .false.
# endif
# if defined SOLVE3D && defined VIS_COEF_3D
   ! &croco_vis_coef_history_fields
   logical :: out_his_visc3d = .false.
# endif
# if defined SOLVE3D && defined DIF_COEF_3D
   ! &croco_dif_coef_history_fields
   logical :: out_his_diff3d = .false.
# endif
# if defined SOLVE3D && defined BIOLOGY && !defined PISCES
   ! &croco_biology_history_fields
   logical :: out_his_hel = .false.
# endif
# if defined SOLVE3D && defined BIO_NChlPZD
   ! &croco_bio_nchlpzd_history_fields
   logical :: out_his_chc = .false.
   logical :: out_his_u10 = .false.
   logical :: out_his_kvo2 = .false.
   logical :: out_his_o2sat = .false.
# endif
# if defined SOLVE3D && defined BIO_BioEBUS
   ! &croco_bio_bioebus_history_fields
   logical :: out_his_aou = .false.
   logical :: out_his_wind10 = .false.
# endif
# if defined SOLVE3D && defined MORPHODYN
   ! &croco_morphodyn_history_fields
   logical :: out_his_hm = .true.
# endif
# if defined AVERAGES && defined SOLVE3D
   ! &croco_auxiliary_averages_fields  (SOLVE3D basics + surface stress)
   logical :: out_avg_rho = .false.
   logical :: out_avg_omega = .false.
   logical :: out_avg_w = .false.
   logical :: out_avg_akv = .false.
   logical :: out_avg_bostr = .false.
   logical :: out_avg_bustr = .false.
   logical :: out_avg_bvstr = .false.
   logical :: out_avg_wstr = .false.
   logical :: out_avg_ustr = .false.
   logical :: out_avg_vstr = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined TEMPERATURE
   ! &croco_temperature_averages_fields
   logical :: out_avg_akt = .false.
   logical :: out_avg_shflx = .false.
   logical :: out_avg_shflx_rsw = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined SALINITY
   ! &croco_salinity_averages_fields
   logical :: out_avg_aks = .false.
   logical :: out_avg_swflx = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined BULK_FLUX
   ! &croco_bulk_flux_averages_fields
   logical :: out_avg_shflx_rlw = .false.
   logical :: out_avg_shflx_lat = .false.
   logical :: out_avg_shflx_sen = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined BHFLUX
   ! &croco_bhflux_averages_fields
   logical :: out_avg_bhflx = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined BWFLUX && defined SALINITY
   ! &croco_bwflux_averages_fields
   logical :: out_avg_bwflx = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && \
   (defined ANA_VMIX||defined LMD_MIXING||\
   defined LMD_SKPP||defined LMD_BKPP||\
   defined GLS_MIXING)
   ! &croco_bvf_averages_fields
   logical :: out_avg_bvf = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && \
   (defined LMD_SKPP||defined GLS_MIXING)
   ! &croco_hbl_averages_fields
   logical :: out_avg_hbl = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined LMD_BKPP
   ! &croco_lmd_bkpp_averages_fields
   logical :: out_avg_hbbl = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined VIS_COEF_3D
   ! &croco_vis_coef_averages_fields
   logical :: out_avg_visc3d = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined DIF_COEF_3D
   ! &croco_dif_coef_averages_fields
   logical :: out_avg_diff3d = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIOLOGY && !defined PISCES
   ! &croco_biology_averages_fields
   logical :: out_avg_hel = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_NChlPZD
   ! &croco_bio_nchlpzd_averages_fields
   logical :: out_avg_chc = .false.
   logical :: out_avg_u10 = .false.
   logical :: out_avg_kvo2 = .false.
   logical :: out_avg_o2sat = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_BioEBUS
   ! &croco_bio_bioebus_averages_fields
   logical :: out_avg_aou = .false.
   logical :: out_avg_wind10 = .false.
# endif
# if defined AVERAGES && defined SOLVE3D && defined MORPHODYN
   ! &croco_morphodyn_averages_fields
   logical :: out_avg_hm = .false.
# endif

END MODULE croco_namelist
