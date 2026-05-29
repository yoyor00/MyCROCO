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
   !! Logfile name

   ! &croco_time_stepping
   real    :: dt = 0.0
   !! Baroclinic time step [in s]
   integer :: ntimes = 0
   !! Number of time-steps required for the simulation
   integer :: ndtfast = 20
   !! Number of barotropic time-steps between each baroclinic time step.
   !! For 2D configurations, ndtfast should be unity
   integer :: ninfo = 1
   !! Number of time-steps between printing of information to standard output

   ! &croco_history
   logical :: ldefhis = .true.
   !! Logical switch used to create the history file.
   !! If TRUE, a new history file is created. If FALSE,
   !! data is appended to an existing history file.
   integer :: nwrt = 72
   !! Number of timesteps between writing of fields into history file.
   integer :: nrpfhis = 0
   !! 0: write several records every NWRT time steps
   !! >0: create more than one file (sequential numbers), NRPFHIS records per file
   !! -1: overwrite record every NWRT time steps
   character(len=180) :: hisname = "CROCO_FILES/croco_his.nc"
   !! Name of history file

   ! &croco_initial
   integer :: nrrec = 1
   !! Switch to indicate start or re-start from a previous solution.
   !! nrrec is the time index of the initial or re-start NetCDF file
   !! assigned for initialization. If nrrec is negative (say nrrec = -1),
   !! the model will start from the most recent time record.
   !! That is, the initialization record is assigned internally.
   character(len=180) :: ininame = "CROCO_FILES/croco_ini.nc"
   !! Name of file containing the initial state.

   ! &croco_restart
   integer :: nrst = 720
   !! Number of time-steps between writing of re-start fields
   integer :: nrpfrst = -1
   !! 0: write several records every NRST time steps
   !! >0: create more than one file (with sequential numbers) and write NRPRST records per file
   !! -1: overwrite record every NRST time steps
   character(len=180) :: rstname = "CROCO_FILES/croco_rst.nc"
   !! Name of restart file.

#ifdef NBQ
   ! &croco_time_stepping_nbq
   integer :: ndtnbq = 1
   real    :: csound_nbq = 1000.0
   real    :: visc2_nbq = 0.01
#endif

#ifdef SOLVE3D
   ! &croco_s_coord
   real :: theta_s = 7.0d0
   real :: theta_b = 2.0d0
   real :: Tcline = 200.0d0
#endif

#ifdef USE_CALENDAR
   ! &croco_use_calendar
   character(len=19) :: start_date = '2000-01-01 00:00:00'
   character(len=19) :: end_date = '2000-02-01 00:00:00'
   real :: dt_his = 1.0
   real :: dt_avg = 6.0
   real :: dt_rst = 12.0
#endif

#ifndef ANA_GRID
   ! &croco_grid
   character(len=180) :: grdname = "CROCO_FILES/croco_grd.nc"
#endif

   ! &croco_forcing
   character(len=180) :: frcname = "CROCO_FILES/croco_frc.nc"
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
#endif

#if !defined ANA_BRY && defined FRC_BRY
   ! &croco_boundary
   character(len=180) :: bry_file = "CROCO_FILES/croco_bry.nc"
#endif

#if defined WKB_WWAVE && !defined ANA_BRY_WKB
   ! &croco_wkb_boundary
   character(len=180) :: brywkb_file = "CROCO_FILES/croco_wkb.nc"
#endif

#ifdef WKB_WWAVE
   ! &croco_wkb_wwave  – primary waves and empirical breaking model parameters
   real :: wkb_amp = 0.25
   !! offshore wave amplitude [m]
   real :: wkb_ang = 190.0
   !! offshore wave angle [deg]
   real :: wkb_prd = 8.0
   !! offshore wave period [s]
   real :: wkb_tide = -2.0
   !! constant offshore water level [m]
   real :: wkb_btg = 1.3
   !! B parameter (breaking)
   real :: wkb_gam = 0.38
   !! gamma parameter (Hrms/h ratio)
#  ifdef WAVE_ROLLER
   ! &croco_wkb_roller  – Svendsen (1984) surface roller model parameters
   real :: wkb_rsb = 0.1
   !! sin(beta) roller dissipation
   real :: wkb_roller = 0.5
   !! breaking contrib to roller: [0,1]
#  endif
#endif

#ifdef WAVE_MAKER
   ! &croco_wave_maker  – offshore wavemaker parameters for wave-resolving simulations
   real :: wmaker_amp = 0.0
   !! RMS wave amplitude [m]
   real :: wmaker_prd = 8.0
   !! peak wave period [s]
   real :: wmaker_dir = 0.0
   !! mean wave angle [deg]
   real :: wmaker_dsp = 0.0
   !! directional spread [deg]
   real :: wmaker_fsp = 3.3
   !! freq. spread (gamma in JONSWAP)
#endif

#ifdef AVERAGES
   ! &croco_averages
   integer :: ntsavg = 1
   integer :: navg = 48
   integer :: nrpfavg = 0
   character(len=180) :: avgname = "CROCO_FILES/croco_avg.nc"
#endif

#if defined ABL1D && !defined XIOS
   ! &croco_abl
   logical :: ldefablhis = .true.
   integer :: nwrtablhis = 36
   integer :: nrpfablhis = 0
   character(len=180) :: ablname = "croco_abl_his.nc"
#  ifdef AVERAGES
   ! &croco_abl_averages
   logical :: ldefablavg = .false.
   integer :: ntsablavg = 1
   integer :: nwrtablavg = 0
   integer :: nrpfablavg = 0
   character(len=180) :: ablname_avg = "croco_abl_avg.nc"
#  endif
#endif

#if defined WAVE_OFFLINE && defined MUSTANG
   ! &croco_wave_offline
   character(len=180) :: wave_file
#endif

#if defined BIOLOGY && defined PISCES
   ! &croco_biology
   character(len=180) :: bioname
#endif

#ifdef BODYFORCE
   ! &croco_bodyforce
   integer :: levsfrc = 1
   !! Deepest level to apply surface stress as a bodyforce
   integer :: levbfrc = 1
   !! Shallowest level to apply bottom stress as a bodyforce
#endif

#if !defined NONLIN_EOS
   ! &croco_lin_eos
   real :: R0 = 1027.0
   !! Background density [kg/m3] used in linear EOS
   real :: T0 = 14.0
   !! Background potential temperature [Celsius]
   real :: S0 = 35.0
   !! Background salinity [PSU]
   real :: Tcoef = 1.7e-4
   !! Thermal expansion coefficient [kg/m3/Celsius]
   real :: Scoef = 7.6e-4
   !! Saline contraction coefficient [kg/m3/PSU]
#endif

#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
   ! &croco_abl_nudg_tra
   real :: ltra_min = 1.0e5
   !! ABL tracer nudging timescale at bottom [s]
   real :: ltra_max = 1.0e4
   !! ABL tracer nudging timescale at top [s]
#endif

#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
   ! &croco_abl_nudg_dyn
   real :: ldyn_min = 1.0e5
   !! ABL dynamics nudging timescale at bottom [s]
   real :: ldyn_max = 1.0e4
   !! ABL dynamics nudging timescale at top [s]
#endif

#ifdef SEDIMENT
   ! &croco_sediments
   character(len=180) :: sedname = "sediment.in"
   !! Sediment model input file
#endif

#ifdef MUSTANG
   ! &croco_sediments_mustang
   character(len=180) :: sedname_must = "paraMUSTANG.txt"
   !! MUSTANG sediment model input file
#endif

#ifdef SUBSTANCE
   ! &croco_substance
   character(len=180) :: subsfilename = "parasubstance.txt"
   !! Substance module input file
#endif

#ifdef OBSTRUCTION
   ! &croco_obstruction
   character(len=180) :: obstname = "obstruction.in"
   !! Obstruction module input file
#endif

END MODULE croco_namelist
