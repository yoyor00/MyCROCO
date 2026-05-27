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

END MODULE croco_namelist
