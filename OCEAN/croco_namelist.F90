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
!            No I/O, no initialisation of derived parameters.
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

#ifdef NBQ
   ! &croco_time_stepping_nbq
   integer :: ndtnbq = 1
   real    :: csound_nbq = 1000.0
   real    :: visc2_nbq = 0.01
#endif

#ifdef USE_CALENDAR
   ! &croco_use_calendar
   character(len=19) :: start_date = '2000-01-01 00:00:00'
   character(len=19) :: end_date = '2000-02-01 00:00:00'
   real :: dt_his = 1.0
   real :: dt_avg = 6.0
   real :: dt_rst = 12.0
#endif

contains

   !---------------------------------------------------------------------
   !  read_nml_fname
   !  Set fname_nml from command-line argument 2 (defaults to croco.nml).
   !---------------------------------------------------------------------
   subroutine read_nml_fname()
      implicit none

      call get_command_argument(2, fname_nml)
      if (len_trim(fname_nml) == 0) fname_nml = 'croco.nml'

#ifdef AGRIF
      if (.Not. Agrif_Root()) then
#  ifdef AGRIF_ADAPTIVE
         fname_nml = trim(fname_nml)//'.1'
#  else
         fname_nml = trim(fname_nml)//'.'//Agrif_Cfixed()
#  endif
      end if
#endif
   end subroutine read_nml_fname

   !---------------------------------------------------------------------
   !  read_nml
   !  Open fname_nml, read all namelists into module variables,
   !  then delegate checks and initialisations to croco_check /
   !  croco_init. Returns ierr /= 0 on any fatal error.
   !
   !  For each namelist section the sequence is:
   !    1. check_nml_presence : detect whether the section exists in the
   !       file (missing optional -> warning; missing mandatory -> fatal)
   !    2. read(..., nml=xxx) : parse the section if present; a non-zero
   !       iostat here means a malformed parameter, which is always fatal
   !    3. check_xxx          : validate the values that were just read
   !    4. init_xxx           : compute derived parameters
   !---------------------------------------------------------------------
   subroutine read_nml(ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
      implicit none
      integer, intent(out) :: ierr

      integer :: nmlunit, ios
      logical :: found

      namelist /croco_title/ title
      namelist /croco_logfile/ logname
      namelist /croco_time_stepping/ dt, ntimes, ndtfast, ninfo
      namelist /croco_history/ ldefhis, nwrt, nrpfhis, hisname
#ifdef NBQ
      namelist /croco_time_stepping_nbq/ ndtnbq, csound_nbq, visc2_nbq
#endif
#ifdef USE_CALENDAR
      namelist /croco_use_calendar/ start_date, end_date, dt_his, dt_avg, dt_rst
#endif

      ierr = 0
      nmlunit = 10

      MPI_master_only write (stdout, *) '*** READING NAMELIST ***'

      ! Open the namelist file; if absent, fall back to default values
      ! for all parameters (non-fatal: the model can still run)
      open (unit=nmlunit, file=trim(fname_nml), status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         MPI_master_only write (stdout, *) &
            'WARNING: namelist file not found: ', trim(fname_nml)
         MPI_master_only write (stdout, *) 'Using default values.'
         return
      end if

      ! --- croco_title (optional) ---
      call check_nml_presence(nmlunit, "croco_title", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_title, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_title (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if

      ! --- croco_time_stepping (mandatory) ---
      call check_nml_presence(nmlunit, "croco_time_stepping", .true., found, ierr)
      if (ierr /= 0) return
      if (found) then
         read (nmlunit, nml=croco_time_stepping, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_time_stepping (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if
      call check_time_stepping(ierr)
      if (ierr /= 0) then
         call fatal_nml_error("croco_time_stepping", nmlunit); return
      end if
      call init_time_stepping()

#ifdef NBQ
      ! --- croco_time_stepping_nbq (mandatory if NBQ) ---
      call check_nml_presence(nmlunit, "croco_time_stepping_nbq", .true., found, ierr)
      if (ierr /= 0) return
      if (found) then
         read (nmlunit, nml=croco_time_stepping_nbq, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_time_stepping_nbq (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if
      call check_time_stepping_nbq(ierr)
      if (ierr /= 0) then
         call fatal_nml_error("croco_time_stepping_nbq", nmlunit); return
      end if
      call init_time_stepping_nbq()
#endif

#ifdef USE_CALENDAR
      ! --- croco_use_calendar (mandatory if USE_CALENDAR) ---
      call check_nml_presence(nmlunit, "croco_use_calendar", .true., found, ierr)
      if (ierr /= 0) return
      if (found) then
         read (nmlunit, nml=croco_use_calendar, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_use_calendar (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if
      call init_calendar()
#endif

      ! --- croco_logfile (optional) ---
      ! Read last so that all warnings above are printed to stdout
      ! before the log file is opened and stdout is redirected.
#ifdef LOGFILE
      call check_nml_presence(nmlunit, "croco_logfile", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_logfile, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_logfile (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if
      call init_logfile(ierr)
      if (ierr /= 0) then
         call fatal_nml_error("croco_logfile", nmlunit); return
      end if
#endif

      ! --- croco_history (mandatory) ---
      call check_nml_presence(nmlunit, "croco_history", .true., found, ierr)
      if (ierr /= 0) return
      if (found) then
         read (nmlunit, nml=croco_history, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_history (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if
      call check_history(ierr)
      if (ierr /= 0) then
         call fatal_nml_error("croco_history", nmlunit); return
      end if
      call init_history(ierr)
      if (ierr /= 0) then
         call fatal_nml_error("croco_history (init)", nmlunit); return
      end if

      close (nmlunit)

      ! Echo all namelists to stdout so the user can verify what was read
      MPI_master_only WRITE (stdout, nml=croco_title)
#ifdef LOGFILE
      MPI_master_only WRITE (stdout, nml=croco_logfile)
#endif
      MPI_master_only WRITE (stdout, nml=croco_time_stepping)
#ifdef NBQ
      MPI_master_only WRITE (stdout, nml=croco_time_stepping_nbq)
#endif
#ifdef USE_CALENDAR
      MPI_master_only WRITE (stdout, nml=croco_use_calendar)
#endif
      MPI_master_only WRITE (stdout, nml=croco_history)

   end subroutine read_nml

   !---------------------------------------------------------------------
   !  check_nml_presence
   !
   !  Scan the namelist file line by line looking for the section header
   !  "&nml_name" (case-insensitive).  On return:
   !    found = .true.  -> section exists; the caller should read it
   !    found = .false. -> section is absent
   !       is_mandatory = .true.  -> fatal error, ierr is incremented
   !       is_mandatory = .false. -> warning only, ierr is unchanged
   !
   !  The file is rewound before and after the search so that subsequent
   !  reads start from the beginning.
   !---------------------------------------------------------------------
   subroutine check_nml_presence(nmlunit, nml_name, is_mandatory, found, ierr)
      use param, ONLY: stdout
      use tools_string, ONLY: to_lowercase
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(in)    :: nmlunit
      character(len=*), intent(in)    :: nml_name
      logical, intent(in)    :: is_mandatory
      logical, intent(out)   :: found
      integer, intent(inout) :: ierr

      character(len=256) :: line
      integer            :: ios

      ! Search for the section header "&nml_name" at the start of a line
      ! (after leading blanks).  Fortran namelist headers are case-
      ! insensitive, so we compare lower-case versions of both strings.
      found = .false.
      rewind (nmlunit)
      do
         read (nmlunit, '(A)', iostat=ios) line
         if (ios /= 0) exit   ! end of file or read error
         if (index(adjustl(to_lowercase(line)), &
                   '&'//to_lowercase(trim(nml_name))) == 1) then
            found = .true.
            exit
         end if
      end do
      rewind (nmlunit)

      if (.not. found) then
         if (is_mandatory) then
            MPI_master_only write (stdout, '(3a)') &
               'Error - mandatory namelist section "', &
               trim(nml_name), '" not found.'
            ierr = ierr + 1
         else
            MPI_master_only write (stdout, '(3a)') &
               'Warning - optional namelist section "', trim(nml_name), &
               '" not found, default values are used.'
         end if
      end if

   end subroutine check_nml_presence

   !---------------------------------------------------------------------
   !  fatal_nml_error
   !
   !  Print a single "aborting" line and close the namelist file.
   !  Detailed error messages are the responsibility of the check_xxx
   !  routines; this helper only provides the final summary line.
   !---------------------------------------------------------------------
   subroutine fatal_nml_error(nml_name, nmlunit)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      character(len=*), intent(in) :: nml_name
      integer, intent(in) :: nmlunit

      MPI_master_only write (stdout, '(3a)') &
         'FATAL: aborting after errors in "', trim(nml_name), '" namelist.'
      close (nmlunit)

   end subroutine fatal_nml_error


   !---------------------------------------------------------------------
   !  init_time_stepping
   !  Derive dtfast from dt and ndtfast.
   !---------------------------------------------------------------------
   subroutine init_time_stepping()
      use scalars, ONLY: dtfast
      implicit none

      dtfast = dt/real(ndtfast)

   end subroutine init_time_stepping

   !---------------------------------------------------------------------
   !  init_history
   !  Adjust hisname for MPI/ENSEMBLE and derive nwrt from dt_his
   !  when USE_CALENDAR is active.
   !---------------------------------------------------------------------
   subroutine init_history(ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode, mynode2, NNODES2
#endif
      implicit none
#ifdef ENSEMBLE
#include "mpi_cpl.h"
#endif
      integer, intent(inout) :: ierr

#if defined MPI && defined PARALLEL_FILES
      call insert_node(hisname, len_trim(hisname), mynode2, NNODES2, ierr)
      if (ierr /= 0) then
         MPI_master_only write (stdout, *) &
            'Error in init_history: insert_node failed for hisname'
         return
      end if
#endif

      hisname = trim(hisname)

#ifdef ENSEMBLE
      hisname = cmember//hisname
#endif

#ifdef USE_CALENDAR
      nwrt = ceiling((dt_his*3600.0)/dt)
#endif

   end subroutine init_history

#ifdef LOGFILE
   !---------------------------------------------------------------------
   !  init_logfile
   !  Open the log file and redirect stdout to it.
   !---------------------------------------------------------------------
   subroutine init_logfile(ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode, mynode2, NNODES2
#endif
      implicit none
#ifdef ENSEMBLE
#include "mpi_cpl.h"
#endif
      integer, intent(inout) :: ierr

      integer :: ios
      character(len=len(logname) + 20) :: logfile_path

#if defined MPI && defined PARALLEL_FILES
      call insert_node(logname, len_trim(logname), mynode2, NNODES2, ierr)
      if (ierr /= 0) then
         write (6, *) 'Error in init_logfile: insert_node failed for logname'
         return
      end if
#endif

      if (mynode == 0) then
#ifndef AGRIF
         stdout = 80
         logfile_path = trim(logname)
#else
         stdout = Agrif_Get_Unit()
         if (.Not. Agrif_Root()) then
            logfile_path = trim(logname)//'.'//Agrif_Cfixed()
         else
            logfile_path = trim(logname)
         end if
#endif

#ifdef ENSEMBLE
         logfile_path = trim(cmember)//trim(logname)
#endif

         open (unit=stdout, file=trim(logfile_path), &
               status='REPLACE', form='formatted', iostat=ios)
         if (ios /= 0) then
            write (6, '(/1x,3A/)') 'ERROR: Cannot open log file: ', trim(logfile_path)
            ierr = ierr + 1
            return
         end if
      else
         stdout = 6
      end if

   end subroutine init_logfile
#endif /* LOGFILE */

#ifdef NBQ
   !---------------------------------------------------------------------
   !  init_time_stepping_nbq
   !  Propagate dtfast into the NBQ common block.
   !  NOTE: the COMMON block is legacy; it should eventually be replaced
   !        by a proper module variable.
   !---------------------------------------------------------------------
   subroutine init_time_stepping_nbq()
      use scalars, ONLY: dtfast
      implicit none

      real :: dtnbq
      common/time_nbq2/dtnbq   ! TODO: replace with module variable

      dtnbq = dtfast
      ndtnbq = 1   ! reset: substep ratio is always 1 at init

   end subroutine init_time_stepping_nbq
#endif

#ifdef USE_CALENDAR
   !---------------------------------------------------------------------
   !  init_calendar
   !  Parse start_date string and populate ncscrum origin variables.
   !  Only active when nrrec == 0 (ANA_INITIAL).
   !---------------------------------------------------------------------
   subroutine init_calendar()
      use scalars, ONLY: nrrec, start_time
      use ncscrum, ONLY: origin_year, origin_month, origin_day, &
                         origin_hour, origin_minute, origin_second, &
                         origin_date, origin_date_in_sec
      implicit none

      real*8, external :: tool_datosec

#ifdef ANA_INITIAL
      if (nrrec == 0) then
         origin_date = start_date
         origin_date_in_sec = tool_datosec(origin_date)
         READ (origin_date(1:4), fmt='(i4)') origin_year
         READ (origin_date(6:7), fmt='(i2)') origin_month
         READ (origin_date(9:10), fmt='(i2)') origin_day
         READ (origin_date(12:13), fmt='(i2)') origin_hour
         READ (origin_date(15:16), fmt='(i2)') origin_minute
         READ (origin_date(18:19), fmt='(i2)') origin_second
         start_time = origin_date_in_sec
      end if
#endif

   end subroutine init_calendar
#endif


   !---------------------------------------------------------------------
   !  check_time_stepping
   !  Validate &croco_time_stepping parameters.
   !---------------------------------------------------------------------
subroutine check_time_stepping(ierr)
   use param, ONLY: stdout, NWEIGHT
#if defined MPI
   use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
   implicit none
   integer, intent(inout) :: ierr

   if (dt == 0.0) then
      MPI_master_only write (stdout, '(a,f10.1)') &
         'Error - Null baroclinic time step dt: ', dt
      ierr = ierr + 1
   end if

   if (ntimes == 0) then
      MPI_master_only write (stdout, '(a,i0)') &
         'Error - Null number of time steps ntimes: ', ntimes
      ierr = ierr + 1
   end if

   if (NWEIGHT < (2*ndtfast - 1)) then
      MPI_master_only write (stdout, '(a,i0)') &
         'Error - Number of 2D time steps (2*ndtfast-1): ', 2*ndtfast - 1
      MPI_master_only write (stdout, '(a,i0)') &
         '        exceeds barotropic weight dimension NWEIGHT: ', NWEIGHT
      ierr = ierr + 1
   end if

end subroutine check_time_stepping

!---------------------------------------------------------------------
!  check_history
!  Validate &croco_history parameters.
!---------------------------------------------------------------------
subroutine check_history(ierr)
   use param, ONLY: stdout
#if defined MPI
   use scalars, ONLY: mynode
#endif
   implicit none
   integer, intent(inout) :: ierr

   if (nwrt <= 0) then
      MPI_master_only write (stdout, '(a,i0)') &
         'Error - History output frequency nwrt must be > 0: ', nwrt
      ierr = ierr + 1
   end if

   if (nrpfhis < -1) then
      MPI_master_only write (stdout, '(a,i0)') &
         'Error - nrpfhis must be >= -1: ', nrpfhis
      ierr = ierr + 1
   end if

end subroutine check_history

#ifdef NBQ
!---------------------------------------------------------------------
!  check_time_stepping_nbq
!  Validate &croco_time_stepping_nbq parameters.
!---------------------------------------------------------------------
subroutine check_time_stepping_nbq(ierr)
   use param, ONLY: stdout
#  if defined MPI
   use scalars, ONLY: mynode
#  endif
   implicit none
   integer, intent(inout) :: ierr

   if (ndtnbq <= 0) then
      MPI_master_only write (stdout, '(a,i0)') &
         'Error - NBQ acoustic substep ratio ndtnbq must be > 0: ', ndtnbq
      ierr = ierr + 1
   end if

   if (csound_nbq > 1500.0) then
      MPI_master_only write (stdout, '(a,f12.4,a)') &
         'Error - NBQ pseudo-acoustic speed csound_nbq = ', csound_nbq, &
         ' must not exceed real acoustic speed (1500 m/s).'
      ierr = ierr + 1
   end if

   if (visc2_nbq < 0.0) then
      MPI_master_only write (stdout, '(a,f12.4,a)') &
         'Error - NBQ bulk viscosity visc2_nbq = ', visc2_nbq, &
         ' must be positive.'
      ierr = ierr + 1
   end if

   ! NOTE: the check csound_nbq > 5*sqrt(g*hmax) requires hmax which
   ! is not yet available at namelist-reading time. It should be moved
   ! to the model initialisation phase once hmax is known.

end subroutine check_time_stepping_nbq
#endif

END MODULE croco_namelist
