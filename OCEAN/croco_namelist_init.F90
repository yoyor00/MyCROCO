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
!  MODULE croco_namelist_init
!
!  Purpose : Initialisation routines derived from namelist parameters.
!            Called only after check_all has returned ierr == 0.
!
!  Public entry points:
!    init_all               – call every init_xxx in order
!
!  Individual init routines are kept public for unit-testing purposes.
!=======================================================================

MODULE croco_namelist_init
   implicit none
   private

   public :: init_all
   public :: init_time_stepping
   public :: init_history
   public :: init_initial
#ifdef LOGFILE
   public :: init_logfile
#endif
#ifdef NBQ
   public :: init_time_stepping_nbq
#endif
#ifdef USE_CALENDAR
   public :: init_calendar
#endif

   integer :: testunit = 40 ! to test availability of input files

contains

   !---------------------------------------------------------------------
   !  init_all
   !  Single entry point that runs every init_xxx in sequence.
   !  Guaranteed to be called with ierr == 0 (check_all passed), but
   !  init_history and init_logfile can still fail on system errors
   !  (e.g. filesystem issues), so ierr is checked between them.
   !---------------------------------------------------------------------
   subroutine init_all(ierr)
      implicit none
      integer, intent(out) :: ierr

      ierr = 0
      call init_time_stepping()
#ifdef NBQ
      call init_time_stepping_nbq()
#endif
#ifdef USE_CALENDAR
      call init_calendar()
#endif
      call init_history(ierr)
      call init_initial(ierr)
#ifdef LOGFILE
      if (ierr == 0) call init_logfile(ierr)
#endif

   end subroutine init_all

   !---------------------------------------------------------------------
   !  init_time_stepping
   !  Derive dtfast from dt and ndtfast.
   !---------------------------------------------------------------------
   subroutine init_time_stepping()
      use croco_namelist, ONLY: dt, ndtfast
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
      use croco_namelist, ONLY: hisname, nwrt
#ifdef USE_CALENDAR
      use croco_namelist, ONLY: dt_his, dt
#endif
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

   !---------------------------------------------------------------------
   !  init_initial
   !  Adjust ininame for MPI/ENSEMBLE
   !  Check its availability (in the case
   !  of analytical initial conditions and nrrec=0 initial conditions are
   !  created internally and no file is needed).
   !---------------------------------------------------------------------
   subroutine init_initial(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: ininame, nrrec
#if defined MPI
      use scalars, ONLY: mynode, mynode2, NNODES2
#endif
      implicit none
#ifdef ENSEMBLE
#include "mpi_cpl.h"
#endif
      integer, intent(inout) :: ierr
      integer :: ios

#ifdef ANA_INITIAL
        if (nrrec.gt.0) then
#endif

#if defined MPI && defined PARALLEL_FILES
      call insert_node(ininame, len_trim(ininame), mynode2, NNODES2, ierr)
      if (ierr /= 0) then
         MPI_master_only write (stdout, *) &
            'Error in init_initial: insert_node failed for ininame'
         return
      end if
#endif

      ininame = trim(ininame)

#ifdef ENSEMBLE
      ininame = cmember//ininame
#endif

   open(testunit, file=trim(ininame), status='old', iostat=ios)

   if (ios == 0) then
      close(testunit)
   else
      MPI_master_only write(stdout, *) &
         'Error: cannot open initial file ', trim(ininame)
      return
   endif
#ifdef ANA_INITIAL
   endif
#endif

   end subroutine init_initial

#ifdef LOGFILE
   !---------------------------------------------------------------------
   !  init_logfile
   !  Open the log file and redirect stdout to it.
   !---------------------------------------------------------------------
   subroutine init_logfile(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: logname
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
   !
   !  NOTE 1: ndtnbq is intentionally forced to 1 here regardless of the
   !          namelist value.  The substep ratio is always 1 at startup;
   !          the namelist variable is reserved for future dynamic use.
   !          This behaviour is documented in the namelist declaration.
   !  NOTE 2: the COMMON block is legacy; it should eventually be replaced
   !          by a proper module variable.
   !---------------------------------------------------------------------
   subroutine init_time_stepping_nbq()
      use croco_namelist, ONLY: ndtnbq
      use scalars, ONLY: dtfast
      implicit none

      real :: dtnbq
      common/time_nbq2/dtnbq   ! TODO: replace with module variable

      dtnbq = dtfast
      ndtnbq = 1   ! reset: substep ratio is always 1 at init

   end subroutine init_time_stepping_nbq
#endif /* NBQ */

#ifdef USE_CALENDAR
   !---------------------------------------------------------------------
   !  init_calendar
   !  Parse start_date string and populate ncscrum origin variables.
   !  Only active when nrrec == 0 (ANA_INITIAL).
   !---------------------------------------------------------------------
   subroutine init_calendar()
      use croco_namelist, ONLY: start_date, nrrec
      use scalars, ONLY: start_time
      use ncscrum, ONLY: origin_year, origin_month, origin_day, &
                         origin_hour, origin_minute, origin_second, &
                         origin_date, origin_date_in_sec
      implicit none

      real(kind=8), external :: tool_datosec

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
#endif /* USE_CALENDAR */

END MODULE croco_namelist_init
