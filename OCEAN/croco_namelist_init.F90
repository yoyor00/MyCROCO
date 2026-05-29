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
   public :: init_restart
#ifndef ANA_GRID
   public :: init_grid
#endif
   public :: init_forcing
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
   public :: init_bulk_forcing
#endif
#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || \
   (defined ZCLIMATOLOGY  && !defined ANA_SSH) || \
   (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
   public :: init_climatology
#endif
#if !defined ANA_BRY && defined FRC_BRY
   public :: init_boundary
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
   public :: init_wkb_boundary
#endif
#ifdef AVERAGES
   public :: init_averages
#endif
#if defined ABL1D && !defined XIOS
   public :: init_abl
#  ifdef AVERAGES
   public :: init_abl_averages
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
   public :: init_wave_offline
#endif
#if defined BIOLOGY && defined PISCES
   public :: init_biology
#endif
#ifdef LOGFILE
   public :: init_logfile
#endif
#ifdef NBQ
   public :: init_time_stepping_nbq
#endif
#ifdef USE_CALENDAR
   public :: init_calendar
#endif

   integer, parameter :: testunit = 40 ! unit used to probe input file availability

contains

   !---------------------------------------------------------------------
   !  init_all
   !  Single entry point that runs every init_xxx in sequence.
   !  Guaranteed to be called with ierr == 0 (check_all passed), but
   !  filesystem errors can still occur in init_history / init_initial /
   !  init_restart / init_logfile, so ierr is checked before init_logfile.
   !---------------------------------------------------------------------
   subroutine init_all(ierr)
      use croco_namelist, ONLY: use_frcname
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
      call init_restart(ierr)
#ifndef ANA_GRID
      call init_grid(ierr)
#endif
      if (use_frcname) then
         call init_forcing(ierr)
      end if
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
      call init_bulk_forcing(ierr)
#endif
#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || \
   (defined ZCLIMATOLOGY  && !defined ANA_SSH) || \
   (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
      call init_climatology(ierr)
#endif
#if !defined ANA_BRY && defined FRC_BRY
      call init_boundary(ierr)
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      call init_wkb_boundary(ierr)
#endif
#ifdef AVERAGES
      call init_averages(ierr)
#endif
#if defined ABL1D && !defined XIOS
      call init_abl(ierr)
#  ifdef AVERAGES
      call init_abl_averages(ierr)
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      call init_wave_offline(ierr)
#endif
#if defined BIOLOGY && defined PISCES
      call init_biology(ierr)
#endif
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
      use croco_namelist, ONLY: hisname, nwrt
#ifdef USE_CALENDAR
      ! Two separate use statements on croco_namelist are intentional:
      ! a Fortran continuation line (&) cannot span a CPP directive,
      ! so the conditional dt_his/dt cannot be merged into the block above.
      use croco_namelist, ONLY: dt_his, dt
#endif
      implicit none
      integer, intent(inout) :: ierr

      call adjust_filename_parallel(hisname, "hisname", ierr)
      call adjust_filename_ensemble(hisname)
      if (ierr /= 0) return

#ifdef USE_CALENDAR
      nwrt = ceiling((dt_his*3600.0)/dt)
#endif

   end subroutine init_history

   !---------------------------------------------------------------------
   !  init_initial
   !  Adjust ininame for MPI/ENSEMBLE and check file availability.
   !  Skipped when ANA_INITIAL is defined and nrrec == 0 (initial
   !  conditions are generated internally; no file needed).
   !---------------------------------------------------------------------
   subroutine init_initial(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: ininame, nrrec
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

#ifdef ANA_INITIAL
      if (nrrec == 0) return
#endif

      call adjust_filename_parallel(ininame, "ininame", ierr)
      call adjust_filename_ensemble(ininame)
      if (ierr /= 0) return

      open (testunit, file=trim(ininame), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open initial file ', trim(ininame)
         ierr = ierr + 1
      end if

   end subroutine init_initial

   !---------------------------------------------------------------------
   !  init_restart
   !  Adjust rstname for MPI/ENSEMBLE and derive nrst from dt_rst
   !  when USE_CALENDAR is active.
   !---------------------------------------------------------------------
   subroutine init_restart(ierr)
      use croco_namelist, ONLY: rstname, nrst
#ifdef USE_CALENDAR
      use croco_namelist, ONLY: dt_rst, dt
#endif
      implicit none
      integer, intent(inout) :: ierr

      call adjust_filename_parallel(rstname, "rstname", ierr)
      call adjust_filename_ensemble(rstname)
      if (ierr /= 0) return

#ifdef USE_CALENDAR
      nrst = ceiling((dt_rst*3600.0)/dt)
#endif

   end subroutine init_restart

#ifndef ANA_GRID
   !---------------------------------------------------------------------
   !  init_grid
   !  Adjust grdname for MPI and check file availability.
   !---------------------------------------------------------------------
   subroutine init_grid(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: grdname
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      call adjust_filename_parallel(grdname, "grdname", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(grdname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open grid file ', trim(grdname)
         ierr = ierr + 1
      end if

   end subroutine init_grid
# endif

   !---------------------------------------------------------------------
   !  init_forcing
   !  Adjust frcname for MPI and check file availability.
   !---------------------------------------------------------------------
   subroutine init_forcing(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: frcname
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      call adjust_filename_parallel(frcname, "frcname", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(frcname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open forcing file ', trim(frcname)
         ierr = ierr + 1
      end if

   end subroutine init_forcing

#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
   !---------------------------------------------------------------------
   !  init_bulk_forcing
   !  Adjust bulkname for MPI and check file availability.
   !---------------------------------------------------------------------
   subroutine init_bulk_forcing(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: bulkname
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      call adjust_filename_parallel(bulkname, "bulkname", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(bulkname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open bulk forcing file ', trim(bulkname)
         ierr = ierr + 1
      end if

   end subroutine init_bulk_forcing
#endif

#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || \
   (defined ZCLIMATOLOGY  && !defined ANA_SSH) || \
   (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
   !---------------------------------------------------------------------
   !  init_climatology
   !  Adjust clmname for MPI and check file availability.
   !  For AGRIF child grids the climatology file is not used.
   !---------------------------------------------------------------------
   subroutine init_climatology(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: clmname
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

#  ifdef AGRIF
      if (.not. Agrif_Root()) return
#  endif

      call adjust_filename_parallel(clmname, "clmname", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(clmname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open climatology file ', trim(clmname)
         ierr = ierr + 1
      end if

   end subroutine init_climatology
#endif

#if defined WAVE_OFFLINE && defined MUSTANG
   !---------------------------------------------------------------------
   !  init_wave_offline
   !  Adjust wave_file for MPI and check file availability.
   !---------------------------------------------------------------------
   subroutine init_wave_offline(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: wave_file
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      call adjust_filename_parallel(wave_file, "wave_file", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(wave_file), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open wave_offline file ', trim(wave_file)
         ierr = ierr + 1
      end if

   end subroutine init_wave_offline
# endif

#if defined BIOLOGY && defined PISCES
   !---------------------------------------------------------------------
   !  init_biology
   !  Adjust bioname for MPI
   !---------------------------------------------------------------------
   subroutine init_biology(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: bioname
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      call adjust_filename_parallel(bioname, "bioname", ierr)
      if (ierr /= 0) return

      ! Notes : no check of availability here

   end subroutine init_biology
# endif

#if !defined ANA_BRY && defined FRC_BRY
   !---------------------------------------------------------------------
   !  init_boundary
   !  Adjust bry_file for MPI and check file availability.
   !  AGRIF child grids do not use a boundary file.
   !---------------------------------------------------------------------
   subroutine init_boundary(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: bry_file
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

#  ifdef AGRIF
      if (.not. Agrif_Root()) return
#  endif

      call adjust_filename_parallel(bry_file, "bry_file", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(bry_file), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open boundary file ', trim(bry_file)
         ierr = ierr + 1
      end if

   end subroutine init_boundary
#endif

#if defined WKB_WWAVE && !defined ANA_BRY_WKB
   !---------------------------------------------------------------------
   !  init_wkb_boundary
   !  Adjust brywkb_file for MPI and check file availability.
   !  AGRIF child grids do not use a WKB boundary file.
   !---------------------------------------------------------------------
   subroutine init_wkb_boundary(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: brywkb_file
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

#  ifdef AGRIF
      if (.not. Agrif_Root()) return
#  endif

      call adjust_filename_parallel(brywkb_file, "brywkb_file", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(brywkb_file), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open WKB boundary file ', trim(brywkb_file)
         ierr = ierr + 1
      end if

   end subroutine init_wkb_boundary
#endif

#ifdef AVERAGES
   !---------------------------------------------------------------------
   !  init_averages
   !  Adjust avgname for MPI/ENSEMBLE and derive navg from dt_avg
   !  when USE_CALENDAR is active.
   !---------------------------------------------------------------------
   subroutine init_averages(ierr)
      use croco_namelist, ONLY: avgname, navg
#  ifdef USE_CALENDAR
      use croco_namelist, ONLY: dt_avg, dt
#  endif
      implicit none
      integer, intent(inout) :: ierr

      call adjust_filename_parallel(avgname, "avgname", ierr)
      call adjust_filename_ensemble(avgname)
      if (ierr /= 0) return

#  ifdef USE_CALENDAR
      navg = ceiling((dt_avg*3600.0)/dt)
#  endif

   end subroutine init_averages
#endif

#if defined ABL1D && !defined XIOS
   !---------------------------------------------------------------------
   !  init_abl
   !  Adjust ablname for MPI/ENSEMBLE.
   !---------------------------------------------------------------------
   subroutine init_abl(ierr)
      use croco_namelist, ONLY: ablname
      implicit none
      integer, intent(inout) :: ierr

      call adjust_filename_parallel(ablname, "ablname", ierr)
      call adjust_filename_ensemble(ablname)

   end subroutine init_abl

#  ifdef AVERAGES
   !---------------------------------------------------------------------
   !  init_abl_averages
   !  Adjust ablname_avg for MPI/ENSEMBLE.
   !  Apply default nwrtablavg = navg when nwrtablavg is 0.
   !---------------------------------------------------------------------
   subroutine init_abl_averages(ierr)
      use croco_namelist, ONLY: ablname_avg, nwrtablavg, navg
      implicit none
      integer, intent(inout) :: ierr

      if (nwrtablavg == 0) nwrtablavg = navg

      call adjust_filename_parallel(ablname_avg, "ablname_avg", ierr)
      call adjust_filename_ensemble(ablname_avg)

   end subroutine init_abl_averages
#  endif
#endif

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

      call adjust_filename_parallel(logname, "logname", ierr)
      if (ierr /= 0) return
      call adjust_filename_ensemble(logname)

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

         open (unit=stdout, file=trim(logfile_path), &
               status='REPLACE', form='formatted', iostat=ios)
         if (ios /= 0) then
            write (6, '(/1x,3A/)') 'ERROR: Cannot open log file: ', trim(logfile_path)
            ierr = ierr + 1
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

   !---------------------------------------------------------------------
   !  adjust_filename_parallel  (private helper)
   !
   !  Apply all standard filename adjustments in one place:
   !    1. insert MPI node suffix  (only if MPI + PARALLEL_FILES)
   !    2. trim trailing blanks
   !    3. prepend ensemble member prefix  (only if ENSEMBLE)
   !
   !  context : caller name used in error messages (e.g. "hisname")
   !---------------------------------------------------------------------
   subroutine adjust_filename_parallel(fname, context, ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode, mynode2, NNODES2
#endif
      implicit none
      character(len=*), intent(inout) :: fname
      character(len=*), intent(in)    :: context
      integer, intent(inout) :: ierr

#if defined MPI && defined PARALLEL_FILES
      call insert_node(fname, len_trim(fname), mynode2, NNODES2, ierr)
      if (ierr /= 0) then
         MPI_master_only write (stdout, '(3a)') &
            'Error in adjust_filename_parallel: insert_node failed for ', &
            trim(context), '.'
         return
      end if
#endif
      fname = trim(fname)
   end subroutine adjust_filename_parallel

   !---------------------------------------------------------------------
   !  adjust_filename_ensemble (private helper)
   !---------------------------------------------------------------------
   subroutine adjust_filename_ensemble(fname)
      implicit none
#ifdef ENSEMBLE
#include "mpi_cpl.h"
#endif
      character(len=*), intent(inout) :: fname

      fname = trim(fname)
#ifdef ENSEMBLE
      fname = cmember//fname
#endif

   end subroutine adjust_filename_ensemble

END MODULE croco_namelist_init
