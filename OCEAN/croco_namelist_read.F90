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
!  MODULE croco_namelist_read
!
!  Purpose : Open the namelist file, read all sections, then call
!            check_all and init_all.
!
!  Public entry points:
!    read_nml_fname  – set fname_nml from command-line argument 2
!    read_nml        – full read / check / init sequence
!
!  Internal helpers:
!    check_nml_presence – locate a section header in the file
!    fatal_nml_error    – print summary abort line and close file
!=======================================================================

MODULE croco_namelist_read
   implicit none
   private

   public :: read_nml_fname
   public :: read_nml

contains

   !---------------------------------------------------------------------
   !  read_nml_fname
   !  Set fname_nml from command-line argument 2 (defaults to croco.nml).
   !---------------------------------------------------------------------
   subroutine read_nml_fname()
      use croco_namelist, ONLY: fname_nml
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
   !  Phase 1 : read all namelist sections from file.
   !  Phase 2 : validate all parameters (check_all).
   !  Phase 3 : compute derived quantities (init_all).
   !
   !  Separating the three phases means:
   !    - all parse errors are caught before any check runs
   !    - all value errors are reported before any init runs
   !    - each phase remains independently readable / testable
   !---------------------------------------------------------------------
   subroutine read_nml(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: fname_nml, &
                                title, logname, &
                                dt, ntimes, ndtfast, ninfo, &
                                ldefhis, nwrt, nrpfhis, hisname
      use croco_namelist_check, ONLY: check_all
      use croco_namelist_init, ONLY: init_all
#ifdef NBQ
      use croco_namelist, ONLY: ndtnbq, csound_nbq, visc2_nbq
#endif
#ifdef SOLVE3D
      use croco_namelist, ONLY: theta_s, theta_b, Tcline
#endif
#ifdef USE_CALENDAR
      use croco_namelist, ONLY: start_date, end_date, &
                                dt_his, dt_avg, dt_rst
#endif
#if defined MPI
      use scalars, ONLY: mynode
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
#ifdef SOLVE3D
      namelist /croco_s_coord/ theta_s, theta_b, Tcline
#endif
#ifdef USE_CALENDAR
      namelist /croco_use_calendar/ start_date, end_date, &
         dt_his, dt_avg, dt_rst
#endif

      ierr = 0
      nmlunit = 10

      MPI_master_only write (stdout, *) '*** READING NAMELIST ***'

      ! Open the namelist file; if absent, fall back to default values
      ! for all parameters (non-fatal: the model can still run).
      open (unit=nmlunit, file=trim(fname_nml), status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         MPI_master_only write (stdout, *) &
            'WARNING: namelist file not found: ', trim(fname_nml)
         MPI_master_only write (stdout, *) 'Using default values.'
         return
      end if

      ! ----------------------------------------------------------------
      ! Phase 1 – read all sections
      ! ----------------------------------------------------------------

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
#endif /* NBQ */

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
#endif /* USE_CALENDAR */

#ifdef SOLVE3D
      ! --- croco_s_coord (optional if SOLVE3D) ---
      call check_nml_presence(nmlunit, "croco_s_coord", .false., found, ierr)
      if (ierr /= 0) return
      if (found) then
         read (nmlunit, nml=croco_s_coord, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_s_coord (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if
#endif /* SOLVE3D */

      ! --- croco_history (optional) ---
      call check_nml_presence(nmlunit, "croco_history", .false., found, ierr)
      if (ierr /= 0) return
      if (found) then
         read (nmlunit, nml=croco_history, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_history (parse error)", nmlunit)
            ierr = ierr + 1; return
         end if
      end if

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
#endif /* LOGFILE */

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
#ifdef SOLVE3D
      MPI_master_only WRITE (stdout, nml=croco_s_coord)
#endif
#ifdef USE_CALENDAR
      MPI_master_only WRITE (stdout, nml=croco_use_calendar)
#endif
      MPI_master_only WRITE (stdout, nml=croco_history)

      ! ----------------------------------------------------------------
      ! Phase 2 – validate all parameter values
      ! ----------------------------------------------------------------
      call check_all(ierr)
      if (ierr /= 0) return

      ! ----------------------------------------------------------------
      ! Phase 3 – compute derived quantities
      ! ----------------------------------------------------------------
      call init_all(ierr)

   end subroutine read_nml

   !---------------------------------------------------------------------
   !  check_nml_presence  (private helper)
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

      found = .false.
      rewind (nmlunit)
      do
         read (nmlunit, '(A)', iostat=ios) line
         if (ios /= 0) exit
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
   !  fatal_nml_error  (private helper)
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

END MODULE croco_namelist_read
