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
!    fatal_nml_error    – print summary abort message (does NOT close
!                         the file; the caller owns the close)
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
   !    - all presence/parse errors are caught before any check runs
   !    - all value errors are reported before any init runs
   !    - each phase remains independently readable / testable
   !---------------------------------------------------------------------
   subroutine read_nml(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: fname_nml, &
                                title, logname, &
                                dt, ntimes, ndtfast, ninfo, &
                                ldefhis, nwrt, nrpfhis, hisname, &
                                nrrec, ininame, &
                                nrst, nrpfrst, rstname, &
                                use_frcname, frcname
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
#ifndef ANA_GRID
      use croco_namelist, ONLY: grdname
#endif
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
      use croco_namelist, ONLY: bulkname
#endif
#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY( & !defined ANA_SSH) || \
       defined M2CLIMATOLOGY(d & !defined ANA_M2CLIMA) || \
                             efined M3CLIMATOLOGY use & !defined ANA_M3CLIMA)
                             croco_namelist, ONLY:clmname
#endif
#if !defined ANA_BRY && defined FRC_BRY
      use croco_namelist, ONLY: bry_file
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      use croco_namelist, ONLY: brywkb_file
#endif
#ifdef AVERAGES
      use croco_namelist, ONLY: ntsavg, navg, nrpfavg, avgname
#endif
#if defined ABL1D && !defined XIOS
      use croco_namelist, ONLY: ldefablhis, nwrtablhis, nrpfablhis, ablname
#  ifdef AVERAGES
      use croco_namelist, ONLY: ldefablavg, ntsablavg, nwrtablavg, nrpfablavg, &
                                ablname_avg
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      use croco_namelist, ONLY: wave_file
#endif
#if defined BIOLOGY && defined PISCES
      use croco_namelist, ONLY: bioname
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
      namelist /croco_initial/ nrrec, ininame
      namelist /croco_restart/ nrst, nrpfrst, rstname
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
#ifndef ANA_GRID
      namelist /croco_grid/ grdname
#endif
      namelist /croco_forcing/ frcname
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
      namelist /croco_bulk_forcing/ bulkname
#endif
#if (  defined TCLIMATOLOGY  && \
      !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY&&\
      !defined ANA_SSH) || \
      (defined M2CLIMATOLOGY&&\
      !defined ANA_M2CLIMA) || \
      (defined M3CLIMATOLOGY&&\
      !defined ANA_M3CLIMA)
      namelist /croco_climatology/ clmname
#endif
#if !defined ANA_BRY && defined FRC_BRY
      namelist /croco_boundary/ bry_file
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      namelist /croco_wkb_boundary/ brywkb_file
#endif
#ifdef AVERAGES
      namelist /croco_averages/ ntsavg, navg, nrpfavg, avgname
#endif
#if defined ABL1D && !defined XIOS
      namelist /croco_abl/ ldefablhis, nwrtablhis, nrpfablhis, ablname
#  ifdef AVERAGES
      namelist /croco_abl_averages/ ldefablavg, ntsablavg, nwrtablavg, &
         nrpfablavg, ablname_avg
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      namelist /croco_wave_offline/ wave_file
#endif
#if defined BIOLOGY && defined PISCES
      namelist /croco_biology/ bioname
#endif
      ierr = 0

      MPI_master_only write (stdout, *) '*** READING NAMELIST ***'

      ! Open the namelist file; if absent, fall back to default values
      ! for all parameters (non-fatal: the model can still run).
      ! newunit= lets the compiler pick a free unit number automatically.
      ! open (newunit=nmlunit, file=trim(fname_nml), status='old', &
      !       action='read', iostat=ios)
      ! AGRIF conv incompatibility with newunit for now...
      ! TODO, update AGRIF conv
      nmlunit = 10
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
      ! Two distinct error strategies:
      !   - presence errors (mandatory section absent): accumulated across
      !     all sections, file closed once at end of phase, single return.
      !   - parse errors (malformed section): immediate close + return,
      !     because continuing to read a corrupt file is meaningless.
      ! In both cases the file is closed by read_nml.
      ! ----------------------------------------------------------------

      ! --- croco_title (optional) ---
      call check_nml_presence(nmlunit, "croco_title", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_title, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_title (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_time_stepping (mandatory) ---
      call check_nml_presence(nmlunit, "croco_time_stepping", .true., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_time_stepping, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_time_stepping (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

#ifdef NBQ
      ! --- croco_time_stepping_nbq (mandatory if NBQ) ---
      call check_nml_presence(nmlunit, "croco_time_stepping_nbq", .true., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_time_stepping_nbq, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_time_stepping_nbq (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif /* NBQ */

#ifdef USE_CALENDAR
      ! --- croco_use_calendar (mandatory if USE_CALENDAR) ---
      call check_nml_presence(nmlunit, "croco_use_calendar", .true., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_use_calendar, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_use_calendar (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif /* USE_CALENDAR */

#ifdef SOLVE3D
      ! --- croco_s_coord (optional if SOLVE3D) ---
      call check_nml_presence(nmlunit, "croco_s_coord", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_s_coord, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_s_coord (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif /* SOLVE3D */

      ! --- croco_history (optional) ---
      call check_nml_presence(nmlunit, "croco_history", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_history, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_history (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_initial (optional) ---
      call check_nml_presence(nmlunit, "croco_initial", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_initial, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_initial (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_restart (optional) ---
      call check_nml_presence(nmlunit, "croco_restart", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_restart, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_restart (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

#ifndef ANA_GRID
      ! --- croco_grid (optional) ---
      call check_nml_presence(nmlunit, "croco_grid", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_grid, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_grid (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif /* ANA_GRID */

      if (use_frcname) then
         ! --- croco_forcing (optional) ---
         call check_nml_presence(nmlunit, "croco_forcing", .false., found, ierr)
         if (found) then
            read (nmlunit, nml=croco_forcing, iostat=ios); rewind (nmlunit)
            if (ios /= 0) then
               call fatal_nml_error("croco_forcing (parse error)")
               ierr = ierr + 1; close (nmlunit); return
            end if
         end if
      end if

#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
      ! --- croco_bulk_forcing (optional) ---
      call check_nml_presence(nmlunit, "croco_bulk_forcing", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_bulk_forcing, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_bulk_forcing (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if (  defined TCLIMATOLOGY  && \
      !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY&&\
      !defined ANA_SSH) || \
      (defined M2CLIMATOLOGY&&\
      !defined ANA_M2CLIMA) || \
      (defined M3CLIMATOLOGY&&\
      !defined ANA_M3CLIMA)
      ! --- croco_climatology (optional) ---
      call check_nml_presence(nmlunit, "croco_climatology", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_climatology, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_climatology (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if defined WAVE_OFFLINE && defined MUSTANG
      ! --- croco_wave_offline (optional) ---
      call check_nml_presence(nmlunit, "croco_wave_offline", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_wave_offline, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_wave_offline (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if defined BIOLOGY && defined PISCES
      ! --- croco_biology (optional) ---
      call check_nml_presence(nmlunit, "croco_biology", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_biology, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_biology (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if !defined ANA_BRY && defined FRC_BRY
      ! --- croco_boundary (optional) ---
      call check_nml_presence(nmlunit, "croco_boundary", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_boundary, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_boundary (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      ! --- croco_wkb_boundary (optional) ---
      call check_nml_presence(nmlunit, "croco_wkb_boundary", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_wkb_boundary, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_wkb_boundary (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef AVERAGES
      ! --- croco_averages (optional) ---
      call check_nml_presence(nmlunit, "croco_averages", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_averages, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_averages (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if defined ABL1D && !defined XIOS
      ! --- croco_abl (optional) ---
      call check_nml_presence(nmlunit, "croco_abl", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_abl, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_abl (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_abl_averages (optional) ---
      call check_nml_presence(nmlunit, "croco_abl_averages", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_abl_averages, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_abl_averages (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif

      ! --- croco_logfile (optional) ---
#ifdef LOGFILE
      call check_nml_presence(nmlunit, "croco_logfile", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_logfile, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_logfile (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif /* LOGFILE */

      ! Single close: the file is no longer needed after this point.
      close (nmlunit)

      ! Presence errors (mandatory sections missing) are fatal: stop here
      ! so that check_all and init_all never run on incomplete data.
      if (ierr /= 0) return

      ! ----------------------------------------------------------------
      ! Phase 2 – validate all parameter values
      ! All check_xxx routines run regardless of individual failures so
      ! the user sees the complete error list in one pass.
      ! ----------------------------------------------------------------
      call check_all(ierr)
      if (ierr /= 0) return

      ! ----------------------------------------------------------------
      ! Phase 3 – compute derived quantities
      ! init logfile change stdout direction
      ! ----------------------------------------------------------------
      call init_all(ierr)
      if (ierr /= 0) return

      ! ----------------------------------------------------------------
      ! Phase 4 – write nml in log
      ! ----------------------------------------------------------------
      ! Echo all namelists to stdout so the user can verify what was read
      ! Note : can not do that in subroutine due to incompatibility of
      ! AGRIF conv with a namelist declaration at module level
      MPI_master_only WRITE (stdout, *) "---------------------"
      MPI_master_only WRITE (stdout, *) "CROCO namelist :"
      MPI_master_only WRITE (stdout, *)
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
      MPI_master_only WRITE (stdout, nml=croco_initial)
      MPI_master_only WRITE (stdout, nml=croco_restart)
#ifndef ANA_GRID
      MPI_master_only WRITE (stdout, nml=croco_grid)
#endif
      if (use_frcname) then
         MPI_master_only WRITE (stdout, nml=croco_forcing)
      end if
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
      MPI_master_only WRITE (stdout, nml=croco_bulk_forcing)
#endif
#if (  defined TCLIMATOLOGY  && \
      !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY&&\
      !defined ANA_SSH) || \
      (defined M2CLIMATOLOGY&&\
      !defined ANA_M2CLIMA) || \
      (defined M3CLIMATOLOGY&&\
      !defined ANA_M3CLIMA)
      MPI_master_only WRITE (stdout, nml=croco_climatology)
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      MPI_master_only WRITE (stdout, nml=croco_wave_offline)
#endif
#if defined BIOLOGY && defined PISCES
      MPI_master_only WRITE (stdout, nml=croco_biology)
#endif
#if !defined ANA_BRY && defined FRC_BRY
      MPI_master_only WRITE (stdout, nml=croco_boundary)
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      MPI_master_only WRITE (stdout, nml=croco_wkb_boundary)
#endif
#ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_averages)
#endif
#if defined ABL1D && !defined XIOS
      MPI_master_only WRITE (stdout, nml=croco_abl)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_abl_averages)
#  endif
#endif

      MPI_master_only WRITE (stdout, *) "End of CROCO namelist"
      MPI_master_only WRITE (stdout, *) "---------------------"
      MPI_master_only WRITE (stdout, *)

   end subroutine read_nml

   !---------------------------------------------------------------------
   !  check_nml_presence  (private helper)
   !
   !  Scan the namelist file line by line looking for "&nml_name"
   !  (case-insensitive).  Rewinds before and after the search.
   !  Does NOT close the file on error; the caller owns the close.
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
   !
   !  Print a single "aborting" line.
   !  Does NOT close the namelist file: the caller is responsible for
   !  closing it, which avoids any risk of double-close.
   !---------------------------------------------------------------------
   subroutine fatal_nml_error(nml_name)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      character(len=*), intent(in) :: nml_name

      MPI_master_only write (stdout, '(3a)') &
         'FATAL: aborting after errors in "', trim(nml_name), '" namelist.'

   end subroutine fatal_nml_error

END MODULE croco_namelist_read
