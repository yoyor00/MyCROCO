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
      use param, ONLY: stdout, NT
      use croco_namelist
      use croco_namelist_check, ONLY: check_all
      use croco_namelist_init, ONLY: init_all
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
#ifdef WKB_WWAVE
      namelist /croco_wkb_wwave/ wkb_amp, wkb_ang, wkb_prd, wkb_tide, wkb_btg, wkb_gam
#  ifdef WAVE_ROLLER
      namelist /croco_wkb_roller/ wkb_rsb, wkb_roller
#  endif
#endif
#ifdef WAVE_MAKER
      namelist /croco_wave_maker/ wmaker_amp, wmaker_prd, wmaker_dir, wmaker_dsp, wmaker_fsp
#endif
#ifdef AVERAGES
      namelist /croco_averages/ ntsavg, navg, nrpfavg, avgname
#endif
#if defined OUTPUTS_SURFACE && !defined XIOS
      namelist /croco_surf/ ldefsurf, nwrtsurf, nrpfsurf, surfname
#  ifdef AVERAGES
      namelist /croco_surf_avg/ ldefsurf_avg, ntssurf_avg, nwrtsurf_avg, &
         nrpfsurf_avg, surfname_avg
#  endif
#endif
#if defined ABL1D && !defined XIOS
      namelist /croco_abl/ ldefablhis, nwrtablhis, nrpfablhis, ablname
#  ifdef AVERAGES
      namelist /croco_abl_averages/ ldefablavg, ntsablavg, nwrtablavg, &
         nrpfablavg, ablname_avg
#  endif
#endif

#if defined DIAGNOSTICS_TS
      namelist /croco_diagnostics_ts/ ldefdia, nwrtdia, nrpfdia, dianame
#  ifdef AVERAGES
      namelist /croco_diag_avg/ ldefdia_avg, ntsdia_avg, nwrtdia_avg, &
         nrpfdia_avg, dianame_avg
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
      namelist /croco_diag_mld_dens/ mld_crit_D, mld_crit_T
#  endif
#endif
#if defined DIAGNOSTICS_UV
      namelist /croco_diagnosticsM/ ldefdiaM, nwrtdiaM, nrpfdiaM, dianameM
#  ifdef AVERAGES
      namelist /croco_diagM_avg/ ldefdiaM_avg, ntsdiaM_avg, nwrtdiaM_avg, &
         nrpfdiaM_avg, dianameM_avg
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      namelist /croco_diags_vrt/ ldefdiags_vrt, nwrtdiags_vrt, nrpfdiags_vrt, &
         diags_vrtname
#  ifdef AVERAGES
      namelist /croco_diags_vrt_avg/ ldefdiags_vrt_avg, ntsdiags_vrt_avg, &
         nwrtdiags_vrt_avg, nrpfdiags_vrt_avg, diags_vrtname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_EK
      namelist /croco_diags_ek/ ldefdiags_ek, nwrtdiags_ek, nrpfdiags_ek, &
         diags_ekname
#  ifdef AVERAGES
      namelist /croco_diags_ek_avg/ ldefdiags_ek_avg, ntsdiags_ek_avg, &
         nwrtdiags_ek_avg, nrpfdiags_ek_avg, diags_ekname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      namelist /croco_diags_pv/ ldefdiags_pv, nwrtdiags_pv, nrpfdiags_pv, &
         diags_pvname
#  ifdef AVERAGES
      namelist /croco_diags_pv_avg/ ldefdiags_pv_avg, ntsdiags_pv_avg, &
         nwrtdiags_pv_avg, nrpfdiags_pv_avg, diags_pvname_avg
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
      namelist /croco_diags_eddy/ ldefdiags_eddy, nwrtdiags_eddy, nrpfdiags_eddy, &
         diags_eddyname
#  ifdef AVERAGES
      namelist /croco_diags_eddy_avg/ ldefdiags_eddy_avg, ntsdiags_eddy_avg, &
         nwrtdiags_eddy_avg, nrpfdiags_eddy_avg, diags_eddyname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      namelist /croco_diagnostics_bio/ ldefdiabio, nwrtdiabio, nrpfdiabio, dianamebio
#  ifdef AVERAGES
      namelist /croco_diagbio_avg/ ldefdiabio_avg, ntsdiabio_avg, nwrtdiabio_avg, &
         nrpfdiabio_avg, dianamebio_avg
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      namelist /croco_wave_offline/ wave_file
#endif
#if defined BIOLOGY && defined PISCES
      namelist /croco_biology/ bioname
#endif
#ifdef BODYFORCE
      namelist /croco_bodyforce/ levsfrc, levbfrc
#endif
#if !defined NONLIN_EOS
      namelist /croco_lin_eos/ R0, T0, S0, Tcoef, Scoef
#endif
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
      namelist /croco_abl_nudg_tra/ ltra_min, ltra_max
#endif
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
      namelist /croco_abl_nudg_dyn/ ldyn_min, ldyn_max
#endif
#ifdef SEDIMENT
      namelist /croco_sediments/ sedname
#endif
#ifdef MUSTANG
      namelist /croco_sediments_mustang/ sedname_must
#endif
#ifdef SUBSTANCE
      namelist /croco_substance/ subsfilename
#endif
#ifdef OBSTRUCTION
      namelist /croco_obstruction/ obstname
#endif
#ifdef XIOS
      namelist /croco_xios_origin_date/ xios_origin_date
#endif
#ifdef ASSIMILATION
      namelist /croco_assimilation/ aparnam, assname
#endif
      namelist /croco_rho0/ rho0
      namelist /croco_bottom_drag/ rdrg, rdrg2, Zobt, Cdb_min, Cdb_max
      namelist /croco_gamma2/ gamma2
      namelist /croco_lateral_visc/ visc2, visc4
#ifdef SOLVE3D
#  ifdef TRACERS
      namelist /croco_tracer_diff2/ tnu2
      namelist /croco_tracer_diff4/ tnu4
#  endif
#  if !defined LMD_MIXING && !defined GLS_MIXING
      namelist /croco_vertical_mixing/ Akv_bak, Akt_bak
#  endif
#endif
      ierr = 0

      ! Allocate tracer arrays with the actual NT size (from use param).
#ifdef SOLVE3D
#  ifdef TRACERS
      if (.not. allocated(tnu2))    then; allocate(tnu2(NT));    tnu2    = 0.0;   endif
      if (.not. allocated(tnu4))    then; allocate(tnu4(NT));    tnu4    = 0.0;   endif
#    if !defined LMD_MIXING && !defined GLS_MIXING
      if (.not. allocated(Akt_bak)) then; allocate(Akt_bak(NT)); Akt_bak = 1.e-6; endif
#    endif
#  endif
#endif

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

#ifdef WKB_WWAVE
      ! --- croco_wkb_wwave (optional) ---
      call check_nml_presence(nmlunit, "croco_wkb_wwave", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_wkb_wwave, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_wkb_wwave (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef WAVE_ROLLER
      ! --- croco_wkb_roller (optional) ---
      call check_nml_presence(nmlunit, "croco_wkb_roller", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_wkb_roller, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_wkb_roller (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif

#ifdef WAVE_MAKER
      ! --- croco_wave_maker (optional) ---
      call check_nml_presence(nmlunit, "croco_wave_maker", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_wave_maker, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_wave_maker (parse error)")
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

#if defined OUTPUTS_SURFACE && !defined XIOS
      ! --- croco_surf (optional) ---
      call check_nml_presence(nmlunit, "croco_surf", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_surf, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_surf (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_surf_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_surf_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_surf_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_surf_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
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

#if defined DIAGNOSTICS_TS
      ! --- croco_diagnostics_ts (optional) ---
      call check_nml_presence(nmlunit, "croco_diagnostics_ts", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diagnostics_ts, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diagnostics_ts (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diag_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diag_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diag_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diag_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
      ! --- croco_diag_mld_dens (optional) ---
      call check_nml_presence(nmlunit, "croco_diag_mld_dens", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diag_mld_dens, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diag_mld_dens (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#if defined DIAGNOSTICS_UV
      ! --- croco_diagnosticsM (optional) ---
      call check_nml_presence(nmlunit, "croco_diagnosticsM", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diagnosticsM, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diagnosticsM (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diagM_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diagM_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diagM_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diagM_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      ! --- croco_diags_vrt (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_vrt", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_vrt, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_vrt (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diags_vrt_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_vrt_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_vrt_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_vrt_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#ifdef DIAGNOSTICS_EK
      ! --- croco_diags_ek (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_ek", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_ek, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_ek (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diags_ek_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_ek_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_ek_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_ek_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      ! --- croco_diags_pv (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_pv", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_pv, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_pv (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diags_pv_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_pv_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_pv_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_pv_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
      ! --- croco_diags_eddy (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_eddy", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_eddy, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_eddy (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diags_eddy_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diags_eddy_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diags_eddy_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diags_eddy_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      ! --- croco_diagnostics_bio (optional) ---
      call check_nml_presence(nmlunit, "croco_diagnostics_bio", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diagnostics_bio, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diagnostics_bio (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  ifdef AVERAGES
      ! --- croco_diagbio_avg (optional) ---
      call check_nml_presence(nmlunit, "croco_diagbio_avg", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_diagbio_avg, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_diagbio_avg (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif
#ifdef BODYFORCE
      ! --- croco_bodyforce (optional) ---
      call check_nml_presence(nmlunit, "croco_bodyforce", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_bodyforce, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_bodyforce (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if !defined NONLIN_EOS
      ! --- croco_lin_eos (optional) ---
      call check_nml_presence(nmlunit, "croco_lin_eos", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_lin_eos, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_lin_eos (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
      ! --- croco_abl_nudg_tra (optional) ---
      call check_nml_presence(nmlunit, "croco_abl_nudg_tra", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_abl_nudg_tra, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_abl_nudg_tra (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
      ! --- croco_abl_nudg_dyn (optional) ---
      call check_nml_presence(nmlunit, "croco_abl_nudg_dyn", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_abl_nudg_dyn, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_abl_nudg_dyn (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef SEDIMENT
      ! --- croco_sediments (optional) ---
      call check_nml_presence(nmlunit, "croco_sediments", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_sediments, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_sediments (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef MUSTANG
      ! --- croco_sediments_mustang (optional) ---
      call check_nml_presence(nmlunit, "croco_sediments_mustang", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_sediments_mustang, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_sediments_mustang (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef SUBSTANCE
      ! --- croco_substance (optional) ---
      call check_nml_presence(nmlunit, "croco_substance", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_substance, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_substance (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef OBSTRUCTION
      ! --- croco_obstruction (optional) ---
      call check_nml_presence(nmlunit, "croco_obstruction", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_obstruction, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_obstruction (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef XIOS
      ! --- croco_xios_origin_date (optional) ---
      call check_nml_presence(nmlunit, "croco_xios_origin_date", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_xios_origin_date, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_xios_origin_date (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

#ifdef ASSIMILATION
      ! --- croco_assimilation (optional) ---
      call check_nml_presence(nmlunit, "croco_assimilation", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_assimilation, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_assimilation (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#endif

      ! --- croco_rho0 (optional) ---
      call check_nml_presence(nmlunit, "croco_rho0", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_rho0, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_rho0 (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_bottom_drag (optional) ---
      call check_nml_presence(nmlunit, "croco_bottom_drag", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_bottom_drag, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_bottom_drag (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_gamma2 (optional) ---
      call check_nml_presence(nmlunit, "croco_gamma2", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_gamma2, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_gamma2 (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_lateral_visc (optional) ---
      call check_nml_presence(nmlunit, "croco_lateral_visc", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_lateral_visc, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_lateral_visc (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

#ifdef SOLVE3D
#  ifdef TRACERS
      ! --- croco_tracer_diff2 (optional) ---
      call check_nml_presence(nmlunit, "croco_tracer_diff2", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_tracer_diff2, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_tracer_diff2 (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if

      ! --- croco_tracer_diff4 (optional) ---
      call check_nml_presence(nmlunit, "croco_tracer_diff4", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_tracer_diff4, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_tracer_diff4 (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif

#  if !defined LMD_MIXING && !defined GLS_MIXING
      ! --- croco_vertical_mixing (optional) ---
      call check_nml_presence(nmlunit, "croco_vertical_mixing", .false., found, ierr)
      if (found) then
         read (nmlunit, nml=croco_vertical_mixing, iostat=ios); rewind (nmlunit)
         if (ios /= 0) then
            call fatal_nml_error("croco_vertical_mixing (parse error)")
            ierr = ierr + 1; close (nmlunit); return
         end if
      end if
#  endif
#endif

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
#ifdef WKB_WWAVE
      MPI_master_only WRITE (stdout, nml=croco_wkb_wwave)
#  ifdef WAVE_ROLLER
      MPI_master_only WRITE (stdout, nml=croco_wkb_roller)
#  endif
#endif
#ifdef WAVE_MAKER
      MPI_master_only WRITE (stdout, nml=croco_wave_maker)
#endif
#ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_averages)
#endif
#if defined OUTPUTS_SURFACE && !defined XIOS
      MPI_master_only WRITE (stdout, nml=croco_surf)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_surf_avg)
#  endif
#endif
#if defined ABL1D && !defined XIOS
      MPI_master_only WRITE (stdout, nml=croco_abl)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_abl_averages)
#  endif
#  if defined ABL_NUDGING && defined ABL_NUDGING_TRA
      MPI_master_only WRITE (stdout, nml=croco_abl_nudg_tra)
#  endif
#  if defined ABL_NUDGING && defined ABL_NUDGING_DYN
      MPI_master_only WRITE (stdout, nml=croco_abl_nudg_dyn)
#  endif
#endif
#if defined DIAGNOSTICS_TS
      MPI_master_only WRITE (stdout, nml=croco_diagnostics_ts)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diag_avg)
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
      MPI_master_only WRITE (stdout, nml=croco_diag_mld_dens)
#  endif
#endif
#if defined DIAGNOSTICS_UV
      MPI_master_only WRITE (stdout, nml=croco_diagnosticsM)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diagM_avg)
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      MPI_master_only WRITE (stdout, nml=croco_diags_vrt)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diags_vrt_avg)
#  endif
#endif
#ifdef DIAGNOSTICS_EK
      MPI_master_only WRITE (stdout, nml=croco_diags_ek)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diags_ek_avg)
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      MPI_master_only WRITE (stdout, nml=croco_diags_pv)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diags_pv_avg)
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
      MPI_master_only WRITE (stdout, nml=croco_diags_eddy)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diags_eddy_avg)
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      MPI_master_only WRITE (stdout, nml=croco_diagnostics_bio)
#  ifdef AVERAGES
      MPI_master_only WRITE (stdout, nml=croco_diagbio_avg)
#  endif
#endif
#ifdef BODYFORCE
      MPI_master_only WRITE (stdout, nml=croco_bodyforce)
#endif
#if !defined NONLIN_EOS
      MPI_master_only WRITE (stdout, nml=croco_lin_eos)
#endif
#ifdef SEDIMENT
      MPI_master_only WRITE (stdout, nml=croco_sediments)
#endif
#ifdef MUSTANG
      MPI_master_only WRITE (stdout, nml=croco_sediments_mustang)
#endif
#ifdef SUBSTANCE
      MPI_master_only WRITE (stdout, nml=croco_substance)
#endif
#ifdef OBSTRUCTION
      MPI_master_only WRITE (stdout, nml=croco_obstruction)
#endif
#ifdef XIOS
      MPI_master_only WRITE (stdout, nml=croco_xios_origin_date)
#endif
#ifdef ASSIMILATION
      MPI_master_only WRITE (stdout, nml=croco_assimilation)
#endif
      MPI_master_only WRITE (stdout, nml=croco_rho0)
      MPI_master_only WRITE (stdout, nml=croco_bottom_drag)
      MPI_master_only WRITE (stdout, nml=croco_gamma2)
      MPI_master_only WRITE (stdout, nml=croco_lateral_visc)
#ifdef SOLVE3D
#  ifdef TRACERS
      MPI_master_only WRITE (stdout, nml=croco_tracer_diff2)
      MPI_master_only WRITE (stdout, nml=croco_tracer_diff4)
#  endif
#  if !defined LMD_MIXING && !defined GLS_MIXING
      MPI_master_only WRITE (stdout, nml=croco_vertical_mixing)
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
