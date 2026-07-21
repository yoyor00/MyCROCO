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
!=======================================================================

MODULE croco_namelist_init
   implicit none
   private

   public :: init_all

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
#ifdef SOLVE3D
      call init_croco_s_coord()
#endif
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
#if defined ANA_PSOURCE
      call init_psource(ierr)
#endif
#if defined T_FRC_BRY     || defined M2_FRC_BRY    || \
      defined M3_FRC_BRY||defined Z_FRC_BRY||\
      defined W_FRC_BRY||defined NBQ_FRC_BRY||\
      defined TCLIMATOLOGY||defined M2CLIMATOLOGY||\
      defined M3CLIMATOLOGY||defined ZCLIMATOLOGY||\
      defined WCLIMATOLOGY||defined NBQCLIMATOLOGY
      call init_nudging()
#endif
#if defined BHFLUX || (defined BWFLUX && defined SALINITY)
      call init_bottom_forcing(ierr)
#endif
      if (use_frcname) then
         call init_forcing(ierr)
      end if
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA && !defined ONLINE
      call init_bulk_forcing(ierr)
#endif
#if (  defined TCLIMATOLOGY  && \
      !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY&&\
      !defined ANA_SSH) || \
      (defined M2CLIMATOLOGY&&\
      !defined ANA_M2CLIMA) || \
      (defined M3CLIMATOLOGY&&\
      !defined ANA_M3CLIMA)
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
#if defined OUTPUTS_SURFACE && !defined XIOS
      call init_surf(ierr)
#  ifdef AVERAGES
      call init_surf_avg(ierr)
#  endif
#endif
#if defined DIAGNOSTICS_TS
      call init_diagnostics_ts(ierr)
#  ifdef AVERAGES
      call init_diag_avg(ierr)
#  endif
#endif
#if defined DIAGNOSTICS_UV
      call init_diagnosticsM(ierr)
#  ifdef AVERAGES
      call init_diagM_avg(ierr)
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      call init_diags_vrt(ierr)
#  ifdef AVERAGES
      call init_diags_vrt_avg(ierr)
#  endif
#endif
#ifdef DIAGNOSTICS_KE
      call init_diags_ek(ierr)
#  ifdef AVERAGES
      call init_diags_ek_avg(ierr)
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      call init_diags_pv(ierr)
#  ifdef AVERAGES
      call init_diags_pv_avg(ierr)
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
#  ifdef AVERAGES
      call init_diags_eddy_avg(ierr)
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      call init_diagnostics_bio(ierr)
#  ifdef AVERAGES
      call init_diagbio_avg(ierr)
#  endif
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
#ifdef SEDIMENT
      call init_sediments(ierr)
#endif
#ifdef MUSTANG
      call init_sediments_mustang(ierr)
#endif
#ifdef SUBSTANCE
      call init_substance(ierr)
#endif
#ifdef OBSTRUCTION
      call init_obstruction(ierr)
#endif
#ifdef XIOS
      call init_xios_origin_date()
#endif
#ifdef ASSIMILATION
      call init_assimilation(ierr)
#endif
#ifdef LOGFILE
      if (ierr == 0) call init_logfile(ierr)
#endif
#ifdef STATIONS
      call init_stations(ierr)
#endif
      call init_primary_history_fields()
#ifdef AVERAGES
      call init_primary_averages_fields()
#endif
#if defined SOLVE3D || (!defined SOLVE3D && defined RIP)
      call init_auxiliary_history_fields()
#endif
#if defined SOLVE3D && defined TEMPERATURE
      call init_temperature_history_fields()
#endif
#if defined SOLVE3D && defined SALINITY
      call init_salinity_history_fields()
#endif
#if defined SOLVE3D && defined BULK_FLUX
      call init_bulk_flux_history_fields()
#endif
#if defined SOLVE3D && defined BHFLUX
      call init_bhflux_history_fields()
#endif
#if defined SOLVE3D && defined BWFLUX && defined SALINITY
      call init_bwflux_history_fields()
#endif
#if defined SOLVE3D && \
      (defined ANA_VMIX||defined LMD_MIXING||\
      defined LMD_SKPP||defined LMD_BKPP||\
      defined GLS_MIXING)
      call init_bvf_history_fields()
#endif
#if defined SOLVE3D && (defined LMD_SKPP || defined GLS_MIXING)
      call init_hbl_history_fields()
#endif
#if defined SOLVE3D && defined LMD_BKPP
      call init_lmd_bkpp_history_fields()
#endif
#if defined SOLVE3D && defined VIS_COEF_3D
      call init_vis_coef_history_fields()
#endif
#if defined SOLVE3D && defined DIF_COEF_3D
      call init_dif_coef_history_fields()
#endif
#if defined SOLVE3D && defined BIOLOGY && !defined PISCES
      call init_biology_history_fields()
#endif
#if defined SOLVE3D && defined BIO_NChlPZD
      call init_bio_nchlpzd_history_fields()
#endif
#if defined SOLVE3D && defined BIO_BioEBUS
      call init_bio_bioebus_history_fields()
#endif
#if defined SOLVE3D && defined MORPHODYN
      call init_morphodyn_history_fields()
#endif
#if defined AVERAGES && defined SOLVE3D
      call init_auxiliary_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined TEMPERATURE
      call init_temperature_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined SALINITY
      call init_salinity_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined BULK_FLUX
      call init_bulk_flux_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined BHFLUX
      call init_bhflux_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined BWFLUX && defined SALINITY
      call init_bwflux_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && \
      (defined ANA_VMIX||defined LMD_MIXING||\
      defined LMD_SKPP||defined LMD_BKPP||\
      defined GLS_MIXING)
      call init_bvf_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && \
      (defined LMD_SKPP||defined GLS_MIXING)
      call init_hbl_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined LMD_BKPP
      call init_lmd_bkpp_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined VIS_COEF_3D
      call init_vis_coef_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined DIF_COEF_3D
      call init_dif_coef_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined BIOLOGY && !defined PISCES
      call init_biology_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined BIO_NChlPZD
      call init_bio_nchlpzd_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined BIO_BioEBUS
      call init_bio_bioebus_averages_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined MORPHODYN
      call init_morphodyn_averages_fields()
#endif
#ifdef DIAGNOSTICS_TS
      call init_diag_ts_fields()
#endif
#ifdef DIAGNOSTICS_UV
      call init_diagM_fields()
#endif
#ifdef DIAGNOSTICS_VRT
      call init_diags_vrt_fields()
#endif
#ifdef DIAGNOSTICS_KE
      call init_diags_ek_fields()
#endif
#ifdef DIAGNOSTICS_PV
      call init_diags_pv_fields()
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS && defined AVERAGES
      call init_diags_eddy_fields()
#endif
#ifdef DIAGNOSTICS_BIO
      call init_diagbio_fields()
#endif
#ifdef STOGEN
      call init_stochastic_history_fields()
#endif
#if defined SOLVE3D && defined GLS_MIXING
      call init_gls_history_fields()
#endif
#if defined AVERAGES && defined SOLVE3D && defined GLS_MIXING
      call init_gls_averages_fields()
#endif
#if defined ABL1D && !defined XIOS
      call init_abl_history_fields()
# ifdef AVERAGES
      call init_abl_averages_fields()
# endif
#endif
#if defined OUTPUTS_SURFACE && !defined XIOS
      call init_surf_history_fields()
# ifdef AVERAGES
      call init_surf_average_fields()
# endif
#endif
#ifdef STATIONS
      call init_station_fields()
#endif
#if defined SOLVE3D && defined SEDIMENT
      call init_sediment_history_fields()
#endif
#ifdef BBL
      call init_bbl_history_fields()
#endif
#ifdef MRL_WCI
      call init_wci_history_fields()
# ifdef AVERAGES
      call init_wci_average_fields()
# endif
#endif
#if defined MRL_WCI || defined OW_COUPLING
      call init_wave_history_fields()
# ifdef AVERAGES
      call init_wave_average_fields()
# endif
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

#ifdef SOLVE3D
   !---------------------------------------------------------------------
   !  init_croco_s_coord
   !  Set Tcline = hc for backward-compatible NetCDF global attribute.
   !---------------------------------------------------------------------
   subroutine init_croco_s_coord()
      use croco_namelist, ONLY: hc, Tcline
      implicit none

      Tcline = hc

   end subroutine init_croco_s_coord
#endif

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

#if defined T_FRC_BRY     || defined M2_FRC_BRY    || \
   defined M3_FRC_BRY||defined Z_FRC_BRY||\
   defined TCLIMATOLOGY||defined M2CLIMATOLOGY||\
   defined M3CLIMATOLOGY||defined ZCLIMATOLOGY
   !---------------------------------------------------------------------
   !  init_nudging
   !  Convert nudging time-scales from days (namelist input) to sec^-1
   !  (expected by the model). Applies only on root when AGRIF is active
   !  without Orlanski OBC schemes, matching the original read_inp.F logic.
   !---------------------------------------------------------------------
   subroutine init_nudging()
      use croco_namelist, ONLY: tauT_in, tauT_out, tauM_in, tauM_out
#  if defined AGRIF && !defined AGRIF_OBC_M2ORLANSKI && \
      !defined AGRIF_OBC_M3ORLANSKI && !defined AGRIF_OBC_TORLANSKI
      use Agrif_Util
#  endif
      implicit none

#  if defined AGRIF && !defined AGRIF_OBC_M2ORLANSKI && \
      !defined AGRIF_OBC_M3ORLANSKI && !defined AGRIF_OBC_TORLANSKI
      if (.not. Agrif_Root()) return
#  endif
      tauT_in = 1.0/(tauT_in*86400.0)
      tauT_out = 1.0/(tauT_out*86400.0)
      tauM_in = 1.0/(tauM_in*86400.0)
      tauM_out = 1.0/(tauM_out*86400.0)

   end subroutine init_nudging
#endif

#if defined BHFLUX || (defined BWFLUX && defined SALINITY)
   !---------------------------------------------------------------------
   !  init_bottom_forcing
   !  Adjust btfname for MPI and check file availability.
   !---------------------------------------------------------------------
   subroutine init_bottom_forcing(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: btfname
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      call adjust_filename_parallel(btfname, "btfname", ierr)
      if (ierr /= 0) return

      open (testunit, file=trim(btfname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open bottom forcing file ', trim(btfname)
         ierr = ierr + 1
      end if

   end subroutine init_bottom_forcing
#endif

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

#if (  defined TCLIMATOLOGY  && \
   !defined ANA_TCLIMA) || \
   (defined ZCLIMATOLOGY&&\
   !defined ANA_SSH) || \
   (defined M2CLIMATOLOGY&&\
   !defined ANA_M2CLIMA) || \
   (defined M3CLIMATOLOGY&&\
   !defined ANA_M3CLIMA)
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

#ifdef SEDIMENT
   !---------------------------------------------------------------------
   !  init_sediments
   !---------------------------------------------------------------------
   subroutine init_sediments(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: sedname
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      open (testunit, file=trim(sedname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open sediment file ', trim(sedname)
         ierr = ierr + 1
      end if

   end subroutine init_sediments
#endif

#ifdef MUSTANG
   !---------------------------------------------------------------------
   !  init_sediments_mustang
   !---------------------------------------------------------------------
   subroutine init_sediments_mustang(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: sedname_must
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      open (testunit, file=trim(sedname_must), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open MUSTANG sediment file ', trim(sedname_must)
         ierr = ierr + 1
      end if

   end subroutine init_sediments_mustang
#endif

#ifdef SUBSTANCE
   !---------------------------------------------------------------------
   !  init_substance
   !---------------------------------------------------------------------
   subroutine init_substance(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: subsfilename
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      open (testunit, file=trim(subsfilename), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open substance file ', trim(subsfilename)
         ierr = ierr + 1
      end if

   end subroutine init_substance
#endif

#ifdef OBSTRUCTION
   !---------------------------------------------------------------------
   !  init_obstruction
   !---------------------------------------------------------------------
   subroutine init_obstruction(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: obstname
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      open (testunit, file=trim(obstname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open obstruction file ', trim(obstname)
         ierr = ierr + 1
      end if

   end subroutine init_obstruction
#endif

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

#if defined OUTPUTS_SURFACE && !defined XIOS
   !---------------------------------------------------------------------
   !  init_surf
   !  Adjust surfname for MPI/ENSEMBLE and apply nwrtsurf default.
   !---------------------------------------------------------------------
   subroutine init_surf(ierr)
      use croco_namelist, ONLY: surfname, nwrtsurf, nwrt
      implicit none
      integer, intent(inout) :: ierr

      call adjust_filename_parallel(surfname, "surfname", ierr)
      call adjust_filename_ensemble(surfname)
      if (ierr /= 0) return
      if (nwrtsurf == 0) nwrtsurf = nwrt

   end subroutine init_surf
#  ifdef AVERAGES
   !---------------------------------------------------------------------
   !  init_surf_avg
   !  Adjust surfname_avg for MPI/ENSEMBLE and apply nwrtsurf_avg default.
   !---------------------------------------------------------------------
   subroutine init_surf_avg(ierr)
      use croco_namelist, ONLY: surfname_avg, nwrtsurf_avg, navg
      implicit none
      integer, intent(inout) :: ierr

      call adjust_filename_parallel(surfname_avg, "surfname_avg", ierr)
      call adjust_filename_ensemble(surfname_avg)
      if (ierr /= 0) return
      if (nwrtsurf_avg == 0) nwrtsurf_avg = navg

   end subroutine init_surf_avg
#  endif
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
   !  NOTE : the COMMON block is legacy; it should eventually be replaced
   !          by a proper module variable.
   !---------------------------------------------------------------------
   subroutine init_time_stepping_nbq()
      use scalars, ONLY: dtfast
      implicit none

      real :: dtnbq
      common/time_nbq2/dtnbq   ! TODO: replace with module variable

      dtnbq = dtfast

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

#ifdef XIOS
   !---------------------------------------------------------------------
   !  init_xios_origin_date
   !  Compute xios_origin_date_in_sec from the date string.
   !---------------------------------------------------------------------
   subroutine init_xios_origin_date()
      use croco_namelist, ONLY: xios_origin_date
      use ncscrum, ONLY: xios_origin_date_in_sec
      implicit none

      real(kind=8), external :: tool_datosec

      xios_origin_date_in_sec = tool_datosec(xios_origin_date)

   end subroutine init_xios_origin_date
#endif

#ifdef ASSIMILATION
   !---------------------------------------------------------------------
   !  init_assimilation
   !  Check that both assimilation files exist.
   !---------------------------------------------------------------------
   subroutine init_assimilation(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: aparnam, assname
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: ios

      open (testunit, file=trim(aparnam), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open assimilation parameters file ', trim(aparnam)
         ierr = ierr + 1
      end if

      open (testunit, file=trim(assname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open assimilation data file ', trim(assname)
         ierr = ierr + 1
      end if

   end subroutine init_assimilation
#endif

#if defined DIAGNOSTICS_TS
   subroutine init_diagnostics_ts(ierr)
      use croco_namelist, ONLY: dianame, nwrtdia, nwrt
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(dianame, "dianame", ierr)
      call adjust_filename_ensemble(dianame)
      if (ierr /= 0) return
      if (nwrtdia == 0) nwrtdia = nwrt
   end subroutine init_diagnostics_ts
#  ifdef AVERAGES
   subroutine init_diag_avg(ierr)
      use croco_namelist, ONLY: dianame_avg, nwrtdia_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(dianame_avg, "dianame_avg", ierr)
      call adjust_filename_ensemble(dianame_avg)
      if (ierr /= 0) return
      if (nwrtdia_avg == 0) nwrtdia_avg = navg
   end subroutine init_diag_avg
#  endif
#endif

#if defined DIAGNOSTICS_UV
   subroutine init_diagnosticsM(ierr)
      use croco_namelist, ONLY: dianameM, nwrtdiaM, nwrt
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(dianameM, "dianameM", ierr)
      call adjust_filename_ensemble(dianameM)
      if (ierr /= 0) return
      if (nwrtdiaM == 0) nwrtdiaM = nwrt
   end subroutine init_diagnosticsM
#  ifdef AVERAGES
   subroutine init_diagM_avg(ierr)
      use croco_namelist, ONLY: dianameM_avg, nwrtdiaM_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(dianameM_avg, "dianameM_avg", ierr)
      call adjust_filename_ensemble(dianameM_avg)
      if (ierr /= 0) return
      if (nwrtdiaM_avg == 0) nwrtdiaM_avg = navg
   end subroutine init_diagM_avg
#  endif
#endif

#ifdef DIAGNOSTICS_VRT
   subroutine init_diags_vrt(ierr)
      use croco_namelist, ONLY: diags_vrtname, nwrtdiags_vrt, nwrt
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_vrtname, "diags_vrtname", ierr)
      call adjust_filename_ensemble(diags_vrtname)
      if (ierr /= 0) return
      if (nwrtdiags_vrt == 0) nwrtdiags_vrt = nwrt
   end subroutine init_diags_vrt
#  ifdef AVERAGES
   subroutine init_diags_vrt_avg(ierr)
      use croco_namelist, ONLY: diags_vrtname_avg, nwrtdiags_vrt_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_vrtname_avg, "diags_vrtname_avg", ierr)
      call adjust_filename_ensemble(diags_vrtname_avg)
      if (ierr /= 0) return
      if (nwrtdiags_vrt_avg == 0) nwrtdiags_vrt_avg = navg
   end subroutine init_diags_vrt_avg
#  endif
#endif

#ifdef DIAGNOSTICS_KE
   subroutine init_diags_ek(ierr)
      use croco_namelist, ONLY: diags_ekname, nwrtdiags_ek, nwrt
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_ekname, "diags_ekname", ierr)
      call adjust_filename_ensemble(diags_ekname)
      if (ierr /= 0) return
      if (nwrtdiags_ek == 0) nwrtdiags_ek = nwrt
   end subroutine init_diags_ek
#  ifdef AVERAGES
   subroutine init_diags_ek_avg(ierr)
      use croco_namelist, ONLY: diags_ekname_avg, nwrtdiags_ek_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_ekname_avg, "diags_ekname_avg", ierr)
      call adjust_filename_ensemble(diags_ekname_avg)
      if (ierr /= 0) return
      if (nwrtdiags_ek_avg == 0) nwrtdiags_ek_avg = navg
   end subroutine init_diags_ek_avg
#  endif
#endif

#ifdef DIAGNOSTICS_PV
   subroutine init_diags_pv(ierr)
      use croco_namelist, ONLY: diags_pvname, nwrtdiags_pv, nwrt
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_pvname, "diags_pvname", ierr)
      call adjust_filename_ensemble(diags_pvname)
      if (ierr /= 0) return
      if (nwrtdiags_pv == 0) nwrtdiags_pv = nwrt
   end subroutine init_diags_pv
#  ifdef AVERAGES
   subroutine init_diags_pv_avg(ierr)
      use croco_namelist, ONLY: diags_pvname_avg, nwrtdiags_pv_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_pvname_avg, "diags_pvname_avg", ierr)
      call adjust_filename_ensemble(diags_pvname_avg)
      if (ierr /= 0) return
      if (nwrtdiags_pv_avg == 0) nwrtdiags_pv_avg = navg
   end subroutine init_diags_pv_avg
#  endif
#endif

#if defined DIAGNOSTICS_EDDY && !defined XIOS
#  ifdef AVERAGES
   subroutine init_diags_eddy_avg(ierr)
      use croco_namelist, ONLY: diags_eddyname_avg, nwrtdiags_eddy_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(diags_eddyname_avg, "diags_eddyname_avg", ierr)
      call adjust_filename_ensemble(diags_eddyname_avg)
      if (ierr /= 0) return
      if (nwrtdiags_eddy_avg == 0) nwrtdiags_eddy_avg = navg
   end subroutine init_diags_eddy_avg
#  endif
#endif

#ifdef DIAGNOSTICS_BIO
   subroutine init_diagnostics_bio(ierr)
      use croco_namelist, ONLY: dianamebio, nwrtdiabio, nwrt
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(dianamebio, "dianamebio", ierr)
      call adjust_filename_ensemble(dianamebio)
      if (ierr /= 0) return
      if (nwrtdiabio == 0) nwrtdiabio = nwrt
   end subroutine init_diagnostics_bio
#  ifdef AVERAGES
   subroutine init_diagbio_avg(ierr)
      use croco_namelist, ONLY: dianamebio_avg, nwrtdiabio_avg, navg
      implicit none
      integer, intent(inout) :: ierr
      call adjust_filename_parallel(dianamebio_avg, "dianamebio_avg", ierr)
      call adjust_filename_ensemble(dianamebio_avg)
      if (ierr /= 0) return
      if (nwrtdiabio_avg == 0) nwrtdiabio_avg = navg
   end subroutine init_diagbio_avg
#  endif
#endif
#ifdef STATIONS
   subroutine init_stations(ierr)
      use croco_namelist, ONLY: staname, staposname
#  ifdef AGRIF
      use croco_namelist, ONLY: ldefsta, nsta, nrpfsta
      use Agrif_Util
#  endif
      implicit none
      integer, intent(inout) :: ierr
#  ifdef AGRIF
      if (.not. Agrif_Root()) return
#  endif
      call adjust_filename_parallel(staname, "staname", ierr)
      call adjust_filename_ensemble(staname)
   end subroutine init_stations
#endif

   !---------------------------------------------------------------------
   !  init_primary_history_fields
   !  Copy out_his_* module variables into the legacy wrthis() array.
   !---------------------------------------------------------------------
   subroutine init_primary_history_fields()
      use croco_namelist, ONLY: out_his_zeta, out_his_ubar, out_his_vbar
#ifdef SOLVE3D
      use croco_namelist, ONLY: out_his_u, out_his_v
# ifdef TRACERS
      use croco_namelist, ONLY: out_his_tracer
# endif
#endif
      use param, ONLY: NT
      use ncscrum, ONLY: wrthis, indxTime, indxZ, indxUb, indxVb
#ifdef SOLVE3D
      use ncscrum, ONLY: indxU, indxV
#endif
      implicit none
      integer :: itrc

      wrthis(indxZ) = out_his_zeta
      wrthis(indxUb) = out_his_ubar
      wrthis(indxVb) = out_his_vbar
      if (out_his_zeta .or. out_his_ubar .or. out_his_vbar) wrthis(indxTime) = .true.
#ifdef SOLVE3D
      wrthis(indxU) = out_his_u
      wrthis(indxV) = out_his_v
      if (out_his_u .or. out_his_v) wrthis(indxTime) = .true.
# ifdef TRACERS
      do itrc = 1, NT
         wrthis(indxV + itrc) = out_his_tracer(itrc)
         if (out_his_tracer(itrc)) wrthis(indxTime) = .true.
      end do
# endif
#endif
   end subroutine init_primary_history_fields

#ifdef AVERAGES
   !---------------------------------------------------------------------
   !  init_primary_averages_fields
   !  Copy out_avg_* module variables into the legacy wrtavg() array.
   !---------------------------------------------------------------------
   subroutine init_primary_averages_fields()
      use croco_namelist, ONLY: out_avg_zeta, out_avg_ubar, out_avg_vbar
# ifdef SOLVE3D
      use croco_namelist, ONLY: out_avg_u, out_avg_v
#  ifdef TRACERS
      use croco_namelist, ONLY: out_avg_tracer
#  endif
# endif
      use param, ONLY: NT
      use ncscrum, ONLY: wrtavg, indxTime, indxZ, indxUb, indxVb
#ifdef SOLVE3D
      use ncscrum, ONLY: indxU, indxV
#endif
      implicit none
      integer :: itrc

      wrtavg(indxZ) = out_avg_zeta
      wrtavg(indxUb) = out_avg_ubar
      wrtavg(indxVb) = out_avg_vbar
      if (out_avg_zeta .or. out_avg_ubar .or. out_avg_vbar) wrtavg(indxTime) = .true.
# ifdef SOLVE3D
      wrtavg(indxU) = out_avg_u
      wrtavg(indxV) = out_avg_v
      if (out_avg_u .or. out_avg_v) wrtavg(indxTime) = .true.
#  ifdef TRACERS
      do itrc = 1, NT
         wrtavg(indxV + itrc) = out_avg_tracer(itrc)
         if (out_avg_tracer(itrc)) wrtavg(indxTime) = .true.
      end do
#  endif
# endif
   end subroutine init_primary_averages_fields
#endif

#if defined SOLVE3D || (!defined SOLVE3D && defined RIP)
   !---------------------------------------------------------------------
   !  init_auxiliary_history_fields
   !  SOLVE3D basics (rho, omega, w, akv) + surface stress (always).
   !---------------------------------------------------------------------
   subroutine init_auxiliary_history_fields()
      use croco_namelist, ONLY: &
         out_his_rho, out_his_omega, out_his_w, out_his_akv, &
         out_his_bostr, out_his_bustr, out_his_bvstr, out_his_wstr, out_his_ustr, out_his_vstr
      use ncscrum, ONLY: wrthis, indxTime, &
                         indxBostr, indxBustr, indxBvstr, &
                         indxWstr, indxUWstr, indxVWstr
#ifdef SOLVE3D
      use ncscrum, ONLY: indxR, indxO, indxW, indxAkv
#endif
      implicit none

      wrthis(indxBostr) = out_his_bostr
      wrthis(indxBustr) = out_his_bustr
      wrthis(indxBvstr) = out_his_bvstr
      wrthis(indxWstr) = out_his_wstr
      wrthis(indxUWstr) = out_his_ustr
      wrthis(indxVWstr) = out_his_vstr
      if (out_his_bostr .or. out_his_bustr .or. out_his_bvstr .or. &
          out_his_wstr .or. out_his_ustr .or. out_his_vstr) wrthis(indxTime) = .true.
#ifdef SOLVE3D
      wrthis(indxR) = out_his_rho
      wrthis(indxO) = out_his_omega
      wrthis(indxW) = out_his_w
      wrthis(indxAkv) = out_his_akv
      if (out_his_rho .or. out_his_omega .or. out_his_w .or. out_his_akv) wrthis(indxTime) = .true.
#endif
   end subroutine init_auxiliary_history_fields
#endif

#if defined SOLVE3D && defined TEMPERATURE
   subroutine init_temperature_history_fields()
      use croco_namelist, ONLY: out_his_akt, out_his_shflx, out_his_shflx_rsw
      use ncscrum, ONLY: wrthis, indxTime, indxAkt, indxShflx, indxShflx_rsw
      implicit none
      wrthis(indxAkt) = out_his_akt
      wrthis(indxShflx) = out_his_shflx
      wrthis(indxShflx_rsw) = out_his_shflx_rsw
      if (out_his_akt .or. out_his_shflx .or. out_his_shflx_rsw) wrthis(indxTime) = .true.
   end subroutine init_temperature_history_fields
#endif

#if defined SOLVE3D && defined SALINITY
   subroutine init_salinity_history_fields()
      use croco_namelist, ONLY: out_his_aks, out_his_swflx
      use ncscrum, ONLY: wrthis, indxTime, indxAks, indxSwflx
      implicit none
      wrthis(indxAks) = out_his_aks
      wrthis(indxSwflx) = out_his_swflx
      if (out_his_aks .or. out_his_swflx) wrthis(indxTime) = .true.
   end subroutine init_salinity_history_fields
#endif

#if defined SOLVE3D && defined BULK_FLUX
   subroutine init_bulk_flux_history_fields()
      use croco_namelist, ONLY: out_his_shflx_rlw, out_his_shflx_lat, out_his_shflx_sen
      use ncscrum, ONLY: wrthis, indxTime, indxShflx_rlw, indxShflx_lat, indxShflx_sen
      implicit none
      wrthis(indxShflx_rlw) = out_his_shflx_rlw
      wrthis(indxShflx_lat) = out_his_shflx_lat
      wrthis(indxShflx_sen) = out_his_shflx_sen
      if (out_his_shflx_rlw .or. out_his_shflx_lat .or. out_his_shflx_sen) wrthis(indxTime) = .true.
   end subroutine init_bulk_flux_history_fields
#endif

#if defined SOLVE3D && defined BHFLUX
   subroutine init_bhflux_history_fields()
      use croco_namelist, ONLY: out_his_bhflx
      use ncscrum, ONLY: wrthis, indxTime, indxBhflx
      implicit none
      wrthis(indxBhflx) = out_his_bhflx
      if (out_his_bhflx) wrthis(indxTime) = .true.
   end subroutine init_bhflux_history_fields
#endif

#if defined SOLVE3D && defined BWFLUX && defined SALINITY
   subroutine init_bwflux_history_fields()
      use croco_namelist, ONLY: out_his_bwflx
      use ncscrum, ONLY: wrthis, indxTime, indxBwflx
      implicit none
      wrthis(indxBwflx) = out_his_bwflx
      if (out_his_bwflx) wrthis(indxTime) = .true.
   end subroutine init_bwflux_history_fields
#endif

#if defined SOLVE3D && (defined ANA_VMIX || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING)
   subroutine init_bvf_history_fields()
      use croco_namelist, ONLY: out_his_bvf
      use ncscrum, ONLY: wrthis, indxTime, indxbvf
      implicit none
      wrthis(indxbvf) = out_his_bvf
      if (out_his_bvf) wrthis(indxTime) = .true.
   end subroutine init_bvf_history_fields
#endif

#if defined SOLVE3D && (defined LMD_SKPP || defined GLS_MIXING)
   subroutine init_hbl_history_fields()
      use croco_namelist, ONLY: out_his_hbl
      use ncscrum, ONLY: wrthis, indxTime, indxHbl
      implicit none
      wrthis(indxHbl) = out_his_hbl
      if (out_his_hbl) wrthis(indxTime) = .true.
   end subroutine init_hbl_history_fields
#endif

#if defined SOLVE3D && defined LMD_BKPP
   subroutine init_lmd_bkpp_history_fields()
      use croco_namelist, ONLY: out_his_hbbl
      use ncscrum, ONLY: wrthis, indxTime, indxHbbl
      implicit none
      wrthis(indxHbbl) = out_his_hbbl
      if (out_his_hbbl) wrthis(indxTime) = .true.
   end subroutine init_lmd_bkpp_history_fields
#endif

#if defined SOLVE3D && defined VIS_COEF_3D
   subroutine init_vis_coef_history_fields()
      use croco_namelist, ONLY: out_his_visc3d
      use ncscrum, ONLY: wrthis, indxTime, indxVisc
      implicit none
      wrthis(indxVisc) = out_his_visc3d
      if (out_his_visc3d) wrthis(indxTime) = .true.
   end subroutine init_vis_coef_history_fields
#endif

#if defined SOLVE3D && defined DIF_COEF_3D
   subroutine init_dif_coef_history_fields()
      use croco_namelist, ONLY: out_his_diff3d
      use ncscrum, ONLY: wrthis, indxTime, indxDiff
      implicit none
      wrthis(indxDiff) = out_his_diff3d
      if (out_his_diff3d) wrthis(indxTime) = .true.
   end subroutine init_dif_coef_history_fields
#endif

#if defined SOLVE3D && defined BIOLOGY && !defined PISCES
   subroutine init_biology_history_fields()
      use croco_namelist, ONLY: out_his_hel
      use ncscrum, ONLY: wrthis, indxTime, indxHel
      implicit none
      wrthis(indxHel) = out_his_hel
      if (out_his_hel) wrthis(indxTime) = .true.
   end subroutine init_biology_history_fields
#endif

#if defined SOLVE3D && defined BIO_NChlPZD
   subroutine init_bio_nchlpzd_history_fields()
      use croco_namelist, ONLY: out_his_chc, out_his_u10, out_his_kvo2, out_his_o2sat
      use ncscrum, ONLY: wrthis, indxTime, indxChC
# ifdef OXYGEN
      use ncscrum, ONLY: indxU10, indxKvO2, indxO2sat
# endif
      implicit none
      wrthis(indxChC) = out_his_chc
      if (out_his_chc) wrthis(indxTime) = .true.
# ifdef OXYGEN
      wrthis(indxU10) = out_his_u10
      wrthis(indxKvO2) = out_his_kvo2
      wrthis(indxO2sat) = out_his_o2sat
      if (out_his_u10 .or. out_his_kvo2 .or. out_his_o2sat) wrthis(indxTime) = .true.
# endif
   end subroutine init_bio_nchlpzd_history_fields
#endif

#if defined SOLVE3D && defined BIO_BioEBUS
   subroutine init_bio_bioebus_history_fields()
      use croco_namelist, ONLY: out_his_aou, out_his_wind10
      use ncscrum, ONLY: wrthis, indxTime, indxAOU, indxWIND10
      implicit none
      wrthis(indxAOU) = out_his_aou
      wrthis(indxWIND10) = out_his_wind10
      if (out_his_aou .or. out_his_wind10) wrthis(indxTime) = .true.
   end subroutine init_bio_bioebus_history_fields
#endif

# if defined SOLVE3D && defined MORPHODYN
   !---------------------------------------------------------------------
   !  init_morphodyn_history_fields
   !  Morphodynamics output field (MORPHODYN).
   !---------------------------------------------------------------------
   subroutine init_morphodyn_history_fields()
      use croco_namelist, ONLY: out_his_hm
      use ncscrum, ONLY: wrthis, indxTime, indxHm
      implicit none

      wrthis(indxHm) = out_his_hm
      if (out_his_hm) wrthis(indxTime) = .true.
   end subroutine init_morphodyn_history_fields
# endif

#if defined AVERAGES && defined SOLVE3D
   !---------------------------------------------------------------------
   !  init_auxiliary_averages_fields
   !  SOLVE3D basics (rho, omega, w, akv) + surface stress.
   !---------------------------------------------------------------------
   subroutine init_auxiliary_averages_fields()
      use croco_namelist, ONLY: &
         out_avg_rho, out_avg_omega, out_avg_w, out_avg_akv, &
         out_avg_bostr, out_avg_bustr, out_avg_bvstr, out_avg_wstr, out_avg_ustr, out_avg_vstr
      use ncscrum, ONLY: wrtavg, indxTime, &
                         indxR, indxO, indxW, indxAkv, &
                         indxBostr, indxBustr, indxBvstr, &
                         indxWstr, indxUWstr, indxVWstr
      implicit none

      wrtavg(indxR) = out_avg_rho
      wrtavg(indxO) = out_avg_omega
      wrtavg(indxW) = out_avg_w
      wrtavg(indxAkv) = out_avg_akv
      if (out_avg_rho .or. out_avg_omega .or. out_avg_w .or. out_avg_akv) wrtavg(indxTime) = .true.
      wrtavg(indxBostr) = out_avg_bostr
      wrtavg(indxBustr) = out_avg_bustr
      wrtavg(indxBvstr) = out_avg_bvstr
      wrtavg(indxWstr) = out_avg_wstr
      wrtavg(indxUWstr) = out_avg_ustr
      wrtavg(indxVWstr) = out_avg_vstr
      if (out_avg_bostr .or. out_avg_bustr .or. out_avg_bvstr .or. &
          out_avg_wstr .or. out_avg_ustr .or. out_avg_vstr) wrtavg(indxTime) = .true.
   end subroutine init_auxiliary_averages_fields
#endif

# if defined AVERAGES && defined SOLVE3D && defined TEMPERATURE
   subroutine init_temperature_averages_fields()
      use croco_namelist, ONLY: out_avg_akt, out_avg_shflx, out_avg_shflx_rsw
      use ncscrum, ONLY: wrtavg, indxTime, indxAkt, indxShflx, indxShflx_rsw
      implicit none
      wrtavg(indxAkt) = out_avg_akt
      wrtavg(indxShflx) = out_avg_shflx
      wrtavg(indxShflx_rsw) = out_avg_shflx_rsw
      if (out_avg_akt .or. out_avg_shflx .or. out_avg_shflx_rsw) wrtavg(indxTime) = .true.
   end subroutine init_temperature_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined SALINITY
   subroutine init_salinity_averages_fields()
      use croco_namelist, ONLY: out_avg_aks, out_avg_swflx
      use ncscrum, ONLY: wrtavg, indxTime, indxAks, indxSwflx
      implicit none
      wrtavg(indxAks) = out_avg_aks
      wrtavg(indxSwflx) = out_avg_swflx
      if (out_avg_aks .or. out_avg_swflx) wrtavg(indxTime) = .true.
   end subroutine init_salinity_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined BULK_FLUX
   subroutine init_bulk_flux_averages_fields()
      use croco_namelist, ONLY: out_avg_shflx_rlw, out_avg_shflx_lat, out_avg_shflx_sen
      use ncscrum, ONLY: wrtavg, indxTime, indxShflx_rlw, indxShflx_lat, indxShflx_sen
      implicit none
      wrtavg(indxShflx_rlw) = out_avg_shflx_rlw
      wrtavg(indxShflx_lat) = out_avg_shflx_lat
      wrtavg(indxShflx_sen) = out_avg_shflx_sen
      if (out_avg_shflx_rlw .or. out_avg_shflx_lat .or. out_avg_shflx_sen) wrtavg(indxTime) = .true.
   end subroutine init_bulk_flux_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined BHFLUX
   subroutine init_bhflux_averages_fields()
      use croco_namelist, ONLY: out_avg_bhflx
      use ncscrum, ONLY: wrtavg, indxTime, indxBhflx
      implicit none
      wrtavg(indxBhflx) = out_avg_bhflx
      if (out_avg_bhflx) wrtavg(indxTime) = .true.
   end subroutine init_bhflux_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined BWFLUX && defined SALINITY
   subroutine init_bwflux_averages_fields()
      use croco_namelist, ONLY: out_avg_bwflx
      use ncscrum, ONLY: wrtavg, indxTime, indxBwflx
      implicit none
      wrtavg(indxBwflx) = out_avg_bwflx
      if (out_avg_bwflx) wrtavg(indxTime) = .true.
   end subroutine init_bwflux_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && \
   (defined ANA_VMIX||defined LMD_MIXING||\
   defined LMD_SKPP||defined LMD_BKPP||\
   defined GLS_MIXING)
   subroutine init_bvf_averages_fields()
      use croco_namelist, ONLY: out_avg_bvf
      use ncscrum, ONLY: wrtavg, indxTime, indxbvf
      implicit none
      wrtavg(indxbvf) = out_avg_bvf
      if (out_avg_bvf) wrtavg(indxTime) = .true.
   end subroutine init_bvf_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && \
   (defined LMD_SKPP||defined GLS_MIXING)
   subroutine init_hbl_averages_fields()
      use croco_namelist, ONLY: out_avg_hbl
      use ncscrum, ONLY: wrtavg, indxTime, indxHbl
      implicit none
      wrtavg(indxHbl) = out_avg_hbl
      if (out_avg_hbl) wrtavg(indxTime) = .true.
   end subroutine init_hbl_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined LMD_BKPP
   subroutine init_lmd_bkpp_averages_fields()
      use croco_namelist, ONLY: out_avg_hbbl
      use ncscrum, ONLY: wrtavg, indxTime, indxHbbl
      implicit none
      wrtavg(indxHbbl) = out_avg_hbbl
      if (out_avg_hbbl) wrtavg(indxTime) = .true.
   end subroutine init_lmd_bkpp_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined VIS_COEF_3D
   subroutine init_vis_coef_averages_fields()
      use croco_namelist, ONLY: out_avg_visc3d
      use ncscrum, ONLY: wrtavg, indxTime, indxVisc
      implicit none
      wrtavg(indxVisc) = out_avg_visc3d
      if (out_avg_visc3d) wrtavg(indxTime) = .true.
   end subroutine init_vis_coef_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined DIF_COEF_3D
   subroutine init_dif_coef_averages_fields()
      use croco_namelist, ONLY: out_avg_diff3d
      use ncscrum, ONLY: wrtavg, indxTime, indxDiff
      implicit none
      wrtavg(indxDiff) = out_avg_diff3d
      if (out_avg_diff3d) wrtavg(indxTime) = .true.
   end subroutine init_dif_coef_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined BIOLOGY && !defined PISCES
   subroutine init_biology_averages_fields()
      use croco_namelist, ONLY: out_avg_hel
      use ncscrum, ONLY: wrtavg, indxTime, indxHel
      implicit none
      wrtavg(indxHel) = out_avg_hel
      if (out_avg_hel) wrtavg(indxTime) = .true.
   end subroutine init_biology_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined BIO_NChlPZD
   subroutine init_bio_nchlpzd_averages_fields()
      use croco_namelist, ONLY: out_avg_chc, out_avg_u10, out_avg_kvo2, out_avg_o2sat
      use ncscrum, ONLY: wrtavg, indxTime, indxChC
#  ifdef OXYGEN
      use ncscrum, ONLY: indxU10, indxKvO2, indxO2sat
#  endif
      implicit none
      wrtavg(indxChC) = out_avg_chc
      if (out_avg_chc) wrtavg(indxTime) = .true.
#  ifdef OXYGEN
      wrtavg(indxU10) = out_avg_u10
      wrtavg(indxKvO2) = out_avg_kvo2
      wrtavg(indxO2sat) = out_avg_o2sat
      if (out_avg_u10 .or. out_avg_kvo2 .or. out_avg_o2sat) wrtavg(indxTime) = .true.
#  endif
   end subroutine init_bio_nchlpzd_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined BIO_BioEBUS
   subroutine init_bio_bioebus_averages_fields()
      use croco_namelist, ONLY: out_avg_aou, out_avg_wind10
      use ncscrum, ONLY: wrtavg, indxTime, indxAOU, indxWIND10
      implicit none
      wrtavg(indxAOU) = out_avg_aou
      wrtavg(indxWIND10) = out_avg_wind10
      if (out_avg_aou .or. out_avg_wind10) wrtavg(indxTime) = .true.
   end subroutine init_bio_bioebus_averages_fields
# endif

# if defined AVERAGES && defined SOLVE3D && defined MORPHODYN
   !---------------------------------------------------------------------
   !  init_morphodyn_averages_fields
   !  Morphodynamics average field (MORPHODYN).
   !---------------------------------------------------------------------
   subroutine init_morphodyn_averages_fields()
      use croco_namelist, ONLY: out_avg_hm
      use ncscrum, ONLY: wrtavg, indxTime, indxHm
      implicit none

      wrtavg(indxHm) = out_avg_hm
      if (out_avg_hm) wrtavg(indxTime) = .true.
   end subroutine init_morphodyn_averages_fields
# endif

#ifdef DIAGNOSTICS_TS
   subroutine init_diag_ts_fields()
      use croco_namelist, ONLY: out_his_dia3D_tracer
# ifdef DIAGNOSTICS_TS_MLD
      use croco_namelist, ONLY: out_his_dia2D_tracer
# endif
# ifdef AVERAGES
      use croco_namelist, ONLY: out_avg_dia3D_tracer
#  ifdef DIAGNOSTICS_TS_MLD
      use croco_namelist, ONLY: out_avg_dia2D_tracer
#  endif
# endif
      use param, ONLY: NT
      use ncscrum, ONLY: wrtdia3D, wrtdia2D
# ifdef AVERAGES
      use ncscrum, ONLY: wrtdia3D_avg, wrtdia2D_avg
# endif
      implicit none
      integer :: itrc
# ifdef TRACERS
      do itrc = 1, NT
         wrtdia3D(itrc) = out_his_dia3D_tracer(itrc)
         if (out_his_dia3D_tracer(itrc)) wrtdia3D(NT + 1) = .true.
      end do
#  ifdef DIAGNOSTICS_TS_MLD
      do itrc = 1, NT
         wrtdia2D(itrc) = out_his_dia2D_tracer(itrc)
         if (out_his_dia2D_tracer(itrc)) wrtdia2D(NT + 1) = .true.
      end do
#  endif
# endif
# if defined AVERAGES && defined TRACERS
      do itrc = 1, NT
         wrtdia3D_avg(itrc) = out_avg_dia3D_tracer(itrc)
         if (out_avg_dia3D_tracer(itrc)) wrtdia3D_avg(NT + 1) = .true.
      end do
#  ifdef DIAGNOSTICS_TS_MLD
      do itrc = 1, NT
         wrtdia2D_avg(itrc) = out_avg_dia2D_tracer(itrc)
         if (out_avg_dia2D_tracer(itrc)) wrtdia2D_avg(NT + 1) = .true.
      end do
#  endif
# endif
   end subroutine init_diag_ts_fields
#endif

#ifdef DIAGNOSTICS_UV
   subroutine init_diagM_fields()
      use croco_namelist, ONLY: out_his_diagM_u, out_his_diagM_v
# ifdef AVERAGES
      use croco_namelist, ONLY: out_avg_diagM_u, out_avg_diagM_v
      use ncscrum, ONLY: wrtdiaM, wrtdiaM_avg
# else
      use ncscrum, ONLY: wrtdiaM
# endif
      implicit none
      wrtdiaM(1) = out_his_diagM_u
      wrtdiaM(2) = out_his_diagM_v
      if (out_his_diagM_u .or. out_his_diagM_v) wrtdiaM(3) = .true.
# ifdef AVERAGES
      wrtdiaM_avg(1) = out_avg_diagM_u
      wrtdiaM_avg(2) = out_avg_diagM_v
      if (out_avg_diagM_u .or. out_avg_diagM_v) wrtdiaM_avg(3) = .true.
# endif
   end subroutine init_diagM_fields
#endif

#ifdef DIAGNOSTICS_VRT
   subroutine init_diags_vrt_fields()
      use croco_namelist, ONLY: out_his_diags_vrt
# ifdef AVERAGES
      use croco_namelist, ONLY: out_avg_diags_vrt
      use ncscrum, ONLY: wrtdiags_vrt, wrtdiags_vrt_avg
# else
      use ncscrum, ONLY: wrtdiags_vrt
# endif
      implicit none
      wrtdiags_vrt(1) = out_his_diags_vrt
      if (out_his_diags_vrt) wrtdiags_vrt(3) = .true.
# ifdef AVERAGES
      wrtdiags_vrt_avg(1) = out_avg_diags_vrt
      if (out_avg_diags_vrt) wrtdiags_vrt_avg(3) = .true.
# endif
   end subroutine init_diags_vrt_fields
#endif

#ifdef DIAGNOSTICS_KE
   subroutine init_diags_ek_fields()
      use croco_namelist, ONLY: out_his_diags_ek
# ifdef AVERAGES
      use croco_namelist, ONLY: out_avg_diags_ek
      use ncscrum, ONLY: wrtdiags_ek, wrtdiags_ek_avg
# else
      use ncscrum, ONLY: wrtdiags_ek
# endif
      implicit none
      wrtdiags_ek(1) = out_his_diags_ek
      if (out_his_diags_ek) wrtdiags_ek(3) = .true.
# ifdef AVERAGES
      wrtdiags_ek_avg(1) = out_avg_diags_ek
      if (out_avg_diags_ek) wrtdiags_ek_avg(3) = .true.
# endif
   end subroutine init_diags_ek_fields
#endif

#ifdef DIAGNOSTICS_PV
   subroutine init_diags_pv_fields()
# ifdef TRACERS
      use croco_namelist, ONLY: out_his_diags_pv_tracer
      use ncscrum, ONLY: wrtdiags_pv
#  ifdef AVERAGES
      use croco_namelist, ONLY: out_avg_diags_pv_tracer
      use ncscrum, ONLY: wrtdiags_pv_avg
#  endif
      use param, ONLY: NT
# endif
      implicit none
      integer :: itrc
# ifdef TRACERS
      do itrc = 1, NT
         wrtdiags_pv(itrc) = out_his_diags_pv_tracer(itrc)
         if (out_his_diags_pv_tracer(itrc)) wrtdiags_pv(NT + 1) = .true.
      end do
#  ifdef AVERAGES
      do itrc = 1, NT
         wrtdiags_pv_avg(itrc) = out_avg_diags_pv_tracer(itrc)
         if (out_avg_diags_pv_tracer(itrc)) wrtdiags_pv_avg(NT + 1) = .true.
      end do
#  endif
# endif
   end subroutine init_diags_pv_fields
#endif

#if defined DIAGNOSTICS_EDDY && !defined XIOS && defined AVERAGES
   subroutine init_diags_eddy_fields()
      use croco_namelist, ONLY: out_avg_diags_eddy
      use ncscrum, ONLY: wrtdiags_eddy_avg
      implicit none
      wrtdiags_eddy_avg(1) = out_avg_diags_eddy
      if (out_avg_diags_eddy) wrtdiags_eddy_avg(3) = .true.
   end subroutine init_diags_eddy_fields
#endif

#ifdef DIAGNOSTICS_BIO
   subroutine init_diagbio_fields()
      use croco_namelist, ONLY: out_his_diagbioFlux, out_his_diagbioVSink
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      use croco_namelist, ONLY: out_his_diagbioGasExc
# endif
# ifdef AVERAGES
      use croco_namelist, ONLY: out_avg_diagbioFlux, out_avg_diagbioVSink
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      use croco_namelist, ONLY: out_avg_diagbioGasExc
#  endif
# endif
      use param, ONLY: NumFluxTerms, NumVSinkTerms
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      use param, ONLY: NumGasExcTerms
# endif
      use ncscrum, ONLY: wrtdiabioFlux, wrtdiabioVSink
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      use ncscrum, ONLY: wrtdiabioGasExc
# endif
# ifdef AVERAGES
      use ncscrum, ONLY: wrtdiabioFlux_avg, wrtdiabioVSink_avg
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      use ncscrum, ONLY: wrtdiabioGasExc_avg
#  endif
# endif
      implicit none
      integer :: iflux
      do iflux = 1, NumFluxTerms
         wrtdiabioFlux(iflux) = out_his_diagbioFlux(iflux)
         if (out_his_diagbioFlux(iflux)) wrtdiabioFlux(NumFluxTerms + 1) = .true.
      end do
      do iflux = 1, NumVSinkTerms
         wrtdiabioVSink(iflux) = out_his_diagbioVSink(iflux)
         if (out_his_diagbioVSink(iflux)) wrtdiabioVSink(NumVSinkTerms + 1) = .true.
      end do
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      do iflux = 1, NumGasExcTerms
         wrtdiabioGasExc(iflux) = out_his_diagbioGasExc(iflux)
         if (out_his_diagbioGasExc(iflux)) wrtdiabioGasExc(NumGasExcTerms + 1) = .true.
      end do
# endif
# ifdef AVERAGES
      do iflux = 1, NumFluxTerms
         wrtdiabioFlux_avg(iflux) = out_avg_diagbioFlux(iflux)
         if (out_avg_diagbioFlux(iflux)) wrtdiabioFlux_avg(NumFluxTerms + 1) = .true.
      end do
      do iflux = 1, NumVSinkTerms
         wrtdiabioVSink_avg(iflux) = out_avg_diagbioVSink(iflux)
         if (out_avg_diagbioVSink(iflux)) wrtdiabioVSink_avg(NumVSinkTerms + 1) = .true.
      end do
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      do iflux = 1, NumGasExcTerms
         wrtdiabioGasExc_avg(iflux) = out_avg_diagbioGasExc(iflux)
         if (out_avg_diagbioGasExc(iflux)) wrtdiabioGasExc_avg(NumGasExcTerms + 1) = .true.
      end do
#  endif
# endif
   end subroutine init_diagbio_fields
#endif

#ifdef STOGEN
   subroutine init_stochastic_history_fields()
      use croco_namelist, ONLY: out_his_xi2d, out_his_xi3d
      use ncscrum, ONLY: wrthis, indxXIsto2d, indxXIsto3d
      use stoalloc, ONLY: ln_hststo2d, ln_hststo3d
      implicit none
      wrthis(indxXIsto2d) = out_his_xi2d
      wrthis(indxXIsto3d) = out_his_xi3d
      ln_hststo2d = out_his_xi2d
      ln_hststo3d = out_his_xi3d
   end subroutine init_stochastic_history_fields
#endif

#if defined SOLVE3D && defined GLS_MIXING
   subroutine init_gls_history_fields()
      use croco_namelist, ONLY: out_his_tke, out_his_gls, out_his_lscale
      use ncscrum, ONLY: wrthis, indxTke, indxGls, indxLsc, indxTime
      implicit none
      wrthis(indxTke) = out_his_tke
      wrthis(indxGls) = out_his_gls
      wrthis(indxLsc) = out_his_lscale
      if (out_his_tke .or. out_his_gls .or. out_his_lscale) wrthis(indxTime) = .true.
   end subroutine init_gls_history_fields
#endif

#if defined AVERAGES && defined SOLVE3D && defined GLS_MIXING
   subroutine init_gls_averages_fields()
      use croco_namelist, ONLY: out_avg_tke, out_avg_gls, out_avg_lscale
      use ncscrum, ONLY: wrtavg, indxAkk, indxAkp, indxTke, indxGls, &
                         indxLsc, indxTime
      implicit none
      wrtavg(indxTke) = out_avg_tke
      wrtavg(indxGls) = out_avg_gls
      wrtavg(indxLsc) = out_avg_lscale
      if (wrtavg(indxAkk) .or. wrtavg(indxAkp) .or. out_avg_tke &
          .or. out_avg_gls .or. out_avg_lscale) wrtavg(indxTime) = .true.
   end subroutine init_gls_averages_fields
#endif

#if defined ABL1D && !defined XIOS
   subroutine init_abl_history_fields()
      use croco_namelist, ONLY: &
         out_his_abl_pu_dta, out_his_abl_pv_dta, out_his_abl_pt_dta, &
         out_his_abl_pq_dta, out_his_abl_pgu_dta, out_his_abl_pgv_dta, &
         out_his_abl_u_abl, out_his_abl_v_abl, out_his_abl_t_abl, out_his_abl_q_abl, &
         out_his_abl_tke_abl, out_his_abl_mxlm_abl, out_his_abl_mxld_abl, &
         out_his_abl_avm_abl, out_his_abl_avt_abl, out_his_abl_ablh_abl, &
         out_his_abl_zr_abl, out_his_abl_zw_abl, out_his_abl_Hzr_abl, out_his_abl_Hzw_abl
      use ncscrum, ONLY: wrtabl, indxTime, &
                         indxabl_pu_dta, indxabl_pv_dta, indxabl_pt_dta, &
                         indxabl_pq_dta, indxabl_pgu_dta, indxabl_pgv_dta, &
                         indxabl_u_abl, indxabl_v_abl, indxabl_t_abl, indxabl_q_abl, &
                         indxabl_tke_abl, indxabl_mxlm_abl, indxabl_mxld_abl, &
                         indxabl_avm_abl, indxabl_avt_abl, indxabl_ablh_abl, &
                         indxabl_zr_abl, indxabl_zw_abl, indxabl_Hzr_abl, indxabl_Hzw_abl
      implicit none
      wrtabl(indxabl_pu_dta) = out_his_abl_pu_dta
      wrtabl(indxabl_pv_dta) = out_his_abl_pv_dta
      wrtabl(indxabl_pt_dta) = out_his_abl_pt_dta
      wrtabl(indxabl_pq_dta) = out_his_abl_pq_dta
      wrtabl(indxabl_pgu_dta) = out_his_abl_pgu_dta
      wrtabl(indxabl_pgv_dta) = out_his_abl_pgv_dta
      wrtabl(indxabl_u_abl) = out_his_abl_u_abl
      wrtabl(indxabl_v_abl) = out_his_abl_v_abl
      wrtabl(indxabl_t_abl) = out_his_abl_t_abl
      wrtabl(indxabl_q_abl) = out_his_abl_q_abl
      wrtabl(indxabl_tke_abl) = out_his_abl_tke_abl
      wrtabl(indxabl_mxlm_abl) = out_his_abl_mxlm_abl
      wrtabl(indxabl_mxld_abl) = out_his_abl_mxld_abl
      wrtabl(indxabl_avm_abl) = out_his_abl_avm_abl
      wrtabl(indxabl_avt_abl) = out_his_abl_avt_abl
      wrtabl(indxabl_ablh_abl) = out_his_abl_ablh_abl
      wrtabl(indxabl_zr_abl) = out_his_abl_zr_abl
      wrtabl(indxabl_zw_abl) = out_his_abl_zw_abl
      wrtabl(indxabl_Hzr_abl) = out_his_abl_Hzr_abl
      wrtabl(indxabl_Hzw_abl) = out_his_abl_Hzw_abl
      if (out_his_abl_pu_dta .or. out_his_abl_pv_dta .or. &
          out_his_abl_pt_dta .or. out_his_abl_pq_dta .or. &
          out_his_abl_pgu_dta .or. out_his_abl_pgv_dta .or. &
          out_his_abl_u_abl .or. out_his_abl_v_abl .or. &
          out_his_abl_t_abl .or. out_his_abl_q_abl .or. &
          out_his_abl_tke_abl .or. out_his_abl_mxlm_abl .or. &
          out_his_abl_mxld_abl .or. out_his_abl_avm_abl .or. &
          out_his_abl_avt_abl .or. out_his_abl_ablh_abl .or. &
          out_his_abl_zr_abl .or. out_his_abl_zw_abl .or. &
          out_his_abl_Hzr_abl .or. out_his_abl_Hzw_abl) wrtabl(indxTime) = .true.
   end subroutine init_abl_history_fields

# ifdef AVERAGES
   subroutine init_abl_averages_fields()
      use croco_namelist, ONLY: &
         out_avg_abl_pu_dta, out_avg_abl_pv_dta, out_avg_abl_pt_dta, &
         out_avg_abl_pq_dta, out_avg_abl_pgu_dta, out_avg_abl_pgv_dta, &
         out_avg_abl_u_abl, out_avg_abl_v_abl, out_avg_abl_t_abl, out_avg_abl_q_abl, &
         out_avg_abl_tke_abl, out_avg_abl_mxlm_abl, out_avg_abl_mxld_abl, &
         out_avg_abl_avm_abl, out_avg_abl_avt_abl, out_avg_abl_ablh_abl, &
         out_avg_abl_zr_abl, out_avg_abl_zw_abl, out_avg_abl_Hzr_abl, out_avg_abl_Hzw_abl
      use ncscrum, ONLY: wrtabl_avg, indxTime, &
                         indxabl_pu_dta, indxabl_pv_dta, indxabl_pt_dta, &
                         indxabl_pq_dta, indxabl_pgu_dta, indxabl_pgv_dta, &
                         indxabl_u_abl, indxabl_v_abl, indxabl_t_abl, indxabl_q_abl, &
                         indxabl_tke_abl, indxabl_mxlm_abl, indxabl_mxld_abl, &
                         indxabl_avm_abl, indxabl_avt_abl, indxabl_ablh_abl, &
                         indxabl_zr_abl, indxabl_zw_abl, indxabl_Hzr_abl, indxabl_Hzw_abl
      implicit none
      wrtabl_avg(indxabl_pu_dta) = out_avg_abl_pu_dta
      wrtabl_avg(indxabl_pv_dta) = out_avg_abl_pv_dta
      wrtabl_avg(indxabl_pt_dta) = out_avg_abl_pt_dta
      wrtabl_avg(indxabl_pq_dta) = out_avg_abl_pq_dta
      wrtabl_avg(indxabl_pgu_dta) = out_avg_abl_pgu_dta
      wrtabl_avg(indxabl_pgv_dta) = out_avg_abl_pgv_dta
      wrtabl_avg(indxabl_u_abl) = out_avg_abl_u_abl
      wrtabl_avg(indxabl_v_abl) = out_avg_abl_v_abl
      wrtabl_avg(indxabl_t_abl) = out_avg_abl_t_abl
      wrtabl_avg(indxabl_q_abl) = out_avg_abl_q_abl
      wrtabl_avg(indxabl_tke_abl) = out_avg_abl_tke_abl
      wrtabl_avg(indxabl_mxlm_abl) = out_avg_abl_mxlm_abl
      wrtabl_avg(indxabl_mxld_abl) = out_avg_abl_mxld_abl
      wrtabl_avg(indxabl_avm_abl) = out_avg_abl_avm_abl
      wrtabl_avg(indxabl_avt_abl) = out_avg_abl_avt_abl
      wrtabl_avg(indxabl_ablh_abl) = out_avg_abl_ablh_abl
      wrtabl_avg(indxabl_zr_abl) = out_avg_abl_zr_abl
      wrtabl_avg(indxabl_zw_abl) = out_avg_abl_zw_abl
      wrtabl_avg(indxabl_Hzr_abl) = out_avg_abl_Hzr_abl
      wrtabl_avg(indxabl_Hzw_abl) = out_avg_abl_Hzw_abl
      if (out_avg_abl_pu_dta .or. out_avg_abl_pv_dta .or. &
          out_avg_abl_pt_dta .or. out_avg_abl_pq_dta .or. &
          out_avg_abl_pgu_dta .or. out_avg_abl_pgv_dta .or. &
          out_avg_abl_u_abl .or. out_avg_abl_v_abl .or. &
          out_avg_abl_t_abl .or. out_avg_abl_q_abl .or. &
          out_avg_abl_tke_abl .or. out_avg_abl_mxlm_abl .or. &
          out_avg_abl_mxld_abl .or. out_avg_abl_avm_abl .or. &
          out_avg_abl_avt_abl .or. out_avg_abl_ablh_abl .or. &
          out_avg_abl_zr_abl .or. out_avg_abl_zw_abl .or. &
          out_avg_abl_Hzr_abl .or. out_avg_abl_Hzw_abl) wrtabl_avg(indxTime) = .true.
   end subroutine init_abl_averages_fields
# endif
#endif

#if defined OUTPUTS_SURFACE && !defined XIOS
   subroutine init_surf_history_fields()
      use croco_namelist, ONLY: out_his_surf
      use ncscrum, ONLY: wrtsurf
      implicit none
      wrtsurf(1) = out_his_surf
      if (out_his_surf) wrtsurf(3) = .true.
   end subroutine init_surf_history_fields

# ifdef AVERAGES
   subroutine init_surf_average_fields()
      use croco_namelist, ONLY: out_avg_surf
      use ncscrum, ONLY: wrtsurf_avg
      implicit none
      wrtsurf_avg(1) = out_avg_surf
      if (out_avg_surf) wrtsurf_avg(3) = .true.
   end subroutine init_surf_average_fields
# endif
#endif

#ifdef STATIONS
   subroutine init_station_fields()
      use croco_namelist, ONLY: sta_grd, sta_temp, sta_salt, sta_rho, sta_vel
#ifdef AGRIF
      use Agrif_Util
#endif
#include "nc_sta.h"
      implicit none
#ifdef AGRIF
      if (.not. Agrif_Root()) return
#endif
      wrtsta(indxstaGrd) = sta_grd
      wrtsta(indxstaTemp) = sta_temp
      wrtsta(indxstaSalt) = sta_salt
      wrtsta(indxstaRho) = sta_rho
      wrtsta(indxstaVel) = sta_vel
   end subroutine init_station_fields
#endif

#if defined SOLVE3D && defined SEDIMENT
   subroutine init_sediment_history_fields()
      use croco_namelist, ONLY: out_his_sed_athk, out_his_sed_bthk, out_his_sed_bpor, &
                                out_his_sed_bfra
# ifdef SUSPLOAD
      use croco_namelist, ONLY: out_his_sed_dflx, out_his_sed_eflx
# endif
# ifdef BEDLOAD
      use croco_namelist, ONLY: out_his_sed_bdlu, out_his_sed_bdlv
# endif
# if defined MIXED_BED || defined COHESIVE_BED
      use croco_namelist, ONLY: out_his_sed_btcr
# endif
# ifdef AVERAGES
      use ncscrum, ONLY: wrtavg
# endif
      use ncscrum, ONLY: wrthis, indxATHK, indxBTHK, indxBPOR, indxBFRA, &
                         indxSed
# ifdef SUSPLOAD
      use ncscrum, ONLY: indxDFLX, indxEFLX
# endif
# ifdef BEDLOAD
      use ncscrum, ONLY: indxBDLU, indxBDLV
# endif
# if defined MIXED_BED || defined COHESIVE_BED
      use ncscrum, ONLY: indxBTCR
# endif
      use param, ONLY: NST
      implicit none
      integer :: itrc

      wrthis(indxATHK) = out_his_sed_athk
      wrthis(indxBTHK) = out_his_sed_bthk
      wrthis(indxBPOR) = out_his_sed_bpor
      do itrc = 1, NST
         wrthis(indxBFRA(itrc)) = out_his_sed_bfra(itrc)
      end do
# ifdef SUSPLOAD
      do itrc = 1, NST
         wrthis(indxDFLX(itrc)) = out_his_sed_dflx(itrc)
         wrthis(indxEFLX(itrc)) = out_his_sed_eflx(itrc)
      end do
# endif
# ifdef BEDLOAD
      do itrc = 1, NST
         wrthis(indxBDLU(itrc)) = out_his_sed_bdlu(itrc)
         wrthis(indxBDLV(itrc)) = out_his_sed_bdlv(itrc)
      end do
# endif
# if defined MIXED_BED || defined COHESIVE_BED
      wrthis(indxBTCR) = out_his_sed_btcr
# endif
# ifdef AVERAGES
      do itrc = indxSed, &
# ifdef SUSPLOAD
# ifdef BEDLOAD
         indxSed + NST + 1 + 4*NST
# else
         indxSed + NST + 1 + 2*NST
# endif
# else
# ifdef BEDLOAD
         indxSed + NST + 1 + 2*NST
# else
         indxSed + NST + 1
# endif
# endif
         wrtavg(itrc) = wrthis(itrc)
      end do
# endif
   end subroutine init_sediment_history_fields
#endif

#ifdef BBL
   subroutine init_bbl_history_fields()
      use croco_namelist, ONLY: out_his_abed, out_his_hripple, out_his_lripple, &
                                out_his_zbnot, out_his_zbapp, out_his_bostrw
      use ncscrum, ONLY: wrthis, indxAbed, indxHrip, indxLrip, &
                         indxZbnot, indxZbapp, indxBostrw
      implicit none
      wrthis(indxAbed) = out_his_abed
      wrthis(indxHrip) = out_his_hripple
      wrthis(indxLrip) = out_his_lripple
      wrthis(indxZbnot) = out_his_zbnot
      wrthis(indxZbapp) = out_his_zbapp
      wrthis(indxBostrw) = out_his_bostrw
   end subroutine init_bbl_history_fields
#endif

#ifdef MRL_WCI
   subroutine init_wci_history_fields()
      use croco_namelist, ONLY: out_his_sup, out_his_ust2d, out_his_vst2d
# ifdef SOLVE3D
      use croco_namelist, ONLY: out_his_ust, out_his_vst, out_his_wst, out_his_akb, &
                                out_his_akw, out_his_kvf, out_his_calp, out_his_kaps
# endif
      use ncscrum, ONLY: wrthis, indxSUP, indxUST2D, indxVST2D
# ifdef SOLVE3D
      use ncscrum, ONLY: indxUST, indxVST, indxWST, indxAkb, indxAkw, &
                         indxKVF, indxCALP, indxKAPS
# endif
      implicit none
      wrthis(indxSUP) = out_his_sup
      wrthis(indxUST2D) = out_his_ust2d
      wrthis(indxVST2D) = out_his_vst2d
# ifdef SOLVE3D
      wrthis(indxUST) = out_his_ust
      wrthis(indxVST) = out_his_vst
      wrthis(indxWST) = out_his_wst
      wrthis(indxAkb) = out_his_akb
      wrthis(indxAkw) = out_his_akw
      wrthis(indxKVF) = out_his_kvf
      wrthis(indxCALP) = out_his_calp
      wrthis(indxKAPS) = out_his_kaps
# endif
   end subroutine init_wci_history_fields

# ifdef AVERAGES
   subroutine init_wci_average_fields()
      use croco_namelist, ONLY: out_avg_sup, out_avg_ust2d, out_avg_vst2d
#  ifdef SOLVE3D
      use croco_namelist, ONLY: out_avg_ust, out_avg_vst, out_avg_wst, out_avg_akb, &
                                out_avg_akw, out_avg_kvf, out_avg_calp, out_avg_kaps
#  endif
      use ncscrum, ONLY: wrtavg, indxSUP, indxUST2D, indxVST2D
#  ifdef SOLVE3D
      use ncscrum, ONLY: indxUST, indxVST, indxWST, indxAkb, indxAkw, &
                         indxKVF, indxCALP, indxKAPS
#  endif
      implicit none
      wrtavg(indxSUP) = out_avg_sup
      wrtavg(indxUST2D) = out_avg_ust2d
      wrtavg(indxVST2D) = out_avg_vst2d
#  ifdef SOLVE3D
      wrtavg(indxUST) = out_avg_ust
      wrtavg(indxVST) = out_avg_vst
      wrtavg(indxWST) = out_avg_wst
      wrtavg(indxAkb) = out_avg_akb
      wrtavg(indxAkw) = out_avg_akw
      wrtavg(indxKVF) = out_avg_kvf
      wrtavg(indxCALP) = out_avg_calp
      wrtavg(indxKAPS) = out_avg_kaps
#  endif
   end subroutine init_wci_average_fields
# endif
#endif

#if defined MRL_WCI || defined OW_COUPLING
   subroutine init_wave_history_fields()
      use croco_namelist, ONLY: out_his_hrm, out_his_frq, out_his_action, out_his_k_xi, &
                                out_his_k_eta, out_his_eps_b, out_his_eps_d, out_his_erol, out_his_eps_r
      use ncscrum, ONLY: wrthis, indxHRM, indxFRQ, indxWAC, indxWKX, &
                         indxWKE, indxEPB, indxEPD, indxWAR, indxEPR
      implicit none
      wrthis(indxHRM) = out_his_hrm
      wrthis(indxFRQ) = out_his_frq
      wrthis(indxWAC) = out_his_action
      wrthis(indxWKX) = out_his_k_xi
      wrthis(indxWKE) = out_his_k_eta
      wrthis(indxEPB) = out_his_eps_b
      wrthis(indxEPD) = out_his_eps_d
      wrthis(indxWAR) = out_his_erol
      wrthis(indxEPR) = out_his_eps_r
   end subroutine init_wave_history_fields

# ifdef AVERAGES
   subroutine init_wave_average_fields()
      use croco_namelist, ONLY: out_avg_hrm, out_avg_frq, out_avg_action, out_avg_k_xi, &
                                out_avg_k_eta, out_avg_eps_b, out_avg_eps_d, out_avg_erol, out_avg_eps_r
      use ncscrum, ONLY: wrtavg, indxHRM, indxFRQ, indxWAC, indxWKX, &
                         indxWKE, indxEPB, indxEPD, indxWAR, indxEPR
      implicit none
      wrtavg(indxHRM) = out_avg_hrm
      wrtavg(indxFRQ) = out_avg_frq
      wrtavg(indxWAC) = out_avg_action
      wrtavg(indxWKX) = out_avg_k_xi
      wrtavg(indxWKE) = out_avg_k_eta
      wrtavg(indxEPB) = out_avg_eps_b
      wrtavg(indxEPD) = out_avg_eps_d
      wrtavg(indxWAR) = out_avg_erol
      wrtavg(indxEPR) = out_avg_eps_r
   end subroutine init_wave_average_fields
# endif
#endif

#if defined ANA_PSOURCE
   !---------------------------------------------------------------------
   !  init_psource
   !  Copy psource_* namelist vars into sources.h COMMON block variables,
   !  verify qbarname file accessibility (PSOURCE_NCFILE), and set up
   !  the MPI tile-relative source indices.
   !---------------------------------------------------------------------
   subroutine init_psource(ierr)
      use param, ONLY: stdout, Msrc
      use ncscrum, ONLY: qbarname
      use croco_namelist, ONLY: psource_Nsrc, psource_Isrc, psource_Jsrc, psource_Dsrc
#  ifndef PSOURCE_NCFILE
      use croco_namelist, ONLY: psource_Qbar
#  endif
#  ifdef PSOURCE_NCFILE
      use croco_namelist, ONLY: psource_qbarname, psource_qbardir
#  endif
#  if defined TRACERS
      use param, ONLY: NT
      use croco_namelist, ONLY: psource_Lsrc, psource_Tsrc0
#  endif
#  ifdef MPI
      use param, ONLY: iminmpi, jminmpi, NNODES
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: is
#  if defined TRACERS
      integer :: itrc
#  endif
      integer :: ios

      integer :: Nsrc
      integer :: Isrc(Msrc), Jsrc(Msrc), Dsrc(Msrc)
      real    :: Qbar(Msrc)
      common/source_Nsrc/Nsrc
      common/source_Isrc/Isrc
      common/source_Jsrc/Jsrc
      common/source_Dsrc/Dsrc
      common/sources_Qbar/Qbar
#  if defined TRACERS
      logical :: Lsrc(Msrc, NT)
      real    :: Tsrc0(Msrc, NT)
      common/source_Lsrc/Lsrc
      common/source_Tsrc0/Tsrc0
#  endif
#  ifdef PSOURCE_NCFILE
      real :: qbardir(Msrc)
      common/source_qbardir/qbardir
#  endif
#  ifdef MPI
      integer :: Isrc_mpi(Msrc, 0:NNODES - 1)
      integer :: Jsrc_mpi(Msrc, 0:NNODES - 1)
      common/source_Isrc_mpi/Isrc_mpi
      common/source_Jsrc_mpi/Jsrc_mpi
#  endif

#  ifdef PSOURCE_NCFILE
      open (testunit, file=trim(psource_qbarname), status='old', iostat=ios)
      if (ios == 0) then
         close (testunit)
      else
         MPI_master_only write (stdout, *) &
            'Error: cannot open runoff file ', trim(psource_qbarname)
         ierr = ierr + 1
         return
      end if
      qbarname = psource_qbarname
#  endif

      Nsrc = psource_Nsrc
      do is = 1, Nsrc
         Isrc(is) = psource_Isrc(is)
         Jsrc(is) = psource_Jsrc(is)
         Dsrc(is) = psource_Dsrc(is)
#  ifndef PSOURCE_NCFILE
         Qbar(is) = psource_Qbar(is)
#  else
         qbardir(is) = psource_qbardir(is)
#  endif
#  if defined TRACERS
         do itrc = 1, NT
            Lsrc(is, itrc) = psource_Lsrc(is, itrc)
            Tsrc0(is, itrc) = psource_Tsrc0(is, itrc)
         end do
#  endif
      end do

#  ifdef MPI
      do is = 1, Nsrc
         Isrc_mpi(is, mynode) = Isrc(is) - iminmpi + 1
         Jsrc_mpi(is, mynode) = Jsrc(is) - jminmpi + 1
      end do
#  endif

   end subroutine init_psource
#endif /* ANA_PSOURCE */

END MODULE croco_namelist_init
