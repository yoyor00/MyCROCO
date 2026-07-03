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
!  MODULE croco_namelist_check
!
!  Purpose : Validation routines for all namelist sections.
!
!  Public entry points:
!    check_all              – call every check_xxx in order, accumulate
!                             all errors before returning
!
!  Individual check routines are not used outside this module and are
!  kept private.
!=======================================================================

MODULE croco_namelist_check
   implicit none
   private

   public :: check_all

contains

   !---------------------------------------------------------------------
   !  check_all
   !  Single entry point that runs every check_xxx in sequence.
   !  All errors are accumulated so that the user sees the full list
   !  in one run rather than stopping at the first failing section.
   !---------------------------------------------------------------------
   subroutine check_all(ierr)
      implicit none
      integer, intent(inout) :: ierr

      call check_time_stepping(ierr)
#ifdef NBQ
      call check_time_stepping_nbq(ierr)
#endif
#ifdef SOLVE3D
      call check_croco_s_coord(ierr)
#endif
      call check_history(ierr)
#ifdef ONLINE
      call check_online(ierr)
#endif
#if defined WKB_WWAVE && defined WAVE_ROLLER
      call check_wkb_roller(ierr)
#endif
#ifdef BODYFORCE
      call check_bodyforce(ierr)
#endif
#ifdef WET_DRY
      call check_wetdry(ierr)
#endif
#ifdef ANA_PSOURCE
      call check_psource(ierr)
#endif

   end subroutine check_all

   !---------------------------------------------------------------------
   !  check_time_stepping
   !---------------------------------------------------------------------
   subroutine check_time_stepping(ierr)
      use param, ONLY: stdout, NWEIGHT
      use croco_namelist, ONLY: dt, ntimes, ndtfast
#if defined MPI
      use scalars, ONLY: mynode
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

#ifdef NBQ
   !---------------------------------------------------------------------
   !  check_time_stepping_nbq
   !---------------------------------------------------------------------
   subroutine check_time_stepping_nbq(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: csound_nbq, visc2_nbq
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr

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
#endif /* NBQ */

#ifdef SOLVE3D
   !---------------------------------------------------------------------
   !  check_nml_croco_s_coord
   !---------------------------------------------------------------------
   subroutine check_croco_s_coord(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: theta_s, theta_b, Tcline
#  if defined MPI
      use scalars, ONLY: mynode
#  endif
      implicit none
      integer, intent(inout) :: ierr

      if (theta_s < 0.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') &
            'Error - S-coord surface control parameter theta_s = ', theta_s, &
            ' must be positive.'
         ierr = ierr + 1
      end if
      if (theta_s > 20.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') &
            'Error - S-coord surface control parameter theta_s = ', theta_s, &
            ' is too large. It should not exceed 20.'
         ierr = ierr + 1
      end if
      if (theta_b < 0.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') &
            'Error - S-coord bottom control parameter theta_b = ', theta_b, &
            ' must be positive.'
         ierr = ierr + 1
      end if
      if (Tcline < 0.d0) then
         MPI_master_only write (stdout, *) &
            'Error - S-coordinate surface/bottom layer width used in', &
            'vertical coordinate stretching Tcline = ', Tcline, &
            ' must be positive.'
         ierr = ierr + 1
      end if

   end subroutine check_croco_s_coord
#endif /* SOLVE3D */

   !---------------------------------------------------------------------
   !  check_history
   !---------------------------------------------------------------------
   subroutine check_history(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: nwrt, nrpfhis
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr

      if (nwrt <= 0) then
         MPI_master_only write (stdout, *) &
            'Error - History output frequency nwrt must be > 0: ', nwrt
         ierr = ierr + 1
      end if

      if (nrpfhis < -1) then
         MPI_master_only write (stdout, *) &
            'Error - nrpfhis must be >= -1: ', nrpfhis
         ierr = ierr + 1
      end if

   end subroutine check_history

#ifdef ONLINE
   !---------------------------------------------------------------------
   !  check_online
   !---------------------------------------------------------------------
   subroutine check_online(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: pathbulk
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr

      if (len_trim(pathbulk) == 0) then
         MPI_master_only write (stdout, '(a)') &
            'Error - pathbulk is empty in &croco_online.'
         ierr = ierr + 1
      end if

   end subroutine check_online
#endif /* ONLINE */

#if defined WKB_WWAVE && defined WAVE_ROLLER
   !---------------------------------------------------------------------
   !  check_wkb_roller
   !---------------------------------------------------------------------
   subroutine check_wkb_roller(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: wkb_roller
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr

      ! if (wkb_roller < 0.d0 .or. wkb_roller > 1.d0) then
      !    MPI_master_only write (stdout, '(a,f12.4,a)') &
      !       'Error - wkb_roller = ', wkb_roller, &
      !       ' must be in [0, 1].'
      !    ierr = ierr + 1
      ! end if
      ! TODO Check this, in SANDBAR, wkb_roller is 1.8..., cf work_item #533
      ! not coherent with comment in read_inp.F : 
      !     # ifdef WAVE_ROLLER
      !           MPI_master_only write(stdout,'(/1x,A,2(/3x,f10.4,2x,A))')
      !      &  'Svenden (1984) surface roller model parameters.',
      !      &   wkb_rsb,   'wkb_rsb     sin(beta) roller dissipation',
      !      &   wkb_roller,'wkb_roller  breaking contrib to roller: [0,1]'
      !     # endif



   end subroutine check_wkb_roller
#endif /* WKB_WWAVE && WAVE_ROLLER */

#ifdef BODYFORCE
   !---------------------------------------------------------------------
   !  check_bodyforce
   !---------------------------------------------------------------------
   subroutine check_bodyforce(ierr)
      use param, ONLY: stdout, N
      use croco_namelist, ONLY: levsfrc, levbfrc
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr

      if (levsfrc < 1 .or. levsfrc > N) then
         MPI_master_only write (stdout, '(a,i4,a)') &
            'Error - levsfrc = ', levsfrc, &
            ' must be in [1, N].'
         ierr = ierr + 1
      end if
      if (levbfrc < 1 .or. levbfrc > N) then
         MPI_master_only write (stdout, '(a,i4,a)') &
            'Error - levbfrc = ', levbfrc, &
            ' must be in [1, N].'
         ierr = ierr + 1
      end if

   end subroutine check_bodyforce
#endif /* BODYFORCE */

#ifdef WET_DRY
   !---------------------------------------------------------------------
   !  check_wetdry
   !---------------------------------------------------------------------
   subroutine check_wetdry(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: D_wetdry
#if defined MRL_WCI && defined WAVE_DRY
      use croco_namelist, ONLY: D_wavedry
#endif
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr

      if (D_wetdry <= 0.0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') &
            'Error - Critical wetting/drying depth D_wetdry = ', D_wetdry, &
            ' must be positive.'
         ierr = ierr + 1
      end if

#if defined MRL_WCI && defined WAVE_DRY
      if (D_wavedry <= D_wetdry) then
         MPI_master_only write (stdout, '(a,f12.4,a,f12.4)') &
            'Error - Wave-drying depth D_wavedry = ', D_wavedry, &
            ' must exceed critical wetting/drying depth D_wetdry = ', D_wetdry
         ierr = ierr + 1
      end if
#endif

   end subroutine check_wetdry
#endif /* WET_DRY */

#ifdef ANA_PSOURCE
   !---------------------------------------------------------------------
   !  check_psource
   !---------------------------------------------------------------------
   subroutine check_psource(ierr)
      use param, ONLY: stdout, Msrc
      use croco_namelist, ONLY: psource_Nsrc
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr

      if (psource_Nsrc < 0 .or. psource_Nsrc > Msrc) then
         MPI_master_only write (stdout, '(a,i4,a,i4,a)') &
            'Error - psource_Nsrc = ', psource_Nsrc, ' must be in [0, Msrc=', Msrc, '].'
         ierr = ierr + 1
      end if
   end subroutine check_psource
#endif /* ANA_PSOURCE */

END MODULE croco_namelist_check
