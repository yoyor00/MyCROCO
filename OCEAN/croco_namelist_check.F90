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
!  Individual check routines are kept public for unit-testing purposes.
!=======================================================================

MODULE croco_namelist_check
   implicit none
   private

   public :: check_all
   public :: check_time_stepping
   public :: check_history
#ifdef NBQ
   public :: check_time_stepping_nbq
#endif
#ifdef SOLVE3D
   public :: check_croco_s_coord
#endif
#if defined WKB_WWAVE && defined WAVE_ROLLER
   public :: check_wkb_roller
#endif
#ifdef BODYFORCE
   public :: check_bodyforce
#endif

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
#if defined WKB_WWAVE && defined WAVE_ROLLER
      call check_wkb_roller(ierr)
#endif
#ifdef BODYFORCE
      call check_bodyforce(ierr)
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

#ifdef NBQ
   !---------------------------------------------------------------------
   !  check_time_stepping_nbq
   !---------------------------------------------------------------------
   subroutine check_time_stepping_nbq(ierr)
      use param, ONLY: stdout
      use croco_namelist, ONLY: ndtnbq, csound_nbq, visc2_nbq
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

      if (wkb_roller < 0.d0 .or. wkb_roller > 1.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') &
            'Error - wkb_roller = ', wkb_roller, &
            ' must be in [0, 1].'
         ierr = ierr + 1
      end if

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

END MODULE croco_namelist_check
