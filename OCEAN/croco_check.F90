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
!  MODULE croco_check
!
!  Purpose : Consistency checks on namelist parameters.
!            Each subroutine writes its own diagnostic messages and
!            increments ierr for every fatal error found.
!            The caller (read_nml) decides what to do with ierr.
!=======================================================================

MODULE croco_check
   implicit none
   private

   public :: check_time_stepping
   public :: check_history
#ifdef NBQ
   public :: check_time_stepping_nbq
#endif

contains

   !---------------------------------------------------------------------
   !  check_time_stepping
   !  Validate &croco_time_stepping parameters.
   !---------------------------------------------------------------------
   subroutine check_time_stepping(ierr)
      use param, ONLY: stdout, NWEIGHT
      use croco_namelist, ONLY: dt, ntimes, ndtfast
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
      use croco_namelist, ONLY: nwrt, nrpfhis
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
#endif

END MODULE croco_check
