MODULE stoalloc

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE stoalloc  ***
   !!
   !! Define and allocate arrays needed:
   !!    for the outout of stochastic fields in history file
   !!    for temporary use by stochastic parameterizations
   !! 
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   sto_alloc_init     : define required arrays
   !!----------------------------------------------------------------------

   USE stoexternal , only : wp
   USE stoarray

   IMPLICIT NONE
   PRIVATE

   PUBLIC sto_alloc_init

   LOGICAL, PUBLIC :: ln_hststo2d = .FALSE.  ! do we output 2d stochastic field
   LOGICAL, PUBLIC :: ln_hststo3d = .FALSE.  ! do we output 3d stochastic field

   ! Working arrays for output in history files
   REAL(wp), PUBLIC, DIMENSION(:,:),   ALLOCATABLE :: sto_xi2d ! output 2d array
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: sto_xi3d ! output 3d array

   ! Working arrays for the stochastic parameterization of vertical mixing
   REAL(wp), PUBLIC, DIMENSION(:,:),   ALLOCATABLE :: sto_p2d ! perturbation of production term
   REAL(wp), PUBLIC, DIMENSION(:,:),   ALLOCATABLE :: sto_b2d ! perturbation of destruction term

   ! Value used to initialize working arrays
   REAL(wp), PARAMETER :: init = 0.

CONTAINS

   SUBROUTINE sto_alloc_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_alloc_init  ***
      !!
      !! Purpose : initialize required arrays
      !!
      !!----------------------------------------------------------------------
      USE param

      LOGICAL :: alloc_xi2d, alloc_xi3d

      ! allocate and initialize arrays for output in history files
      alloc_xi2d = ln_stogen .AND. ln_hststo2d  ! allocate if output of stochastic array is requested
      alloc_xi2d = alloc_xi2d .OR. (ln_stogen.AND.ln_stobulk)    ! needed in any case in stobulk
      alloc_xi2d = alloc_xi2d .OR. (ln_stogen.AND.ln_stostress)  ! needed in any case in stostress

      alloc_xi3d = ln_stogen .AND. ln_hststo3d  ! allocate if output of stochastic array is requested

      IF (alloc_xi2d) THEN
        ! allocate 2D array for output
        allocate(sto_xi2d(GLOBAL_2D_ARRAY))
        sto_xi2d(:,:) = init
      ENDIF

      IF (alloc_xi3d) THEN
        ! allocate 2D array for output
        allocate(sto_xi3d(GLOBAL_2D_ARRAY,N))
        sto_xi3d(:,:,:) = init
      ENDIF

      ! allocate working arrays for vertical mixing
      IF (ln_stogen.AND.ln_stogls) THEN
        allocate(sto_p2d(GLOBAL_2D_ARRAY))
        allocate(sto_b2d(GLOBAL_2D_ARRAY))
      ENDIF

   END SUBROUTINE sto_alloc_init

   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stoalloc
