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
   USE stobulk
   USE stogls

   IMPLICIT NONE
   PRIVATE

   PUBLIC sto_alloc_init

   LOGICAL, PUBLIC :: ln_hststo2d = .FALSE.  ! do we output 2d stochastic field
   LOGICAL, PUBLIC :: ln_hststo3d = .FALSE.  ! do we output 3d stochastic field

   ! Working arrays for output in history files
   REAL(wp), PUBLIC, DIMENSION(:,:),   POINTER :: sto_xi2d ! output 2d array
   REAL(wp), PUBLIC, DIMENSION(:,:,:), POINTER :: sto_xi3d ! output 3d array

   ! Working arrays for the stochastic parameterization of bulk fluxes
   REAL(wp), PUBLIC, DIMENSION(:,:),   POINTER :: sto_bulk_cd ! perturbation of Cd

   ! Working arrays for the stochastic parameterization of the surface stress (non-bulk)
   REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, TARGET :: sto_stress_factor ! perturbation factor

   ! Working arrays for the stochastic parameterization of vertical mixing
   REAL(wp), PUBLIC, DIMENSION(:,:),   POINTER :: sto_gls_s2d ! perturbation of production term
   REAL(wp), PUBLIC, DIMENSION(:,:),   POINTER :: sto_gls_b2d ! perturbation of destruction term
   REAL(wp), PUBLIC, DIMENSION(:,:,:), POINTER :: sto_gls_z3d ! perturbation of destruction term

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

      LOGICAL :: output_xi2d, output_xi3d

      ! define working arrays used for bulk fluxes
      IF (ln_stogen.AND.ln_stobulk) THEN
        sto_bulk_cd(GLOBAL_2D_ARRAY) => sto2d(:,:,jpidxlast2d,stofields(jstobulk_cd)%index)
      ENDIF

      ! define working arrays used for the stochastic surface stress
      IF (ln_stogen.AND.ln_stostress) THEN
        allocate(sto_stress_factor(GLOBAL_2D_ARRAY))
        sto_stress_factor(:,:) = init
      ENDIF

      ! define working arrays used for vertical mixing
      IF (ln_stogen.AND.ln_stogls) THEN
        IF (ln_Sprod) sto_gls_s2d(GLOBAL_2D_ARRAY) => sto2d(:,:,jpidxlast2d,stofields(jstogls_s)%index)
        IF (ln_Bprod) sto_gls_b2d(GLOBAL_2D_ARRAY) => sto2d(:,:,jpidxlast2d,stofields(jstogls_b)%index)
        IF (ln_zlevs) sto_gls_z3d(GLOBAL_2D_ARRAY,1:N) => sto3d(:,:,:,jpidxlast3d,stofields(jstogls_z)%index)
      ENDIF

      ! allocate and initialize arrays for output in history files
      output_xi2d = ln_stogen .AND. ln_hststo2d  ! allocate if output of stochastic array is requested
      output_xi3d = ln_stogen .AND. ln_hststo3d  ! allocate if output of stochastic array is requested

      ! decide which field to output as xi2d (according to namelist)
      IF (output_xi2d) THEN
        SELECT CASE (cn_xi2d)
        CASE('bulk_cd')
          IF (ln_stobulk)             sto_xi2d(GLOBAL_2D_ARRAY) => sto_bulk_cd(GLOBAL_2D_ARRAY)
        CASE('stress')
          IF (ln_stostress)           sto_xi2d(GLOBAL_2D_ARRAY) => sto_stress_factor(GLOBAL_2D_ARRAY)
        CASE('gls_Sprod')
          IF (ln_stogls.AND.ln_Sprod) sto_xi2d(GLOBAL_2D_ARRAY) => sto_gls_s2d(GLOBAL_2D_ARRAY)
        CASE('gls_Bprod')
          IF (ln_stogls.AND.ln_Bprod) sto_xi2d(GLOBAL_2D_ARRAY) => sto_gls_b2d(GLOBAL_2D_ARRAY)
        END SELECT
        IF (.NOT.associated(sto_xi2d)) THEN
          STOP 'No valid array associated to requested xi2d output'
        ENDIF
      ENDIF

      ! decide which field to output as xi3d (according to namelist)
      IF (output_xi3d) THEN
        SELECT CASE (cn_xi3d)
        CASE('gls_zlevs')
          IF (ln_stogls.AND.ln_zlevs) sto_xi3d(GLOBAL_2D_ARRAY,1:N) => sto_gls_z3d(GLOBAL_2D_ARRAY,1:N)
        END SELECT
        IF (.NOT.associated(sto_xi3d)) THEN
          STOP 'No valid array associated to requested xi3d output'
        ENDIF
      ENDIF

   END SUBROUTINE sto_alloc_init

   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stoalloc
