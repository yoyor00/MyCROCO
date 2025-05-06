MODULE stomod

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE stomod  ***
   !!
   !! Implement stochastic parameterizations in the model
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   sto_mod          : apply parameterizations at each model time step
   !!   sto_mod_init     : initialize stochastic parameterizations
   !!   sto_mod_finalize : finalize stochastic parameterizations
   !!----------------------------------------------------------------------

   USE stoarray
   USE stopar
   USE storst
   USE stoalloc
   USE stobulk
   USE stostress
   USE stogls
   USE stoics

   IMPLICIT NONE
   PRIVATE

   !INTEGER, PARAMETER :: jpstomax=100   ! maximum number of stochastic arrays

   PUBLIC sto_mod, sto_mod_init, sto_mod_finalize

CONTAINS

   SUBROUTINE sto_mod(kt)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod  ***
      !!
      !! Purpose : apply parameterizations at each model time step
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      IF (ln_stogen) THEN

        ! Update stochastic fields
        CALL sto_par(kt)

        ! Apply every stochastic parameterization
        ! Here we just include the stochastic parameterization
        ! of the bulk formulation for the air-sea fluxes
        IF (ln_stobulk) CALL sto_bulk(kt)
        IF (ln_stogls)  CALL sto_gls(kt)

      ENDIF

   END SUBROUTINE sto_mod


   SUBROUTINE sto_mod_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod_init  ***
      !!
      !! Purpose : initialize stochastic parameterizations
      !!
      !!----------------------------------------------------------------------

      ! Read general parameters from namelist
      CALL sto_param_init

      IF (ln_stogen) THEN

        ! Request maximum number of stochastic arrays
        CALL sto_array_request_size(jpstomax)

        ! Initialization of the various stochastic schemes
        ! (including requests for stochastic arrays using sto_array_request_new)
        ! There must be one such routine for each stochastic scheme.
        ! Here we just include the stochastic parameterization
        ! of the bulk formulation for the air-sea fluxes
        IF (ln_stobulk)   CALL sto_bulk_init
        IF (ln_stostress) CALL sto_stress_init
        IF (ln_stogls)    CALL sto_gls_init
        IF (ln_stoics)    CALL sto_ics_init

        ! Initialize stochastic arrays
        CALL sto_array_init

        ! Initialize required temporary arrays for output
        CALL sto_alloc_init

        ! Initialize time iteration of stochastic arrays
        CALL sto_par_init

      ENDIF

   END SUBROUTINE sto_mod_init


   SUBROUTINE sto_mod_finalize
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod_finalize  ***
      !!
      !! Purpose : finalize stochastic parameterizations
      !!
      !!----------------------------------------------------------------------

      IF (ln_stogen) THEN

        ! write final restart file
        CALL sto_rst_write

      ENDIF

   END SUBROUTINE sto_mod_finalize

   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stomod
