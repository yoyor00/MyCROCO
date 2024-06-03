MODULE stomod
   !!======================================================================
   !!                       ***  MODULE stomod  ***
   !!
   !! Implement stochastic parameterizations in the model
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   sto_mod      : apply parameterizations at each model time step
   !!   sto_mod_init : initialize stochastic parameterizations
   !!----------------------------------------------------------------------

   USE stoarray
   USE stopar
   USE stobulk

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: jpstomax=100   ! maximum number of stochastic arrays

   PUBLIC sto_mod, sto_mod_init

CONTAINS

   SUBROUTINE sto_mod(kt)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod  ***
      !!
      !! Purpose : apply parameterizations at each model time step
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index


      ! Update stochastic fields
      CALL sto_par(kt)

      ! Apply every stochastic parameterization
      ! Here we just include the stochastic parameterization
      ! of the bulk formulation for the air-sea fluxes
      CALL sto_bulk(kt)

   END SUBROUTINE sto_mod


   SUBROUTINE sto_mod_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod_init  ***
      !!
      !! Purpose : initialize stochastic parameterizations
      !!
      !!----------------------------------------------------------------------
      ! Request maximum number of stochastic arrays
      CALL sto_array_request_size(jpstomax)

      ! Initialization of the various stochastic schemes
      ! (including requests for stochastic arrays using sto_array_request_new)
      ! There must be one such routine for each stochastic scheme.
      ! Here we just include the stochastic parameterization
      ! of the bulk formulation for the air-sea fluxes
      CALL sto_bulk_init

      ! Initialize stochastic arrays
      CALL sto_array_init

      ! Initialize time iteration of stochastic arrays
      CALL sto_par_init

   END SUBROUTINE sto_mod_init

   !!======================================================================
END MODULE stomod
