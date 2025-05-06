MODULE stoics

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE stoics  ***
   !!
   !! Purpose : Stochastic parameterization of initial conditions
   !!======================================================================
   USE stoexternal , only : wp, lwm, lwp, numnam_ref, numond, ctl_nam, &
                          & jpi, jpj, jpk
   USE stoarray

   IMPLICIT NONE
   PRIVATE

   ! Index of stochastic field used for the drag coefficient
   INTEGER, SAVE :: jstoics

   ! Parameters of stochastic fields
   ! (default values are replaced by values read in namelist)
   REAL(wp), SAVE :: std  = 0.0001   ! standard deviation of the multiplicative noise

   PUBLIC sto_ics_init, sto_ics

CONTAINS

   SUBROUTINE sto_ics_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_ics_init  ***
      !!
      !! This routine is called at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------

      ! Request index for a new stochastic array
      CALL sto_array_request_new(jstoics)

      ! Set features of the requested stochastic field from parameters
      ! 1. time structure
      stofields(jstoics)%type_t='white'
      stofields(jstoics)%type_z='white'
      stofields(jstoics)%std=std

   END SUBROUTINE sto_ics_init


   SUBROUTINE sto_ics ( ic, stoxi ) 
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_ics  ***
      !!
      !! This routine implements perturbation initial conditions.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1:jpi,1:jpj,1:jpk), INTENT(inout) :: ic
      REAL(wp), DIMENSION(1:jpi,1:jpj,1:jpk), INTENT(inout) :: stoxi

      stoxi(:,:,:) = stofields(jstoics)%sto3d(:,:,:)
      ic(:,:,:) = ic(:,:,:) * (1 + stoxi(:,:,:))

   END SUBROUTINE sto_ics

   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stoics
