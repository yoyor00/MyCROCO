MODULE stoics

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE stoics  ***
   !!
   !! Purpose : Stochastic parameterization of initial conditions
   !!======================================================================
   USE stoexternal , only : wp, lwp, numnam_ref, numond, ctl_nam, &
                          & jpi, jpj, jpk
   USE stoarray
   

   IMPLICIT NONE
   PRIVATE

   ! Index of stochastic field used for the drag coefficient
   INTEGER, PUBLIC, SAVE :: jstoics

   ! Parameters of stochastic fields
   ! (default values are replaced by values read in namelist)
   REAL(wp), SAVE  :: std  = 0.0001              ! standard deviation of the multiplicative noise
   LOGICAL, PUBLIC :: ln_strict_stoics = .TRUE. ! TRUE: apply to cold start only (i.e. ntstart.eq.1). FALSE: apply to restart as well.
 
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

      ! Read namelist block corresponding to this stochastic scheme
      CALL read_parameters

      ! Request index for a new stochastic array
      CALL sto_array_request_new(jstoics)

      ! Set features of the requested stochastic field from parameters
      ! 1. time structure
      stofields(jstoics)%type_t='constant'
      ! 2. space structure (horizontal)
      ! stofields(jstoics)%type_xy='white' !default, see stoarray.F90
      ! 2.1 space structure (vertical)
      stofields(jstoics)%type_z='white'
      ! 3. distribution parameters (std, marginal, ...)
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

   SUBROUTINE read_parameters
      !!----------------------------------------------------------------------
      !!                  ***  routine read_parameters  ***
      !!
      !! ** Purpose :   Read parameters for this stochastic module
      !!
      !!----------------------------------------------------------------------

      ! Namelist with parameters for this stochastic module
      NAMELIST/namsto_ics/ ln_strict_stoics 
      !!----------------------------------------------------------------------
      INTEGER  ::   ios                            ! Local integer output status for namelist read

      ! Read namsto_ics namelist
      REWIND( numnam_ref )
      READ  ( numnam_ref, namsto_ics, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsto_ics in reference namelist', lwp )

   END SUBROUTINE read_parameters


   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stoics
