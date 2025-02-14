MODULE stostress

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE stostress  ***
   !!
   !! Purpose : Stochastic parameterization of the surface stress
   !!           (outside from bulk formulation, see stobulk instead)
   !!======================================================================
   USE stoexternal , only : wp, lwm, lwp, numnam_ref, numnam_cfg, numond, ctl_nam, &
                          & jpi, jpj, stodt
   USE stoarray

   IMPLICIT NONE
   PRIVATE

   ! Index of stochastic field used for the drag coefficient
   INTEGER, SAVE :: jstostress

   ! Parameters of stochastic fields
   ! (default values are replaced by values read in namelist)
   REAL(wp), SAVE :: std  = 0.01  ! standard deviation of the multiplicative noise
   REAL(wp), SAVE :: tcor = 0.0   ! time correlation (in days)
   INTEGER,  SAVE :: npasses = 20 ! number of passes of the horizontal Laplacian filter
   INTEGER,  SAVE :: arorder = 1  ! order of autoregressive process
   INTEGER,  SAVE :: nupdate = 1  ! update frequency of autoregressive process (in time steps)

   PUBLIC sto_stress_init, sto_stress

CONTAINS

   SUBROUTINE sto_stress_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_stress_init  ***
      !!
      !! This routine is called at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------

      ! Read namelist block corresponding to this stochastic scheme
      CALL read_parameters

      ! Request index for a new stochastic array
      CALL sto_array_request_new(jstostress)

      ! Convert tcor parameter from days to time steps
      tcor = tcor * 86400. / ( stodt * nupdate )

      ! Set features of the requested stochastic field from parameters
      ! 1. time structure
      stofields(jstostress)%type_t='arn'
      stofields(jstostress)%corr_t=tcor
      stofields(jstostress)%nar_order=arorder
      stofields(jstostress)%nar_update=nupdate
      ! 2. space structure (with diffusive operator)
      stofields(jstostress)%type_xy='diffusive'
      stofields(jstostress)%diff_passes=npasses
      stofields(jstostress)%diff_type=1

   END SUBROUTINE sto_stress_init


   SUBROUTINE sto_stress ( suvstr, struct_fcn )
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_stress  ***
      !!
      !! This routine implements perturbation to wind stress.
      !! Only applied to zonal component for now, 
      !! i.e. for application to the BASIN test case
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1:jpi,1:jpj), INTENT(inout) :: suvstr
      REAL(wp), DIMENSION(1:jpi,1:jpj), INTENT(inout) :: struct_fcn

      struct_fcn(:,:) = struct_fcn(:,:) * stofields(jstostress)%sto2d(:,:)
      suvstr(:,:) = suvstr(:,:) * (1 + struct_fcn(:,:))

   END SUBROUTINE sto_stress

   SUBROUTINE read_parameters
      !!----------------------------------------------------------------------
      !!                  ***  routine read_parameters  ***
      !!
      !! ** Purpose :   Read parameters for this stochastic module
      !!
      !!----------------------------------------------------------------------

      ! Namelist with parameters for this stochastic module
      NAMELIST/namsto_bulk/ tcor, std, npasses, arorder, nupdate
      !!----------------------------------------------------------------------
      INTEGER  ::   ios                            ! Local integer output status for namelist read

      ! Read namsto_bulk namelist
      REWIND( numnam_ref )
      READ  ( numnam_ref, namsto_bulk, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsto_bulk in reference namelist', lwp )

      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namsto_bulk, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsto_bulk in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsto_bulk )

   END SUBROUTINE read_parameters



   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stostress
