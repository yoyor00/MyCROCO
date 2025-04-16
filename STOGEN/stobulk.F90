MODULE stobulk

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE stobulk  ***
   !!
   !! Purpose : Stochastic parameterization of the bulk formulation
   !!           for the air-sea fluxes
   !!======================================================================
   USE stoexternal , only : wp, lwm, lwp, numnam_ref, numond, ctl_nam, &
                          & jpi, jpj, stodt, ishift_priv, jshift_priv
   USE stoarray

   IMPLICIT NONE
   PRIVATE

   ! Index of stochastic field used for the drag coefficient
   INTEGER, SAVE :: jstobulk_cd

   ! Parameters of stochastic fields
   ! (default values are replaced by values read in namelist)
   REAL(wp), SAVE :: std  = 0.1   ! standard deviation of the multiplicative noise
   REAL(wp), SAVE :: tcor = 10.0   ! time correlation (in days)
   INTEGER,  SAVE :: npasses = 50 ! number of passes of the horizontal Laplacian filter
   INTEGER,  SAVE :: arorder = 1  ! order of autoregressive process
   INTEGER,  SAVE :: nupdate = 1  ! update frequency of autoregressive process (in time steps)

   PUBLIC sto_bulk, sto_bulk_init, sto_bulk_perturb

CONTAINS

   SUBROUTINE sto_bulk(kt)
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_bulk  ***
      !!
      !! This routine is called at every time step
      !! to make appropriate use of the stochastic fields
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      ! Use arrays (for instance with index jstobulk_cd) as:
      ! stofields(jstobulk_cd)%sto2d(:,:)

      ! Here the effet is directly include in bulk_flux.F
      ! using stofields(jstobulk_cd)%sto2d from stoarray

   END SUBROUTINE sto_bulk


   SUBROUTINE sto_bulk_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_bulk_init  ***
      !!
      !! This routine is calle at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------

      ! Read namelist block corresponding to this stochastic scheme
      CALL read_parameters

      ! Request index for a new stochastic array
      CALL sto_array_request_new(jstobulk_cd)

      ! Convert tcor parameter from days to time steps
      tcor = tcor * 86400. / ( stodt * nupdate )

      ! Set features of the requested stochastic field from parameters
      ! 1. time structure
      stofields(jstobulk_cd)%type_t='arn'
      stofields(jstobulk_cd)%corr_t=tcor
      stofields(jstobulk_cd)%nar_order=arorder
      stofields(jstobulk_cd)%nar_update=nupdate
      ! 2. space structure (with diffusive operator)
      stofields(jstobulk_cd)%type_xy='diffusive'
      stofields(jstobulk_cd)%diff_passes=npasses
      stofields(jstobulk_cd)%diff_type=1
      ! An alternative would be to use the kernel approach
      ! stofields(jstobulk_cd)%type_xy='kernel'
      ! stofields(jstobulk_cd)%corr_xy=10.
      ! 3. modified marginal distribution (here lognormal, with user-defined std)
      stofields(jstobulk_cd)%type_variate='lognormal'
      stofields(jstobulk_cd)%ave=1.0
      stofields(jstobulk_cd)%std=std

   END SUBROUTINE sto_bulk_init


   SUBROUTINE sto_bulk_perturb ( array, array_type )
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_bulk_perturb ***
      !!
      !! This routine implements perturbation of bulk coeffcients
      !! (currently only for Cd but other cases can be added)
      !! 
      !! array : array with bulk coefficient to perturb
      !! array_type : type of array (default=GLOBAL_2D_ARRAY),
      !!              option='private', for PRIVATE_2D_SCRATCH_ARRAY
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: array
      CHARACTER(len=*), INTENT(in), OPTIONAL :: array_type

      INTEGER :: array_type_int, jpiarray, jpjarray

      jpiarray=SIZE(array,1)
      jpjarray=SIZE(array,2)

      array_type_int = 0
      IF (PRESENT(array_type)) THEN
        IF (array_type(1:4) == 'priv') array_type_int = 1
      ENDIF

      IF (array_type_int == 1 ) THEN
        array(:,:) = array(:,:) * stofields(jstobulk_cd)%sto2d(ishift_priv:,jshift_priv:)
      ELSE
        IF (jpiarray.NE.jpi) STOP 'inconsistent array size in sto_bulk'
        IF (jpjarray.NE.jpj) STOP 'inconsistent array size in sto_bulk'
        array(:,:) = array(:,:) * stofields(jstobulk_cd)%sto2d(:,:)
      ENDIF

   END SUBROUTINE sto_bulk_perturb


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

   END SUBROUTINE read_parameters

   !!======================================================================

#endif /* if defined STOGEN */

END MODULE stobulk
