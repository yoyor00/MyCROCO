MODULE stogls

#include "cppdefs.h"
# if defined STOGEN
 
    !!======================================================================
    !!                       ***  MODULE stogls  ***
    !!
    !! Purpose : Stochastic parameterization of the gls mixing
    !!======================================================================
    USE stoexternal , only : wp, lwm, lwp, numnam_ref, numond, ctl_nam, &
                          & jpi, jpj, stodt
    USE stoarray

    IMPLICIT NONE
    PRIVATE

    ! Index of stochastic field used for the production/destruction terms
    INTEGER, SAVE :: jstogls_p
    INTEGER, SAVE :: jstogls_b

    ! Parameters of stochastic fields
    ! (default values are replaced by values read in namelist)
    REAL(wp), SAVE :: std  = 0.1   ! standard deviation of the multiplicative noise
    REAL(wp), SAVE :: tcor = 10.0   ! time correlation (in days)
    INTEGER,  SAVE :: npasses = 50 ! number of passes of the horizontal Laplacian filter
    INTEGER,  SAVE :: arorder = 1  ! order of autoregressive process
    INTEGER,  SAVE :: nupdate = 1  ! update frequency of autoregressive process (in time steps)

    PUBLIC sto_gls, sto_gls_init, sto_gls_pb

CONTAINS

    SUBROUTINE sto_gls(kt)
        !!----------------------------------------------------------------------
        !!
        !!                     ***  ROUTINE sto_gls  ***
        !!
        !! This routine is called at every time step
        !! to make appropriate use of the stochastic fields
        !!
        !!---------------------------------------------------------------------- 
        INTEGER, INTENT( in ) :: kt

    END SUBROUTINE sto_gls


    SUBROUTINE sto_gls_init
        !!----------------------------------------------------------------------
        !!
        !!                     ***  ROUTINE sto_gls_init  ***
        !!
        !! This routine is called at initialization time
        !! to request stochastic field with appropriate features
        !!
        !!----------------------------------------------------------------------

        ! Read namelist block corresponding to this stochastic scheme
        CALL read_parameters

        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstogls_p)
        CALL sto_array_request_new(jstogls_b)

        ! Convert tcor parameter from days to time steps
        tcor = tcor * 86400. / ( stodt * nupdate )

        ! Set features of the requested stochastic field from parameters
        ! 1. time structure
        stofields(jstogls_p)%type_t='arn'
        stofields(jstogls_b)%type_t='arn'
        stofields(jstogls_p)%corr_t=tcor
        stofields(jstogls_b)%corr_t=tcor
        stofields(jstogls_p)%nar_order=arorder
        stofields(jstogls_b)%nar_order=arorder
        stofields(jstogls_p)%nar_update=nupdate
        stofields(jstogls_b)%nar_update=nupdate
        ! 2. space structure (with diffusive operator)
        stofields(jstogls_p)%type_xy='diffusive'
        stofields(jstogls_b)%type_xy='diffusive'
        stofields(jstogls_p)%diff_passes=npasses
        stofields(jstogls_b)%diff_passes=npasses
        stofields(jstogls_p)%diff_type=1
        stofields(jstogls_b)%diff_type=1
        ! An alternative would be to use the kernel approach
        ! stofields(jstogls_pb)%type_xy='kernel'
        ! stofields(jstogls_pb)%corr_xy=10.
        ! 3. modified marginal distribution (here lognormal, with user-defined std)
        stofields(jstogls_p)%type_variate='lognormal'
        stofields(jstogls_b)%type_variate='lognormal'
        stofields(jstogls_p)%ave=1.0
        stofields(jstogls_b)%ave=1.0
        stofields(jstogls_p)%std=std
        stofields(jstogls_b)%std=std

    END SUBROUTINE sto_gls_init


    SUBROUTINE sto_gls_pb( Sprod , Bprod )
        !!----------------------------------------------------------------------
        !!
        !!                     ***  ROUTINE sto_gls_PB  ***
        !!
        !! This routine implements perturbation of the production/destruction terms
        !!
        !!----------------------------------------------------------------------
        REAL(wp), DIMENSION(1:jpi,1:jpj), INTENT(inout) :: Sprod
        REAL(wp), DIMENSION(1:jpi,1:jpj), INTENT(inout) :: Bprod

        Sprod(:,:) = Sprod(:,:) * stofields(jstogls_p)%sto2d(:,:)
        Bprod(:,:) = Bprod(:,:) * stofields(jstogls_b)%sto2d(:,:)

    END SUBROUTINE sto_gls_pb


    SUBROUTINE read_parameters
        !!----------------------------------------------------------------------
        !!                  ***  routine read_parameters  ***
        !!
        !! ** Purpose :   Read parameters for this stochastic module
        !!
        !!----------------------------------------------------------------------

        ! Namelist with parameters for this stochastic module
        NAMELIST/namsto_gls/ tcor, std, npasses, arorder, nupdate
        !!----------------------------------------------------------------------
        INTEGER  ::   ios                            ! Local integer output status for namelist read

        ! Read namsto_gls namelist
        REWIND( numnam_ref )
        READ  ( numnam_ref, namsto_gls, IOSTAT = ios, ERR = 901)
901     IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsto_gls in reference namelist', lwp )

    END SUBROUTINE read_parameters

    !!======================================================================
    
#endif /* if defined STOGEN */

END MODULE stogls
