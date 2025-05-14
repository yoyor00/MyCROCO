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
    INTEGER, PUBLIC, SAVE :: jstogls_s
    INTEGER, PUBLIC, SAVE :: jstogls_b
    INTEGER, PUBLIC, SAVE :: jstogls_z

    ! Parameters of stochastic fields
    ! (default values are replaced by values read in namelist)
    LOGICAL, PUBLIC, SAVE :: ln_Sprod =.TRUE.  ! perturb TKE production term
    LOGICAL, PUBLIC, SAVE :: ln_Bprod =.TRUE.  ! perturb TKE destruction term
    LOGICAL, PUBLIC, SAVE :: ln_zlevs =.FALSE. ! location uncertainty in Sprod and Bprod
    REAL(wp), SAVE :: std  = 0.1   ! standard deviation of the multiplicative noise
    REAL(wp), SAVE :: tcor = 10.0  ! time correlation (in days)
    INTEGER,  SAVE :: npasses = 50 ! number of passes of the horizontal Laplacian filter
    INTEGER,  SAVE :: arorder = 1  ! order of autoregressive process
    INTEGER,  SAVE :: nupdate = 1  ! update frequency of autoregressive process (in time steps)

    PUBLIC sto_gls, sto_gls_init

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
        IF (ln_Sprod) CALL sto_array_request_new(jstogls_s)
        IF (ln_Bprod) CALL sto_array_request_new(jstogls_b)
        IF (ln_zlevs) CALL sto_array_request_new(jstogls_z)

        ! Convert tcor parameter from days to time steps
        tcor = tcor * 86400. / ( stodt * nupdate )

        ! Set features of the requested stochastic field from parameters
        IF (ln_Sprod) THEN
          ! 1. time structure
          stofields(jstogls_s)%type_t='arn'
          stofields(jstogls_s)%corr_t=tcor
          stofields(jstogls_s)%nar_order=arorder
          stofields(jstogls_s)%nar_update=nupdate
          ! 2. space structure (with diffusive operator)
          stofields(jstogls_s)%type_xy='diffusive'
          stofields(jstogls_s)%diff_passes=npasses
          stofields(jstogls_s)%diff_type=1
          ! An alternative would be to use the kernel approach
          ! stofields(jstogls_s)%type_xy='kernel'
          ! stofields(jstogls_s)%corr_xy=10.
          ! 3. modified marginal distribution (here lognormal, with user-defined std)
          stofields(jstogls_s)%type_variate='lognormal'
          stofields(jstogls_s)%ave=1.0
          stofields(jstogls_s)%std=std
        ENDIF

        IF (ln_Bprod) THEN
          ! 1. time structure
          stofields(jstogls_b)%type_t='arn'
          stofields(jstogls_b)%corr_t=tcor
          stofields(jstogls_b)%nar_order=arorder
          stofields(jstogls_b)%nar_update=nupdate
          ! 2. space structure (with diffusive operator)
          stofields(jstogls_b)%type_xy='diffusive'
          stofields(jstogls_b)%diff_passes=npasses
          stofields(jstogls_b)%diff_type=1
          ! An alternative would be to use the kernel approach
          ! stofields(jstogls_b)%type_xy='kernel'
          ! stofields(jstogls_b)%corr_xy=10.
          ! 3. modified marginal distribution (here lognormal, with user-defined std)
          stofields(jstogls_b)%type_variate='lognormal'
          stofields(jstogls_b)%ave=1.0
          stofields(jstogls_b)%std=std
        ENDIF

        IF (ln_zlevs) THEN
          stofields(jstogls_z)%dim=3
          ! 1. time structure
          stofields(jstogls_z)%type_t='arn'
          stofields(jstogls_z)%corr_t=tcor
          stofields(jstogls_z)%nar_order=arorder
          stofields(jstogls_z)%nar_update=nupdate
          ! 2. space structure (with diffusive operator)
          stofields(jstogls_z)%type_xy='diffusive'
          stofields(jstogls_z)%diff_passes=npasses
          stofields(jstogls_z)%diff_type=1
          ! An alternative would be to use the kernel approach
          ! stofields(jstogls_z)%type_xy='kernel'
          ! stofields(jstogls_z)%corr_xy=10.
          ! 3. modified marginal distribution (here lognormal, with user-defined std)
          stofields(jstogls_z)%type_variate='lognormal'
          stofields(jstogls_z)%ave=1.0
          stofields(jstogls_z)%std=std
        ENDIF

    END SUBROUTINE sto_gls_init


    SUBROUTINE read_parameters
        !!----------------------------------------------------------------------
        !!                  ***  routine read_parameters  ***
        !!
        !! ** Purpose :   Read parameters for this stochastic module
        !!
        !!----------------------------------------------------------------------

        ! Namelist with parameters for this stochastic module
        NAMELIST/namsto_gls/ ln_Sprod, ln_Bprod, ln_zlevs, tcor, std, npasses, arorder, nupdate
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
