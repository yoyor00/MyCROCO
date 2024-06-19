MODULE stobulk
   !!======================================================================
   !!                       ***  MODULE stobulk  ***
   !!
   !! Purpose : Stochastic parameterization of the bulk formulation
   !!           for the air-sea fluxes
   !!======================================================================
   USE stoexternal
   USE stoarray

   IMPLICIT NONE
   PRIVATE

   INTEGER, PUBLIC :: jstobulk_cd  ! index of stochastic field used for the drag coefficient

   PUBLIC sto_bulk, sto_bulk_init, sto_bulk_cd

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
      ! -> get parameters

      ! Request index for a new stochastic array
      CALL sto_array_request_new(jstobulk_cd)

      ! Set features of the requested stochastic field from parameters
      ! 1. time structure
      stofields(jstobulk_cd)%type_t='arn'
      stofields(jstobulk_cd)%corr_t=5.0
      stofields(jstobulk_cd)%nar_order=2
      stofields(jstobulk_cd)%nar_update=5
      ! 2. space structure (with diffusive operator)
      stofields(jstobulk_cd)%type_xy='diffusive'
      stofields(jstobulk_cd)%diff_passes=50
      stofields(jstobulk_cd)%diff_type=1  ! option 1 would require the mask
      ! 3. modified marginal distribution (here lognormal, with 30% std)
      stofields(jstobulk_cd)%type_variate='lognormal'
      stofields(jstobulk_cd)%ave=1.0
      stofields(jstobulk_cd)%std=0.3

   END SUBROUTINE sto_bulk_init


   SUBROUTINE sto_bulk_cd ( suvstr )
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_bulk_cd  ***
      !!
      !! This routine implements perturbation of the cd coeffcient
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1:jpi,1:jpj), INTENT(inout) :: suvstr

      print *, 'stogen bounds x:',narea-1,lbound(suvstr,1),ubound(suvstr,1)
      print *, 'stogen bounds y:',narea-1,lbound(suvstr,2),ubound(suvstr,2)
      print *, 'stogen first:',narea-1,suvstr(lbound(suvstr,1),lbound(suvstr,2))

      suvstr(:,:) = suvstr(:,:) * stofields(jstobulk_cd)%sto2d(:,:)

      print *, 'stogen after:',narea-1,suvstr(lbound(suvstr,1),lbound(suvstr,2))

   END SUBROUTINE sto_bulk_cd

   !!======================================================================
END MODULE stobulk
