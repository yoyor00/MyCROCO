#include "cppdefs.h"

MODULE p2zlim
   !!======================================================================
   !!                         ***  MODULE p2zlim  ***
   !! TOP :   Computes the nutrient limitation terms of phytoplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p2z_lim        :   Compute the nutrients limitation terms 
   !!   p2z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES variables
!   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p2z_lim           ! called in p4zbio.F90 
   PUBLIC p2z_lim_init      ! called in trcsms_pisces.F90 
   PUBLIC p2z_lim_alloc     ! called in trcini_pisces.F90

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concnno3    !:  NO3, PO4 half saturation   
   REAL(wp), PUBLIC ::  concbno3    !:  NO3 half saturation  for bacteria 
   REAL(wp), PUBLIC ::  concnfer    !:  Iron half saturation for nanophyto 
   REAL(wp), PUBLIC ::  concbfe     !:  Fe half saturation for bacteria 
   REAL(wp), PUBLIC ::  xsizephy    !:  Minimum size criteria for nanophyto
   REAL(wp), PUBLIC ::  xsizern     !:  Size ratio for nanophytoplankton
   REAL(wp), PUBLIC ::  xkdoc       !:  2nd half-sat. of DOC remineralization  
   REAL(wp), PUBLIC ::  caco3r      !:  mean rainratio 
   REAL(wp), PUBLIC ::  oxymin      !:  half saturation constant for anoxia

   !!* Phytoplankton limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanono3   !: Nanophyto limitation by NO3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphy    !: Nutrient limitation term of nanophytoplankton
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbac    !: Bacterial limitation term
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbacl   !: Bacterial limitation term
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: Nanophyto limitation by Iron

   LOGICAL  :: l_dia
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zlim.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_lim( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!                for the unique phytoplankton species 
      !!
      !! ** Method  : - Limitation is computed according to Monod formalism
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zconcn, zconcn2, z1_trbphy, zconc0n, zconcnf, zlim1, zlim2, zlim3
      REAL(wp) ::   ztem1, ztem2, zetot1, zetot2
      REAL(wp) ::   zferlim, zno3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  &
           & l_dia = iom_use( "LNnut" ) .OR. iom_use( "LNFe" ) .OR. iom_use( "xfracal" )

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE

               ! Tuning of the iron concentration to a minimum level that is set to the detection limit
               !-------------------------------------
               zno3    = trb(ji,jj,K,jpno3) / 40.e-6
               zferlim = MAX( 3e-11 * zno3 * zno3, 5e-12 )
               zferlim = MIN( zferlim, 7e-11 )
               trb(ji,jj,K,jpfer) = MAX( trb(ji,jj,K,jpfer), zferlim )

               ! Computation of a variable Ks for iron on diatoms taking into account
               ! that increasing biomass is made of generally bigger cells
               !------------------------------------------------
               zconcn   = MAX( 0.e0 , trb(ji,jj,K,jpphy) - xsizephy )
               zconcn2  = trb(ji,jj,K,jpphy) - zconcn
               z1_trbphy   = 1. / ( trb(ji,jj,K,jpphy) + rtrn )

               zconcnf = MAX( concnfer, ( zconcn2 * concnfer + concnfer * xsizern * zconcn ) * z1_trbphy )
               zconc0n = MAX( concnno3, ( zconcn2 * concnno3 + concnno3 * xsizern * zconcn ) * z1_trbphy )
               
               ! Michaelis-Menten Limitation term for nutrients Small bacteria
               ! -------------------------------------------------------------
               zlim1  = trb(ji,jj,K,jpno3) / ( concbno3 + trb(ji,jj,K,jpno3) )
               zlim2  = trb(ji,jj,K,jpfer) / ( concbfe + trb(ji,jj,K,jpfer) )
               zlim3  = trb(ji,jj,K,jpdoc) / ( xkdoc   + trb(ji,jj,K,jpdoc) )
               xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2 )
               xlimbac (ji,jj,jk) = xlimbacl(ji,jj,jk) * zlim3

               ! Michaelis-Menten Limitation term for nutrients Small flagellates
               ! -----------------------------------------------
               ! Limitation of nanophytoplankton growth
               xnanono3(ji,jj,jk) = trb(ji,jj,K,jpno3) / ( zconc0n + trb(ji,jj,K,jpno3) )
               xlimnfe (ji,jj,jk) = trb(ji,jj,K,jpfer) / ( zconcnf + trb(ji,jj,K,jpfer) ) 
               xlimphy (ji,jj,jk) = MIN( xlimnfe(ji,jj,jk), xnanono3(ji,jj,jk) )
               !
           END DO
         END DO
      END DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! --------------------------------------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ztem1  = MAX( 0., tsn(ji,jj,K,jp_tem) )
               ztem2  = tsn(ji,jj,K,jp_tem) - 10.
               zetot1 = MAX( 0., etot_ndcy(ji,jj,jk) - 1.) / ( 4. + etot_ndcy(ji,jj,jk) )
               zetot2 = 30. / ( 30. + etot_ndcy(ji,jj,jk) )

               xfracal(ji,jj,jk) = caco3r * xlimphy(ji,jj,jk)                  &
                  &                       * ztem1 / ( 0.1 + ztem1 )                     &
                  &                       * MAX( 1., trb(ji,jj,K,jpphy) * 1.e6 / 2. )  &
                  &                       * zetot1 * zetot2               &
                  &                       * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
                  &                       * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
               xfracal(ji,jj,jk) = MIN( 0.8 , xfracal(ji,jj,jk) )
               xfracal(ji,jj,jk) = MAX( 0.02, xfracal(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trb(ji,jj,K,jpoxy) )    &
                  &                                / ( oxymin + trb(ji,jj,K,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
               !
               ! denitrification factor computed from NO3 levels
               nitrfac2(ji,jj,jk) = MAX( 0.e0,       ( 1.E-6 - trb(ji,jj,K,jpno3) )  &
                  &                                / ( 1.E-6 + trb(ji,jj,K,jpno3) ) )
               nitrfac2(ji,jj,jk) = MIN( 1., nitrfac2(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         IF( l_dia ) THEN
            ALLOCATE( zw3d(GLOBAL_2D_ARRAY,1:jpk) )   ;   zw3d(:,:,:) = 0.
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = xlimphy(ji,jj,jk) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "LNnut", zw3d )  ! Nutrient limitation term
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = xlimnfe(ji,jj,jk) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "LNFe", zw3d )  ! Iron limitation term
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = xfracal(ji,jj,jk) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "xfracal", zw3d )  ! Calcifiers
            DEALLOCATE( zw3d )
         ENDIF
      ENDIF
      !
   END SUBROUTINE p2z_lim


   SUBROUTINE p2z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of the nutrient limitation parameters
      !!
      !! ** Method  :   Read the namp2zlim namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp2zlim
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer

      ! Namelist block
      NAMELIST/namp2zlim/concnno3, concbno3, concnfer, concbfe, &
         &               xsizephy, xsizern,  xkdoc,  &
         &               caco3r, oxymin
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p2z_lim_init : initialization of nutrient limitations'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampislim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp2zlim, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp2zlim in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampislim in configuration namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_cfg, namp2zlim, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp2zlim in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp2zlim )

      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp2zlim'
         WRITE(numout,*) '      NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '      NO3 half saturation of phyto             concnno3  = ', concnno3
         WRITE(numout,*) '      Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '      Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '      half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '      size ratio for phytoplankton             xsizern   = ', xsizern
         WRITE(numout,*) '      Minimum size criteria for phyto          xsizephy  = ', xsizephy
         WRITE(numout,*) '      mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '      halk saturation constant for anoxia      oxymin    =' , oxymin
      ENDIF
      !
      nitrfac (:,:,:) = 0.0
      nitrfac2(:,:,:) = 0.0
      !
   END SUBROUTINE p2z_lim_init


   INTEGER FUNCTION p2z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !! 
      !            Allocation of the arrays used in this module
      !!----------------------------------------------------------------------

      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xnanono3(jpi,jpj,jpk), xlimphy(jpi,jpj,jpk),       &
         &      xlimnfe (jpi,jpj,jpk), xlimbac (jpi,jpj,jpk),       &
         &      xlimbacl(jpi,jpj,jpk),                       STAT=p2z_lim_alloc )
 
      !
      IF( p2z_lim_alloc /= 0 )  CALL ctl_stop( 'p2z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p2z_lim_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_lim                   ! Empty routine
   END SUBROUTINE p2z_lim
#endif
   !!======================================================================
END MODULE p2zlim
