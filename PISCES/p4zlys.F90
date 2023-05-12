#include "cppdefs.h"

MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_lys        :   Compute the CaCO3 dissolution 
   !!   p4z_lys_init   :   Read the namelist parameters
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  Chemical model
!   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_lys         ! called in trcsms_pisces.F90
   PUBLIC   p4z_lys         ! called in trcsms_pisces.F90
   PUBLIC   p4z_lys_init    ! called in trcsms_pisces.F90

#include "ocean2pisces.h90"

   REAL(wp), PUBLIC ::   kdca   !: diss. rate constant calcite
   REAL(wp), PUBLIC ::   nca    !: order of reaction for calcite dissolution

   REAL(wp) ::   calcon = 1.03E-2   ! mean calcite concentration [Ca2+] in sea water [mole/kg solution]
   LOGICAL  :: l_dia

 
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zlys.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p2z_lys( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      !
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zdispot, zfact, zcalcon, zdepexp, zdissol
      REAL(wp) ::   zomegaca, zexcess, zexcess0, zkd, zwsbio
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zco3, zhinit, zhi, zcaco3, ztra
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  &
           & l_dia = iom_use( "PH" ) .OR. iom_use( "CO3" ) .OR. iom_use( "CO3sat" )

      zco3(:,:,:) = 0.
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zhinit(ji,jj,jk) = hi(ji,jj,jk) * 1000. / ( rhop(ji,jj,K) + rtrn )
            END DO
         END DO
      END DO
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general( zhinit, zhi )

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zco3(ji,jj,jk) = trb(ji,jj,K,jpdic) * ak13(ji,jj,jk) * ak23(ji,jj,jk) / (zhi(ji,jj,jk)**2   &
                  &             + ak13(ji,jj,jk) * zhi(ji,jj,jk) + ak13(ji,jj,jk) * ak23(ji,jj,jk) + rtrn )
               hi  (ji,jj,jk) = zhi(ji,jj,jk) * rhop(ji,jj,K) / 1000.
            END DO
         END DO
      END DO

      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE

               ! DEVIATION OF [CO3--] FROM SATURATION VALUE
               ! Salinity dependance in zomegaca and divide by rhop/1000 to have good units
               zcalcon  = calcon * ( salinprac(ji,jj,jk) / 35.0 )
               zfact    = rhop(ji,jj,K) / 1000.0
               zomegaca = ( zcalcon * zco3(ji,jj,jk) ) / ( aksp(ji,jj,jk) * zfact + rtrn )

               ! SET DEGREE OF UNDER-/SUPERSATURATION
               excess(ji,jj,jk) = 1. - zomegaca
               zexcess0 = MAX( 0., excess(ji,jj,jk) )

              IF( zomegaca < 0.8 ) THEN
                zexcess = zexcess0**nca
               ! AMOUNT CACO3 THAT RE-ENTERS SOLUTION
                zdispot = kdca * zexcess
             ELSE
                zkd = kdca * 0.2**(nca - 0.2)
                zexcess = zexcess0**0.2
                zdispot = zkd * zexcess
             ENDIF
        
             !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
             !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
             ztra(ji,jj,jk)  = zdispot / rmtss ! calcite dissolution
            END DO
         END DO
      END DO
      !
      DO jj = JRANGE
         DO ji = IRANGE
            zcaco3(ji,jj,1) = prodcal(ji,jj,1) * rfact2r / ( wsbio4(ji,jj,1) / e3t_n(ji,jj,KSURF) / rday + ztra(ji,jj,1) )
         END DO
      END DO
      DO jk = KRANGEL
         DO jj = JRANGE
            DO ji = IRANGE
                zdissol = 0.0
                zwsbio = wsbio4(ji,jj,jk) / rday
                IF( tmask(ji,jj,1) == 1. ) THEN
                   IF( ztra(ji,jj,jk) == 0.0 ) THEN
                      zcaco3(ji,jj,jk) = zcaco3(ji,jj,jk-1) + prodcal(ji,jj,jk) * rfact2r / zwsbio * e3t_n(ji,jj,K)
                   ELSE
                      zdepexp = exp( - ztra(ji,jj,jk) * e3t_n(ji,jj,K) / zwsbio )
                      zcaco3(ji,jj,jk) = prodcal(ji,jj,jk) * rfact2r / ztra(ji,jj,jk)  &
                               & * (1.0 - zdepexp ) + zcaco3(ji,jj,jk-1) * zdepexp
                      zdissol = prodcal(ji,jj,jk) * e3t_n(ji,jj,K) + prodcal(ji,jj,jk)   &
                         &      * zwsbio / ztra(ji,jj,jk) * ( zdepexp - 1.0 )   &
                         &      + zwsbio * zcaco3(ji,jj,jk-1) * ( 1.0 - zdepexp ) * rfact2
                      zdissol = zdissol / e3t_n(ji,jj,K)
                   ENDIF
                   tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2. * zdissol
                   tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +      zdissol
                ENDIF
            END DO
         END DO
      END DO
      !
      DO jj = JRANGE
         DO ji = IRANGE
            sinkcalb(ji,jj) = wsbio4(ji,jj,ikt) * zcaco3(ji,jj,ikt) * rfact2 / rday
         ENDDO
      ENDDO
      !
# if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_hi    )  = -1. * LOG10( MAX( hi(ji,jj,jk), rtrn) ) * tmask(ji,jj,jk)  ! PH
               trc3d(ji,jj,K,jp_co3   )  = zco3(ji,jj,jk)     * 1.e+3 * tmask(ji,jj,jk)  ! Ion carbonate
               trc3d(ji,jj,K,jp_co3sat)  = aksp(ji,jj,jk) / calcon * 1.e+3 * tmask(ji,jj,jk)
            ENDDO
         ENDDO
      ENDDO
# endif

      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         IF( l_dia ) THEN
            ALLOCATE( zw3d(GLOBAL_2D_ARRAY,1:jpk) )   ;   zw3d(:,:,:) = 0.
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = -1. * LOG10( MAX( hi(ji,jj,jk), rtrn) ) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "PH", zw3d(:,:,:) )
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = zco3(ji,jj,jk) * 1.e+3 * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "CO3", zw3d(:,:,:) )
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = aksp(ji,jj,jk) / calcon * 1.e+3 * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "CO3sat", zw3d(:,:,:) )
            DEALLOCATE( zw3d ) 
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc( charout, ltra='tra')
!       CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p2z_lys

   SUBROUTINE p4z_lys( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      !
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zdispot, zfact, zcalcon
      REAL(wp) ::   zomegaca, zexcess, zexcess0
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zco3, zcaldiss, zhinit, zhi, zco3sat
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  &
           & l_dia = iom_use( "PH" )     .OR. iom_use( "CO3" ) .OR.  &
           &         iom_use( "CO3sat" ) .OR. iom_use( "DCAL" )
      !
      zco3    (:,:,:) = 0.
      zcaldiss(:,:,:) = 0.
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zhinit(ji,jj,jk) = hi(ji,jj,jk) * 1000. / ( rhop(ji,jj,K) + rtrn )
            END DO
         END DO
      END DO
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general( zhinit, zhi )

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zco3(ji,jj,jk) = trb(ji,jj,K,jpdic) * ak13(ji,jj,jk) * ak23(ji,jj,jk) / (zhi(ji,jj,jk)**2   &
                  &             + ak13(ji,jj,jk) * zhi(ji,jj,jk) + ak13(ji,jj,jk) * ak23(ji,jj,jk) + rtrn )
               hi  (ji,jj,jk) = zhi(ji,jj,jk) * rhop(ji,jj,K) / 1000.
            END DO
         END DO
      END DO

      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE

               ! DEVIATION OF [CO3--] FROM SATURATION VALUE
               ! Salinity dependance in zomegaca and divide by rhop/1000 to have good units
               zcalcon  = calcon * ( salinprac(ji,jj,jk) / 35.0 )
               zfact    = rhop(ji,jj,K) / 1000.0
               zomegaca = ( zcalcon * zco3(ji,jj,jk) ) / ( aksp(ji,jj,jk) * zfact + rtrn )
               zco3sat(ji,jj,jk) = aksp(ji,jj,jk) * zfact / ( zcalcon + rtrn )

               ! SET DEGREE OF UNDER-/SUPERSATURATION
               excess(ji,jj,jk) = 1. - zomegaca
               zexcess0 = MAX( 0., excess(ji,jj,jk) )
               zexcess  = zexcess0**nca

               ! AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
               !       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
               !       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
               zdispot = kdca * zexcess * trb(ji,jj,K,jpcal)
              !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
              !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
              zcaldiss(ji,jj,jk)  = zdispot * rfact2 / rmtss ! calcite dissolution
              !
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2. * zcaldiss(ji,jj,jk)
              tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) -      zcaldiss(ji,jj,jk)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +      zcaldiss(ji,jj,jk)
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
                    zw3d(ji,jj,jk ) = -1. * LOG10( MAX( hi(ji,jj,jk), rtrn) ) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "PH", zw3d(:,:,:) )
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = zco3(ji,jj,jk) * 1.e+3 * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "CO3", zw3d(:,:,:) )
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = aksp(ji,jj,jk) / calcon * 1.e+3 * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "CO3sat", zw3d(:,:,:) )
            !
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                    zw3d(ji,jj,jk ) = zcaldiss(ji,jj,jk) * 1.e+3 * rfact2r * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put( "DCAL", zw3d(:,:,:) )
            DEALLOCATE( zw3d ) 
         ENDIF
      ENDIF
      !
# if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_hi    )  = -1. * LOG10( MAX( hi(ji,jj,jk), rtrn) ) * tmask(ji,jj,jk)  ! PH
               trc3d(ji,jj,K,jp_co3   )  = zco3(ji,jj,jk)     * 1.e+3 * tmask(ji,jj,jk)  ! Ion carbonate
               trc3d(ji,jj,K,jp_co3sat)  = zco3sat(ji,jj,jk)  * 1.e+3* tmask(ji,jj,jk)
            ENDDO
         ENDDO
      ENDDO
# endif

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc( charout, ltra='tra')
!       CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p4z_lys


   SUBROUTINE p4z_lys_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lys_init  ***
      !!
      !! ** Purpose :   Initialization of CaCO3 dissolution parameters
      !!
      !! ** Method  :   Read the nampiscal namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiscal
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/nampiscal/ kdca, nca
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_lys_init : initialization of CaCO3 dissolution'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampiscal in reference namelist : Pisces CaCO3 dissolution
      READ  ( numnatp_ref, nampiscal, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampiscal in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampiscal in configuration namelist : Pisces CaCO3 dissolution
      READ  ( numnatp_cfg, nampiscal, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampiscal in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampiscal )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : nampiscal'
         WRITE(numout,*) '      diss. rate constant calcite (per month)        kdca =', kdca
         WRITE(numout,*) '      order of reaction for calcite dissolution      nca  =', nca
      ENDIF
      !
   END SUBROUTINE p4z_lys_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_lys: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_lys
#endif 

   !!======================================================================
END MODULE p4zlys
