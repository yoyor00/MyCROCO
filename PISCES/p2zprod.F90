#include "cppdefs.h"

MODULE p2zprod
   !!======================================================================
   !!                         ***  MODULE p2zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p2z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p2z_prod_init  : Initialization of the parameters for growth
   !!   p2z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p2zlim          ! Co-limitations of differents nutrients
!  USE prtctl_trc      ! print control for debugging
!  USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_prod         ! called in p2zbio.F90
   PUBLIC   p2z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p2z_prod_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   REAL(wp), PUBLIC ::   pislopen     !:
   REAL(wp), PUBLIC ::   xadap        !:
   REAL(wp), PUBLIC ::   excretn      !:
   REAL(wp), PUBLIC ::   bresp        !:
   REAL(wp), PUBLIC ::   chlcnm       !:
   REAL(wp), PUBLIC ::   chlcmin      !:

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zprod.F90 11117 2019-06-17 08:50:02Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zpislopen
      REAL(wp) ::   zmxltst, zmxlday
      REAL(wp) ::   zrum, zcodel, zargu, zval, chlcnm_n
      REAL(wp) ::   zfact, zmsk
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zw2d, zstrn
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zw3d
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprmax, zprbio, zprorcan
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zmxl_fac, zmxl_chl
      !
      !!---------------------------------------------------------------------
      !
      !  Allocate temporary workspace
      !
      zprorcan(:,:,:) = 0. ;  zprbio  (:,:,:) = 0. 
      zmxl_fac(:,:,:) = 0. ; zmxl_chl(:,:,:) = 0. 

      ! Computation of the optimal production
      zprmax(:,:,:) = 0.8 * r1_rday * tgfunc(:,:,:)

      ! compute the day length depending on latitude and the day
      zrum = FLOAT( nday_year - 80 ) / nyear_len
      zcodel = ASIN(  SIN( zrum * rpi * 2. ) * SIN( rad * 23.5 )  )

      ! day length in hours
      zstrn(:,:) = 0.
      DO jj = JRANGE
         DO ji = IRANGE
            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
            zargu = MAX( -1., MIN(  1., zargu ) )
            zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
         END DO
      END DO

      ! Impact of the day duration and light intermittency on phytoplankton growth
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  zval = MAX( 1., zstrn(ji,jj) )
                  IF( gdept_n(ji,jj,K) <= hmld(ji,jj) ) THEN
                     zval = zval * MIN(1., heup_01(ji,jj) / ( hmld(ji,jj) + rtrn ))
                  ENDIF
                  zmxl_chl(ji,jj,jk) = zval / 24.
                  zmxl_fac(ji,jj,jk) = 1.5 * zval / ( 12. + zval )
               ENDIF
            END DO
         END DO
      END DO

      zprbio(:,:,:) = zprmax(:,:,:) * zmxl_fac(:,:,:)

      ! Maximum light intensity
      WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24.

      ! Computation of the P-I slope for nanos and diatoms
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                   ! Computation of production function for Carbon
                   !  ---------------------------------------------
                   zpislopen = pislopen * thetanano(ji,jj,jk) &
                      &      / ( zprbio(ji,jj,jk) * rday * xlimphy(ji,jj,jk) + rtrn )
                   zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) )  )
               ENDIF
            END DO
         END DO
      END DO

      !  Computation of a proxy of the N/C ratio
      !  ---------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
                zval = xnanono3(ji,jj,jk) * zprmax(ji,jj,jk) / ( zprbio(ji,jj,jk) + rtrn )
                quotan(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval )

               ! Diagnostic Chl/C ratio according to Geider et al. (1997)
               ! --------------------------------------------------------
               thetanano(ji,jj,jk) = chlcnm / ( 1.0 + pislopen * chlcnm * emoy(ji,jj,jk)   &
                &                  / ( 2.0 * zprmax(ji,jj,jk) * zmxl_fac(ji,jj,jk) * xlimphy(ji,jj,jk) * rday + rtrn ) )
               thetanano(ji,jj,jk) = MAX( chlcmin, thetanano(ji,jj,jk) )
            END DO
         END DO
      END DO

      !  Mixed-layer effect on production 
      !  Sea-ice effect on production

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO

      ! Computation of the various production terms 
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. (C)
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,K,jpphy) * rfact2
                  !
               ENDIF
            END DO
         END DO
      END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
              IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                 tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zprorcan(ji,jj,jk)
                 tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zprorcan(ji,jj,jk) * excretn
                 tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zprorcan(ji,jj,jk) 
                 !
                 tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zprorcan(ji,jj,jk) * ferat3 * texcretn
                 tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) 
                 tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zprorcan(ji,jj,jk)
              ENDIF
           END DO
        END DO
     END DO
     !
#if defined key_iomput
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "TPP" )  )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprorcan(:,:,jk) * zfact * tmask(:,:,jk)  ! primary production by nanophyto
              END DO
              CALL iom_put( "TPP"  , zw3d )
              !
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprmax(:,:,jk) * tmask(:,:,jk)   ! Maximum growth rate
              END DO
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprbio(:,:,jk) * xlimphy(:,:,jk) * tmask(:,:,jk)  ! Realized growth rate for nanophyto
              END DO
              CALL iom_put( "MuN"  , zw3d )
              !
          ENDIF
          IF( iom_use( "LNlight" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprbio (:,:,jk) / (zprmax(:,:,jk) + rtrn) * tmask(:,:,jk) ! light limitation term
              END DO
              CALL iom_put( "LNlight"  , zw3d )
              !
          ENDIF
          IF( iom_use( "INTPP" ) ) THEN   
             zw2d(:,:) = 0.
             DO jk = KRANGE
                zw2d(:,:) = zw2d(:,:) + zprorcan(:,:,jk) * e3t_n(:,:,K) * zfact * tmask(:,:,jk) ! vert. integrated pp
             ENDDO
             CALL iom_put( "INTPP" , zw2d )
          ENDIF
          !
       ENDIF
     ENDIF
#endif

#if defined key_trc_diaadd 
      !   Supplementary diagnostics
     zfact = 1.e3 * rfact2r
     DO jk = KRANGE
        DO jj = JRANGE
          DO ji = IRANGE
             zmsk = zfact * tmask(ji,jj,K)
             trc3d(ji,jj,K,jp_pphy  )  = zprorcan(ji,jj,jk) * zmsk  ! primary production by nanophyto
             trc3d(ji,jj,K,jp_pnew  )  = thetanano(ji,jj,jk) * zmsk ! new primary production by nanophyto
             trc3d(ji,jj,K,jp_pnewo2)  = ( o2ut + o2nit ) * zprorcan(ji,jj,jk) * zmsk   ! Oxygen production by the New Produc.
         END DO
        END DO
      END DO
#endif

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
     ENDIF
      !
   END SUBROUTINE p2z_prod


   SUBROUTINE p2z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp2zprod/ pislopen, bresp, excretn,  &
         &                 chlcnm, chlcmin
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'p2z_prod_init : phytoplankton growth'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp2zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp2zprod in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp2zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp2zprod in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp2zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp2zprod'
         WRITE(numout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '      Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
      ENDIF
      !
      r1_rday   = 1. / rday 
      texcretn  = 1. - excretn
      !
   END SUBROUTINE p2z_prod_init


   INTEGER FUNCTION p2z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(PRIV_3D_BIOARRAY), STAT = p2z_prod_alloc )
      !
      IF( p2z_prod_alloc /= 0 ) CALL ctl_warn( 'p2z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p2z_prod_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_prod                    ! Empty routine
   END SUBROUTINE p2z_prod
#endif


   !!======================================================================
END MODULE p2zprod
