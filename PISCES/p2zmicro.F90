#include "cppdefs.h"

MODULE p2zmicro
   !!======================================================================
   !!                         ***  MODULE p2zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
#if defined key_pisces
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p2z_micro      : Compute the sources/sinks for microzooplankton
   !!   p2z_micro_init : Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p2zlim          ! Co-limitations
   USE p2zprod         ! production
!   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_micro         ! called in p2zbio.F90
   PUBLIC   p2z_micro_init    ! called in trcsms_pisces.F90

   !!* Substitution
#  include "top_substitute.h90"
#  include "ocean2pisces.h90"

   REAL(wp), PUBLIC ::   part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::   xprefc      !: microzoo preference for POC 
   REAL(wp), PUBLIC ::   xprefn      !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::   xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::   mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::   grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::   xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::   unass       !: Non-assimilated part of food
   REAL(wp), PUBLIC ::   sigma1      !: Fraction of microzoo excretion as DOM 
   REAL(wp), PUBLIC ::   epsher      !: growth efficiency for grazing 1 
   REAL(wp), PUBLIC ::   epshermin   !: minimum growth efficiency for grazing 1

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zmicro.F90 10374 2018-12-06 09:49:35Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   knt   ! ??? 
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze  , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zbeta
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasrat, zgrasratn
      REAL(wp) :: zgrazp, zgrazm, zgrazsd
      REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zgrazing, zfezoo
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zw3d
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zcompaz = MAX( ( trb(ji,jj,K,jpzoo) - 1.e-9 ), 0.e0 )
               zfact   = xstep * tgfunc2(ji,jj,jk) * zcompaz

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz = resrat * zfact * trb(ji,jj,K,jpzoo)   &
                  &   / ( xkmort + trb(ji,jj,K,jpzoo) )       &
                  &   + resrat * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation.
               !  ---------------------------------------------------------------
               ztortz = mzrat * 1.e6 * zfact * trb(ji,jj,K,jpzoo) * (1. - nitrfac(ji,jj,jk))

               zcompaph  = MAX( ( trb(ji,jj,K,jpphy) - xthreshphy ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,K,jppoc) - xthreshpoc ), 0.e0 )
               
               !     Microzooplankton grazing
               !     ------------------------
               zfood     = xprefn * zcompaph + xprefc * zcompapoc
               zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
               zdenom    = zfoodlim / ( xkgraz + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze    = grazrat * xstep * tgfunc2(ji,jj,jk)    &
               &           * trb(ji,jj,K,jpzoo) * (1. - nitrfac(ji,jj,jk))

               zgrazp    = zgraze  * xprefn * zcompaph  * zdenom2 
               zgrazm    = zgraze  * xprefc * zcompapoc * zdenom2 

               !
               zgraztotc = zgrazp  + zgrazm 
               zgraztotn = zgrazp * quotan(ji,jj,jk) + zgrazm 

               ! Grazing by microzooplankton
               zgrazing(ji,jj,jk) = zgraztotc

               !    Various remineralization and excretion terms
               !    --------------------------------------------
               zgrasratn = ( zgraztotn + rtrn ) / ( zgraztotc + rtrn )
               zepshert  =  MIN( 1., zgrasratn)
               zbeta     = MAX(0., (epsher - epshermin) )
               zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
               zepsherv  = zepsherf * zepshert 

               zgrarem   = zgraztotc * ( 1. - zepsherv - unass )
               zgrapoc   = zgraztotc * unass

               !  Update of the TRA arrays
               !  ------------------------
               zgrarsig  = zgrarem * sigma1
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem - zgrarsig
               !
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer
               zfezoo(ji,jj,jk)    = zgrafer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc
               prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + zgrapoc
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               zmortz = ztortz + zrespz
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zmortz + zepsherv * zgraztotc 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazp
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz - zgrazm
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortz
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazm
               !
               ! calcite production
               zprcaca = xfracal(ji,jj,jk) * zgrazp
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
            END DO
         END DO
      END DO
      !
#if defined key_iomput
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
           ALLOCATE( zw3d(PRIV_3D_BIOARRAY) )
           IF( iom_use( "GRAZ1" ) ) THEN
              zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !  Total grazing of phyto by zooplankton
              CALL iom_put( "GRAZ1", zw3d )
           ENDIF
           IF( iom_use( "FEZOO" ) ) THEN
              zw3d(:,:,:) = zfezoo(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
              CALL iom_put( "FEZOO", zw3d )
           ENDIF
           DEALLOCATE( zw3d )
         ENDIF
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_grapoc) = zgrazing(ji,jj,jk) * 1.e+3 * rfact2r * tmask(ji,jj,jk) !  grazing of phyto by microzoo
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_mico2) = zgrazing(ji,jj,jk) * ( 1.- epsher - unass ) &
                  &                      * (-o2ut) * sigma1 * 1.e+3 * rfact2r * tmask(ji,jj,jk)   ! o2 consumption by Microzoo
            END DO
         END DO
      END DO
#endif
      !
      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p2z_micro


   SUBROUTINE p2z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp2zzoo/ part, grazrat, resrat, mzrat, xprefn, xprefc, &
         &                xthreshphy,  xthreshpoc, &
         &                xthresh, xkgraz, epsher, epshermin, sigma1, unass
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p2z_micro_init : Initialization of microzooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampiszoo in reference namelist : Pisces microzooplankton
      READ  ( numnatp_ref, namp2zzoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp2zzoo in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampiszoo in configuration namelist : Pisces microzooplankton
      READ  ( numnatp_cfg, namp2zzoo, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp2zzoo in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp2zzoo )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp2zzoo'
         WRITE(numout,*) '      part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '      microzoo preference for POC                     xprefc      =', xprefc
         WRITE(numout,*) '      microzoo preference for nano                    xprefn      =', xprefn
         WRITE(numout,*) '      nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) '      poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(numout,*) '      feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(numout,*) '      exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(numout,*) '      microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(numout,*) '      maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(numout,*) '      non assimilated fraction of P by microzoo       unass       =', unass
         WRITE(numout,*) '      Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(numout,*) '      Minimum efficicency of microzoo growth          epshermin   =', epshermin
         WRITE(numout,*) '      Fraction of microzoo excretion as DOM           sigma1      =', sigma1
         WRITE(numout,*) '      half sturation constant for grazing 1           xkgraz      =', xkgraz
      ENDIF
      !
   END SUBROUTINE p2z_micro_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_micro                    ! Empty routine
   END SUBROUTINE p2z_micro
#endif 

   !!======================================================================
END MODULE p2zmicro
