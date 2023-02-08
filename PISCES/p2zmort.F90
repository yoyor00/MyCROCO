#include "cppdefs.h"

MODULE p2zmort
   !!======================================================================
   !!                         ***  MODULE p2zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p2z_mort       : Compute the mortality terms for phytoplankton
   !!   p2z_mort_init  : Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p2zprod         ! Primary productivity 
   USE p2zlim          ! Phytoplankton limitation terms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_mort    
   PUBLIC   p2z_mort_init    

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   REAL(wp), PUBLIC ::   wchl     !:
   REAL(wp), PUBLIC ::   mprat    !:

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zmort.F90 10227 2018-10-25 14:42:24Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsizerat, zcompaph
      REAL(wp) ::   zprcaca
      REAL(wp) ::   ztortp , zrespp , zmortp 
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      prodcal(:,:,:) = 0.   ! calcite production variable set to zero
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zcompaph = MAX( ( trb(ji,jj,K,jpphy) - 1e-8 ), 0.e0 )
               !     When highly limited by macronutrients, very small cells 
               !     dominate the community. As a consequence, aggregation
               !     due to turbulence is negligible. Mortality is also set
               !     to 0
               zsizerat = MIN(1., MAX( 0., (quotan(ji,jj,jk) - 0.2) / 0.3) ) * trb(ji,jj,K,jpphy)
               !     Squared mortality of Phyto similar to a sedimentation term during
               !     blooms (Doney et al. 1996)
               zrespp = wchl * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * zsizerat 

               !     Phytoplankton mortality. This mortality loss is slightly
               !     increased when nutrients are limiting phytoplankton growth
               !     as observed for instance in case of iron limitation.
               ztortp = mprat * xstep * zcompaph / ( xkmort + trb(ji,jj,K,jpphy) ) * zsizerat

               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks

               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               !
               zprcaca = xfracal(ji,jj,jk) * zmortp
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) +  zmortp
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
   END SUBROUTINE p2z_mort

   SUBROUTINE p2z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!              called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp2zmort/ wchl, mprat
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p2z_mort_init : Initialization of phytoplankton mortality parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampismort in reference namelist : Pisces phytoplankton
      READ  ( numnatp_ref, namp2zmort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp2zmort in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampismort in configuration namelist : Pisces phytoplankton
      READ  ( numnatp_cfg, namp2zmort, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp2zmort in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp2zmort )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp2zmort'
         WRITE(numout,*) '      quadratic mortality of phytoplankton        wchl   =', wchl
         WRITE(numout,*) '      phytoplankton mortality rate                mprat  =', mprat
      ENDIF
      !
   END SUBROUTINE p2z_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_mort                    ! Empty routine
   END SUBROUTINE p2z_mort
#endif 

   !!======================================================================
END MODULE p2zmort
