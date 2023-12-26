#include "cppdefs.h"

MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined XIOS && defined key_pisces 
   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE sms_pisces  ! PISCES variables
   USE trc         ! passive tracers common variables 
   USE sedwri
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

   !! * Substitutions
#include "ocean2pisces.h90"
#include "do_loop_substitute.h90"


CONTAINS

   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)   :: cltra
      INTEGER              :: ji, jj, jk, jn
      REAL(wp), DIMENSION(A2D(0),jpk) :: ztra
      !!---------------------------------------------------------------------
 
      IF( ln_sediment )  CALL sed_wri

      DO jn = jp_pcs0, jp_pcs1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         DO_3D( 0, 0, 0, 0, 1, jpk)
            ztra(ji,jj,jk) = MAX( 0., tr(ji,jj,jk,jn,nnew) )
         END_3D   
         CALL iom_put( cltra, ztra )
      END DO

# if defined key_trc_diaadd
      CALL iom_put( "Cflx"    , trc2d(:,:,jp_flxco2)   )
      CALL iom_put( "Oflx"    , trc2d(:,:,jp_flxo2)    )
      CALL iom_put( "Kg"      , trc2d(:,:,jp_kgco2)    )
      CALL iom_put( "Dpco2"   , trc2d(:,:,jp_dpco2)    )
      CALL iom_put( "EPC100"  , trc2d(:,:,jp_sinkco2)  )
      CALL iom_put( "Heup"    , trc2d(:,:,jp_heup)     )
      CALL iom_put( "INTNFIX" , trc2d(:,:,jp_nfix)     )
      CALL iom_put( "No3dep"  , trc2d(:,:,jp_no3dep)   ) 
#if ! defined key_pisces_light
      CALL iom_put( "EPFE100" , trc2d(:,:,jp_sinkfer)  )
      CALL iom_put( "EPSI100" , trc2d(:,:,jp_sinksil)  )
      CALL iom_put( "EPCAL100", trc2d(:,:,jp_sinkcal)  )
      CALL iom_put( "Nh4dep"  , trc2d(:,:,jp_nh4dep)   ) 
      CALL iom_put( "Sildep"  , trc2d(:,:,jp_sildep)   ) 
      CALL iom_put( "Po4dep"  , trc2d(:,:,jp_po4dep)   ) 
#endif
      !
      CALL iom_put( "PH"      , trc3d(:,:,:,jp_hi)     )
      CALL iom_put( "CO3"     , trc3d(:,:,:,jp_co3)    )
      CALL iom_put( "CO3sat"  , trc3d(:,:,:,jp_co3sat) )
      CALL iom_put( "PAR"     , trc3d(:,:,:,jp_etot)   )
      CALL iom_put( "PPPHY"   , trc3d(:,:,:,jp_pphy)   )
      CALL iom_put( "GRAZ1" , trc3d(:,:,:,jp_grapoc)  )
      CALL iom_put( "MicroZo2", trc3d(:,:,:,jp_mico2 )  )
      CALL iom_put( "Remino2" , trc3d(:,:,:,jp_remino2))
      CALL iom_put( "Nfixo2"  , trc3d(:,:,:,jp_nfixo2) )
      CALL iom_put( "Irondep" , trc3d(:,:,:,jp_irondep)  ) 
      CALL iom_put( "Ironsed" , trc3d(:,:,:,jp_ironsed) )
#if  defined key_pisces_light
      CALL iom_put( "Thetanano" , trc3d(:,:,:,jp_pnew) )
#else
      CALL iom_put( "PPNEWo2" , trc3d(:,:,:,jp_pnewo2) )
      CALL iom_put( "PPNEWN"  , trc3d(:,:,:,jp_pnew)   )
      CALL iom_put( "PPPHY2"  , trc3d(:,:,:,jp_pphy2)  )
      CALL iom_put( "PPNEWD"  , trc3d(:,:,:,jp_pnew2)  )
      CALL iom_put( "PBSi"    , trc3d(:,:,:,jp_pbsi)   )
      CALL iom_put( "PFeN"    , trc3d(:,:,:,jp_pfen)   )
      CALL iom_put( "PFeD"    , trc3d(:,:,:,jp_pfed)   )
      CALL iom_put( "PPRego2" , trc3d(:,:,:,jp_prego2) )
      CALL iom_put( "GRAZ2"   , trc3d(:,:,:,jp_grapoc2) )
      CALL iom_put( "MesoZo2" , trc3d(:,:,:,jp_meso2 )  )
      CALL iom_put( "Nitrifo2", trc3d(:,:,:,jp_nitrifo2) )
#endif
#endif
      !
   END SUBROUTINE trc_wri_pisces
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_pisces
CONTAINS
   SUBROUTINE trc_wri_pisces                     ! Empty routine  
   END SUBROUTINE trc_wri_pisces
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_pisces.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_pisces
