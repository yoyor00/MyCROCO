
MODULE phytoplankton

#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"
    
    USE module_BIOLink
    USE comBIOLink
    USE COMBLOOM

    USE comBIOLink_helping
    USE shared !, only: effetchaleur, effetturbidite, epn, l_BLOOM_opt2



      IMPLICIT NONE


      PUBLIC :: phytoplankton_dynamics, storage_diag_phytoplankton, get_Tchl

 
      REAL(KIND=rsh) :: effetazotediat
      REAL(KIND=rsh) :: effetazotedino, effetphosphoredino, effetazotenano, effetphosphorenano
      REAL(KIND=rsh) :: effetlumierediat,effetlumieredino,effetlumierenano
      REAL(KIND=rsh) :: effetsilice,effetphosphorediat, effetselnutdiat
      REAL(KIND=rsh) :: diatmorteau, dinomorteau, nanomorteau


      contains

          SUBROUTINE phytoplankton_dynamics(i,j,k,c,dc)

              INTEGER, intent(in)                 :: i,j,k
              REAL(KIND=rsh),DIMENSION(nv_state), intent(in)  :: c
              REAL(KIND=rsh),DIMENSION(nv_state), intent(inout)  :: dc

              REAL(KIND=rsh) :: fluxrelatifsurf,fluxrelatiffond

              REAL(KIND=rsh) :: fractionno3diat, fractionno3dino,fractionno3nano
              REAL(KIND=rsh) :: fractionnh4diat, fractionnh4dino,fractionnh4nano
              !REAL(KIND=rsh) :: effetsilice,effetphosphorediat

              REAL(KIND=rsh) :: rationdiat,excretdiat
              REAL(KIND=rsh) :: effetnitdiat,effetamdiat,diat_iksmith
              REAL(KIND=rsh) :: pslimdiat

              REAL(KIND=rsh) :: rationdino,excretdino
              REAL(KIND=rsh) :: effetnitdino,effetamdino,effetselnutdino
              REAL(KIND=rsh) :: pslimdino,dino_iksmith

              REAL(KIND=rsh) :: rationnano,excretnano
              REAL(KIND=rsh) :: effetnitnano,effetamnano,effetselnutnano
              REAL(KIND=rsh) :: pslimnano,nano_iksmith
              REAL(KIND=rsh) :: oresphyto

              

              IF(l_physadaptation) THEN
                  diat_iksmith=p_diat_iksmith*(1.0_rsh+0.5_rsh*SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLINK-100))) &
                      *max(0.5_rsh,(1.0_rsh-effetturbidite))
                  dino_iksmith=p_dino_iksmith*max(0.5_rsh,(1.0_rsh-effetturbidite))
                  nano_iksmith=p_nano_iksmith*max(0.5_rsh,(1.0_rsh-effetturbidite))


              ELSE
                  diat_iksmith=p_diat_iksmith
                  dino_iksmith=p_dino_iksmith
                  nano_iksmith=p_nano_iksmith
              ENDIF

              fluxrelatifsurf=PAR_top_layer(k,i,j)/(diat_iksmith+0.0000000001_rsh)
              fluxrelatiffond=PAR_top_layer(k-1,i,j)/(diat_iksmith+0.0000000001_rsh)
              effetlumierediat=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*    &
                  log((fluxrelatifsurf+sqrt(1.0_rsh+fluxrelatifsurf* &
                  fluxrelatifsurf))/(fluxrelatiffond +               &
                  sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))

              fluxrelatifsurf=PAR_top_layer(k,i,j)/(dino_iksmith+0.0000000001_rsh)
              fluxrelatiffond=PAR_top_layer(k-1,i,j)/(dino_iksmith+0.0000000001_rsh)
              effetlumieredino=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*    &
                  log((fluxrelatifsurf+sqrt(1.0_rsh+fluxrelatifsurf* &
                  fluxrelatifsurf))/(fluxrelatiffond +               &
                  Sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))

              fluxrelatifsurf=PAR_top_layer(k,i,j)/(nano_iksmith+0.0000000001_rsh)
              fluxrelatiffond=PAR_top_layer(k-1,i,j)/(nano_iksmith+0.0000000001_rsh)
              effetlumierenano=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*    &
                  log((fluxrelatifsurf+sqrt(1.0_rsh+fluxrelatifsurf* &
                  fluxrelatifsurf))/(fluxrelatiffond +               &
                  sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))


              effetnitdiat=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_diat_kNO3+ &
                  (c(iv_nutr_NH4)*p_diat_kNO3/p_diat_kNH4))
              effetnitdino=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_dino_kNO3+ &
                  (c(iv_nutr_NH4)*p_dino_kNO3/p_dino_kNH4))
              effetnitnano=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_nano_kNO3+ &
                  (c(iv_nutr_NH4)*p_nano_kNO3/p_nano_kNH4))


              effetamdiat=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_diat_kNH4+  &
                  (c(iv_nutr_NO3)*p_diat_kNH4/p_diat_kNO3))
              effetamdino=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_dino_kNH4+  &
                  (c(iv_nutr_NO3)*p_dino_kNH4/p_dino_kNO3))
              effetamnano=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_nano_kNH4+  &
                  (c(iv_nutr_NO3)*p_nano_kNH4/p_nano_kNO3))

              effetazotediat=effetnitdiat+effetamdiat
              effetazotedino=effetnitdino+effetamdino
              effetazotenano=effetnitnano+effetamnano

              effetsilice=c(iv_nutr_SiOH)/(c(iv_nutr_SiOH)+p_diat_kSi)

              effetphosphorediat=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_diat_kPO4)
              effetphosphoredino=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_dino_kPO4)
              effetphosphorenano=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_nano_kPO4)

              effetselnutdiat=min(effetazotediat,effetsilice,effetphosphorediat)
              effetselnutdino=min(effetazotedino,effetphosphoredino)
              effetselnutnano=min(effetazotenano,effetphosphorenano)

              pslimdiat=min(effetlumierediat,effetselnutdiat)
              pslimdino=min(effetlumieredino,effetselnutdino)
              pslimnano=min(effetlumierenano,effetselnutnano)

              rationdiat=p_diat_mumax*effetchaleur*pslimdiat
              rationdino=p_dino_mumax*effetchaleur*pslimdino * max(0.0_rsh,(1.0_rsh-ECT_kij/p_dino_thhold_ect))
              IF(l_BLOOM_opt2) rationdino=p_dino_mumax*effetchaleur*pslimdino
              rationnano=p_nano_mumax*effetchaleur*pslimnano

              excretdiat=0.0_rsh
              excretdino=0.0_rsh
              excretnano=0.0_rsh

              diatmorteau=p_diat_mort*effetchaleur
              dinomorteau=p_dino_mort*effetchaleur
              nanomorteau=p_nano_mort*effetchaleur

              IF(c(iv_phyto_diat_N).le.p_diat_thhold_mort) diatmorteau=0.0_rsh
              IF(c(iv_phyto_dino_N).le.p_dino_thhold_mort) dinomorteau=0.0_rsh
              IF (c(iv_phyto_nano_N).lt.p_nano_thhold_mort) nanomorteau=0.0_rsh

              IF (effetazotediat.gt.0.000001_rsh) THEN
                  fractionno3diat=effetnitdiat/effetazotediat
                  fractionnh4diat=effetamdiat/effetazotediat
              ELSE
                  fractionno3diat=0.0_rsh
                  fractionnh4diat=0.0_rsh
              ENDIF
              IF (effetazotedino.gt.0.000001_rsh) THEN
                  fractionno3dino=effetnitdino/effetazotedino
                  fractionnh4dino=effetamdino/effetazotedino
              ELSE
                  fractionno3dino=0.0_rsh
                  fractionnh4dino=0.0_rsh
              ENDIF
              IF (effetazotenano.gt.0.000001_rsh) THEN
                  fractionno3nano=effetnitnano/effetazotenano
                  fractionnh4nano=effetamnano/effetazotenano
              ELSE
                  fractionno3nano=0.0_rsh
                  fractionnh4nano=0.0_rsh
              ENDIF

              oresphyto=0.0_rsh
              IF (c(iv_oxygen).gt.0.2_rsh) then
                  oresphyto = p_phyto_resp*( (1.0_rsh- effetlumierediat)*c(iv_phyto_diat_N) &
                      + (1.0_rsh-effetlumieredino)*c(iv_phyto_dino_N)  &
                      + (1.0_rsh-effetlumierenano)*c(iv_phyto_nano_N) )
              ENDIF

              !++++++++++++++++++++++++++++++++++++++++++++++++++
              ! EQUATIONS D EVOLUTION  ! concentration en jour-1
              !++++++++++++++++++++++++++++++++++++++++++++++++++

              ! Evolution des biomasses:
              !------------------------
              dc(iv_phyto_nano_N) = dc(iv_phyto_nano_N)    &
                  + c(iv_phyto_nano_N)*(rationnano-excretnano-nanomorteau) 

              dc(iv_phyto_diat_N) = dc(iv_phyto_diat_N)    &
                  + c(iv_phyto_diat_N)*(rationdiat-diatmorteau-excretdiat)

              dc(iv_phyto_dino_N) = dc(iv_phyto_dino_N)    &
                  + (rationdino-dinomorteau-excretdino)*c(iv_phyto_dino_N)

              dc(iv_detr_N)  = dc(iv_detr_N)    &
                  + (diatmorteau+excretdiat)*c(iv_phyto_diat_N)

              dc(iv_detr_Si) = dc(iv_detr_Si)  &
                  + rapsiaz*(diatmorteau+excretdiat)*c(iv_phyto_diat_N) 

              dc(iv_detr_P)  = dc(iv_detr_P)   &
                  + (diatmorteau+excretdiat)*c(iv_phyto_diat_N)*rappaz


!              IF(l_BLOOM_opt2) THEN
#if defined key_BLOOM_opt2

                  dc(iv_diss_N)    = dc(iv_diss_N)             &
                      + dinomorteau *c(iv_phyto_dino_N)        &
                      + nanomorteau *c(iv_phyto_nano_N)

                  dc(iv_diss_P)  = dc(iv_diss_P)               &
                      + nanomorteau*c(iv_phyto_nano_N)*rappaz  &
                      + dinomorteau*c(iv_phyto_dino_N)*rappaz

                  ! Production azotée
                  dc(iv_phyto_nano_pp) = dc(iv_phyto_nano_pp)  &
                      + 14.e-3*rationnano*c(iv_phyto_nano_N)*epn

                  dc(iv_phyto_diat_pp) = dc(iv_phyto_diat_pp)   &
                      + 14.e-3*rationdiat*c(iv_phyto_diat_N)*epn

                  dc(iv_phyto_dino_pp) = dc(iv_phyto_dino_pp)   &
                      + 14.e-3*rationdino*c(iv_phyto_dino_N)*epn
                  

              !ELSE
#else

                  dc(iv_detr_N) = dc(iv_detr_N)                     &
                      + (dinomorteau+excretdino)*c(iv_phyto_dino_N) &
                      + (nanomorteau+excretnano)*c(iv_phyto_nano_N) 

                  dc(iv_detr_P)  = dc(iv_detr_P)   &
                      + rappaz*((dinomorteau+excretdino)*c(iv_phyto_dino_N) &
                      + (nanomorteau+excretnano)*c(iv_phyto_nano_N))



                  ! Production azotée
                  dc(iv_phyto_nano_pp) = dc(iv_phyto_nano_pp)          &
                      + 12.e-3*p_phyto_CNratio*rationnano*c(iv_phyto_nano_N)*epn

                  dc(iv_phyto_diat_pp) = dc(iv_phyto_diat_pp)          &
                      + 12.e-3*p_phyto_CNratio*rationdiat*c(iv_phyto_diat_N)*epn

                  dc(iv_phyto_dino_pp) = dc(iv_phyto_dino_pp)          &
                      + 12.e-3*p_phyto_CNratio*rationdino*c(iv_phyto_dino_N)*epn

              !ENDIF
#endif         

              ! Evolution de l ammonium
              ! -----------------------

              dc(iv_nutr_NH4)= dc(iv_nutr_NH4)                     &
                  - fractionnh4diat*rationdiat*c(iv_phyto_diat_N)  &
                  -fractionnh4dino*rationdino*c(iv_phyto_dino_N)   &
                  -fractionnh4nano*rationnano*c(iv_phyto_nano_N)          

              ! Evolution du nitrate
              ! --------------------

              dc(iv_nutr_NO3)= dc(iv_nutr_NO3)                    &
                  -fractionno3diat*rationdiat*c(iv_phyto_diat_N) &
                  -fractionno3dino*rationdino*c(iv_phyto_dino_N) &
                  -fractionno3nano*rationnano*c(iv_phyto_nano_N)

              ! Evolution de la silice dissoute
              ! -------------------------------

              dc(iv_nutr_SiOH)= dc(iv_nutr_SiOH)  &
                  -rapsiaz*rationdiat*c(iv_phyto_diat_N)

              ! Evolution du phosphate dissous
              !-------------------------------
              dc(iv_nutr_PO4) = dc(iv_nutr_PO4)           &
                  - rappaz*(rationdiat*c(iv_phyto_diat_N) &
                  + rationdino*c(iv_phyto_dino_N)         &
                  + rationnano*c(iv_phyto_nano_N))

              ! Evolution de l'oxygène dissous
              !-------------------------------
              dc(iv_oxygen) = dc(iv_oxygen)         &
                  + ratio_mgO_to_mumolN_anabol *    &
                  ( rationdiat*c(iv_phyto_diat_N)   &
                  + rationdino*c(iv_phyto_dino_N)   &
                  + rationnano*c(iv_phyto_nano_N) ) &
                  - ratio_mgO_to_mumolN_catabol * effetchaleur * oresphyto


            ! sinking velocity
!#if ! defined key_BLOOM_opt2

!                  CALL diatom_sinking_velocity(i,j,k)
!#endif


          END SUBROUTINE phytoplankton_dynamics


          function get_Tchl(i,j,k,c) result(Tchl)
            implicit none
            INTEGER, intent(in) :: i,j,k
            REAL(KIND=rsh), DIMENSION(nv_state), intent(in) :: c
            REAL(KIND=rsh) :: Tchl, fact_phyto_ChlNratio

            IF(l_ChlNratio_var) THEN
                  IF((CURRENT_TIME-TIME_BEGIN) > 345600.00_rlg) THEN
                      ! Cas ou on transforme l azote en chloro par un rapport Chloro/N dependant de l extinction
                      !   lumineuse moyennee sur les 4 jours precedents  et de la limitation en azote
                      fact_phyto_ChlNratio=p_phyto_ChlNratiomax*extinction_ave4d(k,i,j)/p_phyto_ChlN_ksmithextinct/ &
                          (sqrt(1+(extinction_ave4d(k,i,j)/p_phyto_ChlN_ksmithextinct)**2))
                      fact_phyto_ChlNratio=fact_phyto_ChlNratio*(0.5_rsh*(1.0_rsh-effetazotediat)+1.5_rsh*effetazotediat)
                  ELSE
                      fact_phyto_ChlNratio=p_phyto_ChlNratiomax*effetturbidite
                  ENDIF
                  IF(fact_phyto_ChlNratio<0.7_rsh) fact_phyto_ChlNratio=0.7_rsh
              ELSE
                  fact_phyto_ChlNratio=p_phyto_ChlNratio
              ENDIF

              Tchl = (c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)) *fact_phyto_ChlNratio

          end function get_Tchl



          !--------------------------------------------------------------------------------------------------
          ! Storing diagnostic variables
          !--------------------------------------------------------------------------------------------------

          SUBROUTINE storage_diag_phytoplankton(i,j,k,c)

              INTEGER, intent(in)                 :: i,j,k
              REAL(KIND=rsh)             ::fact_phyto_ChlNratio
              REAL(KIND=rsh),DIMENSION(nv_state), intent(in)  :: c


              effetlumiere_day_diat(k,i,j) = effetlumiere_day_diat(k,i,j) + effetlumierediat*dtbio
              diag_3d_wat(irk_diag(id_diat_limlight),k,i,j)=effetlumierediat

              effetlumiere_day_dino(k,i,j)=effetlumiere_day_dino(k,i,j)+effetlumieredino*dtbio
              diag_3d_wat(irk_diag(id_dino_limlight),k,i,j)=effetlumieredino

              effetlumiere_day_nano(k,i,j)=effetlumiere_day_nano(k,i,j)+effetlumierenano*dtbio
              diag_3d_wat(irk_diag(id_nano_limlight),k,i,j)=effetlumierenano

              diag_3d_wat(irk_diag(id_diat_limN),k,i,j)=effetazotediat
              diag_3d_wat(irk_diag(id_diat_limSi),k,i,j)=effetsilice
              diag_3d_wat(irk_diag(id_diat_limP),k,i,j)=effetphosphorediat

              diag_3d_wat(irk_diag(id_dino_limN),k,i,j)=effetazotedino
              diag_3d_wat(irk_diag(id_dino_limP),k,i,j)=effetphosphoredino

              diag_3d_wat(irk_diag(id_nano_limN),k,i,j)=effetazotenano
              diag_3d_wat(irk_diag(id_nano_limP),k,i,j)=effetphosphorenano

              ! chlorophylle [mgChl/m3]
              
              
              diag_3d_wat(irk_diag(id_totalchl),k,i,j)= get_Tchl(i,j,k,c) 



          END SUBROUTINE storage_diag_phytoplankton





END MODULE phytoplankton





