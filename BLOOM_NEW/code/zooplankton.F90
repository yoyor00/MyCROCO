
MODULE zooplankton

#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"

    USE module_BIOLink
    USE comBIOLink
    USE comBIOLink_helping
    USE COMBLOOM

    USE shared !, only: effetchaleur,effetturbidite, l_BLOOM_opt2

      IMPLICIT NONE

      PUBLIC :: zooplankton_dynamics

      REAL(KIND=rsh) :: rationmesozoo, assimilmesozoo,excretionmesozoo,txmortmesozoo
      REAL(KIND=rsh) :: rationmicrozoo,excretionmicrozoo, assimilmicrozoo, txmortmicrozoo
      REAL(KIND=rsh) :: broumicrozoodet,broumicrozoodiat, broumicrozoonano, broumicrozoodino
      REAL(KIND=rsh) :: broumesozoodiat,broumesozoodino,broumesozoomicrozoo
      REAL(KIND=rsh) :: rNC


      contains

          SUBROUTINE zooplankton_dynamics(i,j,k,c,dc)

              INTEGER, intent(in)                 :: i,j,k
              REAL(KIND=rsh),DIMENSION(nv_state), intent(in)  :: c
              REAL(KIND=rsh),DIMENSION(nv_state), intent(inout) :: dc
              !REAL(KIND=rsh) , intent(in) :: effetchaleur,effetturbidite

              REAL(KIND=rsh) :: totproiebrut,totproie,fmesozoo
              REAL(KIND=rsh) :: captdiatmeso,captdino,captmicrozoo,captmicrozooN
              REAL(KIND=rsh) :: meso_kivlev, p_mesz_thrN_turb

              REAL(KIND=rsh) :: fmicrozoo,captnano,captdet,captdinomicro
              REAL(KIND=rsh) :: captdiatmicro
              REAL(KIND=rsh) :: totproie_thr, orespzoo ! Mokrane


              ! Mesozooplankton

              captdiatmeso=c(iv_phyto_diat_N)*p_mesz_captdiat    !en microMolN.l-1
              captdino=c(iv_phyto_dino_N)*p_mesz_captdino        !en microMolN.l-1
              captmicrozoo=c(iv_zoo_micr_N)*p_mesz_captmicz      !en microMolN.l-1
              rNC=1.0_rsh
              captmicrozooN=captmicrozoo*rNC

              p_mesz_thrN_turb=p_mesz_thrN*effetturbidite
              meso_kivlev=p_mesz_kivlev

              ! Microzooplankton

              captnano=c(iv_phyto_nano_N)*p_micz_captnano        !en microMolN.l-1
              captdet=min(c(iv_detr_N),c(iv_detr_P)*p_phyto_NPratio)*p_micz_captdet   !en microMolN.l-1
              captdinomicro=c(iv_phyto_dino_N)*p_micz_captdino   !en microMolN.l-1
              captdiatmicro=c(iv_phyto_diat_N)*p_micz_captdiat   !en microMolN.l-1
              totproie_thr = 0.0_rsh   ! Mokrane : new parameter

              IF(l_physadaptation)  THEN
                  captdiatmicro= captdiatmicro * (0.1_rsh+0.9_rsh*min(200.0_rsh,BATHY_H0(i,j))/200.0_rsh)
              ENDIF
              

              IF(l_BLOOM_opt2) THEN  ! Mokrane : Dans ce cas le zoo est exprimé en carbone d'où sa conversion en azote
                  !Mesozooplankton
                  rNC=1.0_rsh/(12.0_rsh*p_zoo_CNratio)                 !en micromolN.microgC-1
                  IF (c(iv_phyto_diat_N) .le. 0.01_rsh) captdiatmeso= 0.0_rsh
                  IF (c(iv_phyto_dino_N) .le. 0.01_rsh) captdino=0.0_rsh
                  IF((c(iv_zoo_micr_N)*rNC) .le. 0.01_rsh) captmicrozoo=0.0_rsh
                  p_mesz_thrN_turb = 0.0_rsh
                  totproie_thr = 0.005_rsh   ! Mokrane : new parameter TODO : declare it


                  ! Microzooplankton
                  IF (c(iv_phyto_nano_N) .le. 0.01_rsh) captnano=0.0_rsh
                  IF (c(iv_detr_N).le.0.01_rsh) captdet=0.0_rsh
                  IF(c(iv_phyto_diat_N) .le. 0.01_rsh) captdiatmicro=0.0_rsh
                  captdinomicro=0.0_rsh
                  p_micz_thrnano = 0.0_rsh

              ENDIF

              ! Meso
              totproiebrut = captdiatmeso + captdino + captmicrozoo*rNC
              totproie=max(0.0_rsh,(totproiebrut - p_mesz_thrN_turb)) ! en microgMolN.l-1
              IF (totproie .gt. totproie_thr) THEN

                  IF(l_physadaptation)  THEN
                      IF (cmes_3dmgl(k,i,j).gt.p_mesz_thhold_mes_kivlev) meso_kivlev=p_mesz_kivlev*min(1.0_rsh,(1.0_rsh-effetturbidite))
                  ENDIF

                  fmesozoo=1.0_rsh-exp(-meso_kivlev*totproie)              ! food lim [0-1]
                  rationmesozoo=p_mesz_mumax*effetchaleur*fmesozoo     ! specif growth rate d-1
                  ! assimilation depend de quantitee ingeree Hofmann & Hambler, 1988 J.M.R 
                  assimilmesozoo=0.3_rsh*(3.0_rsh-0.67_rsh*fmesozoo)       ! ass rate s.u.
                  broumesozoodiat=rationmesozoo*captdiatmeso/totproiebrut ! d-1
                  broumesozoodino=rationmesozoo*captdino/totproiebrut ! d-1
                  broumesozoomicrozoo=rationmesozoo*captmicrozooN/totproiebrut ! d-1
              ELSE
                  totproiebrut=0.0_rsh
                  totproie=0.0_rsh
                  fmesozoo=0.0_rsh
                  assimilmesozoo=0.0_rsh
                  rationmesozoo=0.0_rsh
                  broumesozoodiat=0.0_rsh
                  broumesozoodino=0.0_rsh
                  broumesozoomicrozoo=0.0_rsh
              ENDIF
              excretionmesozoo=p_mesz_excret*effetchaleur*fmesozoo ! specif excret rate d-1
              txmortmesozoo=(p_mesz_mort1+p_mesz_mort2*c(iv_zoo_meso_N))*effetchaleur ! specif mort rate d-1

              IF (c(iv_zoo_meso_N).lt.p_mesz_thhold_mort)  txmortmesozoo=0.0_rsh

              ! Micro
              totproiebrut=captdiatmicro+captnano+captdet+captdinomicro
              totproie=max(0.0_rsh,(totproiebrut-p_micz_thrnano)) ! en microgMolN.l-1
              IF (totproie .gt. totproie_thr) THEN

                  fmicrozoo=totproie/(p_micz_kgraz+totproie)   ! food lim [0-1]
                  rationmicrozoo=p_micz_mumax*effetchaleur*fmicrozoo  ! specif growth rate d-1
                  assimilmicrozoo=0.3_rsh*(3.0_rsh-0.67_rsh*fmicrozoo)        ! ass rate s.u.

                  IF(l_BLOOM_opt2) assimilmicrozoo=p_micz_assim

                  broumicrozoonano=rationmicrozoo*captnano/totproiebrut ! d-1
                  broumicrozoodet=rationmicrozoo*captdet/totproiebrut ! d-1
                  broumicrozoodiat=rationmicrozoo*captdiatmicro/totproiebrut ! d-1
                  broumicrozoodino=rationmicrozoo*captdinomicro/totproiebrut ! d-1

              ELSE
                  fmicrozoo=0.0_rsh
                  assimilmicrozoo=0.0_rsh
                  rationmicrozoo=0.0_rsh
                  broumicrozoonano=0.0_rsh
                  broumicrozoodet=0.0_rsh
                  broumicrozoodiat=0.0_rsh
                  broumicrozoodino=0.0_rsh
              ENDIF
              excretionmicrozoo=p_micz_excret*effetchaleur*fmicrozoo ! specif excret rate d-1
              txmortmicrozoo=(p_micz_mort+p_mesz_mort2*c(iv_zoo_micr_N))*effetchaleur  !specif excret rate d-1

              IF (c(iv_zoo_micr_N).lt.p_micz_thhold_mort) txmortmicrozoo=0.0_rsh

              orespzoo=p_zoo_resp*(c(iv_zoo_meso_N) + c(iv_zoo_micr_N))
              IF (c(iv_oxygen).lt.0.2_rsh) THEN
                  txmortmicrozoo=100.0_rsh*txmortmicrozoo
                  txmortmesozoo=100.0_rsh*txmortmesozoo
                  orespzoo = 0._rsh
              ENDIF

              !++++++++++++++++++++++++++++++++++++++++++++++++++
              ! EQUATIONS D EVOLUTION  ! concentration en jour-1
              !++++++++++++++++++++++++++++++++++++++++++++++++++

              dc(iv_zoo_meso_N) = dc(iv_zoo_meso_N)     &
                  + (rationmesozoo*assimilmesozoo-excretionmesozoo-txmortmesozoo)*c(iv_zoo_meso_N)

              dc(iv_zoo_micr_N) = dc(iv_zoo_micr_N)     &
                  + c(iv_zoo_micr_N)*(rationmicrozoo*assimilmicrozoo-excretionmicrozoo) &
                  - txmortmicrozoo*c(iv_zoo_micr_N)     &
                  - broumesozoomicrozoo*c(iv_zoo_meso_N)

              dc(iv_phyto_nano_N) = dc(iv_phyto_nano_N) &
                  - broumicrozoonano*c(iv_zoo_micr_N)*rNC

              dc(iv_phyto_diat_N) = dc(iv_phyto_diat_N) &
                  - broumesozoodiat*c(iv_zoo_meso_N)*rNC &
                  - broumicrozoodiat*c(iv_zoo_micr_N)*rNC

              dc(iv_phyto_dino_N) = dc(iv_phyto_dino_N) &
                  - broumesozoodino*c(iv_zoo_meso_N)*rNC &
                  - broumicrozoodino*c(iv_zoo_micr_N)*rNC

              dc(iv_nutr_NH4) = dc(iv_nutr_NH4)               &
                  + (excretionmesozoo*c(iv_zoo_meso_N)        &
                  + excretionmicrozoo*c(iv_zoo_micr_N))*rNC

              dc(iv_nutr_PO4) = dc(iv_nutr_PO4)                &
                  + rappaz*((excretionmesozoo*c(iv_zoo_meso_N) &
                  + excretionmicrozoo*c(iv_zoo_micr_N))*rNC)

              dc(iv_detr_N) = dc(iv_detr_N)         &
                  + (txmortmesozoo+(1.0_rsh-assimilmesozoo)*rationmesozoo)*c(iv_zoo_meso_N)*rNC   &
                  + (txmortmicrozoo-broumicrozoodet)*c(iv_zoo_micr_N)*rNC

              dc(iv_detr_Si) = dc(iv_detr_Si)       &
                  + broumesozoodiat*c(iv_zoo_meso_N)*rNC*rapsiaz

              dc(iv_detr_P) = dc(iv_detr_P)         &
                  - broumicrozoodet*c(iv_zoo_micr_N)*rNC*rappaz


#if defined key_BLOOM_opt2
                !IF(l_BLOOM_opt2) THEN
                  dc(iv_nutr_NH4)= dc(iv_nutr_NH4)  &
                      + (1.0_rsh-assimilmicrozoo)*rationmicrozoo*(1.0_rsh-p_micz_diss)*c(iv_zoo_micr_N)*rNC
                  dc(iv_nutr_PO4) = dc(iv_nutr_PO4) &
                      + (1.0_rsh-assimilmicrozoo)*rationmicrozoo*(1.0_rsh-p_micz_diss)*c(iv_zoo_micr_N)*rappaz*rNC
                  dc(iv_detr_P) = dc(iv_detr_P)         &
                      + txmortmesozoo*c(iv_zoo_meso_N)*rNC*rappaz &
                      + (1.0_rsh-assimilmesozoo)*rationmesozoo*c(iv_zoo_meso_N)*rNC*rappaz &
                      + txmortmicrozoo*c(iv_zoo_micr_N)*rNC*rappaz
                  dc(iv_diss_N) = dc(iv_diss_N)         &
                      + (1.0_rsh-assimilmicrozoo)*rationmicrozoo*p_micz_diss*c(iv_zoo_micr_N)*rNC  
                  dc(iv_diss_Si)= dc(iv_diss_Si)        &
                      + broumicrozoodiat*c(iv_zoo_micr_N)*rNC*rapsiaz
                  dc(iv_diss_P) = dc(iv_diss_P)         &
                      + (1.0_rsh-assimilmicrozoo)*rationmicrozoo*p_micz_diss*c(iv_zoo_micr_N)*rNC*rappaz


              !ELSE
#else
                  dc(iv_detr_N) = dc(iv_detr_N)         &
                      + (1.0_rsh-assimilmicrozoo)*rationmicrozoo*c(iv_zoo_micr_N)*rNC
                  dc(iv_detr_Si) = dc(iv_detr_Si)       &
                      + broumicrozoodiat*c(iv_zoo_micr_N)*rNC*rapsiaz
                  dc(iv_detr_P) = dc(iv_detr_P)         &
                      + rappaz * rNC * ((txmortmesozoo+(1.0_rsh-assimilmesozoo)*rationmesozoo)*c(iv_zoo_meso_N) &
                      + (txmortmicrozoo+(1.0_rsh-assimilmicrozoo)*rationmicrozoo)*c(iv_zoo_micr_N))
              !ENDIF
#endif

              

                ! Evolution of dissolved oxygen
                !------------------------------
                dc(iv_oxygen) = dc(iv_oxygen) &
                    - orespzoo * ratio_mgO_to_mumolN_catabol * effetchaleur
                



          END SUBROUTINE zooplankton_dynamics

END MODULE zooplankton



