MODULE SINKING_VELOCITY

        !&E--------------------------------------------------------------------
        !&E                 ***  incellwat_bloom_settling  ***
        !&E
        !&E ** Purpose : estimate settling velocity for biologic variable (detrital..)
        !&E              a l interieur d une maille d eau
        !&E
        !&E       !  2010-03    (B. Thouvenin issued fro bloomdynzwat ) Original code
        !&E       !  2019-07    (B.Thouvenin ) update
        !&E       !  2026-05    (M.Belharet ) Modularization
        !&E
        !&E      use from general modele : temper,sali, WATER_DENSITY  si ca existe, (h0 mais on pourrait mettre htot)
        !&E      use from BIOLink  :
        !&E
        !&E      use from basic bloom modele : c, epn, l_SNeffect_settle, effetselnutdiat, l_phyzoodeteffect_settle,
        !&E                                    diatmorteau, dinomorteau, nanomorteau, rationmicrozoo, assimilmicrozoo, txmortmicrozoo
        !&E                                    rationmesozoo, assimilmesozoo, txmortmesozoo, rationmicrozoo, assimilmicrozoo, txmortmicrozoo
        !&E      use from options bloom modules : effetselnutpsnz, effetselnutphaeocystiscolo, effetlumierephaeocystis
        !&E                                       psnzmorteau, kareniamorteau, phaeocystismortcolo, phaeocystislysecolo, phaeocystismortcell
        !&E      OUTPUT : ws3
        !&E---------------------------------------------------------------------
#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"
    
        USE module_BIOLink
        USE comBIOLink
        USE COMBLOOM

        USE comBIOLink_helping
        USE shared
        USE phytoplankton , only : effetselnutdiat, diatmorteau, dinomorteau, nanomorteau
        USE zooplankton   , only:  rationmicrozoo, assimilmicrozoo, txmortmicrozoo, rationmesozoo , &
                                assimilmesozoo, txmortmesozoo


        IMPLICIT NONE

        PRIVATE

        PUBLIC :: estimate_sinking_velocity, storage_diag_sinking_velocity

       contains
              SUBROUTINE estimate_sinking_velocity(i,j,k,c)
                
                INTEGER, intent(in)                 :: i,j,k
                REAL(KIND=rsh), DIMENSION(nv_state), INTENT(in)          :: c
!                REAL(KIND=rsh), ALLOCATABLE, INTENT(inout) :: ws3(:,:,:,:)
                REAL(KIND=rsh)                             :: modulationsedpardenseau,rdet,wchutedet
                REAL(KIND=rsh)                             :: phytodetritus,zoodetritus,effetchutedia,abovesinkingrate
                INTEGER                                    :: ivp 
                

                IF (l_SNeffect_settle) THEN
                        !  Phytoplancton
                        !===============
                        ws3(k,iv_phyto_diat_N,i,j) = ws_free_max(iv_phyto_diat_N)
                        effetchutedia = MAX(0.0_rsh,effetselnutdiat)**0.2_rsh
                        ws3(k,iv_phyto_diat_N,i,j) = (ws_free_min(iv_phyto_diat_N)*effetchutedia+      &
                                ws3(k,iv_phyto_diat_N,i,j)*(1.0_rsh-effetchutedia))
                        
                END IF

                ! Materiel detritique
                !======================
                ! modulation en fonction du mesozooplancton

                IF(l_phyzoodeteffect_settle .and. c(iv_zoo_meso_N)>0.0_rsh) THEN
                        ! PHYTODETRITUS
                        phytodetritus=diatmorteau*c(iv_phyto_diat_N)  &
                                +dinomorteau*c(iv_phyto_dino_N)  &
                                +nanomorteau*c(iv_phyto_nano_N)  
                                
                        ! ZOODETRITUS
                        zoodetritus=(rationmesozoo*(1.0_rsh-assimilmesozoo)+txmortmesozoo)*c(iv_zoo_meso_N)   &
                                        +(rationmicrozoo*(1.0_rsh-assimilmicrozoo)+txmortmicrozoo)*c(iv_zoo_micr_N)
                        IF (zoodetritus > 0.0_rsh) THEN
                                rdet=phytodetritus/(zoodetritus+epsilon_BIOLink)
                                wchutedet=p_detzoo_wsed*(1.0_rsh/(rdet+1.0_rsh))+      &
                                        p_detphy_wsed*(1.0_rsh-(1.0_rsh/(rdet+1.0_rsh)))
                        ELSE
                                wchutedet=p_detphy_wsed
                        ENDIF
                ELSE
                        wchutedet=p_detphy_wsed
                ENDIF

                ws3(k,iv_detr_N,i,j)=wchutedet
                ws3(k,iv_detr_Si,i,j)=wchutedet
                ws3(k,iv_detr_P,i,j)=wchutedet

                ! bornage de la vitesse de chute des variables particulaires en chaque maille
                ! mise a zero aux limites ouvertes
                ! ---------------------------------------------------------------------------
                DO ivp=1,nvp
                        ws3(k,ivp,i,j)=sign(MIN(0.95_rlg*epn/dtbio,REAL(ABS(ws3(k,ivp,i,j)),rlg)),ws3(k,ivp,i,j))
                END DO

              
              END SUBROUTINE estimate_sinking_velocity 

              SUBROUTINE storage_diag_sinking_velocity(i,j,k)
                      
                      INTEGER, intent(in)                 :: i,j,k

                      diag_3d_wat(irk_diag(id_diatsettling),k,i,j)=ws3(k,iv_phyto_diat_N,i,j)*86400.0_rsh
                      diag_3d_wat(irk_diag(id_detsettling),k,i,j)=ws3(k,iv_detr_N,i,j)*86400.0_rsh

              END SUBROUTINE storage_diag_sinking_velocity




END MODULE SINKING_VELOCITY
