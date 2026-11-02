
MODULE chemistry
#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"

    USE module_BIOLink
    USE comBIOLink
    USE comBIOLink_helping
    USE COMBLOOM
    USE shared 

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: chemistry_processes

    REAL(KIND=rsh) :: desorpeau,adsorpeau
    REAL(KIND=rsh) :: dissolMOPeau
    REAL(KIND=rsh) :: fractionDET_N, fractionDET_P, fractionDET_Si, &
        fractionDIAT_N, fractionDIAT_P, fractionDIAT_Si
    REAL(KIND=rsh) :: xnitrifeau,reminazdeteau,reminpdeteau,dissolsiliceeau


    contains
        SUBROUTINE chemistry_processes(i,j,k,c,dc) 

            INTEGER, intent(in)                 :: i,j,k
            REAL(KIND=rsh),DIMENSION(nv_state), intent(in)  :: c
            REAL(KIND=rsh),DIMENSION(nv_state), intent(inout)  :: dc

            REAL(KIND=rsh) :: varads, cmaxdesorp, saturmesp,cmaxadsorp
            REAL(KIND=rsh) :: p_det_fragm


            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ! Remineralisation de la matiere detritique dans l eau
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            reminazdeteau=p_N_remin*effetchaleur
            xnitrifeau=p_nitrif*effetchaleur
            IF(l_BLOOM_insed) xnitrifeau=p_nitrif/10*effetchaleur*c(iv_oxygen)/(c(iv_oxygen)+1.0)  ! nitrification much slower in wat than in sed
            reminpdeteau=p_P_remin*effetchaleur
            dissolsiliceeau=p_BSi_dissEau*effetchaleur
            IF(l_BLOOM_opt2) dissolMOPeau=p_det_fragm*effetchaleur

            IF(c(iv_oxygen).lt.0.2_rsh) THEN
                reminazdeteau=0.0_rsh
                xnitrifeau=0.0_rsh
                reminpdeteau=0.0_rsh
            ENDIF

            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ! Adsorption et desorption du phosphore dans l eau
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++

            ! adsorption
            !p_P_adsormaxspim = max(5.0_rsh,p_P_adsormaxspim*exp(-0.35*c(iv_spim)))

            saturmesp=max(0.0_rsh,(p_P_adsormaxspim*c(iv_spim)-c(iv_nutr_Pads))) !c(iv_spim) en g/l
            adsorpeau=p_P_adsor*saturmesp
            cmaxadsorp=adsorpeau*c(iv_nutr_PO4)*dtbiojour
            IF (cmaxadsorp.ge.c(iv_nutr_PO4))  adsorpeau=0.9_rsh/dtbiojour


            ! desorption
            varads=0.0_rsh
            IF(c(iv_spim).gt.0.000001_rsh) varads=min(1.0_rsh,c(iv_nutr_Pads)/(p_P_adsormaxspim*c(iv_spim)))
            desorpeau=p_P_desor*varads
            cmaxdesorp=desorpeau*c(iv_nutr_Pads)*dtbiojour
            IF (cmaxdesorp.ge.c(iv_nutr_Pads)) desorpeau=0.9_rsh/dtbiojour

            IF(l_BLOOM_opt2) THEN
                fractionDIAT_N=0.0_rsh
                fractionDIAT_P=0.0_rsh
                fractionDIAT_Si=0.0_rsh
                fractionDET_N=0.0_rsh
                fractionDET_Si=0.0_rsh
                fractionDET_P=0.0_rsh

                IF (k.eq.1) THEN
                    fractionDIAT_N=c(iv_phyto_diat_N)*(min((ws3(2,iv_phyto_diat_N,i,j)*dtbio/epn),1.0_rsh))
                    fractionDIAT_P=fractionDIAT_N*rappaz
                    fractionDIAT_Si=fractionDIAT_N*rapsiaz
                    fractionDET_N=c(iv_detr_N)*(min((ws3(2,iv_detr_N,i,j)*dtbio/epn),1.0_rsh))
                    fractionDET_Si=c(iv_detr_Si)*(min((ws3(2,iv_detr_Si,i,j)*dtbio/epn),1.0_rsh))
                    fractionDET_P=c(iv_detr_P)*(min((ws3(2,iv_detr_P,i,j)*dtbio/epn),1.0_rsh))
                ENDIF

            ENDIF

            !++++++++++++++++++++++++++++++++++++++++++++++++++
            ! EQUATIONS D EVOLUTION  ! concentration en jour-1
            !++++++++++++++++++++++++++++++++++++++++++++++++++

            dc(iv_nutr_NO3)= dc(iv_nutr_NO3)       &
                + xnitrifeau*c(iv_nutr_NH4)

            dc(iv_nutr_NH4) = dc(iv_nutr_NH4) &
                - xnitrifeau*c(iv_nutr_NH4) 

            dc(iv_nutr_PO4) = dc(iv_nutr_PO4) &
                + desorpeau*c(iv_nutr_Pads)   &
                - adsorpeau*c(iv_nutr_PO4)    

            dc(iv_nutr_Pads) = dc(iv_nutr_Pads) &
                + adsorpeau*c(iv_nutr_PO4)      &
                - desorpeau*c(iv_nutr_Pads)



            !IF(l_BLOOM_opt2) THEN
#if defined key_BLOOM_opt2
                dc(iv_nutr_NH4)  = dc(iv_nutr_NH4)   + reminazdeteau*c(iv_diss_N)
                dc(iv_nutr_SiOH) = dc(iv_nutr_SiOH)  + dissolsiliceeau*c(iv_diss_Si)
                dc(iv_nutr_PO4)  = dc(iv_nutr_PO4)   + reminpdeteau*c(iv_diss_P)
                dc(iv_detr_N)    = dc(iv_detr_N)  - dissolMOPeau*c(iv_detr_N) - fractionDET_N
                dc(iv_detr_Si)   = dc(iv_detr_Si) - dissolMOPeau*c(iv_detr_Si)- fractionDET_Si
                dc(iv_detr_P)    = dc(iv_detr_P)  - dissolMOPeau*c(iv_detr_P) - fractionDET_P
                dc(iv_diss_N)    = dc(iv_diss_N)             &
                    - reminazdeteau*c(iv_diss_N)             &
                    + dissolMOPeau*c(iv_detr_N)              &
                    + 0.25* (fractionDIAT_N + fractionDET_N) 
                dc(iv_diss_fond_Nitr) = dc(iv_diss_fond_Nitr) &
                    + 0.75 * (fractionDIAT_N + fractionDET_N)
                dc(iv_diss_Si) = dc(iv_diss_Si)       &
                    - dissolsiliceeau*c(iv_diss_Si)   &
                    + dissolMOPeau*c(iv_detr_Si)      &
                    + fractionDET_Si + fractionDIAT_Si
                dc(iv_diss_P) = dc(iv_diss_P)         &
                    - reminpdeteau*c(iv_diss_P)       &
                    + dissolMOPeau*c(iv_detr_P)       &
                    + fractionDET_P + fractionDIAT_P

                dc(iv_phyto_diat_N) = dc(iv_phyto_diat_N) &
                    - fractionDIAT_N
#else
            !ELSE
                dc(iv_nutr_NH4)  = dc(iv_nutr_NH4)   + reminazdeteau*c(iv_detr_N)
                dc(iv_nutr_SiOH) = dc(iv_nutr_SiOH)  + dissolsiliceeau*c(iv_detr_Si)
                dc(iv_nutr_PO4)  = dc(iv_nutr_PO4)    + reminpdeteau*c(iv_detr_P)
                dc(iv_detr_N)    = dc(iv_detr_N)     - reminazdeteau*c(iv_detr_N)
                dc(iv_detr_Si)   = dc(iv_detr_Si)    - dissolsiliceeau*c(iv_detr_Si)
                dc(iv_detr_P)    = dc(iv_detr_P)     - reminpdeteau*c(iv_detr_P)
            !ENDIF
#endif

            ! Evolution of dissolved oxygen
            !------------------------------
            dc(iv_oxygen) = dc(iv_oxygen) &
                - reminazdeteau*c(iv_detr_N)*ratio_mgO_to_mumolN_catabol &
                - xnitrifeau*c(iv_nutr_NH4)*0.064_rsh


        END SUBROUTINE chemistry_processes



END MODULE chemistry

