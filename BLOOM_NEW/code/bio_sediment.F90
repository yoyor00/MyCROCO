MODULE bio_sediment

#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"

    USE module_BIOLink
    USE comBIOLink
    USE comBIOLink_helping
    USE COMBLOOM
#if defined MUSTANG
   USE comMUSTANG
#endif


    IMPLICIT NONE

    PRIVATE

    PUBLIC :: sediment_bio_dynamics

    contains
        SUBROUTINE sediment_bio_dynamics(i,j,k,dtbiojour, zmiddle, cw_bott)

            INTEGER, intent(in)  :: i, j, k
            REAL(KIND=rsh), INTENT(in) :: dtbiojour
            REAL(KIND=rsh), INTENT(inout) :: zmiddle
            REAL(KIND=rsh),DIMENSION(nv_state) :: dcw_filtbent, dcdt, cs 
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: cw_bott
            REAL(KIND=rsh) :: foncsinus, ksmin, ksmax,txfiltbenth_max
            REAL(KIND=rsh) :: flim1_O2, flim2_O2, flim3_O2, glim1_O2, glim2_O2,flim1_NO3, glim1_NO3, Sfliminv
            REAL(KIND=rsh) :: temper, effetchaleur, porosite_inv

            ksmin=ksmi(i,j)
            ksmax=ksma(i,j)

            txfiltbenth_max=p_txfiltbenthmax


            ! INITIALIZE THE FLUXES TO 0
            dcdt(:)=0.0_rsh

            ! transfer 3D var cv_sed into 1D var c
            cs(1:nv_state) = 0.5_rsh*( cv_sed(1:nv_state,k,i,j) + ABS(cv_sed(1:nv_state,k,i,j)) )
            cs(1:nv_state) = cs(1:nv_state)/unit_modif_mudbio_N2dw(irk_fil(1:nv_state))

            porosite_inv=1._rsh/poro(k,i,j)

            !diag_3d_sed(id_porosite_sed,k,i,j)=poro(k,i,j) ! TODO

            ! ztop=-dzs(ksmax,i,j) Mokrane : pas utilisé

            zmiddle=zmiddle+dzs(k+1,i,j)/2._rsh+dzs(k,i,j)/2._rsh
            if(k==ksmax) zmiddle=dzs(k,i,j)/2


            dcw_filtbent(:)=0.0_rsh
            IF(k == ksma(i,j) .AND. htot(i,j) > RESIDUAL_THICKNESS_WAT .AND. txfiltbenth_max .NE. 0.0_rsh) CALL filtr_on_wat(i,j,k,dtbiojour,cs,dcdt, dcw_filtbent, cw_bott, txfiltbenth_max)

            !$$-- Temperature limitation --$$

            temper=cv_sed(-1,k,i,j)

            effetchaleur=exp(p_T_effect*temper)

            !$$-- Oxygen limitation & inhibition --$$ 

            flim1_O2 = 0.0_rsh
            flim2_O2 = 0.0_rsh
            flim3_O2 = 0.0_rsh
            glim1_O2 = 1.0_rsh
            glim2_O2 = 1.0_rsh

            IF(cs(iv_oxygen).gt.0.0_rsh) THEN
                flim1_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_reminO2)
                flim2_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_nit)
                flim3_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_reoxyd)
                glim1_O2 = 1._rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_remin0O2)
                glim2_O2 = 1._rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_denit)
            ENDIF

            !$$-- Nitrate limitation & inhibition --$$
            flim1_NO3 = 0.0_rsh
            glim1_NO3 = 1.0_rsh
            IF(cs(iv_nutr_NO3).gt.0.0_rsh) THEN
                flim1_NO3 = cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kNO3_reminssO2)
                glim1_NO3 = 1._rsh - cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kiNO3_remin0O2)
            ENDIF

            !$$-- Total limitation --$$

            Sfliminv = 1.0_rsh / (flim1_O2 + flim1_NO3 * glim2_O2 + glim1_NO3 * glim1_O2)


            CALL nitrogen_cycle(i,j,k,dtbiojour,cs,dcdt,zmiddle,flim1_O2,flim2_O2,flim3_O2,glim1_O2,glim2_O2,flim1_NO3,glim1_NO3,Sfliminv,effetchaleur,porosite_inv)
            CALL phosphorus_cycle(i,j,k,dtbiojour,cs,dcdt,zmiddle,flim1_O2,glim1_O2,glim2_O2,flim1_NO3,glim1_NO3,Sfliminv,effetchaleur,porosite_inv)
            CALL silicate_cycle(i,j,k,dtbiojour,cs,dcdt,zmiddle,temper,effetchaleur,porosite_inv)
            CALL MO_aging(k,dtbiojour,cs,dcdt,ksmin)

            CALL update_state(i,j,k,dtbiojour,cs,dcdt,dcw_filtbent,zmiddle,ksmax)



        END SUBROUTINE sediment_bio_dynamics


        SUBROUTINE filtr_on_wat(i,j,k,dtbiojour,cs, dcdt, dcw_filtbent, cw_bott, txfiltbenth_max)

            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(inout)  :: dcdt, dcw_filtbent
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: cw_bott, cs

            INTEGER, intent(in)  :: i, j, k
            REAL(KIND=rsh), INTENT(in) :: dtbiojour, txfiltbenth_max
            REAL(KIND=rsh) :: epn_bott, factlim, txfiltbenth_wat, filtr_benth, txfiltbenth

            epn_bott=epn_bottom_MUSTANG(i,j)
#ifdef key_MARS
            IF (htot(i,j) < hm ) epn_bott=htot(i,j)
#endif      

            txfiltbenth=MAX(-0.0333_rsh*cmes_3dmgl(1,i,j)+1,0.0_rsh)*txfiltbenth_max !de 0 mg/L a 30 mg/L
            ! conversion par metre3 eau
            txfiltbenth_wat=txfiltbenth/epn_bott
            ! resolution explicite pour cw_bott ==> limitation pas de temps
            factlim=MIN(1.0_rsh,1.0_rsh/dtbiojour/(txfiltbenth_wat+epsilon_BIOLink))
            IF(factlim .NE. 1.0_rsh) write(*,*)'filt',(CURRENT_TIME-TIME_BEGIN)/3600._rsh, &
                i,j,k,epn_bott,factlim,txfiltbenth_wat,dtbiojour

            txfiltbenth_wat=txfiltbenth_wat*factlim
            txfiltbenth=txfiltbenth*factlim/dzs(k,i,j)
            filtr_benth= txfiltbenth*(cw_bott(iv_phyto_diat_N)+cw_bott(iv_phyto_dino_N)+cw_bott(iv_phyto_nano_N))

            dcw_filtbent(iv_phyto_diat_N)=  -txfiltbenth_wat*cw_bott(iv_phyto_diat_N)
            dcw_filtbent(iv_phyto_dino_N)=  -txfiltbenth_wat*cw_bott(iv_phyto_dino_N)
            dcw_filtbent(iv_phyto_nano_N)=  -txfiltbenth_wat*cw_bott(iv_phyto_nano_N)

            dcdt(iv_detr_P) = dcdt(iv_detr_P) &
                + filtr_benth * rappaz * dtbiojour

            dcdt(iv_detr_Si) = dcdt(iv_detr_Si) &
                +rapsiaz * txfiltbenth * cw_bott(iv_phyto_diat_N) * dtbiojour

            diag_3d_sed(id_filtr_benth,k,i,j)=filtr_benth * dtbiojour * dzs(k,i,j)

          

        END SUBROUTINE filtr_on_wat

        SUBROUTINE nitrogen_cycle(i,j,k,dtbiojour,cs,dcdt,zmiddle,flim1_O2,flim2_O2,flim3_O2,glim1_O2,glim2_O2,flim1_NO3,glim1_NO3,Sfliminv,effetchaleur,porosite_inv)

            INTEGER, intent(in)  :: i, j, k
            REAL(KIND=rsh), INTENT(in) :: dtbiojour
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  ::cs

            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(inout)  :: dcdt
            REAL(KIND=rsh), INTENT(in) :: zmiddle,flim1_O2,flim2_O2,flim3_O2,glim1_O2,glim2_O2,flim1_NO3,glim1_NO3,Sfliminv,effetchaleur,porosite_inv

            REAL(rsh) :: xtmp, F_remin_aerN, flimz, F_reminR_aerN
            REAL(rsh) :: F_nitrif, F_remin_NO3, F_reminR_NO3
            REAL(rsh) :: F_denit, F_denitR, F_dnra, F_dnraR
            REAL(rsh) :: F_oxyd_ODU, F_solid_ODU
            REAL(rsh) :: F_remin_anaerN, F_reminR_anaerN

           !$$-- Aerobic Mineralisation --$$
           xtmp = effetchaleur * flim1_O2 * dtbiojour
           F_remin_aerN = p_N_remin * xtmp * Sfliminv ! Labile  
           flimz = max(p_xflimz,exp(-p_k_remin*zmiddle))
           F_reminR_aerN = p_N_reminR * flimz * xtmp * Sfliminv ! Refractaire

           !$$-- Nitrification (NH4 --> NO3) --$$
           F_nitrif = p_nitrif * effetchaleur * flim2_O2 * dtbiojour

           !$$-- Oxydation by NO3 --$$
           xtmp = effetchaleur * flim1_NO3 * glim2_O2 * dtbiojour
           F_remin_NO3 = p_N_remin * xtmp * Sfliminv ! Labile
           F_reminR_NO3 = p_N_reminR * flimz * xtmp * Sfliminv ! Refractaire

           !$$-- Denitrification: MO --> N2 using NO3 --$$
           F_denit = p_DNO3_denit * F_remin_NO3
           F_denitR = p_DNO3_denit * F_reminR_NO3

           !$$-- Dissimilatory Nitrate Reduction to Ammonium : MO --> NH4 using NO3 --$$ 
           F_dnra  = (1.0_rsh - p_DNO3_denit)  * F_remin_NO3
           F_dnraR = (1.0_rsh - p_DNO3_denit) * F_reminR_NO3

           !$$-- ODU oxydation and solidification --$$
           F_oxyd_ODU = p_ODU_oxy * flim3_O2 * dtbiojour
           F_solid_ODU = p_ODU_precip * dtbiojour

           !$$-- Anaerobic Remineralization --$$
           xtmp = effetchaleur * glim1_O2 * glim1_NO3 * dtbiojour
           F_remin_anaerN = p_N_remin * xtmp * Sfliminv
           F_reminR_anaerN = p_N_reminR *flimz * xtmp * Sfliminv

           !-------- EVOLUTION ----------
           dcdt(iv_oxygen) = dcdt(iv_oxygen) &
               + (-F_oxyd_ODU * cs(iv_ODU)   &
                  -F_remin_aerN * cs(iv_detr_N) * porosite_inv * p_GO2_NOrg    &
                  -F_reminR_aerN * cs(iv_detrR_N) * porosite_inv * p_GO2_NOrgR &
                  -F_nitrif * cs(iv_nutr_NH4) * p_GO2_NH4) &
                  * 0.032_rsh ! transformation en mg/L 

           dcdt(iv_detr_N) = dcdt(iv_detr_N) &
               -(F_remin_aerN + F_remin_anaerN + F_remin_NO3) * cs(iv_detr_N) 

           dcdt(iv_detrR_N) = dcdt(iv_detrR_N) &
               -(F_reminR_aerN + F_reminR_anaerN + F_reminR_NO3) * cs(iv_detrR_N)

           dcdt(iv_nutr_NH4) = dcdt(iv_nutr_NH4) &
               + (F_remin_aerN + F_remin_anaerN + F_dnra) * cs(iv_detr_N) * porosite_inv     &
               + (F_reminR_aerN + F_reminR_anaerN + F_dnraR) * cs(iv_detrR_N) * porosite_inv &
               - F_nitrif * cs(iv_nutr_NH4)

           dcdt(iv_nutr_NO3) = dcdt(iv_nutr_NO3) &
               + F_nitrif * cs(iv_nutr_NH4)      &
               - (F_denit + F_dnra) * cs(iv_detr_N) * p_GNO3_Norg * porosite_inv &
               - (F_denitR + F_dnraR) * cs(iv_detrR_N) * p_GNO3_NorgR * porosite_inv

           dcdt(iv_ODU) = dcdt(iv_ODU) &
               + F_remin_anaerN * cs(iv_detr_N) * p_GODU_NOrg * porosite_inv   &
               + F_reminR_anaerN * cs(iv_detrR_N) * p_GODU_NOrgR *porosite_inv &
               -(F_oxyd_ODU+F_solid_ODU)*cs(iv_ODU)

           !----------- Diagnostic -----------------
           diag_3d_sed(id_remin_aerN,k,i,j)=(F_remin_aerN * cs(iv_detr_N) + F_reminR_aerN * cs(iv_detrR_N)) * dzs(k,i,j)
           diag_3d_sed(id_remin_anaerN,k,i,j)=(F_remin_anaerN * cs(iv_detr_N) + F_reminR_anaerN * cs(iv_detrR_N)) * dzs(k,i,j)
           diag_3d_sed(id_remin_nitrateN,k,i,j)=(F_remin_NO3 * cs(iv_detr_N) + F_reminR_NO3 * cs(iv_detrR_N)) * dzs(k,i,j)
           diag_3d_sed(id_nitrif,k,i,j)=F_nitrif * cs(iv_nutr_NH4) * dzs(k,i,j)
           diag_3d_sed(id_remin_dnraN,k,i,j)=F_dnra * cs(iv_detr_N) * porosite_inv * dzs(k,i,j) &
               +F_dnraR * cs(iv_detrR_N) * porosite_inv * dzs(k,i,j)
           diag_3d_sed(id_remin_denitN,k,i,j)=F_denit * cs(iv_detr_N) * p_GNO3_Norg * porosite_inv * dzs(k,i,j)  &
               +F_denitR * cs(iv_detrR_N) * p_GNO3_NorgR * porosite_inv * dzs(k,i,j)
           diag_3d_sed(id_oxyd_solid_ODU,k,i,j)=(F_oxyd_ODU+F_solid_ODU)*cs(iv_ODU)*dzs(k,i,j)



        END SUBROUTINE nitrogen_cycle


        SUBROUTINE phosphorus_cycle(i,j,k,dtbiojour,cs,dcdt,zmiddle,flim1_O2,glim1_O2,glim2_O2,flim1_NO3,glim1_NO3,Sfliminv,effetchaleur,porosite_inv)

            INTEGER, intent(in)  :: i, j, k
            REAL(KIND=rsh), INTENT(in) :: dtbiojour
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: cs

            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(inout)  :: dcdt
            REAL(KIND=rsh), INTENT(in) :: zmiddle,flim1_O2,glim1_O2,glim2_O2,flim1_NO3,glim1_NO3,Sfliminv,effetchaleur,porosite_inv

            REAL(rsh) :: xtmp, flim4_O2, glim3_O2, flim2_NO3, glim2_NO3
            REAL(rsh) :: F_remin_aerP, F_reminR_aerP, F_precPFE_O2
            REAL(rsh) :: F_remin_NO3_P, F_reminR_NO3_P, F_precPFE_NO3
            REAL(rsh) :: F_remin_anaerP, F_reminR_anaerP, F_dissolPFE
            REAL(rsh) :: adsormax, F_adsorP, F_desorP, flimz


            flim4_O2 = 0.0_rsh
            glim3_O2 = 1.0_rsh
            IF(cs(iv_oxygen).gt.0.0_rsh) THEN
                flim4_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_precPFe)
                glim3_O2 = 1.0_rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_dissPFe)
            ENDIF

            flim2_NO3 = 0.0_rsh
            glim2_NO3 = 1.0_rsh
            IF(cs(iv_nutr_NO3).gt.0.0_rsh) THEN
                flim2_NO3 = cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kNO3_precPFe)
                glim2_NO3 = 1.0_rsh - cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kiNO3_dissPFe)
            ENDIF

            ! aerobiose
            xtmp = effetchaleur * flim1_O2 * dtbiojour
            flimz = max(p_xflimz,exp(-p_k_remin*zmiddle))
            F_remin_aerP = p_P_remin * xtmp * Sfliminv
            F_reminR_aerP = p_P_reminR * flimz * xtmp * Sfliminv
            F_precPFE_O2 = p_P_precFeO2 * flim4_O2 * dtbiojour

            ! sub-oxique (NO3 comme oxydant)
            xtmp = effetchaleur * flim1_NO3 * glim2_O2 * dtbiojour
            F_remin_NO3_P = p_P_remin * xtmp * Sfliminv
            F_reminR_NO3_P = p_P_reminR * flimz * xtmp * Sfliminv
            F_precPFE_NO3 = p_P_precFeNO3 * glim3_O2 * flim2_NO3 * dtbiojour

            ! anaerobiose
            xtmp = effetchaleur * glim1_O2 * glim1_NO3 * dtbiojour
            F_remin_anaerP = p_P_speedup_reminanaer * p_P_remin * xtmp * Sfliminv
            F_reminR_anaerP = p_P_speedup_reminanaer * p_P_reminR * flimz * xtmp * Sfliminv
            F_dissolPFE = p_P_dissFe * glim3_O2 * glim2_NO3 * flimz * dtbiojour

            ! adsorption / desorption sur mes
            adsormax = p_P_adsormaxsed * cs(iv_spim)
            F_adsorP = p_P_adsor * max(0.0_rsh, (adsormax - cs(iv_nutr_Pads))) * dtbiojour
            IF(cs(iv_oxygen).le.0.2_rsh) F_adsorP = F_adsorP/100
            
            F_desorP = 0.0_rsh
            IF(adsormax.gt.0.0_rsh) F_desorP = p_P_desor * (cs(iv_nutr_Pads) / adsormax) * dtbiojour

            !--------- EVOLUTION -----------------
            dcdt(iv_nutr_PO4) = dcdt(iv_nutr_PO4) &
                + (F_remin_aerP + F_remin_anaerP + F_remin_NO3_P) * cs(iv_detr_P) * porosite_inv     &
                + (F_reminR_aerP + F_reminR_anaerP + F_reminR_NO3_P) * cs(iv_detrR_P) * porosite_inv &
                + (F_desorP * cs(iv_nutr_Pads) * porosite_inv)                                       &
                + (F_dissolPFE * cs(iv_PFe) * porosite_inv)                                          &
                - (F_adsorP + F_precPFE_O2 + F_precPFE_NO3)*cs(iv_nutr_PO4)

            dcdt(iv_detr_P) = dcdt(iv_detr_P) &
                - (F_remin_aerP + F_remin_anaerP + F_remin_NO3_P) * cs(iv_detr_P)

            dcdt(iv_detrR_P) = dcdt(iv_detrR_P) &
                - (F_reminR_aerP + F_reminR_anaerP + F_reminR_NO3_P) * cs(iv_detrR_P)

            dcdt(iv_PFe) = dcdt(iv_PFe)  &
                + (F_precPFE_O2 + F_precPFE_NO3) * cs(iv_nutr_PO4) * poro(k,i,j) &
                - F_dissolPFE * cs(iv_PFe)

            dcdt(iv_nutr_Pads) = dcdt(iv_nutr_Pads) &
                + F_adsorP * cs(iv_nutr_PO4) * poro(k,i,j) &
                - F_desorP * cs(iv_nutr_Pads)

            dcdt(iv_oxygen)= dcdt(iv_oxygen)  &
                - F_precPFE_O2 * cs(iv_nutr_PO4) * p_GO2_PFe * 0.032_rsh

            dcdt(iv_nutr_NO3)=dcdt(iv_nutr_NO3) &
                - F_precPFE_NO3 * cs(iv_nutr_PO4) * p_GNO3_PFe

            !----------- Diagnostic -----------------
            diag_3d_sed(id_remin_aerP,k,i,j)=(F_remin_aerP * cs(iv_detr_P) + F_reminR_aerP * cs(iv_detrR_P)) * dzs(k,i,j)
            diag_3d_sed(id_remin_anaerP,k,i,j)=(F_remin_anaerP * cs(iv_detr_P) + F_reminR_anaerP * cs(iv_detrR_P)) * dzs(k,i,j)
            diag_3d_sed(id_remin_nitrateP,k,i,j)=(F_remin_NO3_P * cs(iv_detr_P) + F_reminR_NO3_P * cs(iv_detrR_P)) * dzs(k,i,j)
            diag_3d_sed(id_dissol_PFe,k,i,j)=(F_dissolPFE * cs(iv_PFe) * dzs(k,i,j))
            diag_3d_sed(id_precipit_P,k,i,j)=(F_precPFE_O2 + F_precPFE_NO3) * cs(iv_nutr_PO4) * poro(k,i,j) * dzs(k,i,j)
            diag_3d_sed(id_adsor_desorb_P,k,i,j)=dcdt(iv_nutr_Pads) * dzs(k,i,j)



        END SUBROUTINE phosphorus_cycle



        SUBROUTINE silicate_cycle(i,j,k,dtbiojour,cs,dcdt, zmiddle, temper,effetchaleur,porosite_inv)

            INTEGER, intent(in)  :: i, j, k
            REAL(KIND=rsh), INTENT(in) :: dtbiojour
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: cs

            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(inout)  :: dcdt
            REAL(KIND=rsh), INTENT(in) :: zmiddle, temper,effetchaleur,porosite_inv

            REAL(rsh) :: glim_silica, effetchaleurSi, F_remin_Si, F_precSi, diatmortsed

            glim_silica = 1-cs(iv_nutr_SiOH)/p_Si_Eq
            IF(glim_silica.lt.0._rsh) glim_silica=0._rsh
            effetchaleurSi = exp(p_T_effectSi*temper)
            F_remin_Si = ((p_BSi_dissSurfSed-p_BSi_dissFondSed) * exp(-p_kSi*zmiddle) + p_BSi_dissFondSed) * &
                effetchaleurSi * glim_silica *dtbiojour

            F_precSi = p_Si_precip *(cs(iv_nutr_SiOH)-p_Si_EqPrec) *dtbiojour
            IF(cs(iv_nutr_SiOH).le.p_Si_EqPrec) F_precSi = 0._rsh


            ! mortalite des diatomees
            diatmortsed=p_diat_mort_sed*effetchaleur*dtbiojour

            !--------- EVOLUTION --------

            dcdt(iv_detr_Si) = dcdt(iv_detr_Si) &
                 - F_remin_Si * cs(iv_detr_Si)  &
                 + diatmortsed * cs(iv_phyto_diat_N) * p_phyto_SiNratio

            dcdt(iv_nutr_SiOH) = dcdt(iv_nutr_SiOH) &
                + F_remin_Si * cs(iv_detr_Si) * porosite_inv &
                - F_precSi
            dcdt(iv_phyto_diat_N) = dcdt(iv_phyto_diat_N) &
                -diatmortsed*cs(iv_phyto_diat_N)

            dcdt(iv_detr_N) = dcdt(iv_detr_N) &
                + diatmortsed * cs(iv_phyto_diat_N)


            dcdt(iv_detr_P) = dcdt(iv_detr_P)                &
                + diatmortsed * cs(iv_phyto_diat_N) * rappaz

            !----------- Diagnostic -----------------
            diag_3d_sed(id_remin_aerSi,k,i,j)=F_remin_Si * cs(iv_detr_Si) * dzs(k,i,j)
            diag_3d_sed(id_precipit_Si,k,i,j)=F_precSi * dzs(k,i,j)
            diag_3d_sed(id_morta_phyto,k,i,j)=dcdt(iv_phyto_diat_N)*dzs(k,i,j)



        END SUBROUTINE silicate_cycle



        SUBROUTINE MO_aging(k,dtbiojour, cs, dcdt, ksmin)
            INTEGER, INTENT(in) :: k
            REAL(KIND=rsh), INTENT(in) :: dtbiojour, ksmin
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: cs
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(inout)  :: dcdt
            REAL(rsh) :: F_transf_MO, F_burried

            F_transf_MO = p_aging_MO * dtbiojour
            F_burried = 0.0_rsh
            IF(k==ksmin) F_burried = p_burried * dtbiojour

            dcdt(iv_detr_N) = dcdt(iv_detr_N) &
                -(F_transf_MO + F_burried) * cs(iv_detr_N) 

            dcdt(iv_detrR_N) = dcdt(iv_detrR_N)   &
                - F_burried   * cs(iv_detrR_N)    &
                + F_transf_MO * cs(iv_detr_N)

            dcdt(iv_detr_P) = dcdt(iv_detr_P)                &
                - (F_transf_MO + F_burried)  * cs(iv_detr_P)

            dcdt(iv_detrR_P) = dcdt(iv_detrR_P) &
                - F_burried * cs(iv_detrR_P)   &
                + F_transf_MO * cs(iv_detr_P)

            dcdt(iv_PFe) = dcdt(iv_PFe)  &
                - F_burried * cs(iv_PFe)

            dcdt(iv_nutr_Pads) = dcdt(iv_nutr_Pads) &
                - F_burried * cs(iv_nutr_Pads)

            dcdt(iv_detr_Si) = dcdt(iv_detr_Si) &
                - F_burried * cs(iv_detr_Si)



        END SUBROUTINE MO_aging

        SUBROUTINE update_state(i,j,k,dtbiojour,cs, dcdt, dcw_filtbent, zmiddle, ksmax)

            INTEGER, INTENT(in) :: i,j,k
            INTEGER :: iv
            REAL(KIND=rsh), INTENT(in) :: dtbiojour
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: cs
            REAL(KIND=rsh),DIMENSION(nv_state), INTENT(in)  :: dcdt, dcw_filtbent
            REAL(KIND=rsh), INTENT(in) :: zmiddle, ksmax

            REAL(rsh) :: cvolp, F_aeration,KO2_aeration,KzO2_aeration,o2sats,tempabs,gst


            ! update cv_sed and transfer from 1D to 3D
            cv_sed(1:nv_adv,k,i,j) = cv_sed(1:nv_adv,k,i,j)+  &
                dcdt(1:nv_adv)*unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))

            ! prevents dissolved_oxygen and nitrate from being negative
            IF(cv_sed(iv_oxygen,k,i,j)<0.0_rsh) cv_sed(iv_oxygen,k,i,j)=0.0_rsh
            IF(cv_sed(iv_nutr_NO3,k,i,j)<0.0_rsh) cv_sed(iv_nutr_NO3,k,i,j)=0.0_rsh

            ! introduction of rearation pour l oxygene
            ! cycle de l oxygene
            ! si pas d eau ==> reaeration par l air en surface
            ! traite apres coup pour tenir compte des fortes consommations

            IF(htot(i,j) < RESIDUAL_THICKNESS_WAT) THEN
                tempabs=cv_sed(-1,ksmax,i,j)+273.15_rsh
                gst=-173.4292_rsh+249.6329_rsh*100.0_rsh/tempabs+143.3483_rsh*log(tempabs/100.0_rsh) &
                    -21.8492_rsh*tempabs/100.0_rsh+cv_sed(0,ksmax,i,j)*(-0.033096_rsh+0.014259_rsh*tempabs/100.0_rsh &
                    -0.0017_rsh*(tempabs/100.0_rsh)**2)
                o2sats=1.429_rsh*exp(gst)
                KO2_aeration=p_KO2sed_aeration
                KzO2_aeration=p_Kzsed_aeration

                F_aeration=KO2_aeration*exp(-KzO2_aeration*zmiddle)*(o2sats-cv_sed(iv_oxygen,k,i,j))*dtbiojour
                IF(id_fluxsed_aeration .ne. 0) diag_3d_sed(id_fluxsed_aeration,k,i,j)=F_aeration*dzs(k,i,j)
                cv_sed(iv_oxygen,k,i,j)=cv_sed(iv_oxygen,k,i,j)+ F_aeration
            ENDIF

            ! if variable bio = var.part.constitutiv ==> update poro
            IF(SUM(unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))) .NE. nv_adv) THEN
                c_sedtot(k,i,j)=0.0_rsh
                cvolp=0.0_rsh
                DO iv=1,nvpc
                    c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                    cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
                ENDDO
                poro(k,i,j)=1.0_rsh-cvolp
                IF(poro(k,i,j) .LE. 0.0_rsh) THEN
                    write(*,*)'reactions_in_sed poro<0 !! ',CURRENT_TIME,i,j,k,poro(k,i,j),cv_sed(1:nvpc,k,i,j)
                ENDIF
            ENDIF

            cw_bottom_MUSTANG(1:nv_adv,i,j)=cw_bottom_MUSTANG(1:nv_adv,i,j)+  &
                dcw_filtbent(1:nv_adv)*dtbiojour*unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))

!TODO : Attention à vérifier ! Ce bloc est normalement soumis à une condition : IF(k==ksma(i,j) .AND. htot(i,j) >
!RESIDUAL_THICKNESS_WAT .AND. txfiltbenth .NE. 0.0_rsh . Cela n'aura
!aucun effet sur la mise à jour de cw_bottom_MUSTANG car dans le cas où elle
!n'est pas valide dcw_filtbent sera égal à 0. La question se pose uniquement
!pour la mise à jour de cv_wat

#ifdef key_MARS
            IF (htot(i,j) < hm) THEN
                cw_bottom_MUSTANG(1:nv_adv,i,j)=cw_bottom_MUSTANG(1:nv_adv,i,j) &
                    + dcw_filtbent(1:nv_adv)*dtbiojour/htot(i,j)*unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))
                DO iv=1,nv_adv
                    cv_wat(iv,:,i,j)=cw_bottom_MUSTANG(iv,i,j) 
                ENDDO
            ELSE
                cv_wat(1:nv_adv,1,i,j)=cw_bottom_MUSTANG(1:nv_adv,i,j)
            ENDIF
#endif



        END SUBROUTINE update_state



END MODULE bio_sediment
