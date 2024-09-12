#ifdef key_gracilaria_GAMELAG

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_gracilaria_GAMELAG  ***
   !&E
   !&E ** Purpose : Quota model for gracilaria
   !&E
   !&E       !  2024-09    (F. Dufois ) Adapted from GAMELAG model written in R

   !&E---------------------------------------------------------------------

          
        IF(c(iv_graci) .gt. 0.0) THEN

              omegaD = (1.-omegaP)
              thetaD = (1.-(thetaP+thetaPO))

              ! Gracilaria accesory functions
              f_I_Graci = 1.-exp(-(I_phyto/Iopt_Graci)) ! Adim
              f_Temp_Graci = 1./(1.+exp(-zeta1_Graci*(temper-Temp_Graci))) ! Adim
        
              ! Graci nutrient quota
              f_QN_Graci = max((c(iv_Graci_N)- Qmin_N_Graci)/(Qmax_N_Graci-Qmin_N_Graci),0.0_rsh) ! Adim
              f_QP_Graci = max((c(iv_Graci_P) - Qmin_P_Graci)/(Qmax_P_Graci-Qmin_P_Graci),0.0_rsh) ! Adim
              
              ! Limitation nutrients
              f_NO3_Graci = c(iv_nutr_NO3)/(c(iv_nutr_NO3)+K_NO3_Graci) ! Adim
              f_NH4_Graci = c(iv_nutr_NH4)/(c(iv_nutr_NH4)+K_NH4_Graci) ! Adim
              f_N_Graci = f_NO3_Graci + f_NH4_Graci ! Adim
              f_PO4_Graci = c(iv_nutr_PO4)/(c(iv_nutr_PO4)+K_PO4_Graci) ! Adim
              
              ! Graci uptake of N & P
              NO3_uptake_Graci = V_NO3_Graci*f_NO3_Graci * (1. - f_QN_Graci) ! mmolN/gdw/d
              NH4_uptake_Graci = V_NH4_Graci*f_NH4_Graci * (1. - f_QN_Graci) ! mmolN/gdw/d
              PO4_uptake_Graci = V_PO4_Graci*f_PO4_Graci * (1. - f_QP_Graci) ! mmolP/gdw/d
              
              ! Lim function for nutrients
              f_SN_Graci = min(f_QN_Graci, f_QP_Graci) ! Adim
              
              mu_Graci = mu_max_Graci * f_I_Graci * f_Temp_Graci * f_SN_Graci ! d-1
              
              ! Graci growth
              growth_Graci = mu_Graci * c(iv_Graci)! gdw/m3/d
                      
              ! Graci death (from Zaldivar et al 2003)

              death_Graci = MaxMort_Graci * (1./(1.+exp(-zeta2_Graci*(temper-MaxTemp_Graci)))) * c(iv_Graci) ! gdw/m3/d
              
              photos_Graci = p_Phi_photo*growth_Graci ! in gO2/m3/d

#ifdef GAMELAG_EXACT
              respGraci = p_R_photo*(RD_Graci*death_Graci+RG_Graci*growth_Graci)*c(iv_Graci_N)*c(iv_oxygen)/(c(iv_oxygen)+0.001)  ! in gO2/m3/d
#else
              respGraci = p_R_photo*(RD_Graci*death_Graci+RG_Graci*growth_Graci)*c(iv_Graci_N)  ! in gO2/m3/d
#endif               

              ! Rate of change of nitrate
              ! --------------------
              dc(iv_nutr_NO3)=dc(iv_nutr_NO3) - NO3_uptake_Graci * c(iv_Graci)
              ! Rate of change of ammonium
              ! ----------------------------------
              dc(iv_nutr_NH4)=dc(iv_nutr_NH4) - NH4_uptake_Graci * c(iv_Graci)

              ! Rate of change of phosphate
              ! ----------------------------------
              dc(iv_nutr_PO4)=dc(iv_nutr_PO4) + thetaPO*death_Graci* c(iv_Graci_P) - PO4_uptake_Graci * c(iv_Graci)

              ! Rate of change of oxygen
              ! --------------------------------------
              dc(iv_oxygen)= dc(iv_oxygen) + photos_Graci - respGraci
              
              ! Rate of change of dissolved organic nitrogen 
              ! --------------------------------------

              dc(iv_diss_detr_N)=dc(iv_diss_detr_N) + p_labi * omegaD * death_Graci * c(iv_Graci_N)

              dc(iv_diss_detrR_N)=dc(iv_diss_detrR_N) + (1.-p_labi) * omegaD * death_Graci * c(iv_Graci_N)

              ! Rate of change of dissolved organic phosphorus 
              ! --------------------------------------
              dc(iv_diss_detr_P)=dc(iv_diss_detr_P) + p_labi * thetaD * death_Graci * c(iv_Graci_P)

              dc(iv_diss_detrR_P)=dc(iv_diss_detrR_P) + (1.-p_labi) * thetaD * death_Graci * c(iv_Graci_P)

              ! Rate of change of particulate organic nitrogen 
              ! --------------------------------------
              
              dc(iv_detr_N)=dc(iv_detr_N)  + p_labi * omegaP * death_Graci * c(iv_Graci_N)

              dc(iv_detrR_N)=dc(iv_detrR_N)  + (1.-p_labi) * omegaP * death_Graci * c(iv_Graci_N)

              ! Rate of change of particulate organic phosphorus 
              ! --------------------------------------

              dc(iv_detr_P)=dc(iv_detr_P)  + p_labi * thetaP * death_Graci * c(iv_Graci_P)

              dc(iv_detrR_P)=dc(iv_detrR_P) + (1.-p_labi) * thetaP * death_Graci * c(iv_Graci_P)

              ! Rate of change of Graci biomass and N/P quotas 
              ! --------------------------------------              
              dc(iv_Graci)= growth_Graci - death_Graci ! in gdw/m3/d

              dc(iv_Graci_N)=NO3_uptake_Graci + NH4_uptake_Graci - mu_Graci*c(iv_Graci_N)  ! in mmolN/gdw/d

              dc(iv_Graci_P)=PO4_uptake_Graci - mu_Graci* c(iv_Graci_P)                   ! in mmolP/gdw/d

        ENDIF

#endif
