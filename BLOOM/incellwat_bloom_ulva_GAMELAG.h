#ifdef key_ulva_GAMELAG

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_ulva_GAMELAG  ***
   !&E
   !&E ** Purpose : Quota model for ulva
   !&E
   !&E       !  2024-09    (F. Dufois ) Adapted from GAMELAG model written in R

   !&E---------------------------------------------------------------------

          
        IF(c(iv_ulva) .gt. 0.0) THEN

              etaD = (1.-etaP)
              phiD = (1.-(phiP+phiPO))

              ! Ulva accesory functions
              f_I_Ulva = 1.-exp(-(I_phyto/Iopt_Ulva)) ! Adim
              f_Temp_Ulva = 1./(1.+exp(-zeta1_Ulva*(temper-Temp_Ulva))) ! Adim
        
              ! Ulva nutrient quota
              f_QN_Ulva = max((c(iv_ulva_N)- Qmin_N_Ulva)/(Qmax_N_Ulva-Qmin_N_Ulva),0.0_rsh) ! Adim
              f_QP_Ulva = max((c(iv_ulva_P) - Qmin_P_Ulva)/(Qmax_P_Ulva-Qmin_P_Ulva),0.0_rsh) ! Adim
              
              ! Limitation nutrients
              f_NO3_Ulva = c(iv_nutr_NO3)/(c(iv_nutr_NO3)+K_NO3_Ulva) ! Adim
              f_NH4_Ulva = c(iv_nutr_NH4)/(c(iv_nutr_NH4)+K_NH4_Ulva) ! Adim
              f_N_Ulva = f_NO3_Ulva + f_NH4_Ulva ! Adim
              f_PO4_Ulva = c(iv_nutr_PO4)/(c(iv_nutr_PO4)+K_PO4_Ulva) ! Adim
              
              ! Ulva uptake of N & P
              NO3_uptake_Ulva = V_NO3_Ulva*f_NO3_Ulva * (1. - f_QN_Ulva) ! mmolN/gdw/d
              NH4_uptake_Ulva = V_NH4_Ulva*f_NH4_Ulva * (1. - f_QN_Ulva) ! mmolN/gdw/d
              PO4_uptake_Ulva = V_PO4_Ulva*f_PO4_Ulva * (1. - f_QP_Ulva) ! mmolP/gdw/d
              
              ! Lim function for nutrients
              f_SN_Ulva = min(f_QN_Ulva, f_QP_Ulva) ! Adim
              
              mu_Ulva = mu_max_Ulva * f_I_Ulva * f_Temp_Ulva * f_SN_Ulva ! d-1
              
              ! Ulva growth
              growth_Ulva = mu_Ulva * c(iv_ulva)! gdw/m3/d
                      
              ! Ulva death (from Zaldivar et al 2003)

              death_Ulva = MaxMort_Ulva * (1./(1.+exp(-zeta2_Ulva*(temper-MaxTemp_Ulva)))) * c(iv_ulva) ! gdw/m3/d
              
              photos_Ulva = p_Phi_photo*growth_Ulva ! in gO2/m3/d

#ifdef GAMELAG_EXACT
              respUlva = p_R_photo*(RD_Ulva*death_Ulva+RG_Ulva*growth_Ulva)*c(iv_ulva_N)*c(iv_oxygen)/(c(iv_oxygen)+0.001)  ! in gO2/m3/d
#else
              respUlva = p_R_photo*(RD_Ulva*death_Ulva+RG_Ulva*growth_Ulva)*c(iv_ulva_N)  ! in gO2/m3/d
#endif               

              ! Rate of change of nitrate
              ! --------------------
              dc(iv_nutr_NO3)=dc(iv_nutr_NO3) - NO3_uptake_Ulva * c(iv_ulva)
              ! Rate of change of ammonium
              ! ----------------------------------
              dc(iv_nutr_NH4)=dc(iv_nutr_NH4) - NH4_uptake_Ulva * c(iv_ulva)

              ! Rate of change of phosphate
              ! ----------------------------------
              dc(iv_nutr_PO4)=dc(iv_nutr_PO4) + phiPO*death_Ulva* c(iv_ulva_P) - PO4_uptake_Ulva * c(iv_ulva)

              ! Rate of change of oxygen
              ! --------------------------------------
              dc(iv_oxygen)= dc(iv_oxygen) + photos_Ulva - respUlva
              
              ! Rate of change of dissolved organic nitrogen 
              ! --------------------------------------

              dc(iv_diss_detr_N)=dc(iv_diss_detr_N) + p_labi * etaD * death_Ulva * c(iv_ulva_N)

              dc(iv_diss_detrR_N)=dc(iv_diss_detrR_N) + (1.-p_labi) * etaD * death_Ulva * c(iv_ulva_N)

              ! Rate of change of dissolved organic phosphorus 
              ! --------------------------------------
              dc(iv_diss_detr_P)=dc(iv_diss_detr_P) + p_labi * phiD * death_Ulva * c(iv_ulva_P)

              dc(iv_diss_detrR_P)=dc(iv_diss_detrR_P) + (1.-p_labi) * phiD * death_Ulva * c(iv_ulva_P)

              ! Rate of change of particulate organic nitrogen 
              ! --------------------------------------
              
              dc(iv_detr_N)=dc(iv_detr_N)  + p_labi * etaP * death_Ulva * c(iv_ulva_N)

              dc(iv_detrR_N)=dc(iv_detrR_N)  + (1.-p_labi) * etaP * death_Ulva * c(iv_ulva_N)

              ! Rate of change of particulate organic phosphorus 
              ! --------------------------------------

              dc(iv_detr_P)=dc(iv_detr_P)  + p_labi * phiP * death_Ulva * c(iv_ulva_P)

              dc(iv_detrR_P)=dc(iv_detrR_P) + (1.-p_labi) * phiP * death_Ulva * c(iv_ulva_P)

              ! Rate of change of ulva biomass and N/P quotas 
              ! --------------------------------------              
              dc(iv_ulva)= growth_Ulva - death_Ulva ! in gdw/m3/d

              dc(iv_ulva_N)=NO3_uptake_Ulva + NH4_uptake_Ulva - mu_Ulva*c(iv_ulva_N)  ! in mmolN/gdw/d

              dc(iv_ulva_P)=PO4_uptake_Ulva - mu_Ulva* c(iv_ulva_P)                   ! in mmolP/gdw/d

        ENDIF

#endif
