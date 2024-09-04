  
#if defined key_oyster_DEB_GAMELAG

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_oysterDEB_GAMELAG  ***
   !&E
   !&E ** Purpose : Calcul des variations dues a la croissance des filtreurs
   !&E              selon le modèle DEB
   !&E              
    !&E       !          (Inspired by S. Petton DEB code) Original code
   !&E       
   !&E      use from general modele : CELL_SURF  
   !&E      use from general modele : 
   !&E      use from BIOLink  : ,dtbiojour, temper, extinction, PAR_top_layer, htot
   !&E                           cmes_3dmgl,diag_3d_wat(irk_diag(id_totalchl)
   !&E      use from basic bloom modele : c, epn,  effetchaleur, nbhuitre (in initdefine), fact_phyto_ChlNratio
   !&E                                    effetselnutdiat
   !&E      OUTPUT :   dc, diag_2d(irk_diag(id_oys_DEB....
   !&E---------------------------------------------------------------------
 
         IF (k.eq.2) THEN 
          
             IF(SUM(c) .ne. 0.0) THEN

              if(PO2_sat .lt. DOcrit) then
              cDO = 0.
              else
              cDO = min(1.,((100.*c(iv_oxygen)/O2_sat_deb)-DOcrit)/((100.*c(iv_oxygen)/O2_sat_deb)-DOcrit+bDO)   &
                                                                 /((100.-DOcrit)/(100.-DOcrit+bDO)))
              endif

              ! Temperature effect
              kTO = exp(Ta/T1-Ta/tempK) * &
              (1 + exp(Tal/T1 - Tal/Tl) + exp(Tah/Th - Tah/T1)) / &
              (1 + exp(Tal/tempK  - Tal/Tl) + exp(Tah/Th - Tah/tempK ))

              ! Food and Assimilation
#ifdef GAMELAG_EXACT
              Xtot_N = alpha_PS * c(iv_phyto_nano_N) + alpha_PL * c(iv_phyto_diat_N)*c(iv_phyto_diat_N)/(c(iv_phyto_diat_N)+0.01) &
              + alpha_ZS * c(iv_zoo_micr_N) + alpha_ZL * c(iv_zoo_meso_N)      &
              + alpha_PON * (c(iv_detr_N)*c(iv_detr_N) /(c(iv_detr_N) +0.0001) &
              + c(iv_detrR_N)*c(iv_detrR_N) /(c(iv_detrR_N) +0.0001)) ! mmolN/m3
#else
              Xtot_N = alpha_PS * c(iv_phyto_nano_N) + alpha_PL * c(iv_phyto_diat_N)  &
              + alpha_ZS * c(iv_zoo_micr_N) + alpha_ZL * c(iv_zoo_meso_N)      &
              + alpha_PON * (c(iv_detr_N) + c(iv_detrR_N)) ! mmolN/m3
#endif
#ifdef GAMELAG_EXACT
              Xtot_P = (alpha_PS * c(iv_phyto_nano_N) + alpha_PL *                  &
              c(iv_phyto_diat_N)*c(iv_phyto_diat_N)/(c(iv_phyto_diat_N)+0.01)) / p_phyto_NPratio  &
              + (alpha_ZS * c(iv_zoo_micr_N) + alpha_ZL * c(iv_zoo_meso_N)) / p_phyto_NPratio   &      ! in bloom same NP ratio for zoo and phyto
              + alpha_POP * (c(iv_detr_P)*c(iv_detr_P) /(c(iv_detr_P) +0.0001) &
              + c(iv_detrR_P)*c(iv_detrR_P) /(c(iv_detrR_P) +0.0001)) ! mmolN/m3              
#else
              Xtot_P = (alpha_PS * c(iv_phyto_nano_N) + alpha_PL * c(iv_phyto_diat_N)) / p_phyto_NPratio  &
              + (alpha_ZS * c(iv_zoo_micr_N) + alpha_ZL * c(iv_zoo_meso_N)) / p_phyto_NPratio   &      ! in bloom same NP ratio for zoo and phyto
              + alpha_POP * (c(iv_detr_P) + c(iv_detrR_P)) ! mmolP/m3
#endif
             ENDIF 

#ifdef key_oysterspat_DEB_GAMELAG
            DO ideb=1,2 
#else
            DO ideb=1,1
#endif 

               if (ideb.eq.1) then
                  iv_ndeb=iv_oysdeb
                  iv_E=iv_oysdeb_E
                  iv_E_V=iv_oysdeb_E_V
                  iv_E_GO=iv_oysdeb_E_GO
                  iv_E_R=iv_oysdeb_E_R
                  XkN=XkN_OYST
                  XkP=XkP_OYST
               else
                  iv_ndeb=iv_spatdeb
                  iv_E=iv_spatdeb_E
                  iv_E_V=iv_spatdeb_E_V
                  iv_E_GO=iv_spatdeb_E_GO
                  iv_E_R=iv_spatdeb_E_R
                  XkN=XkN_SPAT
                  XkP=XkP_SPAT
               endif

             IF(c(iv_ndeb) .ne. 0.0) THEN

              ! DEB model (adapted from S. Petton Fortran Script to match the GAMELAG DEB formulation)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!
              E=c(iv_E)
              E_V=c(iv_E_V)
              E_GO=c(iv_E_GO)
              E_R=c(iv_E_R)
              oyster_density = c(iv_ndeb) / (epn * CELL_SURF(i,j)) ! number of oyster per m3

              ! Volume
              V = E_V / (muV*d_v)  ! Here we divide by muV*d_V instead of Eg because the cost of growth is substracted from E_V and converted into NH4/PO4 release 
              V_2_3 = V**(2.0/3.0)
              ! Longueur
              L = V**(1.0/3.0)


              f_nut = min(Xtot_N / (Xtot_N + XkN), Xtot_P / (Xtot_P + XkP)) 

              ! Ingestion
              pX= cDO * Jxm * f_nut * V_2_3 * kTO * T_im

              !Assimilation
              pA =  KappaX * pX

              ! Metabolism : Mobilisation
              pC1 = (E / V )/(Eg + Kappa * (E/V)) * ((Eg * KappaX * Jxm  / Em * V_2_3) + Pm_deb  * V) * kTO

              ! Somatic Maintenance
              pMt = Pm_deb * V * kTO    

              pM1 = min(pMt, Kappa * pC1) ! fraction of somatic maintenance directly supported by reserve mobilization (J/d)  / In GAMELAG but Not in S. Petton Fortran DEB

              ! Growth
              pG = max(Kappa * pC1 - pM1 , 0.0)   

              Vp=(deltam * Lp)**3.0
              ! Maturity Maintenance
              pJ = min(Vp,V)*(1.-Kappa)/Kappa * Pm_deb * kTO

              ! Reproduction
              if (V .ge. Vp) then
              pR = max((1.-Kappa)*pC1 - pJ,0.0)
              pC2 = E_R* ((KappaX*Jxm*kTO)/(Em*V**(1.0/3.0))+Pm_deb*kTO/Eg)*(1.0-Kappa*E/(Eg*V+Kappa*E))
              else
              pR = 0.0
              pC2 = 0.0
              endif

              pM2 = max(min((pMt-pM1), pc2),0.)   ! Pm2>0, exists only if energy deficit available for structural maintenance

              ! Energy available for gamete production
              pGo = pC2 - pM2

              !! In S. Petton Fortran Deb 
              !! Shrinking if 'energy density' < length
              !! e = E / V / Em
              !! l = L / L_m
              !! pL2=0.0 ! Lysis of the gonade
              !! pL1=0.0 ! Extreme case : decrease of the structure
              !! pL1 and pL2 defined only with conditions (if (e<l))
              pL2 = max(min((pMt - pM1 - pM2) / Ygo, E_GO), 0.0) !Pl2>0: If structural maintenance deficit persists: gametes resorption energy intake
              pL1 = max((pMt - pM1 - pM2 - Ygo * pL2)/yv, 0.0) !Pl1>0: If still in deficit for structural maintenance, energy is drawn from the structure: structural lysis

              ! Dry Flesh Mass
              DFM = (E + E_R)/muE + d_v * V + dgo *E_Go/Eggo

              ! Storage
              dc(iv_E) = pA - pC1 
              ! Energy for structure
              dc(iv_E_V) = (muV * d_v /Eg) * pG  - pL1  

              ! Energy for reproduction
              dc(iv_E_R) = pR - pC2
              ! Energy allocated to gonades
              dc(iv_E_GO)= kappaGo * pGo - pL2 

              ! Spawning  
              ! Different from GAMELAG
              spaw = 0.0
              GSR = (dgo *E_Go/Eggo) / DFM
              if( (GSR .gt. GSR_spawn) .and. (tempK .gt. T_spawn) ) then
              spaw = E_Go
              c(iv_E_GO) = 0.
              dc(iv_E_GO) = 0.
              endif

              !!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! DEB to BGC formulation
              !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
              ! Food stoichiometry
              NP_X = Xtot_N / Xtot_P
              ! Feaces + Pseudo-feces
              Oyster_feces_N = ((1 - KappaX) * pX / C_oyster) * oyster_density ! mmolN/m3/d
              Oyster_feces_P = ((1 - KappaX) * pX / (C_oyster * NP_X)) * oyster_density ! mmolP/m3/d
              ! Nutrients regeneration
              ! Excess N regenerated to NH4
              if (NP_X > NP_oyster) then
              N_regen_oyster = ((pA / C_oyster) * (1 - NP_oyster * (1/NP_X))) ! mmolN/d
              else
              N_regen_oyster = 0
              endif  
              ! Excess P regenerated to PO4
              if (NP_X <= NP_oyster) then
              P_regen_oyster = (pA / C_oyster * ((1/NP_X) - (1/NP_oyster))) !mmolP/d
              else
              P_regen_oyster = 0
              endif


              ! Grazing

              G_oyster_PS = (((1/C_oyster) * alpha_PS * c(iv_phyto_nano_N) * pX) / Xtot_N) &
              * oyster_density ! mmolN/m3/d
#ifdef GAMELAG_EXACT
              G_oyster_PL = (((1/C_oyster) * alpha_PL * c(iv_phyto_diat_N) * pX) / Xtot_N) &
              *c(iv_phyto_diat_N)/(c(iv_phyto_diat_N)+0.01) * oyster_density ! mmolN/m3/d
#else
              G_oyster_PL = (((1/C_oyster) * alpha_PL * c(iv_phyto_diat_N) * pX) / Xtot_N) &
              * oyster_density ! mmolN/m3/d
#endif
              G_oyster_ZS = (((1/C_oyster) * alpha_ZS * c(iv_zoo_micr_N) * pX) / Xtot_N) &
              * oyster_density ! mmolN/m3/d
              G_oyster_ZL = (((1/C_oyster) * alpha_ZL * c(iv_zoo_meso_N) * pX) / Xtot_N) &
              * oyster_density ! mmolN/m3/d
#ifdef GAMELAG_EXACT              
              G_oyster_PON = (((1/C_oyster) * alpha_PON * (c(iv_detr_N)*c(iv_detr_N)/(c(iv_detr_N)+0.0001)      &
                        + c(iv_detrR_N)*c(iv_detrR_N)/(c(iv_detrR_N)+0.0001)) * pX) / Xtot_N)    &
                         * oyster_density ! mmolN/m3/d
  
              G_oyster_POP = (((1/C_oyster) * alpha_POP * (c(iv_detr_P)*c(iv_detr_P) /(c(iv_detr_P) +0.0001) &
              + c(iv_detrR_P)*c(iv_detrR_P) /(c(iv_detrR_P) +0.0001) ) * pX) / Xtot_P) &
              * oyster_density ! mmolP/m3/d
#else
              G_oyster_PON = (((1/C_oyster) * alpha_PON * (c(iv_detr_N) + c(iv_detrR_N)) * pX) / Xtot_N) &
              * oyster_density ! mmolN/m3/d
              G_oyster_POP = (((1/C_oyster) * alpha_POP * (c(iv_detr_P) + c(iv_detrR_P)) * pX) / Xtot_P) &
              * oyster_density ! mmolP/m3/d
#endif

              ! T. Guyondet comment in GAMELAG:  Modif connect OYST mortality with POM dynamics (23/10/2022)
              if (PO2_sat > DOcrit_MR) then
              NMOR_oyster =  oyster_mortality*(E + V * Eg + E_GO + E_R)/muE_N * oyster_density ! mmolN/m3/d
              PMOR_oyster =  oyster_mortality*(E + V * Eg + E_GO + E_R)/muE_N/NP_oyster * oyster_density ! mmolP/m3/d

              ! Rate of change of oyster number
              dc(iv_ndeb)=dc(iv_ndeb) - c(iv_ndeb) * oyster_mortality

              else
              NMOR_oyster =  OMR_DO*oyster_mortality*(E + V * Eg + E_GO + E_R)/muE_N * oyster_density ! mmolN/m3/d
              PMOR_oyster =  OMR_DO*oyster_mortality*(E + V * Eg + E_GO + E_R)/muE_N/NP_oyster * oyster_density ! mmolP/m3/d

              ! Rate of change of oyster number
              dc(iv_ndeb)=dc(iv_ndeb) - c(iv_ndeb) * OMR_DO * oyster_mortality

              endif

              NH4_oyster = (N_regen_oyster +                             &   
              (pMt + pJ) / muE_N +              &     ! cout de la maintenance, energie dissipée transformée en mmolN
              (pG * (Eg - d_v * muV) / Eg) / muE_N +  &   ! cout de la croissance
              (pGO * (Eggo - dgo * muGo) / Eggo) / muE_N )  &  ! cout fabrication des oeufs (1-KappaGo) ! why dont we compute KappaGO=dgo*muGo/Eggo
              * oyster_density ! mmolP/m3/d

              PO4_oyster = (P_regen_oyster +                             &     
              (pMt + pJ) / muE_N / NP_oyster +            &       ! cout de la maintenance
              (pG * (Eg - d_v * muV) / Eg) / muE_N / NP_oyster +  &   ! cout de la croissance
              (pGO * (Eggo - dgo * muGo) / Eggo) / muE_N / NP_oyster) & ! cout fabrication des oeufs
              * oyster_density ! mmolP/m3/d

              ! T. Guyondet comment in GAMELAG: Flux becomes null at low DO
              if (c(iv_oxygen) >= p_O2_Threshold) then
              respOyst = p_R_photo*NH4_oyster  ! gO2/m3/d
              else
              respOyst = 0.
              endif


              ! Rate of change of ammonium
              ! ----------------------------------
              dc(iv_nutr_NH4)=dc(iv_nutr_NH4) + NH4_oyster
              ! Rate of change of phosphate
              ! ----------------------------------
              dc(iv_nutr_PO4)=dc(iv_nutr_PO4) + PO4_oyster

              ! Rate of change of oxygen
              ! --------------------------------------
              dc(iv_oxygen)=dc(iv_oxygen) - respOyst


              ! Rate of change of dissolved organic nitrogen 
              ! --------------------------------------

              dc(iv_diss_detr_N)=dc(iv_diss_detr_N)    &
              + p_labi*(1.-epsOyst)*(1.-OMR_Dpo)*NMOR_oyster

              dc(iv_diss_detrR_N)=dc(iv_diss_detrR_N)    &
              + (1.-p_labi)*(1.-epsOyst)*(1.-OMR_Dpo)*NMOR_oyster

              ! Rate of change of dissolved organic phosphorus 
              ! --------------------------------------
              dc(iv_diss_detr_P)=dc(iv_diss_detr_P)    &
              + p_labi*(1.-epsOyst)*(1.-OMR_Dpo)*PMOR_oyster

              dc(iv_diss_detrR_P)=dc(iv_diss_detrR_P)    &
              + (1.-p_labi)*(1.-epsOyst)*(1.-OMR_Dpo)*PMOR_oyster


              ! Rate of change of particulate organic nitrogen 
              ! --------------------------------------
              dc(iv_detr_N)=dc(iv_detr_N)    &
              + p_labi*((1.-BioDpo)*Oyster_feces_N - G_oyster_PON + epsOyst*(1.-OMR_Dpo)*NMOR_oyster)

              dc(iv_detrR_N)=dc(iv_detrR_N)    &
              + (1.-p_labi)*((1.-BioDpo)*Oyster_feces_N - G_oyster_PON + epsOyst*(1.-OMR_Dpo)*NMOR_oyster)

              ! Rate of change of particulate organic phosphorus 
              ! --------------------------------------

              dc(iv_detr_P)=dc(iv_detr_P)    &
              + p_labi*((1.-BioDpo)*Oyster_feces_P - G_oyster_POP + epsOyst*(1.-OMR_Dpo)*PMOR_oyster)

              dc(iv_detrR_P)=dc(iv_detrR_P)    &
              + (1.-p_labi)*((1.-BioDpo)*Oyster_feces_P - G_oyster_POP + epsOyst*(1.-OMR_Dpo)*PMOR_oyster)


              ! Rate of change of small phyto
              ! -----------------------------------------------------------------
              dc(iv_phyto_nano_N)=dc(iv_phyto_nano_N) - G_oyster_PS

              ! Rate of change of large phyto
              ! -----------------------------------------------------------------
              dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N) - G_oyster_PL

              ! Rate of change of small zoo
              ! -----------------------------------------------------------------
              dc(iv_zoo_micr_N)=dc(iv_zoo_micr_N) - G_oyster_ZS

              ! Rate of change of large zoo
              ! -----------------------------------------------------------------
              dc(iv_zoo_meso_N)=dc(iv_zoo_meso_N) - G_oyster_ZL

 
              ! Direct fluxes to sediment
              flx_w2s_sum_CROCO(i,j,iv_detr_N+itsubs1-1)=flx_w2s_sum_CROCO(i,j,iv_detr_N+itsubs1-1)         &
                                +(BioDpo*Oyster_feces_N+OMR_Dpo*NMOR_oyster)*dtbiojour*thicklayerW_C(k,i,j)

              flx_w2s_sum_CROCO(i,j,iv_detr_P+itsubs1-1)=flx_w2s_sum_CROCO(i,j,iv_detr_P+itsubs1-1)          &
#ifdef GAMELAG_EXACT            
                                +(BioDpo*Oyster_feces_P)*dtbiojour*thicklayerW_C(k,i,j)
#else
                                +(BioDpo*Oyster_feces_P+OMR_Dpo*PMOR_oyster)*dtbiojour*thicklayerW_C(k,i,j)
#endif                                

#    

             ELSE  ! nbhuitre > 0 
              
              c(iv_E)=0.0_rsh
              c(iv_E_V)=0.0_rsh
              c(iv_E_GO)=0.0_rsh
              c(iv_E_R)=0.0_rsh

            ENDIF


             

            ENDDO
       ENDIF 
   !!======================================================================
#endif
