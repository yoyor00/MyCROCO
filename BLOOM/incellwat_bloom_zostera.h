#ifdef key_zostera


   !&E-------------------------------------------------------------------
   !&E                 ***   incellwat_zostera  ***
   !&E
   !&E ** Purpose : Calcul des variables concernant les herbiers de zosteres
   !&E              
   !&E       !  2013-09    (M. Plus)
   !&E       !
   !&E
   !&E      use from general modele : temper,  
   !&E      use from BIOLink  : PAR_top_layer,jjulien_BIOLink, dtbiojour
   !&E                           
   !&E      use from basic bloom modele : c, epn,epnsed, ibbenth
   !&E      OUTPUT :  diag_2d(irk_diag(id_zost_prod)),diag_3d_wat( et id_zost...) + dc
   !&E---------------------------------------------------------------------
   !! * Modules used

          ! quota en N et P des graines
          ! ---------------------------
          SNquota = 0.5_rsh*(p_zost_LNquotamax-p_zost_LNquotamin)
          SPquota = 0.5_rsh*(p_zost_LPquotamax-p_zost_LPquotamin) 
 
          IF (ibbenth.eq.1) THEN
 
            IF(c(iv_zost_LB).gt.1.e-10 .and. c(iv_zost_RB).gt.1e-10) THEN

               ! Effets de la temperature sur les zosteres (sd)
               ! ----------------------------------------------
               effetchaleurzprod = exp(p_zost_T_prod*temper)
               effetchaleurzrespf = exp(p_zost_T_respf*temper)
               effetchaleurzrespr = exp(p_zost_T_respr*temper)
               effetchaleurzmort = exp(p_zost_T_mort*(temper-20.0_rsh))
               effetchaleurzrecr = exp(p_zost_T_recr*(temper-5.0_rsh))

               ! Effet de la lumiere ("big leaf" model, Tchamitchian 1993)
               ! ---------------------------------------------------------
               LAI = p_zost_klai*c(iv_zost_LB)
               Qcan = p_parradratio*PAR_top_layer(k,i,j)
               !effetlumzost = tanh(Qcan*(1-exp(-p_zost_leafabscoef*LAI))/p_zost_Ik)
               effetlumzost = tanh(Qcan/p_zost_Ik)
     
               alumfond = Qcan*exp(-p_zost_leafabscoef*LAI)

               ! Effet des sels nutritifs sur l absorption (sd)
               ! ----------------------------------------------
               LNquota = c(iv_zost_LN)/c(iv_zost_LB)
               LNsat = MAX((LNquota - p_zost_LNquotamin)/(p_zost_LNquotamax - p_zost_LNquotamin),0.0_rsh)
               IF(LNsat.lt.0.99_rsh .AND. c(iv_nutr_NO3).gt.0.5_rsh) THEN 
                  limabsLno3 = c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_zost_KNO3L)*(1-LNsat)
                ELSE
                  limabsLno3 = 0.0_rsh 
                ENDIF
                IF(LNsat.lt.0.99_rsh .AND. c(iv_nutr_NH4).gt.0.1_rsh) THEN 
                  limabsLnh4 = c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_zost_KNH4L)*(1-LNsat)
                ELSE
                  limabsLnh4 = 0.0_rsh 
                ENDIF
                LPquota = c(iv_zost_LP)/c(iv_zost_LB)
                LPsat = MAX((LPquota - p_zost_LPquotamin)/(p_zost_LPquotamax - p_zost_LPquotamin),0.0_rsh)
                IF(LPsat.lt.0.99_rsh .AND. c(iv_nutr_PO4).gt.0.01_rsh) THEN
                  limabsLpo4 = c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_zost_KPL)*(1-LPsat)
                ELSE
                  limabsLpo4 = 0.0_rsh
                ENDIF      
                RNquota = c(iv_zost_RN)/c(iv_zost_RB)
                RNsat = MAX((RNquota - p_zost_RNquotamin)/(p_zost_RNquotamax - p_zost_RNquotamin),0.0_rsh)
                RPquota = c(iv_zost_RP)/c(iv_zost_RB)
                RPsat = MAX((RPquota - p_zost_RPquotamin)/(p_zost_RPquotamax - p_zost_RPquotamin),0.0_rsh)
                IF(RNsat.lt.0.99_rsh .AND. c(iv_benth_NH4).gt.1.0_rsh) THEN 
                  limabsRnh4 = c(iv_benth_NH4)/(c(iv_benth_NH4)+p_zost_KNH4R)*(1-RNsat)**p_zost_delta1
                ELSE
                  limabsRnh4 = 0.0_rsh
                ENDIF
                IF(RPsat.lt.0.99_rsh .AND. c(iv_benth_PO4).gt.0.1_rsh) THEN
                  limabsRpo4 = c(iv_benth_PO4)/(c(iv_benth_PO4)+p_zost_KPR)*(1-RPsat)**p_zost_delta1
                ELSE
                  limabsRpo4 = 0.0_rsh
                ENDIF

                ! Effet des quotas en N et P sur la production (sd)
                ! -------------------------------------------------
                limprodLN = LNsat**p_zost_delta2
                limprodLP = LPsat**p_zost_delta2

                ! Photosynthese (gO2/mmolC/j) et respiration (gO2/mmolC/j)
                ! ----------------------------------------------------------
                Gross_Prod_zost = p_zost_maxprod0*effetchaleurzprod*effetlumzost*MIN(limprodLN,limprodLP)

                Leaf_Resp_zost = p_zost_respf0*effetchaleurzrespf
                Root_Resp_zost = p_zost_respr0*effetchaleurzrespr

                ! Croissance des feuilles, Leaf_growth, et des rhizomes, Root_growth (j-1)
                ! --------------------------------------------------------
                Net_Prod_zost = Gross_Prod_zost - (Leaf_Resp_zost + Root_Resp_zost)
                Net_Prod_zost = Net_Prod_zost*1000/(32*p_zost_qphotos)   ! transformation de g02/mmolC/j en mmolC/mmolC/j
                Kzost = p_zost_K
                IF(Net_Prod_zost.gt.0.0_rsh) THEN
                  Root_growth = Net_Prod_zost*Kzost
                  Leaf_growth = Net_Prod_zost - Root_growth
                ELSE
                  Root_growth = 0.0_rsh
                  Leaf_growth = 0.0_rsh
                ENDIF

                ! Absorption des sels nutritifs (molN/molC/j)
                ! -------------------------------------------
                Labs_nh4 = p_zost_LNnh4maxabs*limabsLnh4
                Labs_no3 = p_zost_LNno3maxabs*limabsLno3
                Labs_po4 = p_zost_LPmaxabs*limabsLpo4*effetlumzost
      
                Rabs_nh4 = p_zost_RNmaxabs*limabsRnh4
                Rabs_po4 = p_zost_RPmaxabs*limabsRpo4
 
                ! transfert de N et P dans la plante (molN ou P/molC/j). Positif sens acropetal
                ! -----------------------------------------------------------------------------
                deltasat_N = RNsat - LNsat
                IF(deltasat_N .gt. 0._rsh .and. LNsat.lt.0.95_rsh) THEN 
                  trans_N_acro = p_zost_transfrate*deltasat_N*c(iv_zost_RN)
                  trans_N_basi = 0.0_rsh
                ELSE IF(deltasat_N .lt. 0._rsh .and. RNsat.lt.0.95_rsh) THEN
                  trans_N_basi = p_zost_transfrate*ABS(deltasat_N)*c(iv_zost_LN)
                  trans_N_acro = 0.0_rsh
                ELSE      
                  trans_N_acro = 0.0_rsh
                  trans_N_basi = 0.0_rsh
                END IF  
         
                deltasat_P = RPsat - LPsat
                IF(deltasat_P .gt. 0._rsh .and. LPsat.lt.0.95_rsh) THEN 
                  trans_P_acro = p_zost_transfrate*deltasat_P*c(iv_zost_RP)
                  trans_P_basi = 0.0_rsh
                ELSE IF(deltasat_P .lt. 0._rsh .and. RPsat.lt.0.95_rsh) THEN
                  trans_P_basi = p_zost_transfrate*ABS(deltasat_P)*c(iv_zost_LP)
                  trans_P_acro = 0.0_rsh
                ELSE
                  trans_P_acro = 0.0_rsh
                  trans_P_basi = 0.0_rsh
                END IF  

                ! conservation de N et P (molN ou P/molC/j) liee au processus de senescence
                ! -------------------------------------------------------------------------
                LNrecla = MAX(p_zost_reclamax*(1-(LNsat)**2),0.0_rsh)
                RNrecla = MAX(p_zost_reclamax*(1-(RNsat)**2),0.0_rsh)
                LPrecla = MAX(p_zost_reclamax*(1-(LPsat)**2),0.0_rsh)
                RPrecla = MAX(p_zost_reclamax*(1-(RPsat)**2),0.0_rsh)
 
                ! Mortalite des feuilles et des rhizomes (j-1)
                ! --------------------------------------------
                IF(WIND_SPEED(i,j).gt.6.0_rsh) THEN
                  effet_vent = (WIND_SPEED(i,j)-6.0_rsh)/10._rsh*exp(-1.2_rsh*TOTAL_WATER_HEIGHT(i,j))*exp(0.1_rsh*LAI)
                ELSE 
                  effet_vent = 0.0_rsh
                ENDIF
                Leaf_mort = p_zost_mort*effetchaleurzmort
                IF(c(iv_zost_LB).gt.1.e-7) THEN
                  Root_mort = Leaf_mort*p_zost_K
                ELSE
                  Root_mort = Leaf_mort*10.0_rsh
                ENDIF
      
                ! Recrutement de nouveaux pieds par reproduction vegetative (j-1) 
                ! ---------------------------------------------------------------
                limselfshad = alumfond/(alumfond+p_zost_KREC1)
                IF(c(iv_zost_RB).gt.10.0_rsh) THEN
                  limRB = c(iv_zost_RB)/(c(iv_zost_RB)+p_zost_KREC2)
                ELSE
                  limRB = 0.0_rsh
                ENDIF
                limRN = RNsat**p_zost_delta2
                limRP = RPsat**p_zost_delta2
                limnutrecr = MIN(limRN,limRP)
                Recruit_rate = p_zost_RECRmax*effetchaleurzrecr*MIN(limselfshad,limRB,limnutrecr)

                ! Effort de reproduction sexuee (s.d.) et production de graines (mmolC/pied) 
                ! --------------------------------------------------------------------------
                IF(c(iv_zost_LB).gt.200.0_rsh) THEN
                  revolangle = 0.2163108_rsh+2*atan(0.9671396_rsh*tan(0.0086_rsh*(jjulien_BIOLink-186))) !Forsythe et al. (1995)
                  declinangle = asin(0.39795_rsh*cos(revolangle))
                  dureej = 24-24/NUMBER_PI*acos((sin(0.833_rsh*NUMBER_PI/180)+sin(44.6667_rsh*NUMBER_PI/180)*sin(declinangle)) &
                         /(cos(44.6667_rsh*NUMBER_PI/180)*cos(declinangle)))
                  limdureej = MAX(dureej-14,0.0_rsh)/1.57_rsh
                  flowering_rate = p_zost_ERSmax*limdureej*effetchaleurzprod
                ELSE
                  flowering_rate = 0.0_rsh
                ENDIF

                ! tests sur quelques variables
                ! ----------------------------
!               ! IF(LNquota.lt.p_zost_LNquotamin .OR. LNquota.gt.p_zost_LNquotamax) THEN
!               !   print*,'Attention LNquota sort des bornes !!!! en',i,j,', LNquota =',LNquota
!               !   print*,'Labs =',(Labs_nh4 + Labs_no3)*c(iv_zost_LB)
!               !   print*,'Recruit_rate =',Recruit_rate*p_zost_SB1*c(iv_zost_D)*c(iv_zost_RN)/c(iv_zost_RB)
!               !   print*,'trans_N_acro, trans_N_basi =',trans_N_acro,trans_N_basi
!               !   print*,'LNrecla =',LNrecla
!               !   print*,'deltasat_N, RNsat, LNsat =',deltasat_N,RNsat,LNsat
!               !   print*,'-(Leaf_mort+effet_vent) =',-(Leaf_mort+effet_vent)*c(iv_zost_LN)
!               !   print*,'Germin_rate=',0.5_rsh*Germin_rate*c(iv_zost_benth_seed)*epnsed*SNquota 
                   IF(LNquota.lt.0.005_rsh .OR. LNquota.gt.0.15_rsh) THEN
                     print*,'Attention LNquota sort des bornes !!!! en',i,j,', LNquota =',LNquota
                     stop
                   ENDIF
!               ! ENDIF
!               ! IF(RNquota.lt.p_zost_RNquotamin .OR. RNquota.gt.p_zost_RNquotamax) THEN
!               !   print*,'Attention RNquota sort des bornes !!!! en ',i,j,', RNquota = ',RNquota
!               !   print*,'limabsRnh4 =',limabsRnh4
!               !   print*,'Recruit_rate =',Recruit_rate*p_zost_SB1*c(iv_zost_D)*c(iv_zost_RN)/c(iv_zost_RB)
!               !   print*,'trans_N_acro, trans_N_basi =',trans_N_acro,trans_N_basi
!               !   print*,'RNrecla =',RNrecla
!               !   print*,'deltasat_N, RNsat, LNsat =',deltasat_N,RNsat,LNsat
!               !  print*,'-Root_mort =',-Root_mort*c(iv_zost_RN)
!               !   print*,'Germin_rate=',0.5_rsh*Germin_rate*c(iv_zost_benth_seed)*epnsed*SNquota
!               !   print*,'LB=',c(iv_zost_LB),'RB=',c(iv_zost_RB)
!               !   print*,'LN=',c(iv_zost_LN),'RN=',c(iv_zost_RN)
                  IF(RNquota.lt.0.01_rsh .OR. RNquota.gt.0.135_rsh) THEN
                    print*,'Attention RNquota sort des bornes !!!! en ',i,j,', RNquota = ',RNquota
                    stop
                  ENDIF
!               ! ENDIF
!               ! IF(LPquota.lt.p_zost_LPquotamin .OR. LPquota.gt.p_zost_LPquotamax) THEN
!               !   print*,'Attention LPquota sort des bornes !!!! en ',i,j,', LPquota = ',LPquota
!               !   print*,'LPsat, limabsLpo4 =',LPsat,limabsLpo4
                  IF(LPquota.lt.0.0005_rsh .OR. LPquota.gt.0.009_rsh) THEN
                    print*,'Attention LPquota sort des bornes !!!! en ',i,j,', LPquota = ',LPquota
                    stop
                  ENDIF
!               ! ENDIF
!               ! IF(RPquota.lt.p_zost_RPquotamin .OR. RPquota.gt.p_zost_RPquotamax) THEN
!               !   print*,'Attention RPquota sort des bornes !!!! en ',i,j,' RPquota = ',RPquota
!               !   print*,'RP :'
!               !   print*,'Rabs, Recruit_rate =',Rabs_po4*c(iv_zost_RB),-Recruit_rate*p_zost_SB1*c(iv_zost_D)*RPquota
!               !   print*,'transacro,transbasi,RPrecla =',-trans_P_acro,trans_P_basi,RPrecla
!               !   print*,'Root_mort =',-(1-RPrecla)*Root_mort*c(iv_zost_RP)
!               !   print*,'RB :'
!               !   print*,'Root_growth, Recruit_rate, Root_mort =',Root_growth*c(iv_zost_LB),  &
!               !          -Recruit_rate*p_zost_SB1*c(iv_zost_D),-Root_mort*c(iv_zost_RB)
                  IF(RPquota.lt.0.0005_rsh .OR. RPquota.gt.0.009_rsh) THEN
                    print*,'Attention RPquota sort des bornes !!!! en ',i,j,' RPquota = ',RPquota
                    stop
                  ENDIF
!               ! ENDIF

             ELSE
                  RNquota = 0.0_rsh
                  RPquota = 0.0_rsh
                  Gross_Prod_zost = 0.0_rsh         
                  Leaf_Resp_zost = 0.0_rsh 
                  Root_Resp_zost = 0.0_rsh 
                  Root_growth = 0.0_rsh
                  Leaf_growth = 0.0_rsh
                  Labs_nh4 = 0.0_rsh
                  Labs_no3 = 0.0_rsh
                  Labs_po4 = 0.0_rsh
                  Rabs_nh4 = 0.0_rsh
                  Rabs_po4 = 0.0_rsh
                  trans_N_acro = 0.0_rsh
                  trans_N_basi = 0.0_rsh
                  trans_P_acro = 0.0_rsh
                  trans_P_basi = 0.0_rsh
                  LNrecla =  0.0_rsh
                  RNrecla =  0.0_rsh
                  LPrecla =  0.0_rsh
                  RPrecla =  0.0_rsh
                  effet_vent = 0.0_rsh
                  Leaf_mort = 0.0_rsh
                  Root_mort = 0.0_rsh
                  Recruit_rate = 0.0_rsh
                  flowering_rate = 0.0_rsh
          
             END IF     ! fin test presence zosteres
     
             IF(c(iv_zost_benth_seed).gt.0.0_rsh) THEN

               ! Mortalite des graines du sediment (j-1)
               ! ---------------------------------------
               IF((p_zost_Smort*c(iv_zost_benth_seed))*dtbiojour .lt. c(iv_zost_benth_seed)) THEN
                  Seed_mort = p_zost_Smort
               ELSE
                  Seed_mort = 1.0_rsh/dtbiojour
               ENDIF

               ! Recrutement de nouveaux pieds par germination (j-1)
               ! ---------------------------------------------------
               LAI = p_zost_klai*c(iv_zost_LB)
               Qcan = p_parradratio*PAR_top_layer(k,i,j)
               alumfond = Qcan*exp(-p_zost_leafabscoef*LAI)
               IF(c(iv_zost_benth_seed) .gt. p_zost_SB0 .and. alumfond .gt. 150._rsh) THEN
                 IF(c(iv_zost_LB) .le. p_zost_KGERM) THEN
                   limLBgerm = 1-(c(iv_zost_LB)/p_zost_KGERM)**2
                 ELSE
                   limLBgerm = 0._rsh
                 END IF
                 effetchaleurgerm = exp(p_zost_T_prod*temper) 
                 Germin_rate = p_zost_GERmax*limLBgerm*effetchaleurgerm
               ELSE
                 Germin_rate = 0.0_rsh
               ENDIF
             ELSE
                 Germin_rate = 0.0_rsh
             ENDIF      ! fin test sur presence de graines benthiques

             IF(c(iv_zost_seed).gt.0.0_rsh) THEN

               ! Mortalite des graines en suspension (j-1)
               ! -----------------------------------------
               IF((p_zost_Smort*c(iv_zost_seed))*dtbiojour .lt. c(iv_zost_seed)) THEN
                 Seed_mort = p_zost_Smort
               ELSE
                 Seed_mort = 1.0_rsh/dtbiojour
               ENDIF
             ELSE
                 Seed_mort = 0.0_rsh
             ENDIF


             !++++++++++++++++++++++++++++++++++++
             ! Processus lies aux zosteres 
             !++++++++++++++++++++++++++++++++++++
    
             ! Evolution du nitrate dissous (mmolN/m3/j)
             ! -----------------------------------------
              dc(iv_nutr_NO3) = dc(iv_nutr_NO3) - Labs_no3*c(iv_zost_LB)/epn                                    
        
             ! Evolution de l ammonium (mmol/m3/j)
             !!!! ??????????  c est quoi cette division par 10 ???????????????
             ! -----------------------------------
              dc(iv_nutr_NH4) = dc(iv_nutr_NH4) - Labs_nh4*c(iv_zost_LB)/epn + reminazdeteau/10.0_rsh*c(iv_detr_zost_N)      
 
             ! Evolution du phosphore dissous (mmol/m3/j)
             ! ------------------------------------------   
              dc(iv_nutr_PO4) = dc(iv_nutr_PO4) - Labs_po4*c(iv_zost_LB)/epn + reminpdeteau/10.0_rsh*c(iv_detr_zost_P)      

             ! Evolution du phosphore dissous benthique (mmol/m3/j)
             ! ------------------------------------------------------   
              dc(iv_benth_PO4) = - Rabs_po4*c(iv_zost_RB)/epnsed

             ! Evolution de l ammonium benthique (mmol/m3/j)
             ! -----------------------------------------------
              dc(iv_benth_NH4) = - Rabs_nh4*c(iv_zost_RB)/epnsed

             ! Evolution des feuilles de zosteres (mmolC/m2/j)
             ! -----------------------------------------------
              dc(iv_zost_LB) = Leaf_growth*c(iv_zost_LB) + Recruit_rate*p_zost_SB1*c(iv_zost_D) - (Leaf_mort+effet_vent)*c(iv_zost_LB) &
                              + 0.5_rsh*Germin_rate*c(iv_zost_benth_seed)
 
             ! Evolution des rhizomes et racines de zosteres (mmolC/m2/j) 
             ! ----------------------------------------------------------
              dc(iv_zost_RB) = Root_growth*c(iv_zost_LB) - Recruit_rate*p_zost_SB1*c(iv_zost_D) - Root_mort*c(iv_zost_RB)  &
                               + 0.5_rsh*Germin_rate*c(iv_zost_benth_seed)
 
             ! Evolution de la densite de pieds de zosteres (/m2/j)
             ! ----------------------------------------------------
              dc(iv_zost_D) = (Recruit_rate - (Leaf_mort+effet_vent))*c(iv_zost_D)                                                &
                              + Germin_rate*c(iv_zost_benth_seed)*epnsed/p_zost_SB0
 
             ! Evolution du pool de N des feuilles de zosteres (mmolN/m2/j)
             ! ------------------------------------------------------------
              dc(iv_zost_LN) = (Labs_nh4 + Labs_no3)*c(iv_zost_LB)                                                 &
                               + Recruit_rate*p_zost_SB1*c(iv_zost_D)*RNquota                                             &
                               + trans_N_acro - trans_N_basi - ((1-LNrecla)*Leaf_mort+effet_vent)*c(iv_zost_LN)           &
                               + 0.5_rsh*Germin_rate*c(iv_zost_benth_seed)*SNquota 
        
             ! Evolution du pool de P des feuilles de zosteres (mmolP/m2/j)
             ! ------------------------------------------------------------
              dc(iv_zost_LP) = Labs_po4*c(iv_zost_LB) + Recruit_rate*p_zost_SB1*c(iv_zost_D)*RPquota                      &
                               + trans_P_acro - trans_P_basi - ((1-LPrecla)*Leaf_mort+effet_vent)*c(iv_zost_LP)           &
                              + 0.5_rsh*Germin_rate*c(iv_zost_benth_seed)*SPquota 

             ! Evolution du pool de N des rhizomes de zosteres (mmolN/m2/j)
             ! -----------------------------------------------------------
              dc(iv_zost_RN) = Rabs_nh4*c(iv_zost_RB) - Recruit_rate*p_zost_SB1*c(iv_zost_D)*RNquota                      &
                               - trans_N_acro + trans_N_basi - (1-RNrecla)*Root_mort*c(iv_zost_RN)                        &
                               + 0.5_rsh*Germin_rate*c(iv_zost_benth_seed)*SNquota

             ! Evolution du pool de P des rhizomes de zosteres (mmolP/m2/j)
             ! ------------------------------------------------------------
              dc(iv_zost_RP) = Rabs_po4*c(iv_zost_RB) - Recruit_rate*p_zost_SB1*c(iv_zost_D)*RPquota                      &
                               - trans_P_acro + trans_P_basi - (1-RPrecla)*Root_mort*c(iv_zost_RP)                        &
                               + 0.5_rsh*Germin_rate*c(iv_zost_benth_seed)*SPquota

             ! Evolution de l azote detritique peu labile dans le sediment mmolN/m2/j
             ! ----------------------------------------------------------------------
              dc(iv_zost_benth_N) = (1-RNrecla)*Root_mort*c(iv_zost_RN) + Seed_mort*c(iv_zost_benth_seed)*SNquota

             ! Evolution du phosphore detritique peu labile dans le sediment mmolP/m2/j
             ! ------------------------------------------------------------------------
              dc(iv_zost_benth_P) = (1-RPrecla)*Root_mort*c(iv_zost_RP) + Seed_mort*c(iv_zost_benth_seed)*SPquota

             ! Evolution des pousses de zosteres derivantes exprimee en mmolN/m3/j
             ! -------------------------------------------------------------------
              dc(iv_detr_zost_N) = (1-flowering_rate)*((1-LNrecla)*Leaf_mort+effet_vent)*c(iv_zost_LN)/epn  &
                                  + Seed_mort*c(iv_zost_seed)*SNquota &
                                  - reminazdeteau/10.0_rsh*c(iv_detr_zost_N)      

             ! Evolution des pousses de zosteres derivantes exprimee en mmolP/m3/j
             ! -------------------------------------------------------------------
              dc(iv_detr_zost_P) = (1-flowering_rate)*((1-LPrecla)*Leaf_mort+effet_vent)*c(iv_zost_LP)/epn &
                                  + Seed_mort*c(iv_zost_seed)*SPquota &
                                  - reminpdeteau/10.0_rsh*c(iv_detr_zost_P)

             ! Evolution des graines derivantes (mmolC/m3/j)
             ! ---------------------------------------------
              dc(iv_zost_seed) = flowering_rate*(Leaf_mort+effet_vent)*c(iv_zost_LB)/epn - Seed_mort*c(iv_zost_seed)

             ! Evolution des graines du sediment (mmolC/m2/j)
             ! ----------------------------------------------
              dc(iv_zost_benth_seed) = -Seed_mort*c(iv_zost_benth_seed) - Germin_rate*c(iv_zost_benth_seed)

             ! Production cumulee des zosteres (gC/m2/j)
             ! -----------------------------------------
             ! dc(iv_zost_pp) = MAX(Net_Prod_zost,0.0_rsh)*c(iv_zost_LB)*12.0_rsh/1000.0_rsh  
    
             ! Pompage de N cumule des feuilles de zosteres (gN/m2/j)
             ! -------------------------------------------------------
             ! dc(iv_zost_LNuptake) = (Labs_nh4 + Labs_no3)*c(iv_zost_LB)*14.0_rsh/1000.0_rsh 
 
             ! Pompage de P cumule des feuilles de zosteres (gP/m2/j)
             ! ------------------------------------------------------
             ! dc(iv_zost_LPuptake) = Labs_po4*c(iv_zost_LB)*31.0_rsh/1000.0_rsh

             ! Pompage de N cumule des racines de zosteres (gN/m2/j)
             ! -------------------------------------------------------
             ! dc(iv_zost_RNuptake) = Rabs_nh4*c(iv_zost_RB)*14.0_rsh/1000.0_rsh  
 
             ! Pompage de P cumule des racines de zosteres (gP/m2/j)
             ! ------------------------------------------------------
             ! dc(iv_zost_RPuptake) = Rabs_po4*c(iv_zost_RB)*31.0_rsh/1000.0_rsh

          ELSE    ! ibbenth = 0

             IF(c(iv_zost_seed).gt.0.0_rsh) THEN

               ! Mortalite des graines en suspension (j-1)
               ! -----------------------------------------
               IF((p_zost_Smort*c(iv_zost_seed))*dtbiojour .lt. c(iv_zost_seed)) THEN
                   Seed_mort = p_zost_Smort
               ELSE
                   Seed_mort = 1.0_rsh/dtbiojour
               ENDIF
             ELSE
               Seed_mort = 0.0_rsh
             ENDIF

             ! Evolution des pousses de zosteres derivantes exprimee en mmolN/m3/j
             ! -------------------------------------------------------------------
             dc(iv_detr_zost_N) = Seed_mort*c(iv_zost_seed)*SNquota - reminazdeteau/10.0_rsh*c(iv_detr_zost_N)

             ! Evolution des pousses de zosteres derivantes exprimee en mmolP/m3/j
             ! -------------------------------------------------------------------
             dc(iv_detr_zost_P) = Seed_mort*c(iv_zost_seed)*SPquota - reminpdeteau/10.0_rsh*c(iv_detr_zost_P)

             ! Evolution des graines derivantes (mmolC/m3/j)
             ! ---------------------------------------------
             dc(iv_zost_seed) = - Seed_mort*c(iv_zost_seed) 
        
             ! Evolution de l ammonium (mmol/m3/j)
             ! -----------------------------------
             dc(iv_nutr_NH4) = dc(iv_nutr_NH4) + reminazdeteau/10.0_rsh*c(iv_detr_zost_N)      
 
             ! Evolution du phosphore dissous (mmol/m3/j)
             ! ------------------------------------------   
             dc(iv_nutr_PO4) = dc(iv_nutr_PO4) + reminpdeteau/10.0_rsh*c(iv_detr_zost_P)      

          END IF  ! fin test sur la maille de fond

        !   IF (i.eq.85 .and. j.eq.161 .and. ibbenth.eq.1) THEN
        !     print*,'h0+xe =',h0(i,j)+xe
        !     print*,'epn =',epn
        !     print*,'Labs_nh4*c(iv_zost_LB)/epn =',Labs_nh4*c(iv_zost_LB)/epn
        !     print*,'Labs_no3*c(iv_zost_LB)/epn =',Labs_no3*c(iv_zost_LB)/epn
        !     print*,'Eabs_no3*c(iv_zost_EB)/p_phyto_CNratio/epn =',Eabs_no3*c(iv_zost_EB)/p_phyto_CNratio/epn
        !     print*,'Rabs_nh4*c(iv_zost_RB)/epn =',Rabs_nh4*c(iv_zost_RB)/epn
        !     print*,'minsed_Nrefrac*c(iv_zost_detRN) =',minsed_Nrefrac*c(iv_zost_detRN)
        !     print*,'EM*c(iv_zost_EB)/p_phyto_CNratio/epn =',EM*c(iv_zost_EB)/p_phyto_CNratio/epn
        !     print*,'Eabs_nh4*c(iv_zost_EB)/p_phyto_CNratio/epn =',Eabs_nh4*c(iv_zost_EB)/p_phyto_CNratio/epn
        !     print*,'Recruit_rate*p_zost_SB1*c(iv_zost_D)*c(iv_zost_RN)/c(iv_zost_RB) =', &
        !             Recruit_rate*p_zost_SB1*c(iv_zost_D)*c(iv_zost_RN)/c(iv_zost_RB)
        !     print*,'trans_N_acro =',trans_N_acro
        !     print*,'trans_N_basi =',trans_N_basi
        !     print*,'Root_mort*c(iv_zost_RN) =',Root_mort*c(iv_zost_RN)
        !     print*,'Leaf_mort*c(iv_zost_LN)/epn =',Leaf_mort*c(iv_zost_LN)/epn
        !     print*,'min_Nrefrac*c(iv_zost_detLN) =',min_Nrefrac*c(iv_zost_detLN)
        !     print*,'*******************************'
        !   ENDIF


            ! ++++++++++++++++++++++++++++++++++++++++++++++++
            !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
            ! ++++++++++++++++++++++++++++++++++++++++++++++++
          IF (ibbenth.eq.1) THEN
             IF (c(iv_zost_LB).gt.1.e-10 .AND. c(iv_zost_RB).gt.1.e-10 .or. c(iv_zost_benth_seed).gt.p_zost_SB0) THEN
        !       diag_3d_wat(irk_diag(id_zost_LBDW),k,i,j)=c(iv_zost_LB)*12/1000/0.37_rsh
        !       diag_3d_wat(irk_diag(id_zost_RBDW),k,i,j)=c(iv_zost_RB)*12/1000/0.34_rsh
               diag_3d_wat(irk_diag(id_zost_Qcan),k,i,j)=Qcan
               diag_3d_wat(irk_diag(id_zost_LAI),k,i,j)=LAI
               diag_3d_wat(irk_diag(id_zost_limlum),k,i,j)=effetlumzost
               diag_3d_wat(irk_diag(id_zost_limprodLN),k,i,j)=limprodLN
               diag_3d_wat(irk_diag(id_zost_limprodLP),k,i,j)=limprodLP
               diag_3d_wat(irk_diag(id_zost_limRN),k,i,j)=limRN
               diag_3d_wat(irk_diag(id_zost_limRP),k,i,j)=limRP
        !       diag_3d_wat(irk_diag(id_zost_effetchaleurzprod),k,i,j)=effetchaleurzprod
        !       diag_3d_wat(irk_diag(id_zost_effetchaleurzrespf),k,i,j)=effetchaleurzrespf
        !       diag_3d_wat(irk_diag(id_zost_effetchaleurzrespr),k,i,j)=effetchaleurzrespr
        !       diag_3d_wat(irk_diag(id_zost_limabsLnh4),k,i,j)=limabsLnh4
        !       diag_3d_wat(irk_diag(id_zost_limabsLno3),k,i,j)=limabsLno3
        !       diag_3d_wat(irk_diag(id_zost_limabsLpo4),k,i,j)=limabsLpo4
        !       diag_3d_wat(irk_diag(id_zost_limabsRnh4),k,i,j)=limabsRnh4
        !       diag_3d_wat(irk_diag(id_zost_limabsRpo4),k,i,j)=limabsRpo4
               diag_3d_wat(irk_diag(id_zost_limselfshad),k,i,j)=limselfshad
               diag_3d_wat(irk_diag(id_zost_limRB),k,i,j)=limRB
        !       diag_3d_wat(irk_diag(id_zost_effetvent),k,i,j)=effet_vent  
        !       diag_3d_wat(irk_diag(id_zost_Germin_rate),k,i,j)=Germin_rate
        !       diag_3d_wat(irk_diag(id_zost_ERS),k,i,j)=flowering_rate
               diag_2d(irk_diag(id_zost_prod),i,j)=MAX(Net_Prod_zost,0.0_rsh)*c(iv_zost_LB)*12.0_rsh/1000.0_rsh
             END IF
          END IF
   
   
   !!======================================================================

#endif

