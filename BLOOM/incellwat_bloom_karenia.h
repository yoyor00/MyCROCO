#ifdef key_karenia

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_karenia  ***
   !&E
   !&E ** Purpose : Calcul des variations dues a la presence de Karenia
   !&E              a l interieur d une maille d eau
   !&E
   !&E       !  2010-03    (A. Vanhoutte, A. Menesguen ) Original code
   !&E
   !&E      use from general modele : temper,sali, WATER_DENSITY  si ca existe
   !&E      use from BIOLink  : PAR_top_layer,extinction, rappaz, ws3, iheure_BIOLink,iminu_BIOLink,isec_BIOLink
   !&E                         l_waterdensity_known ,thicklayerW_C 
   !&E      use from basic bloom modele : c, dtbiojour, epn,ibbenth, txfiltbenthij,, effetchaleur
   !&E      OUTPUT : effetlumiere_day_karenia, diag_3d_wat( id_totalchl et id_karenia...) + dc
   !&E---------------------------------------------------------------------

          deltadensite=0.0_rsh
#if ! defined key_MARS
          IF(l_waterdensity_known) THEN
             IF (k.lt.NB_LAYER_WAT) &
                deltadensite=abs(WATER_DENSITY_kij-WATER_DENSITY_kp1ij)/((epn+thicklayerW_C(k+1,i,j))/2.0_rsh)
          ELSE
#endif
             ! calcul car non calcule par le modele hydro
             IF (k.lt.NB_LAYER_WAT) THEN
               ! Calcul de la difference de densite avec la couche superieure
               ! density of pure water at temperature temper
               row = 999.842594_rlg + (  6.793952d-02 + ( - 9.095290d-03 +                &
                     (1.001685d-04 - 1.120083d-06*temper + 6.536332d-09*temper*temper)*temper )*temper  )*temper
               ! density of sea water at salinity sali and temperature temper
               rost = row + sali * (   0.824493_rlg + (  -4.0899d-03 + ( 7.6438d-05 +                   &
                     (-8.2467d-07+5.3875d-09*temper)*temper )*temper)*temper) +             &
                     sali**1.5_rlg * ( -5.72466d-03 + (1.0227d-04-1.6546d-06*temper)*temper) +  &
                    sali*sali*4.8314d-04
               deltadensite=rost
               salisup=max(0.0_rsh,SAL_BIOLink(k+1,i,j))
               row = 999.842594_rlg + (  6.793952d-02 + ( - 9.095290d-03 + (1.001685d-04 - 1.120083d-06*TEMP_BIOLink(k+1,i,j)  &
                                       + 6.536332d-09*TEMP_BIOLink(k+1,i,j)*TEMP_BIOLink(k+1,i,j))*TEMP_BIOLink(k+1,i,j) )       &
                                                     *TEMP_BIOLink(k+1,i,j)  )*TEMP_BIOLink(k+1,i,j)
               rost = row + salisup * (   0.824493_rlg + (  -4.0899d-03 + ( 7.6438d-05 +                   &
                                        (-8.2467d-07+5.3875d-09*TEMP_BIOLink(k+1,i,j))*TEMP_BIOLink(k+1,i,j) )          &
                                                 *TEMP_BIOLink(k+1,i,j))*TEMP_BIOLink(k+1,i,j)) + salisup**1.5_rlg *    &
                                         ( -5.72466d-03 + (1.0227d-04-1.6546d-06*TEMP_BIOLink(k+1,i,j))*TEMP_BIOLink(k+1,i,j)) +  &
                                            salisup*salisup*4.8314d-04 
               deltadensite=abs(deltadensite-rost)/((epn+thicklayerW_C(k+1,i,j))/2.0_rsh)
             ENDIF
#if ! defined key_MARS
          ENDIF
#endif
           
          ! effetchaleur : courbe de Blanchard et al. (1996)
          effetchaleurkarenia=0.0_rsh
          if (temper.lt.p_karenia_templethal)  effetchaleurkarenia=              &
            ((p_karenia_templethal-temper)/(p_karenia_templethal-p_karenia_tempopt))**p_karenia_beta   &
            *exp(-p_karenia_beta*((p_karenia_templethal-temper)/(p_karenia_templethal-p_karenia_tempopt)-1))

          ! Effetlumiere : Formulation de Smith integree sur la profondeur de la boite
          ! ............................................................

          fluxrelatifsurf=PAR_top_layer(k,i,j)/p_karenia_iksmith
          fluxrelatiffond=PAR_top_layer(k-1,i,j)/p_karenia_iksmith
          effetlumierekarenia=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*log((fluxrelatifsurf+         &
                       sqrt(1.0_rsh+fluxrelatifsurf*fluxrelatifsurf))/(fluxrelatiffond+           &
                       sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))
          IF ((effetlumierekarenia.lt.0.0_rsh).or.(effetlumierekarenia.gt.1.0_rsh)) &
                  print*,'effetlumierekar=',effetlumierekarenia,'k=',k,'PAR_top_layer=',PAR_top_layer(k,i,j)
   
          ! effet des nutriments
          ! ............................................................

          effetno3karenia=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_karenia_kNO3+(c(iv_nutr_NH4)*p_karenia_kNO3/p_karenia_kNH4))
          !effetnh4karenia=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_karenia_kNH4+(c(iv_nutr_NO3)*p_karenia_kNH4/p_karenia_kNO3))
          effetnh4karenia=max(0.0_rsh,c(iv_nutr_NH4)-0.3_rsh)/   &
               (max(0.0_rsh,c(iv_nutr_NH4)-0.3_rsh)+p_karenia_kNH4+(c(iv_nutr_NO3)*p_karenia_kNH4/p_karenia_kNO3))

          ! modele sans quota azote
          effetazotekarenia=effetno3karenia+effetnh4karenia
          IF (((c(iv_nutr_NO3).gt.1.e-8_rsh).or.(c(iv_nutr_NH4).gt.1.e-8_rsh)).and.(effetazotekarenia.gt.0.0_rsh)) THEN
                  fractionno3karenia=effetno3karenia/effetazotekarenia
                  fractionnh4karenia=effetnh4karenia/effetazotekarenia
          ELSE
                  fractionno3karenia=0.0_rsh
                  fractionnh4karenia=0.0_rsh
          ENDIF
              
          ! modele avec quota azote mais utilise les variables du sans quota 
          uptakeN=0.0_rsh
          effetazotekarenia=0.0_rsh
          maxquotaN=4.0_rsh*p_karenia_qminN
          IF (c(iv_phyto_karenia_C) .gt. 1.e-10) THEN
             quotaN=c(iv_phyto_karenia_N)/c(iv_phyto_karenia_C)
             IF (quotaN.ge.p_karenia_qminN) THEN
                effetazotekarenia=(quotaN-p_karenia_qminN)/(p_karenia_kqN+quotaN-p_karenia_qminN)/  &
                 ((maxquotaN-p_karenia_qminN)/(p_karenia_kqN+maxquotaN-p_karenia_qminN))
             ENDIF
             IF (quotaN.lt.maxquotaN) THEN
                uptakeN=min((p_karenia_vmaxN*(effetno3karenia+effetnh4karenia)) &
                     ,(0.9_rsh*(c(iv_nutr_NH4)*fractionnh4karenia+c(iv_nutr_NO3)*fractionno3karenia)/dtbiojour) &
                     /c(iv_phyto_karenia_C))
                uptakeN=min(uptakeN,((maxquotaN-quotaN)/dtbiojour))
             ENDIF
          ENDIF

          ! modele sans quota phosphore
          !effetphosphorekarenia=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_karenia_kPO4)     

          ! modele avec quota phosphore
          uptakeP=0.0_rsh
          effetphosphorekarenia=0.0_rsh
          maxquotaP=4.0_rsh*p_karenia_qminP
          IF (c(iv_phyto_karenia_C) .gt. 1.e-10) THEN
              quotaP=c(iv_phyto_karenia_P)/c(iv_phyto_karenia_C)
              IF (quotaP.gt.p_karenia_qminP) THEN
                  effetphosphorekarenia=(quotaP-p_karenia_qminP)/(p_karenia_kqP+quotaP-p_karenia_qminP) /  &
                               ((maxquotaP-p_karenia_qminP)/(p_karenia_kqP+maxquotaP-p_karenia_qminP))
              ENDIF
              IF (quotaP.lt.maxquotaP) THEN
                 uptakeP=min(p_karenia_vmaxP*c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_karenia_kPO4),0.9_rsh*c(iv_nutr_PO4)/dtbiojour &
                             /c(iv_phyto_karenia_C))
                 uptakeP=min(uptakeP,((maxquotaP-quotaP)/dtbiojour))
              ENDIF
          ENDIF
     
          effetselnutkarenia=min(effetazotekarenia,effetphosphorekarenia)
   
          pslimkarenia=min(effetlumierekarenia,effetselnutkarenia)
          rationkarenia=p_karenia_mumax*(0.15_rsh+0.85_rsh*effetchaleurkarenia)*pslimkarenia

          ! taux max de croissance diminue dans les zones a forte ECT(croissance nulle pour ect.ge.0.0005 m2.s-2)
          ! -----------------------------------------
          rationkarenia=rationkarenia*max(0.0_rsh,(1.0_rsh-ECT_kij/0.0005_rsh))

          ! excretion and mortality
          ! -----------------------------------------
          excretkarenia=0.0_rsh
          kareniamorteau=p_karenia_mort*effetchaleur
          IF (c(iv_phyto_karenia_N).lt.p_dino_thhold_mort) kareniamorteau=0.0_rsh
   
          ! Evolution de l ammonium  (micromoles/L N)
          ! -----------------------------------------
          dc(iv_nutr_NH4)=dc(iv_nutr_NH4)-fractionnh4karenia*uptakeN*c(iv_phyto_karenia_C)
          !   dc(iv_nutr_NH4)=dc(iv_nutr_NH4)-fractionnh4karenia*rationkarenia*c(iv_phyto_karenia_N)

          ! Evolution du nitrate  (micromoles/L N)
          ! --------------------------------------   
          dc(iv_nutr_NO3)=dc(iv_nutr_NO3)-fractionno3karenia*uptakeN*c(iv_phyto_karenia_C)
          !  dc(iv_nutr_NO3)=dc(iv_nutr_NO3)-fractionno3karenia*rationkarenia*c(iv_phyto_karenia_N)
      
          ! Evolution du phosphate dissous  (micromoles/L P)
          ! ------------------------------------------------   
          dc(iv_nutr_PO4)=dc(iv_nutr_PO4)- uptakeP*c(iv_phyto_karenia_C)   
           !dc(iv_nutr_PO4)=dc(iv_nutr_PO4)-rationkarenia*c(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)
           !dc(iv_nutr_PO4)=dc(iv_nutr_PO4)-rationkarenia*rappaz*c(iv_phyto_karenia_N)

          ! Evolution du carbone des Karenia (micromoles/L C)
          ! --------------------------------
          IF (c(iv_phyto_karenia_N) .gt. 1.e-10) THEN
             dc(iv_phyto_karenia_C)=c(iv_phyto_karenia_C)     &
                                *(rationkarenia-kareniamorteau-excretkarenia-txfiltbenthij*ibbenth/epn)      &
                              -(broumesozookarenia*c(iv_zoo_meso_N)+broumicrozookarenia*c(iv_zoo_micr_N))  &
                                *c(iv_phyto_karenia_C)/c(iv_phyto_karenia_N)
          ENDIF
       
          ! Evolution de l azote des Karenia (micromoles/L N)
          ! -------------------------------------------------
          dc(iv_phyto_karenia_N)=c(iv_phyto_karenia_C)*uptakeN      &
                          -(kareniamorteau+excretkarenia+txfiltbenthij*ibbenth/epn)*c(iv_phyto_karenia_N) &
                          -broumesozookarenia*c(iv_zoo_meso_N)-broumicrozookarenia*c(iv_zoo_micr_N)
 
          ! Evolution du phosphore des Karenia (micromoles/L P)
          ! ---------------------------------------------------
          IF (c(iv_phyto_karenia_N) .gt. 1.e-10) THEN
             dc(iv_phyto_karenia_P)=c(iv_phyto_karenia_C)*uptakeP   &
                            -(kareniamorteau+excretkarenia+txfiltbenthij*ibbenth/epn)*c(iv_phyto_karenia_P) &
                            - (broumesozookarenia*c(iv_zoo_meso_N)+broumicrozookarenia*c(iv_zoo_micr_N))  &
                                *c(iv_phyto_karenia_P)/c(iv_phyto_karenia_N)
          ENDIF 

         ! Production carbonee des Karenia dans la couche k, cumulee depuis le 1er janvier
         ! -------------------------------------------------------------------------------
         dc(iv_phyto_karenia_pp)=12.e-3_rsh*rationkarenia*c(iv_phyto_karenia_C)*epn
           !  dc(iv_phyto_karenia_pp)=12.e-3_rsh*p_phyto_CNratio*rationkarenia*c(iv_phyto_karenia_N)*epn

         ! Evolution de l azote detritique  (micromoles/L N)
         ! -------------------------------------------------   
         dc(iv_detr_N)=dc(iv_detr_N)+ (kareniamorteau+excretkarenia)*c(iv_phyto_karenia_N)
#if ! defined key_benthos
         dc(iv_detr_N)=dc(iv_detr_N)+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_N)
#endif

         ! Evolution du phosphore detritique  (micromoles/L P)    
         ! ---------------------------------------------------
         !dc(iv_detr_P)=dc(iv_detr_P)+c(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)*kareniamorteau  
         !IF (c(iv_phyto_karenia_N).gt.0.0_rsh) dc(iv_detr_P)=dc(iv_detr_P) &
         !   + c(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)*broumesozookarenia*(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)  &
         !        *(1.0_rsh/c(iv_phyto_karenia_N)/(p_phyto_NPratio*p_phyto_CNratio))  &
         !       -rappaz*broumesozookarenia*(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)
         dc(iv_detr_P)=dc(iv_detr_P)+ kareniamorteau*c(iv_phyto_karenia_P)
         IF (c(iv_phyto_karenia_N).gt.0.0_rsh) &
           dc(iv_detr_P)=dc(iv_detr_P) + excretkarenia*c(iv_phyto_karenia_P)      &
                 + broumesozookarenia*c(iv_zoo_meso_N)*(c(iv_phyto_karenia_P)/c(iv_phyto_karenia_N)-rappaz*assimilmesozoo) &
                 + broumicrozookarenia*c(iv_zoo_micr_N)*(c(iv_phyto_karenia_P)/c(iv_phyto_karenia_N)-rappaz*assimilmicrozoo)

#if ! defined key_benthos
        !dc(iv_detr_P)=dc(iv_detr_P)+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)
        ! dc(iv_detr_P)=dc(iv_detr_P)+rappaz*txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_N)
        dc(iv_detr_P)=dc(iv_detr_P)+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_P)
#endif

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! calcul des vitesses de chute de Karenia  
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ws3(k,iv_phyto_karenia_N,i,j)=0.0_rsh

        !   Vitesse sinusoidale   
        !  version ete 2014 : vitesse de chute des Karenia mise a zero
        !IF (kz(k,i,j) >= 0.006_rsh) ws3(k,iv_phyto_karenia_N,i,j)=0.0003*SIN(2.0_rsh*3.14159_rsh/24.0_rsh*(iheure_BIOLink-18.0_rsh))
        !IF (k==1) ws3(k,iv_phyto_karenia_N,i,j)=min(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))
        !IF (k==NB_LAYER_WAT) ws3(k,iv_phyto_karenia_N,i,j)=max(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))
        ! version Alain Octobre 2015
        !ws3(k,iv_phyto_karenia_N,i,j)=0.001*SIN(2.0_rsh*3.14159_rsh/24.0_rsh*(iheure_BIOLink+9.0_rsh))
        !IF (k==1) ws3(k,iv_phyto_karenia_N,i,j)=min(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))
        !IF (k==NB_LAYER_WAT) ws3(k,iv_phyto_karenia_N,i,j)=max(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))
        ! version 2019 : Vitesse sinusoidale (vers le bas durant la nuit)
        IF (((effetselnutkarenia.lt.0.95_rsh).and.((iheure_BIOLink.lt.6).or.(iheure_BIOLink.gt.18))).or.  &
            ((effetlumierekarenia.gt.0.05_rsh).and.((iheure_BIOLink.ge.6).and.(iheure_BIOLink.le.18))))   &
            ws3(k,iv_phyto_karenia_N,i,j)=-p_karenia_wmax*SIN(2.0_rsh*3.14159_rsh/24.0_rsh*(iheure_BIOLink-2.0_rsh))
        IF (k.lt.NB_LAYER_WAT) THEN
              IF (TEMP_BIOLink(k+1,i,j).gt.(p_karenia_tempopt+4.0_rsh)) & 
                ws3(k,iv_phyto_karenia_N,i,j)=max(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))
        ENDIF
        IF (k==1) ws3(k,iv_phyto_karenia_N,i,j)=min(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))
        IF (k==NB_LAYER_WAT) ws3(k,iv_phyto_karenia_N,i,j)=max(0.0_rsh,ws3(k,iv_phyto_karenia_N,i,j))

        ! bornage de la vitesse de chute des variables particulaires en chaque maille
        ! ---------------------------------------------------------------------------
        ws3(k,iv_phyto_karenia_N,i,j)=sign(MIN(0.95_rlg*epn/dtbio,REAL(ABS(ws3(k,iv_phyto_karenia_N,i,j)),rlg)),  &
                                          ws3(k,iv_phyto_karenia_N,i,j))
   
        ws3(k,iv_phyto_karenia_C,i,j)=ws3(k,iv_phyto_karenia_N,i,j)
        ws3(k,iv_phyto_karenia_P,i,j)=ws3(k,iv_phyto_karenia_N,i,j)
#ifdef key_N_tracer
        DO iso=1,nb_source_tracerN
           ws3(k,iv_phyto_karenia_tra_N(iso),i,j)=ws3(k,iv_phyto_karenia_N,i,j)
#ifdef key_age_tracer
           ws3(k,iv_phyto_karenia_age_tra_N(iso),i,j)=ws3(k,iv_phyto_karenia_N,i,j)
#endif
        ENDDO
#endif	
#ifdef key_P_tracer
        DO iso=1,nb_source_tracerP
          ws3(k,iv_phyto_karenia_tra_P(iso),i,j)=ws3(k,iv_phyto_karenia_N,i,j)
#ifdef key_age_tracer
          ws3(k,iv_phyto_karenia_age_tra_P(iso),i,j)=ws3(k,iv_phyto_karenia_N,i,j)
#endif
        ENDDO
#endif	
        ! ++++++++++++++++++++++++++++++++++++++++++++++++
        !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
        ! ++++++++++++++++++++++++++++++++++++++++++++++++
   
        ! effet lumiere: somme a chaque pas de temps = moyenne a minuit pour memoriser la variable diagnostique
        effetlumiere_day_karenia(k,i,j)=effetlumiere_day_karenia(k,i,j)+effetlumierekarenia*dtbio
        IF(iheure_BIOLink==0 .and. iminu_BIOLink ==0 .and. isec_BIOLink <= dtbio) THEN
            diag_3d_wat(irk_diag(id_karenia_limlight),k,i,j)=effetlumiere_day_karenia(k,i,j)/86400.0_rsh
            effetlumiere_day_karenia(k,i,j)=0.0_rsh
        ENDIF

       diag_3d_wat(irk_diag(id_karenia_limN),k,i,j)=effetazotekarenia
       diag_3d_wat(irk_diag(id_karenia_limP),k,i,j)=effetphosphorekarenia

       !   chlorophylle [mgCChl/m3]  
       diag_3d_wat(irk_diag(id_totalchl),k,i,j)=diag_3d_wat(irk_diag(id_totalchl),k,i,j)+c(iv_phyto_karenia_C)*0.5


#endif

 
