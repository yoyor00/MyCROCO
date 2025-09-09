#ifdef key_psnz

   !&E---------------------------------------------------------------------
   !&E                 ***   incellwat_pseudonitzschia  ***
   !&E
   !&E ** Purpose : Calcul des variations dues a la presence de Pseudo-nitzschia
   !&E              a l interieur d une maille d eau
   !&E
   !&E       !  2010-02    (C. Penard A. Menesguen ) Original code
   !&E
   !&E      use from general modele : temper,sali,
   !&E      use from BIOLink  : extinction, rappaz, ws3, iheure_BIOLink,iminu_BIOLink,isec_BIOLink 
   !&E      use from basic bloom modele : effetturbidite, c, dtbiojour, epn,ibbenth, txfiltbenthij, l_SNeffect_settle
   !&E                                    effetchaleur
   !&E      if key_physadaptation : jjulien_BIOLink,
   !&E      if key_MANGAbio : , lat2d
   !&E      OUTPUT : effetlumiere_day_psnz, diag_3d_wat( id_totalchl et id_psnz...) + dc
   !&E---------------------------------------------------------------------

  
          ! Effet de la temperature : courbe de Blanchard et al. (1996)
          effetchaleurpsnz=0.0_rsh   
          IF (temper.lt.p_psnz_templethal) &
             effetchaleurpsnz=((p_psnz_templethal-temper)/(p_psnz_templethal-p_psnz_tempopt))**p_psnz_beta &
                   *exp(-p_psnz_beta*((p_psnz_templethal-temper)/(p_psnz_templethal-p_psnz_tempopt)-1))           

          ! Effet de la salinite: courbe de Gauss(moyenne: p_psnz_saliopt, ecart-type: p_psnz_salisigma)
           effetsalinitepsnz=exp(-((sali-30.0_rsh)/10.0_rsh)**2/2.0_rsh)
   
          ! Effetlumiere : Formulation de Smith integree sur la profondeur de la boite
          ! ............................................................
#ifdef key_physadaptation
#ifdef key_MANGAbio
          psnz_iksmith=p_psnz_iksmith*(1.0_rsh+0.5_rsh*SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLink-100)))  &
                *max(0.5_rsh,(1.0_rsh-effetturbidite))*                                   &
                (2.0_rsh*(53.0_rsh-lat2d(i,j))+0.5_rsh*(lat2d(i,j)-43.0_rsh))/10.0_rsh
#else
          psnz_iksmith=p_psnz_iksmith*(1.0_rsh+0.5_rsh*SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLink-100)))  &
                *max(0.5_rsh,(1.0_rsh-effetturbidite))              
#endif
#else
          psnz_iksmith=p_psnz_iksmith
#endif
          fluxrelatifsurf=PAR_top_layer(k,i,j)/psnz_iksmith
          fluxrelatiffond=PAR_top_layer(k-1,i,j)/psnz_iksmith
          effetlumierepsnz=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*log((fluxrelatifsurf+         &
             sqrt(1.0_rsh+fluxrelatifsurf*fluxrelatifsurf))/(fluxrelatiffond+          &
             sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))

          ! Effet des nutriments

          effetno3psnz=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_psnz_kNO3+(c(iv_nutr_NH4)*p_psnz_kNO3/p_psnz_kNH4))
          effetnh4psnz=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_psnz_kNH4+(c(iv_nutr_NO3)*p_psnz_kNH4/p_psnz_kNO3))
          effetazotepsnz=effetno3psnz+effetnh4psnz
          IF (((c(iv_nutr_NO3).gt.1.e-8_rsh).or.(c(iv_nutr_NH4).gt.1.e-8_rsh)).and.(effetazotepsnz.gt.0.0_rsh)) THEN
              fractionno3psnz=effetno3psnz/effetazotepsnz
              fractionnh4psnz=effetnh4psnz/effetazotepsnz
          ELSE
              fractionno3psnz=0.0_rsh
              fractionnh4psnz=0.0_rsh
          ENDIF

          uptakeSi=0.0_rsh
          quotamaxsurminSi=2.0_rsh 
          effetsilicepsnz=0.0_rsh
          IF (c(iv_phyto_psnz_N) .gt. 1.e-10) THEN
             quotaSi=c(iv_phyto_psnz_Si)/c(iv_phyto_psnz_N)
             IF (quotaSi.ge.p_psnz_qminSi) &
                   effetsilicepsnz=(quotaSi-p_psnz_qminSi)/(p_psnz_kqSi+quotaSi-p_psnz_qminSi)/  &
                       ((quotamaxsurminSi-1)/(quotamaxsurminSi-1+p_psnz_kqSi/p_psnz_qminSi))
             IF (quotaSi.lt.(quotamaxsurminSi*p_psnz_qminSi))   &
                    uptakeSi=min((p_psnz_vmaxSi*Si/(p_psnz_kSi+Si)),((0.9_rsh*Si/dtbiojour+dc(iv_nutr_SiOH))/c(iv_phyto_psnz_N)))

          ENDIF
          effetphosphorepsnz=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_psnz_kPO4)

          effetselnutpsnz=min(effetazotepsnz,min(effetsilicepsnz,effetphosphorepsnz))
   
          pslimpsnz=min(effetlumierepsnz,effetselnutpsnz)
          rationpsnz=p_psnz_mumax*(0.15_rsh+0.85_rsh*effetchaleurpsnz)*effetsalinitepsnz*pslimpsnz

          excretpsnz=0.0_rsh
          psnzmorteau=p_psnz_mort*effetchaleur
          !psnzmorteau=0.3_rsh*(1-effetselnutpsnz**0.2)
          IF (c(iv_phyto_psnz_N).lt.(p_diat_thhold_mort/2.0_rsh)) psnzmorteau=0.0_rsh


 
          ! Source d acide domoique
          ! --------------------------------
          prodacidedomo=-0.01_rsh*c(iv_phyto_psnz_ad)
          IF (effetsilicepsnz.lt.p_psnz_thhold_Si) &
             !   Formulation en tout ou rien par rapport au seuil     
              prodacidedomo=p_psnz_prod_domoic*max(c(iv_phyto_psnz_N),0.0_rsh)
             !   Formulation proportionnelle a l effet limitant de la silice     
          !    prodacidedomo=p_psnz_prod_domoic*max(c(iv_phyto_psnz_N),0.0_rsh)  &
          !                   *(p_psnz_thhold_Si-effetsilicepsnz)/p_psnz_thhold_Si

    
          ! Evolution de l ammonium  (micromoles/L N)
          ! -----------------------------------------
          dc(iv_nutr_NH4)=dc(iv_nutr_NH4)-fractionnh4psnz*rationpsnz*c(iv_phyto_psnz_N)
      
          ! Evolution du nitrate  (micromoles/L N)
          ! --------------------------------------   
          dc(iv_nutr_NO3)=dc(iv_nutr_NO3)-fractionno3psnz*rationpsnz*c(iv_phyto_psnz_N)
      
          ! Evolution du phosphate dissous  (micromoles/L P)
          ! ------------------------------------------------   
          dc(iv_nutr_PO4)=dc(iv_nutr_PO4)-rationpsnz*c(iv_phyto_psnz_N)*rappaz
 
          ! Evolution du silicate  (micromoles/L Si)
          ! --------------------------------------   
          dc(iv_nutr_SiOH)=dc(iv_nutr_SiOH)-uptakeSi*c(iv_phyto_psnz_N)  
      
          ! Evolution de l azote des Pseudo-nitzschia
          ! -----------------------------------------
   
          dc(iv_phyto_psnz_N)=c(iv_phyto_psnz_N)*(rationpsnz-psnzmorteau-excretpsnz-txfiltbenthij*ibbenth/epn)  &
                          -broumesozoopsnz*c(iv_zoo_meso_N) -broumicrozoopsnz*c(iv_zoo_micr_N)

          ! Evolution du silicium des Pseudo-nitzschia (micromoles/L N)
          ! -----------------------------------------------------------
          IF (c(iv_phyto_psnz_N) .gt. 1.e-10) THEN 
             dc(iv_phyto_psnz_Si)=c(iv_phyto_psnz_N)*uptakeSi-(psnzmorteau+excretpsnz+txfiltbenthij*ibbenth/epn)*c(iv_phyto_psnz_Si) &
                            -(broumesozoopsnz*c(iv_zoo_meso_N)+broumicrozoopsnz*c(iv_zoo_micr_N))*quotaSi
          ENDIF

          ! Evolution de l acide domoique dans l ensemble eau+cellules (microgrammes/L)
          ! ----------------------------------------------------------
          dc(iv_phyto_psnz_ad)=prodacidedomo-p_psnz_decay_domoic*c(iv_phyto_psnz_ad)
  
          ! Production carbonee des Pseudo-nitzschia dans la couche k, cumulee depuis le 1er janvier
          ! ----------------------------------------------------------------------------------------
          dc(iv_phyto_psnz_pp)=12.e-3_rsh*p_phyto_CNratio*rationpsnz*c(iv_phyto_psnz_N)*epn
            
          ! Evolution de l azote detritique  (micromoles/L N)
          ! -------------------------------------------------   
          dc(iv_detr_N)=dc(iv_detr_N)+ (psnzmorteau+excretpsnz)*c(iv_phyto_psnz_N)
#if ! defined key_benthos
          dc(iv_detr_N)=dc(iv_detr_N)+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_N)
#endif
      
          ! Evolution de la silice particulaire  (micromoles/L Si)
          ! ------------------------------------------------------
          dc(iv_detr_Si)=dc(iv_detr_Si)+ (psnzmorteau+excretpsnz)*c(iv_phyto_psnz_Si)
          IF (c(iv_phyto_psnz_N) .gt. 1.e-10) THEN 
             dc(iv_detr_Si)=dc(iv_detr_Si)+ (broumesozoopsnz*c(iv_zoo_meso_N)+broumicrozoopsnz*c(iv_zoo_micr_N))*quotaSi
          ENDIF  
#if ! defined key_benthos
          dc(iv_detr_Si)=dc(iv_detr_Si)+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_Si)
#endif    
          ! Evolution du phosphore detritique  (micromoles/L P)
          ! ---------------------------------------------------
          dc(iv_detr_P)=dc(iv_detr_P)+rappaz*(psnzmorteau+excretpsnz)*c(iv_phyto_psnz_N)
#if ! defined key_benthos
          dc(iv_detr_P)=dc(iv_detr_P)+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_N)*rappaz
#endif

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! calcul des vitesses de chute de Pseudo-nitzschia si choix 
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(l_SNeffect_settle) THEN
             effetchutepsnz=effetselnutpsnz**0.2_rsh
             ws3(k,iv_phyto_psnz_N,i,j)=(ws_free_min(iv_phyto_psnz_N)*effetchutepsnz+      &
                                  ws_free_max(iv_phyto_psnz_N)*(1.-effetchutepsnz))
             ws3(k,iv_phyto_psnz_Si,i,j)=ws3(k,iv_phyto_psnz_N,i,j)
             ws3(k,iv_phyto_psnz_ad,i,j)=ws3(k,iv_phyto_psnz_N,i,j)

#ifdef key_N_tracer
            DO iso=1,nb_source_tracerN
               ws3(k,iv_phyto_psnz_tra_N(iso),i,j)=ws3(k,iv_phyto_psnz_N,i,j)
#ifdef key_age_tracer
               ws3(k,iv_phyto_psnz_age_tra_N(iso),i,j)=ws3(k,iv_phyto_psnz_N,i,j)
#endif
            ENDDO
#endif	
#ifdef key_P_tracer
            DO iso=1,nb_source_tracerP
               ws3(k,iv_phyto_psnz_tra_P(iso),i,j)=ws3(k,iv_phyto_psnz_N,i,j)
#ifdef key_age_tracer
               ws3(k,iv_phyto_psnz_age_tra_P(iso),i,j)=ws3(k,iv_phyto_psnz_N,i,j)
#endif
            ENDDO
#endif	

          ENDIF

          ! ++++++++++++++++++++++++++++++++++++++++++++++++
          !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
          ! ++++++++++++++++++++++++++++++++++++++++++++++++
  
          ! Effet lumiere : somme a chaque pas de temps = moyenne a minuit pour memoriser la variable diagnostique
          effetlumiere_day_psnz(k,i,j)=effetlumiere_day_psnz(k,i,j)+effetlumierepsnz*dtbio
          IF(iheure_BIOLink==0 .and. iminu_BIOLink ==0 .and. isec_BIOLink <= dtbio) THEN
             diag_3d_wat(irk_diag(id_psnz_limlight),k,i,j)=effetlumiere_day_psnz(k,i,j)/86400.0_rsh
             effetlumiere_day_psnz(k,i,j)=0.0_rsh
          ENDIF

          diag_3d_wat(irk_diag(id_psnz_limN),k,i,j)=effetazotepsnz
          diag_3d_wat(irk_diag(id_psnz_limSi),k,i,j)=effetsilicepsnz
          diag_3d_wat(irk_diag(id_psnz_limP),k,i,j)=effetphosphorepsnz

          !   chlorophylle [mgCChl/m3]   
          diag_3d_wat(irk_diag(id_totalchl),k,i,j)=diag_3d_wat(irk_diag(id_totalchl),k,i,j)+ &
                                                 c(iv_phyto_psnz_N)*fact_phyto_ChlNratio
   !!======================================================================
#endif

