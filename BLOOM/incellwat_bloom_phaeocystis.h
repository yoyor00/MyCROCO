#ifdef key_phaeocystis

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_phaeocystis  ***
   !&E
   !&E ** Purpose : Calcul des variations dues a la presence de phaeocystis
   !&E              a l interieur d une maille d eau
   !&E
   !&E       !  2010-03    (A. Vanhoutte, A. Menesguen ) Original code
   !&E       !  2019-07    (A. Menesguen ) update
   !&E
   !&E      use from general modele : temper,sali, 
   !&E      use from BIOLink  : extinction, rappaz, iheure_BIOLink,iminu_BIOLink,isec_BIOLink
   !&E                           
   !&E      use from basic bloom modele : c, dtbiojour, epn,ibbenth, txfiltbenthij, effetturbidite,effetchaleur
   !&E      OUTPUT : effetlumiere_day_phaeocystis, diag_3d_wat( id_totalchl et id_phaeocystis...) + dc
   !&E---------------------------------------------------------------------

      
          !  Effet de la temperature  : courbe de Blanchard et al. (1996)
          effetchaleurphaeocystis=0.0_rsh   
          IF (temper.lt.p_phaeocystis_templethal)  effetchaleurphaeocystis=  &
               ((p_phaeocystis_templethal-temper)/(p_phaeocystis_templethal-p_phaeocystis_tempopt))**p_phaeocystis_beta &
                  *exp(-p_phaeocystis_beta*((p_phaeocystis_templethal-temper)/(p_phaeocystis_templethal-p_phaeocystis_tempopt)-1))

          ! Effet de la salinite (Peperzak, L. (2002). The wax and wane of Phaeocystis globosa blooms)
          effetsalinitephaeo=(sali-19.0_rsh)/(0.012*(sali-19.0_rsh)**2+0.75*(sali-19.0_rsh)+2.36_rsh)
          if (sali.lt.19.0_rsh) effetsalinitephaeo=0.0_rsh

          ! Effetlumiere : Formulation de Smith integree sur la profondeur de la boite
          ! ................................................................................   
#ifdef key_physadaptation
          ! Ik varie selon la turbidite (Amelioration des resultats sur les cotes belges)
          phaeocystis_iksmith=p_phaeocystis_iksmith*(1.0_rsh-effetturbidite)
          !  test par precaution au cas ou 1-effeturbidite=0.
          IF ( phaeocystis_iksmith == 0.0_rsh) THEN
             fluxrelatifsurf=0.0_rsh
             fluxrelatiffond=0.0_rsh
          ELSE
             fluxrelatifsurf=PAR_top_layer(k,i,j)/phaeocystis_iksmith
             fluxrelatiffond=PAR_top_layer(k-1,i,j)/p_phaeocystis_iksmith
          ENDIF
#else
          fluxrelatifsurf=PAR_top_layer(k,i,j)/p_phaeocystis_iksmith
          fluxrelatiffond=PAR_top_layer(k-1,i,j)/p_phaeocystis_iksmith
#endif
          effetlumierephaeocystis=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*log((fluxrelatifsurf+         &
              sqrt(1.0_rsh+fluxrelatifsurf*fluxrelatifsurf))/(fluxrelatiffond+           &
              sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))

          ! Effet des nutriments sur les colonies
          ! ................................................................................   
          effetno3phaeocystiscolo=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_phaeocystis_kNO3colo+ &
                           (c(iv_nutr_NH4)*p_phaeocystis_kNO3colo/p_phaeocystis_kNH4))
          effetnh4phaeocystiscolo=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+  &
                           p_phaeocystis_kNH4+(c(iv_nutr_NO3)*p_phaeocystis_kNH4/p_phaeocystis_kNO3colo))
          effetazotephaeocystiscolo=effetno3phaeocystiscolo+effetnh4phaeocystiscolo
          IF (((c(iv_nutr_NO3).gt.1.e-8_rsh).or.(c(iv_nutr_NH4).gt.1.e-8_rsh)).and. &
                                           (effetazotephaeocystiscolo.gt.0.0_rsh)) THEN
            fractionno3phaeocystiscolo=effetno3phaeocystiscolo/effetazotephaeocystiscolo
            fractionnh4phaeocystiscolo=effetnh4phaeocystiscolo/effetazotephaeocystiscolo
          ELSE
            fractionno3phaeocystiscolo=0.0_rsh
            fractionnh4phaeocystiscolo=0.0_rsh
          ENDIF

          effetphosphorephaeocystiscolo=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_phaeocystis_kPO4colo)
          effetselnutphaeocystiscolo=min(effetazotephaeocystiscolo,effetphosphorephaeocystiscolo)

          ! Effet des nutriments sur les cellules libres

          effetno3phaeocystiscell=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_phaeocystis_kNO3cell+  &
                           (c(iv_nutr_NH4)*p_phaeocystis_kNO3cell/p_phaeocystis_kNH4))
          effetnh4phaeocystiscell=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_phaeocystis_kNH4+  &
                           (c(iv_nutr_NO3)*p_phaeocystis_kNH4/p_phaeocystis_kNO3cell))
          effetazotephaeocystiscell=effetno3phaeocystiscell+effetnh4phaeocystiscell
          IF (((c(iv_nutr_NO3).gt.1.e-8_rsh).or.(c(iv_nutr_NH4).gt.1.e-8_rsh)).and.(effetazotephaeocystiscell.gt.0.0_rsh)) THEN
              fractionno3phaeocystiscell=effetno3phaeocystiscell/effetazotephaeocystiscell
              fractionnh4phaeocystiscell=effetnh4phaeocystiscell/effetazotephaeocystiscell
          ELSE
             fractionno3phaeocystiscell=0.0_rsh
             fractionnh4phaeocystiscell=0.0_rsh
          ENDIF

          effetphosphorephaeocystiscell=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_phaeocystis_kPO4cell)

          effetselnutphaeocystiscell=min(effetazotephaeocystiscell,effetphosphorephaeocystiscell)  

          ! taux de croissance global   
          ! ............................ 
          pslimphaeocystiscolo=min(effetlumierephaeocystis,effetselnutphaeocystiscolo)
          rationphaeocystiscolo=p_phaeocystis_mumaxcolo*(0.15_rsh+0.85_rsh*effetchaleurphaeocystis)  &
                              *pslimphaeocystiscolo*effetsalinitephaeo

          pslimphaeocystiscell=min(effetlumierephaeocystis,effetselnutphaeocystiscell)
          rationphaeocystiscell=p_phaeocystis_mumaxcell*(0.15_rsh+0.85_rsh*effetchaleurphaeocystis)  &
                              *pslimphaeocystiscell*effetsalinitephaeo
  
          !formation de colonies
          ! ............................ 
          initcolonie=0.0_rsh
          IF(c(iv_phyto_phaeocystis_cell_N).gt.p_phaeocystis_coloseuil) initcolonie=p_phaeocystis_coloinit

          ! mortalite du 1er ordre
          phaeocystismortcolo=p_phaeocystis_mort*effetchaleur
          ! ou ????
          !phaeocystismortcolo=p_phaeocystis_mort*effetchaleurphaeocystis
          IF (c(iv_phyto_phaeocystis_colo_N).lt.0.01_rsh) phaeocystismortcolo=0.0_rsh
          phaeocystismortcell=p_phaeocystis_mort*effetchaleur
          ! ou ??   
          !phaeocystismortcell=p_phaeocystis_mort*effetchaleurphaeocystis
          IF (c(iv_phyto_phaeocystis_cell_N).lt.0.01_rsh) phaeocystismortcell=0.0_rsh
   
          !excretion
          ! ............................ 
          !excretphaeocystis=p_phyto_resp*effetchaleur*(1.0_rsh-effetlumierephaeocystis)
          excretphaeocystis=0.0_rsh 
   
          ! rupture des colonies trop mucilagineuses (contenu en biomasse vivante faible)
          ! ................................................................................   
          Ctomucusratio=c(iv_phyto_phaeocystis_colo_N)*p_phyto_CNratio*12.0_rsh/1000.0_rsh  &
                         /(c(iv_phyto_phaeocystis_mucus)+epsilon_BIOLink)
          !if ((i.eq.55).and.(j.eq.59).and.(k.eq.30)) print*,'Ctomucusratio= ',Ctomucusratio
          IF (Ctomucusratio.le.0.1_rsh)  phaeocystismortcolo=phaeocystismortcolo*3.0_rsh
          
          ! lyse des colonies senescentes
          ! ............................ 
          phaeocystislysecolo=p_phaeocystis_lyse*effetchaleur*(1.0_rsh-effetselnutphaeocystiscolo)

          ! production de mucus en fonction du rapport C/mucus   
          ! ................................................................................   
          mucusprod=0.0_rsh
          IF(Ctomucusratio.gt.0.1_rsh) THEN
                 IF(rationphaeocystiscolo.gt.1.15_rsh) mucusprod=rationphaeocystiscolo*1.2_rsh
          ENDIF
          IF(effetselnutphaeocystiscolo.lt.0.4_rsh) mucusprod=mucusprod*2.0_rsh

          ! Evolution de l ammonium  (micromoles/L N)
          ! -----------------------------------------
          dc(iv_nutr_NH4)=dc(iv_nutr_NH4)-fractionnh4phaeocystiscolo*rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)  &
                                     -fractionnh4phaeocystiscell*rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N)
          !  IF (dc(iv_nutr_NH4).lt.(-c(iv_nutr_NH4)/dtbiojour)) print*,'NH4phaeo  i=',i,'j=',j,'k=',k
      
          ! Evolution du nitrate  (micromoles/L N)
          ! --------------------------------------   
          dc(iv_nutr_NO3)=dc(iv_nutr_NO3)-fractionno3phaeocystiscolo*rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)  &
                                     -fractionno3phaeocystiscell*rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N)
          IF (dc(iv_nutr_NO3).lt.(-c(iv_nutr_NO3)/dtbiojour)) print*,'NO3 initial ',c(iv_nutr_NO3),' epuise en:  i=',i,'j=',j,'k=',k
      
          ! Evolution du phosphate dissous  (micromoles/L P)
          ! ------------------------------------------------   
          dc(iv_nutr_PO4)=dc(iv_nutr_PO4)-(rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)+  &
                                       rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))*rappaz
   
          ! Evolution de l azote des phaeocystis coloniaux
          ! ----------------------------------------------
          dc(iv_phyto_phaeocystis_colo_N)=c(iv_phyto_phaeocystis_colo_N)*  &
                                  (rationphaeocystiscolo-phaeocystismortcolo-excretphaeocystis-phaeocystislysecolo)  &
                                 -broumesozoophaeocolo*c(iv_zoo_meso_N)+c(iv_phyto_phaeocystis_cell_N)*initcolonie

          ! Evolution de l azote des phaeocystis unicellulaires
          ! ---------------------------------------------------
          dc(iv_phyto_phaeocystis_cell_N)=c(iv_phyto_phaeocystis_cell_N)*(rationphaeocystiscell &
                                      -phaeocystismortcell-excretphaeocystis-txfiltbenthij*ibbenth/epn)  &
                             -broumicrozoophaeocell*c(iv_zoo_micr_N) - c(iv_phyto_phaeocystis_cell_N)*initcolonie
        
          ! Production carbonee des phaeocystis dans la couche k, cumulee depuis le 1er janvier
          ! -----------------------------------------------------------------------------------    
          dc(iv_phyto_phaeocystis_pp)=12.e-3_rsh*p_phyto_CNratio*epn*(rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)  &
                                  +rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))

          ! Production de mucus
          ! -------------------     
          dc(iv_phyto_phaeocystis_mucus)=mucusprod*c(iv_phyto_phaeocystis_colo_N) -  &
                                     p_phaeocystis_mucus_decay*c(iv_phyto_phaeocystis_mucus)

      
          ! Evolution de l azote detritique  (micromoles/L N)
          ! -------------------------------------------------   
          dc(iv_detr_N)=dc(iv_detr_N)+ (phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_N)+ &
                                    phaeocystismortcell*c(iv_phyto_phaeocystis_cell_N) + &
                                    excretphaeocystis*(c(iv_phyto_phaeocystis_colo_N)+c(iv_phyto_phaeocystis_cell_N))
#if ! defined key_benthos
          dc(iv_detr_N)=dc(iv_detr_N)+txfiltbenthij*ibbenth/epn*c(iv_phyto_phaeocystis_cell_N)
#endif
          ! Evolution du phosphore detritique  (micromoles/L P)
          ! ---------------------------------------------------
          dc(iv_detr_P)=dc(iv_detr_P)+rappaz*((phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_N)+ &
                                   phaeocystismortcell*c(iv_phyto_phaeocystis_cell_N)+ &
                                   excretphaeocystis*(c(iv_phyto_phaeocystis_colo_N)+c(iv_phyto_phaeocystis_cell_N)))
#if ! defined key_benthos
          dc(iv_detr_P)=dc(iv_detr_P)+txfiltbenthij*ibbenth/epn*c(iv_phyto_phaeocystis_cell_N)*rappaz
#endif

          ! ++++++++++++++++++++++++++++++++++++++++++++++++
          !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
          ! ++++++++++++++++++++++++++++++++++++++++++++++++
   
          ! Effet lumiere : somme a chaque pas de temps = moyenne a minuit pour memoriser la variable diagnostique
          effetlumiere_day_phaeocystis(k,i,j)=effetlumiere_day_phaeocystis(k,i,j)+effetlumierephaeocystis*dtbio
          IF(iheure_BIOLink==0 .and. iminu_BIOLink ==0 .and. isec_BIOLink <= dtbio) THEN
             diag_3d_wat(irk_diag(id_phaeocystis_limlight),k,i,j)=effetlumiere_day_phaeocystis(k,i,j)/86400.0_rsh
             effetlumiere_day_phaeocystis(k,i,j)=0.0_rsh
          ENDIF

         diag_3d_wat(irk_diag(id_phaeocystis_limN),k,i,j)=effetazotephaeocystiscolo
         diag_3d_wat(irk_diag(id_phaeocystis_limP),k,i,j)=effetphosphorephaeocystiscolo
         !   chlorophylle [mgCChl/m3]
         diag_3d_wat(irk_diag(id_totalchl),k,i,j)=diag_3d_wat(irk_diag(id_totalchl),k,i,j)+ &
                                          (c(iv_phyto_phaeocystis_colo_N)+c(iv_phyto_phaeocystis_cell_N))*fact_phyto_ChlNratio  

   !!======================================================================     
#endif


 
