#if defined key_P_tracer
 
   !&E---------------------------------------------------------------------
   !&E                 ***   incellwat_bloom_phosphoretracer  ***
   !&E
   !&E ** Purpose : Tracage dans l ecosysteme du phosphore emis par une ou plusieurs source 
   !&E              
   !&E       !  2010-03    (A. Menesguen ) Original code
   !&E
   !&E      use from general modele : name_varbio
   !&E      use from BIOLink  : 
   !&E      use from basic bloom modele : c, epn,ibbenth,rappaz, txfiltbenthij, rNC
   !&E                reminazdeteau,fractionnh4diat,rationdiat,fractionnh4dino,rationdino
   !&E               fractionnh4nano,rationnano ,xnitrifeau, excretionmesozoo,excretionmicrozoo
   !&E               fractionno3diat,rationdiat,fractionno3dino,rationdino,fractionno3nano,rationnano  
   !&E               assimilmicrozoo, broumicrozoonano,broumicrozoodiat,broumicrozoodet,excretionmicrozoo,txmortmicrozoo,broumesozoomicrozoo
   !&E               assimilmesozoo,broumesozoodiat,broumesozoodino,broumesozoomicrozoo,excretionmesozoo,txmortmesozoo
   !&E               diatmorteau,dinomorteau, nanomorteau   
   !&E               desorpeau,adsorpeau,reminpdeteau        
   !&E                                 + pour options   
   !&E      OUTPUT :   dc, diag_3d(irk_diag(id_sign...
   !&E---------------------------------------------------------------------
            

          seuil_P_tracer=seuil_N_tracer/p_phyto_NPratio
          seuil_P_age_tracer=seuil_N_age_tracer/p_phyto_NPratio
          seuil_P_tracer_diag=seuil_N_tracer_diag/p_phyto_NPratio*0.5_rsh
          seuil_P_age_tracer_phytot=seuil_N_age_tracer_phytot


          ! Calcul de la fraction de traceur dans les variables phophorees  (s.u)
          ! ------------------------------------------------------------
          DO iso=1,nb_source_tracerP
             DO ivtra=1,nb_var_tracerP 
               iv_tra=iv_tracer_P(ivtra,iso)
               signature_P(iv_tra)=0.0_rsh
               iv_sign=iv_signed_P(ivtra,iso)
               IF (c(iv_sign).gt.seuil_P_tracer) signature_P(iv_tra)=c(iv_tra)/c(iv_sign)
               !IF (signature_P(iv_tra).gt.1000.0_rsh) print*,'im=',im,'jm=',jm,'k=',k, &
               !             'signature_P',TRIM(ADJUSTL(ADJUSTR(name_varbio(isubs_tracer_P(ivtra,iso))))),signature_P(iv_tra)

#ifdef key_age_tracer
              ! Calcul des ages du traceur dans les variables phosphorees  (jours)
              ! ----------------------------------------------------------------
               ivage=iv_age_P(ivtra,iso)
               age_P(ivage)=0.0_rsh
               IF (signature_P(iv_tra) .GT. seuil_P_age_tracer) age_P(ivage)=c(ivage)/c(iv_tra)
#endif
             ENDDO
   
       
             ! Evolution de variables traceurs marques (micromoles/L P)
             ! -------------------------------------------------
             ! Evolution du phosphate marque dissous
             dc(iv_nutr_PO4_tra_P(iso))=desorpeau*c(iv_nutr_Pads_tra_P(iso))            &
                     -adsorpeau*c(iv_nutr_PO4_tra_P(iso))             &
                     -rappaz*(rationdiat*c(iv_phyto_diat_N)+rationdino*c(iv_phyto_dino_N)+  &
                              rationnano*c(iv_phyto_nano_N))*signature_P(iv_nutr_PO4_tra_P(iso)) &
                     +excretionmesozoo*c(iv_zoo_meso_tra_P(iso))+excretionmicrozoo*c(iv_zoo_micr_tra_P(iso))    &
                     +reminpdeteau*c(iv_detr_tra_P(iso))
 
#ifdef key_karenia
             dc(iv_nutr_PO4_tra_P(iso))=dc(iv_nutr_PO4_tra_P(iso))- &
                     rationkarenia*c(iv_phyto_karenia_C)*signature_P(iv_nutr_PO4_tra_P(iso))/(p_phyto_NPratio*p_phyto_CNratio)
#endif
#ifdef key_phaeocystis
             dc(iv_nutr_PO4_tra_P(iso))=dc(iv_nutr_PO4_tra_P(iso))-(rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)    &
                    +rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))*signature_P(iv_nutr_PO4_tra_P(iso))*rappaz
#endif
#ifdef key_psnz
             dc(iv_nutr_PO4_tra_P(iso))=dc(iv_nutr_PO4_tra_P(iso))-rationpsnz*c(iv_phyto_psnz_N)*signature_P(iv_nutr_PO4_tra_P(iso))*rappaz
#endif
#ifdef key_benthos
             dc(iv_nutr_PO4_tra_P(iso))=dc(iv_nutr_PO4_tra_P(iso))+reminpdeteau*p_reminbenth*c(iv_benth_tra_P(iso))*ibbenth/epn
#endif

#ifdef key_age_tracer
             dc(iv_nutr_PO4_age_tra_P(iso))=c(iv_nutr_PO4_tra_P(iso)) &
                     +desorpeau*c(iv_nutr_Pads_age_tra_P(iso))-adsorpeau*c(iv_nutr_PO4_age_tra_P(iso))      &
                     -rappaz*(rationdiat*c(iv_phyto_diat_N)+rationdino*c(iv_phyto_dino_N)+rationnano*c(iv_phyto_nano_N))  &
                             *signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso)) &
                     +excretionmesozoo*c(iv_zoo_meso_age_tra_P(iso))+excretionmicrozoo*c(iv_zoo_micr_age_tra_P(iso))       &
                     +reminpdeteau*c(iv_detr_age_tra_P(iso))
#ifdef key_karenia
             dc(iv_nutr_PO4_age_tra_P(iso))=dc(iv_nutr_PO4_age_tra_P(iso))-rationkarenia*c(iv_phyto_karenia_C) &
                     *signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso))/(p_phyto_NPratio*p_phyto_CNratio)   
#endif
#ifdef key_phaeocystis
             dc(iv_nutr_PO4_age_tra_P(iso))=dc(iv_nutr_PO4_age_tra_P(iso))-(rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)    &
                       +rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))  &
                       *signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso))*rappaz
#endif
#ifdef key_psnz
             dc(iv_nutr_PO4_age_tra_P(iso))=dc(iv_nutr_PO4_age_tra_P(iso))  &
                       -rationpsnz*c(iv_phyto_psnz_N)*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso))*rappaz
#endif
#ifdef key_benthos
             dc(iv_nutr_PO4_age_tra_P(iso))=dc(iv_nutr_PO4_age_tra_P(iso))+reminpdeteau*p_reminbenth*c(iv_benth_age_tra_P(iso))*ibbenth/epn
#endif
#endif
   
             ! Evolution du phosphore marque adsorbe
             ! --------------------------------------
             dc(iv_nutr_Pads_tra_P(iso))=adsorpeau*c(iv_nutr_PO4_tra_P(iso))-desorpeau*c(iv_nutr_Pads_tra_P(iso))

#ifdef key_age_tracer
             dc(iv_nutr_Pads_age_tra_P(iso))=adsorpeau*c(iv_nutr_PO4_age_tra_P(iso))-desorpeau*c(iv_nutr_Pads_age_tra_P(iso))
#endif
   
             ! Evolution du phosphore marque via l azote marque des diatomees  (micromoles/L P)
             ! -------------------------------------------------------------
             dc(iv_phyto_diat_tra_P(iso))=rationdiat*c(iv_phyto_diat_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))  &
                     -(diatmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_diat_tra_P(iso))    &
                     -(broumesozoodiat*c(iv_zoo_meso_N)+broumicrozoodiat*c(iv_zoo_micr_N)) &
                          *rappaz*signature_P(iv_phyto_diat_tra_P(iso))

#ifdef key_age_tracer
              !     dc(iv_phyto_diat_age_tra_P(iso))=c(iv_phyto_diat_age_tra_P(iso))  &
             dc(iv_phyto_diat_age_tra_P(iso))=c(iv_phyto_diat_tra_P(iso))  &
                      +rationdiat*c(iv_phyto_diat_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso)) &
                     -(diatmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_diat_age_tra_P(iso))    &
                     -(broumesozoodiat*c(iv_zoo_meso_N)+broumicrozoodiat*c(iv_zoo_micr_N))  &
                           *rappaz*signature_P(iv_phyto_diat_tra_P(iso))*age_P(iv_phyto_diat_age_tra_P(iso))
#endif
   
             ! Evolution du phosphore marque via l azote marque des dinoflagelles  (micromoles/L P)
             ! -----------------------------------------------------------------
             dc(iv_phyto_dino_tra_P(iso))=rationdino*c(iv_phyto_dino_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))  &
                     -(dinomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_dino_tra_P(iso))  &
                       -broumesozoodino*c(iv_zoo_meso_N)*signature_P(iv_phyto_dino_tra_P(iso))*rappaz

#ifdef key_age_tracer  
             dc(iv_phyto_dino_age_tra_P(iso))=c(iv_phyto_dino_tra_P(iso))  &
                       +rationdino*c(iv_phyto_dino_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso)) &
                       -(dinomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_dino_age_tra_P(iso))   &
                       -broumesozoodino*c(iv_zoo_meso_N) &
                             *rappaz*signature_P(iv_phyto_dino_tra_P(iso))*age_P(iv_phyto_dino_age_tra_P(iso))
#endif
   
             ! Evolution du phosphore marque via l azote marque des nanoflagelles  (micromoles/L P)
             ! -----------------------------------------------------------------
             dc(iv_phyto_nano_tra_P(iso))=rationnano*c(iv_phyto_nano_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))  & 
                   -(nanomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_nano_tra_P(iso)) &
                    -broumicrozoonano*c(iv_zoo_micr_N)*rappaz*signature_P(iv_phyto_nano_tra_P(iso))

#ifdef key_age_tracer  
             dc(iv_phyto_nano_age_tra_P(iso))=c(iv_phyto_nano_tra_P(iso))  &
                       +rationnano*c(iv_phyto_nano_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso)) &
                        -(nanomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_nano_age_tra_P(iso))   &
                        -broumicrozoonano*c(iv_zoo_micr_N)  &
                             *rappaz*signature_P(iv_phyto_nano_tra_P(iso))*age_P(iv_phyto_nano_age_tra_P(iso))
#endif
      
#ifdef key_karenia

             ! Evolution du phosphore marque via l azote marque de Karenia mikimotoi (micromoles/L P)
             ! ---------------------------------------------------------------   
             dc(iv_phyto_karenia_tra_P(iso))= uptakeP*c(iv_phyto_karenia_C)*signature_P(iv_nutr_PO4_tra_P(iso))  &
                           -(kareniamorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_karenia_tra_P(iso))  &
                            -broumesozookarenia*c(iv_zoo_meso_N)*rappaz*signature_P(iv_phyto_karenia_tra_P(iso))
#ifdef key_age_tracer
             dc(iv_phyto_karenia_age_tra_P(iso))=c(iv_phyto_karenia_tra_P(iso))   &
                      +uptakeP*c(iv_phyto_karenia_C)*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso)) &
                      -(kareniamorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_karenia_age_tra_P(iso))   &
                      -broumesozookarenia*c(iv_zoo_meso_N)  &
                           *rappaz*signature_P(iv_phyto_karenia_tra_P(iso))*age_P(iv_phyto_karenia_age_tra_P(iso))
#endif
#endif

#ifdef key_phaeocystis

             ! Evolution du phosphore marque via l azote marque des colonies de Phaeocystis globosa (micromoles/L P)
             ! ----------------------------------------------------------------------------------
             dc(iv_phyto_phaeocystis_colo_tra_P(iso))=rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N) &
                                                 *rappaz*signature_P(iv_nutr_PO4_tra_P(iso))  & 
                                -(phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_tra_P(iso))   &
                                  -broumesozoophaeocolo*c(iv_zoo_meso_N)  &
                                                  *rappaz*signature_P(iv_phyto_phaeocystis_colo_tra_P(iso)) &
                                +initcolonie*c(iv_phyto_phaeocystis_cell_tra_P(iso))
#ifdef key_age_tracer
             dc(iv_phyto_phaeocystis_colo_age_tra_P(iso))=c(iv_phyto_phaeocystis_colo_tra_P(iso))  &
                      +rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N) &
                                               *rappaz*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso))  & 
                      -(phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_age_tra_P(iso))   &
                      -broumesozoophaeocolo*c(iv_zoo_meso_N)*rappaz*signature_P(iv_phyto_phaeocystis_colo_tra_P(iso))   &
                         *age_P(iv_phyto_phaeocystis_colo_age_tra_P(iso))+initcolonie*c(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif

             ! Evolution du phophore marque via l azote marque des cellules libres de Phaeocystis globosa (micromoles/L P)
             ! -----------------------------------------------------------------------------------------
             dc(iv_phyto_phaeocystis_cell_tra_P(iso))=rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N) &
                                                             *rappaz*signature_P(iv_nutr_PO4_tra_P(iso))    & 
                     -(phaeocystismortcell+txfiltbenthij*ibbenth/epn)*c(iv_phyto_phaeocystis_cell_tra_P(iso))   &
                     -broumicrozoophaeocell*c(iv_zoo_micr_N)*rappaz*signature_P(iv_phyto_phaeocystis_cell_tra_P(iso)) &
                     -initcolonie*c(iv_phyto_phaeocystis_cell_tra_P(iso))
#ifdef key_age_tracer
            dc(iv_phyto_phaeocystis_cell_age_tra_P(iso))=c(iv_phyto_phaeocystis_cell_tra_P(iso))  &
                      +rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N) &
                                           *rappaz*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso))    & 
                      -(phaeocystismortcell+txfiltbenthij*ibbenth/epn)*c(iv_phyto_phaeocystis_cell_age_tra_P(iso))   &
                      -broumicrozoophaeocell*c(iv_zoo_micr_N)*rappaz*signature_P(iv_phyto_phaeocystis_cell_tra_P(iso))    &
                       *age_P(iv_phyto_phaeocystis_cell_age_tra_P(iso))-initcolonie*c(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif
#endif

#ifdef key_psnz

             ! Evolution du phophore marque via l azote marque des Pseudo-nitzschia (micromoles/L P)
             ! -------------------------------------------------------------------
             dc(iv_phyto_psnz_tra_P(iso))=rationpsnz*c(iv_phyto_psnz_N)*rappaz*signature_P(iv_nutr_PO4_tra_P(iso))  &
                     -(psnzmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_psnz_tra_P(iso))   &
                     -(broumicrozoopsnz*c(iv_zoo_micr_N)+broumesozoopsnz*c(iv_zoo_meso_N)) &
                                  *rappaz*signature_P(iv_phyto_psnz_tra_P(iso))
#ifdef key_age_tracer
             dc(iv_phyto_psnz_age_tra_P(iso))=c(iv_phyto_psnz_tra_P(iso))  &
                       +rationpsnz*c(iv_phyto_psnz_N) &
                          *rappaz*signature_P(iv_nutr_PO4_tra_P(iso))*age_P(iv_nutr_PO4_age_tra_P(iso))  &
                      -(psnzmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_psnz_age_tra_P(iso))  &
                      -(broumicrozoopsnz*c(iv_zoo_micr_N)+broumesozoopsnz*c(iv_zoo_meso_N)) &
                          *rappaz*signature_P(iv_phyto_psnz_tra_P(iso))*age_P(iv_phyto_psnz_age_tra_P(iso))
#endif
#endif

             ! Evolution du phophore marque via l azote marque du microzooplancton (micromoles/L P)
             ! ------------------------------------------------------------------  
             dc(iv_zoo_micr_tra_P(iso))=assimilmicrozoo*c(iv_zoo_micr_N) &
                             *rappaz*(broumicrozoonano*signature_P(iv_phyto_nano_tra_P(iso))  &
                                     +broumicrozoodiat*signature_P(iv_phyto_diat_tra_P(iso))  &
                                    + broumicrozoodet*signature_P(iv_detr_tra_P(iso)))         &
                                   -(excretionmicrozoo+txmortmicrozoo)*c(iv_zoo_micr_tra_P(iso))   &
                                  - broumesozoomicrozoo*c(iv_zoo_meso_N)*rappaz*signature_P(iv_zoo_micr_tra_P(iso))

#ifdef key_phaeocystis
             dc(iv_zoo_micr_tra_P(iso))=dc(iv_zoo_micr_tra_P(iso)) + assimilmicrozoo*c(iv_zoo_micr_N) &
                                  *rappaz*broumicrozoophaeocell*signature_P(iv_phyto_phaeocystis_cell_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_micr_tra_P(iso))=dc(iv_zoo_micr_tra_P(iso)) + assimilmicrozoo*c(iv_zoo_micr_N) &
                                  *rappaz*broumicrozoopsnz*signature_P(iv_phyto_psnz_tra_P(iso))
#endif

#ifdef key_age_tracer
             dc(iv_zoo_micr_age_tra_P(iso))=c(iv_zoo_micr_tra_P(iso)) +assimilmicrozoo*c(iv_zoo_micr_N) &
                       *rappaz*(broumicrozoonano*signature_P(iv_phyto_nano_tra_P(iso))*age_P(iv_phyto_nano_age_tra_P(iso))  &
                               +broumicrozoodiat*signature_P(iv_phyto_diat_tra_P(iso))*age_P(iv_phyto_diat_age_tra_P(iso))  &
                               +broumicrozoodet*signature_P(iv_detr_tra_P(iso))*age_P(iv_detr_age_tra_P(iso)))   &
                       -(excretionmicrozoo+txmortmicrozoo)*c(iv_zoo_micr_age_tra_P(iso))  &
                       - broumesozoomicrozoo*c(iv_zoo_meso_N)* &
                              rappaz*signature_P(iv_zoo_micr_tra_P(iso))*age_P(iv_zoo_micr_age_tra_P(iso))
#ifdef key_phaeocystis
              dc(iv_zoo_micr_age_tra_P(iso))=dc(iv_zoo_micr_age_tra_P(iso)) + assimilmicrozoo*c(iv_zoo_micr_N)  &
                   *rappaz*broumicrozoophaeocell*signature_P(iv_phyto_phaeocystis_cell_tra_P(iso)) &
                                   *age_P(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif
#ifdef key_psnz
              dc(iv_zoo_micr_age_tra_P(iso))=dc(iv_zoo_micr_age_tra_P(iso))+ assimilmicrozoo*c(iv_zoo_micr_N)  &
                   *rappaz*broumicrozoopsnz*signature_P(iv_phyto_psnz_tra_P(iso))*age_P(iv_phyto_psnz_age_tra_P(iso))
#endif
#endif

             ! Evolution du phophore marque via l azote marque du mesozooplancton (micromoles/L P)
             ! -----------------------------------------------------------------
             dc(iv_zoo_meso_tra_P(iso))=assimilmesozoo*c(iv_zoo_meso_N) &
                   *rappaz*(broumesozoodiat*signature_P(iv_phyto_diat_tra_P(iso))  &
                           +broumesozoodino*signature_P(iv_phyto_dino_tra_P(iso))  &
                           +broumesozoomicrozoo*signature_P(iv_zoo_micr_tra_P(iso)))   &
                          -(excretionmesozoo+txmortmesozoo)*c(iv_zoo_meso_tra_P(iso))
#ifdef key_phaeocystis
             dc(iv_zoo_meso_tra_P(iso))=dc(iv_zoo_meso_tra_P(iso)) +assimilmesozoo*c(iv_zoo_meso_N) &
                   *rappaz*broumesozoophaeocolo*signature_P(iv_phyto_phaeocystis_colo_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_meso_tra_P(iso))=dc(iv_zoo_meso_tra_P(iso))+assimilmesozoo*c(iv_zoo_meso_N) &
                   *rappaz*broumesozoopsnz*signature_P(iv_phyto_psnz_tra_P(iso))
#endif
#ifdef key_karenia
             dc(iv_zoo_meso_tra_P(iso))=dc(iv_zoo_meso_tra_P(iso))+assimilmesozoo*c(iv_zoo_meso_N) &
                   *rappaz*broumesozookarenia*signature_P(iv_phyto_karenia_tra_P(iso))
#endif

#ifdef key_age_tracer
             dc(iv_zoo_meso_age_tra_P(iso))=c(iv_zoo_meso_tra_P(iso)) +assimilmesozoo*c(iv_zoo_meso_N) &
                   *rappaz*(broumesozoodiat*signature_P(iv_phyto_diat_tra_P(iso))*age_P(iv_phyto_diat_age_tra_P(iso))   &
                           +broumesozoodino*signature_P(iv_phyto_dino_tra_P(iso))*age_P(iv_phyto_dino_age_tra_P(iso))  &
                           +broumesozoomicrozoo*signature_P(iv_zoo_micr_tra_P(iso))*age_P(iv_zoo_micr_age_tra_P(iso)))   &
                           -(excretionmesozoo+txmortmesozoo)*c(iv_zoo_meso_age_tra_P(iso))
#ifdef key_phaeocystis
             dc(iv_zoo_meso_age_tra_P(iso))=dc(iv_zoo_meso_age_tra_P(iso)) +assimilmesozoo*c(iv_zoo_meso_N) &
                   *rappaz*broumesozoophaeocolo*signature_P(iv_phyto_phaeocystis_colo_tra_P(iso)) &
                            *age_P(iv_phyto_phaeocystis_colo_age_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_meso_age_tra_P(iso))=dc(iv_zoo_meso_age_tra_P(iso)) +assimilmesozoo*c(iv_zoo_meso_N)  &
                   *rappaz*broumesozoopsnz*signature_P(iv_phyto_psnz_tra_P(iso))*age_P(iv_phyto_psnz_age_tra_P(iso))
#endif
#ifdef key_karenia
             dc(iv_zoo_meso_age_tra_P(iso))=dc(iv_zoo_meso_age_tra_P(iso)) +assimilmesozoo*c(iv_zoo_meso_N)   &
                   *rappaz*broumesozookarenia*signature_P(iv_phyto_karenia_tra_P(iso))*age_P(iv_phyto_karenia_age_tra_P(iso))
#endif
#endif
   
             ! Evolution du phosphore detritique marque (micromoles/L P)
             ! ---------------------------------------------------------
             dc(iv_detr_tra_P(iso))=diatmorteau*c(iv_phyto_diat_tra_P(iso))                       &
                -broumicrozoodet*c(iv_zoo_micr_N)*signature_P(iv_detr_tra_P(iso))*rNC*rappaz          &
                 -reminpdeteau*c(iv_detr_tra_P(iso))                                                   &
                +(dinomorteau*c(iv_phyto_dino_tra_P(iso))+nanomorteau*c(iv_phyto_nano_tra_P(iso)))      &
                +((txmortmesozoo+(1.0_rsh-assimilmesozoo)*rationmesozoo)*c(iv_zoo_meso_N)*  &
                     (broumesozoodiat*signature_P(iv_phyto_diat_tra_P(iso)) &
                     +broumesozoodino*signature_P(iv_phyto_dino_tra_P(iso))  &
                     +broumesozoomicrozoo*signature_P(iv_zoo_micr_tra_P(iso)))    &
                + (txmortmicrozoo+(1.0_rsh-assimilmicrozoo)*rationmicrozoo)*c(iv_zoo_micr_N)* &
                     (broumicrozoonano*signature_P(iv_phyto_nano_tra_P(iso)) &
                     +broumicrozoodiat*signature_P(iv_phyto_diat_tra_P(iso)) &
                     +broumicrozoodet*signature_P(iv_detr_tra_P(iso))))*rNC*rappaz
                                                     
#ifdef key_phaeocystis
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+ (1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N) &
               *rappaz*broumesozoophaeocolo*signature_P(iv_phyto_phaeocystis_colo_tra_P(iso))  &
                                                  + (1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)  &
                *rappaz*broumicrozoophaeocell*signature_P(iv_phyto_phaeocystis_cell_tra_P(iso))  &
                   + (phaeocystislysecolo+phaeocystismortcolo)*c(iv_phyto_phaeocystis_colo_tra_P(iso))+ &
                      phaeocystismortcell*c(iv_phyto_phaeocystis_cell_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+((1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*rappaz*broumesozoopsnz+       &
                                                   +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*rappaz*broumicrozoopsnz)  &
                              *signature_P(iv_phyto_psnz_tra_P(iso))+ psnzmorteau*c(iv_phyto_psnz_tra_P(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N) &
                              *rappaz*broumesozookarenia*signature_P(iv_phyto_karenia_tra_P(iso))       &
                              + kareniamorteau*c(iv_phyto_karenia_tra_P(iso))
#endif
#ifdef key_benthos
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))-(-detresusppho*signature_P(iv_benth_tra_P(iso))  &
                               +vitessechutedet*c(iv_detr_tra_P(iso)))/epn
#endif
#if ! defined key_benthos
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+  &
                        txfiltbenthij*ibbenth/epn*(c(iv_phyto_diat_tra_P(iso))+c(iv_phyto_dino_tra_P(iso))+c(iv_phyto_nano_tra_P(iso)))
#ifdef key_phaeocystis
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+ &
                        txfiltbenthij*ibbenth/epn*c(iv_phyto_phaeocystis_cell_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_tra_P(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_tra_P(iso))=dc(iv_detr_tra_P(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_tra_P(iso)) 
#endif      
#endif

#ifdef key_age_tracer
              dc(iv_detr_age_tra_P(iso))=c(iv_detr_tra_P(iso))  &
               +diatmorteau*c(iv_phyto_diat_age_tra_P(iso))    &
                -broumicrozoodet*c(iv_zoo_micr_N)*signature_P(iv_detr_tra_P(iso))*age_P(iv_detr_age_tra_P(iso))*rNC*rappaz   &
                -reminpdeteau*c(iv_detr_age_tra_P(iso))                                                   &
                +(dinomorteau*c(iv_phyto_dino_age_tra_P(iso))+nanomorteau*c(iv_phyto_nano_age_tra_P(iso)))      &
                 +((txmortmesozoo+(1.0_rsh-assimilmesozoo)*rationmesozoo)*c(iv_zoo_meso_N)*  &
                      (broumesozoodiat*signature_P(iv_phyto_diat_tra_P(iso))*age_P(iv_phyto_diat_age_tra_P(iso))  &
                      +broumesozoodino*signature_P(iv_phyto_dino_tra_P(iso))*age_P(iv_phyto_dino_age_tra_P(iso))   &
                      +broumesozoomicrozoo*signature_P(iv_zoo_micr_tra_P(iso))*age_P(iv_zoo_micr_age_tra_P(iso)))  &   
                 +(txmortmicrozoo+(1.0_rsh-assimilmicrozoo)*rationmicrozoo)*c(iv_zoo_micr_N)*  &
                      (broumicrozoonano*signature_P(iv_phyto_nano_tra_P(iso))*age_P(iv_phyto_nano_age_tra_P(iso))  &
                      +broumicrozoodiat*signature_P(iv_phyto_diat_tra_P(iso))*age_P(iv_phyto_diat_age_tra_P(iso))  &
                      +broumicrozoodet*signature_P(iv_detr_tra_P(iso))*age_P(iv_detr_age_tra_P(iso)))) &
                       *rNC*rappaz
                                 
#ifdef key_phaeocystis
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))   &
                 +(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*rappaz &
                  *broumesozoophaeocolo*signature_P(iv_phyto_phaeocystis_colo_tra_P(iso)) &
                            *age_P(iv_phyto_phaeocystis_colo_age_tra_P(iso))  &
                 +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*rappaz*broumicrozoophaeocell  &
                  *signature_P(iv_phyto_phaeocystis_cell_tra_P(iso))*age_P(iv_phyto_phaeocystis_cell_age_tra_P(iso))  & 
                 + (phaeocystislysecolo+phaeocystismortcolo)*c(iv_phyto_phaeocystis_colo_age_tra_P(iso))+ &
                  phaeocystismortcell*c(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))  &
                +((1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*rappaz*broumesozoopsnz+       &
                +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*rappaz*broumicrozoopsnz)  &
                   *signature_P(iv_phyto_psnz_tra_P(iso))*age_P(iv_phyto_psnz_age_tra_P(iso))       &
                + psnzmorteau*c(iv_phyto_psnz_age_tra_P(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))  &
                +(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*rappaz*broumesozookarenia &
                    *signature_P(iv_phyto_karenia_tra_P(iso))*age_P(iv_phyto_karenia_age_tra_P(iso))  &
                + kareniamorteau*c(iv_phyto_karenia_age_tra_P(iso))
#endif
#ifdef key_benthos
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))-  &
                 (-detresusppho*signature_P(iv_benth_tra_P(iso))*age_P(iv_benth_age_tra_P(iso))  &
                 +vitessechutedet*c(iv_detr_age_tra_P(iso)))/epn
#endif
#if ! defined key_benthos
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso)) +txfiltbenthij*ibbenth/epn* &
                 (c(iv_phyto_diat_age_tra_P(iso))+c(iv_phyto_dino_age_tra_P(iso))+c(iv_phyto_nano_age_tra_P(iso)))
#ifdef key_phaeocystis
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_age_tra_P(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_age_tra_P(iso))=dc(iv_detr_age_tra_P(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_age_tra_P(iso)) 
#endif      
#endif
#endif     

#ifdef key_benthos
             ! Evolution du phosphore organique benthique marque  (mmol.m-2.jour-1)
             ! -------------------------------------------------------------------
             dc(iv_benth_tra_P(iso))=(-detresusppho*signature_P(iv_benth_tra_P(iso))+vitessechutedet*c(iv_detr_tra_P(iso))) &
                    -(reminpdeteau*p_reminbenth+p_burial)*c(iv_benth_tra_P(iso))  &
                    +txfiltbenthij*(c(iv_phyto_diat_tra_P(iso))+c(iv_phyto_dino_tra_P(iso))+c(iv_phyto_nano_tra_P(iso)))

#if defined key_psnz
             dc(iv_benth_tra_P(iso))=dc(iv_benth_tra_P(iso))+txfiltbenthij*ibbenth*c(iv_phyto_psnz_tra_P(iso))
#endif
#if defined key_karenia
             dc(iv_benth_tra_P(iso))=dc(iv_benth_tra_P(iso))+txfiltbenthij*ibbenth*c(iv_phyto_karenia_tra_P(iso))
#endif
#ifdef key_phaeocystis
             dc(iv_benth_tra_P(iso))=dc(iv_benth_tra_P(iso))+txfiltbenthij*ibbenth*c(iv_phyto_phaeocystis_cell_tra_P(iso))
#endif

#ifdef key_age_tracer
             dc(iv_benth_age_tra_P(iso))=c(iv_benth_tra_P(iso))+ &
                     (-detresusppho*signature_P(iv_benth_tra_P(iso))*age_P(iv_benth_age_tra_P(iso))  &
                       +vitessechutedet*c(iv_detr_age_tra_P(iso))) &
                     -(reminpdeteau*p_reminbenth+p_burial)*c(iv_benth_age_tra_P(iso))  &
                    +txfiltbenthij*(c(iv_phyto_diat_age_tra_P(iso))+c(iv_phyto_dino_age_tra_P(iso))+c(iv_phyto_nano_age_tra_P(iso)))
#if defined key_psnz
             dc(iv_benth_age_tra_P(iso))=dc(iv_benth_age_tra_P(iso))+txfiltbenthij*ibbenth*c(iv_phyto_psnz_age_tra_P(iso))
#endif
#if defined key_karenia
             dc(iv_benth_age_tra_P(iso))=dc(iv_benth_age_tra_P(iso))+txfiltbenthij*ibbenth*c(iv_phyto_karenia_age_tra_P(iso))
#endif
#ifdef key_phaeocystis
             dc(iv_benth_age_tra_P(iso))=dc(iv_benth_age_tra_P(iso))+txfiltbenthij*ibbenth*c(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif
#ifdef key_ulvas
             dc(iv_benth_age_tra_P(iso))=dc(iv_benth_age_tra_P(iso))+ulvebenthmortsed*c(iv_ulv_benth_age_tra_P(iso))*ibbenth/0.014_rsh
#endif
#endif     
#endif
   
             ! ++++++++++++++++++++++++++++++++++++++++++++++++
             !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
             ! ++++++++++++++++++++++++++++++++++++++++++++++++  
             DO ivtra=1,nb_var_tracerP 
                iv_tra=iv_tracer_P(ivtra,iso)
                iv_sign=iv_signed_P(ivtra,iso)
                id_sign=id_tracer_signP(ivtra,iso)
                diag_3d_wat(id_sign,k,i,j)=0.0_rsh
                IF (c(iv_sign).gt.seuil_P_tracer_diag) diag_3d_wat(id_sign,k,i,j)=signature_P(iv_tra)
            ENDDO
#ifdef key_age_tracer  
            id_sign=id_sign+2
#else
            id_sign=id_sign+1
#endif
            diag_3d_wat(id_sign,k,i,j)=c(iv_phyto_diat_tra_P(iso))+c(iv_phyto_dino_tra_P(iso))+c(iv_phyto_nano_tra_P(iso))
            Pphytototal=(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))*rappaz

#ifdef key_karenia
            diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_karenia_tra_P(iso))
            Pphytototal=Pphytototal+c(iv_phyto_karenia_N)*rappaz
#endif
#ifdef key_phaeocystis
            diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j) &
                           +c(iv_phyto_phaeocystis_colo_tra_P(iso))+c(iv_phyto_phaeocystis_cell_tra_P(iso))
            Pphytototal=Pphytototal+(c(iv_phyto_phaeocystis_colo_N)+c(iv_phyto_phaeocystis_cell_N))*rappaz
#endif
#ifdef key_psnz
            diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_psnz_tra_P(iso))
            Pphytototal=Pphytototal+c(iv_phyto_psnz_N)*rappaz
#endif
            phytototal_tra_P=diag_3d_wat(id_sign,k,i,j)
            IF (Pphytototal.gt.seuil_P_tracer_diag) THEN
              diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)/(Pphytototal)
            ELSE
              diag_3d_wat(id_sign,k,i,j)=0.0_rsh
            ENDIF
   
#ifdef key_age_tracer
            age_max=(CURRENT_TIME-REAL(tdeb_tracerP,rlg))/86400.0_rlg
            DO ivtra=1,nb_var_tracerP 
                iv_tra=iv_tracer_P(ivtra,iso)
                iv_sign=iv_signed_P(ivtra,iso)
                id_sign=id_tracer_ageP(ivtra,iso)
                ivage=iv_age_P(ivtra,iso)
                !IF (signature_P(iv_tra) .GT. 0.02_rsh) THEN
                IF (age_P(ivage)==0.0_rsh) THEN
                   diag_3d_wat(id_sign,k,i,j)=-1.0_rsh
                ELSE IF (age_P(ivage) > age_max) THEN
                   diag_3d_wat(id_sign,k,i,j)=-1.0_rsh
                ELSE
                   diag_3d_wat(id_sign,k,i,j)=age_P(ivage)
                ENDIF
            ENDDO  
            id_sign=id_sign+2
            diag_3d_wat(id_sign,k,i,j)=c(iv_phyto_diat_age_tra_P(iso))+c(iv_phyto_dino_age_tra_P(iso))+c(iv_phyto_nano_age_tra_P(iso))
 
#ifdef key_karenia
            diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_karenia_age_tra_P(iso))
#endif
#ifdef key_phaeocystis
            diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j) &
                        +c(iv_phyto_phaeocystis_colo_age_tra_P(iso))+c(iv_phyto_phaeocystis_cell_age_tra_P(iso))
#endif
#ifdef key_psnz
            diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_psnz_age_tra_P(iso))
#endif
            IF (diag_3d_wat(id_sign-1,k,i,j) > seuil_P_age_tracer_phytot) THEN
               diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)/phytototal_tra_P  
            ELSE
               diag_3d_wat(id_sign,k,i,j)= -1.0_rsh  
            ENDIF
#endif

          ENDDO
         
         
!========================================
#endif

