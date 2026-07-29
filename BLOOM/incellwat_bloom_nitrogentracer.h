#if defined key_N_tracer
 
   !&E---------------------------------------------------------------------
   !&E                 ***   incellwat_bloom_nitrogentracer  ***
   !&E
   !&E ** Purpose : Tracage dans l ecosysteme de l azote emis par une ou plusieurs source 
   !&E              
   !&E       !  2010-03    (A. Menesguen ) Original code
   !&E
   !&E      use from general modele : name_var
   !&E      use from BIOLink  : 
   !&E      use from basic bloom modele : c, epn,ibbenth, txfiltbenthij
   !&E                reminazdeteau,fractionnh4diat,rationdiat,fractionnh4dino,rationdino
   !&E               fractionnh4nano,rationnano ,xnitrifeau, excretionmesozoo,excretionmicrozoo
   !&E               fractionno3diat,rationdiat,fractionno3dino,rationdino,fractionno3nano,rationnano  
   !&E               assimilmicrozoo, broumicrozoonano,broumicrozoodiat,broumicrozoodet,excretionmicrozoo,txmortmicrozoo,broumesozoomicrozoo
   !&E               assimilmesozoo,broumesozoodiat,broumesozoodino,broumesozoomicrozoo,excretionmesozoo,txmortmesozoo
   !&E               diatmorteau,dinomorteau, nanomorteau           
   !&E                                 + pour options   
   !&E      OUTPUT :   dc, diag_3d(irk_diag(id_sign...
   !&E---------------------------------------------------------------------
            
 
          ! Calcul de la fraction de traceur dans les variables azoteees  (s.u)
          ! ------------------------------------------------------------
          DO iso=1,nb_source_tracerN
             DO ivtra=1,nb_var_tracerN 
               iv_tra=iv_tracer_N(ivtra,iso)
               signature_N(iv_tra)=0.0_rsh
               iv_sign=iv_signed_N(ivtra,iso)
               IF (c(iv_sign).gt.seuil_N_tracer) signature_N(iv_tra)=c(iv_tra)/c(iv_sign)
               IF (signature_N(iv_tra).gt.1000.0_rsh)  print*,'i=',i,'j=',j,'k=',k,'signature_N', &
                             TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_N(ivtra,iso))))),signature_N(iv_tra)

#ifdef key_age_tracer
          ! Calcul des ages du traceur dans les variables azoteees  (jours)
          ! ----------------------------------------------------------------
               ivage=iv_age_N(ivtra,iso)
               age_N(ivage)=0.0_rsh
               IF (signature_N(iv_tra) .GT. seuil_N_age_tracer) age_N(ivage)=c(ivage)/c(iv_tra)
#endif
             ENDDO
   
       
             ! Evolution de variables traceurs marques (micromoles/L N)
             ! -------------------------------------------------
             ! ammonium
             dc(iv_nutr_NH4_tra_N(iso))=reminazdeteau*c(iv_detr_tra_N(iso))    &
                      -(fractionnh4diat*rationdiat*c(iv_phyto_diat_N)+fractionnh4dino*rationdino*c(iv_phyto_dino_N)+ &
                        fractionnh4nano*rationnano*c(iv_phyto_nano_N))   &
                      *signature_N(iv_nutr_NH4_tra_N(iso))-xnitrifeau*c(iv_nutr_NH4_tra_N(iso))+ &
                        excretionmesozoo*c(iv_zoo_meso_tra_N(iso)) +excretionmicrozoo*c(iv_zoo_micr_tra_N(iso))
 
#ifdef key_karenia
             dc(iv_nutr_NH4_tra_N(iso))=dc(iv_nutr_NH4_tra_N(iso))- &
                                fractionnh4karenia*uptakeN*c(iv_phyto_karenia_C)*signature_N(iv_nutr_NH4_tra_N(iso))   
#endif
#ifdef key_phaeocystis
             dc(iv_nutr_NH4_tra_N(iso))=dc(iv_nutr_NH4_tra_N(iso))- &
                                (fractionnh4phaeocystiscolo*rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)  &
                                 +fractionnh4phaeocystiscell*rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N)) &
                                    *signature_N(iv_nutr_NH4_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_nutr_NH4_tra_N(iso))=dc(iv_nutr_NH4_tra_N(iso))- &
                                fractionnh4psnz*rationpsnz*c(iv_phyto_psnz_N)*signature_N(iv_nutr_NH4_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_nutr_NH4_tra_N(iso))=dc(iv_nutr_NH4_tra_N(iso))- &
                                fractionnh4*signature_N(iv_nutr_NH4_tra_N(iso))*(pompageazoteulve*c(iv_ulv_susdrywght)   &
                                +pompageazoteulvebenth*c(iv_ulv_benthdrywght)*partulvessurface*ibbenth/epn)
#endif
#ifdef key_benthos
             dc(iv_nutr_NH4_tra_N(iso))=dc(iv_nutr_NH4_tra_N(iso))+reminazdeteau*p_reminbenth*c(iv_benth_tra_N(iso))*ibbenth/epn
#endif

#ifdef key_age_tracer
             dc(iv_nutr_NH4_age_tra_N(iso))=c(iv_nutr_NH4_tra_N(iso))+reminazdeteau*c(iv_detr_age_tra_N(iso))      &
                 -(fractionnh4diat*rationdiat*c(iv_phyto_diat_N)+  &
                   fractionnh4dino*rationdino*c(iv_phyto_dino_N)+fractionnh4nano*rationnano*c(iv_phyto_nano_N))      &
                  *signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))-xnitrifeau*c(iv_nutr_NH4_age_tra_N(iso)) &
                  +excretionmesozoo*c(iv_zoo_meso_age_tra_N(iso))+excretionmicrozoo*c(iv_zoo_micr_age_tra_N(iso))	  
#ifdef key_karenia
             dc(iv_nutr_NH4_age_tra_N(iso))=dc(iv_nutr_NH4_age_tra_N(iso))- fractionnh4karenia*uptakeN*   &
                             c(iv_phyto_karenia_C)*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))
#endif
#ifdef key_phaeocystis
             dc(iv_nutr_NH4_age_tra_N(iso))=dc(iv_nutr_NH4_age_tra_N(iso))- &
                         (fractionnh4phaeocystiscolo*rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)    &
                         +fractionnh4phaeocystiscell*rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))*  &
                                     signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_nutr_NH4_age_tra_N(iso))=dc(iv_nutr_NH4_age_tra_N(iso))    &
                         -fractionnh4psnz*rationpsnz*c(iv_phyto_psnz_N)*signature_N(iv_nutr_NH4_tra_N(iso))* &
                                     age_N(iv_nutr_NH4_age_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_nutr_NH4_age_tra_N(iso))=dc(iv_nutr_NH4_age_tra_N(iso))  &
                        -fractionnh4*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))* &
                           (pompageazoteulve*c(iv_ulv_susdrywght)           &
                        +pompageazoteulvebenth*c(iv_ulv_benthdrywght)*partulvessurface*ibbenth/epn)
#endif
#ifdef key_benthos
             dc(iv_nutr_NH4_age_tra_N(iso))=dc(iv_nutr_NH4_age_tra_N(iso))+  &
                           reminazdeteau*p_reminbenth*c(iv_benth_age_tra_N(iso))*ibbenth/epn
#endif
#endif
   
          ! Evolution du nitrate marque (micromoles/L N)
          ! ----------------------------------------------
             dc(iv_nutr_NO3_tra_N(iso))=xnitrifeau*c(iv_nutr_NH4_tra_N(iso))- &
                                        (fractionno3diat*rationdiat*c(iv_phyto_diat_N)+  &
                                         fractionno3dino*rationdino*c(iv_phyto_dino_N)+  &
                                         fractionno3nano*rationnano*c(iv_phyto_nano_N))  &
                                               *signature_N(iv_nutr_NO3_tra_N(iso))

#ifdef key_karenia
             dc(iv_nutr_NO3_tra_N(iso))=dc(iv_nutr_NO3_tra_N(iso))-  &
                                          fractionno3karenia*uptakeN*c(iv_phyto_karenia_C)  &
                                               *signature_N(iv_nutr_NO3_tra_N(iso)) 
#endif
#ifdef key_phaeocystis
             dc(iv_nutr_NO3_tra_N(iso))=dc(iv_nutr_NO3_tra_N(iso))-  &
                                 (fractionno3phaeocystiscolo*rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)  &
                                 +fractionno3phaeocystiscell*rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))  &
                                               *signature_N(iv_nutr_NO3_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_nutr_NO3_tra_N(iso))=dc(iv_nutr_NO3_tra_N(iso))-  &
                                 fractionno3psnz*rationpsnz*c(iv_phyto_psnz_N)*signature_N(iv_nutr_NO3_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_nutr_NO3_tra_N(iso))=dc(iv_nutr_NO3_tra_N(iso))-  &
                                  fractionno3*signature_N(iv_nutr_NO3_tra_N(iso))*(pompageazoteulve*c(iv_ulv_susdrywght)   &
                                   +pompageazoteulvebenth*c(iv_ulv_benthdrywght)*partulvessurface*ibbenth/epn)
#endif

#ifdef key_age_tracer
             dc(iv_nutr_NO3_age_tra_N(iso))=c(iv_nutr_NO3_tra_N(iso))+xnitrifeau*c(iv_nutr_NH4_age_tra_N(iso))  &
                                 -(fractionno3diat*rationdiat*c(iv_phyto_diat_N)+  &
                                   fractionno3dino*rationdino*c(iv_phyto_dino_N)+  &
                                   fractionno3nano*rationnano*c(iv_phyto_nano_N))   &
                                   *signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso))
#ifdef key_karenia
             dc(iv_nutr_NO3_age_tra_N(iso))=dc(iv_nutr_NO3_age_tra_N(iso))     &
                                    - fractionno3karenia*uptakeN*c(iv_phyto_karenia_C)  &
                                   *signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)) 
#endif
#ifdef key_phaeocystis
             dc(iv_nutr_NO3_age_tra_N(iso))=dc(iv_nutr_NO3_age_tra_N(iso))    &
                - (fractionno3phaeocystiscolo*rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)      &
                   +fractionno3phaeocystiscell*rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))    &
                                    *signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_nutr_NO3_age_tra_N(iso))=dc(iv_nutr_NO3_age_tra_N(iso))   &
                                      -fractionno3psnz*rationpsnz*c(iv_phyto_psnz_N)  &
                                    *signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_nutr_NO3_age_tra_N(iso))=dc(iv_nutr_NO3_age_tra_N(iso))   &
                                      -fractionno3*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso))  &
                                        *(pompageazoteulve*c(iv_ulv_susdrywght)                   &
                                         +pompageazoteulvebenth*c(iv_ulv_benthdrywght)*partulvessurface*ibbenth/epn)
#endif
#endif
   
             ! Evolution de l azote marque des diatomees  (micromoles/L N)
             ! -------------------------------------------------------------
             dc(iv_phyto_diat_tra_N(iso))=rationdiat*c(iv_phyto_diat_N)*  &
                     (fractionnh4diat*signature_N(iv_nutr_NH4_tra_N(iso))+  &
                      fractionno3diat*signature_N(iv_nutr_NO3_tra_N(iso)))  &
                     -(diatmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_diat_tra_N(iso))    &
                     -(broumesozoodiat*c(iv_zoo_meso_N)+broumicrozoodiat*c(iv_zoo_micr_N))  &
                                       *signature_N(iv_phyto_diat_tra_N(iso))

#ifdef key_age_tracer
             dc(iv_phyto_diat_age_tra_N(iso))=c(iv_phyto_diat_tra_N(iso)) +rationdiat*c(iv_phyto_diat_N)* &
                      (fractionnh4diat*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))   &
                       +fractionno3diat*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))  &
                     -(diatmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_diat_age_tra_N(iso))    &
                     -(broumesozoodiat*c(iv_zoo_meso_N)+broumicrozoodiat*c(iv_zoo_micr_N))   &
                        *signature_N(iv_phyto_diat_tra_N(iso))*age_N(iv_phyto_diat_age_tra_N(iso))
#endif
   
             ! Evolution de l azote marque des dinoflagelles  (micromoles/L N)
             ! -----------------------------------------------------------------
             dc(iv_phyto_dino_tra_N(iso))=rationdino*c(iv_phyto_dino_N)*  &
                      (fractionnh4dino*signature_N(iv_nutr_NH4_tra_N(iso))+fractionno3dino*signature_N(iv_nutr_NO3_tra_N(iso)))  &
                      -(dinomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_dino_tra_N(iso))-broumesozoodino*c(iv_zoo_meso_N)  &
                             *signature_N(iv_phyto_dino_tra_N(iso))

#ifdef key_age_tracer  
             dc(iv_phyto_dino_age_tra_N(iso))=c(iv_phyto_dino_tra_N(iso)) +rationdino*c(iv_phyto_dino_N)* &
                       (fractionnh4dino*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))   &
                       +fractionno3dino*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))  &
                       -(dinomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_dino_age_tra_N(iso))   &
                        -broumesozoodino*c(iv_zoo_meso_N)*signature_N(iv_phyto_dino_tra_N(iso))*age_N(iv_phyto_dino_age_tra_N(iso))
#endif
   
            ! Evolution de l azote marque des nanoflagelles  (micromoles/L N)
             ! -----------------------------------------------------------------
             dc(iv_phyto_nano_tra_N(iso))=rationnano*c(iv_phyto_nano_N)*  &
                         (fractionnh4nano*signature_N(iv_nutr_NH4_tra_N(iso))+   &
                          fractionno3nano*signature_N(iv_nutr_NO3_tra_N(iso)))  &
                         -(nanomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_nano_tra_N(iso)) &
                              -broumicrozoonano*c(iv_zoo_micr_N)*signature_N(iv_phyto_nano_tra_N(iso))

#ifdef key_age_tracer  
             dc(iv_phyto_nano_age_tra_N(iso))=c(iv_phyto_nano_tra_N(iso))+rationnano*c(iv_phyto_nano_N)*  &
                       (fractionnh4nano*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))   &
                       +fractionno3nano*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))  &
                       -(nanomorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_nano_age_tra_N(iso))   &
                       -broumicrozoonano*c(iv_zoo_micr_N)*signature_N(iv_phyto_nano_tra_N(iso))*age_N(iv_phyto_nano_age_tra_N(iso))
#endif
      
#ifdef key_karenia

             ! Evolution de l azote marque de Karenia mikimotoi (cellules/L)
             ! ---------------------------------------------------------------   
             dc(iv_phyto_karenia_tra_N(iso))= uptakeN*c(iv_phyto_karenia_C)*  &
                                       (fractionnh4karenia*signature_N(iv_nutr_NH4_tra_N(iso))  &
                                       +fractionno3karenia*signature_N(iv_nutr_NO3_tra_N(iso)))  &
                                       -(kareniamorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_karenia_tra_N(iso))  &
                                        -broumesozookarenia*c(iv_zoo_meso_N)*signature_N(iv_phyto_karenia_tra_N(iso))
#ifdef key_age_tracer
             dc(iv_phyto_karenia_age_tra_N(iso))=c(iv_phyto_karenia_tra_N(iso)) +uptakeN*c(iv_phyto_karenia_C)* &
                      (fractionnh4karenia*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))   &
                      +fractionno3karenia*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))  &
                      -(kareniamorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_karenia_age_tra_N(iso))   &
                      -broumesozookarenia*c(iv_zoo_meso_N)*signature_N(iv_phyto_karenia_tra_N(iso))  &
                                         *age_N(iv_phyto_karenia_age_tra_N(iso))
#endif
#endif

#ifdef key_phaeocystis

             ! Evolution de l azote marque des colonies de Phaeocystis globosa (micromoles/L N)
             ! ----------------------------------------------------------------------------------
             dc(iv_phyto_phaeocystis_colo_tra_N(iso))=rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)* &
                     (fractionnh4phaeocystiscolo*signature_N(iv_nutr_NH4_tra_N(iso))  & 
                     +fractionno3phaeocystiscolo*signature_N(iv_nutr_NO3_tra_N(iso)))- &
                     (phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_tra_N(iso))   &
                     -broumesozoophaeocolo*c(iv_zoo_meso_N)*signature_N(iv_phyto_phaeocystis_colo_tra_N(iso))+ &
                      initcolonie*c(iv_phyto_phaeocystis_cell_tra_N(iso))
#ifdef key_age_tracer
             dc(iv_phyto_phaeocystis_colo_age_tra_N(iso))=c(iv_phyto_phaeocystis_colo_tra_N(iso))  &
                      +rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N)* &
                       (fractionnh4phaeocystiscolo*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))  & 
                       +fractionno3phaeocystiscolo*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))  &
                      -(phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_age_tra_N(iso))   &
                      -broumesozoophaeocolo*c(iv_zoo_meso_N)*signature_N(iv_phyto_phaeocystis_colo_tra_N(iso))   &
                   *age_N(iv_phyto_phaeocystis_colo_age_tra_N(iso))+initcolonie*c(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif

             ! Evolution de l azote marque des cellules libres de Phaeocystis globosa (micromoles/L N)
             ! -----------------------------------------------------------------------------------------
             dc(iv_phyto_phaeocystis_cell_tra_N(iso))=rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N)*  &
                         (fractionnh4phaeocystiscell*signature_N(iv_nutr_NH4_tra_N(iso))    & 
                         +fractionno3phaeocystiscell*signature_N(iv_nutr_NO3_tra_N(iso)))-  &
                         (phaeocystismortcell+txfiltbenthij*ibbenth/epn)*c(iv_phyto_phaeocystis_cell_tra_N(iso))   &
                         -broumicrozoophaeocell*c(iv_zoo_micr_N)*signature_N(iv_phyto_phaeocystis_cell_tra_N(iso)) &
                          -initcolonie*c(iv_phyto_phaeocystis_cell_tra_N(iso))
#ifdef key_age_tracer
            dc(iv_phyto_phaeocystis_cell_age_tra_N(iso))=c(iv_phyto_phaeocystis_cell_tra_N(iso))  &
                           +rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N)*  &
                          (fractionnh4phaeocystiscell*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))  & 
                          +fractionno3phaeocystiscell*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso))) &
                      -(phaeocystismortcell+txfiltbenthij*ibbenth/epn)*c(iv_phyto_phaeocystis_cell_age_tra_N(iso))   &
                      -broumicrozoophaeocell*c(iv_zoo_micr_N)*signature_N(iv_phyto_phaeocystis_cell_tra_N(iso))    &
                        *age_N(iv_phyto_phaeocystis_cell_age_tra_N(iso))-initcolonie*c(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif
#endif

#ifdef key_psnz

             ! Evolution de l azote marque des Pseudo-nitzschia (micromoles/L N)
             ! -------------------------------------------------------------------
             dc(iv_phyto_psnz_tra_N(iso))=rationpsnz*c(iv_phyto_psnz_N)*  &
                                    (fractionnh4psnz*signature_N(iv_nutr_NH4_tra_N(iso))+  &
                                     fractionno3psnz*signature_N(iv_nutr_NO3_tra_N(iso)))  &
                                    -(psnzmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_psnz_tra_N(iso))   &
                                     -(broumicrozoopsnz*c(iv_zoo_micr_N)+broumesozoopsnz*c(iv_zoo_meso_N))  &
                                                 *signature_N(iv_phyto_psnz_tra_N(iso))
#ifdef key_age_tracer
             dc(iv_phyto_psnz_age_tra_N(iso))=c(iv_phyto_psnz_tra_N(iso)) +rationpsnz*c(iv_phyto_psnz_N)* &
                       (fractionnh4psnz*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))  &
                        +fractionno3psnz*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))  &
                      -(psnzmorteau+txfiltbenthij*ibbenth/epn)*c(iv_phyto_psnz_age_tra_N(iso))  &
                      -(broumicrozoopsnz*c(iv_zoo_micr_N)+broumesozoopsnz*c(iv_zoo_meso_N))  &
                                       *signature_N(iv_phyto_psnz_tra_N(iso))*age_N(iv_phyto_psnz_age_tra_N(iso))
#endif
#endif

#ifdef key_ulvas
             ! Evolution de l azote marque des ulves en suspension [gN.m-3]
             ! --------------------------------------------------------------
             ! le 0.014 permet de passer de micromolN.l-1 en mgN.l-1 (ou gN.m-3)
             dc(iv_ulv_tra_N(iso))=0.014_rsh*pompageazoteulve*c(iv_ulv_susdrywght)*  &
                            (fractionnh4*signature_N(iv_nutr_NH4_tra_N(iso))  &
                            +fractionno3*signature_N(iv_nutr_NO3_tra_N(iso)))    &
                              -ulvemorteau*c(iv_ulv_tra_N(iso))+(ulveresuspazote*signature_N(iv_ulv_benth_tra_N(iso)) &
                              -vitessechuteulve*c(iv_ulv_tra_N(iso)))/epn*ibbenth
#ifdef key_age_tracer
             dc(iv_ulv_sus_age_tra_N(iso))=c(iv_ulv_sus_tra_N(iso)) +0.014_rsh*pompageazoteulve*c(iv_ulv_susdrywght)  &
                     *(fractionnh4*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))   &
                      +fractionno3*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso)))    &
                      -ulvemorteau*c(iv_ulv_sus_age_tra_N(iso))  &
                      +(ulveresuspazote*signature_N(iv_ulv_benth_tra_N(iso))*age_N(iv_ulv_benth_age_tra_N(iso)) &
                      -vitessechuteulve*c(iv_ulv_sus_age_tra_N(iso)))/epn*ibbenth
#endif

             ! Evolution de l azote marque des ulves en depot [gN.m-2]
             ! ---------------------------------------------------------
             ! le 0.014 permet de passer de micromolN.l-1 en mgN.l-1 (ou gN.m-3)
             dc(iv_ulv_benth_tra_N(iso))=(0.014_rsh*pompageazoteulvebenth*c(iv_ulv_benthdrywght)*partulvessurface* &
                      (fractionnh4*signature_N(iv_nutr_NH4_tra_N(iso))    &
                      +fractionno3*signature_N(iv_nutr_NO3_tra_N(iso)))          &
                      -ulvebenthmortsed*c(iv_ulv_benth_tra_N(iso))  &
                      -ulveresuspazote*signature_N(iv_ulv_benth_tra_N(iso))+vitessechuteulve*c(iv_ulv_sus_tra_N(iso)))*ibbenth  
#ifdef key_age_tracer
             dc(iv_ulv_benth_age_tra_N(iso))=c(iv_ulv_benth_tra_N(iso))   &
                    +(0.014_rsh*pompageazoteulvebenth*c(iv_ulv_benthdrywght)*partulvessurface* &
                         (fractionnh4*signature_N(iv_nutr_NH4_tra_N(iso))*age_N(iv_nutr_NH4_age_tra_N(iso))  &
                         +fractionno3*signature_N(iv_nutr_NO3_tra_N(iso))*age_N(iv_nutr_NO3_age_tra_N(iso))) &
                        -ulvebenthmortsed*c(iv_ulv_benth_age_tra_N(iso))  &
                         -ulveresuspazote*signature_N(iv_ulv_benth_tra_N(iso))*age_N(iv_ulv_benth_age_tra_N(iso)) &
                     +vitessechuteulve*c(iv_ulv_sus_age_tra_N(iso)))*ibbenth 
#endif
#endif

             ! Evolution de l azote marque du microzooplancton (micromoles/L N)
             ! ------------------------------------------------------------------  
             dc(iv_zoo_micr_tra_N(iso))=assimilmicrozoo*c(iv_zoo_micr_N)*(broumicrozoonano*signature_N(iv_phyto_nano_tra_N(iso))  &
                                                                      +broumicrozoodiat*signature_N(iv_phyto_diat_tra_N(iso))  &
                                                                     + broumicrozoodet*signature_N(iv_detr_tra_N(iso)))        &
                                   -(excretionmicrozoo+txmortmicrozoo)*c(iv_zoo_micr_tra_N(iso))   &
                                    - broumesozoomicrozoo*c(iv_zoo_meso_N)*signature_N(iv_zoo_micr_tra_N(iso))

#ifdef key_phaeocystis
             dc(iv_zoo_micr_tra_N(iso))=dc(iv_zoo_micr_tra_N(iso)) +  &
                          assimilmicrozoo*c(iv_zoo_micr_N)*broumicrozoophaeocell*signature_N(iv_phyto_phaeocystis_cell_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_micr_tra_N(iso))=dc(iv_zoo_micr_tra_N(iso)) + &
                          assimilmicrozoo*c(iv_zoo_micr_N)*broumicrozoopsnz*signature_N(iv_phyto_psnz_tra_N(iso))
#endif

#ifdef key_age_tracer
             dc(iv_zoo_micr_age_tra_N(iso))=c(iv_zoo_micr_tra_N(iso)) +assimilmicrozoo*c(iv_zoo_micr_N)* &
                       (broumicrozoonano*signature_N(iv_phyto_nano_tra_N(iso))*age_N(iv_phyto_nano_age_tra_N(iso))   &
                        +broumicrozoodiat*signature_N(iv_phyto_diat_tra_N(iso))*age_N(iv_phyto_diat_age_tra_N(iso))  &
                        +broumicrozoodet*signature_N(iv_detr_tra_N(iso))*age_N(iv_detr_age_tra_N(iso)))   &
                       -(excretionmicrozoo+txmortmicrozoo)*c(iv_zoo_micr_age_tra_N(iso))  &
                       - broumesozoomicrozoo*c(iv_zoo_meso_N)*signature_N(iv_zoo_micr_tra_N(iso))*age_N(iv_zoo_micr_age_tra_N(iso))
#ifdef key_phaeocystis
             dc(iv_zoo_micr_age_tra_N(iso))=dc(iv_zoo_micr_age_tra_N(iso)) &
                        + assimilmicrozoo*c(iv_zoo_micr_N)*broumicrozoophaeocell  &
                          *signature_N(iv_phyto_phaeocystis_cell_tra_N(iso))*age_N(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_micr_age_tra_N(iso))=dc(iv_zoo_micr_age_tra_N(iso))  &
                   + assimilmicrozoo*c(iv_zoo_micr_N)*broumicrozoopsnz &
                          *signature_N(iv_phyto_psnz_tra_N(iso))*age_N(iv_phyto_psnz_age_tra_N(iso))
#endif
#endif

             ! Evolution de l azote marque du mesozooplancton (micromoles/L N)
             ! -----------------------------------------------------------------
             dc(iv_zoo_meso_tra_N(iso))=assimilmesozoo*c(iv_zoo_meso_N)*(broumesozoodiat*signature_N(iv_phyto_diat_tra_N(iso))  &
                                                                  +broumesozoodino*signature_N(iv_phyto_dino_tra_N(iso))  &
                                                                  +broumesozoomicrozoo*signature_N(iv_zoo_micr_tra_N(iso))) &
                                     -(excretionmesozoo+txmortmesozoo)*c(iv_zoo_meso_tra_N(iso))
#ifdef key_phaeocystis
             dc(iv_zoo_meso_tra_N(iso))=dc(iv_zoo_meso_tra_N(iso))+ &
                         assimilmesozoo*c(iv_zoo_meso_N)*broumesozoophaeocolo*signature_N(iv_phyto_phaeocystis_colo_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_meso_tra_N(iso))=dc(iv_zoo_meso_tra_N(iso))+ &
                         assimilmesozoo*c(iv_zoo_meso_N)*broumesozoopsnz*signature_N(iv_phyto_psnz_tra_N(iso))
#endif
#ifdef key_karenia
             dc(iv_zoo_meso_tra_N(iso))=dc(iv_zoo_meso_tra_N(iso))+  &
                         assimilmesozoo*c(iv_zoo_meso_N)*broumesozookarenia*signature_N(iv_phyto_karenia_tra_N(iso))
#endif

#ifdef key_age_tracer
             dc(iv_zoo_meso_age_tra_N(iso))=c(iv_zoo_meso_tra_N(iso)) +assimilmesozoo*c(iv_zoo_meso_N)* &
                        (broumesozoodiat*signature_N(iv_phyto_diat_tra_N(iso))*age_N(iv_phyto_diat_age_tra_N(iso))   &
                        +broumesozoodino*signature_N(iv_phyto_dino_tra_N(iso))*age_N(iv_phyto_dino_age_tra_N(iso))  &
                        +broumesozoomicrozoo*signature_N(iv_zoo_micr_tra_N(iso))*age_N(iv_zoo_micr_age_tra_N(iso)))   &
                         -(excretionmesozoo+txmortmesozoo)*c(iv_zoo_meso_age_tra_N(iso))
#ifdef key_phaeocystis
             dc(iv_zoo_meso_age_tra_N(iso))=dc(iv_zoo_meso_age_tra_N(iso))  &
                             +assimilmesozoo*c(iv_zoo_meso_N)*broumesozoophaeocolo  &
                             *signature_N(iv_phyto_phaeocystis_colo_tra_N(iso))*age_N(iv_phyto_phaeocystis_colo_age_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_zoo_meso_age_tra_N(iso))=dc(iv_zoo_meso_age_tra_N(iso))   &
                                  +assimilmesozoo*c(iv_zoo_meso_N)*broumesozoopsnz  &
                                  *signature_N(iv_phyto_psnz_tra_N(iso))*age_N(iv_phyto_psnz_age_tra_N(iso))
#endif
#ifdef key_karenia
             dc(iv_zoo_meso_age_tra_N(iso))=dc(iv_zoo_meso_age_tra_N(iso))   &
                                     +assimilmesozoo*c(iv_zoo_meso_N)*broumesozookarenia  &
                                      *signature_N(iv_phyto_karenia_tra_N(iso))*age_N(iv_phyto_karenia_age_tra_N(iso))
#endif
#endif
   
             ! Evolution de l azote detritique marque (micromoles/L N)
             ! ---------------------------------------------------------
             dc(iv_detr_tra_N(iso))=diatmorteau*c(iv_phyto_diat_tra_N(iso))+dinomorteau*c(iv_phyto_dino_tra_N(iso))+ &
                              nanomorteau*c(iv_phyto_nano_tra_N(iso))    &           
                             -reminazdeteau*c(iv_detr_tra_N(iso))+txmortmesozoo*c(iv_zoo_meso_tra_N(iso))  &
                                                                 +txmortmicrozoo*c(iv_zoo_micr_tra_N(iso)) &               
                    +(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*  &
                             (broumesozoodiat*signature_N(iv_phyto_diat_tra_N(iso))+  &
                              broumesozoodino*signature_N(iv_phyto_dino_tra_N(iso))+   &
                              broumesozoomicrozoo*signature_N(iv_zoo_micr_tra_N(iso)))   &
                    +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)* &
                              (broumicrozoonano*signature_N(iv_phyto_nano_tra_N(iso))+  &
                               broumicrozoodiat*signature_N(iv_phyto_diat_tra_N(iso))  &
                               +broumicrozoodet*signature_N(iv_detr_tra_N(iso)))-broumicrozoodet*c(iv_zoo_micr_N) &
                                         *signature_N(iv_detr_tra_N(iso)) 

#ifdef key_phaeocystis
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso)) &
                             + (1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*broumesozoophaeocolo  &
                                                    *signature_N(iv_phyto_phaeocystis_colo_tra_N(iso))  &
                             + (1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*broumicrozoophaeocell  &
                                                    *signature_N(iv_phyto_phaeocystis_cell_tra_N(iso))  &
                              + (phaeocystislysecolo+phaeocystismortcolo)*c(iv_phyto_phaeocystis_colo_tra_N(iso))+ &
                                 phaeocystismortcell*c(iv_phyto_phaeocystis_cell_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+ &
                             ((1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*broumesozoopsnz+       &
                             +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*broumicrozoopsnz)  &
                                     *signature_N(iv_phyto_psnz_tra_N(iso))+ psnzmorteau*c(iv_phyto_psnz_tra_N(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+ &
                              (1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*broumesozookarenia  &
                                     *signature_N(iv_phyto_karenia_tra_N(iso))+ kareniamorteau*c(iv_phyto_karenia_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+ulvemorteau*c(iv_ulv_sus_tra_N(iso))/0.014_rsh
#endif
#ifdef key_benthos
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+(detresuspazote*signature_N(iv_benth_tra_N(iso)) &
                                     -vitessechutedet*c(iv_detr_tra_N(iso)))*ibbenth/epn
#endif
#if ! defined key_benthos
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+txfiltbenthij*ibbenth/epn* &
                              (c(iv_phyto_diat_tra_N(iso))+c(iv_phyto_dino_tra_N(iso))+c(iv_phyto_nano_tra_N(iso)))
#ifdef key_phaeocystis
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_phaeocystis_cell_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_tra_N(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_tra_N(iso)) 
#endif      
#ifdef key_ulvas
             dc(iv_detr_tra_N(iso))=dc(iv_detr_tra_N(iso))+ulvebenthmortsed*c(iv_ulv_benth_tra_N(iso))*ibbenth/epn/0.014_rsh
#endif
#endif

#ifdef key_age_tracer
             dc(iv_detr_age_tra_N(iso))=c(iv_detr_tra_N(iso))  &
               +diatmorteau*c(iv_phyto_diat_age_tra_N(iso))+dinomorteau*c(iv_phyto_dino_age_tra_N(iso)) &
                                                           +nanomorteau*c(iv_phyto_nano_age_tra_N(iso))    &           
               -reminazdeteau*c(iv_detr_age_tra_N(iso)) &
               + txmortmesozoo*c(iv_zoo_meso_age_tra_N(iso))+txmortmicrozoo*c(iv_zoo_micr_age_tra_N(iso))       &  
               +(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*(broumesozoodiat*signature_N(iv_phyto_diat_tra_N(iso)) &
                                                            *age_N(iv_phyto_diat_age_tra_N(iso))  &
                                                          +broumesozoodino*signature_N(iv_phyto_dino_tra_N(iso)) &
                                                            *age_N(iv_phyto_dino_age_tra_N(iso))   &
                                                          +broumesozoomicrozoo*signature_N(iv_zoo_micr_tra_N(iso)) &
                                                            *age_N(iv_zoo_micr_age_tra_N(iso)))  &
               +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*(broumicrozoonano*signature_N(iv_phyto_nano_tra_N(iso)) &
                                                            *age_N(iv_phyto_nano_age_tra_N(iso))  &
                                                           +broumicrozoodiat*signature_N(iv_phyto_diat_tra_N(iso)) &
                                                            *age_N(iv_phyto_diat_age_tra_N(iso))  &
                                                           +broumicrozoodet*signature_N(iv_detr_tra_N(iso)) &
                                                            *age_N(iv_detr_age_tra_N(iso)))              &
               -broumicrozoodet*c(iv_zoo_micr_N)*signature_N(iv_detr_tra_N(iso))*age_N(iv_detr_age_tra_N(iso))                                                     
#ifdef key_phaeocystis
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))  &
               +(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*broumesozoophaeocolo &
                        *signature_N(iv_phyto_phaeocystis_colo_tra_N(iso))*age_N(iv_phyto_phaeocystis_colo_age_tra_N(iso)) &
               +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*broumicrozoophaeocell &
                        *signature_N(iv_phyto_phaeocystis_cell_tra_N(iso))*age_N(iv_phyto_phaeocystis_cell_age_tra_N(iso)) &
               +(phaeocystislysecolo+phaeocystismortcolo)*c(iv_phyto_phaeocystis_colo_age_tra_N(iso))+  &
                 phaeocystismortcell*c(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))  &
              +((1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*broumesozoopsnz+       &
              +(1.0_rsh-assimilmicrozoo)*c(iv_zoo_micr_N)*broumicrozoopsnz)   &
                          *signature_N(iv_phyto_psnz_tra_N(iso))*age_N(iv_phyto_psnz_age_tra_N(iso))       &
              + psnzmorteau*c(iv_phyto_psnz_age_tra_N(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))  &
              +(1.0_rsh-assimilmesozoo)*c(iv_zoo_meso_N)*broumesozookarenia  &
                          *signature_N(iv_phyto_karenia_tra_N(iso))*age_N(iv_phyto_karenia_age_tra_N(iso))  &
              + kareniamorteau*c(iv_phyto_karenia_age_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))+ulvemorteau*c(iv_ulv_sus_age_tra_N(iso))/0.014_rsh
#endif

#ifdef key_benthos
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))+  &
                            (detresuspazote*signature_N(iv_benth_tra_N(iso))*age_N(iv_benth_age_tra_N(iso))  &
                                 -vitessechutedet*c(iv_detr_age_tra_N(iso)))*ibbenth/epn
#endif
#if ! defined key_benthos
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso)) +txfiltbenthij*ibbenth/epn* &
                           (c(iv_phyto_diat_age_tra_N(iso))+c(iv_phyto_dino_age_tra_N(iso))+c(iv_phyto_nano_age_tra_N(iso)))
#ifdef key_phaeocystiso
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif
#ifdef key_psnz
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_psnz_age_tra_N(iso))
#endif	 
#ifdef key_karenia
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))+txfiltbenthij*ibbenth/epn*c(iv_phyto_karenia_age_tra_N(iso)) 
#endif      
#ifdef key_ulvas
             dc(iv_detr_age_tra_N(iso))=dc(iv_detr_age_tra_N(iso))+ulvebenthmortsed*c(iv_ulv_benth_age_tra_N(iso))*ibbenth/epn/0.014_rsh
#endif
#endif
#endif     

#ifdef key_benthos
   ! Evolution de l azote organique benthique marque  (mmol.m-2.jour-1)
   ! -------------------------------------------------------------------
             dc(iv_benth_tra_N(iso))=ibbenth*(vitessechutedet*c(iv_detr_tra_N(iso))-detresuspazote*signature_N(iv_benth_tra_N(iso)) &
                    -(reminazdeteau*p_reminbenth+p_burial)*c(iv_benth_tra_N(iso))   &
                    +txfiltbenthij*(c(iv_phyto_diat_tra_N(iso))+c(iv_phyto_dino_tra_N(iso))+c(iv_phyto_nano_tra_N(iso))))
#if defined key_psnz
             dc(iv_benth_tra_N(iso))=dc(iv_benth_tra_N(iso))+txfiltbenthij*ibbenth*c(iv_phyto_psnz_tra_N(iso))
#endif
#if defined key_karenia
             dc(iv_benth_tra_N(iso))=dc(iv_benth_tra_N(iso))+txfiltbenthij*ibbenth*c(iv_phyto_karenia_tra_N(iso))
#endif
#ifdef key_phaeocystis
             dc(iv_benth_tra_N(iso))=dc(iv_benth_tra_N(iso))+txfiltbenthij*ibbenth*c(iv_phyto_phaeocystis_cell_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_benth_tra_N(iso))=dc(iv_benth_tra_N(iso))+ulvebenthmortsed*c(iv_ulv_benth_tra_N(iso))*ibbenth/0.014_rsh
#endif

#ifdef key_age_tracer
             dc(iv_benth_age_tra_N(iso))=c(iv_benth_tra_N(iso))   &
                     +ibbenth*(vitessechutedet*c(iv_detr_age_tra_N(iso))  &
                     -detresuspazote*signature_N(iv_benth_tra_N(iso))*age_N(iv_benth_age_tra_N(iso)) &
                     -(reminazdeteau*p_reminbenth+p_burial)*c(iv_benth_age_tra_N(iso))   &
                     +txfiltbenthij*(c(iv_phyto_diat_age_tra_N(iso))+c(iv_phyto_dino_age_tra_N(iso))+c(iv_phyto_nano_age_tra_N(iso))))
#if defined key_psnz
             dc(iv_benth_age_tra_N(iso))=dc(iv_benth_age_tra_N(iso))+txfiltbenthij*ibbenth*c(iv_phyto_psnz_age_tra_N(iso))
#endif
#if defined key_karenia
             dc(iv_benth_age_tra_N(iso))=dc(iv_benth_age_tra_N(iso))+txfiltbenthij*ibbenth*c(iv_phyto_karenia_age_tra_N(iso))
#endif
#ifdef key_phaeocystis
             dc(iv_benth_age_tra_N(iso))=dc(iv_benth_age_tra_N(iso))+txfiltbenthij*ibbenth*c(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif
#ifdef key_ulvas
             dc(iv_benth_age_tra_N(iso))=dc(iv_benth_age_tra_N(iso))+ulvebenthmortsed*c(iv_ulv_benth_age_tra_N(iso))*ibbenth/0.014_rsh
#endif
#endif     
#endif
   
             ! ++++++++++++++++++++++++++++++++++++++++++++++++
             !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
             ! ++++++++++++++++++++++++++++++++++++++++++++++++  
             DO ivtra=1,nb_var_tracerN 
                iv_tra=iv_tracer_N(ivtra,iso)
                iv_sign=iv_signed_N(ivtra,iso)
                id_sign=id_tracer_signN(ivtra,iso)
                diag_3d_wat(id_sign,k,i,j)=0.0_rsh
                IF (c(iv_sign).gt.seuil_N_tracer_diag) diag_3d_wat(id_sign,k,i,j)=signature_N(iv_tra)
             ENDDO
#ifdef key_age_tracer  
             id_sign=id_sign+2
#else
             id_sign=id_sign+1
#endif
             diag_3d_wat(id_sign,k,i,j)=c(iv_phyto_diat_tra_N(iso))+c(iv_phyto_dino_tra_N(iso))+c(iv_phyto_nano_tra_N(iso))
             Nphytototal=c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)

#ifdef key_karenia
             diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_karenia_tra_N(iso))
             Nphytototal=Nphytototal+c(iv_phyto_karenia_N)
#endif
#ifdef key_phaeocystis
             diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+ &
                                  c(iv_phyto_phaeocystis_colo_tra_N(iso))+c(iv_phyto_phaeocystis_cell_tra_N(iso))
             Nphytototal=Nphytototal+c(iv_phyto_phaeocystis_colo_N)+c(iv_phyto_phaeocystis_cell_N)
#endif
#ifdef key_psnz
             diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_psnz_tra_N(iso))
             Nphytototal=Nphytototal+c(iv_phyto_psnz_N)
#endif
             phytototal_tra_N=diag_3d_wat(id_sign,k,i,j)
             IF (Nphytototal.gt.seuil_N_tracer_diag) THEN
               diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)/(Nphytototal)
             ELSE
               diag_3d_wat(id_sign,k,i,j)=0.0_rsh
             ENDIF
   
#ifdef key_age_tracer
             age_max=(CURRENT_TIME-tdeb_tracerN)/86400.0_rlg
             DO ivtra=1,nb_var_tracerN 
                 iv_tra=iv_tracer_N(ivtra,iso)
                 iv_sign=iv_signed_N(ivtra,iso)
                 id_sign=id_tracer_ageN(ivtra,iso)
                 ivage=iv_age_N(ivtra,iso)
                 !IF (signature_N(iv_tra) .GT. seuil_N_age_tracer_phytot) THEN
                 IF (age_N(ivage)==0.0_rsh) THEN
                    diag_3d_wat(id_sign,k,i,j)=-1.0_rsh
                    ! age superieur a 3 ans non pris en compte
                 ELSE IF (age_N(ivage) > age_max .OR. age_N(ivage) > 1095.0_rsh) THEN
                    diag_3d_wat(id_sign,k,i,j)=-1.0_rsh
                 ELSE
                    diag_3d_wat(id_sign,k,i,j)=age_N(ivage)
                 ENDIF
             ENDDO  
             id_sign=id_sign+2
             diag_3d_wat(id_sign,k,i,j)=c(iv_phyto_diat_age_tra_N(iso))+c(iv_phyto_dino_age_tra_N(iso))+c(iv_phyto_nano_age_tra_N(iso))
 
#ifdef key_karenia
             diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_karenia_age_tra_N(iso))
#endif
#ifdef key_phaeocystis
             diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+  &
                                   c(iv_phyto_phaeocystis_colo_age_tra_N(iso))+c(iv_phyto_phaeocystis_cell_age_tra_N(iso))
#endif
#ifdef key_psnz
             diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)+c(iv_phyto_psnz_age_tra_N(iso))
#endif
             !if (phytototal_tra_N.gt.1.e-4_rsh)  THEN
             IF (diag_3d_wat(id_sign-1,k,i,j) > seuil_N_age_tracer_phytot) THEN
                diag_3d_wat(id_sign,k,i,j)=diag_3d_wat(id_sign,k,i,j)/phytototal_tra_N  
             ELSE
                diag_3d_wat(id_sign,k,i,j)= -1.0_rsh  
             ENDIF
#endif

           ENDDO  ! iso (sources)
         
!   ===================================================================================         
#endif
