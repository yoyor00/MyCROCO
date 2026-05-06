#ifdef key_oxygen

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_oxygen  ***
   !&E
   !&E ** Purpose : Calcul de la variable oxygene
   !&E              
   !&E       !  2010-02    (A.Chapelle) Original code
   !&E
   !&E      use from general modele : WIND_SPEED
   !&E      use from BIOLink  : 
   !&E                           
   !&E      use from basic bloom modele : c, epn,ibbenth,ibsurf,  effetchaleur,ratio_mgO_to_mumolN_anabol
   !&E               rationdiat,rationdino ,rationnano ,rationpsnz ,rationkarenia , rationphaeocystiscolo, rationphaeocystiscell
   !&E               effetlumierediat,effetlumieredino,effetlumierenano, ..
   !&E                ratio_mgO_to_mumolN_catabol,reminazdeteau,xnitrifeau, temper,tempabs
   !&E      OUTPUT :  diag_3d_wat(irk_diag(id_oxy_sat  + dc
   !&E---------------------------------------------------------------------
 

          !++++++++++++++++++++++++++++++++++++
          ! Processus lies a l oxygene dissous 
          !++++++++++++++++++++++++++++++++++++

          ! production et consommation d oxygene dans l eau
          ! -----------------------------------------------
          ophotos=ratio_mgO_to_mumolN_anabol*(rationdiat*c(iv_phyto_diat_N)+ &
                                         rationdino*c(iv_phyto_dino_N)+rationnano*c(iv_phyto_nano_N))

#ifdef key_psnz  
          ophotos=ophotos+ratio_mgO_to_mumolN_anabol*rationpsnz*c(iv_phyto_psnz_N)
#endif
#ifdef key_karenia  
          ophotos=ophotos+ratio_mgO_to_mumolN_anabol*rationkarenia*c(iv_phyto_karenia_N)
#endif
#ifdef key_phaeocystis  
          ophotos=ophotos+ratio_mgO_to_mumolN_anabol*(rationphaeocystiscolo*c(iv_phyto_phaeocystis_colo_N) &
                  +rationphaeocystiscell*c(iv_phyto_phaeocystis_cell_N))
#endif
#ifdef key_ulvas
          ! le 0.44 permet de passer de gPoidsSec.l-1 en gC.l-1 
          ! le 0.001*32/12 permet de passer de gC.l-1 en mgO2.l-1
          ophotos=ophotos+0.44*0.001*32/12(rationulve*c(iv_ulv_susdrywght)+ &
                                    rationulvebenth*c(iv_ulv_susdrywght)*ibbenth*partulvessurface/epn)
#endif     
#ifdef key_zostera
          IF(c(iv_zost_LB)>0.0_rsh) THEN
            ophotos = ophotos + Gross_Prod_zost*c(iv_zost_LB)/epn
          ENDIF
#endif

           oresphyto=0.0_rsh
           orespzoo=0.0_rsh
           IF (c(iv_oxygen).gt.0.2_rsh) then  
              oresphyto=p_phyto_resp*((1.0_rsh-effetlumierediat)*c(iv_phyto_diat_N) &
                         + (1.0_rsh-effetlumieredino)*c(iv_phyto_dino_N)  &
                         +(1.0_rsh-effetlumierenano)*c(iv_phyto_nano_N))
#ifdef key_psnz
              oresphyto=oresphyto+  &
                         p_phyto_resp*(1.0_rsh-effetlumierepsnz)*c(iv_phyto_psnz_N)
#endif
#ifdef key_karenia
              oresphyto=oresphyto+  &
                         p_phyto_resp*(1.0_rsh-effetlumierekarenia)*c(iv_phyto_karenia_N)

#endif
#ifdef key_phaeocystis
              oresphyto=oresphyto+ p_phyto_resp*(1.0_rsh-effetlumierephaeocystis)* &
                        (c(iv_phyto_phaeocystis_colo_N)+ c(iv_phyto_phaeocystis_cell_N))
#endif 
              oresphyto=oresphyto*ratio_mgO_to_mumolN_catabol*effetchaleur
              orespzoo=p_zoo_resp*(c(iv_zoo_meso_N)+ &
                            c(iv_zoo_micr_N))*ratio_mgO_to_mumolN_catabol*effetchaleur
          ENDIF
          oremineau=reminazdeteau*c(iv_detr_N)*ratio_mgO_to_mumolN_catabol
#ifdef key_benthos   
          oremineau=oremineau+(reminazdeteau*p_reminbenth)*c(iv_benth_N)*ratio_mgO_to_mumolN_catabol*ibbenth/epn
#endif  
          onitrifeau=xnitrifeau*c(iv_nutr_NH4)*0.064_rsh
          perteo2=oresphyto+orespzoo+oremineau+onitrifeau
#ifdef key_zostera
          IF(c(iv_zost_LB)>0.0_rsh) THEN
            perteo2 = perteo2 + Leaf_Resp_zost*c(iv_zost_LB)/epn + Root_Resp_zost*c(iv_zost_RB)/epn
          ENDIF
#endif
    
          ! reoxygenation
          ! -------------
          gst=-173.4292_rsh+249.6329_rsh*100.0_rsh/tempabs+143.3483_rsh*log(tempabs/100.0_rsh)    &
                -21.8492_rsh*tempabs/100.0_rsh+sali*(-0.033096_rsh+0.014259_rsh*tempabs/100.0_rsh &
                -0.0017_rsh*(tempabs/100.0_rsh)**2)
          o2sat=1.429_rsh*exp(gst)

          !Formule de Riley et Stephan (1988)
          oair=(0.64_rsh+0.0256_rsh*(WIND_SPEED(i,j)/0.447_rsh)**2)*(o2sat-c(iv_oxygen))/epn

          !Formule de Wanninkhof (1992)
          !schmidt=(1800.6_rsh-120.1_rsh*temper+3.7818_rsh*temper**2-0.047608_rsh*temper**3)*(1.0_rsh+0.00314_rsh*sali)
          !le facteur 0.31 cm/h est converti en 24*0.31/100=0.0744 m/jour
          !oair=0.0744_rsh*(WIND_SPEED(i,j)**2)/sqrt(schmidt/660.0_rsh)*(o2sat-c(iv_oxygen))/epn

          !Formule de Wanninkhof & McGillis (1999)
          !le facteur 0.0283 cm/h est converti en 24*0.0283/100=0.00679 m/jour
          !oair=0.00679_rsh*(WIND_SPEED(i,j)**3)/sqrt(schmidt/660.0_rsh)*(o2sat-c(iv_oxygen))/epn

          ! Evolution de l oxygene dissous
          ! ------------------------------
          dc(iv_oxygen)=ophotos+oair*ibsurf-perteo2

          ! pourcentage de saturation de l oxygene
          diag_3d_wat(irk_diag(id_oxy_sat),k,i,j)=100.0_rsh*c(iv_oxygen)/o2sat

          ! Evolution du nitrate (si denitrification)
          ! ------------------------------------------
          !  IF((c(iv_oxygen)/o2sat).le.0.05_rsh) THEN
          !  dc(iv_nutr_NO3)=dc(iv_nutr_NO3)-((-0.2_rsh/o2sat)*c(iv_oxygen)+0.2_rsh)*c(iv_nutr_NO3)*ibbenth
          !  END IF
   

   !!======================================================================
#endif

