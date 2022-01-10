#if ! defined key_BLOOM_opt2

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_settling  ***
   !&E
   !&E ** Purpose : estimate settling velocity for biologic variable (detrital..)
   !&E              a l interieur d une maille d eau
   !&E
   !&E       !  2010-03    (B. Thouvenin issued fro bloomdynzwat ) Original code
   !&E       !  2019-07    (B.Thouvenin ) update
   !&E
   !&E      use from general modele : temper,sali, WATER_DENSITY  si ca existe, (h0 mais on pourrait mettre htot) 
   !&E      use from BIOLink  :
   !&E                           
   !&E      use from basic bloom modele : c, epn, l_SNeffect_settle, effetselnutdiat, l_phyzoodeteffect_settle,
   !&E                                    diatmorteau, dinomorteau, nanomorteau, rationmicrozoo, assimilmicrozoo, txmortmicrozoo
   !&E                                    rationmesozoo, assimilmesozoo, txmortmesozoo, rationmicrozoo, assimilmicrozoo, txmortmicrozoo
   !&E      use from options bloom modules : effetselnutpsnz, effetselnutphaeocystiscolo, effetlumierephaeocystis
   !&E                                       psnzmorteau, kareniamorteau, phaeocystismortcolo, phaeocystislysecolo, phaeocystismortcell
   !&E      OUTPUT : ws3
   !&E---------------------------------------------------------------------


         IF(l_SNeffect_settle) THEN

           !  Phytoplancton
           !===============
           ws3(k,iv_phyto_diat_N,i,j)=ws_free_max(iv_phyto_diat_N)
!#ifdef key_physadaptation
!             ! Modulation par la profondeur (grosses diats cotieres, petites oceaniques)
!              attention test multiplie par 0 donc non operant  donc mis en commentaire
!              ws3(k,iv_phyto_diat_N,i,j)=max(ws_free_min(iv_phyto_diat_N), ws_free_max(iv_phyto_diat_N)*  &
!                      !!!!    (1.0+0.0*min(1.0_rsh,max(0.0_rsh,(5000.0_rsh-h0(i,j))/4900.0_rsh))))
!                          (1.0+0.0*min(1.0_rsh,max(0.0_rsh,(5000.0_rsh-htot(i,j))/4900.0_rsh))))
!#endif
               effetchutedia=MAX(0.0_rsh,effetselnutdiat)**0.2_rsh
               ws3(k,iv_phyto_diat_N,i,j)=(ws_free_min(iv_phyto_diat_N)*effetchutedia+      &
                                  ws3(k,iv_phyto_diat_N,i,j)*(1.0_rsh-effetchutedia))

#ifdef key_psnz
               effetchutedia=MAX(0.0_rsh,effetselnutpsnz)**0.2_rsh
               ws3(k,iv_phyto_psnz_N,i,j)=(ws_free_min(iv_phyto_psnz_N)*effetchutedia+      &
                                  ws_free_max(iv_phyto_psnz_N)*(1.0_rsh-effetchutedia)) 
               ws3(k,iv_phyto_psnz_Si,i,j)=ws3(k,iv_phyto_psnz_N,i,j)
               ws3(k,iv_phyto_psnz_ad,i,j)=ws3(k,iv_phyto_psnz_N,i,j)
#endif  

#ifdef key_phaeocystis
               ws3(k,iv_phyto_phaeocystis_colo_N,i,j)=0.00001_rsh*(1.0_rsh+ &
                                                     5.0_rsh*(1.0_rsh-effetselnutphaeocystiscolo))
               !Ajout d une flottabilite positive en presence de lumiere
               ws3(k,iv_phyto_phaeocystis_colo_N,i,j)=ws3(k,iv_phyto_phaeocystis_colo_N,i,j)  &
                                                      +0.00005_rsh*effetlumierephaeocystis
               ws3(k,iv_phyto_phaeocystis_mucus,i,j)=ws3(k,iv_phyto_phaeocystis_colo_N,i,j)
#endif                            

          ENDIF

#ifdef key_physadaptation
          modulationsedpardenseau=1.0_rsh
          ! Modulation de la vitesse de chute du materiel particulaire en 
          ! fonction de la densite de l eau de mer
          ! -------------------------------------------------------------
          ! la loi de Stokes permet de calculer la vitesse limite de chute dans une eau a T,S
          ! en fonction de la vitesse de chute connue dans une eau a Tref,Sref par:
          ! v(T,S)=v(Tref,Sref) * [1 + (1 - ro(T,S)/ro(Tref,Sref))/(roparticule-ro(Tref,Sref)]
          ! roparticule=masse volumique de la particule (=1.05 kg.m-3 pour le phytoplancton)
          ! ro(Tref,Sref)=1.026 kg.m-3 pour l eau de mer a 15 degres C et 35 pour mille
          ! density of pure water at temperature temper
#if ! defined key_MARS
          ! IF(l_waterdensity_known) THEN
          !       modulationsedpardenseau=1.0_rsh+(1.0_rsh-(rost/ro3515))/(1.05_rsh-ro3515)abs(WATER_DENSITY_kij-WATER_DENSITY_kp1ij)/((epn+thicklayerW_C(k+1,i,j))/2.0_rsh)
          !      modulationsedpardenseau=1.0_rsh+(1.0_rsh-WATER_DENSITY_kij/1025.97275386538_rsh)/(1.05_rsh-1025.97275386538_rsh)
          ! ELSE
#endif
 
          ! row = 999.842594_rlg + (  6.793952d-02 + ( - 9.095290d-03 +                &
          !         (1.001685d-04 - 1.120083d-06*temper + 6.536332d-09*temper*temper)*temper )*temper  )*temper
          ! density of sea water at salinity sali and temperature temper
          !rost = row + sali * (   0.824493_rlg + (  -4.0899d-03 + ( 7.6438d-05 +                   &
          !           (-8.2467d-07+5.3875d-09*temper)*temper )*temper)*temper) +             &
          !           sali**1.5_rlg * ( -5.72466d-03 + (1.0227d-04-1.6546d-06*temper)*temper) +  &
          !           sali*sali*4.8314d-04
          ! modulationsedpardenseau=1.0_rsh+(1.0_rsh-REAL(rost/1025.97275386538_rlg,rsh))/(1.05_rsh-1025.97275386538_rsh)
#if ! defined key_MARS
          ! ENDIF
#endif
          !ws3(k,iv_phyto_diat_N,i,j)=ws3(k,iv_phyto_diat_N,i,j)*modulationsedpardenseau
          !!ws3(k,iv_phyto_dino_N,i,j)=ws3(k,iv_phyto_dino_N,i,j)*modulationsedpardenseau

          !IF (modulationsedpardenseau .ne.1.0_rsh) THEN

#ifdef key_psnz
          !   ws3(k,iv_phyto_psnz_N,i,j)=ws3(k,iv_phyto_psnz_N,i,j)*modulationsedpardenseau
          !   ws3(k,iv_phyto_psnz_Si,i,j)=ws3(k,iv_phyto_psnz_N,i,j)
#endif	
#ifdef key_karenia	
          !   ws3(k,iv_phyto_karenia_N,i,j)=ws3(k,iv_phyto_karenia_N,i,j)*modulationsedpardenseau
          !   ws3(k,iv_phyto_karenia_C,i,j)=ws3(k,iv_phyto_karenia_N,i,j)
          !   ws3(k,iv_phyto_karenia_P,i,j)=ws3(k,iv_phyto_karenia_N,i,j)
#endif	
#ifdef key_phaeocystis
          !   ws3(k,iv_phyto_phaeocystis_colo_N,i,j)=ws3(k,iv_phyto_phaeocystis_colo_N,i,j)*modulationsedpardenseau
          !   ws3(k,iv_phyto_phaeocystis_mucus,i,j)=ws3(k,iv_phyto_phaeocystis_colo_N,i,j)
#endif

          !ENDIF  ! fin if modulationsedpardenseau different de 1
#endif

          ! Materiel detritique
          !======================
          ! modulation en fonction du mesozooplancton
          IF(l_phyzoodeteffect_settle .and. c(iv_zoo_meso_N)>0.0_rsh) THEN

           ! PHYTODETRITUS
            phytodetritus=diatmorteau*c(iv_phyto_diat_N)  &
                     +dinomorteau*c(iv_phyto_dino_N)  &
                     +nanomorteau*c(iv_phyto_nano_N)  &
            +(rationmicrozoo*(1.0_rsh-assimilmicrozoo)+txmortmicrozoo)*c(iv_zoo_micr_N)

#ifdef key_psnz
            phytodetritus=phytodetritus+psnzmorteau*c(iv_phyto_psnz_N)
#endif
#ifdef key_karenia
            phytodetritus=phytodetritus+kareniamorteau*c(iv_phyto_karenia_N)
#endif
#ifdef key_phaeocystis
            phytodetritus=phytodetritus+(phaeocystismortcolo+phaeocystislysecolo)*c(iv_phyto_phaeocystis_colo_N) &
                         +c(iv_phyto_phaeocystis_cell_N)*phaeocystismortcell
#endif

 
            ! ZOODETRITUS

            zoodetritus=(rationmesozoo*(1.0_rsh-assimilmesozoo)+txmortmesozoo)*c(iv_zoo_meso_N)   &
                  +(rationmicrozoo*(1.0_rsh-assimilmicrozoo)+txmortmicrozoo)*c(iv_zoo_micr_N)
            IF (zoodetritus > 0.0_rsh) THEN
               rdet=phytodetritus/(zoodetritus+epsilon_BIOLink)       
               wchutedet=p_detzoo_wsed*(1.0_rsh/(rdet+1.0_rsh))+          &
                      p_detphy_wsed*(1.0_rsh-(1.0_rsh/(rdet+1.0_rsh)))
            ELSE
               wchutedet=p_detphy_wsed
            ENDIF
          ELSE
            wchutedet=p_detphy_wsed
          ENDIF

#ifdef key_MANGAbio
          ! Acceleration en cours de descente par formation de neige
          IF (k.lt.kmax) THEN
            abovesinkingrate=0.0_rsh
            DO k=k+1,kmax
              abovesinkingrate=abovesinkingrate+diag_3d_wat(irk_diag(id_detsettling),k,i,j)
            ENDDO
            abovesinkingrate=max((abovesinkingrate/(kmax-k)/86400.0_rsh),wchutedet)
            wchutedet=0.1_rsh*wchutedet+0.9_rsh*abovesinkingrate
          ENDIF
#endif
          ws3(k,iv_detr_N,i,j)=wchutedet
          ws3(k,iv_detr_Si,i,j)=wchutedet
          ws3(k,iv_detr_P,i,j)=wchutedet

#if defined key_N_tracer
          DO iso=1,nb_source_tracerN
            DO ivtra=1,nb_var_tracerN 
             iv_tra=iv_tracer_N(ivtra,iso)
             IF(iv_tra <= nvp) THEN
                iv_sign=iv_signed_N(ivtra,iso)
                ws3(k,iv_tra,i,j)=ws3(k,iv_sign,i,j)
#ifdef key_age_tracer
                ivage=iv_age_N(ivtra,iso)
                ws3(k,ivage,i,j)=ws3(k,iv_sign,i,j)
#endif
             ENDIF
            ENDDO
          ENDDO
#endif
#if defined key_P_tracer
          DO iso=1,nb_source_tracerP
            DO ivtra=1,nb_var_tracerP 
             iv_tra=iv_tracer_P(ivtra,iso)
             IF(iv_tra <= nvp) THEN
                iv_sign=iv_signed_P(ivtra,iso)
                ws3(k,iv_tra,i,j)=ws3(k,iv_sign,i,j)
#ifdef key_age_tracer
                ivage=iv_age_P(ivtra,iso)
                ws3(k,ivage,i,j)=ws3(k,iv_sign,i,j)
#endif
             ENDIF
            ENDDO
          ENDDO
#endif
  
          ! bornage de la vitesse de chute des variables particulaires en chaque maille
          ! mise a zero aux limites ouvertes
          ! ---------------------------------------------------------------------------
#ifdef key_MANGAbio
          IF ((l_obc_west .AND. im < imin+obc_width+1) .OR. (l_obc_east .AND. im > imax-obc_width) .OR. &
              (l_obc_south .AND. jm < jmin+obc_width+1) .OR. (l_obc_north .AND. jm > jmax-obc_width)) THEN
               ws3(k,:,i,j)=0.0_rsh
          ELSE
#endif
               DO ivp=1,nvp
                  ws3(k,ivp,i,j)=sign(MIN(0.95_rlg*epn/dtbio,REAL(ABS(ws3(k,ivp,i,j)),rlg)),ws3(k,ivp,i,j))
               END DO
#ifdef key_MANGAbio
          ENDIF
#endif

          !!======================================================================
#endif 

 
