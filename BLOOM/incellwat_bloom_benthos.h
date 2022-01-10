#if defined key_benthos
   

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom benthos  ***
   !&E
   !&E ** Purpose : Calcul du stockage du detritique dans 3 ou 4 variables benthiques
   !&E              
   !&E       !  2010-02    (A. Menesguen, M. Dussauze) Original code
   !&E
   !&E      use from general modele : htot,uz,vz,  
   !&E      use from general modele : si key_MANGAbio  ubr_vague,vbr_vague, LATITUDE,LONGITUDE
   !&E      use from BIOLink  : ,dtbio, dtbiojour, ws3
   !&E                           
   !&E      use from basic bloom modele : c, epn,epnsed, rapsiaz, rappaz, effetchaleur,txfiltbenthij
   !&E      OUTPUT :   dc
   !&E---------------------------------------------------------------------
  
           IF (ibbenth.eq.1) THEN

          !++++++++++++++++++++++++++++++++++++
          ! Processus lies au benthos 
          !++++++++++++++++++++++++++++++++++++
#if defined key_MARS
            IF(htot(i,j).le.hm) THEN
               ufond=0.5_rsh*sqrt((u(i-1,j)+u(i,j))**2.+(v(i,j)+v(i,j-1))**2.)
               ! ifirst,ilast,jfirst,jlast sont tels qu on ne peut que avoir i>1 ou j>1
              !IF(i > 1 .and. j > 1) ufond=0.5_rsh*sqrt((u(i-1,j)+u(i,j))**2.+(v(i,j)+v(i,j-1))**2.)
              !IF(i==1 .and. j > 1) ufond=0.5_rsh*sqrt((u(i,j)+u(i,j))**2.+(v(i,j)+v(i,j-1))**2.)
              !IF(i > 1 .and. j==1) ufond=0.5_rsh*sqrt((u(i-1,j)+u(i,j))**2.+(v(i,j)+v(i,j))**2.)
              !IF(i==1 .and. j==1) ufond=sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            ELSE
#endif
              ufond=BOTTOM_CURRENT_ij
#if defined key_MARS
              ! precautions a programmer si besoin
              !IF(i > 1 .and. j > 1) ufond=0.5_rsh*sqrt((U_BOTTOM_BIOLink(i-1,j)+U_BOTTOM_BIOLink(i,j))**2.+(V_BOTTOM_BIOLink(i,j)+V_BOTTOM_BIOLink(i,j-1))**2.)
              !IF(i==1 .and. j > 1) ufond=0.5_rsh*sqrt((U_BOTTOM_BIOLink(i,j)+U_BOTTOM_BIOLink(i,j))**2.+(V_BOTTOM_BIOLink(i,j)+V_BOTTOM_BIOLink(i,j-1))**2.)
              !IF(i > 1 .and. j==1) ufond=0.5_rsh*sqrt((U_BOTTOM_BIOLink(i-1,j)+U_BOTTOM_BIOLink(i,j))**2.+(V_BOTTOM_BIOLink(i,j)+V_BOTTOM_BIOLink(i,j))**2.)
              !IF(i==1 .and. j==1) ufond=sqrt(U_BOTTOM_BIOLink(i,j)*U_BOTTOM_BIOLink(i,j)+V_BOTTOM_BIOLink(i,j)*V_BOTTOM_BIOLink(i,j))

            ENDIF
#endif
    
#if defined key_MANGAbio && defined key_MANGAbiovague
            ! specifique MANGA pour remettre en suspension le detritique benthique quand il y a des vagues
            !! pas d interpolation  pour simplifier. on prend la valeur a la date inf

            ubr_interp=ubr_vague(i,j,idateinf)
            vbr_interp=vbr_vague(i,j,idateinf)

            IF (htot(i,j).le.200.0_rsh) THEN
             ufond_vague=sqrt(ubr_interp*ubr_interp+vbr_interp*vbr_interp)
            ELSE
             ufond_vague=0.0_rsh
            END IF
            ufond=ufond+ufond_vague
#endif

            reminbenth=p_reminbenth
            erodvitcrit=p_erodvitcrit

#if defined key_MANGAbio
            !!!ESTUAIRE DE LA LOIRE ACCUMULATION ARTIFICIELLE DE DETRITIQUE PAS DE PRISE EN COMPTE
            IF(LATITUDE(i,j) .ge. 47.11 .and. LATITUDE(i,j) .le. 47.4 .and. &
             LONGITUDE(i,j) .ge. -2.23 .and. LONGITUDE(i,j) .le. -2.02) THEN
                 erodvitcrit=100.0_rsh
                 reminbenth=0.0_rsh
            ENDIF
#endif

            !detresusp=(p_erodflux*max(0.0_rsh,(ufond/erodvitcrit-1.0_rsh)))
            detresusp=(p_erodflux*max(0.0_rsh,(ufond**2/erodvitcrit**2-1.0_rsh)))


            IF (c(iv_benth_N) > 1.e-4) THEN 
              detresuspazote=detresusp
              IF (detresuspazote*dtbio > c(iv_benth_N)) THEN
                  detresuspazote=0.95_rsh*c(iv_benth_N)/dtbio
              END IF
             detresuspazote=detresuspazote*86400.0_rsh
            ELSE
             detresuspazote=0.0_rsh
            END IF
   
            IF (c(iv_benth_Si) > 1.e-4) THEN
              detresuspsil=detresusp*rapsiaz
              IF (detresuspsil*dtbio > c(iv_benth_Si)) THEN 
                 detresuspsil=0.95_rsh*c(iv_benth_Si)/dtbio
              END IF 
              detresuspsil=detresuspsil*86400.0_rsh
            ELSE
              detresuspsil=0.0_rsh
            END IF
     
            IF (c(iv_benth_P) > 1.e-4) THEN
             detresusppho=detresusp*rappaz
             IF (detresusppho*dtbio > c(iv_benth_P)) THEN
                detresusppho=0.95_rsh*c(iv_benth_P)/dtbio
             END IF     
             detresusppho=detresusppho*86400.0_rsh
            ELSE
             detresusppho=0.0_rsh
            END IF

#if defined key_diatbenth || defined key_zostera

            IF (c(iv_benth_phyto_diat_N) > 1.e-15) THEN
             diatresusp=detresusp
             IF (diatresusp*dtbio > c(iv_benth_phyto_diat_N)) THEN
               diatresusp=c(iv_benth_phyto_diat_N)/dtbio
             END IF        
             diatresusp=diatresusp*86400._rsh
            ELSE
             diatresusp=0._rsh
            END IF
#endif
#if defined key_zostera 
            IF (c(iv_zost_benth_seed) > 1.e-15) THEN
             seedresusp=detresusp
             IF (seedresusp*dtbio > c(iv_zost_benth_seed)) THEN
                seedresusp=c(iv_zost_benth_seed)/dtbio
             END IF        
             seedresusp=seedresusp*86400._rsh
            ELSE
             seedresusp=0._rsh
            END IF

            IF (c(iv_zost_benth_N) > 1.e-15) THEN
             detNzostresusp=detresusp
             IF (detNzostresusp*dtbio > c(iv_zost_benth_N)) THEN
               detNzostresusp=c(iv_zost_benth_N)/dtbio
             END IF        
             detNzostresusp=detNzostresusp*86400._rsh
            ELSE
             detNzostresusp=0._rsh
            END IF

            IF (c(iv_zost_benth_P) > 1.e-15) THEN
             detPzostresusp=detresusp
             IF (detPzostresusp*dtbio > c(iv_zost_benth_P)) THEN
               detPzostresusp=c(iv_zost_benth_P)/dtbio
             END IF        
             detPzostresusp=detPzostresusp*86400._rsh
            ELSE
             detPzostresusp=0._rsh
            END IF
#endif   

            !vitesse de chute en m.jour-1
            ! --------------------------------------------
            ! ********ATTENTION deux choix possibles *********
            !  : choix en vitesse au carre : (1-u^2/(p_depovitcrit)^2)
            !  : ou choix en vitesse       :  (1-u/p_depovitcrit)
            factvitessmax=max(0.,(1.0_rsh-ufond/p_depovitcrit))
           !factvitessmax=max(0.,(1.0_rsh-ufond**2./p_depovitcrit**2.))
            vitessechutedet=min(0.95_rsh*epn,dtbio*ws3(1,iv_detr_N,i,j)*factvitessmax)/dtbiojour
            vitessechutediat=min(0.95_rsh*epn,dtbio*ws3(1,iv_phyto_diat_N,i,j)*factvitessmax)/dtbiojour

#if defined key_zostera
            vitessechuteseed=min(0.95_rsh*epn,dtbio*ws3(1,iv_zost_seed,i,j)*factvitessmax)/dtbiojour
            vitessechutezostdet=min(0.95_rsh*epn,dtbio*ws3(1,iv_detr_zost_N,i,j)*factvitessmax)/dtbiojour
#endif
          
#if defined key_NPbenth || defined key_zostera
            ! diffusion a l interface eau-sed (mmole/m2/j)
            ! --------------------------------------------
            IF(epn > 0.01) THEN
              diffusion_PO4 = p_kzf/deltaz*(c(iv_benth_PO4) - c(iv_nutr_PO4))*86400.0_rsh ! positif vers la colonne d eau
              diffusion_NH4 = p_kzf/deltaz*(c(iv_benth_NH4) - c(iv_nutr_NH4))*86400.0_rsh
            ELSE
              diffusion_PO4 = 0.0_rsh
              diffusion_NH4 = 0.0_rsh
            ENDIF
#endif
 
#if defined key_diatbenth || defined key_zostera
            diatmortsed=effetchaleur*p_diat_mort*0.3_rsh

          ! Evolution de l azote organique benthique en mmol.m-2.jour-1
          ! ----------------------------------------
            dc(iv_benth_N)=(-detresuspazote+vitessechutedet*c(iv_detr_N)) &
                    -(reminazdeteau*reminbenth+p_burial)*c(iv_benth_N)   &
                    +diatmortsed*c(iv_benth_phyto_diat_N)         &
                    +txfiltbenthij*(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))
            
#else
            ! Evolution de l azote organique benthique en mmol.m-2.jour-1
            ! ----------------------------------------
            dc(iv_benth_N)=(-detresuspazote+vitessechutedet*c(iv_detr_N)+vitessechutediat*c(iv_phyto_diat_N)) &
                    -(reminazdeteau*reminbenth+p_burial)*c(iv_benth_N)   &
                    +txfiltbenthij*(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))
#endif
   
#if defined key_psnz
            dc(iv_benth_N)=dc(iv_benth_N)+txfiltbenthij*c(iv_phyto_psnz_N)
#endif
#if defined key_karenia
            dc(iv_benth_N)=dc(iv_benth_N)+txfiltbenthij*c(iv_phyto_karenia_N)
#endif
#ifdef key_phaeocystis
            dc(iv_benth_N)=dc(iv_benth_N)+txfiltbenthij*c(iv_phyto_phaeocystis_cell_N)
#endif
#ifdef key_ulvas
            dc(iv_benth_N)=dc(iv_benth_N)+ulvebenthmortsed*c(iv_ulv_benth_N)/0.014_rsh
#endif
                               
            ! Evolution du silicium organique benthique en mmol.m-2.jour-1
            ! -----------------------------------------
            dc(iv_benth_Si)=-detresuspsil+vitessechutedet*c(iv_detr_Si) &
                    -(dissolsiliceeau*reminbenth+p_burial)*c(iv_benth_Si) &
                    +rapsiaz*txfiltbenthij*c(iv_phyto_diat_N)
#if defined key_diatbenth || defined key_zostera
            dc(iv_benth_Si)=dc(iv_benth_Si)+diatmortsed*c(iv_benth_phyto_diat_N)*p_phyto_SiNratio
#else
            dc(iv_benth_Si)=dc(iv_benth_Si)+vitessechutediat*c(iv_phyto_diat_N)*rapsiaz
#endif
#if defined key_psnz
            dc(iv_benth_Si)=dc(iv_benth_Si)+txfiltbenthij*c(iv_phyto_psnz_Si)
#endif
 
            ! Evolution du phosphore organique benthique en mmol.m-2.jour-1
            ! ------------------------------------------
            dc(iv_benth_P)=(-detresusppho+vitessechutedet*c(iv_detr_P)) &
                    -(reminpdeteau*reminbenth+p_burial)*c(iv_benth_P)  &
                    +txfiltbenthij*(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))*rappaz
#if defined key_diatbenth || defined key_zostera
            dc(iv_benth_P)=dc(iv_benth_P)+diatmortsed*c(iv_benth_phyto_diat_N)/p_phyto_NPratio 
#else
            dc(iv_benth_P)=dc(iv_benth_P)+vitessechutediat*c(iv_phyto_diat_N)*rappaz
#endif
#if defined key_psnz
            dc(iv_benth_P)=dc(iv_benth_P)+txfiltbenthij*c(iv_phyto_psnz_N)*rappaz
#endif
#if defined key_karenia
            !   dc(iv_benth_P)=dc(iv_benth_P)+txfiltbenthij*c(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)
            !    dc(iv_benth_P)=dc(iv_benth_P)+txfiltbenthij*c(iv_phyto_karenia_N)*rappaz
            dc(iv_benth_P)=dc(iv_benth_P)+txfiltbenthij*c(iv_phyto_karenia_P)
#endif
#ifdef key_phaeocystis
            dc(iv_benth_P)=dc(iv_benth_P)+txfiltbenthij*c(iv_phyto_phaeocystis_cell_N)*rappaz
#endif
#ifdef key_ulvas
            dc(iv_benth_P)=dc(iv_benth_P)+ulvebenthmortsed*c(iv_ulv_benth_P)/0.031_rsh
#endif

            ! Evolution de l azote detritique en micromol.l-1.jour-1
            ! -------------------------------
            dc(iv_detr_N)=dc(iv_detr_N)-(-detresuspazote+vitessechutedet*c(iv_detr_N))/epn
  
            ! Evolution de la silice detritique en micromol.l-1.jour-1
            ! ---------------------------------
            dc(iv_detr_Si)=dc(iv_detr_Si)-(-detresuspsil+vitessechutedet*c(iv_detr_Si))/epn
   
            ! Evolution du phosphore detritique en micromol.l-1.jour-1
            ! ------------------------------
            dc(iv_detr_P)=dc(iv_detr_P)-(-detresusppho+vitessechutedet*c(iv_detr_P))/epn
      
            ! Evolution du silicium dissous
            ! ------------------------------
            dc(iv_nutr_SiOH)=dc(iv_nutr_SiOH)+(dissolsiliceeau*reminbenth)*c(iv_benth_Si)/epn
   
#if defined key_diatbenth || defined key_zostera
            ! Evolution des diatomes benthiques en mmol.m-2.jour-1
            ! ------------------------------
            dc(iv_benth_phyto_diat_N)=(-diatresusp+vitessechutediat*c(iv_phyto_diat_N))  &
                              -diatmortsed*c(iv_benth_phyto_diat_N)

            ! Evolution des diatomes pelagiques
            ! --------------------------------
            dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)-(vitessechutediat*c(iv_phyto_diat_N)-diatresusp)/epn
#else
            dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)-vitessechutediat*c(iv_phyto_diat_N)/epn
#endif
#if defined key_NPbenth
            ! Evolution du phosphore dissous benthique (mmol/m3ei/j)
            ! ----------------------------------------------------
            dc(iv_benth_PO4)=dc(iv_benth_PO4)+((reminpdeteau*p_reminbenth)*c(iv_benth_P)  &
                           -diffusion_PO4)/epnsed

            ! Evolution de l ammonium benthique (mmol/m3ei/j)
            ! -----------------------------------------------
            dc(iv_benth_NH4)=dc(iv_benth_NH4)+((reminazdeteau*p_reminbenth)*c(iv_benth_N)  &
                     -diffusion_NH4)/epnsed
    
            ! Evolution du phosphore dissous
            ! ------------------------------   
            dc(iv_nutr_PO4)=dc(iv_nutr_PO4)+diffusion_PO4/epn

            ! Evolution de l ammonium
            ! -----------------------
            dc(iv_nutr_NH4)=dc(iv_nutr_NH4)+diffusion_NH4/epn

#elif defined key_zostera
            ! Evolution de l azote detritique peu labile (mmol/m2/j)
            ! ------------------------------------------------------
            dc(iv_zost_benth_N)=dc(iv_zost_benth_N)-detNzostresusp+vitessechutezostdet*c(iv_detr_zost_N)  &
                              -(reminazdeteau*reminbenth/10._rsh)*c(iv_zost_benth_N)

            ! Evolution du phosphore detritique peu labile (mmol/m2/j)
            ! ------------------------------------------------------
            dc(iv_zost_benth_P)=dc(iv_zost_benth_P)-detPzostresusp+vitessechutezostdet*c(iv_detr_zost_P)  &
                        -(reminpdeteau*reminbenth/10._rsh)*c(iv_zost_benth_P)
            ! Evolution du phosphore dissous benthique (mmol/m3ei/j)
            ! ----------------------------------------------------
            dc(iv_benth_PO4)=dc(iv_benth_PO4)+((reminpdeteau*p_reminbenth)*c(iv_benth_P)  &
                     +(reminpdeteau*p_reminbenth/10._rsh)*c(iv_zost_benth_P)  &
                     -diffusion_PO4)/epnsed

            ! Evolution de l ammonium benthique (mmol/m3ei/j)
            ! -----------------------------------------------
            dc(iv_benth_NH4)=dc(iv_benth_NH4)+((reminazdeteau*p_reminbenth)*c(iv_benth_N)  &
                     +(reminazdeteau*p_reminbenth/10._rsh)*c(iv_zost_benth_N)  &
                     -diffusion_NH4)/epnsed
   
            ! Evolution des graines de zosteres derivantes (mmolC/m3/j)
            ! -----------------------------------------------------------
            dc(iv_zost_seed)=dc(iv_zost_seed)-(-seedresusp+vitessechuteseed*c(iv_zost_seed))/epn

            ! Evolution des graines de zosteres benthiques (mmolC/m2/j)
            ! -----------------------------------------------------------
            dc(iv_zost_benth_seed)=dc(iv_zost_benth_seed)-seedresusp+vitessechuteseed*c(iv_zost_seed)

            ! Evolution du phosphore dissous
            ! ------------------------------   
            dc(iv_nutr_PO4)=dc(iv_nutr_PO4)+diffusion_PO4/epn

            ! Evolution de l ammonium
            ! -----------------------
            dc(iv_nutr_NH4)=dc(iv_nutr_NH4)+diffusion_NH4/epn

#else
            ! Evolution du phosphore dissous
            ! ------------------------------   
            dc(iv_nutr_PO4)=dc(iv_nutr_PO4)+(reminpdeteau*reminbenth)*c(iv_benth_P)/epn
  
            ! Evolution de l ammonium
            ! ------------------------------
            dc(iv_nutr_NH4)=dc(iv_nutr_NH4)+(reminazdeteau*reminbenth)*c(iv_benth_N)/epn
#endif

#ifdef key_MANGAbio
            ! Evolution de la MES : perte definitive dans le sediment
            ! ATTENTION, ca piege toute la MES en amont des estuaires si pas de contrainte de profondeur
            ! modification par rapport a version Alain 2019 : utilisation de htot au lieu de h0 pour eviter le transfert de h0
            ! -------------------------------------------------------  
             !IF (h0(i,j).gt.15.0_rsh) dc(iv_spim)=dc(iv_spim)-amin1(0.95_rsh*epn,dtbio*ws3(1,iv_spim,i,j)*factvitessmax)/dtbiojour*c(iv_spim)/epn
             IF (htot(i,j).gt.15.0_rsh) dc(iv_spim)=dc(iv_spim)-amin1(0.95_rsh*epn,dtbio*ws3(1,iv_spim,i,j)*factvitessmax)/dtbiojour*c(iv_spim)/epn
#endif
          ENDIF    ! ibbenth=1
   !!======================================================================
#endif
