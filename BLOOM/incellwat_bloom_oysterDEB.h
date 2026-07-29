  
#if defined key_oyster_DEB

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_oysterDEB  ***
   !&E
   !&E ** Purpose : Calcul des variations dues a la croissance des filtreurs
   !&E              selon le modèle DEB
   !&E              
    !&E       !          (S. Pouvreau, K. Grangere) Original code
   !&E       
   !&E      use from general modele : CELL_SURF  
   !&E      use from general modele : 
   !&E      use from BIOLink  : ,dtbiojour, temper, extinction, PAR_top_layer, htot
   !&E                           cmes_3dmgl,diag_3d_wat(irk_diag(id_totalchl)
   !&E      use from basic bloom modele : c, epn,  effetchaleur, nbhuitre (in initdefine), fact_phyto_ChlNratio
   !&E                                    effetselnutdiat
   !&E      OUTPUT :   dc, diag_2d(irk_diag(id_oys_DEB....
   !&E---------------------------------------------------------------------
 
          IF (ibbenth.eq.1) THEN 
          
             ! Parametres des entrees
             ! -----------------------
             !  espece=espece consideree
             !  chloro=concentration en chlorophylle dans l eau (microg/l)
             !  temp=temperature de l eau (degrC)
             !  NRJ_reserves=Quantitee d energie des reserves (variable d etat Joules)
             !  NRJ_structures=Quantitee d energie des structures (variable d etat Joules)
             !  NRJ_gonades=Quantitee d energie des gonades (variable d etat Joules)
             !
             ! Parametres des sorties
             ! ----------------------
             !  pa=taux d assimilation (J j-1)
             !  pc=flux cataboliques (J j-1)
             !  pg=flux vers la croissance (J j-1)
             !  amaig=taux d ammaigrissement (J j-1)
             !  pr=taux de croissance des goandes (J j-1)
             !  ponte=flux vers la ponte (J j-1)
             !  lyse=lyse des gonades (J j-1)
             !  DWres=poids sec des reserves
             !  DWstruct=poids sec des structures
             !  DWgon=poids sec des gonades
             !  DWtot=poids sec total
   
             tempref=20.

             !!SR utilisation du C comme unite de nourriture
             !        p_oysDEB_kchl=150
             !  Utilisation d une relation entre p_oysDEB_kchl et chloro moyenne fond (Desclaux et al. 2014)
             !         p_oysDEB_kchl=104.091*chlorofondmoy(i,j)+52.26
             !         p_oysDEB_kchl=0.89*(104.091*chlorofondmoy(i,j)+52.26)+0.09*(74.42*tempfondmax(i,j)-1084.54)


             IF(nbhuitre(i,j) .ne. 0.0 .and. htot(i,j) > hautable(i,j)) THEN
    
                untier=1.0_rsh/3.0_rsh
                deuxtier=2.0_rsh/3.0_rsh

                ! calcul de la chlorophylle  
                chloro=diag_3d_wat(irk_diag(id_totalchl),k,i,j)/fact_phyto_ChlNratio
                MES=cmes_3dmgl(k,i,j)/1000_rsh

                ! Calcul du rapport Chloro/carbone
                ChlC=0.003+0.0154*exp(0.05*temper)*(exp(-0.059*            &
                     (((PAR_top_layer(k,i,j)/1951000)*86400)/(EXTINCTION_RAD(k,i,j)  &
                      *epn))*(1-exp(-(EXTINCTION_RAD(k,i,j)*epn)))))*effetselnutdiat
     
                !!SR utilisation du C comme unite de nourriture
                chloro=chloro/ChlC
      
                !!SR    NRJ_reservesH=c(iv_oysdeb_res)
                !!    NRJ_structuresH=c(iv_oysdeb_str)
                !!    NRJ_gonadesH=c(iv_oysdeb_gon)
   
                NRJ_reserves=max(0.,c(iv_oysdeb_res))
                NRJ_structures=max(0.,c(iv_oysdeb_str))
                NRJ_gonades=max(0.,c(iv_oysdeb_gon))

                !SR verification de transfert
                ! print*,'DEB1',i,j,c(iv_oysdeb_res),c(iv_oysdeb_str),c(iv_oysdeb_gon)

                ! Maximum surface ingestion rate {Pxm}
                !! Correction temperature avec formulation Y. Bourles 2009 : ingestion
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                !!SR       pxm = p_oysDEB_pxm0 * exp(p_oysDEB_TA*((1./(273+20))-(1./(273+temp))))/   &
                !!SR       (1+exp(p_oysDEB_TAL*(1/(273+temp)-1/p_oysDEB_TL))+exp(p_oysDEB_TAH*(1/p_oysDEB_THing-1/     &
                !!SR       (273+temp))))    
                facttemp=exp(p_oysDEB_TA*((1./(273+tempref))-(1./(273+temper)))) *   &
                           (1+exp(p_oysDEB_TAL*(1/(273+tempref)-1/p_oysDEB_TL))+exp(p_oysDEB_TAH*(1/p_oysDEB_THing-1/     &
                               (273+tempref)))) /                                           &
                           (1+exp(p_oysDEB_TAL*(1/(273+temper)-1/p_oysDEB_TL))+exp(p_oysDEB_TAH*(1/p_oysDEB_THing-1/     &
                               (273+temper))))

                pxm = p_oysDEB_pxm0 * facttemp

                ! Maximum surface assimilation rate {Pam}
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                pam = p_oysDEB_ae * pxm

                ! Reponse fonctionnelle
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                ! xk variable avec MES (ATTENTION MES en g par litre) !!!!!!!!!
                p_oysDEB_kchl=p_oysDEB_kchl*(1+(MES/0.06))

                ! forcage par chla mesuree - p_oysDEB_kchl=coef de demi saturation
                fchl = chloro/(chloro+p_oysDEB_kchl) 
       
                vol=(NRJ_structures/p_oysDEB_Eg)+epsilon_BIOLink

                ! Assimilation rate Pa
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                 pa = pam * fchl * (vol**deuxtier)       
                !
                ! Volume cost of maintenance [Pm]
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                !! Correction temperature avec formulation Y. Bourles 2009 : respiration
                !!SR        pm = p_oysDEB_pm0 * exp(p_oysDEB_TA*((1./(273+20))-(1./(273+temp))))/   &
                !!       (1+exp(p_oysDEB_TAL*(1/(273+temp)-1/p_oysDEB_TL))+exp(p_oysDEB_TAH*(1/p_oysDEB_THresp-1/   &
                !!       (273+temp))))

                facttemp=exp(p_oysDEB_TA*((1./(273+tempref))-(1./(273+temper))))   &
                    *(1+exp(p_oysDEB_TAL*(1/(273+tempref)-1/p_oysDEB_TL))+exp(p_oysDEB_TAH*(1/p_oysDEB_THresp-1/     &
                    (273+tempref))))                                            &
                    /(1+exp(p_oysDEB_TAL*(1/(273+temper)-1/p_oysDEB_TL))+exp(p_oysDEB_TAH*(1/p_oysDEB_THresp-1/     &
                    (273+temper))))

                pm = p_oysDEB_pm0 * facttemp
     
                !  Catabolic flux Pc
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                pc = ((p_oysDEB_Eg/p_oysDEB_Em)*pam*(vol**(-untier)) + pm) /  &
                           ((p_oysDEB_Kappa/vol)+(p_oysDEB_Eg/(NRJ_reserves+epsilon)))
      
                !   Flux vers la croissance
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                pg = p_oysDEB_Kappa*pc - pm*vol
                IF (pg .le. 0)  pg=0.
        
                !   Maturite et reproduction
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                pj = min((p_oysDEB_Vp*((1-p_oysDEB_Kappa)/p_oysDEB_Kappa)*pm),(vol*((1-p_oysDEB_Kappa)/p_oysDEB_Kappa)*pm))
                
                IF (vol .gt. p_oysDEB_Vp) THEN
                   ! Maintenance developpement
                   ! Developpement - repro
                   pr=(1-p_oysDEB_Kappa) * pc - pj
                ELSE
                   pr = 0.
                ENDIF
                IF(pr.lt.0.) THEN
                   IF (NRJ_gonades .le. 0.) THEN
                       pr=0.
                   ELSE
                       pr=-min(-pr,NRJ_gonades/dtbiojour)
                   ENDIF
                ENDIF
       
                !  amaigrissement
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                seuil = (vol**untier)*pm/(p_oysDEB_Kappa*pam)
                energ = (NRJ_reserves/vol)/p_oysDEB_Em
                !       if ((energ .le. seuil) .and. (NRJ_gonades .eq. 0.)) then
                   IF ((energ .le. seuil) .and. (NRJ_gonades .le. 0.)) THEN
                           amaig = -p_oysDEB_Kappa*pc + pm*vol
                           amaigmax=max(0.,NRJ_structures/dtbiojour)
                           amaig=min(amaig, amaigmax)
                   ELSE
                      amaig = 0.0_rsh
                      amaigmax = 0.0_rsh
                   ENDIF
       
                !  lyse
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                !       if (energ .le. seuil) then
                IF ((energ .le. seuil).and. (NRJ_gonades.gt.0.)) THEN
                   lyse =  ((pm*vol)-p_oysDEB_Kappa*pc)/p_oysDEB_kR
                   lysemax=max(0.,NRJ_gonades/dtbiojour)
                   lyse = min(lyse,lysemax)
                ELSE
                   lyse = 0.0_rsh
                   lysemax = 0.0_rsh
                ENDIF
                !
                !   variables etats poids sec
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                DWres = NRJ_reserves/p_oysDEB_muE
                DWstruct = 0.15*vol
                DWgon = p_oysDEB_kR*NRJ_gonades/p_oysDEB_muE
                DWtot = DWres + DWstruct + DWgon
                Lgtot=vol**untier/p_oysDEB_shape

                !  sauv taux filt
                txfilt=pxm * (vol**deuxtier)/p_oysDEB_kchl/24           &
                !SR uptake  /((4.189*11.4)/ChlC/1000.)/(DWtot+epsilon)
                       /((4.189*11.4)/1000.)/(DWtot+epsilon)
                txing=pxm * (vol**deuxtier)/(chloro+p_oysDEB_kchl)/24   &
                !SR uptake  /((4.189*11.4)/ChlC/1000.)/(DWtot+epsilon)  
                       /((4.189*11.4)/1000.)/(DWtot+epsilon) 

                !   ponte - huitre - ER=rapport gonado-somatique (%)
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                ER = 100.*(p_oysDEB_kR*NRJ_gonades/p_oysDEB_muE)/(DWtot+epsilon)
                !!SR        if ((ER .gt. p_oysDEB_ERlim) .and. (temp .gt. p_oysDEB_Tseuilponte)) then
                IF ((ER .gt. p_oysDEB_ERlim) .and. (temper .gt. p_oysDEB_Tseuilponte)) THEN
                   ponte=NRJ_gonades/dtbiojour
                ELSE
                   ponte=0.0_rsh
                ENDIF

                !!!! pour huitre creuse

                !  Filtration population
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                FR= pxm * (chloro/p_oysDEB_kchl)*(vol**deuxtier)
                FRpop=FR*nbhuitre(i,j)/1000/epn/CELL_SURF(i,j)/p_phyto_ChlNratio/((4.189*11.4)/ChlC/1000.)        

                !  pseudo-feces
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                PF=FR- pxm * fchl * (vol**deuxtier)
                PFpop=PF*nbhuitre(i,j)/1000./CELL_SURF(i,j)/p_phyto_ChlNratio/((4.189*11.4)/ChlC/1000.)

                ! ingestion 
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!       px=filtration * pxm * fchl * (vol**deuxtier)
                px= pxm * fchl * (vol**deuxtier)
                pxpop=px*nbhuitre(i,j)/1000./epn/CELL_SURF(i,j)/p_phyto_ChlNratio/2./((4.189*11.4)/ChlC/1000.)
 
                !  Feces
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                Feces= pxm * fchl * (vol**deuxtier)-pa
                Fpop=Feces*nbhuitre(i,j)/1000./CELL_SURF(i,j)/p_phyto_ChlNratio/((4.189*11.4)/ChlC/1000.)

                 !!SR!!       if (im==75.and.jm==25) print*,'FRpop, PFpop, Fpop', FRpop, PFpop, Fpop


                !!!!!EVOLUTION DES VARIABLES D ETAT
                !!!=====================================

                !!SR       if ((nbhui1(i,j).ne.0.0_rsh).and.(htot(i,j).gt.0.5_rsh)) then

                !  Evolution de l energie des reserves
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                dc(iv_oysdeb_res)=(pa-pc)
                !	  write(*,*)'dc(iv_oysdeb_res)=',i,j,dc(iv_oysdeb_res)

                !  Evolution de l energie des structures
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                dc(iv_oysdeb_str)=(pg-amaig)
                !	  write(*,*)'dc(iv_oysdeb_str)=',i,j,dc(iv_oysdeb_str)

                !  Evolution de l energie des gonades
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                dc(iv_oysdeb_gon)=(pr-ponte-lyse)
               
                !	  if((NRJ_gonades.le.0).and.(dc(iv_oysdeb_gon).lt.0)) dc(iv_oysdeb_gon)=0.
                !	  if((NRJ_gonades.gt.0).and.(dc(iv_oysdeb_gon).lt.0)) then
                !	     dc(iv_oysdeb_gon)=-min(-dc(iv_oysdeb_gon),NRJ_gonades/dtbiojour)
                !	  endif
                !	  write(*,*)'dc(iv_oysdeb_gon)=',i,j,dc(iv_oysdeb_gon)
                !  if(ibbenth.ne.0) print*,'DCGNOADES',dc(iv_oysdeb_gon), pr,ponte, lyse,c(iv_oysdeb_gon),c(iv_oysdeb_res),c(iv_oysdeb_str)
                !
                !!!!!!VARIABLES DIAGNOSTIQUES
                !!!=====================================
     
                diag_2d(irk_diag(id_oysDEB_FRpop),i,j)= FRpop !FILTRATION HUITRE
                diag_2d(irk_diag(id_oysDEB_PFpop),i,j)= (PFpop/epn) !PSEUDOFECES HUITRE
                diag_2d(irk_diag(id_oysDEB_Fpop),i,j)= (Fpop/epn) !FECES HUITRE
                diag_2d(irk_diag(id_oysDEB_DWtot),i,j)= DWtot ! poids sec total HUITRE
                diag_2d(irk_diag(id_oysDEB_kchlvble),i,j)= p_oysDEB_kchl !coefficient de demi-saturation chloro variable
                !   diag_2d(irk_diag(id_oysDEB_kchlvble),k,i,j)= k_chlvble !coefficient de demi-saturation chloro variable
                diag_2d(irk_diag(id_oysDEB_ER),i,j)= ER  !rapport GONADO-SOMATIQUE HUITRE
                diag_2d(irk_diag(id_oysDEB_filtrate),i,j)= txfilt  !taux de filtration HUITRE
                diag_2d(irk_diag(id_oysDEB_Wtot),i,j)= (DWtot/0.09/0.2) ! poids total
                diag_2d(irk_diag(id_oysDEB_Istress),i,j)= ((vol*pm+pj)/(pa+epsilon)) !INDICE DE STRESS
                diag_2d(irk_diag(id_oysDEB_Lgtot),i,j)= Lgtot  !Longueur totale HUITRE
                diag_2d(irk_diag(id_oysDEB_ChlC),i,j)= ChlC  !rapport carbone/chloro
                diag_2d(irk_diag(id_oysDEB_nbhuit),i,j)= nbhuitre(i,j)  !nombre d huitre


                !+++++++++++++++++++++++++++++++++++++++++++++++
                ! EQUATIONS D EVOLUTION  ! concentration en jour-1
                !+++++++++++++++++++++++++++++++++++++++++++++++
                   
                ! Evolution de l azote des diatomees
                ! ----------------------------------
                dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)   &
                   -FRpop*c(iv_phyto_diat_N)/(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)+1e-10)

                !!!!prise ne compte d une production microphytobenthique
                !!if (htot(i,j) > 1.5_rsh .and. htot(i,j) < 2.0_rsh .and. j.ge.33 .and. couvmicrophyto(i,j)==1) then
                !!dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)+350.0_rsh
                !!endif
                !!if (htot(i,j) > 1.5_rsh .and. htot(i,j)< 2.0_rsh .and. j.lt. 33 .and. couvmicrophyto(i,j)==1) then
                !!dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)+200.0_rsh
                !!endif
                !!!if (htot(i,j)> 1.5_rsh .and. htot(i,j)<2.0_rsh .and. j.lt.32 .and. i.lt.57 .and. couvchloro(i,j)==1) then
                !!if (htot(i,j)>1.5_rsh .and. htot(i,j)<2.0_rsh .and. j.lt.24 .and. i.lt.118 .and. couvmicrophyto(i,j).gt.0.0_rsh) then
                !!dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)+10.0_rsh
                !!endif

   
                ! Evolution de l azote des dinoflagelles
                ! --------------------------------------
                dc(iv_phyto_dino_N)=dc(iv_phyto_dino_N)   &
                    -FRpop*c(iv_phyto_dino_N)/(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)+1e-10)

   
                ! Evolution de l azote detritique benthique
                ! --------------------------------------
                !dc(iv_detr_N)=dc(iv_detr_N)+(egestotNhui/epn-consoNdet2)
                !!    dc(iv_detr_N)=dc(iv_detr_N)+(PFpop+Fpop)   !&  !!!!SR
                dc(iv_benth_N)=dc(iv_benth_N)+(PFpop+Fpop)
                  !!    attention aux unites!!!!!
                  !!              -reminazdeteau*0.1*c(iv_detr_N)


                ! Evolution du phosphore detritique benthique
                ! ---------------------------------
                 !dc(iv_detr_P)=dc(iv_detr_P)+(egestotPhui/epn-consoPdet2) 
                 !!   dc(iv_detr_P)=dc(iv_detr_P)+(PFpop+Fpop)/16.0_rsh
                dc(iv_benth_P)=dc(iv_benth_P)+(PFpop+Fpop)/16.0_rsh
                         
                ! Evolution de la silice detritique benthique
                ! ---------------------------------
                dc(iv_benth_Si)=dc(iv_benth_Si)     &
                    +rapsiaz*FRpop*c(iv_phyto_diat_N)/(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)+1e-10)*epn


                ! Evolution de l azote du nano_pico_phytoplancton  (micromolN.l-1)
                ! -----------------------------------------------------------------
                dc(iv_phyto_nano_N)=dc(iv_phyto_nano_N)  &
                    -FRpop*c(iv_phyto_nano_N)/(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)+1e-10)
    
                ! Evolution de l ammonium
                ! ----------------------------------------------------------------- 
                !! excrethuit=1./1000.*24.*PDS*nbhuitre/volmaille
    
                !   dc(11)=reminazdeteau*c(7)-fractionnh4*(rationdiat*c(5)+
                !  $  rationdino*c(6))-xnitrifeau*c(11)+
                ! $  excretionzoo*(c(12)/12.)*(bioparam(58)/bioparam(59))
                ! $  +fluxbenthN(i,j)/epn
                ! $  +reminazdeteau*0.1*c(16)/epn
                !c     $  +reminazdeteau*c(16)/epn
                !     &  +(excretmoul+excrethuit+excretcrep)

            ELSE  ! nbhuitre > 0 et hautable
                dc(iv_oysdeb_res)=0.0_rsh
                dc(iv_oysdeb_str)=0.0_rsh
                dc(iv_oysdeb_gon)=0.0_rsh

            ENDIF

          ENDIF    ! ibbenth=1
   !!======================================================================
#endif
