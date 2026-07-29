  
#if defined key_oyster_SFG

   !&E---------------------------------------------------------------------
   !&E                 ***  incellwat_bloom_oysterSFG  ***
   !&E
   !&E ** Purpose : Calcul des variations dues a la croissance de l huitre creuse
   !&E              
   !&E       !          (L. Barille) Original code
   !&E       
   !&E      use from general modele : CELL_SURF  
   !&E      use from general modele : 
   !&E      use from BIOLink  : ,dtbio, jjulien_BIOLink, ,, temper
   !&E                           imois_BIOLink,iheure_BIOLink,iminu_BIOLink,isec_BIOLink,cmes_3dmgl,diag_3d_wat(irk_diag(id_totalchl)
   !&E      use from basic bloom modele : c, epn,  effetchaleur, nbhuitre (in initdefine)
   !&E      OUTPUT :   dc, diag_2d(irk_diag(id_oys_filt, ,id_oys_orgafring, id_oys_absorg, id_oys_resp
   !&E---------------------------------------------------------------------
 
          IF (ibbenth.eq.1) THEN

            tjour=jjulien_BIOLink+iheure_BIOLink/24.0_rlg+iminu_BIOLink/1440.0_rlg+isec_BIOLink/86400.0_rlg

            ! Filtration - filt2 (L.h-1.ind-1):
            ! Fonction du poids de l individu, de la temperature et des MES (effet colmatage)

            IF(nbhuitre(i,j) .ne. 0) THEN
                volmaille=epn*CELL_SURF(i,j)

                susp=cmes_3dmgl(k,i,j)
                colmat=min(0.0_rsh,(p_oys_colmat_thr-susp))

                effettemper=p_oys_paramfiltr*(temper-p_oys_tempopt_filtr)**2
                IF (susp > p_oys_sestfiltr_thr) THEN
                     filtmes=-p_oys_kfilt*susp+p_oys_yfilt
                ELSE
                     filtmes=p_oys_stdfiltr
                ENDIF
                filt=filtmes-effettemper
                IF (filt < 0.0_rsh)  filt=0.0_rsh
        
                filt1=filt * c(iv_oys_so)**p_oys_allfiltr
                filt2=filt1*exp(p_oys_colmparam*colmat)   !!!*p_oys_subratio

                IF(htot(i,j) <= 0.6_rsh)  filt2=0.0_rsh

                diag_2d(irk_diag(id_oys_filt),i,j)=filt2

                !!          concentration en chlorophylle  en µg/L
                !          cphyto=(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))*p_phyto_ChlNratio
                !         cphyto=(abs(cphyto)+cphyto)*0.5
                ! Essai avec la valeur de chlorophylle calculee precedemment
                cphyto=diag_3d_wat(irk_diag(id_totalchl),k,i,j)  !! ATTENTION c est exprime en chlorophylle

                !!        c(iv_detr_N) en umolN/L --> µmol de carbone --> cdetri en mg/L

                cdetri=(c(iv_detr_N)/16+c(iv_detr_P))*106.0_rsh*12.0_rsh/1000.0_rsh/p_oys_CtoDWratio
                cdetri=(abs(cdetri)+cdetri)*0.5
                IF (susp <= p_oys_retmin_thr) THEN
                      retention=-p_oys_kreten*susp+p_oys_yreten
                ELSE
                      retention=p_oys_retmin
                ENDIF

               !! concentration en MOP  en mg/L          
               IF(imois_BIOLink .GE. 1 .AND. imois_BIOLink .LT. 6) THEN
                   cmop=(0.088_rsh*cphyto+0.52_rsh)*1.1_rsh
               ENDIF
        
               IF(imois_BIOLink .GE. 6 .AND. imois_BIOLink <= 12) THEN
                  cmop=(0.088_rsh*cphyto+0.52_rsh)*1.1_rsh
               ENDIF
              
               !! Consommation phyto - en ugChla.h-1.ind-1          =umolNdiat.h-1.ind-1 :
               !         consochla=filt2*cphyto 

               !!! Consommation MOP - en mg.h-1
               consorg=filt2*cmop

        
               !! Consommation detritique et inorganique en mg.h-1.ind-1
               consodet=filt2*cdetri*retention
               consoi=filt2*susp*retention
        

               !!ce qui est consomme disparait du systeme --> cf dc(iv_phyto_diat_N) et dc(iv_detr_N)
               !! consochla en umolNdiat.h-1.ind-1 --> consochla2 en j-1
               consochla2=filt2*nbhuitre(i,j)*24.0_rsh/1000.0_rsh/volmaille

               !! consodet en mg.h-1.ind-1 --> consodet2 en mg.j-1
               consodet2=consodet*nbhuitre(i,j)*24
        
               !! consoNdet2 et consoPdet2 en µmol.l-1.jour-1           
               consoNdet2=consodet2*p_oys_CtoDWratio /12.0_rsh*16.0_rsh/106.0_rsh/volmaille
               consoPdet2=consodet2*p_oys_CtoDWratio /12.0_rsh/106.0_rsh/volmaille           

               !! Efficacite de selection de la Chla, du detritique et de l inorganique
               val1=min(0.0_rsh,(p_oys_prodpseufec_thr-susp))

               PForg=p_oys_prodpseufeco_lev*(1.0_rsh-exp(p_oys_pseufeco_exp*val1))

               PFi=p_oys_prodpseufeci_lev*(1.0_rsh-exp(p_oys_pseufeci_exp*val1))
        
               IF(PForg > 1.0_rsh)   PForg=1.0_rsh
               IF(PForg < 0.0_rsh)   PForg=0.0_rsh
               IF(PFi > 1.0_rsh)  PFi=1.0_rsh
               IF(PFi < 0.0_rsh)  PFi=0.0_rsh


               !!Ingestion individuelle de la Chla, du detrique et de l inorganique
               cingesorg=(1.0_rsh-PForg)*consorg
               cingesi=(1.0_rsh-PFi)*consoi
               cingestot=cingesorg+cingesi

               !! Production de PF :
               !! pseuchla en ugChla.ind-1.h-1 = umolNdiat.h-1.ind-1
               !! pseudet en mg.ind-1.h-1
               pseuorg=max(0.0_rsh,consorg-cingesorg)
               pseudi=max(0.0_rsh,consoi-cingesi)
        
               !! Absorption de la partie ingeree et production de feces :
               !!feceschla en ugChla.ind-1.h-1 = umolNdiat.h-1.ind-1
               !!fecesdet en mg.ind-1.h-1
               !!! ORGAFRING entre 0 et 100 
               ! print*,'cingesorg cingestot',cingesorg,cingestot       
               ORGAFRING=((cingesorg)/(cingestot+epsilon))*100
               diag_2d(irk_diag(id_oys_orgafring),i,j)=ORGAFRING
        
               !!! EAorg entre 0 et 1
               EAorg=p_oys_absmaxmop*(1.0_rsh-exp(-p_oys_absdet*ORGAFRING))
               IF(EAorg < 0.0_rsh)   EAorg=0.0_rsh
               IF(EAorg > 1.0_rsh)   EAorg=1.0_rsh
        
               !!!Energie absorbee
               ABSorg1=(EAorg*cingesorg/1000.0_rsh)!!*p_oys_subratio

               !!!energie totale absorbee en g MOP.jour-1   
               ABSorg=ABSorg1*24.0_rsh!!!*p_oys_subratio
        
               diag_2d(irk_diag(id_oys_absorg),i,j)=ABSorg
        
               !!!!production de feces en mg MOP.h-1.ind-1
               !!! note Benedicte : dans la doc de Morgan (rapport Gerricao, il y a ici un rapport 1000 comme pour ABSorg ??????????????????????????????
               fecesorg=((1.0_rsh-EAorg)*cingesorg)!!*p_oys_subratio
        
               !! Egestion totale de N detritique :
               !! pseuorg et fecesorg en mg.ind-1.h-1 --> psfdet2 en mmolNdet.j-1.m-2
               psfdet2=(pseuorg+fecesorg)*nbhuitre(i,j)*24.0_rsh
        
               psfNdet2=psfdet2*p_oys_CtoDWratio /12.0_rsh*16.0_rsh/106.0_rsh/CELL_SURF(i,j)
               psfPdet2=psfdet2*p_oys_CtoDWratio /12.0_rsh/106.0_rsh/CELL_SURF(i,j)
        
               egestotNhui=psfNdet2
               egestotPhui=psfPdet2


               !!!!!!RESPIRATION!!!!!!!!!!

               !! Reformulation de la loi de Bougrier en g de MOP.jour-1

               RESP1=(0.2_rsh*exp(0.1_rsh*temper))*(c(iv_oys_so)**p_oys_allres)*0.7_rsh/1000.0_rsh!!*p_oys_subratio/1000.0_rsh
               RESP=RESP1*24.0_rsh!!!*p_oys_subratio

               diag_2d(irk_diag(id_oys_resp),i,j)=RESP

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!REPARTITION DE L ENERGIE ET CROISSANCE!!!!!!!
               !!initialisation de REPOSCOQUILLE ET REPOSGAMETE
               IF (tpostpontecoq(i,j)==1.0_rlg)  REPOSCOQUILLE=0
               IF (tpostpontegam(i,j)==1.0_rlg)  REPOSGAMETE=0
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               !!!!!pas en periode post-ponte!!!!!!!!!!!!!!!!!!
               IF (tpostpontecoq(i,j) < 0.0_rlg)  REPOSCOQUILLE=0
               IF (tpostpontegam(i,j) < 0.0_rlg)  REPOSGAMETE=0
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               !!!!!periode post-ponte de 15 jours durant laquelle la vitellogenese est impossible
               !!!!!ARRET post-ponte de la croissance coquillere pendant 7 jours

               IF (tpostpontegam(i,j) > 0.0_rlg .AND. (tjour-tpostpontegam(i,j)) <= 15.0_rlg) THEN
                     REPOSGAMETE=1
               ENDIF
               IF (tpostpontecoq(i,j) > 0.0_rlg .AND. (tjour-tpostpontecoq(i,j)) <= 7.0_rlg) THEN
                   REPOSCOQUILLE=1
               ENDIF
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               IF (tpostpontecoq(i,j) > 0.0_rlg .AND. (tjour-tpostpontecoq(i,j)) > 7.0_rlg) THEN
                   REPOSCOQUILLE=0
                   tpostpontecoq(i,j)=-1.0_rlg
               ENDIF

               IF (tpostpontegam(i,j) > 0.0_rlg .AND. (tjour-tpostpontegam(i,j)) > 15.0_rlg) THEN
                  REPOSGAMETE=0
                  tpostpontegam(i,j)=-1.0_rlg
               ENDIF

               repsom=c(iv_oys_go)/(c(iv_oys_so)+epsilon)

               !!!GAMETO=2 periode de vitellogenese
               !!!GAMETO=1 hors periode de vitellogenese
               IF(temper > 8.0_rsh .AND. repsom < 0.6_rsh .AND. REPOSGAMETE==0) THEN
                  GAMETO=2
               ELSE
                  GAMETO=1
               ENDIF

               !!!DECLENCHEMENT DE LA PONTE
               IF(temper > 19_rsh .AND. repsom > 0.35_rsh .AND. REPOSGAMETE==0) THEN
                  SEUILPONTE=1
                  tpostpontecoq(i,j)=tjour
                  tpostpontegam(i,j)=tjour
               ELSE
                  SEUILPONTE=0
               ENDIF

               !IF (ibbenth==1) THEN
                  diag_2d(irk_diag(id_oys_tcoq),i,j)=tpostpontecoq(i,j)
                  diag_2d(irk_diag(id_oys_tgam),i,j)=tpostpontegam(i,j)
               !ENDIF

               !!!PONTE
               IF (SEUILPONTE==1) THEN
                  PONTE=(0.6_rsh*c(iv_oys_go)*3600.0_rsh*24.0_rsh)/dtbio
               ELSE
                  PONTE=0.0_rsh
               ENDIF

               !IF (ibbenth==1) THEN
                  diag_2d(irk_diag(id_oys_ponte),i,j)=diag_2d(irk_diag(id_oys_ponte),i,j)+PONTE*dtbio/(3600.0_rsh*24.0_rsh)
               !ENDIF

               !!!REPARTITION ENERGIE!!!

               !!!!!!!!!!!!!!!!!!!
               !!! COQUILLE 
               !!!!!!!!!!!!!!!!!!!
               IF (REPOSCOQUILLE == 1) THEN
                 COQGAM=0.0_rsh
               ELSE
                 COQGAM=1/GAMETO
               ENDIF

               IF ((ABSORG-RESP) > 0.0_rsh) THEN
                 ORGCOQUILLE=(0.05_rsh*ABSORG)*COQGAM
               ELSE
                 ORGCOQUILLE=0.0_rsh
               ENDIF

               GAINCOQUILLE=ORGCOQUILLE*100.0_rsh

               dc(iv_oys_co)=GAINCOQUILLE

               !!!!!!!!!!!!!!!!!!!
               !!!   SOMA
               !!!!!!!!!!!!!!!!!!
               BILANCOQ=ABSORG-ORGCOQUILLE

               IF (c(iv_oys_so) < 1.5_rsh) THEN
                 MAXSOMA=0.0293_rsh*c(iv_oys_so)
               ELSE
                 MAXSOMA=0.044_rsh
               ENDIF
               BILANSOMA=min(MAXSOMA,BILANCOQ)

               IF (GAMETO==2 .OR. repsom < 0.2_rsh) THEN
                 GAINSOMA=BILANSOMA-RESP
               ELSE IF (GAMETO==1) THEN
                 GAINSOMA=BILANSOMA
               ENDIF


               !!!!!!!!!!!!!!!!!!
               !!! REPROD
               !!!!!!!!!!!!!!!!!!!

               IF (GAMETO==2 .OR. repsom < 0.2_rsh) THEN
                 GAINREPROD=BILANCOQ-BILANSOMA
               ELSE IF (GAMETO==1) THEN
                 GAINREPROD=BILANCOQ-(BILANSOMA+RESP)
               ENDIF

               dc(iv_oys_so)=GAINSOMA
               dc(iv_oys_go)=(GAINREPROD-PONTE)

               !!!!!!!!!!!!
               diag_2d(irk_diag(id_oys_gaincoq),i,j)=GAINCOQUILLE
               diag_2d(irk_diag(id_oys_gainsoma),i,j)=GAINSOMA
               diag_2d(irk_diag(id_oys_gainrepr),i,j)=GAINREPROD
               diag_2d(irk_diag(id_oys_bilansoma),i,j)=BILANSOMA
               diag_2d(irk_diag(id_oys_bilancoq),i,j)=BILANCOQ

               diag_2d(irk_diag(id_oys_gameto),i,j)=GAMETO
               diag_2d(irk_diag(id_oys_seuilponte),i,j)=SEUILPONTE
               diag_2d(irk_diag(id_oys_reposcoq),i,j)=REPOSCOQUILLE
               diag_2d(irk_diag(id_oys_reposgam),i,j)=REPOSGAMETE



               !+++++++++++++++++++++++++++++++++++++++++++++++
               ! EQUATIONS D EVOLUTION  ! concentration en jour-1
               !+++++++++++++++++++++++++++++++++++++++++++++++

               ! Evolution de l azote des diatomees
               ! ----------------------------------
               dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)-consochla2*c(iv_phyto_diat_N)


               ! Evolution de l azote des dinoflagelles
               ! --------------------------------------
               dc(iv_phyto_dino_N)=dc(iv_phyto_dino_N)-consochla2*c(iv_phyto_dino_N)
      
               ! Evolution de l azote detritique
               ! --------------------------------------
               dc(iv_detr_N)=dc(iv_detr_N)+(egestotNhui/epn-consoNdet2)

               ! Evolution du phosphore detritique
               ! ---------------------------------
               dc(iv_detr_P)=dc(iv_detr_P)+(egestotPhui/epn-consoPdet2)
             
               ! Evolution de l azote du nano_pico_phytoplancton  (micromolN.l-1)
               ! -----------------------------------------------------------------
               dc(iv_phyto_nano_N)=dc(iv_phyto_nano_N)-consochla2*c(iv_phyto_nano_N)


            ENDIF  ! nbhuitre


          ENDIF    ! ibbenth=1
   !!======================================================================
#endif
