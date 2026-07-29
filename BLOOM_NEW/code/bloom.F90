  MODULE bloom

#include "cppdefs.h"

#if defined SUBSTANCE && defined BLOOM

  !!======================================================================
   !!                   ***  MODULE  BLOOM ***
   !!
   !! Biologic dynamics:  - all routines which are independant of host hydro model
   !!======================================================================

   USE module_BIOLink

#include "coupleur_define_BIOLink.h"

   !! * Modules used
   USE comBIOLink
   USE comBIOLink_physics
   USE comBIOLink_helping
   USE COMBLOOM

#if defined MUSTANG
   USE comMUSTANG
#endif
   USE comsubstance

   USE shared
   USE phytoplankton, only: phytoplankton_dynamics, storage_diag_phytoplankton, get_Tchl
   USE zooplankton, only: zooplankton_dynamics
   USE chemistry, only: chemistry_processes
   USE oxygen, only: oxygen_dynamics
   USE bio_sediment, only: sediment_bio_dynamics
   USE sinking_velocity, only: estimate_sinking_velocity, storage_diag_sinking_velocity
   USE mes , only : get_total_MES

 
   IMPLICIT NONE

   !! * Accessibility
   PUBLIC bloom_sksc_wat,bloom_eval_diag2d,bloom_SPMtot_Chla,bloom_extinction_avg          
#if defined MUSTANG &&  defined key_BLOOM_insed
   PUBLIC bloom_reactions_in_sed
#endif


 CONTAINS

   !!======================================================================
  SUBROUTINE bloom_SPMtot_Chla(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_SPMtot_Chla ***
   !&E                     - 
   !&E
   !&E ** Purpose : estimate concentrations in water column of SPMtot (MES+POM) in mg.L-1 and Chla (mgChl.L-1)
   !&E              in order to evaluate further radiation extinction
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-12 (B. Thouvenin) 
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

   !! * Declaration Subroutine

    !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast

   !! * Local declarations
    INTEGER                  ::  i,j,k,kmaxmod,iv        ! loop indexes
    
   !!----------------------------------------------------------------------
   !! * Executable part


!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
     DO i=ifirst,ilast
       IF (htot(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN
            
           ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
           kmaxmod=NB_LAYER_WAT

           DO k=1,kmaxmod


             ! Total MES in mg.l-1

             cmes_3dmgl(k,i,j)= get_total_MES(i,j,k, cvadv_wat_pos(:,k,i,j))

             !chlorophyl concentration in mug.l-1

             BIOLink_chloro(k,i,j)=get_Tchl(i,j,k,cvadv_wat_pos(:,k,i,j)) 


           ENDDO  ! loop on k       
       ENDIF ! if water
     ENDDO       
   ENDDO 
!$OMP END DO
      
  END SUBROUTINE bloom_SPMtot_Chla

   !!======================================================================
  SUBROUTINE bloom_extinction_avg(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_extinction_avg ***
   !&E                     - 
   !&E
   !&E ** Purpose : estimate extinction average on 4 days
   !&E              
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-12 (B. Thouvenin) 
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

   !! * Declaration Subroutine

    !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast
   
   !! * Local declarations
    INTEGER                  ::  i,j,k,kmaxmod       ! loop indexes
    INTEGER                  ::  iday,jhour ,numday_extinction,numhour_extinction       
    REAL(KIND=rsh)           ::  fact_phyto_ChlNratio   
   !!----------------------------------------------------------------------
   !! * Executable part

           ! ==================================================================
           ! implementation du cumul pour le calcul de l extinction moyenne horaire puis sur 4 jours
           !  t_cum_extinctionh a ete initialise dans BIOLink_alloc au debut de la simulation
           ! ==================================================================
!$OMP SINGLE
   IF(iheure_BIOLINK .ne. ihour_previous) THEN
         !write(*,*)'t_cum_extinctionh ',t_cum_extinctionh,dtbio,iminu_BIOLINK,isec_BIOLINK
         ihour_previous=iheure_BIOLINK
         t_cum_extinctionh=0.0_rsh
   ELSE
         t_cum_extinctionh=t_cum_extinctionh+dtbio
   ENDIF
!$OMP END SINGLE

   numday_extinction=mod(jjulien_BIOLINK,4)+1
   numhour_extinction=iheure_BIOLINK+1


!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
     DO i=ifirst,ilast
       IF (htot(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN
            
               ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
          kmaxmod=NB_LAYER_WAT

         !  memorization of extinction_aveh and extinction_ave4d, extinction_tab
          diag_3d_wat(irk_diag(id_extinctioncoeff),:,i,j)=EXTINCTION_RAD(:,i,j)

         IF(l_ChlNratio_var) THEN
           ! - calcul de l extinction moyenne sur 4 jours  --------------------------------------------
           DO k=LOOPK_SUBSURF_TO_BOTTOM_WAT   ! kmaxmod-1,1,-1
             extinction_aveh(k,i,j)=extinction_aveh(k,i,j)+EXTINCTION_RAD(k,i,j)*dtbio
             IF(t_cum_extinctionh >= (3600._rsh-dtbio) .AND. t_cum_extinctionh > 0.0_rsh) THEN
                extinction_aveh(k,i,j)=extinction_aveh(k,i,j)/t_cum_extinctionh
                extinction_tab(numday_extinction,numhour_extinction,k,i,j)=extinction_aveh(k,i,j)
                extinction_ave4d(k,i,j)=0.0_rsh
                DO iday=1,4
                DO jhour=1,24
                   extinction_ave4d(k,i,j)=extinction_ave4d(k,i,j)+extinction_tab(iday,jhour,k,i,j)
                END DO
                END DO
                extinction_ave4d(k,i,j)=extinction_ave4d(k,i,j)/96.0_rsh
             ENDIF
           ENDDO
           DO k=kmaxmod+1,NB_LAYER_WAT
             EXTINCTION_RAD(k,i,j)=EXTINCTION_RAD(1,i,j)
             extinction_aveh(k,i,j)=extinction_aveh(1,i,j)
             IF(t_cum_extinctionh >= (3600._rsh-dtbio)) THEN
                extinction_tab(numday_extinction,numhour_extinction,k,i,j)= &
                   extinction_tab(numday_extinction,numhour_extinction,1,i,j)
                extinction_ave4d(k,i,j)=extinction_ave4d(1,i,j)
             ENDIF
           ENDDO
        ENDIF ! if l_ChlNratio_var
       ENDIF ! if water
     ENDDO       
   ENDDO 
!$OMP END DO
      
  END SUBROUTINE bloom_extinction_avg


   !!======================================================================
   
  SUBROUTINE bloom_sksc_wat(ifirst,ilast,jfirst,jlast, WIND_SPEED )
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_sksc_wat  ***
   !&E
   !&E              bloomgic dynamics:  - estimate sink and source terms for each variables
   !&E                     - update a few variable
   !&E
   !&E ** Purpose : Estimate variations of concentrations in the water.
   !&E              Estimate new concentrations of fixed variables.
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
  !&E       !          (A. Menesguen P. Cugier, M. Sourisseau ) Original code
   !&E       !  2008-01 (M. Sourisseau) evolutions
   !&E       !  2008-05 (M. Sourisseau) link between microzoo and diatoms
   !&E       !  2012-10  (B. Thouvenin) fusion des versions bloom et mise a jour V906
   !&E       !  2013     (B.Thouvenin) mutiples modifications ane adjustments for OMP and other corrections 
   !&E       !  2013-09 (B.Thouvenin) integration of module tracer_N (A. Menesguen)
   !&E       !  2014-04 (B.Thouvenin) integration of module tracer_P (A. Menesguen) + update module Benthos
   !&E       !  2014-05 (B.Thouvenin) integration of key_zoo_prod (A. MENESGUEN)
   !&E       !  2014-05 (B.Thouvenin) update for option 1 (+ key_diatbenth) (A. MENESGUEN,  M. DUSSAUZ)
   !&E       !  2014-06 (B.Thouvenin) integration of zostere module  (+ key_zostera) (M. PLUS)
   !&E       !  2014-07 (A. Menesguen, M. Dussauz) effetturbidite get as an argument
   !&E       !  2015-03 (B. Thouvenin) integration of OysterDEB module (P. CUGIER)
   !&E       !  2015-03 (B. Thouvenin) integration of OysterSFG module (modele Scope for Growth de Barille, 1997) (A. MENESGUEN)
   !&E       !  2019-02  (B.Thouvenin ) : removal of key parsub_newdt, ws3max - local fractionnary step for particulate vertical transport -
   !&E       !  2019-08   (B.Thouvenin ) : adaptation for transfert
   !&E       !  2026-06  (M. BELHARET) : Modularization
   !&E---------------------------------------------------------------------
   !! * Modules used
   
   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast   !,kmax
   REAL(KIND=rsh),DIMENSION(ARRAY_WINDSPEED),INTENT(IN)       :: WIND_SPEED

   !! * Local declarations
   INTEGER                :: i,j,k,iv,kkmax,ivp
   REAL(KIND=rsh),DIMENSION(nv_state)  :: c,dc

   INTEGER :: ibsurf
   

#ifdef MUSTANG
   REAL(KIND=rsh)  :: tocdpe,depo
#endif


!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
    DO i=ifirst,ilast
      IF (htot(i,j) > RESIDUAL_THICKNESS_WAT)  THEN
            
       kkmax=NB_LAYER_WAT
       !!!!!!!!!!!!!!!!!!!!!!
       !!    loop on k     !!
       !!!!!!!!!!!!!!!!!!!!!!
       DO k=1,kkmax
       
         ! initialize to 0
         dc(1:nv_state)=0.0_rsh


         ! test pour identification des mailles de surface
         ! -----------------------------------------------
         ibsurf=0
         IF (k.eq.kkmax) ibsurf=1

         IF (kkmax > 1) THEN
           ! 3D to 1D for variable WATER_CONCENTRATION
          DO iv=1,nv_adv
              c(iv)=cvadv_wat_pos(iv,k,i,j)/unit_modif_mudbio_N2dw(irk_fil(iv))
          ENDDO
          DO iv=1,nv_fix
              c(nv_adv+iv)=cvfix_wat_pos(iv,k,i,j)/unit_modif_mudbio_N2dw(irk_fil(iv))
          ENDDO
         ENDIF
           
   ! test on salinity (if < p_sali_thhold_bio ==> no bio sink-sources)
         sali=max(0.0_rsh,SAL_BIOLink(k,i,j))

         IF(sali < p_sali_thhold_bio) THEN
         
          dc(:)=0.0_rsh
          
         ELSE

             CALL shared_parameters(i,j,k)

             ! DYNAMICS
             ! ............................................................

             CALL phytoplankton_dynamics(i,j,k,c,dc)
             CALL storage_diag_phytoplankton(i,j,k,c)
             CALL zooplankton_dynamics(i,j,k,c,dc)
             CALL chemistry_processes(i,j,k,c,dc)
             CALL oxygen_dynamics(i,j,k,c,dc,ibsurf)
#if ! defined key_BLOOM_opt2
             CALL estimate_sinking_velocity(i,j,k,c)
             CALL storage_diag_sinking_velocity(i,j,k)
#endif



          !**********************************************
          ! Verification of conservativity at one  point
          !**********************************************
#ifdef key_BIOLink_verif_conserv
          IF ((i .eq. i_BIOLink_verif) .and. (j .eq. j_BIOLink_verif) .AND. &
             mod(CURRENT_TIME,dt_conserv_BIO) <=TRANSPORT_TIME_STEP .and. (k == 1 .or. k == kkmax)) THEN
           !bilan dc et c pour un point de grille
              CALL verif_bilan_NP(k,i,j,dc,c)
          ENDIF
#endif
 
         ENDIF  !fin du test sur la salinite (si trop faible, pas de reaction bloomgique)

           ! 1D to 3D or variable BIO_SINKSOURCES(nv_state,NB_LAYER_WAT,limin:limax,ljmin:ljmax)

           !! cas ou on a fait un tableau expres avec bon ordre d indices
           !  BIO_SKSC_BIOLink (STATE_VAR_INDEXkij)=dc(1:nv_state)/86400.0_rsh  &
           !                         *unit_modif_mudbio_N2dw(irk_fil(1:nv_state))
           !! cas ou on garde le meme tableau qua dans modele hydro meme si pas bon ordre 
           ! evite de creer un nouveau tableau ce qui prend de la place
           BIO_SINKSOURCES(ADV_VAR_INDEXkij)=dc(1:nv_adv)/86400.0_rsh  &
                                 *unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))

           BIO_SKSC_FIX(FIXED_SKSC_INDEXkij)=dc(nv_adv+1:nv_adv+nv_fix)/86400.0_rsh  &
                                 *unit_modif_mudbio_N2dw(irk_fil(nv_adv+1:nv_adv+nv_fix))

       
       END DO  ! fin do k

       IF (kkmax == 1) THEN 
           DO k=2,NB_LAYER_WAT
             BIO_SINKSOURCES (ADV_VAR_INDEXkij)=0.0_rsh
             BIO_SKSC_FIX(FIXED_SKSC_INDEXkij)=0.0_rsh
           ENDDO
       ENDIF 

!!  rajout pour si on utilise MUSTANG qui doit connaitre la tendance au depot avec
!!  les nouvelles vitesses de chute mis a jour 
#if defined MUSTANG && ! defined key_Pconstitonly_insed
       DO ivp=nvpc+1,nvp
            tocdpe=tocd(ivp)+0.00001_rsh
            flx_w2s(ivp,i,j)=ws3(1,ivp,i,j)                   &
                          *MAX(0.0_rsh,1.0_rsh-(tauskin(i,j)/tocdpe)) &
                          *tocd(ivp)/tocdpe
            depo=flx_w2s(ivp,i,j)*cw_bottom_MUSTANG(ivp,i,j)
            IF(irkm_var_assoc(ivp) > 0 ) THEN
              IF (flx_w2s(irkm_var_assoc(ivp),i,j)==0.0_rsh) flx_w2s(ivp,i,j)=0.0_rsh
!            ELSE
!              IF(depo .LE. epsdep_MUSTANG)flx_w2s(ivp,i,j)=0.0_rsh
            END IF
       ENDDO
#endif

     END IF  ! H>RESIDUAL_THICKNESS_WAT
   END DO
   END DO
!$OMP END DO

  END SUBROUTINE bloom_sksc_wat

   !!======================================================================
#ifdef key_BIOLink_verif_conserv
     SUBROUTINE verif_bilan_NP(k,i,j,dc,c)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE verif_bilan_NP  ***
   !&E
   !&E ** Purpose : VERIFICATION conservativite bilan de matiere en N et P
   !&E              en mmolN/m3 en une maille k,i,j
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !         Benedicte issu de la version r1023 de Alain (MANGA4)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   
   !! * Arguments
   INTEGER, INTENT(in)                       :: k,i,j      ! ocean spatial-step index
   REAL(KIND=rsh),DIMENSION(:),INTENT(in)    :: dc,c

   ! test de conservativite de N et P
   ! --------------------------------
   REAL(KIND=rsh) :: somdcN,somdcP

   
!     Conservativite en azote     
!     -----------------------------------------------------------------------------
#ifdef key_BLOOM_opt2
    somdcN=dc(iv_nutr_NH4)+dc(iv_nutr_NO3)+dc(iv_phyto_diat_N)+dc(iv_phyto_dino_N)+dc(iv_phyto_nano_N)+ &
           dc(iv_detr_N)+dc(iv_diss_N)+dc(iv_diss_fond_Nitr)+rNC*(dc(iv_zoo_meso_N)+dc(iv_zoo_micr_N))
#else
    somdcN=dc(iv_nutr_NH4)+dc(iv_nutr_NO3)+dc(iv_phyto_diat_N)+dc(iv_phyto_dino_N)+dc(iv_phyto_nano_N)+ &
           dc(iv_detr_N)+dc(iv_zoo_meso_N)+dc(iv_zoo_micr_N) 

#endif

        IF(abs(somdcN) > 1.0e-3) THEN
          write(*,*)'Probleme conservativite N - Voir fichier conserv_BIOLink...'
          write(iscreenlog_conserv,*)'____________________________________________________'
          write(iscreenlog_conserv,*)'!!!! PB conservativite N: somdcN =',somdcN
          write(iscreenlog_conserv,*)'!!!! maille i, j, k:',i,j,k
          write(iscreenlog_conserv,*)'dc: NH4,NO3,diat,dino,nano,detr ',dc(iv_nutr_NH4),dc(iv_nutr_NO3),dc(iv_phyto_diat_N), &
                                                       dc(iv_phyto_dino_N),dc(iv_phyto_nano_N),dc(iv_detr_N)
          write(iscreenlog_conserv,*)'dc: mesozoo,microzoo',dc(iv_zoo_meso_N),dc(iv_zoo_micr_N)
!          stop
        ENDIF


!     Conservativite en phosphore  
!     ---------------------------------------------------------------------------------
#ifdef key_BLOOM_opt2
        somdcP=dc(iv_nutr_PO4)+dc(iv_nutr_Pads)+dc(iv_detr_P)+dc(iv_diss_P)+ &
              (dc(iv_phyto_diat_N)+dc(iv_phyto_dino_N)+dc(iv_phyto_nano_N)+ &
              +rNC*(dc(iv_zoo_meso_N)+dc(iv_zoo_micr_N)))*rappaz 
#else
        somdcP=dc(iv_nutr_PO4)+dc(iv_nutr_Pads)+dc(iv_detr_P)+ &
              (dc(iv_phyto_diat_N)+dc(iv_phyto_dino_N)+dc(iv_phyto_nano_N)+ dc(iv_zoo_meso_N)+dc(iv_zoo_micr_N))*rappaz 

#endif

        IF(abs(somdcP) > 1.0e-3) THEN
          write(*,*)'Probleme conservativite P - Voir fichier conserv_BIOLink...'
          write(iscreenlog_conserv,*)'____________________________________________________'
          write(iscreenlog_conserv,*)'!!!! PB conservativite P: somdcP =',somdcP
          write(iscreenlog_conserv,*)'!!!! maille i, j, k:',i,j,k
          write(iscreenlog_conserv,*)'dc: PO4,Pads,diat,dino,nano,detr=   ',dc(iv_nutr_PO4),dc(iv_nutr_Pads), &
                        dc(iv_phyto_diat_N)*rappaz,dc(iv_phyto_dino_N)*rappaz,dc(iv_phyto_nano_N)*rappaz,dc(iv_detr_P)
          write(iscreenlog_conserv,*)'dc: mesozoo,microzoo=  ',dc(iv_zoo_meso_N)*rappaz,dc(iv_zoo_micr_N)*rappaz
        ENDIF
     
   END SUBROUTINE verif_bilan_NP
#endif
 
   !!======================================================================
   
  SUBROUTINE bloom_eval_diag2d(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_eval_diag2d  ***
   !&E---------------------------------------------------------------------

   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast  !,kmax

   !! * Local declarations
   INTEGER                  ::  i,j,k         ! loop indexes


   REAL(KIND=rsh)   :: dzw

                            ! diagnoses
                            ! ---------
   INTEGER                  :: kmaxdiat,kmaxdino,kmaxnano,kminO2,kdcdz
   REAL(KIND=rsh)           :: xmaxdiat,datemaxdiat,xmaxdino,datemaxdino, &
                               xmaxnano,datemaxnano,xminO2,dateminO2,     &
                               cumuldiat,cumuldino,cumulnano,&
                               xmaxchloro,kmaxchloro,datemaxchloro,d,voleau
   REAL(KIND=rlg)           :: dcdzmax,dcdz

                            ! variables temperature, salinite
                            ! -------------------------------
   REAL(KIND=rsh)           :: temper,tempabs,sali
   REAL(KIND=rlg)           :: row,rost
   REAL(KIND=rlg),DIMENSION(NB_LAYER_WAT)  :: density



!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
   

    DO i=ifirst,ilast

#if ! defined key_BLOOM_opt2
       !variable diagnostique des max et des dat de max remises a zero tous les mois dans l option 1
       IF ((imois_BIOLINK.eq.1).and.(ijour_BIOLINK.eq.1).and.(iheure_BIOLINK.eq.0)) THEN
            diag_2d(irk_diag(id_diat_max),i,j)=0.0_rsh
            diag_2d(irk_diag(id_dino_max),i,j)=0.0_rsh
            diag_2d(irk_diag(id_nano_max),i,j)=0.0_rsh

            diag_2d(irk_diag(id_diat_datemax),i,j)=0.0_rsh
            diag_2d(irk_diag(id_dino_datemax),i,j)=0.0_rsh
            diag_2d(irk_diag(id_nano_datemax),i,j)=0.0_rsh

            diag_2d(irk_diag(id_diat_depthmax),i,j)=0.0_rsh
            diag_2d(irk_diag(id_dino_depthmax),i,j)=0.0_rsh
            diag_2d(irk_diag(id_nano_depthmax),i,j)=0.0_rsh
       ENDIF
#endif

      IF (htot(i,j) > RESIDUAL_THICKNESS_WAT)  THEN

 
   ! diatomees
   ! ---------
         kmaxdiat=1
         xmaxdiat=0.0_rsh
         datemaxdiat=1
         DO k=1,NB_LAYER_WAT
           IF (cvadv_wat_pos(iv_phyto_diat_N,k,i,j) .GE. xmaxdiat) THEN
             xmaxdiat=cvadv_wat_pos(iv_phyto_diat_N,k,i,j)
             kmaxdiat=k
             datemaxdiat=jjulien_BIOLINK
           END IF
         END DO

         IF (xmaxdiat.GT.diag_2d(irk_diag(id_diat_max),i,j)) THEN
           diag_2d(irk_diag(id_diat_max),i,j)=xmaxdiat
           diag_2d(irk_diag(id_diat_datemax),i,j)=datemaxdiat
           diag_2d(irk_diag(id_diat_depthmax),i,j)=-(thicklayerW_W(kmaxdiat,i,j)-WATER_ELEV_ij)
         END IF

   ! dinoflagelles
   ! -------------
         xmaxdino=0.0_rsh
         kmaxdino=1
         datemaxdino=1
         DO k=1,NB_LAYER_WAT
           IF (cvadv_wat_pos(iv_phyto_dino_N,k,i,j).GE.xmaxdino) THEN
             xmaxdino=cvadv_wat_pos(iv_phyto_dino_N,k,i,j)
             kmaxdino=k
             datemaxdino=jjulien_BIOLINK
           END IF
         END DO

         IF (xmaxdino.GT.diag_2d(irk_diag(id_dino_max),i,j)) THEN
           diag_2d(irk_diag(id_dino_max),i,j)=xmaxdino
           diag_2d(irk_diag(id_dino_datemax),i,j)=datemaxdino
           diag_2d(irk_diag(id_dino_depthmax),i,j)=-(thicklayerW_W(kmaxdino,i,j)-WATER_ELEV_ij)
         END IF

   ! nanoflagelles
   ! -------------
         xmaxnano=0.0_rsh
         kmaxnano=1
         datemaxnano=1
         DO k=1,NB_LAYER_WAT
           IF (cvadv_wat_pos(iv_phyto_nano_N,k,i,j).GE.xmaxnano) THEN
             xmaxnano=cvadv_wat_pos(iv_phyto_nano_N,k,i,j)
             kmaxnano=k
             datemaxnano=jjulien_BIOLINK
           END IF
         END DO

         IF (xmaxnano.GT.diag_2d(irk_diag(id_nano_max),i,j)) THEN
           diag_2d(irk_diag(id_nano_max),i,j)=xmaxnano
           diag_2d(irk_diag(id_nano_datemax),i,j)=datemaxnano
           diag_2d(irk_diag(id_nano_depthmax),i,j)=- (thicklayerW_W(kmaxnano,i,j)-WATER_ELEV_ij)
         END IF


   ! calcul de l intensite et profondeur de l halocline
   ! --------------------------------------------------
         dcdzmax=0.0_rlg
         kdcdz=NB_LAYER_WAT
         DO k=1,NB_LAYER_WAT-1
            dcdz=(SAL_BIOLink(k,i,j)-SAL_BIOLink(k+1,i,j))/thicklayerW_W(k,i,j)
            dcdz=ABS(dcdz)
           IF (dcdz > dcdzmax) THEN
             kdcdz=k+1
             dcdzmax=dcdz
           END IF
         END DO
         diag_2d(irk_diag(id_gradsali_max),i,j)=dcdzmax
         IF (dcdzmax.lt.0.05) then
          diag_2d(irk_diag(id_gradsali_depthmax),i,j)= BATHY_H0(i,j)
         ELSE 
          diag_2d(irk_diag(id_gradsali_depthmax),i,j)=-(thicklayerW_W(kdcdz,i,j)-WATER_ELEV_ij)
         ENDIF

         IF(id_gradtemp_max .NE. 0 ) THEN
   ! calcul de l intensite et profondeur de la thermocline
   ! -----------------------------------------------------
           dcdzmax=0.0_rlg
           kdcdz=0
           DO k=1,NB_LAYER_WAT-1
            dcdz=ABS(TEMP_BIOLink(k,i,j)-TEMP_BIOLink(k+1,i,j))/thicklayerW_W(k,i,j)
            IF (dcdz > dcdzmax) THEN
               kdcdz=k+1
               dcdzmax=dcdz
            END IF
           END DO
           diag_2d(irk_diag(id_gradtemp_max),i,j)=dcdzmax
           IF (dcdzmax.lt.0.025) then
            diag_2d(irk_diag(id_gradtemp_depthmax),i,j)= BATHY_H0(i,j)
           ELSE
            diag_2d(irk_diag(id_gradtemp_depthmax),i,j)=-(thicklayerW_W(kdcdz,i,j)-WATER_ELEV_ij)
           ENDIF

         ELSE 
   ! calcul de l intensite et profondeur de la pycnocline
   ! ----------------------------------------------------
         
           ! Calcul de la densite 
           DO k=1,NB_LAYER_WAT
             temper=TEMP_BIOLink(k,i,j)
             sali=SAL_BIOLink(k,i,j)
             ! density of pure water at temperature temper
             row = 999.842594_rlg + (  6.793952d-02 + ( - 9.095290d-03 +                &
                   (1.001685d-04 - 1.120083d-06*temper + 6.536332d-09*temper*temper)*temper )  &
                                        *temper  )*temper
             ! density of sea water at salinity sali and temperature temper
             rost = row + sali * (   0.824493_rlg + (  -4.0899d-03 + ( 7.6438d-05 +                   &
                (-8.2467d-07+5.3875d-09*temper)*temper )*temper)*temper) +             &
                sali**1.5_rlg * ( -5.72466d-03 + (1.0227d-04-1.6546d-06*temper)*temper) +  &
                sali*sali*4.8314d-04
             density(k)=rost
           ENDDO
     
           dcdzmax=0.0_rlg
           kdcdz=NB_LAYER_WAT
           DO k=1,NB_LAYER_WAT-1
              dcdz=ABS(density(k)-density(k+1))/thicklayerW_W(k,i,j)
              IF (dcdz > dcdzmax) THEN
                kdcdz=k+1
                dcdzmax=dcdz
              END IF
           END DO

           diag_2d(irk_diag(id_graddens_max),i,j)=dcdzmax
         ! diag_2d(irk_diag(id_graddens_max),i,j)=int(10000.0_rsh*diag_2d(irk_diag(id_graddens_max),i,j))+dcdzmax/10000.0_rsh
           IF (dcdzmax.lt.0.02) then
              diag_2d(irk_diag(id_graddens_depthmax),i,j)= BATHY_H0(i,j)
              ! diag_2d(irk_diag(id_graddens_depthmax),i,j)=int(10000.0_rsh*diag_2d(irk_diag(id_graddens_depthmax),i,j))+BATHY_H0(i,j)/10000.0_rsh
           ELSE 
              diag_2d(irk_diag(id_graddens_depthmax),i,j)=-(thicklayerW_W(kdcdz,i,j)-WATER_ELEV_ij)
              ! diag_2d(irk_diag(id_gradtemp_depthmax),i,j)=int(10000.0_rsh*diag_2d(irk_diag(id_gradtemp_depthmax),i,j))-(f_zw(BATHY_H0(i,j),WATER_ELEV_ij,kdcdz,i,j)-WATER_ELEV_ij)/10000.0_rsh
           ENDIF

         ENDIF
         
   ! niveau du gradient vertical de temperature max pour injection
   ! des concentrations surface-fond en limite ouest et sud
   !   
   ! calcul de la production cumulee integree sur la verticale
   ! ATTENTION : production exprimee en carbone pour l option 1 et en azote pour l option 2 
   ! ---------------------------------------------------------------------------------------
#if ! defined key_BLOOM_opt2
        !production cumulee depuis le 1er janvier dans l option 1
         IF ((imois_BIOLINK.eq.1).and.(ijour_BIOLINK.eq.1).and.(iheure_BIOLINK.eq.0)) THEN
          DO k=1,NB_LAYER_WAT
            !FIXED_VAR_CONC(ivfix_cumulprod_first:ivfix_cumulprod_last,k,i,j)=0.0_rsh
            FIXED_VAR_CONC(FIXED_VAR_INDEX_CUMULPROD)=0.0_rsh
            cvfix_wat_pos(1:nv_fix,k,i,j)=0.0_rsh

          END DO
         ENDIF
#endif
         cumuldiat=0.0_rsh
         cumuldino=0.0_rsh
         cumulnano=0.0_rsh

         DO k=1,NB_LAYER_WAT
           cumuldiat=cumuldiat+cvfix_wat_pos(iv_phyto_diat_pp-nv_adv,k,i,j)
           cumuldino=cumuldino+cvfix_wat_pos(iv_phyto_dino_pp-nv_adv,k,i,j)
           cumulnano=cumulnano+cvfix_wat_pos(iv_phyto_nano_pp-nv_adv,k,i,j)
         END DO
         diag_2d(irk_diag(id_diat_columnprod),i,j)=cumuldiat
         diag_2d(irk_diag(id_dino_columnprod),i,j)=cumuldino
         diag_2d(irk_diag(id_nano_columnprod),i,j)=cumulnano

#if ! defined key_BLOOM_opt2
   !   production cumulee dans le temps et sur la verticale de tous les types phytos   
         diag_2d(irk_diag(id_columnprodtotal),i,j)=cumuldiat+cumuldino+cumulnano
#endif                                             
           
     END IF  ! H>RESIDUAL_THICKNESS_WAT
   END DO
   END DO
!$OMP END DO
      
  END SUBROUTINE bloom_eval_diag2d
!!============================================================================================
#if defined key_BLOOM_insed  && ! defined key_BLOOM_opt2
  SUBROUTINE bloom_reactions_in_sed(ifirst,ilast,jfirst,jlast,dt_true)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_reactions_in_sed  ***
   !&E
   !&E ** Purpose : Estimate transformations inside sediments
   !&E              Estimate new concentrations in sediments (cv_sed)
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !          (A. Menesguen) Original code
   !&E       !  2006-11 (V. Garnier) module - mise en forme
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   
!
!!========

   !! * Arguments
   INTEGER, INTENT(IN)            :: ifirst,ilast,jfirst,jlast                           
   REAL(KIND=rlg),INTENT(IN)      :: dt_true

   !! * Local declarations
   INTEGER                :: i,j,k,iv,ksmin,ksmax
   REAL(KIND=rsh)         :: dtbiojour,  zmiddle
   REAL(KIND=rsh),DIMENSION(nv_state)  :: cw_bott 

   !!----------------------------------------------------------------------
   !! * Executable part
   
   dtbiojour=dt_true/86400._rsh

!$OMP DO SCHEDULE(RUNTIME)

     DO j=jfirst,jlast
      DO i=ifirst,ilast

       IF (BATHY_H0(i,j) > -valmanq .AND. ksma(i,j) > 0)  THEN  ! on n est pas a terre et il y a du sediment
       
        ksmin=ksmi(i,j)
        ksmax=ksma(i,j)
        
     ! memorisation des concentrations locales dans l eau surnageante en gardant valeurs seulement positives
        IF (htot(i,j)> RESIDUAL_THICKNESS_WAT) THEN  ! test sur hauteur d eau
            cw_bott(1:nv_state)=0.5_rsh*( cw_bottom_MUSTANG(1:nv_state,i,j)+   &
                                ABS(cw_bottom_MUSTANG(1:nv_state,i,j)) ) 
                        
        ELSE
          !   il n y a pas d eau
          cw_bott(1:nv_state)=0.0_rsh
          
        ENDIF
         
        ! conversion des concentrations dans l eau si changement d unite 
        ! si matiere organique supposee constitutive du sediment
        cw_bott(1:nv_state)= cw_bott(1:nv_state)/unit_modif_mudbio_N2dw(irk_fil(1:nv_state))                            
         

        ! calcul des flux reactifs pour chaque variable et chaque processus dans chaque couche
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zmiddle=0.0_rsh
        DO k=ksmax,ksmin,-1
          ! no bio transformation if very small layer
          IF (dzs(k,i,j) > dzsmin) THEN

              CALL sediment_bio_dynamics(i,j,k,dtbiojour, zmiddle, cw_bott)

           ENDIF ! test dzs 
         END DO  ! boucle ks
       END IF   ! test a terre
     END DO
   END DO

!$OMP END DO

  END SUBROUTINE bloom_reactions_in_sed
#endif


#endif
 END MODULE
