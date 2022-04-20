  MODULE bloom

#include "cppdefs.h"

#if defined SUBSTANCE && defined BLOOM

  !!======================================================================
   !!                   ***  MODULE  BLOOM ***
   !!
   !! Biologic dynamics:  - all routines which are independant of host hydro model
   !!======================================================================

#ifdef key_MARS
#include "toolcpp.h"
#include "coupleur_dimhydro_BIOLink.h"
#if defined MUSTANG &&  defined key_BLOOM_insed
#include "../SEDIM/coupleur_define_MUSTANG.h"
#endif
#else
!!! CROCO
!#ifdef MUSTANG
!   USE module_MUSTANG
!#endif
   USE module_BIOLink
#endif

#include "coupleur_define_BIOLink.h"

   !! * Modules used
   USE comBIOLink
   USE comBIOLink_physics
   USE comBIOLink_helping

#if defined MUSTANG
   USE comMUSTANG
#endif
#if ! defined key_MARS
   USE comsubstance
#endif

 
   IMPLICIT NONE

   !! * Accessibility
   PUBLIC bloom_sksc_wat,bloom_eval_diag2d,bloom_SPMtot_Chla,bloom_extinction_avg          
#if defined MUSTANG &&  defined key_BLOOM_insed
   PUBLIC bloom_reactions_in_sed
#endif
#if defined key_MANGAbio && defined key_MANGAbiovague
   PUBLIC bloom_wavefile_MANGAbio
#endif
   !! * Shared module variables

 !  
   !! * Private variables

 CONTAINS

   !!======================================================================
  SUBROUTINE bloom_SPMtot_Chla(ifirst,ilast,jfirst,jlast          &
#if defined key_MANGAbio && defined key_MANGAbiovague
                         ,forcSPMk                                &
#endif 
#if defined key_messat
                         ,forcSPM                                 &
#endif
                         )

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
#if defined key_messat
   REAL(KIND=rsh), DIMENSION(PROC_IN_ARRAY),INTENT(IN)        :: forcSPM 
#endif  
#if defined key_MANGAbio && defined key_MANGAbiovague
   REAL(KIND=rsh), DIMENSION(NB_LAYER_WAT,PROC_IN_ARRAY),INTENT(IN)      :: forcSPMk
#endif 

   
   !! * Local declarations
    INTEGER                  ::  i,j,k,kmaxmod,iv        ! loop indexes
    REAL(KIND=rsh)                                  :: cpom,fact_phyto_ChlNratio
   !!----------------------------------------------------------------------
   !! * Executable part

!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
#ifdef key_MARS
     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
     DO i=ifirst,ilast
#endif
       IF (htot(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN
            
#ifdef key_MARS
           IF(htot(i,j) < hm ) THEN
               kmaxmod=1
           ELSE
               kmaxmod=NB_LAYER_WAT
           ENDIF

#else
               ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
               kmaxmod=NB_LAYER_WAT
#endif

           DO k=1,kmaxmod
             cmes_3dmgl(k,i,j)=0.0_rsh
             DO iv=1,nvpc
               cmes_3dmgl(k,i,j)=cmes_3dmgl(k,i,j)+cvadv_wat_pos(iv,k,i,j)*1000.0_rsh  !cmes_3dmgl passe en mg/l
             END DO
             diag_3d_wat(irk_diag(id_spm_total),k,i,j)=cmes_3dmgl(k,i,j)

             cpom=0.0_rsh

#if defined key_messat
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !! forcing SPMtot with satellite data
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined key_BLOOM_opt2
             IF (forcSPM(i,j)*1000.0_rsh .GT. cmes_3dmgl(k,i,j)) THEN  !forcage en g/l
                diag_3d_wat(irk_diag(id_spim_satused),k,i,j)=1.0_rsh
             ELSE
                diag_3d_wat(irk_diag(id_spim_satused),k,i,j)=0.0_rsh
             END IF
#endif
             cpom=forcSPM(i,j)*1000.0_rsh ! mg/L
#endif


#if defined key_turbclim && defined key_daily_climato_kpar
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !! forcing SPMtot with climato satellite data
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             cpom=mes_sat(i,j)
#endif

             !       Cas ou on ne tient pas compte des MES simulees  ==>
             !cmes_3dmgl(k,i,j)=cpom               !M.E.S en mg/l
             !                 OU
             !       Cas ou on tient compte des MES simulees  ==>
             cmes_3dmgl(k,i,j)=max(cmes_3dmgl(k,i,j),cpom)     !M.E.S en mg/l


#if defined key_MANGAbio && defined key_MANGAbiovague
             cmes_3dmgl(k,i,j)=0.5_rsh*(forcSPMk(k,i,j)+ABS(forcSPMk(k,i,j)))     !Eq en mg/l
             diag_3d_wat(irk_diag(id_spm_total),k,i,j)=cmes_3dmgl(k,i,j)
#endif


             ! chlorophyl concentration in mug.l-1
             ! rate Chloro/N dependent on light extinction
             IF(l_ChlNratio_var) THEN
               IF((CURRENT_TIME-TIME_BEGIN)>345600.0_rlg) THEN  !debut calcul apres 4 jours initiaux
                  fact_phyto_ChlNratio=max(0.70_rsh,p_phyto_ChlNratiomax*(extinction_ave4d(k,i,j)/p_phyto_ChlN_ksmithextinct)/ &
                              (sqrt(1+(extinction_ave4d(k,i,j)/p_phyto_ChlN_ksmithextinct)**2)))
               ELSE
                  fact_phyto_ChlNratio=max(0.70_rsh,p_phyto_ChlNratiomax*(EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct)/ &
                              (sqrt(1+(EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct)**2)))
               ENDIF
             ELSE
                  fact_phyto_ChlNratio=p_phyto_ChlNratio
             ENDIF
             BIOLink_chloro(k,i,j)=(cvadv_wat_pos(iv_phyto_diat_N,k,i,j)+ cvadv_wat_pos(iv_phyto_dino_N,k,i,j)  &
                          +cvadv_wat_pos(iv_phyto_nano_N,k,i,j) &
#ifdef key_psnz
                          +cvadv_wat_pos(iv_phyto_psnz_N,k,i,j)   &
#endif
#ifdef key_karenia
                          +cvadv_wat_pos(iv_phyto_karenia_N,k,i,j)   &
#endif
#ifdef key_phaeocystis
                          +cvadv_wat_pos(iv_phyto_phaeocystis_cell_N,k,i,j)+ cvadv_wat_pos(iv_phyto_phaeocystis_colo_N,k,i,j)  &
#endif
                          +epsilon_BIOLink)*fact_phyto_ChlNratio

              !BIOLink_chloro(k,i,j)=(ABS(BIOLink_chloro(k,i,j))+BIOLink_chloro(k,i,j))*0.5_rsh

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
#ifdef key_MARS
     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
     DO i=ifirst,ilast
#endif
       IF (htot(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN
            
#ifdef key_MARS
           IF(htot(i,j) < hm ) THEN
               kmaxmod=1
           ELSE
               kmaxmod=NB_LAYER_WAT
           ENDIF

#else
               ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
               kmaxmod=NB_LAYER_WAT
#endif

         !  memorization of extinction_aveh and extinction_ave4d, extinction_tab
          diag_3d_wat(irk_diag(id_extinctioncoeff),:,i,j)=EXTINCTION_RAD(:,i,j)

         IF(l_ChlNratio_var) THEN
           ! - calcul de l extinction moyenne sur 4 jours  --------------------------------------------
           DO k=LOOPK_SUBSURF_TO_BOTTOM_WAT   ! kmaxmod-1,1,-1
             extinction_aveh(k,i,j)=extinction_aveh(k,i,j)+EXTINCTION_RAD(k,i,j)*dtbio
             !if(k==10 .and.i==20 .and.j==20)write(*,*)'extinction_aveh',extinction_aveh(k,i,j),EXTINCTION_RAD(k,i,j)
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
                extinction_tab(numday_extinction,numhour_extinction,k,i,j)=extinction_tab(numday_extinction,numhour_extinction,1,i,j)
                extinction_ave4d(k,i,j)=extinction_ave4d(1,i,j)
             ENDIF
           ENDDO
         ENDIF
       ENDIF ! if water
     ENDDO       
   ENDDO 
!$OMP END DO
      
  END SUBROUTINE bloom_extinction_avg


   !!======================================================================
   
  SUBROUTINE bloom_sksc_wat(ifirst,ilast,jfirst,jlast                             &
#if defined BLOOM
#if defined GLS_MIXING
                          ,CIN_TURBULENT_ENERGY                                    &
#endif
#if defined key_zostera || defined key_oxygen  
                          ,WIND_SPEED                                              &       
#endif
#if (defined key_oyster_SFG || defined key_oyster_DEB) &&  defined key_MARS
                         ,CELL_SURF                                                &
#endif
#endif
                          )
                         !TIME_BEGIN )
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
   !&E      !  2019-08   (B.Thouvenin ) : adaptation for transfert
   !&E---------------------------------------------------------------------
   !! * Modules used
   
   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast   !,kmax
   !REAL(KIND=rlg),INTENT(IN)                                  :: TIME_BEGIN
#if defined BLOOM
#if defined GLS_MIXING
   REAL(KIND=rsh),DIMENSION(ARRAY_CIN_TURB_ENERGY),INTENT(IN) :: CIN_TURBULENT_ENERGY
#endif
#if defined key_zostera || defined key_oxygen  
   REAL(KIND=rsh),DIMENSION(ARRAY_WINDSPEED),INTENT(IN)       :: WIND_SPEED
#endif
#if (defined key_oyster_SFG || defined key_oyster_DEB) &&  defined key_MARS
   REAL(KIND=rsh),DIMENSION(ARRAY_CELL_SURF),INTENT(IN)          :: CELL_SURF
#endif   
#endif


   !! * Local declarations
   INTEGER                :: i,j,k,iv,kkmax,ivp
   REAL(KIND=rsh),DIMENSION(nv_state)  :: c,dc
   REAL(KIND=rsh)      :: dtbiojour,foncsinus,txfiltbenth

   INTEGER :: ibbenth,ibsurf
   
   ! variables temperature, salinite, epaisseur couche
   ! -------------------------------
   REAL(KIND=rsh) :: temper,tempabs,sali,epn
   REAL(KIND=rsh) :: effetchaleur,effetturbidite
      
   ! variables lumiere
   ! ------------------------
   REAL(KIND=rsh) :: fluxrelatifsurf,fluxrelatiffond
   REAL(KIND=rsh) :: effetlumierediat,effetlumieredino,effetlumierenano
   
   ! variables phytoplancton
   ! -----------------------
   REAL(KIND=rsh) :: fact_phyto_ChlNratio
   REAL(KIND=rsh) :: fractionno3diat, fractionno3dino,fractionno3nano
   REAL(KIND=rsh) :: fractionnh4diat, fractionnh4dino,fractionnh4nano


   ! variables diatomees
   ! -------------------
   REAL(KIND=rsh) :: rationdiat,diatmorteau,effetselnutdiat,excretdiat
   REAL(KIND=rsh) :: effetnitdiat,effetamdiat,diat_iksmith
   REAL(KIND=rsh) :: effetazotediat,effetsilice,effetphosphorediat
   REAL(KIND=rsh) :: pslimdiat,Si

   ! variables dinoflagelles
   ! -----------------------
   REAL(KIND=rsh) :: rationdino,dinomorteau,excretdino
   REAL(KIND=rsh) :: effetnitdino,effetamdino,effetselnutdino
   REAL(KIND=rsh) :: pslimdino,effetazotedino,dino_iksmith
   REAL(KIND=rsh) :: effetphosphoredino

   ! variables nano_pico_phytoplancton
   ! ---------------------------------
   REAL(KIND=rsh) :: rationnano,nanomorteau,excretnano
   REAL(KIND=rsh) :: effetnitnano,effetamnano,effetazotenano,effetselnutnano
   REAL(KIND=rsh) :: effetphosphorenano,pslimnano,nano_iksmith

   ! phtoplanctons toxiques
   ! ---------------------------------
#ifdef key_karenia
   REAL(KIND=rsh) :: rationkarenia,kareniamorteau,fractionno3karenia,fractionnh4karenia
   REAL(KIND=rsh) :: effetlumierekarenia,effetselnutkarenia,uptakeN,uptakeP
#endif
#ifdef key_psnz
   REAL(KIND=rsh) :: rationpsnz,psnzmorteau,fractionno3psnz,fractionnh4psnz,effetlumierepsnz,effetselnutpsnz
#endif
#ifdef key_phaeocystis
   REAL(KIND=rsh) :: rationphaeocystiscell,fractionno3phaeocystiscell,fractionnh4phaeocystiscell,phaeocystismortcell
   REAL(KIND=rsh) :: rationphaeocystiscolo,fractionno3phaeocystiscolo,fractionnh4phaeocystiscolo,phaeocystismortcolo
   REAL(KIND=rsh) :: initcolonie,effetlumierephaeocystis,phaeocystislysecolo,effetselnutphaeocystiscolo   
   REAL(KIND=rsh) :: effetselnutphaeocystiscell
#endif

   ! variables mesozooplancton
   ! -------------------------   
   REAL(KIND=rsh) :: totproiebrut,totproie,fmesozoo
   REAL(KIND=rsh) :: captdiatmeso,captdino,captmicrozoo,captmicrozooN
   REAL(KIND=rsh) :: meso_kivlev, rationmesozoo, rNC
   REAL(KIND=rsh) :: broumesozoodiat,broumesozoodino,broumesozoomicrozoo
   REAL(KIND=rsh) :: txmortmesozoo,assimilmesozoo,excretionmesozoo 
     
#ifdef key_phaeocystis
   REAL(KIND=rsh) :: broumesozoophaeocolo
#endif
#ifdef key_karenia
   REAL(KIND=rsh) :: broumesozookarenia
#endif
#ifdef key_psnz
   REAL(KIND=rsh) :: broumesozoopsnz
#endif

   ! variables microzooplancton
   ! --------------------------
   REAL(KIND=rsh) :: fmicrozoo,captnano,captdet,captdinomicro
   REAL(KIND=rsh) :: captdiatmicro,broumicrozoodino 
    REAL(KIND=rsh) :: rationmicrozoo,excretionmicrozoo
   REAL(KIND=rsh) :: txmortmicrozoo,broumicrozoonano
   REAL(KIND=rsh) :: broumicrozoodet,broumicrozoodiat,assimilmicrozoo
#ifdef key_phaeocystis
   REAL(KIND=rsh) :: broumicrozoophaeocell
   REAL(KIND=rsh) :: captphaeocell,captphaeocolo
#endif
#ifdef key_psnz
   REAL(KIND=rsh) :: broumicrozoopsnz,captpsnz
#endif       
#ifdef key_karenia
   REAL(KIND=rsh) :: broumicrozookarenia,captkarenia
#endif

 
#ifdef key_BLOOM_opt2
   REAL(KIND=rsh) :: dissolMOPeau  
#endif

#ifdef key_ulvas
   ! variables ulves
   ! ---------------
   REAL(KIND=rsh) :: ulvemorteau,ulvebenthmortsed,pompageazoteulve,pompageazoteulvebenth   
   REAL(KIND=rsh) :: rationulve,rationulvebenth,partulvessurface,vitessechuteulve,ulveresuspazote
   REAL(KIND=rsh) :: fractionno3,fractionnh4
#endif

   ! parametrisation broutage par benthos
   ! ------------------------------------
   REAL(kind=rsh) :: detresuspazote,detresusppho,vitessechutedet,vitessechutediat
   
   ! variables mineralisation
   ! ------------------------
   REAL(KIND=rsh) :: xnitrifeau,reminazdeteau,reminpdeteau,dissolsiliceeau 

   ! parametrisation broutage par benthos
   ! ------------------------------------
   REAL(KIND=rsh)           :: txfiltbenthij

   ! processus d adsorption/desorption
   ! ---------------------------------
   REAL(KIND=rsh) :: desorpeau,adsorpeau
   REAL(KIND=rsh) :: varads,cmaxdesorp,saturmesp,cmaxadsorp
   
#if ! defined key_BLOOM_opt2
   REAL(KIND=rsh)           :: p_mesz_thrN_turb
#endif
#if defined key_BLOOM_opt2
   ! sedimentation mat detr azotee derniere couche
   REAL(KIND=rsh) :: fractionDET_N
   REAL(KIND=rsh) :: fractionDET_P
   REAL(KIND=rsh) :: fractionDET_Si
   REAL(KIND=rsh) :: fractionDIAT_N
   REAL(KIND=rsh) :: fractionDIAT_P
   REAL(KIND=rsh) :: fractionDIAT_Si
#endif
#ifdef key_oyster_benthos
!   filtration benthos par huitres
!   -------------------------------
   REAL(KIND=rsh)           :: cmes1,colmat_huitre,effettemp_huitre,filtmes,filt
#endif   

#if defined key_zostera
   ! "mailles" de sediment 
   ! ---------------------
   REAL(KIND=rsh),PARAMETER :: epnsed=0.4_rsh ! epaisseur (m) de la "maille de sediment" (bricole)
   REAL(KIND=rsh),PARAMETER :: deltaz = 0.3_rsh    ! distance (m) fixe entre la maille de fond et le "sediment" (bricole)
   REAL(KIND=rsh) :: Gross_Prod_zost,Leaf_Resp_zost,Root_Resp_zost
#endif

   ! --------------------------------------------
   ! variables used in include associated modules
   ! --------------------------------------------

#ifdef key_psnz
!  variables used in include incellwat_bloom_pseudocitzschia.h
!  ----------------------------------------------   
   REAL(KIND=rsh) :: effetchutepsnz
   REAL(KIND=rsh) :: effetchaleurpsnz,effetsalinitepsnz,effetno3psnz,effetnh4psnz
   REAL(KIND=rsh) :: effetazotepsnz,effetsilicepsnz,effetphosphorepsnz,pslimpsnz
   REAL(KIND=rsh) :: quotaSi,uptakeSi,prodacidedomo,effetlumierepsnz_24h
   REAL(KIND=rsh) :: psnz_iksmith,excretpsnz,quotamaxsurminSi
#endif   

#ifdef key_karenia
!  variables used in include incellwat_bloom_karenia.h
!  ----------------------------------------------   
   REAL(KIND=rsh) :: deltadensite,salisup,deltahaut,deltabas
   REAL(KIND=rlg) :: row,rost
   REAL(KIND=rsh) :: effetchaleurkarenia,effetno3karenia,effetnh4karenia,quotaN,quotaP
   REAL(KIND=rsh) :: effetazotekarenia,effetphosphorekarenia,pslimkarenia,alpha
   REAL(KIND=rsh) :: maxquotaN,maxquotaP,excretkarenia
#endif   

 
#ifdef key_phaeocystis
   ! variables used in include incellwat_bloom_phaeocystis
   ! ------------------------------------------------------
   REAL(KIND=rsh) :: f
   INTEGER        :: npas
   REAL(KIND=rsh) :: effetchaleurphaeocystis,excretphaeocystis,phaeocystis_iksmith
   REAL(KIND=rsh) :: effetno3phaeocystiscolo,effetnh4phaeocystiscolo,effetazotephaeocystiscolo,effetphosphorephaeocystiscolo, &
                     pslimphaeocystiscolo,effetno3phaeocystiscell,effetnh4phaeocystiscell,effetazotephaeocystiscell, & 
                     effetphosphorephaeocystiscell,pslimphaeocystiscell
   REAL(KIND=rsh) :: mucusprod,Ctomucusratio,effetsalinitephaeo
#endif   

#ifdef key_benthos
  ! variables used in include incellwat_bloom_benthos
   ! ---------------------------------------------------
   REAL(KIND=rsh) :: ufond,detresuspsil,detresusp
   REAL(KIND=rsh) :: erodvitcrit,reminbenth
   REAL(KIND=rsh) :: diatmortsed,factvitessmax
#if defined key_diatbenth || defined key_zostera
   REAL(KIND=rsh) :: diatresusp
#endif
#if defined key_MANGAbio && defined key_MANGAbiovague
   REAL(KIND=rsh)           :: ubr_interp,vbr_interp,ufond_vague
#endif
#if defined key_NPbenth
   REAL(KIND=rsh)           :: diffusion_PO4,diffusion_NH4
   REAL(KIND=rsh),PARAMETER :: p_kzf = 0.00000005_rsh      ! coef. de diffusion du dissous eau/sed (m2.s-1)

   ! "mailles" de sediment 
   ! ---------------------
   REAL(KIND=rsh),PARAMETER :: epnsed=0.4_rsh ! epaisseur (m) de la "maille de sediment" (bricole)
   REAL(KIND=rsh),PARAMETER :: deltaz = 0.3_rsh    ! distance (m) fixe entre la maille de fond et le "sediment" (bricole)
#endif
#if defined key_zostera
   REAL(KIND=rsh) :: vitessechuteseed,vitessechutezostdet,seedresusp,detNzostresusp,detPzostresusp
   REAL(KIND=rsh),PARAMETER :: p_kzf = 0.00000005_rsh      ! coef. de diffusion du dissous eau/sed (m2.s-1)
   REAL(KIND=rsh) :: diffusion_PO4,diffusion_NH4
#endif
#endif

#ifdef key_zostera
   ! variables used in include incellwat_bloom_zostera
   ! ---------------------------------------------------
   REAL(KIND=rsh) :: LAI,Qcan,alumfond,vitvent
   REAL(KIND=rsh) :: effetchaleurzprod,effetchaleurzrespf,effetchaleurzrespr,effetchaleurzmort
   REAL(KIND=rsh) :: effetchaleurzrecr,effet_vent
   REAL(KIND=rsh) :: LNquota,RNquota,LNsat,RNsat,limabsLnh4,limabsLno3,limabsRnh4
   REAL(KIND=rsh) :: limprodLN,limprodLP,limnutrecr,limRN,limRP
   REAL(KIND=rsh) :: LPquota,RPquota,lPsat,RPsat,limabsLpo4,limabsRpo4
   REAL(KIND=rsh) :: effetlumzost,Net_Prod_zost,Root_growth,Leaf_growth,Kzost
   REAL(KIND=rsh) :: Labs_nh4,Labs_no3,Rabs_nh4,Labs_po4,Rabs_po4
   REAL(KIND=rsh) :: deltasat_N,deltasat_P,trans_N_acro,trans_N_basi,trans_P_acro,trans_P_basi
   REAL(KIND=rsh) :: LNrecla,LPrecla,RNrecla,RPrecla
   REAL(KIND=rsh) :: Root_mort,Leaf_mort
   REAL(KIND=rsh) :: revolangle,declinangle,dureej,limdureej,limselfshad,limRB,Recruit_rate
   REAL(KIND=rsh) :: flowering_rate,SNquota,SPquota,Seed_mort,limLBgerm,effetchaleurgerm,Germin_rate
#endif 

#ifdef key_oxygen  
   ! variables used in include incellwat_bloom_oxygen
   ! ---------------------------------------------------
   REAL(KIND=rsh) :: ophotos,oresphyto,orespzoo,oremineau,onitrifeau, &
                     gst,o2sat,oair,perteo2,schmidt
   REAL(KIND=rsh) :: effetnutdiat,effetnutdino,effetnutnano, &
                     effetnutphaeocell,effetnutphaeocolo,effetnutkarenia,effetnutpsnz

   ! le ratio theorique permettant de passer de micromolN.l-1 en mgO2.l-1 vaut :
   ! ratio catabolique :  0.001*6.625*32=0.212  
   !                       si NH4 est le carburant azote (vaut pour la respiration et la remineralisation)
   ! ration anabolique :  0.001*6.625*(1+1.5*16/106)*32=0.260  
   !                       si NO3 est le carburant azote (Neumann,2000 :modele Baltique)
   !                         (vaut pour la photosynthese)
   REAL(KIND=rsh), PARAMETER :: ratio_mgO_to_mumolN_catabol=0.212, &
                                ratio_mgO_to_mumolN_anabol=0.260
#endif   
  
#ifdef key_oyster_SFG  
   ! variables used in include incellwat_bloom_oysterSGF
   ! ---------------------------------------------------
   REAL(KIND=rlg) :: tjour
   REAL(KIND=rsh) :: colmat,effettemper,filtmes,filt,filt1,filt2
   REAL(KIND=rsh) :: cphyto,cdetri,retention,consochla,consodet,consoi
   REAL(KIND=rsh) :: consochla2,consodet2,consoNdet2,consoPdet2,volmaille
   REAL(KIND=rsh) :: val1,PFi,cingesi,susp,repsom,susp0,a,b
   REAL(KIND=rsh) :: cingestot,pseudet,ORGAFRING
   REAL(KIND=rsh) :: ER,BILAN,MAXSOMA,BILANSOM,BILANREP,EAT,pseudi
   REAL(KIND=rsh) :: psfchla2,psfdet2,psfNdet2,psfPdet2,egestotNhui,egestotPhui
   REAL(KIND=rsh) :: cmop,consorg,pseuorg,EAorg,ABSorg,fecesorg,ABSorg1
   REAL(KIND=rsh) :: ORGCOQUILLE,GAINSOMA,GAINREPROD,GAINCOQUILLE,PONTE,GAINSOMA2,GAINREPROD2
   REAL(KIND=rsh) :: BILANCOQ,BILANSOMA,RESP1,RESP,PForg,COQGAM,CINGESORG
   INTEGER        :: GAMETO,REPOSGAMETE,REPOSCOQUILLE,SEUILPONTE
#endif

#ifdef key_oyster_DEB  
   ! variables used in include incellwat_bloom_oysterDEB
   ! ---------------------------------------------------
   REAL(KIND=rsh) :: facttemp,tempref
   REAL(KIND=rsh) :: NRJ_reserves,chloro
   REAL(KIND=rsh) :: NRJ_structures,NRJ_gonades,pa,pc,pg
   REAL(KIND=rsh) :: amaig,pr,ponte,lyse,DWres,DWstruct,DWgon,DWtot
   REAL(KIND=rsh) :: Lgtot,ER,FRpop,PFpop,Fpop,pxpop,MES
   REAL(KIND=rsh) :: k_chlvble,cohnbmoule,pm,pj,txing,txfilt,ChlC,vol
   REAL(KIND=rsh) :: filtration,Tempcorr,pxm,pam,fchl
   REAL(KIND=rsh) :: energ,seuil,pam_max
   REAL(KIND=rsh) :: untier,deuxtier,FR,PF,Feces
   REAL(KIND=rsh) :: lysemax,amaigmax,px
#endif

#ifdef key_N_tracer  
   ! variables used in include incellwat_bloom_nitrogentracer
   ! ---------------------------------------------------
   INTEGER                                   :: ivtra,iv_tra,iv_sign,id_sign,ivage,iso
   REAL(KIND=rlg),DIMENSION(nv_state)        :: signature_N,age_N 
   REAL(KIND=rlg)                            :: age_max 
   REAL(KIND=rsh)                            :: Nphytototal,phytototal_tra_N 
   REAL (KIND=rsh), PARAMETER                :: seuil_N_tracer=1.e-6_rsh,seuil_N_age_tracer=0.01_rsh
   REAL (KIND=rsh), PARAMETER                :: seuil_N_tracer_diag=seuil_N_tracer*500.0_rsh,seuil_N_age_tracer_phytot=0.02_rsh
#endif

#ifdef key_P_tracer  
   ! variables used in include incellwat_bloom_phosphoretracer
   ! ---------------------------------------------------
   INTEGER                                   :: ivtra,iv_tra,iv_sign,id_sign,ivage,iso
   REAL(KIND=rlg),DIMENSION(nv_state)        :: signature_P,age_P 
   REAL(KIND=rlg)                            :: age_max 
   REAL(KIND=rsh)                            :: Pphytototal,phytototal_tra_P
   REAL (KIND=rsh), PARAMETER                :: seuil_N_tracer=1.e-6_rsh,seuil_N_age_tracer=0.01_rsh
   REAL (KIND=rsh), PARAMETER                :: seuil_N_tracer_diag=5.e-4_rsh,seuil_N_age_tracer_phytot=0.02_rsh 
   REAL (KIND=rsh)                           :: seuil_P_tracer,seuil_P_age_tracer,seuil_P_tracer_diag,seuil_P_age_tracer_phytot
#endif

#if ! defined key_BLOOM_opt2
  ! variables used in include incellwat_bloom_settling
   ! ----------------------------------------
   REAL(KIND=rsh)                            :: modulationsedpardenseau,rdet,wchutedet
   REAL(KIND=rsh)                            :: phytodetritus,zoodetritus,effetchutedia,abovesinkingrate
#endif   

#ifdef MUSTANG
   REAL(KIND=rsh)  :: tocdpe,depo
#endif

  !! * Executable part

   dtbiojour=REAL(dtbio/86400.0_rlg,rsh)
   foncsinus=(1.0_rsh+SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLINK-125)))/2.0_rsh
   
#if defined MUSTANG && (! defined key_Pconstitonly_insed || defined key_BLOOM_insed)
   txfiltbenth=0.0_rsh
#elif ! defined key_BLOOM_opt2
   IF(l_filtbenthsinus) THEN
     ! Introduction d une pression de Broutage, d apres Marie SAVINA
     ! -------------------------------------------------------------  
     ! txfiltbenth = taux sinusoidal de filtration du benthos en m3/j/m2
     txfiltbenth=(p_txfiltbenthmax*(0.3_rsh+0.7_rsh*foncsinus))
   ELSE
     txfiltbenth=p_txfiltbenthmax
   ENDIF
#else
  txfiltbenth=0.0_rsh
#endif

!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
#ifdef key_MARS
    DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
    DO i=ifirst,ilast
#endif
      IF (htot(i,j) > RESIDUAL_THICKNESS_WAT)  THEN
            
#ifdef key_MARS
       IF(htot(i,j) < hm) THEN
          ! 3D to 1D for variable cv_wat(nb_var,NB_LAYER_WAT,limin:limax,ljmin:ljmax)
          DO iv=1,nv_adv
              c(iv)=SUM(cvadv_wat_pos(iv,:,i,j)*dsigu(:))
          ENDDO
          DO iv=1,nv_fix
              c(nv_adv+iv)=0.0_rsh
              DO k=1,NB_LAYER_WAT
                IF(cvfix_wat_pos(iv,k,i,j) > c(nv_adv+iv))c(nv_adv+iv)=cvfix_wat_pos(iv,k,i,j)
              ENDDO
          ENDDO
#ifdef key_benthic
          DO iv=1,nv_bent
              c(nv_adv+nv_fix+iv)=cv_bent_pos(iv,i,j)
          ENDDO
#endif
          !c(1:nv_state)=0.5_rsh*( c(1:nv_state)+   &
          !                     ABS(c(1:nv_state)) )
          c(1:nv_state)=c(1:nv_state)/unit_modif_mudbio_N2dw(irk_fil(1:nv_state))
          kkmax=1
       ELSE
          kkmax=NB_LAYER_WAT
       ENDIF
#else
       kkmax=NB_LAYER_WAT
#endif
       !!!!!!!!!!!!!!!!!!!!!!
       !!    loop on k     !!
       !!!!!!!!!!!!!!!!!!!!!!
       DO k=1,kkmax
       
         ! initialize to 0
         dc(1:nv_state)=0.0_rsh

         ! test pour identification des mailles de fond
         ! --------------------------------------------
         ibbenth=0
         IF (k.eq.1) ibbenth=1

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
#ifdef key_benthic
          DO iv=1,nv_bent
              c(nv_adv+nv_fix+iv)=cv_bent_pos(iv,i,j)/unit_modif_mudbio_N2dw(irk_fil(iv))
          ENDDO
#endif
         ENDIF
           
   ! test on salinity (if < p_sali_thhold_bio ==> no bio sink-sources)
         sali=max(0.0_rsh,SAL_BIOLink(k,i,j))
         IF(sali < p_sali_thhold_bio) THEN
         
          dc(:)=0.0_rsh
          
         ELSE
         
          epn=thicklayerW_C(k,i,j)

          ! temperature & effetchaleur
          tempabs=TEMP_BIOLink(k,i,j)+273.15_rsh
          temper=max(0.0_rsh,TEMP_BIOLink(k,i,j))
   ! 
          IF (temper > 30.0_rsh) THEN
             write(*,*) '!!!!!tempwat>30 !!!!!!! (forcee a 30degres)',ijour_BIOLINK,imois_BIOLINK,iheure_BIOLINK,k,i,j,temper
             temper=30.0_rsh
             tempabs=temper+273.15_rsh                 
          ENDIF
          effetchaleur=exp(p_T_effect*temper)
#ifdef key_physadaptation
          effetturbidite=EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct/   &
                    (sqrt(1+(EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct)**2))
#else
          IF(l_ChlNratio_var) THEN
            effetturbidite=EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct/   &
                     (sqrt(1+(EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct)**2))
          ELSE
            effetturbidite=1.0_rsh
          ENDIF
#endif


#if defined MUSTANG && (! defined key_Pconstitonly_insed || defined key_BLOOM_insed)
          txfiltbenthij=0.0_rsh
#elif ! defined key_BLOOM_opt2  
          IF(l_filtbenthmes) THEN
!         !  txfiltbenth fonction des MES (filtration baisse avec la turbidite)
            txfiltbenthij=MAX(-0.0333_rsh*cmes_3dmgl(k,i,j)+1,0.0_rsh)*txfiltbenth  !de 0 mg/L a 30 mg/L
          ELSE
            txfiltbenthij=txfiltbenth
          ENDIF
#else
          txfiltbenthij=txfiltbenth
#endif

!        Filtrage benthique par huitres/moules/crepidules
!         ==================================================
          IF (ibbenth.eq.1) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!     HUITRES     !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_oyster_benthos
           IF(nbhuitre(i,j).ne.0.0_rsh) THEN

              cmes1=c(iv_spim)*1000.0_rsh  !c passe en mg/l
              colmat_huitre=min(0.,(p_huitre_thr_colmat-cmes1))        
              effettemp_huitre=p_huitre_temp1*(temper-p_huitre_temp_opt)**2
              IF (cmes1.gt.p_huitre_thr_mes) THEN
                filtmes=p_huitre_loi_mes1*cmes1+p_huitre_loi_mes2
              ELSE
                filtmes=p_huitre_filt_std
              ENDIF
              filt=filtmes-effettemp_huitre
              filt=filt*(p_huitre_DW**p_huitre_exp_allom)
              filt=filt*exp(p_huitre_loi_colmat*colmat_huitre)
              filt=max(filt,0.)
              txfiltbenthij=filt*nbhuitre(i,j)*24.0_rsh*1.e-3/CELL_SURF(i,j)
              write(*,*)'i,j,filt,txfiltbenthij=',i,j,filt,txfiltbenthij
            ENDIF
#endif
            IF(id_benthos_txf .NE. 0)diag_2d(irk_diag(id_benthos_txf),i,j)=txfiltbenthij/epn
          ENDIF

#ifdef key_ulvas
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Croissance et mortalite du phytoplancton de base
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF (c(iv_nutr_NO3) .ne. 0.0_rsh .or. c(iv_nutr_NH4) .ne. 0.0_rsh) THEN
           fractionno3=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+c(iv_nutr_NH4))
           fractionnh4=c(iv_nutr_NH4)/(c(iv_nutr_NO3)+c(iv_nutr_NH4))
          ELSE
           fractionno3=0.0_rsh
           fractionnh4=0.0_rsh
          END IF
#endif

    ! DIATOMEES
    ! ---------
          ! Effetlumiere : Formulation de Smith integree sur la profondeur de la boite
          ! ............................................................
#ifdef key_physadaptation
#ifdef key_MANGAbio
          ! effet de la latitude en plus
          diat_iksmith=p_diat_iksmith*(1.0_rsh+0.5_rsh*SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLINK-100)))  &
                *max(0.5_rsh,(1.0_rsh-effetturbidite ))*                                   &
                (2.0_rsh*(53.0_rsh-lat2d(i,j))+0.5_rsh*(lat2d(i,j)-43.0_rsh))/10.0_rsh
#else
          ! Cas ou Ik varie selon la turbidite et la saison
          diat_iksmith=p_diat_iksmith*(1.0_rsh+0.5_rsh*SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLINK-100)))  &
                *max(0.5_rsh,(1.0_rsh-effetturbidite))  
#endif
#else
          diat_iksmith=p_diat_iksmith
#endif
          fluxrelatifsurf=PAR_top_layer(k,i,j)/(diat_iksmith+0.0000000001_rsh)
          fluxrelatiffond=PAR_top_layer(k-1,i,j)/(diat_iksmith+0.0000000001_rsh)
          effetlumierediat=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*                     &
                  log((fluxrelatifsurf+sqrt(1.0_rsh+fluxrelatifsurf*                  &
          fluxrelatifsurf))/(fluxrelatiffond+                                 &
                  sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))
          
!          IF ((effetlumierediat.lt.0.0_rsh).or.(effetlumierediat.gt.1.0_rsh)) &
          effetnitdiat=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_diat_kNO3+               &
                   (c(iv_nutr_NH4)*p_diat_kNO3/p_diat_kNH4))
          effetamdiat=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_diat_kNH4+                &
                   (c(iv_nutr_NO3)*p_diat_kNH4/p_diat_kNO3))
          effetazotediat=effetnitdiat+effetamdiat
#if ! defined key_BLOOM_opt2
          Si=max(0.0_rsh,c(iv_nutr_SiOH))
          effetsilice=Si/(Si+p_diat_kSi)
#else
          effetsilice=c(iv_nutr_SiOH)/(c(iv_nutr_SiOH)+p_diat_kSi)
#endif

          effetphosphorediat=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_diat_kPO4)

          effetselnutdiat=min(effetazotediat,effetsilice,effetphosphorediat)

          pslimdiat=min(effetlumierediat,effetselnutdiat)
          rationdiat=p_diat_mumax*effetchaleur*pslimdiat

          !excretdiat=p_phyto_resp*effetchaleur*(1.0_rsh-effetlumierediat)   
          excretdiat=0.0_rsh
   
!#if defined key_BLOOM_opt2
          diatmorteau=p_diat_mort*effetchaleur
!#else
!         diatmorteau=p_diat_mortmax*(1-effetselnutdiat**0.2)
!#endif
          IF(c(iv_phyto_diat_N).le.p_diat_thhold_mort) diatmorteau=0.0_rsh


   ! DINOFLAGELLES
   ! -------------

          ! Effetlumiere : Formulation de Smith integree sur la profondeur de la boite
          ! ............................................................
#ifdef key_physadaptation
          ! Cas ou Ik varie selon la turbidite  
          dino_iksmith=p_dino_iksmith*max(0.5_rsh,(1.0_rsh-effetturbidite))
#else
          dino_iksmith=p_dino_iksmith
#endif
          fluxrelatifsurf=PAR_top_layer(k,i,j)/(dino_iksmith+0.0000000001_rsh)
          fluxrelatiffond=PAR_top_layer(k-1,i,j)/(dino_iksmith+0.0000000001_rsh)
          effetlumieredino=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*                     &
                  log((fluxrelatifsurf+sqrt(1.0_rsh+fluxrelatifsurf*                  &
          fluxrelatifsurf))/(fluxrelatiffond+                                 &
                  sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))

!          IF ((effetlumieredino.lt.0.0_rsh).or.(effetlumieredino.gt.1.0_rsh)) &

          effetnitdino=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_dino_kNO3+                               &
                   (c(iv_nutr_NH4)*p_dino_kNO3/p_dino_kNH4))
          effetamdino=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_dino_kNH4+                                &
                   (c(iv_nutr_NO3)*p_dino_kNH4/p_dino_kNO3))
          effetazotedino=effetnitdino+effetamdino 
          effetphosphoredino=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_dino_kPO4)
          effetselnutdino=min(effetazotedino,effetphosphoredino)
          pslimdino=min(effetlumieredino,effetselnutdino)
          rationdino=p_dino_mumax*effetchaleur*pslimdino

#if ! defined key_BLOOM_opt2
          ! taux max de croissance diminue dans les zones a forte ECT(croissance nulle pour ect>=p_dino_thhold_ect m2.s-2)
          rationdino=rationdino*max(0.0_rsh,(1.0_rsh-ECT_kij/p_dino_thhold_ect))
#endif

          !excretdino=p_phyto_resp*effetchaleur*(1.0_rsh-effetlumieredino)   
          excretdino=0.0_rsh
          dinomorteau=p_dino_mort*effetchaleur

          IF(c(iv_phyto_dino_N).le.p_dino_thhold_mort) dinomorteau=0.0_rsh

   ! NANOPICOPLANCTON
   ! ----------------
          ! Effetlumiere : Formulation de Smith integree sur la profondeur de la boite
          ! ............................................................
#ifdef key_physadaptation
          nano_iksmith=p_nano_iksmith*max(0.5_rsh,(1.0_rsh-effetturbidite))
#else
          nano_iksmith=p_nano_iksmith
#endif
          fluxrelatifsurf=PAR_top_layer(k,i,j)/(nano_iksmith+0.0000000001_rsh)
          fluxrelatiffond=PAR_top_layer(k-1,i,j)/(nano_iksmith+0.0000000001_rsh)
          effetlumierenano=1.0_rsh/epn/EXTINCTION_RAD(k,i,j)*                     &
                  log((fluxrelatifsurf+sqrt(1.0_rsh+fluxrelatifsurf*                  &
            fluxrelatifsurf))/(fluxrelatiffond+                                 &
                  sqrt(1.0_rsh+fluxrelatiffond*fluxrelatiffond)))

!          IF ((effetlumierenano.lt.0.0_rsh).or.(effetlumierenano.gt.1.0_rsh)) &


          effetnitnano=c(iv_nutr_NO3)/(c(iv_nutr_NO3)+p_nano_kNO3+                                 &
                     (c(iv_nutr_NH4)*p_nano_kNO3/p_nano_kNH4))
          effetamnano=c(iv_nutr_NH4)/(c(iv_nutr_NH4)+p_nano_kNH4+                                &
                     (c(iv_nutr_NO3)*p_nano_kNH4/p_nano_kNO3))
          effetazotenano=effetnitnano+effetamnano

          effetphosphorenano=c(iv_nutr_PO4)/(c(iv_nutr_PO4)+p_nano_kPO4)
          effetselnutnano=min(effetazotenano,effetphosphorenano)
          pslimnano=min(effetlumierenano,effetselnutnano)
          rationnano=p_nano_mumax*effetchaleur*pslimnano

          IF (effetazotediat.gt.0.000001_rsh) THEN
           fractionno3diat=effetnitdiat/effetazotediat
           fractionnh4diat=effetamdiat/effetazotediat
          ELSE
           fractionno3diat=0.0_rsh
           fractionnh4diat=0.0_rsh
          ENDIF
          IF (effetazotedino.gt.0.000001_rsh) THEN
           fractionno3dino=effetnitdino/effetazotedino
           fractionnh4dino=effetamdino/effetazotedino
          ELSE
           fractionno3dino=0.0_rsh
           fractionnh4dino=0.0_rsh
          ENDIF
          IF (effetazotenano.gt.0.000001_rsh) THEN
           fractionno3nano=effetnitnano/effetazotenano
           fractionnh4nano=effetamnano/effetazotenano
          ELSE
           fractionno3nano=0.0_rsh
           fractionnh4nano=0.0_rsh
          ENDIF

          excretnano=0.0_rsh
          nanomorteau=p_nano_mort*effetchaleur
          IF (c(iv_phyto_nano_N).lt.p_nano_thhold_mort) nanomorteau=0.0_rsh

   !+++++++++++++++++++++++++++++++++++++++++++++
   ! Croissance et mortalite du zooplancton
   !+++++++++++++++++++++++++++++++++++++++++++++

   ! MESOZOOPLANCTON
   ! ----------------
          ! concentrations de proies visibles 
#if defined key_BLOOM_opt2
          ! limitation excretion par rapport au broutage (virg 17sept2011) entre autre sur les mcrozoo

          IF (c(iv_phyto_diat_N) .gt. 0.01_rsh) THEN
            captdiatmeso=c(iv_phyto_diat_N)*p_mesz_captdiat        !en microMolN.l-1
          ELSE
            captdiatmeso=0.0_rsh
          ENDIF
          IF (c(iv_phyto_dino_N) .gt. 0.01_rsh) THEN 
            captdino=c(iv_phyto_dino_N)*p_mesz_captdino        !en microMolN.l-1
          ELSE
            captdino=0.0_rsh
          ENDIF
          rNC=1.0_rsh/(12.0_rsh*p_zoo_CNratio)                 !en micromolN.microgC-1
          ! dans l option 2 le zoo est exprime en carbone
          IF((c(iv_zoo_micr_N)*rNC) .gt. 0.01_rsh) THEN
            captmicrozoo=c(iv_zoo_micr_N)*p_mesz_captmicz      !en microMolN.l-1
          ELSE
            captmicrozoo=0.0_rsh
          ENDIF
          captmicrozooN=captmicrozoo*rNC                       ! en microMolN.l-1
          totproiebrut=captdiatmeso+captdino+captmicrozooN
#else
          ! dans l option 1 le zoo est exprime en azote
          captdiatmeso=c(iv_phyto_diat_N)*p_mesz_captdiat    !en microMolN.l-1
          captdino=c(iv_phyto_dino_N)*p_mesz_captdino        !en microMolN.l-1
          captmicrozoo=c(iv_zoo_micr_N)*p_mesz_captmicz      !en microMolN.l-1
          rNC=1.0_rsh
          captmicrozooN=captmicrozoo*rNC                     ! en microMolN.l-1
          totproiebrut=captdiatmeso+captdino+captmicrozoo
          p_mesz_thrN_turb=p_mesz_thrN*effetturbidite
#endif

#ifdef key_psnz
          captpsnz=0.0_rsh                                        
          IF (c(iv_phyto_psnz_N) .gt. 0.1_rsh) captpsnz=c(iv_phyto_psnz_N)*p_mesz_captpsnz
          totproiebrut=totproiebrut+captpsnz
#endif
#ifdef key_phaeocystis   
          captphaeocolo=c(iv_phyto_phaeocystis_colo_N)*p_mesz_captphaeocolo   ! microMolN.l-1
          totproiebrut=totproiebrut+captphaeocolo
#endif
#ifdef key_karenia
          captkarenia=c(iv_phyto_karenia_N)*p_mesz_captkarenia
          totproiebrut=totproiebrut+captkarenia
#endif


!          IF (captdiatmeso.ne.0.0_rsh .or. captdino.ne.0.0_rsh .or. captmicrozooN.ne.0.0_rsh) THEN
#if defined key_BLOOM_opt2
          totproie=totproiebrut
          IF (totproiebrut .gt.0.005_rsh) THEN
#else
          totproie=max(0.0_rsh,(totproiebrut-p_mesz_thrN_turb)) ! en microgMolN.l-1
          IF (totproie .ne. 0.0_rsh) THEN
#endif

       ! Cas ou le parametre d Ivlev varie selon la turbidite
#ifdef key_physadaptation
            IF (cmes_3dmgl(k,i,j).le.p_mesz_thhold_mes_kivlev) THEN
              meso_kivlev=p_mesz_kivlev
            ELSE
             meso_kivlev=p_mesz_kivlev*min(1.0_rsh,(1.0_rsh-effetturbidite))
            ENDIF
#else
            meso_kivlev=p_mesz_kivlev
#endif

            fmesozoo=1.0_rsh-exp(-meso_kivlev*totproie)              ! food lim [0-1]
            rationmesozoo=p_mesz_mumax*effetchaleur*fmesozoo     ! specif growth rate d-1
#if defined key_BLOOM_opt2
!           Cas ou le taux d assimilation varie selon le taux de repletion  
!           assimilation depend de quantitee ingeree Hofmann & Hambler, 1988 J.M.R  
            assimilmesozoo=0.3_rsh*(3.0_rsh-0.67_rsh*fmesozoo)       ! ass rate s.u.
#else
            assimilmesozoo=0.3_rsh*(3.0_rsh-0.67_rsh*fmesozoo)       ! ass rate s.u.
#endif
            broumesozoodiat=rationmesozoo*captdiatmeso/totproiebrut                 ! d-1
            broumesozoodino=rationmesozoo*captdino/totproiebrut                 ! d-1
            ! captmicrozooN=captmicrozoo dans option 1 : converti en azote dans l option2
            broumesozoomicrozoo=rationmesozoo*captmicrozooN/totproiebrut         ! d-1
#ifdef key_psnz   
            broumesozoopsnz=rationmesozoo*captpsnz/totproiebrut                 ! d-1   
#endif
#ifdef key_karenia   
            broumesozookarenia=rationmesozoo*captkarenia/totproiebrut                 ! d-1   
#endif
#ifdef key_phaeocystis   
            broumesozoophaeocolo=rationmesozoo*captphaeocolo/totproiebrut       ! d-1   
#endif

          ELSE
            totproiebrut=0.0_rsh
            totproie=0.0_rsh
            fmesozoo=0.0_rsh
            assimilmesozoo=0.0_rsh
            rationmesozoo=0.0_rsh
            broumesozoodiat=0.0_rsh
            broumesozoodino=0.0_rsh
            broumesozoomicrozoo=0.0_rsh
#ifdef key_psnz   
            broumesozoopsnz=0.0_rsh   
#endif
#ifdef key_karenia   
            broumesozookarenia=0.0_rsh                 ! d-1   
#endif
#ifdef key_phaeocystis   
            broumesozoophaeocolo=0.0_rsh   
#endif
          ENDIF

          excretionmesozoo=p_mesz_excret*effetchaleur*fmesozoo ! specif excret rate d-1
          txmortmesozoo=(p_mesz_mort1+p_mesz_mort2*c(iv_zoo_meso_N))*effetchaleur ! specif mort rate d-1

          IF (c(iv_zoo_meso_N).lt.p_mesz_thhold_mort)  txmortmesozoo=0.0_rsh


   ! MICROZOOPLANCTON
   ! ----------------
          ! concentrations de proies visibles 
#if defined key_BLOOM_opt2
          ! limitation excretion par rapport au broutage (virg 17sept2011) entre autre sur les mcrozoo
          IF (c(iv_phyto_nano_N) .gt. 0.01_rsh) THEN 
            captnano=c(iv_phyto_nano_N)*p_micz_captnano                     !en microMolN.l-1
          ELSE
            captnano=0.0_rsh
          ENDIF
          IF (c(iv_detr_N).gt.0.01_rsh) THEN
            captdet=min(c(iv_detr_N),c(iv_detr_P)*p_phyto_NPratio)*p_micz_captdet     !en microMolN.l-1
          ELSE
           captdet=0.0_rsh 
          ENDIF
          IF(c(iv_phyto_diat_N) .gt. 0.01_rsh) THEN
            captdiatmicro=c(iv_phyto_diat_N)*p_micz_captdiat    ! microMolN.l-1
          ELSE
           captdiatmicro=0.0_rsh
          ENDIF
          captdinomicro=0.0_rsh
          totproiebrut=(captnano+captdet+captdiatmicro)
#else
          captnano=c(iv_phyto_nano_N)*p_micz_captnano                               !en microMolN.l-1
          captdet=min(c(iv_detr_N),c(iv_detr_P)*p_phyto_NPratio)*p_micz_captdet     !en microMolN.l-1
          captdinomicro=c(iv_phyto_dino_N)*p_micz_captdino                          !en microMolN.l-1

#if defined key_physadaptation
          ! Modulation de la capturabilite des diatomees par la profondeur (grosses especes cotieres, petites oceaniques)
          captdiatmicro=c(iv_phyto_diat_N)*p_micz_captdiat*(0.1_rsh+0.9_rsh*min(200.0_rsh,BATHY_H0(i,j))/200.0_rsh)
#else
          captdiatmicro=c(iv_phyto_diat_N)*p_micz_captdiat                          ! microMolN.l-1
#endif	        
          totproiebrut=captdiatmicro+captnano+captdet+captdinomicro
#endif

#ifdef key_karenia
          captkarenia=c(iv_phyto_karenia_N)*p_micz_captkarenia
          totproiebrut=totproiebrut+captkarenia
#endif

#ifdef key_psnz
#if defined key_BLOOM_opt2
          captpsnz=c(iv_phyto_psnz_N)*p_micz_captnano
#else
          IF(c(iv_phyto_psnz_N) .gt. 0.1_rsh) THEN
            captpsnz=c(iv_phyto_psnz_N)*p_mesz_captpsnz
          ELSE
            captpsnz=0.0_rsh
          ENDIF    
#endif
          totproiebrut=totproiebrut+captpsnz
#endif
#ifdef key_phaeocystis   
          IF(c(iv_phyto_phaeocystis_cell_N) .gt. 0.01_rsh) THEN
            captphaeocell=c(iv_phyto_phaeocystis_cell_N)*p_micz_captnano               ! microMolN.l-1
          ELSE
            captphaeocell=0.0_rsh
          ENDIF   
          totproiebrut=totproiebrut+captphaeocell
#endif
#if defined key_BLOOM_opt2
          totproie=totproiebrut
          IF (totproiebrut .gt.0.005_rsh) THEN
#else
          totproie=max(0.0_rsh,(totproiebrut-p_micz_thrnano)) ! en microgMolN.l-1
          IF (totproie .ne. 0.0_rsh) THEN
#endif
             fmicrozoo=totproie/(p_micz_kgraz+totproie)                              ! food lim [0-1]
             rationmicrozoo=p_micz_mumax*effetchaleur*fmicrozoo      ! specif growth rate d-1
#if defined key_BLOOM_opt2
!             version  assimilation constante
             assimilmicrozoo=p_micz_assim
#else
!             assimilation depend de la quantite ingeree Hofmann & Hambler, 1988 J.M.R  
             assimilmicrozoo=0.3_rsh*(3.0_rsh-0.67_rsh*fmicrozoo)        ! ass rate s.u.
#endif
             broumicrozoonano=rationmicrozoo*captnano/totproiebrut                   ! d-1
             broumicrozoodet=rationmicrozoo*captdet/totproiebrut                     ! d-1
             broumicrozoodiat=rationmicrozoo*captdiatmicro/totproiebrut                   ! d-1
             broumicrozoodino=rationmicrozoo*captdinomicro/totproiebrut                   ! d-1
#ifdef key_karenia 
             broumicrozookarenia=rationmicrozoo*captkarenia/totproiebrut           ! d-1 
#endif
#ifdef key_psnz   
             broumicrozoopsnz=rationmicrozoo*captpsnz/totproiebrut                 ! d-1   
#endif
#ifdef key_phaeocystis   
             broumicrozoophaeocell=rationmicrozoo*captphaeocell/totproiebrut         ! d-1   
#endif
          ELSE
             fmicrozoo=0.0_rsh
             assimilmicrozoo=0.0_rsh
             rationmicrozoo=0.0_rsh
             broumicrozoonano=0.0_rsh
             broumicrozoodet=0.0_rsh
             broumicrozoodiat=0.0_rsh
             broumicrozoodino=0.0_rsh
#ifdef key_karenia
             broumicrozookarenia=0.0_rsh
#endif
#ifdef key_psnz   
             broumicrozoopsnz=0.0_rsh   
#endif
#ifdef key_phaeocystis   
             broumicrozoophaeocell=0.0_rsh   
#endif
          ENDIF

          excretionmicrozoo=p_micz_excret*effetchaleur*fmicrozoo ! specif excret rate d-1
          txmortmicrozoo=(p_micz_mort+p_mesz_mort2*c(iv_zoo_micr_N))*effetchaleur! specif excret rate d-1

          IF (c(iv_zoo_micr_N).lt.p_micz_thhold_mort) txmortmicrozoo=0.0_rsh 

#if defined key_oxygen
          IF (c(iv_oxygen).lt.0.2_rsh) THEN
            txmortmicrozoo=100.0_rsh*txmortmicrozoo
            txmortmesozoo=100.0_rsh*txmortmesozoo  
          ENDIF
#endif
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Remineralisation de la matiere detritique dans l eau
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          reminazdeteau=p_N_remin*effetchaleur
          xnitrifeau=p_nitrif*effetchaleur

          reminpdeteau=p_P_remin*effetchaleur

          dissolsiliceeau=p_BSi_dissEau*effetchaleur
   
#if defined key_BLOOM_opt2
          dissolMOPeau=p_det_fragm*effetchaleur
#endif
#if defined key_oxygen
          IF (c(iv_oxygen).lt.0.2_rsh) THEN
            reminazdeteau=0.0_rsh
            xnitrifeau=0.0_rsh
            reminpdeteau=0.0_rsh  
          ENDIF
#endif


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Adsorption et desorption du phosphore dans l eau
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++

          ! adsorption
          saturmesp=max(0.0_rsh,(p_P_adsormaxspim*c(iv_spim)-c(iv_nutr_Pads)))  !c(iv_spim) en g/l
          adsorpeau=p_P_adsor*saturmesp
          cmaxadsorp=adsorpeau*c(iv_nutr_PO4)*dtbiojour
          IF (cmaxadsorp.ge.c(iv_nutr_PO4))  adsorpeau=0.9_rsh/dtbiojour

          ! desorption
          IF(c(iv_spim).gt.0.000001_rsh) THEN
            varads=min(1.0_rsh,c(iv_nutr_Pads)/(p_P_adsormaxspim*c(iv_spim)))
            !varads=c(iv_nutr_Pads)/(p_P_adsormaxspim*c(iv_spim))
          ELSE
            varads=0.0_rsh
          END IF
          desorpeau=p_P_desor*varads
          cmaxdesorp=desorpeau*c(iv_nutr_Pads)*dtbiojour
          IF (cmaxdesorp.ge.c(iv_nutr_Pads)) desorpeau=0.9_rsh/dtbiojour


#if defined key_BLOOM_opt2
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !Calcule fraction de la matiere detr azotee qui va jusqu au fond
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
          IF (k.eq.1) THEN 
             fractionDIAT_N=c(iv_phyto_diat_N)*(min((ws3(2,iv_phyto_diat_N,i,j)*dtbio/epn),1.0_rsh))
             fractionDIAT_P=fractionDIAT_N*rappaz
             fractionDIAT_Si=fractionDIAT_N*rapsiaz
             fractionDET_N=c(iv_detr_N)*(min((ws3(2,iv_detr_N,i,j)*dtbio/epn),1.0_rsh))
             fractionDET_Si=c(iv_detr_Si)*(min((ws3(2,iv_detr_Si,i,j)*dtbio/epn),1.0_rsh))
             fractionDET_P=c(iv_detr_P)*(min((ws3(2,iv_detr_P,i,j)*dtbio/epn),1.0_rsh))
          ELSE
             fractionDIAT_N=0.0_rsh
             fractionDIAT_P=0.0_rsh
             fractionDIAT_Si=0.0_rsh
             fractionDET_N=0.0_rsh
             fractionDET_Si=0.0_rsh
             fractionDET_P=0.0_rsh
          ENDIF
#endif

   !++++++++++++++++++++++++++++++++++++++++++++++++++
   ! EQUATIONS D EVOLUTION  ! concentration en jour-1
   !++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Evolution de l ammonium
   ! -----------------------

          dc(iv_nutr_NH4)=-fractionnh4diat*rationdiat*c(iv_phyto_diat_N)  &
                   -fractionnh4dino*rationdino*c(iv_phyto_dino_N)          &
                   -fractionnh4nano*rationnano*c(iv_phyto_nano_N)          &
                   -xnitrifeau*c(iv_nutr_NH4)                              &
                   +(excretionmesozoo*c(iv_zoo_meso_N)+excretionmicrozoo*c(iv_zoo_micr_N))*rNC &
#if defined key_BLOOM_opt2
                   +reminazdeteau*c(iv_diss_N)                             &                 
                   +(1.0_rsh-assimilmicrozoo)*rationmicrozoo*(1.0_rsh-p_micz_diss)*c(iv_zoo_micr_N)*rNC
#else
                   +reminazdeteau*c(iv_detr_N)
#endif


   ! Evolution du nitrate
   ! --------------------
          dc(iv_nutr_NO3)=xnitrifeau*c(iv_nutr_NH4)                  &
                   -fractionno3diat*rationdiat*c(iv_phyto_diat_N)    &
                   -fractionno3dino*rationdino*c(iv_phyto_dino_N)    &
                   -fractionno3nano*rationnano*c(iv_phyto_nano_N)
   ! Evolution de la silice dissoute
   ! -------------------------------
          dc(iv_nutr_SiOH)= -rapsiaz*rationdiat*c(iv_phyto_diat_N)   &
#if defined key_BLOOM_opt2
                    +dissolsiliceeau*c(iv_diss_Si)
#else
                    +dissolsiliceeau*c(iv_detr_Si)
#endif

   ! Evolution du phosphate dissous
   ! ------------------------------
          dc(iv_nutr_PO4)=desorpeau*c(iv_nutr_Pads)     &
                  -adsorpeau*c(iv_nutr_PO4)             &
                  -rappaz*(rationdiat*c(iv_phyto_diat_N)+rationdino*c(iv_phyto_dino_N)+rationnano*c(iv_phyto_nano_N)) &
                  +rappaz*((excretionmesozoo*c(iv_zoo_meso_N)+excretionmicrozoo*c(iv_zoo_micr_N))*rNC)                &
#if defined key_BLOOM_opt2
                  +reminpdeteau*c(iv_diss_P)            &
                  +(1.0_rsh-assimilmicrozoo)*rationmicrozoo*(1.0_rsh-p_micz_diss)*c(iv_zoo_micr_N)*rappaz*rNC
#else
                  +reminpdeteau*c(iv_detr_P)
#endif

   ! Evolution du phosphate adsorbe
   ! --------
          dc(iv_nutr_Pads)=adsorpeau*c(iv_nutr_PO4)-desorpeau*c(iv_nutr_Pads)

   ! Evolution de l azote du nano_pico_phytoplancton  (micromolN.l-1)
   ! -----------------------------------------------------------------
          dc(iv_phyto_nano_N)=c(iv_phyto_nano_N)*(rationnano-excretnano-nanomorteau)-broumicrozoonano*c(iv_zoo_micr_N)*rNC
#if ! defined key_BLOOM_opt2 && (! defined MUSTANG || defined key_Pconstitonly_insed) 
          !  pas de prise en compte du filtrage benthique ici si MUSTANG et Bio dans sedim (sera avec la bio dans les sediments)
          dc(iv_phyto_nano_N)=dc(iv_phyto_nano_N)-txfiltbenthij*ibbenth/epn*c(iv_phyto_nano_N)
#endif

   ! Evolution de l azote des diatomees
   ! ----------------------------------
          dc(iv_phyto_diat_N)=c(iv_phyto_diat_N)*(rationdiat-diatmorteau-excretdiat) &
                      -broumesozoodiat*c(iv_zoo_meso_N)*rNC                          &
                      -broumicrozoodiat*c(iv_zoo_micr_N)*rNC 
#if defined key_BLOOM_opt2
          dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)-fractionDIAT_N
#elif (! defined MUSTANG || defined key_Pconstitonly_insed) 
          !  pas de prise en compte du filtrage benthique ici si MUSTANG et Bio dans sedim (sera avec la bio dans les sediments)
          dc(iv_phyto_diat_N)=dc(iv_phyto_diat_N)-c(iv_phyto_diat_N)*txfiltbenthij*ibbenth/epn
#endif

   ! Evolution de l azote des dinoflagelles
   ! --------------------------------------
          dc(iv_phyto_dino_N)=(rationdino-dinomorteau-excretdino)*c(iv_phyto_dino_N)  &
                       -broumesozoodino*c(iv_zoo_meso_N)*rNC                          &
                       -broumicrozoodino*c(iv_zoo_micr_N)*rNC
#if ! defined key_BLOOM_opt2 && (! defined MUSTANG || defined key_Pconstitonly_insed) 
          !  pas de prise en compte du filtrage benthique ici si MUSTANG et Bio dans sedim (sera avec la bio dans les sediments)
          dc(iv_phyto_dino_N)=dc(iv_phyto_dino_N)-txfiltbenthij*ibbenth/epn*c(iv_phyto_dino_N)
#endif

   ! Evolution de l azote detritique
   ! -------------------------------
          dc(iv_detr_N)=(diatmorteau+excretdiat)*c(iv_phyto_diat_N)                                 &
                +(txmortmesozoo+(1.0_rsh-assimilmesozoo)*rationmesozoo)*c(iv_zoo_meso_N)*rNC        &
                +(txmortmicrozoo-broumicrozoodet)*c(iv_zoo_micr_N)*rNC                              &
#if defined key_BLOOM_opt2
                -dissolMOPeau*c(iv_detr_N)                                                          &
                -fractionDET_N
#else
                 +(dinomorteau+excretdino)*c(iv_phyto_dino_N)                                       &
                 +(nanomorteau+excretnano)*c(iv_phyto_nano_N)                                       &
                 -reminazdeteau*c(iv_detr_N)                                                        &
                 +(1.0_rsh-assimilmicrozoo)*rationmicrozoo*c(iv_zoo_micr_N)*rNC
#if ! defined key_benthos && ! defined MUSTANG
          dc(iv_detr_N)=dc(iv_detr_N)+txfiltbenthij*ibbenth/epn*(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))
#endif
#endif


   ! Evolution de la silice particulaire
   ! -----------------------------------
          dc(iv_detr_Si)=rapsiaz*(diatmorteau+excretdiat)*c(iv_phyto_diat_N)         &
                 +broumesozoodiat*c(iv_zoo_meso_N)*rNC*rapsiaz   &
#if defined key_BLOOM_opt2
                 -dissolMOPeau*c(iv_detr_Si)                     &
                 -fractionDET_Si
#else
                 -dissolsiliceeau*c(iv_detr_Si)                  &
                 +broumicrozoodiat*c(iv_zoo_micr_N)*rNC*rapsiaz
#if ! defined key_benthos && ! defined MUSTANG
          dc(iv_detr_Si)=dc(iv_detr_Si)+rapsiaz*txfiltbenthij*ibbenth/epn*c(iv_phyto_diat_N)
#endif  
#endif                     

   ! Evolution du phosphore detritique
   ! ---------------------------------
          dc(iv_detr_P)=(diatmorteau+excretdiat)*c(iv_phyto_diat_N)*rappaz                    &
                -broumicrozoodet*c(iv_zoo_micr_N)*rNC*rappaz                                  &
#if defined key_BLOOM_opt2
                +txmortmesozoo*c(iv_zoo_meso_N)*rNC*rappaz                                    &
                +(1.0_rsh-assimilmesozoo)*rationmesozoo*c(iv_zoo_meso_N)*rNC*rappaz           &
                +txmortmicrozoo*c(iv_zoo_micr_N)*rNC*rappaz                                   &
                -dissolMOPeau*c(iv_detr_P)                                                    &
                -fractionDET_P 
#else
                -reminpdeteau*c(iv_detr_P)                                                   &
                +rappaz*((dinomorteau+excretdino)*c(iv_phyto_dino_N)                         &
                        +(nanomorteau+excretnano)*c(iv_phyto_nano_N)                         &
                    +rNC*((txmortmesozoo+(1.0_rsh-assimilmesozoo)*rationmesozoo)*c(iv_zoo_meso_N)    &
                         +(txmortmicrozoo+(1.0_rsh-assimilmicrozoo)*rationmicrozoo)*c(iv_zoo_micr_N))) 
#if ! defined key_benthos && ! defined MUSTANG
          dc(iv_detr_P)=dc(iv_detr_P)+txfiltbenthij*ibbenth/epn*(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N))*rappaz
#endif
#endif                                                                  

   ! Evolution de l azote du mesozooplancton (microgC/l for option 2, microgN/l for option 1)
   ! -------------------------------------------------
          dc(iv_zoo_meso_N)=(rationmesozoo*assimilmesozoo-excretionmesozoo-txmortmesozoo)*c(iv_zoo_meso_N)

   ! Evolution de l azote du microzooplancton  (microgC.l-1 for option 2, microgN/l for option 1)
   ! -----------------------------------------------------------
          dc(iv_zoo_micr_N)=c(iv_zoo_micr_N)*(rationmicrozoo*assimilmicrozoo-excretionmicrozoo)       &
                     -txmortmicrozoo*c(iv_zoo_micr_N)                                          &
                     -broumesozoomicrozoo*c(iv_zoo_meso_N)

   ! Evolution de la masse de MES  (g.l-1)
   ! -------------------------------------
          dc(iv_spim)=0.0_rsh

#if defined key_BLOOM_opt2
   ! Evolution de l azote small particulaire
   ! ---------------------------------------
          dc(iv_diss_N)=(1.0_rsh-assimilmicrozoo)*rationmicrozoo*p_micz_diss*c(iv_zoo_micr_N)*rNC &
                -reminazdeteau*c(iv_diss_N)        &
                +dissolMOPeau*c(iv_detr_N)         &
                +dinomorteau*c(iv_phyto_dino_N)    &
                +nanomorteau*c(iv_phyto_nano_N)    &
                +0.25*fractionDIAT_N               &
                +0.25*fractionDET_N


   ! Evolution du N2 accumule dans la derniere couche, perdu pour le cycle N
   ! ------------------------------------------------------------------------- 
          dc(iv_diss_fond_Nitr)=0.75*fractionDIAT_N   & 
                        +0.75*fractionDET_N
                        

   ! Evolution de la silice small particulaire
   ! ------------------------------------------
          dc(iv_diss_Si)=-dissolsiliceeau*c(iv_diss_Si) &
                 +dissolMOPeau*c(iv_detr_Si)            &
                 +broumicrozoodiat*c(iv_zoo_micr_N)*rNC*rapsiaz &
                 +fractionDET_Si                        &
                 +fractionDIAT_Si

   ! Evolution du phosphore small particulaire
   ! ------------------------------------------
          dc(iv_diss_P)=-reminpdeteau*c(iv_diss_P)     &
                +(1.0_rsh-assimilmicrozoo)*rationmicrozoo*p_micz_diss*c(iv_zoo_micr_N)*rNC*rappaz &
                +nanomorteau*c(iv_phyto_nano_N)*rappaz &
                +dinomorteau*c(iv_phyto_dino_N)*rappaz &
                +dissolMOPeau*c(iv_detr_P)             &
                +fractionDET_P                         &
                +fractionDIAT_P
#endif

   ! Production carbonee cumulee des nanoflagelles de la maille (i,j,k)
   ! ---------------------------------------------------------------------
#if defined key_BLOOM_opt2
          ! Production azotee (option 2)
          dc(iv_phyto_nano_pp)=14.e-3*rationnano*c(iv_phyto_nano_N)*epn
#else
          ! Production carbonee (option 1) 
          dc(iv_phyto_nano_pp)=epn*12.e-3*p_phyto_CNratio*rationnano*c(iv_phyto_nano_N)
#endif

   ! Production cumulee des diatomees de la maille (i,j,k)
   ! --------------------------------------------------------
#if defined key_BLOOM_opt2
          ! Production azotee (option 2)
          dc(iv_phyto_diat_pp)=14.e-3*rationdiat*c(iv_phyto_diat_N)*epn
#else
          ! Production carbonee (option 1) 
          dc(iv_phyto_diat_pp)=epn*12.e-3*p_phyto_CNratio*rationdiat*c(iv_phyto_diat_N)
#endif

   ! Production carbonee cumulee des dinoflagelles de la maille (i,j,k)
   ! ---------------------------------------------------------------------
#if defined key_BLOOM_opt2
          ! Production azotee (option 2)
          dc(iv_phyto_dino_pp)=14.e-3*rationdino*c(iv_phyto_dino_N)*epn
#else
          ! Production carbonee (option 1) 
          dc(iv_phyto_dino_pp)=epn*12.e-3*p_phyto_CNratio*rationdino*c(iv_phyto_dino_N)
#endif

   ! Production cumulee du zooplancton de la maille (i,j,k)
   ! ---------------------------------------------------------
#if defined key_zoo_prod
          ! Production carbonee (option 1) 
          dc(iv_zoo_micr_ps)=epn*12.e-3*p_zoo_CNratio*rationmicrozoo*assimilmicrozoo*c(iv_zoo_micr_N)
          dc(iv_zoo_meso_ps)=epn*12.e-3*p_zoo_CNratio*rationmesozoo*assimilmesozoo*c(iv_zoo_meso_N)
#endif

  ! ++++++++++++++++++++++++++++++++++++++++++++++++
  !   ENREGISTREMENT VARIABLES DIAGNOSTIQUES
  ! ++++++++++++++++++++++++++++++++++++++++++++++++

         ! effet lumiere DIAT : somme a chaque pas de temps = moyenne a minuit pour memoriser la variable diagnostique
          effetlumiere_day_diat(k,i,j)=effetlumiere_day_diat(k,i,j)+effetlumierediat*dtbio
          diag_3d_wat(irk_diag(id_diat_limlight),k,i,j)=effetlumierediat

!          IF(iheure_BIOLINK==0 .and. iminu_BIOLINK ==0 .and. isec_BIOLINK <= dtbio) THEN
!             diag_3d_wat(irk_diag(id_diat_limlight),k,i,j)=effetlumiere_day_diat(k,i,j)/86400.0_rsh
!             effetlumiere_day_diat(k,i,j)=0.0_rsh
!          ENDIF

          diag_3d_wat(irk_diag(id_diat_limN),k,i,j)=effetazotediat
          diag_3d_wat(irk_diag(id_diat_limSi),k,i,j)=effetsilice
          diag_3d_wat(irk_diag(id_diat_limP),k,i,j)=effetphosphorediat

          ! effet lumiere DINO : somme a chaque pas de temps = moyenne a minuit pour memoriser la variable diagnostique
          effetlumiere_day_dino(k,i,j)=effetlumiere_day_dino(k,i,j)+effetlumieredino*dtbio
          diag_3d_wat(irk_diag(id_dino_limlight),k,i,j)=effetlumieredino

!          IF(iheure_BIOLINK==0 .and. iminu_BIOLINK ==0 .and. isec_BIOLINK <= dtbio) THEN
!             diag_3d_wat(irk_diag(id_dino_limlight),k,i,j)=effetlumiere_day_dino(k,i,j)/86400.0_rsh
!             effetlumiere_day_dino(k,i,j)=0.0_rsh
!          ENDIF
          diag_3d_wat(irk_diag(id_dino_limN),k,i,j)=effetazotedino
          diag_3d_wat(irk_diag(id_dino_limP),k,i,j)=effetphosphoredino

          ! effet lumiere NANO : somme a chaque pas de temps = moyenne a minuit pour memoriser la variable diagnostique
          effetlumiere_day_nano(k,i,j)=effetlumiere_day_nano(k,i,j)+effetlumierenano*dtbio
          diag_3d_wat(irk_diag(id_nano_limlight),k,i,j)=effetlumierenano
!          IF(iheure_BIOLINK==0 .and. iminu_BIOLINK ==0 .and. isec_BIOLINK <= dtbio) THEN
!             diag_3d_wat(irk_diag(id_nano_limlight),k,i,j)=effetlumiere_day_nano(k,i,j)/86400.0_rsh
!             effetlumiere_day_nano(k,i,j)=0.0_rsh
!          ENDIF

          diag_3d_wat(irk_diag(id_nano_limN),k,i,j)=effetazotenano
          diag_3d_wat(irk_diag(id_nano_limP),k,i,j)=effetphosphorenano

          ! chlorophylle [mgChl/m3]
          !  Cas ou on stocke simplement l azote phyto total
          !fact_phyto_ChlNratio=1.0_rsh
          IF(l_ChlNratio_var) THEN
              IF((CURRENT_TIME-TIME_BEGIN) > 345600.00_rlg) THEN
                 ! Cas ou on transforme l azote en chloro par un rapport Chloro/N dependant de l extinction
                 !   lumineuse moyennee sur les 4 jours precedents  et de la limitation en azote
                 fact_phyto_ChlNratio=p_phyto_ChlNratiomax*extinction_ave4d(k,i,j)/p_phyto_ChlN_ksmithextinct/  &
                             (sqrt(1+(extinction_ave4d(k,i,j)/p_phyto_ChlN_ksmithextinct)**2))
                 fact_phyto_ChlNratio=fact_phyto_ChlNratio*(0.5_rsh*(1.0_rsh-effetazotediat)+1.5_rsh*effetazotediat)   
              ELSE
                 fact_phyto_ChlNratio=p_phyto_ChlNratiomax*effetturbidite    
              ENDIF
              IF(fact_phyto_ChlNratio<0.7_rsh) fact_phyto_ChlNratio=0.7_rsh

          ELSE
              fact_phyto_ChlNratio=p_phyto_ChlNratio
          ENDIF
          diag_3d_wat(irk_diag(id_totalchl),k,i,j)=(c(iv_phyto_diat_N)+c(iv_phyto_dino_N)+c(iv_phyto_nano_N)) &
                                                *fact_phyto_ChlNratio


#if defined key_psnz
   ! Evolution des variables liees a Pseudo-nitzschia
   ! corrige et calcule les dc , evalue total chloro, variables diagnostiques et vitesses de chute pseudo_nitzschia
   ! ------------------------------------------------
   
#include "incellwat_bloom_pseudonitzschia.h"

#endif
#if defined key_karenia
   ! Evolution des variables liees a Karenia
   ! corrige et calcule les dc , evalue total chloro, variables diagnostiques et vitesses de chute karenia
   ! ---------------------------------------

#include "incellwat_bloom_karenia.h"

#endif
#if defined key_phaeocystis
   ! Evolution des variables liees a Phaeocystis
   ! corrige et calcule les dc , evalue total chloro, variables diagnostiques 
   ! -------------------------------------------

#include "incellwat_bloom_phaeocystis.h"

#endif

!   vitesses de chute des diatomees et du detritique enregistrees a la fin de la subroutine
#if ! defined key_BLOOM_opt2
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! modification des vitesses de chute de chaque variable bio si option 1 
   ! ws3 des diatomees seront modifiees si l_SNeffect_settle=true (effet des sels nuts sur chute des diatomees)
   !                                    si key_physadaptation (modulation par la profondeur)
   ! ws3 du materiel detritique seront modifiees si l_phyzoodeteffect_settle=true (modulation en fonction de l origine)
   !                                    si key_physadaptation (modulation par la profondeur)        

#include "incellwat_bloom_settling.h"

#endif        

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! enregistrement de la variable diagnostique : vitesse de chute
   !   vitesses de chute des diatomees et du detritique
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          diag_3d_wat(irk_diag(id_diatsettling),k,i,j)=ws3(k,iv_phyto_diat_N,i,j)*86400.0_rsh
          diag_3d_wat(irk_diag(id_detsettling),k,i,j)=ws3(k,iv_detr_N,i,j)*86400.0_rsh

#if defined key_ulvas
   ! Evolution des ulves  	ATTENTION arguments a lister
   ! ATTENTION NON inclus DANS LA NOUVELLE VERSION
   ! -------------------
          !CALL incellwat_ulvas(k,i,j,xe,c,dc)
#endif
#if defined key_zostera
   ! Evolution des variables liees aux herbiers de zosteres  atetention a besoin de pi
   ! corrige et calcule les dc , evalue total chloro, variables diagnostiques 
   ! ------------------------------------------------------

#include "incellwat_bloom_zostera.h"

#endif

#if defined key_benthos
   ! Evolution des 3 variables N,Si,P dans le sediment
   ! ------------------------------------------------
#include "incellwat_bloom_benthos.h"

        !  IF (ibbenth.eq.1) call incellwat_benthos(k,i,j,xe,c,dc)
#endif
#if defined key_suspensivores
   ! Evolution de la variable generique  suspensivore
   ! ------------------------------------------------
         ! CALL incellwat_suspensivores(k,i,j,xe,c,dc)
#endif
#if defined key_benthos_gener
   ! Evolution de 4 variables benthiques generiques : Bacterie+Meifaune, Deposivore-Herbivore, suspensivore, Carnivore
   ! ATTENTION NON OPERATIONNEL DANS LA NOUVELLE VERSION
   ! -------------------------------------------------------
          !CALL incellwat_benthos_generique(k,i,j,xe,c,dc)
#endif

#if defined key_microtracers
   ! Evolution des variables liees aux micro-traceurs
   ! ATTENTION NON inclus DANS LA NOUVELLE VERSION
   ! ------------------------------------------------
          !CALL incellwat_microtracers(k,i,j,xe,c,dc)
#endif
#if defined key_oxygen
   ! Evolution de la variable oxygene
   ! --------------------------------
#include "incellwat_bloom_oxygen.h"
 
#endif
#if defined key_oyster_SFG
   ! Evolution des huitres selon le modele SFG (Scope for growth) de Barille (1997)
   ! -------------------------------------------------------------------------------
#include "incellwat_bloom_oysterSFG.h"

#endif
#if defined key_oyster_DEB
   ! Evolution des huitres selon le modele DEB
   ! -------------------------------------------------------
#include "incellwat_bloom_oysterDEB.h"
   
#endif

#if defined key_N_tracer
   !   Tracage de l azote 
   !   --------------------
#include "incellwat_bloom_nitrogentracer.h"
 
#endif

#if defined key_P_tracer
   !   Tracage du phosphore 
   !     ATTENTION probleme possible pour karenia
   !   --------------------
#include "incellwat_bloom_phosphoretracer.h"

#endif

#if defined key_BIOLink_verif_conserv || defined key_MANGAbio
          !**********************************************
          ! Verification of conservativity at one  point
          !**********************************************
#ifdef key_BIOLink_verif_conserv
          IF ((i .eq. i_BIOLink_verif) .and. (j .eq. j_BIOLink_verif) .AND. &
             mod(CURRENT_TIME,dt_conserv_BIO) <=TRANSPORT_TIME_STEP .and. (k == 1 .or. k == kkmax)) THEN
#elif  key_MANGAbio
          IF((i == 70 .and. j == 180 .and. (k == 1 .or. k == kkmax)) .or. &
            (i == 145 .and. j == 175 .and. k == kkmax) .AND. &
             mod(CURRENT_TIME,dt_conserv_BIO) <=TRANSPORT_TIME_STEP) THEN
#endif            
           !bilan dc et c pour un point de grille
              CALL verif_bilan_NP(k,i,j,dc,c)
          ENDIF
#endif
 
! fin du test sur la salinite (si trop faible, pas de reaction bloomgique)
         ENDIF

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
                          *MAX(0.0_rsh,1.0_rsh-(tenfon(i,j)/tocdpe)) &
                          *tocd(ivp)/tocdpe
            depo=flx_w2s(ivp,i,j)*cw_bottom_MUSTANG(ivp,i,j)
            IF(irkm_var_assoc(ivp) > 0 ) THEN
              IF (flx_w2s(irkm_var_assoc(ivp),i,j)==0.0_rsh) flx_w2s(ivp,i,j)=0.0_rsh
            ELSE
              IF(depo .LE. epsdep_MUSTANG)flx_w2s(ivp,i,j)=0.0_rsh
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

#if defined key_benthos
        somdcN=somdcN+(dc(iv_benth_N)+p_burial*c(iv_benth_N))/thicklayerW_C(k,i,j)   
#endif
#if defined key_psnz
        somdcN=somdcN+dc(iv_phyto_psnz_N)   
#endif
#if defined key_karenia
        somdcN=somdcN+dc(iv_phyto_karenia_N)
#endif
#if defined key_phaeocystis
        somdcN=somdcN+dc(iv_phyto_phaeocystis_colo_N)+dc(iv_phyto_phaeocystis_cell_N) 
#endif
#if defined key_ulvas
        somdcN=somdcN+dc(iv_ulv_sus_N)+ dc(iv_ulv_benth_N)/thicklayerW_C(k,i,j)/0.014  
#endif
#endif

        IF(abs(somdcN) > 1.0e-3) THEN
          write(*,*)'Probleme conservativite N - Voir fichier conserv_BIOLink...'
          write(iscreenlog_conserv,*)'____________________________________________________'
          write(iscreenlog_conserv,*)'!!!! PB conservativite N: somdcN =',somdcN
          write(iscreenlog_conserv,*)'!!!! maille i, j, k:',i,j,k
          write(iscreenlog_conserv,*)'dc: NH4,NO3,diat,dino,nano,detr ',dc(iv_nutr_NH4),dc(iv_nutr_NO3),dc(iv_phyto_diat_N), &
                                                       dc(iv_phyto_dino_N),dc(iv_phyto_nano_N),dc(iv_detr_N)
          write(iscreenlog_conserv,*)'dc: mesozoo,microzoo',dc(iv_zoo_meso_N),dc(iv_zoo_micr_N)
#if defined key_benthos
          write(iscreenlog_conserv,*)'dc: benthos=   ',(dc(iv_benth_N)+p_burial*c(iv_benth_N))/thicklayerW_C(k,i,j)   
#endif
#if defined key_psnz
          write(iscreenlog_conserv,*)'dc: psnz=   ',dc(iv_phyto_psnz_N)   
#endif
#if defined key_karenia
          write(iscreenlog_conserv,*)'dc: karenia=   ',dc(iv_phyto_karenia_N)
#endif
#if defined key_phaeocystis
          write(iscreenlog_conserv,*)'dc: phaeocystis_colo, phaeocystis_cell=   ',dc(iv_phyto_phaeocystis_colo_N),dc(iv_phyto_phaeocystis_cell_N)
#endif
#if defined key_ulvas
          write(iscreenlog_conserv,*)'dc: suspended ulvas, benthic ulvas=   ',dc(iv_ulv_sus_N),dc(iv_ulv_benth_N)/thicklayerW_C(k,i,j)/0.014   
#endif
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

#if defined key_benthos
        somdcP=somdcP+(dc(iv_benth_P)+p_burial*c(iv_benth_P))/thicklayerW_C(k,i,j)   
#endif
#if defined key_psnz
        somdcP=somdcP+dc(iv_phyto_psnz_N)*rappaz   
#endif
#if defined key_karenia
!        somdcP=somdcP+dc(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)
!        somdcP=somdcP+dc(iv_phyto_karenia_N)*rappaz
        somdcP=somdcP+dc(iv_phyto_karenia_P)
#endif
#if defined key_phaeocystis
        somdcP=somdcP+(dc(iv_phyto_phaeocystis_colo_N)+dc(iv_phyto_phaeocystis_cell_N))*rappaz 
#endif
#if defined key_ulvas
        somdcP=somdcP+dc(iv_ulv_sus_P)+ dc(iv_ulv_benth_P)/thicklayerW_C(k,i,j)/0.031  
#endif
#endif

        IF(abs(somdcP) > 1.0e-3) THEN
          write(*,*)'Probleme conservativite P - Voir fichier conserv_BIOLink...'
          write(iscreenlog_conserv,*)'____________________________________________________'
          write(iscreenlog_conserv,*)'!!!! PB conservativite P: somdcP =',somdcP
          write(iscreenlog_conserv,*)'!!!! maille i, j, k:',i,j,k
          write(iscreenlog_conserv,*)'dc: PO4,Pads,diat,dino,nano,detr=   ',dc(iv_nutr_PO4),dc(iv_nutr_Pads), &
                        dc(iv_phyto_diat_N)*rappaz,dc(iv_phyto_dino_N)*rappaz,dc(iv_phyto_nano_N)*rappaz,dc(iv_detr_P)
          write(iscreenlog_conserv,*)'dc: mesozoo,microzoo=  ',dc(iv_zoo_meso_N)*rappaz,dc(iv_zoo_micr_N)*rappaz
#if defined key_benthos
          write(iscreenlog_conserv,*)'dc: benthos=   ',(dc(iv_benth_P)+p_burial*c(iv_benth_P))/thicklayerW_C(k,i,j)   
#endif
#if defined key_psnz
          write(iscreenlog_conserv,*)'dc: psnz=   ',dc(iv_phyto_psnz_N)*rappaz   
#endif
#if defined key_karenia
!          write(iscreenlog_conserv,*)'dc: karenia=   ',dc(iv_phyto_karenia_C)/(p_phyto_NPratio*p_phyto_CNratio)
!          write(iscreenlog_conserv,*)'dc: karenia=   ',dc(iv_phyto_karenia_N)*rappaz
          write(iscreenlog_conserv,*)'dc: karenia=   ',dc(iv_phyto_karenia_P) 
#endif
#if defined key_phaeocystis
          write(iscreenlog_conserv,*)'dc: phaeocystis_colo, phaeocystis_cell=   ',  &
                       dc(iv_phyto_phaeocystis_colo_N)*rappaz,dc(iv_phyto_phaeocystis_cell_N)*rappaz 
#endif
#if defined key_ulvas
          write(iscreenlog_conserv,*)'dc: suspended ulvas, benthic ulvas=   ',dc(iv_ulv_sus_P),dc(iv_ulv_benth_P)/thicklayerW_C(k,i,j)/0.031   
#endif
         ! write(iscreenlog_conserv,*) 'assimilmesozoo=',assimilmesozoo,'broumesozoodiat=',broumesozoodiat,'rationmesozoo=',rationmesozoo  
         ! write(iscreenlog_conserv,*) 'broumesozoodino=',broumesozoodino,'broumesozoomicrozoo=',broumesozoomicrozoo 
!          stop
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

#ifdef key_psnz
                            ! variables Pseudo-nitzschia
                            ! --------------------------
   REAL(KIND=rsh)           :: xmaxpsnz,datemaxpsnz,cumulpsnz
   INTEGER                  :: kmaxpsnz
#endif
#ifdef key_karenia
                            ! variables Karenia
                            ! -----------------
   REAL(KIND=rsh)           :: xmaxkarenia,datemaxkarenia,cumulkarenia
   INTEGER                  :: kmaxkarenia
#endif
#ifdef key_phaeocystis
                            ! variables Phaeocystis
                            ! ---------------------
   REAL(KIND=rsh)           :: xmaxphaeocystis,datemaxphaeocystis,cumulphaeocystis   
   INTEGER                  :: kmaxphaeocystis
#endif
#ifdef key_zoo_prod
                            ! variables production secondaire zooplanctonique
                            ! -----------------------------------------------
   REAL(KIND=rsh)           :: cumulpszoo_micr,cumulpszoo_meso   
#endif


!$OMP DO SCHEDULE(RUNTIME)
   DO j=jfirst,jlast
   

#ifdef key_MARS
    DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
    DO i=ifirst,ilast
#endif

#if ! defined key_BLOOM_opt2
       !variable diagnostique des max et des dat de max remises a zero tous les mois dans l option 1
       IF ((imois_BIOLINK.eq.1).and.(ijour_BIOLINK.eq.1).and.(iheure_BIOLINK.eq.0)) THEN
            diag_2d(irk_diag(id_diat_max),i,j)=0.0_rsh
            diag_2d(irk_diag(id_dino_max),i,j)=0.0_rsh
            diag_2d(irk_diag(id_nano_max),i,j)=0.0_rsh
#ifdef key_psnz
            diag_2d(irk_diag(id_psnz_max),i,j)=0.0_rsh
#endif
#ifdef key_karenia
            diag_2d(irk_diag(id_karenia_max),i,j)=0.0_rsh
#endif
#ifdef key_phaeocystis
            diag_2d(irk_diag(id_phaeocystis_max),i,j)=0.0_rsh
#endif

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

#ifdef key_psnz
   ! PseudoNitzschia
   ! -------------
         xmaxpsnz=0.0_rsh
         kmaxpsnz=1
         datemaxpsnz=1
         DO k=1,NB_LAYER_WAT
           IF (cvadv_wat_pos(iv_phyto_psnz_N,k,i,j).GE.xmaxpsnz) THEN
             xmaxpsnz=cvadv_wat_pos(iv_phyto_psnz_N,k,i,j)
             kmaxpsnz=k
             datemaxpsnz=jjulien_BIOLINK
           END IF
         END DO

         IF (xmaxpsnz.GT.diag_2d(irk_diag(id_psnz_max),i,j)) THEN
           diag_2d(irk_diag(id_psnz_max),i,j)=xmaxpsnz
           diag_2d(irk_diag(id_psnz_datemax),i,j)=datemaxpsnz
           diag_2d(irk_diag(id_psnz_depthmax),i,j)=- (thicklayerW_W(kmaxpsnz,i,j)-WATER_ELEV_ij)
         END IF
#endif
#ifdef key_karenia
   ! Karenia
   ! -------
         xmaxkarenia=0.0_rsh
         kmaxkarenia=1
         datemaxkarenia=1
         DO k=1,NB_LAYER_WAT
           IF (cvadv_wat_pos(iv_phyto_karenia_N,k,i,j).GE.xmaxkarenia) THEN
             xmaxkarenia=cvadv_wat_pos(iv_phyto_karenia_N,k,i,j)
             kmaxkarenia=k
             datemaxkarenia=jjulien_BIOLINK
           END IF
         END DO

         IF (xmaxkarenia.GT.diag_2d(irk_diag(id_karenia_max),i,j)) THEN
           diag_2d(irk_diag(id_karenia_max),i,j)=xmaxkarenia
           diag_2d(irk_diag(id_karenia_datemax),i,j)=datemaxkarenia
           diag_2d(irk_diag(id_karenia_depthmax),i,j)=- (thicklayerW_W(kmaxkarenia,i,j)-WATER_ELEV_ij)
         END IF
#endif
#ifdef key_phaeocystis
   ! Phaeocystis
   ! -----------
         xmaxphaeocystis=0.0_rsh
         kmaxphaeocystis=1
         datemaxphaeocystis=1
         DO k=1,NB_LAYER_WAT
           IF (cvadv_wat_pos(iv_phyto_phaeocystis_colo_N,k,i,j).GE.xmaxphaeocystis) THEN
             xmaxphaeocystis=cvadv_wat_pos(iv_phyto_phaeocystis_colo_N,k,i,j)
             kmaxphaeocystis=k
             datemaxphaeocystis=jjulien_BIOLINK
           END IF
         END DO

         IF (xmaxphaeocystis.GT.diag_2d(irk_diag(id_phaeocystis_max),i,j)) THEN
           diag_2d(irk_diag(id_phaeocystis_max),i,j)=xmaxphaeocystis
           diag_2d(irk_diag(id_phaeocystis_datemax),i,j)=datemaxphaeocystis
           diag_2d(irk_diag(id_phaeocystis_depthmax),i,j)=- (thicklayerW_W(kmaxphaeocystis,i,j)-WATER_ELEV_ij)
         END IF
#endif

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
   !   !!!! ??? est ce que ca sert vraiment ???
   ! -------------------------------------------------------------
#ifdef key_MARS
         IF (ljmin==jmin .AND. j.EQ.jb(i)) obc_cv_ktherm_s(i)=kdcdz
         IF (limin==imin .AND. i.EQ.ig(j)) obc_cv_ktherm_w(j)=kdcdz
#endif

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
#ifdef key_psnz
         cumulpsnz=0.0_rsh
#endif
#ifdef key_karenia
         cumulkarenia=0.0_rsh
#endif
#ifdef key_phaeocystis
         cumulphaeocystis=0.0_rsh
#endif
#ifdef key_zoo_prod
         cumulpszoo_micr=0.0_rsh
         cumulpszoo_meso=0.0_rsh
#endif
         DO k=1,NB_LAYER_WAT
           cumuldiat=cumuldiat+cvfix_wat_pos(iv_phyto_diat_pp-nv_adv,k,i,j)
           cumuldino=cumuldino+cvfix_wat_pos(iv_phyto_dino_pp-nv_adv,k,i,j)
           cumulnano=cumulnano+cvfix_wat_pos(iv_phyto_nano_pp-nv_adv,k,i,j)
#ifdef key_psnz
           cumulpsnz=cumulpsnz+cvfix_wat_pos(iv_phyto_psnz_pp-nv_adv,k,i,j)
#endif
#ifdef key_karenia
           cumulkarenia=cumulkarenia+cvfix_wat_pos(iv_phyto_karenia_pp-nv_adv,k,i,j)
#endif
#ifdef key_phaeocystis
           cumulphaeocystis=cumulphaeocystis+cvfix_wat_pos(iv_phyto_phaeocystis_pp-nv_adv,k,i,j)
#endif
#ifdef key_zoo_prod
           cumulpszoo_micr=cumulpszoo_micr+cvfix_wat_pos(iv_zoo_micr_ps-nv_adv,k,i,j)
           cumulpszoo_meso=cumulpszoo_meso+cvfix_wat_pos(iv_zoo_meso_ps-nv_adv,k,i,j)
#endif
         END DO
         diag_2d(irk_diag(id_diat_columnprod),i,j)=cumuldiat
         diag_2d(irk_diag(id_dino_columnprod),i,j)=cumuldino
         diag_2d(irk_diag(id_nano_columnprod),i,j)=cumulnano
#ifdef key_psnz
         diag_2d(irk_diag(id_psnz_columnprod),i,j)=cumulpsnz
#endif
#ifdef key_karenia
         diag_2d(irk_diag(id_karenia_columnprod),i,j)=cumulkarenia
#endif
#ifdef key_phaeocystis
         diag_2d(irk_diag(id_phaeocystis_columnprod),i,j)=cumulphaeocystis
#endif
#ifdef key_zoo_prod
         diag_2d(irk_diag(id_zoo_micr_columnprod),i,j)=cumulpszoo_micr
         diag_2d(irk_diag(id_zoo_meso_columnprod),i,j)=cumulpszoo_meso
#endif
#if ! defined key_BLOOM_opt2
   !   production cumulee dans le temps et sur la verticale de tous les types phytos   
         diag_2d(irk_diag(id_columnprodtotal),i,j)=cumuldiat+cumuldino+cumulnano
#ifdef key_psnz
         diag_2d(irk_diag(id_columnprodtotal),i,j)=diag_2d(irk_diag(id_columnprodtotal),i,j) + cumulpsnz
#endif 
#ifdef key_karenia
         diag_2d(irk_diag(id_columnprodtotal),i,j)=diag_2d(irk_diag(id_columnprodtotal),i,j) + cumulkarenia
#endif  
#ifdef key_phaeocystis
         diag_2d(irk_diag(id_columnprodtotal),i,j)=diag_2d(irk_diag(id_columnprodtotal),i,j) + cumulphaeocystis
#endif 
#ifdef key_zoo_prod
         diag_2d(irk_diag(id_zoo_columnprodtotal),i,j)=cumulpszoo_micr+cumulpszoo_meso
#endif                        
#endif                                             
           
     END IF  ! H>RESIDUAL_THICKNESS_WAT
   END DO
   END DO
!$OMP END DO
      
  END SUBROUTINE bloom_eval_diag2d
!!============================================================================================
#if defined key_BLOOM_insed && defined key_oxygen && ! defined key_BLOOM_opt2
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
   REAL(KIND=rsh)         :: porosite_inv,dtbiojour,txfiltbenth,txfiltbenth_max,foncsinus, &
                              cvolp,txfiltbenth_wat,epn_bott,temper,effetchaleur, &
                              factlim,filtr_benth,diatmortsed,xtmp,ztop, &
                              zmiddle,effetchaleurSi
   REAL(KIND=rsh)         :: flim1_O2,flim2_O2,flim1_NO3,glim1_NO3,Sfliminv, SfliminvR, &
            glim1_O2,F_remin_aerN,F_remin_anaerN,F_reminR_aerN,F_reminR_anaerN, flimz,&
            F_reminR_NO3,F_denitR,F_dnraR,F_remin_NO3,F_denit,F_dnra,F_nitrif, F_burried, &
            F_oxyd_ODU,F_solid_ODU,flim3_O2,flim4_O2,flim2_NO3,glim2_O2,glim3_O2,glim2_NO3, &
            F_remin_aerP,F_remin_anaerP,F_reminR_aerP,F_reminR_anaerP,  &
            F_remin_NO3_P,F_reminR_NO3_P,adsormax,F_adsorP,F_desorP,F_transf_MO,  &
            F_precPFE_O2,F_precPFE_NO3,F_dissolPFE,  &
!            glim4_O2,F_remin_aerS,F_remin_anaerS,F_precSi, &
            glim_silica,F_remin_Si,F_precSi,     &
            cvO2new,cvO2old,Ftot_conso2_remin,Ftot_conso2_reminR,reduc_aer, &
            F_aeration,KO2_aeration,KzO2_aeration,o2sats,tempabs,gst
   REAL(KIND=rsh),DIMENSION(nv_state)  :: cs,cw_bott,dcdt,dcw_filtbent

   !!----------------------------------------------------------------------
   !! * Executable part
   
   dtbiojour=dt_true/86400._rsh

#ifdef key_MARS
   IF(l_filtbenthsinus) THEN
     ! Introduction d une pression de Broutage, d apres Marie SAVINA
     ! -------------------------------------------------------------  
     ! txfiltbenth = taux sinusoidal de filtration du benthos en m3/j/m2
     foncsinus=(1.0_rsh+SIN(2.0_rsh*3.14159_rsh/365.0_rsh*(jjulien_BIOLink-125)))/2.0_rsh
     txfiltbenth_max=(p_txfiltbenthmax*(0.3_rsh+0.7_rsh*foncsinus))
   ELSE
     txfiltbenth_max=p_txfiltbenthmax
   ENDIF
#else
     txfiltbenth_max=p_txfiltbenthmax
#endif


!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,iv,ksmin,ksmax,cs,cw_bott,dcdt,dcw_filtbent, &
!$OMP txfiltbenth_wat,cvolp,porosite_inv,foncsinus,zmiddle,txfiltbenth, flimz, &
!$OMP epn_bott,temper,effetchaleur,factlim,filtr_benth,diatmortsed,xtmp,ztop,  &
!$OMP flim1_O2,flim2_O2,flim1_NO3,glim1_NO3,Sfliminv,glim2_O2,effetchaleurSi, &
!$OMP glim1_O2,F_remin_aerN,F_remin_anaerN,F_reminR_aerN,F_reminR_anaerN, &
!$OMP F_reminR_NO3,F_denitR,F_dnraR,F_remin_NO3,F_denit,F_dnra,F_nitrif,F_burried, &
!$OMP F_oxyd_ODU,F_solid_ODU,flim3_O2,flim4_O2,flim2_NO3,glim3_O2,glim2_NO3, &
!$OMP F_remin_aerP,F_remin_anaerP,F_reminR_aerP,F_reminR_anaerP,SfliminvR,  &
!$OMP F_remin_NO3_P,F_reminR_NO3_P,adsormax,F_adsorP,F_desorP,F_transf_MO,  &
!$OMP F_precPFE_O2,F_precPFE_NO3,F_dissolPFE,F_precSi,cvO2new,cvO2old,      &
!$OMP Ftot_conso2_remin,Ftot_conso2_reminR,reduc_aer,glim_silica,F_remin_Si, &
!$OMP F_aeration,KO2_aeration,KzO2_aeration,o2sats,tempabs,gst)
     DO j=jfirst,jlast
#ifdef key_MARS
      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1) 
#else
      DO i=ifirst,ilast
#endif

       IF (BATHY_H0(i,j) > -valmanq .AND. ksma(i,j) > 0)  THEN  ! on n est pas a terre et il y a du sediment
       
        ksmin=ksmi(i,j)
        ksmax=ksma(i,j)
        
     ! memorisation des concentrations locales dans l eau surnageante en gardant valeurs seulement positives
        IF (htot(i,j)> RESIDUAL_THICKNESS_WAT) THEN  ! test sur hauteur d eau
#ifdef key_MARS
          IF (htot(i,j) < hm ) THEN  
          ! cas  : faible profondeur < hm : en 2D in MARS
            epn_bott=htot(i,j)
            DO iv=1,nv_state
              cw_bott(iv)=SUM(cv_wat(iv,:,i,j)*dsigu(:))
              cw_bott(iv)=0.5_rsh*(cw_bott(iv)+ABS(cw_bott(iv)))
            ENDDO
          
          ELSE IF (htot(i,j) > hm) THEN
          !  profondeur > hm
            epn_bott=epn_bottom_MUSTANG(i,j)
#else
            epn_bott=epn_bottom_MUSTANG(i,j)
#endif 
            cw_bott(1:nv_state)=0.5_rsh*( cw_bottom_MUSTANG(1:nv_state,i,j)+   &
                                ABS(cw_bottom_MUSTANG(1:nv_state,i,j)) ) 
#ifdef key_MARS
          ENDIF
#endif
          o2sats=0.0_rsh
          KO2_aeration=0.0_rsh
          KzO2_aeration=0.0_rsh
                        
        ELSE
          !   il n y a pas d eau
          cw_bott(1:nv_state)=0.0_rsh
          
          ! si y a pas d eau reaeration par l air pour evolution de l oxygen en surface
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          tempabs=cv_sed(-1,ksmax,i,j)+273.15_rsh
          gst=-173.4292_rsh+249.6329_rsh*100.0_rsh/tempabs+143.3483_rsh*log(tempabs/100.0_rsh)    &
           -21.8492_rsh*tempabs/100.0_rsh+cv_sed(0,ksmax,i,j)*(-0.033096_rsh+0.014259_rsh*tempabs/100.0_rsh &
           -0.0017_rsh*(tempabs/100.0_rsh)**2)
          o2sats=1.429_rsh*exp(gst)
          KO2_aeration=p_KO2sed_aeration
          KzO2_aeration=p_Kzsed_aeration
        ENDIF
         
        ! conversion des concentrations dans l eau si changement d unite 
        ! si matiere organique supposee constitutive du sediment
        cw_bott(1:nv_state)= cw_bott(1:nv_state)/unit_modif_mudbio_N2dw(irk_fil(1:nv_state))                            
         

        ! calcul des flux reactifs pour chaque variable et chaque processus dans chaque couche
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ztop=-dzs(ksmax,i,j)
        zmiddle=0.0_rsh
        DO k=ksmax,ksmin,-1
          ! no bio transformation if very small layer
          IF (dzs(k,i,j) > dzsmin) THEN
             ! transfer 3D var cv_sed into 1D var c
             cs(1:nv_state)=0.5_rsh*( cv_sed(1:nv_state,k,i,j)+   &
                                  ABS(cv_sed(1:nv_state,k,i,j)) )
             cs(1:nv_state)= cs(1:nv_state)/unit_modif_mudbio_N2dw(irk_fil(1:nv_state))
             porosite_inv=1._rsh/poro(k,i,j) 
             diag_3d_sed(id_porosite_sed,k,i,j)=poro(k,i,j)
             ztop=ztop+dzs(k,i,j)
             if(k==ksmax) then
               zmiddle=dzs(k,i,j)/2
             else
               zmiddle=zmiddle+dzs(k+1,i,j)/2._rsh+dzs(k,i,j)/2._rsh
             endif
             ! initialize to 0
             !fluxbio_tot(1:nv_state,k)=0.0_rsh
             dcdt(:)=0.0_rsh
             temper=cv_sed(-1,k,i,j)
             !IF (temper.gt.30.) THEN
             !write(*,*)'temp sed>30:',k,i,j,temper
             !    temper=15.
             !END IF
             effetchaleur=exp(p_T_effect*temper)
             
             ! Pression de broutage (filtrage / broutage ) (Marina Savina)
             ! -----------------------------------------------------------
             filtr_benth=0.0_rsh
             txfiltbenth=0.0_rsh
             txfiltbenth_wat=0.0_rsh
             dcw_filtbent(:)=0.0_rsh
             IF(k == ksma(i,j) .AND. htot(i,j) > RESIDUAL_THICKNESS_WAT .AND. txfiltbenth_max .NE. 0.0_rsh) THEN

                 ! and txfiltbenth fonction des MES (filtration baisse avec la turbidité)
                 txfiltbenth=MAX(-0.0333_rsh*cmes_3dmgl(1,i,j)+1,0.0_rsh)*txfiltbenth_max  !de 0 mg/L a 30 mg/L
                 ! conversion par metre3 eau
                 txfiltbenth_wat=txfiltbenth/epn_bott
                 ! resolution explicite pour cw_bott ==> limitation pas de temps
                 factlim=MIN(1.0_rsh,1.0_rsh/dtbiojour/(txfiltbenth_wat+epsilon_BIOLink))
                 IF(factlim .NE. 1.0_rsh) write(*,*)'filt',(CURRENT_TIME-TIME_BEGIN)/3600._rsh, &
                                        i,j,k,epn_bott,factlim,txfiltbenth_wat,dtbiojour
                 txfiltbenth_wat=txfiltbenth_wat*factlim
                 txfiltbenth=txfiltbenth*factlim/dzs(k,i,j)
                 filtr_benth= txfiltbenth*(cw_bott(iv_phyto_diat_N)+cw_bott(iv_phyto_dino_N)+cw_bott(iv_phyto_nano_N))
                 ! filtrage sur les concentrations dans l eau
                 dcw_filtbent(iv_phyto_diat_N)=-txfiltbenth_wat*cw_bott(iv_phyto_diat_N)
                 dcw_filtbent(iv_phyto_dino_N)=-txfiltbenth_wat*cw_bott(iv_phyto_dino_N)
                 dcw_filtbent(iv_phyto_nano_N)=-txfiltbenth_wat*cw_bott(iv_phyto_nano_N)
#if defined key_psnz
                 dcw_filtbent(iv_phyto_psnz_N)=-txfiltbenth_wat*cw_bott(iv_phyto_psnz_N)
                 dcw_filtbent(iv_phyto_psnz_Si)=-txfiltbenth_wat*cw_bott(iv_phyto_psnz_Si)
#endif
#if defined key_karenia
                 dcw_filtbent(iv_phyto_karenia_N)=-txfiltbenth_wat*cw_bott(iv_phyto_karenia_N)
                 dcw_filtbent(iv_phyto_karenia_C)=-txfiltbenth_wat*cw_bott(iv_phyto_karenia_C)
#endif
#ifdef key_phaeocystis
                 dcw_filtbent(iv_phyto_phaeocystis_cell_N)=-txfiltbenth_wat*cw_bott(iv_phyto_phaeocystis_cell_N)
#endif
#ifdef key_ulvas
                 dcw_filtbent(iv_ulv_benth_N)=-ulvebenthmortsed*cs(iv_ulv_benth_N)/0.014_rsh
                 dcw_filtbent(iv_ulv_benth_P)=-ulvebenthmortsed*cs(iv_ulv_benth_P)/0.014_rsh
#endif
             ENDIF
            
             ! Minéralisations
             ! -------------------
             ! cycle de l azote
               ! fonctions limitantes
               cvO2old=cs(iv_oxygen)
               IF(cs(iv_oxygen).gt.0.0001_rsh) THEN
                 flim1_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_reminO2)
                 flim2_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_nit)
                 flim3_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_reoxyd)
               ELSE
                 flim1_O2 = 0.0_rsh
                 flim2_O2 = 0.0_rsh
                 flim3_O2 = 0.0_rsh
               ENDIF
               IF(cs(iv_nutr_NO3).gt.0.001_rsh) THEN
                 flim1_NO3 = cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kNO3_reminssO2)
               ELSE
                 flim1_NO3 = 0.0_rsh
               ENDIF
               glim1_O2 = 1._rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_remin0O2)
               glim2_O2 = 1._rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_denit)
               glim1_NO3 = 1._rsh - cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kiNO3_remin0O2)
               Sfliminv = 1._rsh / (flim1_O2 + flim1_NO3 * glim2_O2 + glim1_NO3 * glim1_O2)
               !IF(flim1_O2 > 0._rsh .OR. flim1_NO3 > 0._rsh) THEN 
               !  SfliminvR = 1._rsh / (flim1_O2 + flim1_NO3 * glim2_O2)
               !ELSE
               !  SfliminvR = 0.0_rsh                ! Pas de remineralisation anaerobie de la MO refractaire
               !ENDIF
               ! aerobiose
               ! flux deja multiplie par dtbiojour pour le calcul ensuite des flux en var diag
               xtmp = effetchaleur * flim1_O2 * dtbiojour
               flimz = max(p_xflimz,exp(-p_k_remin*zmiddle))
               F_remin_aerN = p_N_remin * xtmp * Sfliminv
               !F_reminR_aerN = p_N_reminR * xtmp * SfliminvR
               F_reminR_aerN = p_N_reminR * flimz * xtmp * Sfliminv
               F_nitrif = p_nitrif * effetchaleur * flim2_O2 * dtbiojour
               !F_nitrif = p_nitrif * flim2_O2 * dtbiojour
  
               ! sub-oxique (NO3 comme oxydant)
               xtmp = effetchaleur * flim1_NO3 * glim2_O2 * dtbiojour
               F_remin_NO3 = p_N_remin * xtmp * Sfliminv
               !F_reminR_NO3 = p_N_reminR * xtmp * SfliminvR
               F_reminR_NO3 = p_N_reminR * flimz * xtmp * Sfliminv
               
               F_denit = p_DNO3_denit * F_remin_NO3
               F_denitR = p_DNO3_denit * F_reminR_NO3
               F_dnra = (1._rsh - p_DNO3_denit) * F_remin_NO3
               F_dnraR = (1._rsh - p_DNO3_denit) * F_reminR_NO3
               !F_dnra = 0.0_rsh 
               !F_dnraR = 0.0_rsh
               !F_denit = 0.0_rsh
               !F_denitR = 0.0_rsh
               
               ! ODU
               F_oxyd_ODU = p_ODU_oxy * flim3_O2 * dtbiojour
               F_solid_ODU = p_ODU_precip * dtbiojour

               ! anaerobiose
               xtmp = effetchaleur * glim1_O2 * glim1_NO3 * dtbiojour
               F_remin_anaerN = p_N_remin * xtmp * Sfliminv
               F_reminR_anaerN = p_N_reminR *flimz * xtmp * Sfliminv 
               !F_reminR_anaerN = 0.0_rsh
            
             ! cycle du phosphore
               ! fonctions limitantes
               IF(cs(iv_oxygen).gt.0.01_rsh) THEN
                 flim4_O2 = cs(iv_oxygen) / (cs(iv_oxygen) + p_kO2_precPFe)
               ELSE
                 flim4_O2 = 0.0_rsh
               ENDIF
               IF(cs(iv_nutr_NO3).gt.0.001_rsh) THEN
                 flim2_NO3 = cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kNO3_precPFe)
               ELSE
                 flim2_NO3 = 0.0_rsh
               ENDIF
               glim3_O2 = 1._rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_dissPFe)
               glim2_NO3 = 1._rsh - cs(iv_nutr_NO3) / (cs(iv_nutr_NO3) + p_kiNO3_dissPFe)
  
               ! aerobiose
               xtmp = effetchaleur * flim1_O2 * dtbiojour
               F_remin_aerP = p_P_remin * xtmp * Sfliminv
               !F_reminR_aerP = p_P_reminR * xtmp * SfliminvR
               F_reminR_aerP = p_P_reminR * flimz * xtmp * Sfliminv
               F_precPFE_O2 = p_P_precFeO2 * flim4_O2 * dtbiojour
 
               ! sub-oxique (NO3 comme oxydant)
               xtmp = effetchaleur * flim1_NO3 * glim2_O2 * dtbiojour
               F_remin_NO3_P = p_P_remin * xtmp * Sfliminv
               !F_reminR_NO3_P = p_P_reminR * xtmp * SfliminvR
               F_reminR_NO3_P = p_P_reminR * flimz * xtmp * Sfliminv
               F_precPFE_NO3 = p_P_precFeNO3 * glim3_O2 * flim2_NO3 * dtbiojour

               ! anaerobiose
               xtmp = effetchaleur * glim1_O2 * glim1_NO3 * dtbiojour
               F_remin_anaerP = p_P_speedup_reminanaer * p_P_remin * xtmp * Sfliminv
               F_reminR_anaerP = p_P_speedup_reminanaer * p_P_reminR * flimz * xtmp * Sfliminv 
               !F_reminR_anaerP = 0.0_rsh
               ! Benedicte : suppression de ce test anachronique et sans doute non generalisable
               !             remplace par la reduction de la vitesse de dissolution du PFe qui diminue avec la profondeur
               !             de la meme maniere que la mineralisation poue tenir compte du fit que le PFe devient de plus en plus
               !             non reactif au fur et a mesure que la matiere organique vieillit et protege le PFe en l entourant ou en lencapsulant
               ! Le Pfe n est alors plus disponible pour sa dissolution (voir V. Mesnage)
               !IF(cs(iv_PFe).gt.(3.0_rsh*c_sedtot(k,i,j))) THEN         !attention donner les seuils en mmol/m3 pour les vars part.
                 F_dissolPFE = p_P_dissFe * glim3_O2 * glim2_NO3 * flimz * dtbiojour
               !ELSE
               !  F_dissolPFE = 0.0_rsh                                  ! il reste toujours une base de PFe dans le sed profond
               !ENDIF            
 
               ! adsorption / desorption sur mes 
               !adsormax = p_P_adsormaxsed * (cs(iv_MES_local) + cs(iv_spim)) * dtbiojour
               adsormax = p_P_adsormaxsed * cs(iv_spim)             
               IF(cs(iv_oxygen).gt.0.2_rsh) THEN
                 F_adsorP = p_P_adsor * max(0.0_rsh, (adsormax - cs(iv_nutr_Pads))) * dtbiojour
               ELSE
                 F_adsorP = p_P_adsor/100 * max(0.0_rsh, (adsormax - cs(iv_nutr_Pads))) * dtbiojour
               ENDIF
               IF(adsormax.gt.0.0_rsh) THEN
                 F_desorP = p_P_desor * (cs(iv_nutr_Pads) / adsormax) * dtbiojour
               ELSE
                 F_desorP = 0.0_rsh
               ENDIF

             ! cycle de la silice (voir Khalil et al., 2007)
               ! fonctions limitantes
               !glim4_O2 = 1._rsh - cs(iv_oxygen) / (cs(iv_oxygen) + p_kiO2_precSi)  
                glim_silica = 1-cs(iv_nutr_SiOH)/p_Si_Eq
                if(glim_silica.lt.0._rsh) glim_silica=0._rsh
                effetchaleurSi = exp(p_T_effectSi*temper)
                F_remin_Si = ((p_BSi_dissSurfSed-p_BSi_dissFondSed) * exp(-p_kSi*zmiddle) + p_BSi_dissFondSed) * &
                             effetchaleurSi * glim_silica *dtbiojour
                if(cs(iv_nutr_SiOH).gt.p_Si_EqPrec) then
                  F_precSi = p_Si_precip * (cs(iv_nutr_SiOH)-p_Si_EqPrec) *dtbiojour
                else
                  F_precSi = 0._rsh
                endif 

               ! aerobiose
               !F_remin_aerS = p_Si_diss * effetchaleur * flim1_O2 * dtbiojour
               
               ! anaerobiose
               !F_remin_anaerS = p_Si_diss/100 * glim1_O2 * effetchaleur * dtbiojour
               !F_precSi = p_Si_precip * glim4_O2 * dtbiojour

                
             ! transfert de MO labile en MO refractaire par vieillissement
               F_transf_MO = p_aging_MO * dtbiojour
               IF(k==ksmin) THEN
                 F_burried = p_burried * dtbiojour
               ELSE
                 F_burried = 0.0_rsh
               ENDIF

             ! cycle de l oxygene
             ! si pas d eau ==> reaeration par l air en surface
     ! vilaine il y a toujours de l eau donc mis en commentaire
             !F_aeration=KO2_aeration*exp(-KzO2_aeration*ztop)*(o2sats-cs(iv_oxygen))*dtbiojour
             !diag_3d_sed(id_fluxsed_aeration,k,i,j)=diag_3d_sed(id_fluxsed_aeration,k,i,j)+F_aeration*dzs(k,i,j)
            
             ! mortalite des diatomees
             ! -----------------------
             diatmortsed=p_diat_mort_sed*effetchaleur*dtbiojour


             ! calcul des dc (*dt)==> variation M/V pendant dt
             ! ***********************************************-
            
                     !!!!!!!!!!!!!!!!!!
                     !cycle de l azote
                     !!!!!!!!!!!!!!!!!!
                         
!               F_oxyd_ODU = 0.0_rsh
!               F_remin_aerN = 0.0_rsh
!               F_reminR_aerN = 0.0_rsh
!               F_nitrif = 0.0_rsh
!               F_precPFE_O2 = 0.0_rsh 

             ! Dissolved oxygen (mg/l) 
             ! -------------------------
!             dcdt(iv_oxygen)= -F_oxyd_ODU *cs(iv_ODU)*porosite_inv                &
             dcdt(iv_oxygen)= (-F_oxyd_ODU * cs(iv_ODU)                                      &
                               -F_remin_aerN * cs(iv_detr_N) * porosite_inv * p_GO2_NOrg     &
                               -F_reminR_aerN * cs(iv_detrR_N) * porosite_inv * p_GO2_NOrgR  &
                               -F_nitrif * cs(iv_nutr_NH4) * p_GO2_NH4)                      &
                              * 0.032_rsh ! transformation en mg/L

!            ................................................................... 
            
             ! Detrital N (micromol/lsed) 
             ! -----------------------
             dcdt(iv_detr_N)=(-F_remin_aerN - F_remin_anaerN - F_remin_NO3 - F_transf_MO - F_burried) * cs(iv_detr_N)      &
                            +diatmortsed * cs(iv_phyto_diat_N)                    &
#if defined key_psnz
                            +txfiltbenth * cw_bott(iv_phyto_psnz_N)*dtbiojour                &
#endif
#if defined key_karenia
                            +txfiltbenth * cw_bott(iv_phyto_karenia_N)*dtbiojour             &
#endif
#ifdef key_phaeocystis
                            +txfiltbenth * cw_bott(iv_phyto_phaeocystis_cell_N) * dtbiojour    &
#endif
#ifdef key_ulvas
                            +ulvebenthmortsed * cs(iv_ulv_benth_N) / 0.014_rsh * dtbiojour       &
#endif

                            +filtr_benth * dtbiojour   ! (broutage et filtre par phytoplanctons de base)
                           
             dcdt(iv_detrR_N)=(-F_reminR_aerN - F_reminR_anaerN - F_reminR_NO3 - F_burried) * cs(iv_detrR_N)  &
                             +F_transf_MO * cs(iv_detr_N)
                      
             ! var diag : flux en M/m2 cumule dans le temps 
             diag_3d_sed(id_remin_aerN,k,i,j)=diag_3d_sed(id_remin_aerN,k,i,j)          &
                                             +(F_remin_aerN * cs(iv_detr_N) + F_reminR_aerN * cs(iv_detrR_N)) * dzs(k,i,j)
             diag_3d_sed(id_remin_anaerN,k,i,j)=diag_3d_sed(id_remin_anaerN,k,i,j)      &
                                               +(F_remin_anaerN * cs(iv_detr_N) + F_reminR_anaerN * cs(iv_detrR_N)) * dzs(k,i,j)
             diag_3d_sed(id_remin_nitrateN,k,i,j)=diag_3d_sed(id_remin_nitrateN,k,i,j)  &
                                                 +(F_remin_NO3 * cs(iv_detr_N) + F_reminR_NO3 * cs(iv_detrR_N)) * dzs(k,i,j)
             diag_3d_sed(id_filtr_benth,k,i,j)=diag_3d_sed(id_filtr_benth,k,i,j)        &
                                              + filtr_benth * dtbiojour * dzs(k,i,j)
!            ................................................................... 
             
             ! Ammonium (micromol/l ei) 
             ! ---------------------
             dcdt(iv_nutr_NH4)=(F_remin_aerN + F_remin_anaerN + F_dnra) * cs(iv_detr_N) * porosite_inv  &
                              +(F_reminR_aerN + F_reminR_anaerN + F_dnraR) * cs(iv_detrR_N) * porosite_inv  &
!#if defined key_zostera
!                  +(F_remin_aerN/10._rsh)*cs(iv_detr_N)*porosite_inv   &
!#endif
                              -F_nitrif * cs(iv_nutr_NH4)

             diag_3d_sed(id_nitrif,k,i,j)=diag_3d_sed(id_nitrif,k,i,j)  &
                                         +F_nitrif * cs(iv_nutr_NH4) * dzs(k,i,j)
             diag_3d_sed(id_remin_drnaN,k,i,j)=diag_3d_sed(id_remin_drnaN,k,i,j)  &
                                              +F_dnra * cs(iv_detr_N) * porosite_inv * dzs(k,i,j) &
                                              +F_dnraR * cs(iv_detrR_N) * porosite_inv * dzs(k,i,j)
             
             !IF(cs(iv_nutr_NH4)>60.0 .AND. k<40 .AND. k>30) THEN
             !  print*,'!!!!!! k=',k,'NH4=',cs(iv_nutr_NH4)
             !  print*,'F_remin_aerN=',F_remin_aerN,'F_reminR_aerN=',F_reminR_aerN
             !  print*,'effetchaleur',effetchaleur,'flim1_O2=',flim1_O2
             !  print*,'F_remin_anaerN=',F_remin_anaerN
             !  print*,'glim1_O2=',glim1_O2,'glim1_NO3=',glim1_NO3
             !  print*,'F_nitrif=',F_nitrif,'flim2_O2=',flim2_O2
             !ENDIF
!            ................................................................... 

             ! Nitrate (micromol/l ei) 
             ! --------------------
             dcdt(iv_nutr_NO3)=F_nitrif * cs(iv_nutr_NH4)                 &
                              -(F_denit + F_dnra) * cs(iv_detr_N) * p_GNO3_Norg * porosite_inv        &
                              -(F_denitR + F_dnraR) * cs(iv_detrR_N) * p_GNO3_NorgR * porosite_inv          
             diag_3d_sed(id_remin_denitN,k,i,j)=diag_3d_sed(id_remin_denitN,k,i,j)  &
                                               +F_denit * cs(iv_detr_N) * p_GNO3_Norg * porosite_inv * dzs(k,i,j)  &
                                               +F_denitR * cs(iv_detrR_N) * p_GNO3_NorgR * porosite_inv * dzs(k,i,j)
!            ................................................................... 

             ! Oxygen Demand Unit (micromol/l ei) 
             ! -------------------------------
             dcdt(iv_ODU)= F_remin_anaerN * cs(iv_detr_N) * p_GODU_NOrg * porosite_inv     &
                         + F_reminR_anaerN * cs(iv_detrR_N) * p_GODU_NOrgR *porosite_inv  &
                         -(F_oxyd_ODU+F_solid_ODU)*cs(iv_ODU)
 
             diag_3d_sed(id_oxyd_solid_ODU,k,i,j)=diag_3d_sed(id_oxyd_solid_ODU,k,i,j)  &
                                                 -(F_oxyd_ODU+F_solid_ODU)*cs(iv_ODU)*dzs(k,i,j)
!            ...................................................................                          

                    !!!!!!!!!!!!!!!!!!!
                    !cycle du phosphore
                    !!!!!!!!!!!!!!!!!!!
            
             ! Phosphate (micromol/l ei) 
             ! ----------------------
             dcdt(iv_nutr_PO4)=(F_remin_aerP + F_remin_anaerP + F_remin_NO3_P) * cs(iv_detr_P) * porosite_inv  &
                              +(F_reminR_aerP + F_reminR_anaerP + F_reminR_NO3_P) * cs(iv_detrR_P) * porosite_inv  &
                              +F_desorP * cs(iv_nutr_Pads) * porosite_inv               &
                              +F_dissolPFE * cs(iv_PFe) * porosite_inv                  &
!#if defined key_zostera
!                    +(F_remin_aerP/10._rsh)*cs(iv_detr_P)*porosite_inv*dtbiojour                &
!#endif
                              -(F_adsorP + F_precPFE_O2 + F_precPFE_NO3)*cs(iv_nutr_PO4)
!            ...................................................................                       
                            
             ! Detrital P (micromol/l sed) 
             ! -----------------------
             dcdt(iv_detr_P)=(-F_remin_aerP - F_remin_anaerP - F_remin_NO3_P - F_transf_MO - F_burried)  * cs(iv_detr_P)     &
                            +diatmortsed * cs(iv_phyto_diat_N) * rappaz                &
#if defined key_psnz
                            +txfiltbenth * cw_bott(iv_phyto_psnz_N) * rappaz * dtbiojour                        &
#endif
#if defined key_karenia
                            +txfiltbenth * cw_bott(iv_phyto_karenia_C) / (p_phyto_NPratio*p_phyto_CNratio) * dtbiojour  &
#endif
#ifdef key_phaeocystis
                            +txfiltbenth * cw_bott(iv_phyto_phaeocystis_cell_N) * rappaz * dtbiojour            &
#endif
#ifdef key_ulvas
                            +ulvebenthmortsed * cw_bott(iv_ulv_benth_P) / 0.031_rsh *dtbiojour                &
#endif 
                            +filtr_benth * rappaz * dtbiojour   ! (broutage et filtre par phytoplanctons de base)        
 
             dcdt(iv_detrR_P)=(-F_reminR_aerP - F_reminR_anaerP - F_reminR_NO3_P - F_burried) * cs(iv_detrR_P)           &  
                              +F_transf_MO * cs(iv_detr_P)
            
             diag_3d_sed(id_remin_aerP,k,i,j)=diag_3d_sed(id_remin_aerP,k,i,j)  &
                       +(F_remin_aerP * cs(iv_detr_P) + F_reminR_aerP * cs(iv_detrR_P)) * dzs(k,i,j)

             diag_3d_sed(id_remin_anaerP,k,i,j)=diag_3d_sed(id_remin_anaerP,k,i,j)  &
                       +(F_remin_anaerP * cs(iv_detr_P) + F_reminR_anaerP * cs(iv_detrR_P)) * dzs(k,i,j)

             diag_3d_sed(id_remin_nitrateP,k,i,j)=diag_3d_sed(id_remin_nitrateP,k,i,j)  &
                       +(F_remin_NO3_P * cs(iv_detr_P) + F_reminR_NO3_P * cs(iv_detrR_P)) * dzs(k,i,j)
!            ...................................................................  
             
             ! Fe bound P (micromol/l sed) 
             ! ---------------------------
             dcdt(iv_PFe)=(F_precPFE_O2 + F_precPFE_NO3) * cs(iv_nutr_PO4) * poro(k,i,j)  &
                         -F_dissolPFE * cs(iv_PFe) - F_burried * cs(iv_PFe)
                           
             diag_3d_sed(id_dissol_PFe,k,i,j)=diag_3d_sed(id_dissol_PFe,k,i,j)+  &
                                             +F_dissolPFE * cs(iv_PFe) * dzs(k,i,j)

             diag_3d_sed(id_precipit_P,k,i,j)=diag_3d_sed(id_precipit_P,k,i,j)+  &
                                             +(F_precPFE_O2 + F_precPFE_NO3) * cs(iv_nutr_PO4) * poro(k,i,j) * dzs(k,i,j)
!            ................................................................... 

             ! P adsorbe (micromol/lsed) 
             ! --------------------------
             dcdt(iv_nutr_Pads)=F_adsorP * cs(iv_nutr_PO4) * poro(k,i,j)              &
                               -F_desorP * cs(iv_nutr_Pads) - F_burried * cs(iv_nutr_Pads)

             diag_3d_sed(id_adsor_desorb_P,k,i,j)=diag_3d_sed(id_adsor_desorb_P,k,i,j)  &
                                                 +dcdt(iv_nutr_Pads) * dzs(k,i,j)
!            ...................................................................  
     
             ! Dissolved oxygen (mg/l) 
             ! -------------------------
             dcdt(iv_oxygen)= dcdt(iv_oxygen)                                          &
                            -F_precPFE_O2 * cs(iv_nutr_PO4) * p_GO2_PFe * 0.032_rsh                                                            
!            ...................................................................       
             
             ! Nitrate (micomol/l) 
             ! -------------------------
             dcdt(iv_nutr_NO3)=dcdt(iv_nutr_NO3)                                       &
                              -F_precPFE_NO3 * cs(iv_nutr_PO4) * p_GNO3_PFe                                                                        
!            ...................................................................       
 
                       !!!!!!!!!!!!!!!!!!!
                       !cycle de la silice
                       !!!!!!!!!!!!!!!!!!!
             
             ! Biogenic silicon  (micromol/lsed) 
             ! ---------------------------------
             !dcdt(iv_detr_Si)=-(F_remin_aerS + F_remin_anaerS) * cs(iv_detr_Si)                            &
             !                +diatmortsed * cs(iv_phyto_diat_N) * p_phyto_SiNratio     &
             dcdt(iv_detr_Si)=- F_remin_Si * cs(iv_detr_Si)  - F_burried * cs(iv_detr_Si)      &
                              + diatmortsed * cs(iv_phyto_diat_N) * p_phyto_SiNratio           &
#if defined key_psnz
                             +txfiltbenth * cw_bott(iv_phyto_psnz_Si) * dtbiojour              &
#endif
                             +rapsiaz * txfiltbenth * cw_bott(iv_phyto_diat_N) * dtbiojour

             diag_3d_sed(id_remin_aerSi,k,i,j)=diag_3d_sed(id_remin_aerSi,k,i,j)  &
!                                              +F_remin_aerS * cs(iv_detr_Si) * dzs(k,i,j)
                                               +F_remin_Si * cs(iv_detr_Si) * dzs(k,i,j)
!            ...................................................................       
             
             ! Dissolved silicate (micromol/l) 
             ! -------------------------------
             !dcdt(iv_nutr_SiOH)=(F_remin_aerS + F_remin_anaerS) * cs(iv_detr_Si) * porosite_inv             &
             !                  -F_precSi * cs(iv_nutr_SiOH)
             dcdt(iv_nutr_SiOH)= F_remin_Si * cs(iv_detr_Si) * porosite_inv             &
                                 - F_precSi

             !diag_3d_sed(id_precipit_Si,k,i,j)=diag_3d_sed(id_precipit_Si,k,i,j)  &
             !                                 +F_precSi * cs(iv_nutr_SiOH) * dzs(k,i,j)
             diag_3d_sed(id_precipit_Si,k,i,j)=diag_3d_sed(id_precipit_Si,k,i,j)  &
                                              +F_precSi * dzs(k,i,j)
 !            ...................................................................       
           
                    !!!!!!!!!!!!!!!!!!!!!!
                    !    phytoplankton  !!
                    !!!!!!!!!!!!!!!!!!!!!!
            
             dcdt(iv_phyto_diat_N)=-diatmortsed*cs(iv_phyto_diat_N)

             diag_3d_sed(id_morta_phyto,k,i,j)=diag_3d_sed(id_morta_phyto,k,i,j)  &
                                      +dcdt(iv_phyto_diat_N)*dzs(k,i,j)
 !            ...................................................................       


!#if defined key_zostera
                   !!!!!!!!!!!!!!!!!!!!
                   !    zostere
                   !!!!!!!!!!!!!!!!!!!

             ! Evolution de l azote detritique peu labile (mmol/m2/j)
!             dcdt(iv_detr_zost_N)=-(F_remin_aerN/10._rsh)*cs(iv_detr_zost_N)
 !            ...................................................................       
             ! Evolution du phosphore detritique peu labile (mmol/m2/j)
!             dcdt(iv_detr_zost_P)=-(F_remin_aerP/10._rsh)*cs(iv_detr_zost_P)
 !            ...................................................................       
                           
!#endif    
        
             ! update cv_sed and transfer from 1D to 3D 
             cv_sed(1:nv_adv,k,i,j)=cv_sed(1:nv_adv,k,i,j)+  &
               dcdt(1:nv_adv)*unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))

             ! introduction of rearation pour l oxygene 
             ! cycle de l oxygene
             ! si pas d eau ==> reaeration par l air en surface
             ! traite apres coup pour tenir compte des fortes consommations
            IF(htot(i,j) < RESIDUAL_THICKNESS_WAT) THEN
               F_aeration=KO2_aeration*exp(-KzO2_aeration*zmiddle)*(o2sats-cv_sed(iv_oxygen,k,i,j))*dtbiojour
               IF(id_fluxsed_aeration .ne. 0)diag_3d_sed(id_fluxsed_aeration,k,i,j)=diag_3d_sed(id_fluxsed_aeration,k,i,j)+F_aeration*dzs(k,i,j)
               cv_sed(iv_oxygen,k,i,j)=cv_sed(iv_oxygen,k,i,j)+ F_aeration
               
            ENDIF
             ! if variable bio = var.part.constitutiv ==> update poro
             IF(SUM(unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))) .NE. nv_adv) THEN
                c_sedtot(k,i,j)=0.0_rsh
                cvolp=0.0_rsh
                DO iv=1,nvpc
                  c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
                ENDDO
                poro(k,i,j)=1.0_rsh-cvolp
                IF(poro(k,i,j) .LE. 0.0_rsh) THEN
                  write(*,*)'reactions_in_sed poro<0 !! ',CURRENT_TIME,i,j,k,poro(k,i,j),cv_sed(1:nvpc,k,i,j)
                ENDIF
             ENDIF
             
            IF(k==ksma(i,j) .AND. htot(i,j) > RESIDUAL_THICKNESS_WAT .AND. txfiltbenth .NE. 0.0_rsh) THEN
               ! filtrage benthique : phytos dans la couche du fond broutes
               !dcdt(1:nv_state,1,i,j)=dcdt(1:nv_state,1,i,j)+dcw_filtbent(1:nv_state)/86400.0_rsh
#ifdef key_MARS
               IF (htot(i,j) < hm) THEN
                  cw_bottom_MUSTANG(1:nv_adv,i,j)=cw_bottom_MUSTANG(1:nv_adv,i,j)+  &
                      dcw_filtbent(1:nv_adv)*dtbiojour/htot(i,j)*unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))  
                  DO iv=1,nv_adv
                     cv_wat(iv,:,i,j)=cw_bottom_MUSTANG(iv,i,j)
                  ENDDO
               ELSE
#endif
                  cw_bottom_MUSTANG(1:nv_adv,i,j)=cw_bottom_MUSTANG(1:nv_adv,i,j)+  &
                      dcw_filtbent(1:nv_adv)*dtbiojour*unit_modif_mudbio_N2dw(irk_fil(1:nv_adv))  
#ifdef key_MARS
                  cv_wat(1:nv_adv,1,i,j)=cw_bottom_MUSTANG(1:nv_adv,i,j)
               ENDIF 
#endif
                
             ENDIF
            ENDIF ! test dzs 
         END DO  ! boucle ks
       END IF   ! test a terre
     END DO
   END DO
!$OMP END DO
   
  END SUBROUTINE bloom_reactions_in_sed
#endif

!!============================================================================================
#if defined key_MANGAbio && defined key_MANGAbiovague
     SUBROUTINE bloom_wavefile_MANGAbio(ifirst,ilast,jfirst,jlast,icall,forcSPMk)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_wavefile_MANGAbio  ***
   !&E
   !&E ** Purpose : special reading wave file for MANGAbio ( module bloom) 
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
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#ifdef key_MARS
   USE ionc4,       ONLY : ionc4_openr,ionc4_read_dimt,ionc4_read_time, &
                           ionc4_read_xyt,ionc4_close
   USE comvars2d,    ONLY : ig,id,jb,jh
#endif

   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast,icall
   REAL(KIND=rsh), DIMENSION(NB_LAYER_WAT,PROC_IN_ARRAY),INTENT(OUT),OPTIONAL  :: forcSPMk


   !! * Local declarations
   INTEGER    :: k,i,j,iv,it
   INTEGER, DIMENSION(400)       :: t_clim_vague_int
   INTEGER                       :: ijour,imois,ian,iheure,iminu,isec
   INTEGER                       :: jjulien,tool_julien
   CHARACTER(LEN=19)             :: date_start,tool_sectodat
   REAL(kind=rlg)                :: tbid
   REAL(kind=riosh), DIMENSION(COMPLETE_ARRAY) :: tab_mes

   !!--------------------------------------------------------------------------
   !! * Executable part
   
  IF (icall == 0 ) THEN
      !! first call - initialization - open and read files


     ! filevague lu dans namelist nammessat in parasubs
     !filevague_hs='../inputs/hs_MANGA.nc'
     !filevague_t02='../inputs/t02_MANGA.nc'
     !filevague_ubr='../inputs/vitesse_orbi_2008_24h_MANGA.nc'
  
   MPI_master_only WRITE(iscreenlog,*)'------------------------------------------------'
   MPI_master_only WRITE(iscreenlog,*)'lecture fichier climato pour vitesse orbitale'
   MPI_master_only WRITE(iscreenlog,*)'     donnees en m.s-1', TRIM(filevague_ubr)

#ifdef key_MARS
!!! lecture dans MARS
!!!!!!!!!!!!!!!!!!!
   CALL ionc4_openr(filevague_ubr)
   idimt_vague=ionc4_read_dimt(filevague_ubr)
  ! print*,'idimt_vague=',idimt_vague   

   DO it=1,idimt_vague
    CALL ionc4_read_time(filevague_ubr,it,t_clim_vague(it))
   END DO

   !!reading ubr
   DO it=1,idimt_vague
     CALL ionc4_read_xyt(filevague_ubr,'ubr',tab_mes,imin,imax,jmin,jmax,it)
     DO j=jfirst,jlast
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
         IF(h0(i,j).GT.-valmanq) THEN
            IF (tab_mes(i,j).LT.valmanq) THEN
            ! print*,'i j tab_mes',i,j,tab_mes(i,j)
             ubr_vague(i,j,it)=tab_mes(i,j)
             IF(ubr_vague(i,j,it) .LT. 0.0_rsh) THEN
              MPI_master_only WRITE(iscreenlog,*)'pb ubr_vague maille ',i,j,it,ubr_vague(i,j,it)
             END IF
            ELSE
             MPI_master_only WRITE(iscreenlog,*)'pas de ubr_vague climato en ',i,j,it,tab_mes(i,j)
            END IF
          END IF
       END DO
     END DO
   ! MPI_master_only WRITE(iscreenlog,*) t_clim_vague(it)
   END DO

  ! CALL ionc4_close(filevague_ubr)


   !CALL ionc4_openr(filevague_ubr)

   !!reading vbr
   DO it=1,idimt_vague
     CALL ionc4_read_xyt(filevague_ubr,'vbr',tab_mes,imin,imax,jmin,jmax,it)
     DO j=jfirst,jlast
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
         IF(h0(i,j).GT.-valmanq) THEN
            IF (tab_mes(i,j).LT.valmanq) THEN
             vbr_vague(i,j,it)=tab_mes(i,j)
             IF(vbr_vague(i,j,it) .LT. 0.0_rsh) THEN
              MPI_master_only WRITE(iscreenlog,*)'pb vbr_vague maille ',i,j,it,vbr_vague(i,j,it)
             END IF
            ELSE
             MPI_master_only WRITE(iscreenlog,*)'pas de vbr_vague climato en ',i,j,it,tab_mes(i,j)
            END IF
          END IF
       END DO
     END DO
   ! MPI_master_only WRITE(iscreenlog,*) t_clim_vague(it)
   END DO

   CALL ionc4_close(filevague_ubr)

   !MPI_master_only WRITE(iscreenlog,*)'------------------------------------------------'
   !MPI_master_only WRITE(iscreenlog,*)'lecture fichier climato pour hauteur significative'
   !MPI_master_only WRITE(iscreenlog,*)'     donnees en m', TRIM(filevague_hs)

   !CALL ionc4_openr(filevague_hs)
   !idimt_vague=ionc4_read_dimt(filevague_hs)

   !DO it=1,idimt_vague
   ! CALL ionc4_read_time(filevague_hs,it,t_clim_vague(it))
   !END DO

   !!reading
   !DO it=1,idimt_vague
   !  CALL ionc4_read_xyt(filevague_hs,'hs',tab_mes,imin,imax,jmin,jmax,it)
   !  DO j=ljmin,ljmax
   !    DO i=MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
   !      IF(h0(i,j).GT.-valmanq) THEN
   !         IF (tab_mes(i,j).LT.valmanq) THEN
   !         ! print*,'i j tab_mes',i,j,tab_mes(i,j)
   !          hs_vague(i,j,it)=tab_mes(i,j)
   !          IF(hs_vague(i,j,it) .LT. 0.0_rsh) THEN
   !           MPI_master_only WRITE(iscreenlog,*)'pb hs_vague maille ',i,j,it,hs_vague(i,j,it)
   !          END IF
   !         ELSE
   !          MPI_master_only WRITE(iscreenlog,*)'pas de hs_vague climato en ',i,j,it,tab_mes(i,j)
   !         END IF
   !       END IF
   !    END DO
   !  END DO
   ! MPI_master_only WRITE(iscreenlog,*) t_clim_vague(it)
   !END DO

   !CALL ionc4_close(filevague_hs)

#endif

#ifdef key_agrif
   IF (TEST_NOT_INITFROMFILE) THEN
     date_start=tool_sectodat(tdebagrif_messat)
     CALL tool_decompdate(date_start,ijour,imois,ian,iheure,iminu,isec)
     jjulien=tool_julien(ijour,imois,ian)-tool_julien(1,1,ian)+1
     date_s_annee=jjulien*24*60*60
   ELSE
     !date_start=tool_sectodat(CURRENT_TIME)
     date_s_annee=jjulien_BIOLINK*24*60*60
   END IF
#else
   !date_start=tool_sectodat(CURRENT_TIME)
   date_s_annee=jjulien_BIOLINK*24*60*60
#endif

   it=1
   idateinf=1
   t_clim_vague_int(1:idimt_vague)=INT(t_clim_vague(1:idimt_vague))
   IF (date_s_annee .LE. t_clim_vague_int(1)) THEN
    idatesup=1
    idateinf=idimt_vague
   ELSE
    DO WHILE ((it .LT. idimt_vague) .AND.      &
             (t_clim_vague_int(it+1) .LT. date_s_annee))
    it=it+1
    idateinf=it
    END DO
    IF (idateinf .EQ. idimt_vague) THEN
      idatesup=1
    ELSE
      idatesup=idateinf+1
    END IF
   END IF

   MPI_master_only WRITE(iscreenlog,*) 'Fin de lecture fichier climato vague'
   MPI_master_only WRITE(iscreenlog,*) it,' champs vague lus'
   MPI_master_only WRITE(iscreenlog,*) ' champs vague interp entre', idateinf,' et ',idatesup
   MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',date_s_annee
   MPI_master_only WRITE(iscreenlog,*)'----------------------------------------------------'
   MPI_master_only WRITE(*,*)''


  ELSE IF (icall == 1) THEN
     !!  during run :: evaluation of waves and SPM a time t

           ! ==================================================================
           ! calcul du facteur d interpolation des vagues lues dans fichiers de climato vagues
           ! ==================================================================

       IF (date_s_annee .NE. jjulien*24*60*60) THEN ! Definit des dates inf et sup pour l interp
        date_s_annee= jjulien*24*60*60              ! de la climato
        IF (idateinf .LT. idatesup) THEN
         IF (date_s_annee .GT. INT(t_clim_vague(idatesup))) THEN
          idateinf=idatesup
          IF (idatesup .NE. idimt_vague) THEN
           idatesup=idatesup+1
          ELSE
           idatesup=1
          END IF
         END IF
        ELSE
         IF ((date_s_annee .GT. INT(t_clim_vague(idatesup))) .AND.   &
               (date_s_annee .LT. INT(t_clim_vague(idateinf)))) THEN
          idateinf=1
          idatesup=idatesup+1
         END IF
        END IF

       END IF !endif on date_s_annee

       ! mis en commentaire car interp n est pas utilise - on prend la valeur a dateinf
       !IF (idatesup .GT. idateinf)  THEN
       ! interp= (date_s_annee -t_clim_vague(idateinf))/        &
       !          (t_clim_vague(idatesup)-t_clim_vague(idateinf))
       !ELSE
       ! IF (date_s_annee .GE. INT(t_clim_vague(idateinf))) THEN
       ! interp= (date_s_annee-t_clim_vague(idateinf))/&
       !           ((366*24*60*60+t_clim_vague(idatesup))-t_clim_vague(idateinf))
       ! ELSE
       ! interp=(date_s_annee+366*24*60*60-t_clim_vague(idateinf))/ &
       !           ((366*24*60*60+t_clim_vague(idatesup))-t_clim_vague(idateinf))
       ! END IF
       !END IF
!
        !print*, 'Simul a la date',date_s_annee,jjulien,imois,ijour,iheure
        !print*, 'interpolation entre ',idateinf,'et',idatesup
        !print*, 'soit, en s, entre ',INT(t_clim_vague(idateinf)),' et ',INT(t_clim_vague(idatesup))
        !print*, 'interp= ',interp

!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast
#ifdef key_MARS
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
       DO i=ifirst,ilast
#endif
         IF (TOTAL_WATER_HEIGHT .LE. RESIDUAL_THICKNESS_WAT) THEN
            forcSPMk(:,i,j)=0.0_rsh
         ELSE

            ! on cree un gradient de turbidite surface-fond du a la 
            !     remise en suspension par les vagues
            ! -----------------------------------------------
            ubr_interp=ubr_vague(i,j,idateinf)!+(ubr_vague(i,j,idatesup)-ubr_vague(i,j,idateinf))*interp
            vbr_interp=vbr_vague(i,j,idateinf)!+(vbr_vague(i,j,idatesup)-vbr_vague(i,j,idateinf))*interp
      
            IF (TOTAL_WATER_HEIGHT.le.200.0_rsh) THEN
               ufond_vague=sqrt(ubr_interp*ubr_interp+vbr_interp*vbr_interp)
            ELSE
              ufond_vague=0.0_rsh
            END IF
            DO k=1,NB_LAYER_WAT
               forcSPMk(k,i,j)=(forcSPMk(k,i,j)*k+max(forcSPMk(k,i,j),(ufond_vague*200.0_rsh)) &
                                        *(NB_LAYER_WAT-k))/NB_LAYER_WAT
            ENDDO
         END IF  !endif on htot
        END DO
      END DO
!$OMP END DO

  ENDIF
  
   END SUBROUTINE bloom_wavefile_MANGAbio
#endif

!!============================================================================================
#endif
 END MODULE
