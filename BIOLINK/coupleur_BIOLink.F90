!
#include "cppdefs.h"
!---------------------------------------------------------------------------
!
                     MODULE coupleur_BIOLink
!
!---------------------------------------------------------------------------
   

#if defined BIOLink

   !&E======================================================================
   !&E                   ***  MODULE  BIOLink  ***
   !&E
   !&E ** Purpose : concerns coupling Biogeochemical ocean/coastal/sediment model BIOLink with hydro code
   !&E Ocean dynamics Bio :  Initializations, reading of files *.dat (rivers, discharges, bio variables...)
   !&E                       + update sources and sinks terms
   !&E
   !&E ** Description :
   !&E     subroutine BIOLink_initialization   ! initialization of BIOLink  - routine called by main
   !&E     subroutine BIOLink_init             ! specifics initializations of BIO modules  - routine called by main 
   !&E     subroutine BIOLink_update           ! evaluation of sinks and sources terms  - routine called by step 
   !&E
   !&E     subroutine BIOLink_read_vardiag     ! reading diagnostics variables - called by BIOLink_initialization
   !&E     subroutine BIOLink_water_column     ! evaluation of total water height and vertical meshes thickness - called by BIOLink_update
   !&E     subroutine BIOLink_convarray        ! conversion of 3D or 4D array from hydro model to BIOLink - called by BIOLink_update
   !&E     subroutine BIOLink_sinking_rate     ! update and limiting sinking rate for each variables - called by BIOLink_update
   !&E     subroutine BIOLink_eval_PAR         ! evaluation of solar radiation extinction, attenuation and PAR - called by BIOLink_update
   !&E
   !&E     subroutine BIOLink2hydro            !  conversion array BIOLink   to hydro host model 
   !&E     subroutine BIOLink_updateconc_BIO   !  change tracers-substances concentrations if needed
   !&E     subroutine BIOLink_substance        ! definition number et type of variables (if key_nosubstmodule )

   !&E     subroutine BIOLink_SPMsat_file       ! special reading SPM satellite file for key_messat

   !&E ** Purpose :  

   !&E
   !&E

   !&E
   !&E   History :
   !&E    !  2019-08 (B. Thouvenin) issued from bloom and peptic for portability adaptation
   !&E
   !&E===================================================================================================================

#include "coupleur_define_BIOLink.h"

#ifdef key_MARS
#include "coupleur_dimhydro_BIOLink.h"
   USE sflxatm,      ONLY : rad
#else
   USE module_BIOLink
   USE comsubstance
#endif
   USE comBIOLink
#if defined MUSTANG
   USE comMUSTANG ,  ONLY : htot
#endif
#if defined ECO3M   
   USE COUPLEUR_PHY_BIO
#endif

 
   IMPLICIT NONE

   !! * Interface
  
   
   !! * Accessibility

   ! functions & routines of this module, called outside :
    PUBLIC BIOLink_initialization,BIOLink_init                ! initialization of BIOLink  - routines called by main.F90
    PUBLIC BIOLink_update                                     ! routine called by step.F90
    PUBLIC BIOLink_read_vardiag


#ifdef key_MARS
    PUBLIC BIOLink_exchgMPI_cvwat
#endif
#if defined key_nosubstmodule
    PUBLIC BIOLink_substance
#endif

   PRIVATE
#if ! defined key_MARS
   !! * Shared or public module variables (variables used by MUSTANG but issued from hydro model or substances module  )
   !! * for MARS_MODEL, these variables are stored in comvar.. or comsubstance which are not the same for other model

#endif
 
 CONTAINS

  !!======================================================================

  SUBROUTINE BIOLink_initialization(icall)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_initialization  ***
  !&E
  !&E ** Purpose : intialization of  module  BIOLink
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : init
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E       !  2019-08  (B.Thouvenin)  creation 
  !&E
  !&E---------------------------------------------------------------------
  !! * Modules used
#if defined PEPTIC
  USE peptic_initdefine, ONLY : peptic_param,peptic_alloc_var
#elif defined BLOOM
  USE bloom_initdefine, ONLY : bloom_param,bloom_init_iv
#elif defined METeOR
  USE meteor_initdefine, ONLY : meteor_param
#endif

  !! * Arguments
  INTEGER, INTENT(IN) :: icall
   
  !! * Local declarations
  INTEGER               :: isubs
  !REAL(KIND=rlg)        :: tool_datosec
   
   ! reading namelist in parapetic.txt or parabloom.txt
   ! in MARS already called in init.F90 (or in some test case)
   !!!!!    To Program BIO    !!!!!!!!!
   !----------------------------------------------------------------

   IF (icall==0) THEN

#if ! defined key_MARS
#if defined PEPTIC
    CALL peptic_param('r')
#elif defined BLOOM
    CALL bloom_param('r')
#elif defined METeOR
    CALL meteor_param('r')
#endif
#endif

   ELSE
#if ! defined key_MARS
     ! in MARS : TIME_BEGIN = tdeb known and declared by hydro model
     ! in CROCO : TIME_BEGIN ???   => TIME_BEGIN is declared in BIOLink
     TIME_BEGIN=CURRENT_TIME  

     ! in MARS already called in subreaddat

     ! allocation of sources and sinks terms  common for modules ECO
     ! (thicklayerW - htot, BIO_SINKSOURCES, PAR_ATTENUATED_RAD
     ! EXTINCTION_RAD, PAR_top_layer_ave,light_integ
     !----------------------------------------------------------------
     CALL BIOLink_alloc()
       
  ! identification of variable numbers or specific allocation 
  !!!!!    To Program BIO    !!!!!!!!!
  !  ----------------------------------------------------------
#ifdef BLOOM
     DO isubs=1,nv_adv
       CALL bloom_init_iv(isubs,standard_name_var(isubs),1)
     END DO
     DO isubs=nv_adv+1,nv_adv+nv_fix
       CALL bloom_init_iv(isubs,standard_name_var_fix(isubs-nv_adv),2)
     END DO
#ifdef key_benthic
     DO isubs=nv_adv+nv_fix+nv_bent
       CALL bloom_init_iv(isubs,standard_name_var_bent(isubs-nv_adv-nv_fix),3)
     END DO
#endif
#endif
#if defined PEPTIC
     CALL peptic_alloc_var
#endif
#if defined ECO3M
     CALL ALLOC_VAR_Eco3M
     nbcallbio = -1
     CALL main_bio(0.0)
#endif
  
  ! writing BIO used parameters  (To Program BIO )
  !-----------------------------------------------
#if defined PEPTIC
     CALL peptic_param('w')
#elif defined BLOOM
     CALL bloom_param('w')
#elif defined METeOR
     CALL meteor_param('w')
#endif

#endif /* ! defined key_MARS */

  ! reading bio diagnostics variables - 
  !  allocation et reperage pour BLOOM
  !-----------------------------------------------
#if defined BLOOM
     CALL BIOLink_read_vardiag
#endif

   ENDIF  ! /* icall */

   PRINT_DBG*, 'END BIOLink_initialization'


  END SUBROUTINE BIOLink_initialization

  !!======================================================================
 
   SUBROUTINE BIOLink_init(ifirst,ilast,jfirst,jlast)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_init  ***
  !&E
  !&E ** Purpose : Local initialisations, reading of special files
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : main
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E       !  2015-07  
  !&E
  !&E---------------------------------------------------------------------
  !! * Modules used
#if defined BLOOM
   USE bloom_initdefine, ONLY : bloom_userinit
#if defined key_MANGAbio && defined key_MANGAbiovague
   USE bloom,            ONLY : bloom_wavefile_MANGAbio
#endif
#endif
#if defined METeOR
   USE meteor_initdefine, ONLY : meteor_read_react
#endif

  IMPLICIT NONE

  !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast   !,NB_LAYER_WAT

  !! * Local declarations
  !REAL(KIND=rlg)        :: tool_datosec
  INTEGER               :: i,j,k,iv
  CHARACTER(LEN=lchain) :: logfilename


  !!----------------------------------------------------------------------
  !! * Executable part

  ! locating the start time of the calculation of bio sources and sinks
  !---------------------------------------------------------------
#ifdef key_MARS
   t_bio=MAX(CURRENT_TIME,tdebsubs)!!+TRANSPORT_TIME_STEP
#else
   t_bio=CURRENT_TIME+TRANSPORT_TIME_STEP
#endif

#if ! defined MUSTANG
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! initialization of sinking_rate 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast
      DO i=ifirst,ilast
        DO iv=1,nvp
          DO k=1,NB_LAYER_WAT
              WS_BIOLink(k,iv,i,j)=(ws_free_min(iv)+ws_free_max(iv))/2.0_rsh
            ENDDO
          ENDDO
        END DO
      END DO
!$OMP ENDDO
#endif

#if defined BIOLink_verif_conserv 
!  ouverture du fichier log pour test sur conservativite
   !---------------------------------------------------------------
     logfilename='conservBIOLink_'//TRIM(suffix_fileres)//'.log'
     !IF_AGRIF (.NOT. agrif_root()) logfilename= 'conservBIOLink_'//TRIM(suffix_fileres)//'_'//TRIM(Agrif_Cfixed())//'.log'
     IF (iscreenlog_conserv/=6) OPEN(unit=iscreenlog_conserv,file=logfilename)
#endif
   
  ! specifiques Initializations (To Program BIO)
  !--------------------------------------------

#if defined BLOOM
   !CALL bloom_userinit(ifirst,ilast,jfirst,jlast,WATER_CONCENTRATION, &
   !                        BATHY_H0,NB_LAYER_WAT,CURRENT_TIME)
   CALL bloom_userinit(ifirst,ilast,jfirst,jlast)
#endif
#if defined PEPTIC
   !CALL peptic_iniFrac_plkt ! not operational
#endif

#if defined METeOR
    ! reading reactions between contaminant species
    CALL meteor_read_react
#endif

   ! Read special files (To Pogram BIO)
   !-----------------------------------

#if defined key_MANGAbio && defined key_MANGAbiovague
   ! initialization ubr et hs
   ubr_vague(:,:,:)=0.0_rsh
   vbr_vague(:,:,:)=0.0_rsh
   hs_vague(:,:,:)=0.0_rsh
   DO j=jfirst,jlast
#ifdef key_MARS
     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
     DO i=ifirst,ilast
#endif
       IF (BATHY_H0(i,j) .EQ. -valmanq) THEN
         ubr_vague(i,j,:)=valmanq
         vbr_vague(i,j,:)=valmanq
         hs_vague(i,j,:)=valmanq
       END IF
     END DO
   END DO
   
   !!! fichier climato des vitesses orbitales de vague
   !!! routine dans le coupleur (dependant du code)
   CALL bloom_wavefile_MANGAbio(ifirst,ilast,jfirst,jlast,0)
#endif


#if defined key_messat
   ! initialization satellite MES
   messat(:,:,:)=0.0_rsh
   DO j=jfirst,jlast
#ifdef key_MARS
     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
     DO i=ifirst,ilast
#endif
       IF (BATHY_H0(i,j) .EQ. -valmanq) THEN
         messat(i,j,:)=valmanq
       END IF
     END DO
   END DO

   IF (l_messat_clim .and. l_messat_obs) THEN
     IF_MPI (MASTER) THEN
       MPI_master_only WRITE(ierrorlog,*)'ATTENTION : Les climato MES Sat et les observations MES'   &
      //' Sat ne peuvent pas etre utilisees simultanement. Voir dans parasubs.txt'
       MPI_master_only WRITE(ierrorlog,*)'The simulation is stopped'
       CALL_MPI MPI_FINALIZE(IERR_MPI)
       STOP
     ENDIF_MPI
   END IF

   ! read  climato MES file
   !----------------------------------------
   CALL BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,0)

#endif

  END SUBROUTINE BIOLink_init

  !!======================================================================

  SUBROUTINE BIOLink_update(ifirst,ilast,jfirst,jlast   &
#if defined key_MARS && (defined key_oyster_SFG || defined key_oyster_DEB)
          , CELL_SURF                                                &
#endif
         )
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_update  ***
  !&E
  !&E ** Purpose : update  at each time step 
  !&E
  !&E            ******* ATTENTION  we are in a session OMP  **************
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : init
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E       !  2019-08  (B.Thouvenin)  creation 
  !&E
  !&E---------------------------------------------------------------------
   !! * Modules used
#if defined PEPTIC
  USE peptic,            ONLY : peptic_sksc_wat, peptic_SPMtot_Chla
#elif defined BLOOM
  USE bloom,            ONLY : bloom_sksc_wat,bloom_eval_diag2d, bloom_SPMtot_Chla, &
                               bloom_extinction_avg
#if defined key_MANGAbio && defined key_MANGAbiovague
  USE bloom,            ONLY : bloom_wavefile_MANGAbio
#endif
#elif defined ECO3M
  USE COUPLEUR_PHY_BIO
#endif
#if defined key_MARS 
#if defined key_turbclim && defined key_daily_climato_kpar
    USE comvars2d,     ONLY : mes_sat
#endif   
#elif defined METeOR
  USE meteor,          ONLY : meteor_sksc_wat,meteor_reac_equi
#else

#endif   
   !! * Declaration Subroutine

   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast
#if defined key_MARS && (defined key_oyster_SFG || defined key_oyster_DEB)
   REAL(KIND=rsh),DIMENSION(ARRAY_CELL_SURF),INTENT(IN)          :: CELL_SURF
#endif
   
   !! * Local declarations
   REAL(KIND=rlg)         :: dt_bio_cor
   INTEGER                :: i,j,k,kmaxmod,iv,itend,tool_julien
   CHARACTER(LEN=19)        :: tool_sectodat,cdate
   REAL(KIND=rsh), DIMENSION(PROC_IN_ARRAY)   :: forcSPM
#if defined key_MANGAbio && defined key_MANGAbiovague
   REAL(KIND=rsh), DIMENSION(NB_LAYER_WAT,PROC_IN_ARRAY)  :: forcSPMk
#endif 

    itend=0
   IF(CURRENT_TIME .GE. t_bio) THEN
   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! EVALUATION of sinks and sources from bio reactions   !!!!!!
      !!!!  UPDATE ALL ECO_TIME_STEP                             !!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP SINGLE
      dt_bio_cor=MAX(ECO_TIME_STEP,TRANSPORT_TIME_STEP)
      BIO_TIME_STEP=dt_bio_cor+(CURRENT_TIME-t_bio)
!$OMP END SINGLE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  time tracking (day, month, year, Julian day, etc.)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cdate = tool_sectodat(CURRENT_TIME)
      CALL tool_decompdate(cdate,ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK, &
                           iheure_BIOLINK,iminu_BIOLINK,isec_BIOLINK) 
      jjulien_BIOLINK=tool_julien(ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK)-tool_julien(1,1,ian_BIOLINK)+1


      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! coupleur BIOLink : evaluation of total water height and vertical meshes thickness
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL BIOLink_water_column(ifirst,ilast,jfirst,jlast)
      
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!  Conversion array in the index order of BIOLink   ==>temp, sal, ws3
     !!      and keep only positive water concentration ==> cvadv_wat_pos
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL BIOLink_convarray(ifirst,ilast,jfirst,jlast     &
#if defined key_nosubstmodule
                  ,WAT_SETTL                                     &
#endif
#ifdef key_MARS
                  ,TEMPERATURE_MOD,SALINITY_MOD                  &
#endif
                  )     

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! limitation of the sinking speed and actualization of the sink speed for plankton
      !!
      !!! module BLOOM  : sinking_rate evaluated during Sources and sinks terms evaluation  
      !!!                 because needed of some terms calculated localy (bloom_settling)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL BIOLink_sinking_rate(ifirst,ilast,jfirst,jlast)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  specific initializations
    !   To Program BIO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined BLOOM && defined key_oyster_DEB
! Reinitialisation des huitres au 1er janvier de chaque annee.
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)
     DO j=jfirst,jlast
      DO i=ifirst,ilast

        IF ((imois_BIOLINK.eq.1).and.(ijour_BIOLINK.eq.1).and.(iheure_BIOLINK.eq.0)) THEN
          IF (nbhuitre(i,j).ne.0.0_rsh) THEN
             k=1
             iv=iv_oysdeb_res
             WAT_CONCFIX_ivkij=50.0_rsh
             iv=iv_oysdeb_gon
             WAT_CONCFIX_ivkij=500.0_rsh
             iv=iv_oysdeb_str
             WAT_CONCFIX_ivkij=310.0_rsh
          ENDIF
        ENDIF
! Mortalite huitres par cohortes
        nbhuitre(i,j)=nbhuitre(i,j)-txmorthuitco1(imois_BIOLINK)*nbhuitre(i,j)*BIO_TIME_STEP/86400.

      ENDDO
     ENDDO
!$OMP END DO
#endif  

#if defined BIOLink_PAR_eval
      !!=================================================================================
      !!!    steps to evaluate PAR via SPM, CHLORO, EXTINCTION...
      !!            DONE BY BIOLink
      !!        otherwise the BIO module take charge PAR estimation
      !!=================================================================================

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  update SPM evaluation from satellite or model  ==> forcSPM
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP SINGLE
      forcSPM(:,:)=0.0_rsh

#if defined key_messat
      ! SPM read from satellite input file
      CALL BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,1,forcSPM=forcSPM)
#else
      ! SPM read from climato file
#if defined key_MARS

#if defined key_turbclim && defined key_daily_climato_kpar
     !! mes_sat lu dans sflxsurf avec le KPAR pour le code hydro
      forcSPM(:,:)=mes_sat(:,:)
#endif

#else
    ! To Program HYDRO
    ! not MARS : define SPM if forced by readed value (satellite or other mordel or ..)
     ! forcSPM(:,:)=
#endif

#endif

#if defined key_MANGAbio && defined key_MANGAbiovague
     !! forced SPM  with wave action which create a surface_bottom gradient
      CALL bloom_wavefile_MANGAbio(ifirst,ilast,jfirst,jlast,1,forcSPMk=forcSPMk)
     
#endif
!$OMP END SINGLE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! concentrations of mineral and organic suspension particles and cholorophylle
      !!!  for extinction and PAR estimation ( To Program BIO )
      !!!   ***** ATTENTION*** unity of chloro depending on bio module ??? 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined PEPTIC    
      CALL peptic_SPMtot_Chla(ifirst,ilast,jfirst,jlast)

#elif defined BLOOM
      CALL bloom_SPMtot_Chla(ifirst,ilast,jfirst,jlast          &
#if defined key_MANGAbio && defined key_MANGAbiovague
                             ,forcSPMk                          &
#endif 
#if defined key_messat
                              ,forcSPM                          &
#endif
                         )

#endif  

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! solar radiation extinction/attenuation / PAR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate)

#if defined BLOOM
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! module BLOOM : extinction average on 4 days
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL bloom_extinction_avg(ifirst,ilast,jfirst,jlast)
#endif

#if defined PEPTIC && defined key_rand_extinc
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! module PEPTIC : abondance verification (not ready)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL peptic_abondance(ifirst,ilast,jfirst,jlast)  
#endif

      !!=================================================================================
#endif

#if defined BLOOM
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! module bloom : evaluation de certaines variables diagnostiques 2D 
      !!! (production primaire cumulee,  phytos max, position et date...)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL  bloom_eval_diag2d(ifirst,ilast,jfirst,jlast)

#endif

#if defined key_oxygen || defined key_zostera || defined METeOR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NEED TO KNOW WIND SPEED
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   in  MARS : known
!   in CROCO known only if BULK_FLUX
#if ! defined key_MARS && ! defined BULK_FLUX

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)
     DO j=jfirst,jlast
     DO i=ifirst,ilast

        ! wind speed (m.s-1) evaluated from the surface stress values
        ! sustr and svstr  (m2.s-2)
        ! with rho_air=1.3  ; CD=0.0014  
        WIND_SPEED(i,j) = sqrt(sqrt( (0.5*(sustr(i,j)+sustr(i+1,j)))**2
     &                       +(0.5*(svstr(i,j)+svstr(i,j+1)))**2)
     &                       *RHOREF/(1.3*0.0014))
!     &                       *rho0/(rho_air*CD))
      ENDDO
      ENDDO
!$OMP END DO

#endif
#endif
 
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !   calcul des sources et puits biologiques   To Program BIO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined PEPTIC
     CALL peptic_sksc_wat(ifirst,ilast,jfirst,jlast)

#elif defined BLOOM

     CALL bloom_sksc_wat(ifirst,ilast,jfirst,jlast                          &
#if defined BLOOM
#if defined GLS_MIXING
                          ,CIN_TURBULENT_ENERGY                                    &
#endif
#if defined key_zostera || defined key_oxygen
                          ,WIND_SPEED                                              &
#endif
#if defined key_MARS && (defined key_oyster_SFG || defined key_oyster_DEB)
                       ,CELL_SURF                                            &
#endif
#endif
                                 ) 
#elif defined METeOR
    IF ( l_treat_re_dc) THEN
        CALL meteor_sksc_wat(ifirst,ilast,jfirst,jlast,WIND_SPEED)
    ENDIF                      
#endif
   

#if ! defined key_MARS 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Conversion array BIOLink   to hydro host model
   !!   for settling velocities if changed by bio module
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL BIOLink2hydro(ifirst,ilast,jfirst,jlast)     
#endif
     
!$OMP SINGLE
     t_bio=t_bio+BIO_TIME_STEP
!$OMP END SINGLE
     itend=1
 
   ENDIF  ! t>t_bio

#if defined BIOLink_UPDATE_CONCBIO 
   ! a chaque pas de temps, mise a jour des concentrations apres transformations bio
   ! si modele hydro ne le prend pas en charge en resolvant les equations de conservation de la masse
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Conversion array BIOLink   to hydro host model
   !!   for settling velocities if changed by bio module
   !!     and resolution of dc/dt=bio_sink_sources if it is not taking into account in transport equations
   !!     then change water concentrations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined METeOR && ! defined BLOOM && ! defined PEPTIC && ! defined ECO3M
    IF ( l_treat_re_dc) THEN
#endif   

      CALL BIOLink_updateconc_BIO(ifirst,ilast,jfirst,jlast)

#if defined METeOR && ! defined BLOOM && ! defined PEPTIC && ! defined ECO3M
    ENDIF 
#endif
#endif



#if defined METeOR
   ! Contaminant - equilibrium or very fast reactions ==> new concentrations
   ! a each time step
   ! -----------------------------------------------------------------------
       IF ( nreacmax > 0 .AND. l_treat_re_eq) THEN
        IF (itend == 0)THEN
         cdate = tool_sectodat(CURRENT_TIME)
         CALL BIOLink_water_column(ifirst,ilast,jfirst,jlast)
         CALL BIOLink_convarray(ifirst,ilast,jfirst,jlast)
#if defined BIOLink_PAR_eval
!$OMP SINGLE
         forcSPM(:,:)=0.0_rsh
!$OMP END SINGLE
         CALL BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate)
#endif
        ENDIF
       
        CALL meteor_reac_equi(ifirst,ilast,jfirst,jlast,WIND_SPEED)
 
       ENDIF
#endif



     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! evolution des variables fixees :
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !  les variables fixees sont calculees directement 
           !  pour les autres variables (celles qui sont advectees), 
           !            soit elles sont calculees directement (dans routine BIOLink2hydro)
           !            soit les termes sources et puits sont pris en compte dans les equations de conservation de la masse
     IF (nv_fix > 0 ) THEN
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)
      DO j=jfirst,jlast
     DO i=ifirst,ilast
      ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
        kmaxmod=NB_LAYER_WAT

          DO k=1,kmaxmod
             FIXED_VAR_CONC(FIXED_VAR_INDEXkij)=FIXED_VAR_CONC(FIXED_VAR_INDEXkij) &
                                       +TRANSPORT_TIME_STEP*BIO_SKSC_FIX(FIXED_SKSC_INDEXkij)

          ENDDO
        ENDDO
       ENDDO
!$OMP END DO
     ENDIF

    PRINT_DBG*, 'END BIOLink_update'


  END SUBROUTINE BIOLink_update

  !!======================================================================
  SUBROUTINE BIOLink_water_column(ifirst,ilast,jfirst,jlast)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_water_column ***
  !&E
  !&E ** Purpose : memorize total water hight and vertical water layer thicknesses 
  !&E
  !&E
  !&E            ******* ATTENTION  we are in a OMP session   **************
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
  !&E---------------------------------------------------------------------
  !! * Modules used
#ifdef key_MARS
   USE toolgeom,     ONLY : f_dzu,f_dzw
   USE comvars2d,    ONLY : ig,id,jb,jh,hm
#else
#endif
   !! * Arguments 
   INTEGER, INTENT(IN)                           :: ifirst,ilast,jfirst,jlast
 

  !! * Local declarations
    INTEGER                  :: i,j,k,kmaxmod

    

  !!----------------------------------------------------------------------
  !! * Executable part

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)
   DO j=jfirst,jlast
#ifdef key_MARS
      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
        IF(j.GE.jb(i)+1 .AND. j .LE. jh(i)-1) THEN
#else
      DO i=ifirst,ilast
      ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
#endif

#if ! defined MUSTANG
          !  htot : total water height
          ! already knwon if MUSTANG
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MARS
          TOTAL_WATER_HEIGHT(i,j)=BATHY_H0(i,j)+WATER_ELEVATION(i,j)
#else
        !! To Program HYDRO
        !! CROCO : 
           TOTAL_WATER_HEIGHT(i,j)=z_w(i,j,N)+h(i,j)
#endif
#endif
            ! thickness of  water layer
            !  thickness_C : centered around C, T, S
            !  thickness_W : centered around Wz, kz, interface
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MARS
            IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN
               kmaxmod=1
            ELSE
 :diff              kmaxmod=NB_LAYER_WAT
            ENDIF

            IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN
                THICKLAYERWC(1,i,j)=TOTAL_WATER_HEIGHT(i,j)
                THICKLAYERWW(1,i,j)=TOTAL_WATER_HEIGHT(i,j)*0.5_rsh
                THICKLAYERWC(2:NB_LAYER_WAT,i,j)=0.0_rsh
                THICKLAYERWW(2:NB_LAYER_WAT,i,j)=THICKLAYERWW(1,i,j)
            ELSE
              DO k=1,kmaxmod          
                THICKLAYERWC(k,i,j)=f_dzu(BATHY_H0(i,j),WATER_ELEVATION(i,j),k,i,j)
                THICKLAYERWW(k,i,j)=f_dzw(BATHY_H0(i,j),WATER_ELEVATION(i,j),k,i,j)
              ENDDO
            ENDIF
        ENDIF
#else
        !! To Program HYDRO
        !! CROCO : 
          DO k=1,NB_LAYER_WAT-1          
                THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1)
                THICKLAYERWW(k,i,j)=z_r(i,j,k+1)-z_r(i,j,k)
          ENDDO
          k=NB_LAYER_WAT
          THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1)
          THICKLAYERWW(k,i,j)=0._rsh
#endif

      ENDDO
   ENDDO
!$OMP END DO

 
  END SUBROUTINE  BIOLink_water_column

   !!======================================================================
  SUBROUTINE BIOLink_convarray(ifirst,ilast,jfirst,jlast  &
#if defined key_nosubstmodule
                  ,WAT_SETTL                                     &
#endif
                        )

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_convarray ***
  !&E
  !&E ** Purpose : conversion of 3D or 4D array from hydro model to BIOLink
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
  !&E---------------------------------------------------------------------
  !! * Modules used
#ifdef key_MARS
   USE toolgeom,     ONLY : f_dzu,f_dzw
   USE comvars2d,    ONLY : ig,id,jb,jh,hm
#endif
   !! * Arguments 
   INTEGER, INTENT(IN)                           :: ifirst,ilast,jfirst,jlast

  !! * Local declarations
    INTEGER                  :: i,j,k,kmaxmod,iv

    

  !!----------------------------------------------------------------------
  !! * Executable part

#if defined PEPTIC || defined BLOOM 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modules BIO PEPTIC and BLOOM : index ordre : iv,k,i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,iv)
     DO j=jfirst,jlast
     
#ifdef key_MARS
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
       DO i=ifirst,ilast
#endif

          kmaxmod=NB_LAYER_WAT
          WATCONCPOS(:,:,i,j)=0.0_rsh
          FIXCONCPOS(:,:,i,j)=0.0_rsh
#ifdef key_benthic
          BENTCONCPOS(:,i,j)=0.0_rsh
#endif
          DO k=1,kmaxmod
#if ! defined key_MARS
            ! conversion index order if not MARS
            TEMP_BIOLink(k,i,j)=TEMPHYDRO_ijk
            SAL_BIOLink(k,i,j)=SALHYDRO_ijk
#endif
#if defined BLOOM && defined key_benthos
            bottom_current(i,j)=BOTTOM_CURRENT_ij
#endif
            ! variables d etat transportees
            DO iv=1,nv_adv
              WATCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCADV_ivkij+ABS(WAT_CONCADV_ivkij))
            END DO
            ! variables fixees (connu via comsubstance)
            DO iv=1,nv_fix
              FIXCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCFIX_ifixkij+ABS(WAT_CONCFIX_ifixkij))
            END DO
!#if ! defined key_MARS
            ! conversion index order if not MARS si on a cree un tableau pour avoir le bon ordre
          !  DO iv=1,nv_adv
          !    BIO_SKSC_ADV(iv,k,i,j)=BIO_SINKSOURCES(ADV_VAR_INDEXkij)
          !  ENDDO
!#endif
          ENDDO
#ifdef key_benthic
            ! variables benthiques (connu via comsubstance)
            DO iv=1,nspb
              BENTCONCPOS(iv,i,j)=0.5_rsh*(WAT_CONCBENT_ivij+ABS(WAT_CONCBENT_ivij))
            END DO
#endif

#ifdef key_oyster_DEB
! Reinitialisation des huitres au 1er janvier de chaque annee.
        IF ((imois_BIOLink.eq.1).and.(ijour_BIOLink.eq.1).and.(iheure_BIOLink.eq.0)) THEN
          IF (nbhuitre(i,j).ne.0.0_rsh) THEN
             k=1
             iv=iv_oysdeb_res
             BENTHIC_CONCENTRATION(BENTH_INDEXij)=50.0_rsh
             iv=iv_oysdeb_gon
             BENTHIC_CONCENTRATION(BENTH_INDEXij)=500.0_rsh
             iv=iv_oysdeb_str
             BENTHIC_CONCENTRATION(BENTH_INDEXij)=310.0_rsh
          ENDIF
        ENDIF
! Mortalite huitres par cohortes
        nbhuitre(i,j)=nbhuitre(i,j)-txmorthuitco1(imois_BIOLink)*nbhuitre(i,j)*dtbio/86400.
#endif  

        END DO
      END DO
!$OMP END DO

#if ! defined key_MARS
!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast   
      DO i=ifirst,ilast
        DO iv=1,nv_adv
          DO k=1,NB_LAYER_WAT
            ! conversion index order if not MARS
              WS_BIOLink(k,iv,i,j)=WAT_SETTL_ivkij
            ENDDO
          ENDDO
        END DO
      END DO
!$OMP END DO
#endif

!!!!!!!!!! END MODULES BIO PEPTIC &&& BLOOM

#elif defined METeOR && ! defined ECO3M
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modules METeOR : index ordre : iv,k,i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,iv)
     DO j=jfirst,jlast
       DO i=ifirst,ilast
          kmaxmod=NB_LAYER_WAT
          WATCONCPOS(:,:,i,j)=0.0_rsh
          FIXCONCPOS(:,:,i,j)=0.0_rsh
          DO k=1,kmaxmod
            ! conversion index order if not MARS
            TEMP_BIOLink(k,i,j)=TEMPHYDRO_ijk
            SAL_BIOLink(k,i,j)=SALHYDRO_ijk
            ! variables d etat transportees
            DO iv=1,nv_adv
              WATCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCADV_ivkij+ABS(WAT_CONCADV_ivkij))
            END DO
            ! variables fixees (connu via comsubstance)
            DO iv=1,nv_fix
              FIXCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCFIX_ifixkij+ABS(WAT_CONCFIX_ifixkij))
            END DO
          ENDDO 
        END DO
      END DO
!$OMP END DO

!!!!!!!!!! END METeOR

#elif defined ECO3M

#endif

   END SUBROUTINE  BIOLink_convarray

   !!======================================================================
   SUBROUTINE BIOLink_sinking_rate(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_sinking_rate  ***
   !&E
   !&E              Biologic dynamics:  - limiting sinking rate for each variables
   !&E
   !&E ** Purpose :  update sinking rate 
   !&E              
   !&E ** Description :  
   !&E   
   !&E  !!!!!!!!!!!! CHANGE IF NEEDED LOOP ORDER                      !!!!!!!!!!!!!!!!
   !&E  !!!!!!!!!!!!   IF WAT_SETTL indexes are in a different order  !!!!!!!!!!!!!!!!
   !&E  !!!!!!!!!!!!  than in MARS (iv,k,i,j)                         !!!!!!!!!!!!!!!!
   !&E
   !&E
   !&E ** Called by :
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-08 (B. Thouvenin) issued from peptic_dynzwat - verti_quota (A. Menesguen, M. Sourrisseau)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
#ifdef key_MARS
   USE comvars2d,    ONLY : ig,id,hm
   USE obccombine, ONLY : l_obc_south, l_obc_north, l_obc_west, l_obc_east
#endif

   !! * Declaration Subroutine

   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast !,kmax
!   REAL(KIND=rsh),INTENT(IN)                                  :: h0fond 
!   REAL(KIND=rlg),INTENT(IN)                                  :: dt_true,CURRENT_TIME
   
   !! * Local declarations
   INTEGER                :: i,j,k,iv
#ifdef PEPTIC
   INTEGER                :: i_plkt,i_quota_loc,ind
#endif

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,iv)
   DO j=jfirst,jlast
#ifdef key_MARS
      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
      DO i=ifirst,ilast
#endif

        IF (TOTAL_WATER_HEIGHT(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN

            ! =================================================================
            ! Calcul par defaut des vitesses de chute  de toute variable particulaire (si pas MUSTANG)
            !  sans modulation ni particularite (si different : calcul dans module BIO)
            ! =================================================================
            ! si MUSTANG : les vitesses de chute sont calculees dans le module MUSTANG
            ! ----------------------------------------------------------------
#if ! defined MUSTANG
           DO iv=1,nvp
             DO k = NB_LAYER_WAT,1,-1
              WS_BIOLink(k,iv,i,j)=(ws_free_min(iv)+ws_free_max(iv))/2.0_rsh
             ENDDO
           ENDDO
#endif
 
           ! limitation of the sinking speed and actualisation of the sink speed for plankton
           !   To Program BIO
           ! ---------------------------------------------------------------------------
#ifdef PEPTIC
           !correction pour les diffs quotas
           DO i_plkt = 1 , bd_fp%nb_plct
             DO i_quota_loc = 1 , plct(i_plkt)%nb_quota
               IF (plct(i_plkt)%is_quota_var(i_quota_loc)== 1) THEN !only advected var
                 iv=plct(i_plkt)%num_mod_mars(i_quota_loc) 
                 DO k = NB_LAYER_WAT,1,-1
                   WS_BIOLink(k,iv,i,j) = plct(i_plkt)%sink_max  !m.d-1
                 ENDDO
               ENDIF
             ENDDO
           ENDDO
#endif

         ! Bornage
         !!!!!!!!!!!!!
           DO iv=nvpc+1,nvp
                 DO k = NB_LAYER_WAT,1,-1
                   WS_BIOLink(k,iv,i,j)=sign(MIN(0.95_rlg*THICKLAYERWC(k,i,j)/TRANSPORT_TIME_STEP,REAL(ABS(WS_BIOLink(k,iv,i,j)),rlg)),WS_BIOLink(k,iv,i,j))
                 ENDDO
           END DO  !iv
        ENDIF
      END DO  !i
    END DO  !j
!$OMP END DO
         

  END SUBROUTINE BIOLink_sinking_rate

   !!======================================================================
#ifdef BIOLink_PAR_eval
  SUBROUTINE BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_eval_PAR ***
  !&E
  !&E ** Purpose : evaluation of solar radiation extinction, attenuation and PAR 
  !&E
  !&E ** Description : issu de verti_quota
  !&E
  !&E ** Called by : 
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E
  !&E---------------------------------------------------------------------
    !! * Modules used


    !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast
   CHARACTER(LEN=19),INTENT(IN)                               :: cdate

    !! * Local declarations
    INTEGER                  ::  i,j,k,iv,kmaxmod        ! loop indexes
#ifdef PEPTIC
    INTEGER                  ::  i_plkt,i_pom,ind,num_cell,i_quota_loc         ! loop indexes
    REAL(KIND=rsh),DIMENSION(bd_fp%nb_plct) :: conc
#endif
    REAL(KIND=rsh),DIMENSION(NB_LAYER_WAT)     :: attenuation



   !Pour quota Chl
   !integration lumiere sur 24h
   INTEGER                  :: iday,jhour,numday_lum,numhour_lum,klum
   INTEGER                  :: numday_extinction,numhour_extinction



!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,attenuation)
   DO j=jfirst,jlast
#ifdef key_MARS
     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
       IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN
               kmaxmod=1
       ELSE
               kmaxmod=NB_LAYER_WAT
       ENDIF

#else
     DO i=ifirst,ilast
      ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
        kmaxmod=NB_LAYER_WAT
#endif

       !diag_3d_wat(:,:,i,j)=0.0_rsh
       PAR_top_layer(:,i,j)=0.0_rsh
       EXTINCTION_RAD(:,i,j)=0.0_rsh
       attenuation(:)=0.0_rsh
#if defined METeOR
       Flimrad_layer(:,i,j)=0.0_rsh
#endif
       !diag_3d_wat(1,NB_LAYER_WAT,i,j)=PAR_top_layer(NB_LAYER_WAT,i,j)

       IF (TOTAL_WATER_HEIGHT(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN

        ! loop from surface to bottom of water column
         DO k = LOOPK_SURF_TO_BOTTOM_WAT   ! kmaxmod,1,-1

           !absorption water--------------------------------------------------------------------------------
           !----------------
           EXTINCTION_RAD(k,i,j) = PARAM_WAT_EXTINCT

           !absorption with SPMtot -------------------------------------
#ifdef PEPTIC
           EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j) + bd_fp%extincspim * cmes_3dmgl(k,i,j)
#elif defined BLOOM || (defined METeOR && ! defined PEPTIC)
           EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j) + p_extincspim * (cmes_3dmgl(k,i,j)+epsilon_BIOLink)
#endif  


           !absorption chlorophylle -----------------------------------------------------------------------
           !------------------------
#ifdef PEPTIC
           IF (BIOLink_chloro(k,i,j) > 1.0e-10_rsh) THEN
#endif

#if defined PEPTIC || defined BLOOM
            !! ATTENTION chloro unity
            ! module BLOOM : in mug/L
            ! module PEPTIC : in mg/L
             EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j) + PARAM_CHLORO1_EXTINCT * ( BIOLink_chloro(k,i,j) ** PARAM_CHLORO2_EXTINCT )
#ifdef PEPTIC
           ENDIF
#endif
#endif
           
  
           !absorption macro algae -----------------------------------------------------------------------
           !----------------------
#if defined BLOOM && defined key_zostera
           ! - absorption due aux zosteres ---------------------------------------
           IF(k==1 .and. FIXCONCPOS(iv_zost_LB-nv_adv,1,i,j).gt.0.0_rsh) THEN
             EXTINCTION_RAD(k,i,j)=EXTINCTION_RAD(k,i,j)+p_zost_leafabscoef*FIXCONCPOS(iv_zost_LB-nv_adv,1,i,j)*p_zost_klai
           ENDIF
#endif
       
           !estimation of attenuation at layer k ----------------------------------------------------------
           !-------------------------------------
           attenuation(k) = EXP( -EXTINCTION_RAD(k,i,j) * THICKLAYERWC(k,i,j))
         ENDDO !k
           
         ! Estimation of PAR at the top of each layer (W/m2) -----------------------------------------------------
         !             To Program BIO
         !---------------------------------------------------
#ifdef PEPTIC
#if defined key_growth_diurne 
         PAR_top_layer(kmaxmod,i,j)=SOLAR_RAD(i,j)*bd_fp%parradratio / RAD_SRFSCALE
#else
         ! PAR_top_layer(k,i,j)=light_ave_daily(i,j)  ! non connu ????????
#endif
#elif defined BLOOM
         PAR_top_layer(kmaxmod,i,j)=SOLAR_RAD(i,j)*p_parradratio / RAD_SRFSCALE 
#elif defined METeOR
         ! here , only attenuation of radiation is required (for reaction which depend on extinction) 
         Flimrad_layer(kmaxmod,i,j)=attenuation(kmaxmod)
#endif

         DO k=LOOPK_SUBSURF_TO_BOTTOM_WAT   ! kmaxmod-1,1,-1
            PAR_top_layer(k,i,j) = PAR_top_layer(ABOVE_K,i,j) * attenuation(ABOVE_K)       
           !diag_3d_wat(1,k,i,j)= PAR_top_layer(k,i,j)
                 !  MPI_master_only WRITE(*,*)'absorp,bd_fp%extincspim,cpom,:',absorp,bd_fp%extincspim,cpom
                 !  MPI_master_only WRITE(*,*)'bd_fp%extincChl1,BIOLink_chloro,bd_fp%extincChl2:',bd_fp%extincChl1,chloro,bd_fp%extincChl2
                 !  MPI_master_only WRITE(*,*) 'bd_fp%extincwat:',bd_fp%extincwat
#if defined METeOR
            Flimrad_layer(k,i,j)=Flimrad_layer(ABOVE_K,i,j)* attenuation(k)
#endif
                 
         ENDDO !k
         k=0 ! at bottom
         PAR_top_layer(k,i,j) = PAR_top_layer(ABOVE_K,i,j) * attenuation(ABOVE_K)

         !             To Program BIO
         ! Estimation of average PAR in each layer,                   --------------------------------------------------
         ! depending on module : for several phytoplancton species or
         !                    or for total phytoplancton (PEPTIC)
         !                    or will be evaluated in the module 
         !     (exemple module BLOOM : light effect depending on effeturbidite etc.. 
         !                             Smith formulation : one for each  phyto but depending on Ikphyto=f(effeturbidite,season..)) 
         !--------------------------------------------------------------------------------------------------------------
#if defined PEPTIC
         DO k=LOOPK_SUBSURF_TO_BOTTOM_WAT   ! kmaxmod-1,1,-1   
            
            ! module PEPTIC : total phyto ( indice 0)  
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
            klum = k - 1
            IF (k == 1) klum = 1
            IF (EXTINCTION_RAD(klum,i,j) /= 0.0_rsh) THEN
!               !passage W.m-2 = j.s-1.m-2 -> µEin.s-1.m-2 -> µEin.d-1.m-2
!              PAR_avg_layer_phyto(1,k,i,j)= PAR_top_layer(k,i,j) / EXTINCTION_RAD(klum,i,j) * ( 1.0_rsh - attenuation(klum)) &
!              / THICKLAYERWC(k,i,j) * 2.02_rsh * 86400_rsh 
#if defined key_growth_diurne
               !passage W.m-2 -> µEin.s-1.m-2 for the PAR fraction ( 400 - 700 nm)
               PAR_avg_layer_phyto(1,k,i,j) = PAR_top_layer(k,i,j) / EXTINCTION_RAD(klum,i,j) *  &  
                                         ( 1.0_rsh - attenuation(klum))/ THICKLAYERWC(k,i,j) * 4.6_rsh *86400.0_rsh 
#endif
            ENDIF
         ENDDO !k
#endif         
     
       ELSE
          PAR_top_layer(:,i,j)=0.0_rsh 
#if defined PEPTIC
          PAR_avg_layer_phyto(1,:,i,j) = 0.0_rsh
#endif
          !diag_3d_wat(1,:,i,j)=PAR_top_layer(:,i,j)
       ENDIF ! d>RESIDUAL_THICKNESS_WAT
     ENDDO !i
   ENDDO !j
!$OMP END DO


#ifdef PEPTIC
   ! module PEPTIC : estimation of average PAR on 24h on top of layer  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    !integration light / 24h ! Warning : Must restart at 00.00  and one day lag between realistic forcing and application

   IF ((iheure_BIOLINK == 0) .and. (iminu_BIOLINK == 0) .and. (isec_BIOLINK <= BIO_TIME_STEP/2._rsh)) then  !small errors for few s
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k)
      DO j = jfirst,jlast
#ifdef key_MARS
        DO i = MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
#else 
        DO i=ifirst,ilast
#endif
 !         light_ave_daily(i,j)=MAX(light_integ(i,j),0.0_rsh)/(86400.0_rsh)*bd_fp%parradratio !units - same than rad but only PAR
 !         light_integ(i,j)=0.0_rsh
          IF ((i==i_BIOLink_verif) .and. (j==j_BIOLink_verif)) THEN
           MPI_master_only WRITE(iscreenlog,*) 'new daily_aver_PAR_W_m2',cdate,BIO_TIME_STEP,PAR_top_layer_day(:,i,j)/86400.0_rsh ! error?*bd_fp%parradratio !light_ave_daily(i,j) 
          ENDIF
          PAR_top_layer_day(:,i,j)=0.0_rsh 
       ENDDO
      ENDDO
!$OMP END DO
   ELSE
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k)
      DO j = jfirst,jlast
#ifdef key_MARS
        DO i = MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
#else 
        DO i=ifirst,ilast
#endif
           PAR_top_layer_day(:,i,j)=PAR_top_layer_day(:,i,j) + PAR_top_layer(:,i,j)*BIO_TIME_STEP
!           IF (igdu(i,j)/=4) THEN
!#if defined key_siggen || defined key_gencoord
!             umoy=SUM(PAR_top_layer(:,i,j)*hzex(:,i,j)*dsigu(:))/hex(i,j)
!#else
!             umoy=SUM(PAR_top_layer(:,i,j)*dsigu(:))
!#endif
!             uz(:,i,j)=uz(:,i,j)-umoy+uv(i,j)
!           ELSE
!             uz(:,i,j)=uv(i,j)
!           END IF
 !         light_integ(i,j)=light_integ(i,j)+max(0.0,SOLAR_RAD(i,j))*BIO_TIME_STEP
        ENDDO
      ENDDO      
!$OMP END DO
    ENDIF

    !IF ( bd_fp%cle_stop ) THEN
    !  cvadv_wat_pos( plct(2)%num_mod_mars(1:plct(2)%nb_quota), :, :, :) = bd_fp%seuil_retirage
    !  cvadv_wat_pos( plct(2)%num_mod_mars(1:plct(2)%nb_quota), :, 2, 2) = 0.0_rsh
    !  bd_fp%cle_stop = .false.
    !ENDIF

#endif


END SUBROUTINE  BIOLink_eval_PAR
#endif

    !!============================================================================== 
 
  SUBROUTINE BIOLink_read_vardiag

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_read_vardiag  ***
   !&E
   !&E ** Purpose : lecture du fichier descriptif des variables diagnostiques
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
   !&E       !  2019-08 (B. Thouvenin) issued from sub_read_vardiag 
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

!! module bloom :
#ifdef BLOOM
   USE bloom_initdefine, ONLY : bloom_init_id
#ifdef key_N_tracer
#ifdef key_MARS
   USE parameters, ONLY : nb_var_tracerN 
#endif
   USE bloom_initdefine, ONLY : bloom_create_vardiagtracer
#endif
#ifdef key_P_tracer
#ifdef key_MARS 
   USE parameters, ONLY : nb_var_tracerP 
#endif
   USE bloom_initdefine, ONLY : bloom_create_vardiagtracer
#endif
#endif
!!!!!!!!!!

   !! * Arguments

   !! * Local declarations
   LOGICAL               :: ex,l_diag_wat,l_diag_sed
   INTEGER               :: eof,isubs,isubs_r,dimvar,it,ind_white,IERR_MPI
   CHARACTER(LEN=lchain) :: namvar_r,long_name_var_r,standard_name_var_r,unitvar_r
   CHARACTER(LEN=5)      :: comment
#if defined key_N_tracer || defined key_P_tracer
   INTEGER               :: is,id,ivtra
#endif


   !!----------------------------------------------------------------------
   !! * Executable part

   ! save into simu.log
   !-------------------
   IF_MPI (MASTER) THEN
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) '**************************************************'
     MPI_master_only WRITE(iscreenlog,*) '**************** VAR_READ_DIAG.F90 ***************'
     MPI_master_only WRITE(iscreenlog,*) '**************************************************'
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) 'file defining usefull diagnostic variables : ',TRIM(filevardiag)
   ENDIF_MPI

   ! Initialize number of diagnostic variables according to their dimensions
   ndiag_1d = 0
   ndiag_2d = 0
   ndiag_3d = 0
   ndiag_3d_wat = 0
   ndiag_3d_sed = 0
   ndiag_2d_sed = 0
   ndiag_tot = 0
   isubs = 0
!#ifdef key_MARS
   l_out_subs_diag=.false.  ! no saving of diagnoses variable by default
   ! in CROCO: on lit dans croco.in : ldefdiabio  qui remplace ??? a VOIR?????????????????
#ifdef MUSTANG
   l_out_subs_diag_sed=.false.  ! no saving of diagnoses variable by default
#endif
!#endif
   ! read total number of diagnostic variables
   eof = 0
   INQUIRE(file=filevardiag,exist=ex)
   IF (ex) THEN
     OPEN(49,file = filevardiag,form='formatted')
     comment='debut'
     DO
       READ(49,'(a)',iostat=eof) comment
       IF (comment.EQ.'*****') EXIT
     END DO
     DO WHILE(eof==0)
       READ(49,*,iostat=eof) ! variable number
       READ(49,'(a)',iostat=eof) ! variable name
       READ(49,'(a)',iostat=eof) !variable long_name
       READ(49,'(a)',iostat=eof) !variable standard_name
       READ(49,'(a)',iostat=eof) !unit
       READ(49,*,iostat=eof)     !variable valid_min value
       READ(49,*,iostat=eof)     !variable valid_max value
       READ(49,*,iostat=eof)     !dimension (1=1D ; 2=2D ; 3=3D)
       READ(49,*,iostat=eof)     !variable in water    (k,i,j)
       READ(49,*,iostat=eof)     !variable in sediment (i,j,k)
       READ(49,*,iostat=eof)     !saving in file
       IF (eof==0) ndiag_tot=ndiag_tot+1
       READ(49,'(a)',iostat=eof)
     END DO
     CLOSE(49)

#ifdef BLOOM
#if defined key_N_tracer
       ndiag_tracerN=0
       DO is=1,nb_source_tracerN
            ndiag_tot=ndiag_tot+nb_var_tracerN+1
            ndiag_tracerN=ndiag_tracerN+nb_var_tracerN+1
#if defined key_age_tracer
            ndiag_tot=ndiag_tot+nb_var_tracerN+1
            ndiag_tracerN=ndiag_tracerN+nb_var_tracerN+1
#endif
       ENDDO
#endif
#if defined key_P_tracer
       ndiag_tracerP=0
       DO is=1,nb_source_tracerP
            ndiag_tot=ndiag_tot+nb_var_tracerP+1
            ndiag_tracerP=ndiag_tracerP+nb_var_tracerP+1
#if defined key_age_tracer
            ndiag_tot=ndiag_tot+nb_var_tracerP+1
            ndiag_tracerP=ndiag_tracerP+nb_var_tracerP+1
#endif
       ENDDO
#endif
#endif

     IF_MPI (MASTER) THEN
        MPI_master_only WRITE(iscreenlog,*) 'Number of diagnostic variables = ',ndiag_tot
     ENDIF_MPI
   ! allocate arrays for diagnostic variables
     ALLOCATE( idimv_r(ndiag_tot) )
     ALLOCATE( l_diagBIOLink_out(ndiag_tot) )
     ALLOCATE( name_vardiag(ndiag_tot) )
     ALLOCATE( long_name_vardiag(ndiag_tot) )
     ALLOCATE( standard_name_vardiag(ndiag_tot) )
     ALLOCATE( unit_vardiag(ndiag_tot) )
     ALLOCATE( valid_min_vardiag(ndiag_tot) )
     ALLOCATE( valid_max_vardiag(ndiag_tot) )
     ALLOCATE( irk_diag(ndiag_tot) )


   ! read diagnostic variables and order them according to their matrix dimensions
     OPEN(49,file = filevardiag,form='formatted')
     comment='debut'
     DO WHILE (comment /= '*****')
       READ(49,'(a)',iostat=eof) comment
     END DO
     DO WHILE(eof==0)
       READ(49,*,iostat=eof) isubs_r
       isubs=isubs+1
       IF (eof==0) THEN
         READ(49,'(a)',iostat=eof) namvar_r
         ind_white=INDEX(namvar_r,' ')
         name_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(namvar_r(1:ind_white))))
         READ(49,'(a)',iostat=eof) long_name_var_r
         ind_white=INDEX(long_name_var_r,' ') 
         long_name_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(long_name_var_r(1:ind_white))))
         READ(49,'(a)',iostat=eof) standard_name_var_r
         ind_white=INDEX(standard_name_var_r,' ')
         standard_name_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(standard_name_var_r(1:ind_white))))
         READ(49,'(a)',iostat=eof) unitvar_r
         ind_white=INDEX(unitvar_r,' ')
         unit_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(unitvar_r(1:ind_white))))
         READ(49,*,iostat=eof) valid_min_vardiag(isubs)
         READ(49,*,iostat=eof) valid_max_vardiag(isubs)
         READ(49,*,iostat=eof) dimvar
         IF (dimvar==1) THEN
           idimv_r(isubs)=1
           ndiag_1d=ndiag_1d+1
         ELSE IF (dimvar==2) THEN
           idimv_r(isubs)=2
           ndiag_2d=ndiag_2d+1
         ELSE IF (dimvar==3) THEN
           idimv_r(isubs)=3
           ndiag_3d=ndiag_3d+1
         END IF
         READ(49,*,iostat=eof) l_diag_wat
         IF (l_diag_wat .and. dimvar==3) ndiag_3d_wat=ndiag_3d_wat+1
         READ(49,*,iostat=eof) l_diag_sed
         IF (l_diag_sed .and. dimvar==2) THEN
           ndiag_2d_sed=ndiag_2d_sed+1
           idimv_r(isubs)=6
         END IF
         IF (l_diag_sed .and. dimvar==3) THEN
           ndiag_3d_sed=ndiag_3d_sed+1
           IF (l_diag_wat) THEN
             idimv_r(isubs)=4
           ELSE
             idimv_r(isubs)=5
           ENDIF
         END IF
         READ(49,*,iostat=eof) l_diagBIOLink_out(isubs)
         READ(49,*,iostat=eof)
         IF (l_diag_sed) THEN
#if ! defined MUSTANG
          IF_MPI (MASTER) THEN
           MPI_master_only WRITE(iscreenlog,*)' '
           MPI_master_only WRITE(iscreenlog,*)' WARNING : diagnostic variable in sediment '
           MPI_master_only WRITE(iscreenlog,*)'           without CPP key MUSTANG '
           MPI_master_only WRITE(iscreenlog,*)' DIAG. VAR. NAME : ',TRIM(name_vardiag(isubs))
          ENDIF_MPI
#else
          IF (l_diagBIOLink_out(isubs)) THEN
             l_out_subs_diag_sed=.true.
          ENDIF
#endif
         ELSE
           IF (l_diagBIOLink_out(isubs))  l_out_subs_diag=.true.
         ENDIF
         
         IF_MPI (MASTER) THEN
           MPI_master_only WRITE(iscreenlog,*)' '
           MPI_master_only WRITE(iscreenlog,*)' DIAG. VAR. NAME : ',TRIM(name_vardiag(isubs))
           MPI_master_only WRITE(iscreenlog,*)' UNIT            : ',TRIM(unit_vardiag(isubs))
           MPI_master_only WRITE(iscreenlog,*)' DIMENSION       : ', dimvar,'D'
           MPI_master_only WRITE(iscreenlog,*)' DIAG SAVE ?     : ', l_diagBIOLink_out(isubs)
           IF (l_diag_wat)THEN
              MPI_master_only WRITE(iscreenlog,*)' DIAGNOSTIC INTO THE SEA'
           ENDIF
           IF (l_diag_sed) THEN
              MPI_master_only WRITE(iscreenlog,*)' DIAGNOSTIC INTO THE SEDIMENT'
           ENDIF
         ENDIF_MPI
       END IF

     END DO
     CLOSE(49)
   ELSE
     MPI_master_only WRITE(iscreenlog,*)' filevardiag', TRIM(filevardiag), ' not found  !! '
   END IF     ! fin du test sur l existence du fichier

#if defined BLOOM && (defined key_N_tracer || defined key_P_tracer)
       call bloom_create_vardiagtracer
#endif

   ! save into simu.log
   !-------------------
   IF_MPI (MASTER) THEN
     MPI_master_only WRITE(iscreenlog,*)' '
     MPI_master_only WRITE(iscreenlog,*)' STOCK OF DIAGNOSTIC VARIABLES'
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables total :',ndiag_tot
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 1 dim :',ndiag_1d
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 2 dim :',ndiag_2d
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 3 dim :',ndiag_3d
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 3d wat:',ndiag_3d_wat
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 2d sed:',ndiag_2d_sed
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 3d sed:',ndiag_3d_sed
   ENDIF_MPI

   IF (ndiag_tot /= ndiag_1d+ndiag_2d+ndiag_3d) THEN
      MPI_master_only PRINT*,'WARNING number of diagnostic variables is incoherent.'
      MPI_master_only PRINT*,'Check in file :',filevardiag
   END IF

#if ! defined MUSTANG
   IF (ndiag_3d_sed /= 0 .OR. ndiag_2d_sed /= 0 ) THEN
      MPI_master_only PRINT*,'WARNING no MUSTANG, you should not have any diagnostic variable for sediment'
      MPI_master_only PRINT*,'simulation stopped'
      CALL_MPI MPI_FINALIZE(IERR_MPI)
      STOP
   END IF
#endif

   ! allocate diagnostic variables
   ALLOCATE( diag_1d(1:ndiag_1d) )
   ALLOCATE( diag_2d(ndiag_1d+1:ndiag_1d+ndiag_2d,PROC_IN_ARRAY) )
   ALLOCATE( diag_3d_wat(ndiag_2d+1:ndiag_2d+ndiag_3d_wat,NB_LAYER_WAT,PROC_IN_ARRAY) )
   diag_1d(:)=0.0_rsh
   diag_2d(:,:,:)=0.0_rsh
   diag_3d_wat(:,:,:,:)=0.0_rsh
#if defined MUSTANG && defined key_BLOOM_insed
   ALLOCATE( diag_3D_sed(ndiag_tot-ndiag_3d_sed+1:ndiag_tot,ksdmin:ksdmax,PROC_IN_ARRAY) ) 
   diag_3d_sed(:,:,:,:)=0.0_rsh
   ALLOCATE( diag_2D_sed(ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d,PROC_IN_ARRAY) )
   diag_2D_sed(:,:,:)=0.0_rsh
#endif

   ! Storage of diagnostic variables within reading order
#if defined BLOOM
   it = 0
   ! dimvar=1 : diag1D, =2 diag2D
   ! dimvar=3 : diag3D in wat only
   ! dimvar=4 : diag 3D in wat and in sed
   ! dimvar=5 : diag3D in sed only
   ! dimvar=6 : diag2D in sed only
#if defined key_BLOOM_insed
   ! dimvar=1 : diag1D, =2 diag2D
   DO dimvar = 1,2
     DO isubs = 1,ndiag_tot
       IF (idimv_r(isubs) == dimvar) THEN
         it=it+1
!         irk_diag(it)=isubs
         irk_diag(isubs)=it
#if defined key_N_tracer
         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#elif defined key_P_tracer
         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#else
         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#endif
       END IF
     END DO
   END DO
   ! dimvar=6 : diag2D in sed only
   dimvar=6
   DO isubs = 1,ndiag_tot
       IF (idimv_r(isubs) == dimvar) THEN
         it=it+1
!         irk_diag(it)=isubs
         irk_diag(isubs)=it
#if defined key_N_tracer && defined BLOOM
         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#elif defined key_P_tracer && defined BLOOM
         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#else
         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#endif
       END IF
   END DO
   ! dimvar=3 : diag3D in wat only
   ! dimvar=4 : diag 3D in wat and in sed
   ! dimvar=5 : diag3D in sed only
   DO dimvar = 3,5
     DO isubs = 1,ndiag_tot
       IF (idimv_r(isubs) == dimvar) THEN
         it=it+1
!         irk_diag(it)=isubs
         irk_diag(isubs)=it
#if defined key_N_tracer && defined BLOOM
         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#elif defined key_P_tracer && defined BLOOM
         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#else
         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#endif
       END IF
     END DO
   END DO

#else
   DO dimvar = 1,6
     DO isubs = 1,ndiag_tot
       IF (idimv_r(isubs) == dimvar) THEN
         it=it+1
!         irk_diag(it)=isubs
         irk_diag(isubs)=it
#if defined key_N_tracer && defined BLOOM
         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#elif defined key_P_tracer && defined BLOOM
         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#else
         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))
#endif
       END IF
     END DO
   END DO
#endif
#endif

#if defined PEPTIC
   it = 0
   DO dimvar = 1,4
    DO isubs = 1,ndiag_tot
       IF (idimv_r(isubs) == dimvar) THEN
         it=it+1
         irk_diag(isubs)=it
       END IF
    ENDDO
   ENDDO
#endif

   IF_MPI (MASTER) THEN
     DO isubs = 1,ndiag_tot
        MPI_master_only WRITE(iscreenlog,*) isubs,TRIM(name_vardiag(isubs)),idimv_r(isubs),irk_diag(isubs)
     END DO
   ENDIF_MPI

  END SUBROUTINE BIOLink_read_vardiag


   !!======================================================================
#if ! defined key_MARS 
  SUBROUTINE BIOLink2hydro(ifirst,ilast,jfirst,jlast  &
#if defined key_nosubstmodule
                               ,WAT_SETTL  &
#endif
                                 )     

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink2hydro ***
  !&E
  !&E ** Purpose : conversion of setling velocities array from BIOLink to hydro model 
  !&E                  if needed
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
  !&E---------------------------------------------------------------------
  !! * Modules used
   !! * Arguments 
   INTEGER, INTENT(IN)                           :: ifirst,ilast,jfirst,jlast
#if defined key_nosubstmodule
   REAL(KIND=rsh),DIMENSION(ARRAY_WAT_SETTL), INTENT(INOUT)   :: WAT_SETTL
#endif

  !! * Local declarations
    INTEGER                  :: i,j,k,iv
    INTEGER                  :: i1,i2,i3,i4
    REAL(KIND=rsh), DIMENSION(ARRAY_WATER_CONC0) :: xnegtr
    REAL(KIND=rsh)           :: ztra
    

  !!----------------------------------------------------------------------
  !! * Executable part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! conversion index order for settling velocities array from BIOLink to hydro host model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast
       DO i=ifirst,ilast
          DO k=1,NB_LAYER_WAT
            ! conversion index order if not MARS
            DO iv=1,nv_adv
              WAT_SETTL_ivkij=WS_BIOLink(k,iv,i,j)
            ENDDO
          ENDDO
        END DO
      END DO
!$OMP END DO


   ! If needed, echange MPI because WATER_CONCENTRATION have changed
      ! To Program HYDRO

   ! (for MARS : CALL BIOLink_exchgMPI_cvwat), but MARS include IO_SINKSURCES in mass conservation equations, so it is not necessary here)
   
END SUBROUTINE  BIOLink2hydro
#endif

   !!======================================================================
#if defined BIOLink_UPDATE_CONCBIO 
  SUBROUTINE BIOLink_updateconc_BIO(ifirst,ilast,jfirst,jlast)     

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_updateconc_BIO ***
  !&E
  !&E ** Purpose :  update concentrations (depending on hydro host model)
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
  !&E---------------------------------------------------------------------
  !! * Modules used
   !! * Arguments 
   INTEGER, INTENT(IN)                           :: ifirst,ilast,jfirst,jlast

  !! * Local declarations
    INTEGER                  :: i,j,k,iv
    INTEGER                  :: i1,i2,i3,i4
    REAL(KIND=rsh), DIMENSION(ARRAY_WATER_CONC0) :: xnegtr
    REAL(KIND=rsh)           :: ztra,sinsksourcesdt
    

  !!----------------------------------------------------------------------
  !! * Executable part


!**********************
!    To Program HYDRO
!**********************
!  2 solutions :
!     - either the hydro model integrates the terms of sources and sinks in the mass conservation equations 
!                   in this case, sources and sinks terms (BIO_SINKSOURCES) must be known by hydro host model
!                      . if the BIO_SINKSOURCES array is declared in the hydro host model (in substance module as in MARS)
!                            it is known in BIOLink module, because we do "USE comsubstance" at the beginning of the module 
!                            number and order of indexes of BIO_SINKSOURCES array are those from hydro model (ARRAY_SINKSOURCES)
!                            and it is the same in BIOLink 
!                               (we choose not to create a new array for the BIO module in order not to use memory and time for the transfer,
!                                even if this is not effective at the time of their evaluation in the BIO module 
!                                 if the order of loops i, j, k, iv is not appropriate)
!                      . if the BIO_SINKSOURCES array is declared in the BIOLink module (in comBIOLink)
!                                with the dimensions from hydro model (ARRAY_SINKSOURCES)
!                                the hydro host model must know the BIO_SINKSOURCES array and must load it when the source and sink terms 
!                                are taken into account in the mass conservation equations. (USE comBIOLink, ONLY : ..)                          
!
!     - either the BIOLink coupler updates the concentrations of the BIO variables in water column, 
!       after having calculated the terms of sources and sinks.
!                   in this case, sources and sinks terms (BIO_SINKSOURCES) are known only by BIOLink module
!                   and resolution of sink_sources equations for each variable is necessary here 
!                   in order to update concentrations modified by bio transformations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CROCO : resolution of sink_sources equations in order to have new concentrations 
!                 modified by bio transformations
!  loop in order of hydro host model because WATER_CONCENTRATION and BIO_SKSC_ADV array index
!           are in the order of hydro model (chnage the order in coupleur_define_BIOLINK.h)
!  code inspired by  PISCES  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIO_SKSC_ADV are in masse/volume/time (second), evaluated avery dt_bio time step
! but  concentrations update is evaluated every hydro time step (TRANSPORT_TIME_STEP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! test if tracers concentrations fall below 0.
         xnegtr(:,:,:) = 1.e+0
!$OMP DO SCHEDULE(RUNTIME)
         DO i1=IRANGE1
           DO i2=IRANGE2
             DO i3=IRANGE3
               DO i4=IRANGE4
                  sinsksourcesdt=BIO_SKSC_ADV(BIOSKSC_INDEX_EQ)*TRANSPORT_TIME_STEP
                  IF( (WATER_CONCENTRATION(WATCONC_INDEX_EQ)+ sinsksourcesdt ) < 0.0_rsh ) THEN
                     ztra            = ABS(  ( WATER_CONCENTRATION(WATCONC_INDEX_EQ) - epsilon_BIOLink ) &
                                               / ( sinsksourcesdt + epsilon_BIOLink ) )
                        xnegtr(i4,i3,i2) = MIN( xnegtr(i4,i3,i2),  ztra )
                  ENDIF
               END DO
             END DO
           END DO
         END DO
!$OMP END DO
         !                                ! where at least 1 tracer concentration becomes negative
         !                                ! 
!$OMP DO SCHEDULE(RUNTIME)
         DO i1=IRANGE1
           DO i2=IRANGE2
             DO i3=IRANGE3
               DO i4=IRANGE4  
                     WATER_CONCENTRATION(WATCONC_INDEX_EQ) = WATER_CONCENTRATION(WATCONC_INDEX_EQ)    &
                         &           + xnegtr(i4,i3,i2) * BIO_SKSC_ADV(BIOSKSC_INDEX_EQ)*TRANSPORT_TIME_STEP    
               END DO
             END DO
           END DO
         END DO
!$OMP END DO
     
   
END SUBROUTINE  BIOLink_updateconc_BIO
#endif


   !!===========================================================================
#ifdef key_MARS
  SUBROUTINE BIOLink_exchgMPI_cvwat

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_exchgMPI_cvwat  ***
   !&E
   !&E ** Purpose :  exchange MPI cv_wat or fixed variables
   !&E                 ATTENTION ON EST EN ZONE PARALLEL OMP
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :  sed_BIOLink_update
   !&E
   !&E ** Reference : 
   !&E
   !&E ** History :
    !&E       !  2015-12  (B.Thouvenin) reorganization of module SEDIMARS=MUSTANG
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
   USE parameters,   ONLY : liminm1,limaxp2,ljminm1,ljmaxp2

   !! * Arguments


  !! * Local declarations

   !!---------------------------------------------------------------------------
   !! * Executable part


! echange MPI

OMPMPI barrier
OMPMPI master
     CALL_MPI ex_i_rsh(-1,2,nv_state*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,cv_wat(1:nv_state,:,liminm1:limaxp2,ljminm1:ljmaxp2))
     CALL_MPI ex_j_rsh(-1,2,nv_state*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,cv_wat(1:nv_state,:,liminm1:limaxp2,ljminm1:ljmaxp2))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(cv_wat)


  END SUBROUTINE BIOLink_exchgMPI_cvwat
#endif
   !!======================================================================
   
#if defined key_nosubstmodule
  SUBROUTINE BIOLink_substance(icall)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE coupl_BIOLink_dim  ***
   !&E
   !&E ** Purpose : If there is no substance module installed in the code, 
   !&E      this routine makes it possible to define the number and the type of 
   !&E      variables (substances) and to classify them in the desired order by MUSTANG
   !&E             (not used in CROCO)
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
   !&E       !  2016-11  (B. Thouvenin) 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER,INTENT(IN)  :: icall

   !! * Local declarations
   INTEGER    :: iv

   !!--------------------------------------------------------------------------
   !! * Executable part

   IF(icall==0) THEN


     ! definition of the number of simulated variables (substances)  

     !==========================================================================
     ! definition here ??? 
     ! or read a data file or a namelist file 
     !==========================================================================
     ! number of particulate variables  type gravel 
      nv_grav=0
     ! number de particulate variables type sand 
      nv_sand=1
     ! number de particulate variables type mud 
      nv_mud=2
     ! number de particulate variables type sorbed on another partculate variable 
      nv_sorb=0
     ! number of dissolved  variables 
      nv_dis=1


      ! initialize the number of variables according to their type
      nvpc=nv_mud+nv_sand+nv_grav
      nvp=nvpc+nv_ncp+nv_sorb
      nv_adv=nvp+nv_dis
      nv_state=nv_adv
      nv_tot=nv_state
    
    ELSE
    !  icall =1 apres sed_alloc 

      ALLOCATE(name_var(1,nv_tot))
      ALLOCATE(typart(nv_state))
      ALLOCATE(typdiss(nv_adv))
      ALLOCATE(itypv(nv_tot)
      ALLOCATE(irk_fil(nv_tot))
      IF (nvp > 0) THEN
        ALLOCATE(ws_free_opt(nvp)) 
        ALLOCATE(ws_hind_opt(nvp)) 
        ALLOCATE(ws_free_para(4,nvp)) 
        ALLOCATE(ws_free_min(nvp)) 
        ALLOCATE(ws_free_max(nvp)) 
        ALLOCATE(ws_hind_para(2,nvp)) 
        ALLOCATE(tocd(nvp))
        ws_free_min(:)=0.0_rsh
        ws_free_max(:)=0.0_rsh
        ws_free_para(:,:)=0.0_rsh
        ws_hind_para(:,:)=0.0_rsh
        ws_free_opt(:)=0
        ws_hind_opt(:)=0
        tocd(:)=0.0_rsh
      ENDIF
      ALLOCATE(l_subs2D(nv_adv))
      ALLOCATE(ws3(ARRAY_WAT_SETTL))
      ALLOCATE(cv_wat(ARRAY_WATER_CONC))

! and info on  variables
    ! variable names (name_var), type (itypv), settling velocities  (ws_free..), 
    ! itypv=1 for gravels ; =2 pour sand ; =3 pour mud, 
    ! itypv=4 for non constitutives particulate variables, 
    ! itypv=5 for sorbed variables  
    ! itypv=6 for dissoolved
    ! for particulate, give diameter, tocd, ros 
    ! for muds and non constit particulates, give parameters to define settling velocities 
    ! initial concentrations  in water and  sediments
    ! for sand : treatment in 2D or 3D (l_subs2D)
    ! for sorbed particulate, give the indices or numbers of the constitutives 
    !                         part. variables on which they are sorbed (irkm_var_assoc)
    ! identification of igrav1,igrav2,isand1,isand2,imud1,imud2
    ! ---------------------------------------------------------
      igrav1=1
      igrav2=nv_grav
      isand1=igrav2+1
      isand2=igrav2+nv_sand
      imud1=isand2+1
      imud2=isand2+nv_mud 

      typart(1:nv_state)=0.
      typart(1:nvpc)=1.
      typdiss(nvp+1:nv_adv)=1.
      typdiss(1:nvp)=0.

      l_subs2D(:)=.false.

     ! we arrange here to give the variables in the right order for MUSTANG 
     ! so that they are well ranked in the table of concentrations.
     ! so :
      DO iv=1,nv_tot
        irk_fil(iv)=iv
        unit_modif_mudbio_N2dw(iv)=1.0_rsh
      ENDDO
  

   END SUBROUTINE BIOLink_substance
#endif

   !!==============================================================================

    !!==============================================================================
#if defined key_messat
     SUBROUTINE BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,icall,forcSPM)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_SPMsat_file  ***
   !&E
   !&E ** Purpose : special reading SPM satellite file for key_messat
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
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast
   INTEGER, INTENT(IN)                                        :: icall
   REAL(KIND=rsh), DIMENSION(PROC_IN_ARRAY),INTENT(OUT),OPTIONAL  :: forcSPM

   !! * Local declarations
   INTEGER    :: k,i,j,iv,it
   INTEGER                       :: ijour,imois,ian,iheure,iminu,isec
   INTEGER                       :: jjulien_init,tool_julien
   CHARACTER(LEN=19)             :: date_start,tool_sectodat
   REAL(kind=rlg)                :: tbid
   REAL(kind=riosh), DIMENSION(COMPLETE_ARRAY) :: tab_mes
   INTEGER, DIMENSION(51)        :: t_clim_messat_int
   CHARACTER(LEN=lchain)         :: nomfic
   INTEGER                       :: IERR_MPI 
   REAL(kind=rlg)                :: torigin,tool_datosec
   !!--------------------------------------------------------------------------
   !! * Executable part
   
 IF (icall == 0) THEN
     !! first call - initialization - open and read files
        
#ifdef key_MARS

   IF (l_messat_clim .EQV. .TRUE.) THEN
     IF_MPI (MASTER) THEN
       MPI_master_only WRITE(iscreenlog,*)'------------------------------------------------'
       MPI_master_only WRITE(iscreenlog,*)'lecture fichier climato pour mes-sat  '
       MPI_master_only WRITE(iscreenlog,*)'     donnees en g.l-1', TRIM(filemessat_clim)
     ENDIF_MPI
     CALL ionc4_openr(filemessat_clim)
     idimt_messat=ionc4_read_dimt(filemessat_clim)
     DO it=1,idimt_messat
       CALL ionc4_read_time(filemessat_clim,it,t_clim_messat(it))
     END DO

     !!reading
     DO it=1,idimt_messat
       CALL ionc4_read_xyt(filemessat_clim,'MES',tab_mes,imin,imax,jmin,jmax,it)
       DO j=jfirst,jlast
#ifdef key_MARS
         DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
         DO i=ifirst,ilast
#endif
          IF(h0(i,j).GT.-valmanq) THEN
            IF (tab_mes(i,j).LT.valmanq) THEN
             messat(i,j,it)=tab_mes(i,j)
             IF(messat(i,j,it) .LT. 0.0_rsh) THEN
                   MPI_master_only WRITE(iscreenlog,*)'pb messat maille ',i,j,it,messat(i,j,it)
             END IF
            ELSE
                   MPI_master_only WRITE(iscreenlog,*)'pas de MES climato en ',i,j,it,tab_mes(i,j)
            END IF
          END IF
         END DO
       END DO
       MPI_master_only WRITE(iscreenlog,*) t_clim_messat(it)
     END DO

     CALL ionc4_close(filemessat_clim)

#ifdef key_agrif
     IF (.NOT. l_initfromfile) THEN
       date_start=tool_sectodat(tdebagrif_messat)
       CALL tool_decompdate(date_start,ijour,imois,ian,iheure,iminu,isec)
       jjulien_init=tool_julien(ijour,imois,ian)-tool_julien(1,1,ian)+1
       date_s_annee_messat=jjulien_init*24*60*60
     ELSE
       !date_start=tool_sectodat(CURRENT_TIME)
       date_s_annee_messat=jjulien_BIOLINK*24*60*60
     END IF
#else
     !date_start=tool_sectodat(CURRENT_TIME)
     date_s_annee_messat=jjulien_BIOLINK*24*60*60
#endif
     it=1
     idateinf_messat=1
     t_clim_messat_int(1:idimt_messat)=INT(t_clim_messat(1:idimt_messat))
     IF (date_s_annee_messat .LE. t_clim_messat_int(1)) THEN
        idatesup_messat=1
        idateinf_messat=idimt_messat
     ELSE
        DO WHILE ((it .LT. idimt_messat) .AND.      &
             (t_clim_messat_int(it+1) .LT. date_s_annee_messat))
          it=it+1
          idateinf_messat=it
        END DO
        IF (idateinf_messat .EQ. idimt_messat) THEN
          idatesup_messat=1
        ELSE
          idatesup_messat=idateinf_messat+1
        END IF
     END IF
     IF_MPI (MASTER) THEN
       MPI_master_only WRITE(iscreenlog,*) 'Fin de lecture fichier climato satellite MES'
       MPI_master_only WRITE(iscreenlog,*) it,' champs MES lus'
       MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
       MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',date_s_annee_messat
       MPI_master_only WRITE(iscreenlog,*)'----------------------------------------------------'
       MPI_master_only WRITE(iscreenlog,*)''
     ENDIF_MPI

   END IF !endif sur l_messat_clim

   IF (l_messat_obs) THEN

    IF_MPI (MASTER) THEN
      MPI_master_only WRITE(iscreenlog,*)'------------------------------------------------'
      MPI_master_only WRITE(iscreenlog,*)'lecture fichier obs moyennes quizaines pour mes-sat  '
      MPI_master_only WRITE(iscreenlog,*)'     donnees en g.l-1', TRIM(filemessat_obs)
    ENDIF_MPI

    CALL ionc4_openr(filemessat_obs)
    idimt_messat_obs=ionc4_read_dimt(filemessat_obs)
!   idimt_messat=26*6  !attention entree en dur de la taille du fichier- a modif
    ALLOCATE(messat_obs(limin:limax,ljmin:ljmax,idimt_messat_obs))
    ALLOCATE(t_obs_messat(idimt_messat_obs))
#if defined key_MANGAbio
   ! include time lag if time origine is not 01/01/1900 in file
    torigin=tool_datosec('01/01/1900 00:00:00')-tool_datosec('01/01/1998 00:00:00')
#endif

    DO it=1,idimt_messat_obs
      CALL ionc4_read_time(filemessat_obs,it,t_obs_messat(it))
    END DO
    t_obs_messat(:)=t_obs_messat(:)-torigin
    MPI_master_only WRITE(*,*)'tobd_messat debut et fin=',t_obs_messat(1),t_obs_messat(idimt_messat_obs)

    !!reading
    DO it=1,idimt_messat_obs
      CALL ionc4_read_xyt(filemessat_obs,'MES_SAT',tab_mes,imin,imax,jmin,jmax,it)
!     CALL ionc4_read_xyt(filemessat_clim,'MES_SAT',tab_mes,limin,limax,ljmin,ljmax,it)
     DO j=jfirst,jlast
#ifdef key_MARS
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
       DO i=ifirst,ilast
#endif

         IF(h0(i,j).GT.-valmanq) THEN
            IF (tab_mes(i,j).LT.valmanq) THEN
             messat_obs(i,j,it)=tab_mes(i,j)
             IF(messat_obs(i,j,it) .LT. 0.0_rsh) THEN
              MPI_master_only WRITE(iscreenlog,*)'pb messat_obs maille ',i,j,it,messat_obs(i,j,it)
             END IF
            ELSE
             MPI_master_only WRITE(iscreenlog,*)'pas de MES climato en ',i,j,it,tab_mes(i,j)
            END IF
          END IF
       END DO
      END DO
      MPI_master_only WRITE(iscreenlog,*) t_obs_messat(it)
    END DO

    CALL ionc4_close(filemessat_obs)

    it=1
    idateinf_messat=1
    tbid=CURRENT_TIME
    IF (l_initfromfile .and. l_init_rtime) then
       CALL ionc4_openr(file_init)
       CALL ionc4_read_time(file_init,1,tbid)
       CALL ionc4_close(file_init)
    END IF
    date_start=tool_sectodat(tbid)
    MPI_master_only WRITE(iscreenlog,*)'date_start=',date_start
    IF (tbid .LE. t_obs_messat(1)) THEN
       IF_MPI (MASTER) THEN
          MPI_master_only WRITE(ierrorlog,*)'la date de depart est anterieure a la premiere date du fichier messat'
          CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
    ELSE
       DO WHILE ((it .LT. idimt_messat_obs) .AND.(t_obs_messat(it+1) .LT. tbid))
          it=it+1
          idateinf_messat=it
       END DO
       IF (idateinf_messat .EQ. idimt_messat_obs) THEN
         IF_MPI (MASTER) THEN
            MPI_master_only WRITE(ierrorlog,*)'la date de depart est posterieure a la derniere date du fichier messat'    
            CALL_MPI MPI_FINALIZE(IERR_MPI)
            STOP
         ENDIF_MPI
       ELSE
         idatesup_messat=idateinf_messat+1
       END IF
    END IF
   
    IF_MPI (MASTER) THEN
     MPI_master_only WRITE(iscreenlog,*) 'Fin de lecture fichier observation quinzaine satellite MES'
     MPI_master_only WRITE(iscreenlog,*) it,' champs MES lus'
     MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
     MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',tbid
     MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', t_obs_messat(idateinf_messat),' et ',t_obs_messat(idatesup_messat)
     MPI_master_only WRITE(iscreenlog,*)'----------------------------------------------------'
     MPI_master_only WRITE(iscreenlog,*)''
     MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',tbid
     MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
     MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', t_obs_messat(idateinf_messat),' et ',t_obs_messat(idatesup_messat)
    ENDIF_MPI

   END IF !endif sur l_messat_obs
#endif


  ELSE IF (icall == 1) THEN
     !!  during run :: evaluation of SPM a time t

           ! ==================================================================
           ! calcul du facteur d interpolation des MES lues dans fichiers de donnees satellitales
           ! ==================================================================

   IF (l_messat_clim) THEN

     IF (date_s_annee_messat .NE. jjulien_BIOLINK*24*60*60) THEN ! Definit des dates inf et sup pour l interp
        date_s_annee_messat= jjulien_BIOLINK*24*60*60              ! de la climato
        IF (idateinf_messat .LT. idatesup_messat) THEN
          IF (date_s_annee_messat .GT. INT(t_clim_messat(idatesup_messat))) THEN
            idateinf_messat=idatesup_messat
            IF (idatesup_messat .NE. idimt_messat) THEN
              idatesup_messat=idatesup_messat+1
            ELSE
              idatesup_messat=1
            END IF
          END IF
        ELSE
          IF ((date_s_annee_messat .GT. INT(t_clim_messat(idatesup_messat))) .AND.   &
               (date_s_annee_messat .LT. INT(t_clim_messat(idateinf_messat)))) THEN
             idateinf_messat=1
             idatesup_messat=idatesup_messat+1
          END IF
        END IF

!        MPI_master_only WRITE(iscreenlog,*) 'Simul a la date',date_s_annee_messat,jjulien_BIOLINK,imois_BIOLINK,ijour_BIOLINK
!        MPI_master_only WRITE(iscreenlog,*) 'interpolation entre ',idateinf_messat,'et',idatesup_messat
!        MPI_master_only WRITE(iscreenlog,*) 'soit, en s, entre ',INT(t_clim_messat(idateinf_messat)), &
!                               ' et ',INT(t_clim_messat(idatesup_messat))

     END IF !endif on date_s_annee_messat

     IF (idatesup_messat .GT. idateinf_messat)  THEN
         interp= (date_s_annee_messat -t_clim_messat(idateinf_messat))/        &
                 (t_clim_messat(idatesup_messat)-t_clim_messat(idateinf_messat))
     ELSE
        IF (date_s_annee_messat .GE. INT(t_clim_messat(idateinf_messat))) THEN
           interp= (date_s_annee_messat-t_clim_messat(idateinf_messat))/&
                  ((366*24*60*60+t_clim_messat(idatesup_messat))-t_clim_messat(idateinf_messat))
        ELSE
           interp=(date_s_annee_messat+366*24*60*60-t_clim_messat(idateinf_messat))/ &
                  ((366*24*60*60+t_clim_messat(idatesup_messat))-t_clim_messat(idateinf_messat))
        END IF
     END IF

   END IF  !endif on l_messat_clim

   IF (l_messat_obs) THEN
!       file_mes_obs_path=TRIM(filemessat_obs)//'/'//CHAR(ian)//'/'//    &
!                         CHAR(ian)//CHAR(imois)//CHAR(ijour)//'.nc'
      IF (t .GT. t_obs_messat(idatesup_messat)) THEN
          idateinf_messat=idatesup_messat
          IF (idatesup_messat .NE. idimt_messat_obs) THEN
              idatesup_messat=idatesup_messat+1
          ELSE
              MPI_master_only WRITE(*,*)'derniere date du fichier messat depassee'
              stop
          END IF
      END IF
!   MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',t
!   MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
!   MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', t_obs_messat(idateinf_messat),' et ',t_obs_messat(idatesup_messat)
      IF (idatesup_messat .GT. idateinf_messat)  THEN
        interp= (t -t_obs_messat(idateinf_messat))/        &
                 (t_obs_messat(idatesup_messat)-t_obs_messat(idateinf_messat))
      END IF
       
   END IF

!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast
#ifdef key_MARS
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#else
       DO i=ifirst,ilast
#endif
         IF (TOTAL_WATER_HEIGHT .LE. RESIDUAL_THICKNESS_WAT) THEN
            forcSPM(i,j)=valmanq
         ELSE

            ! calcul du forcage MES par donnees satellitales
            ! -----------------------------------------------
            IF (l_messat_obs) THEN
               forcSPM(i,j)=messat_obs(i,j,idateinf_messat)+                                &
               (messat_obs(i,j,idatesup_messat)-messat_obs(i,j,idateinf_messat))*interp

            END IF !endifon l_messat_obs

            IF (l_messat_clim) THEN
              forcSPM(i,j)=messat(i,j,idateinf_messat)+                                &
                 (messat(i,j,idatesup_messat)-messat(i,j,idateinf_messat))*interp

            END IF  !endif on l_messat_clim
         END IF  !endif on htot
        END DO
      END DO
!$OMP END DO
  
  ENDIF

   END SUBROUTINE BIOLink_SPMsat_file
#endif
   !!======================================================================

#endif

END MODULE coupleur_BIOLink

