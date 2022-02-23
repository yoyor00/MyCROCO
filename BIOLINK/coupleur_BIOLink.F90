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
   !&E ** Purpose : concerns coupling Biogeochemical ocean/coastal/sediment 
   !&E              model BIOLink with hydro code
   !&E Ocean dynamics Bio :  Initializations, reading of files *.dat 
   !&E                       (rivers, discharges, bio variables...)
   !&E                       + update sources and sinks terms
   !&E
   !&E ** Description :
   !&E     subroutine BIOLink_initialization   ! initialization of BIOLink  - 
   !&E                                          routine called by main
   !&E
   !&E     subroutine BIOLink_init             ! specifics initializations of 
   !&E                                          BIO modules  - routine called 
   !&E                                          by main 
   !&E                                                                        
   !&E     subroutine BIOLink_alloc            ! Allocates different tables 
   !&E                                           used by BIOLink 
   !&E     subroutine BIOLink_update           ! evaluation of sinks and sources 
   !&E                                          terms  - routine called by step
   !&E                
   !&E     subroutine BIOLink_convarray        ! conversion of 3D or 4D array from 
   !&E                                           hydro model to BIOLink - 
   !&E                                           called by BIOLink_update
   !&E
   !&E     subroutine BIOLink2hydro            !  conversion array BIOLink   
   !&E                                               to hydro host model 
   !&E
   !&E     subroutine BIOLink_updateconc_BIO   !  change tracers-substances 
   !&E                                            concentrations if needed
   !&E
   !&E   History :
   !&E    !  2019-08 (B. Thouvenin) issued from bloom and peptic for portability adaptation
   !&E
   !&E=========================================================================================

!*****************************************************************************************!
!                                                                                         !
!                                      Interface                                          !
!                                                                                         !
!*****************************************************************************************!

#include "coupleur_define_BIOLink.h" ! Equivalence of names between the hydro and the different
                                     ! Tracer and biological model. Also contains a function
                                     ! For the allocation of the main tables

#ifdef key_MARS
#include "coupleur_dimhydro_BIOLink.h"
   USE sflxatm,      ONLY : rad
#else
   USE module_BIOLink ! script that groups all the files of BIOLink together and allows the 
                      ! access to their functions/subroutines/variables
   USE comsubstance   ! Access to the functions/variables of SUBSTANCE ( from the MUSTANG
                      ! sediment model)
#endif /* key_MARS */

   USE comBIOLink     ! Common variables of the BIOLink coupler

   USE coupleur_BIOLink_helping ! For the helping functions of BIOLink

   USE coupleur_BIOLink_physics ! For the physics related functions of 
                                ! BIOLink
#if defined MUSTANG
   USE comMUSTANG ,  ONLY : htot ! Height of the water column
#endif /*MUSTANG*/

#if defined ECO3M   
   USE COUPLEUR_PHY_BIO ! Internal functions of coupling of ECO3M
#endif /* ECO3M */

 
   IMPLICIT NONE
  
!*****************************************************************************************!
!                                                                                         !
!                                  Accessibility                                          !
!                                                                                         !
!*****************************************************************************************!

     !====================================================================
     ! Internal BIOLink functions and routines
     !====================================================================

    PUBLIC BIOLink_initialization ! initialization of BIOLink  
                                  ! - called by main
 
    PUBLIC BIOLink_init           ! Initialization of some of BIOLink tables 
                                  ! - called by main
                                                                           
    PUBLIC BIOLink_alloc          ! Allocation of different tables
    PUBLIC BIOLink_update         ! Update of BIOLink and biology model variables
                                  ! - called by main

#ifdef key_MARS
    PUBLIC BIOLink_exchgMPI_cvwat ! I do not know
#endif


!*****************************************************************************************!
!                                                                                         !
!                         Routines and functions                                          !
!                                                                                         !
!*****************************************************************************************!

 CONTAINS

  !!======================================================================

  SUBROUTINE BIOLink_initialization(icall)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_initialization  ***
  !&E
  !&E ** Purpose : intialization of the biological/tracer model
  !&E
  !&E ** Description : Works in two calls. The first call, with the 
  !&E                  icall = 1 reads the parameter files of the different
  !&E                  models. The second calls use internal functions of 
  !&E                  those models to initialize their variables
  !&E
  !&E ** Called by : main.F90
  !&E
  !&E ** External calls : peptic_param,peptic_alloc_var 
  !&E                     from peptic_initdefine, 
  !&E                     bloom_param,bloom_init_iv 
  !&E                     from peptic_initdefine or 
  !&E                     meteor_param from meteor_initdefine
  !&E
  !&E ** Reference : I do not know
  !&E
  !&E ** History :
  !&E       !  2019-08  (B.Thouvenin) :  creation 
  !&E       !  2022-02  (G. Koenig)   : Commenting 
  !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined PEPTIC
  USE peptic_initdefine, ONLY : peptic_param,peptic_alloc_var
  ! Functions for the reading of parameter files and allocation of variables
  ! of the PEPTIC biological model
#elif defined BLOOM
  USE bloom_initdefine, ONLY : bloom_param,bloom_init_iv
  ! Function for the reading of parameter files and allocation of variables
  ! of the BLOOM biological model
#elif defined METeOR
  USE meteor_initdefine, ONLY : meteor_param
  ! Function for the reading of parameter files of METeOR contaminant model
#endif

     !====================================================================
     ! External arguments
     !====================================================================

  INTEGER, INTENT(IN) :: icall ! To differentiate the two calls to the function
                               ! The first call reads the parameter files
                               ! and the second initialize the tracer/biological
                               ! model 

     !====================================================================
     ! Local declarations of variables
     !====================================================================
  
  INTEGER               :: isubs
   
     !====================================================================
     ! Execution of the function
     !====================================================================
 
   IF (icall==0) THEN ! First call of the function, we read the parameter
                      ! files

#if ! defined key_MARS
#  if defined PEPTIC
    CALL peptic_param('r')

#  elif defined BLOOM
    CALL bloom_param('r')

#  elif defined METeOR
    CALL meteor_param('r')

#  endif /* PEPTIC/BLOOM/METeOR */
#endif /* key_MARS */

   ELSE ! Second call of the function, now we initialize the variable of the 
        ! Biology/Tracer models

#if ! defined key_MARS

     TIME_BEGIN=CURRENT_TIME ! Time begin is the beginning time for BIOLink and 
                             ! biological/tracer models, while current time
                             ! is the time from the hydro model

     CALL BIOLink_alloc() ! Allocation of the variables used by BLOOM :
                          ! Source and sink terms,
                          ! Height of the water column
                          ! and light attenuations related variables
       

! Initialization of the tracer/fixed variables in biological/tracer models
! with their names and characteristics

#  ifdef BLOOM 
     DO isubs=1,nv_adv ! Counter for tracer variables
       CALL bloom_init_iv(isubs,standard_name_var(isubs),1)
     END DO

     DO isubs=nv_adv+1,nv_adv+nv_fix !Counter for fixed variables
       CALL bloom_init_iv(isubs,standard_name_var_fix(isubs-nv_adv),2)
     END DO

#    ifdef key_benthic
     DO isubs=nv_adv+nv_fix+nv_bent ! Counter for benthic variables
       CALL bloom_init_iv(isubs,standard_name_var_bent(isubs-nv_adv-nv_fix),3)
     END DO
#    endif /* key_benthic */
#  endif /* BLOOM */

#  if defined PEPTIC
     CALL peptic_alloc_var ! For PEPTIC all the variables 
                           ! are declared in the same function
#  endif /* PEPTIC */

#  if defined ECO3M 
     CALL ALLOC_VAR_Eco3M ! 
     nbcallbio = -1 ! Internal variable of ECO3M
     CALL main_bio(0.0)
#  endif /* ECO3M */
  
  ! Writing of the biological parameters used in an external file
#  if defined PEPTIC
     CALL peptic_param('w')
#  elif defined BLOOM
     CALL bloom_param('w')
#  elif defined METeOR
     CALL meteor_param('w')
#  endif /* PEPTIC/BLOOM/METeOR */

#endif /* ! defined key_MARS */

  ! Allocation of diagnostic variables for BLOOM
#if defined BLOOM
     CALL BIOLink_read_vardiag
#endif /* BLOOM */

   ENDIF ! /* icall parameter */

   MPI_master_only PRINT*, 'END of BIOLink_initialization'


  END SUBROUTINE BIOLink_initialization

  SUBROUTINE BIOLink_init(ifirst,ilast,jfirst,jlast)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_init  ***
  !&E
  !&E ** Purpose : Initialization of tables of variables in biological 
  !&E              models and of table for the helping functions 
  !&E
  !&E ** Description : First it sets the time at the initial time of 
  !&E                  the hydro model + one slow time step, then it
  !&E                  initializes the settling velocities and the
  !&E                  conservation log file. Then it uses biological/
  !&E                  tracer routines to initialize the table of those
  !&E                  models and finally it initializes the tables of 
  !&E                  the helping functions.
  !&E
  !&E ** Called by : main
  !&E
  !&E ** External calls : bloom_userinit of bloom_initdefine,
  !&E                     bloom_wavefile_MANGAbio from bloom,
  !&E                     meteor_read_react from meteor_initdefine
  !&E
  !&E ** Reference : I do not know
  !&E
  !&E ** History :
  !&E       !  2015-07  : Creation by B. Thouvenin I guess 
  !&E       !  2022-02  : Commenting (G. Koenig)
  !&E
  !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined BLOOM
   USE bloom_initdefine, ONLY : bloom_userinit ! Initialization of the tables
                                               ! used to store variables of BLOOM
#  if defined key_MANGAbio && defined key_MANGAbiovague
   USE bloom,            ONLY : bloom_wavefile_MANGAbio ! Reading of a file
                                                        ! With orbital velocity of 
                                                        ! Waves for the Manche/Golfe de 
                                                        ! Gascogne model
#  endif /* key_MANGAbio && key_MANGAbiovague */
#endif /* BLOOM */

#if defined METeOR
   USE meteor_initdefine, ONLY : meteor_read_react      ! I do not know
#endif /* METeOR */

  IMPLICIT NONE

     !====================================================================
     ! External arguments
     !====================================================================


  !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst,ilast,jfirst,jlast ! Horizontal indexes for 
                                                     ! the navigation in the 
                                                     ! table.
     !====================================================================
     ! Local declarations of variables
     !====================================================================
                                                              
  INTEGER               :: i,j,k,iv ! Counter variables
  CHARACTER(LEN=lchain) :: logfilename ! Name of the file to write the log

     !====================================================================
     ! Execution of the function
     !====================================================================

     !************** Determination of the time steps *********************!
     
#ifdef key_MARS
   t_bio=MAX(CURRENT_TIME,tdebsubs) 
#else
   t_bio=CURRENT_TIME+TRANSPORT_TIME_STEP ! The biology time steps is updated 
                                          ! From the hydrodynamical model 
                                          ! With a transport time step
#endif /* key_MARS */

     !********************* Sinking velocity *****************************!

#if ! defined MUSTANG

     ! Initialization of the sinking rate for all the tracers at all depth and position
     ! The initialization is taken at the mean between the maximum and minimum sinking depth
     DO j=jfirst,jlast
      DO i=ifirst,ilast
        DO iv=1,nvp
          DO k=1,NB_LAYER_WAT
              WS_BIOLink(k,iv,i,j)=(ws_free_min(iv)+ws_free_max(iv))/2.0_rsh
            ENDDO
          ENDDO
        END DO
      END DO

#endif /* MUSTANG */

     !****************** Verification of conservation *********************!

#if defined BIOLink_verif_conserv 
     
      logfilename='conservBIOLink_'//TRIM(suffix_fileres)//'.log' ! Name
                                                                  ! for cons
                                                                  ! -servati
                                                                  ! on file
      IF (iscreenlog_conserv/=6) OPEN(unit=iscreenlog_conserv,file=logfilename)
      ! Opening of the log file to test the conservativity
      ! I think iscreenlog_conserv is initialized at 66 and never changed
      ! So this should always happen

#endif /* BIOLink_verif_conserv */

#if defined BLOOM

      CALL bloom_userinit(ifirst,ilast,jfirst,jlast) ! Initialization of tab
                                                     ! -les in the biological
                                                     ! /tracer models
#elif defined METeOR

      CALL meteor_read_react ! reading reactions between contaminant species

#endif /* BLOOM/METeOR */

     !****************** Wave orbital velocities *********************!

#if defined key_MANGAbio && defined key_MANGAbiovague
   
      ubr_vague(:,:,:)=0.0_rsh ! Initialization of 
      vbr_vague(:,:,:)=0.0_rsh ! wave orbital 
      hs_vague(:,:,:)=0.0_rsh  ! velocities and height

   DO j=jfirst,jlast
#  ifdef key_MARS

     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)

#  else

     DO i=ifirst,ilast

#  endif /* key_MARS */

       IF (BATHY_H0(i,j) .EQ. -valmanq) THEN
         ubr_vague(i,j,:)=valmanq ! Initialization to NaN value 
         vbr_vague(i,j,:)=valmanq ! in case of land
         hs_vague(i,j,:)=valmanq  ! cell
       END IF
     END DO
   END DO
   
   CALL bloom_wavefile_MANGAbio(ifirst,ilast,jfirst,jlast,0) ! Reading
                                                             ! wave climato
                                                             ! -logy file
#endif /* key_MANGAbio && key_MANGAbiovague */

     !*********** Satellital measurement of suspended matter ****************!

#if defined key_messat
   ! initialization of the suspended matter measured by satellite
   messat(:,:,:)=0.0_rsh
   DO j=jfirst,jlast
#  ifdef key_MARS
     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
#  else
     DO i=ifirst,ilast
#  endif /* key_MARS */
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

   CALL BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,0) ! Reading of sus-
                                                         ! pended matter
                                                         ! file
#endif /* key_messat */

  END SUBROUTINE BIOLink_init

  !!======================================================================
  SUBROUTINE BIOLink_alloc()

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_alloc ***
  !&E
  !&E ** Purpose : allocation of common variables tables. Those variables 
  !&E              were declared earlier but their dimensions depends on 
  !&E              the reading of some parametrization files that are read  
  !&E              later in the CROCO model. Once those files are read, we 
  !&E              can allocate the tables their final shapes.
  !&E                            
  !&E
  !&E ** Description : The entire function consist of a succession of 
  !&E                  memory allocation followed by value initialization 
  !&E                  of the tables to dummy values.
  !&E
  !&E ** Called by : coupleur_BIOLink in the second call of the subroutine 
  !&E                   "BIOLink_initialization"
  !&E
  !&E ** External calls :  All the tables are declared as public and 
  !&E                      therefore are accessible even if they are not
  !&E                      given as argument of the subroutine. 
  !&E
  !&E ** Reference : None
  !&E
  !&E ** History : Not sure. Probably created by B. Thouvenin.
  !&E
  !&E---------------------------------------------------------------------      

     !*************************************************************************!
     !*************************************************************************!
     !****************** Variables for the biological or **********************!
     !*********************** chemical model **********************************!
     !*************************************************************************!
     !*************************************************************************!

     !===================================================================
     ! Table of positive concentrations for tracer variables
     !===================================================================

#  if ! defined ECO3M  
      ALLOCATE(WATCONCPOS(nv_adv,NB_LAYER_WAT,PROC_IN_ARRAY))
      WATCONCPOS(:,:,:,:)=0.0_rsh
#  endif /* ECO3M */

     !===================================================================
     ! Table of positive concentrations for fixed variables
     !===================================================================

      ALLOCATE(FIXCONCPOS(nv_fix,NB_LAYER_WAT,PROC_IN_ARRAY))
      FIXCONCPOS(:,:,:,:)=0.0_rsh

     !===================================================================
     ! Table of positive concentrations for benthic variables
     !===================================================================

#  ifdef key_benthic
      ALLOCATE(BENTCONCPOS(nv_bent,PROC_IN_ARRAY))
      BENTCONCPOS(:,:,:,:)=0.0_rsh
#  endif /* key_benthic */ 

     !===================================================================
     ! Table of source and sink terms for tracer variables
     !===================================================================

      ALLOCATE( BIO_SINKSOURCES(ARRAY_SINKSOURCES))
      BIO_SINKSOURCES(:,:,:,:)=0.0_rsh

     !===================================================================
     ! Table of sink and source terms for fixed variables
     !===================================================================

      ALLOCATE(BIO_SKSC_FIX (ARRAY_FIXED_SKSC))
      BIO_SKSC_FIX (:,:,:,:)=0.0_rsh

     !*************************************************************************!
     !*************************************************************************!
     !******************** Variables from the hydro model *********************!
     !*************************************************************************!
     !*************************************************************************! 

     !====================================================================
     ! Height of the water column tables
     !====================================================================

#  if ! defined MUSTANG
        ALLOCATE( TOTAL_WATER_HEIGHT(PROC_IN_ARRAY_m2p2))
#  endif /* MUSTANG */
        
        ALLOCATE( THICKLAYERWC(NB_LAYER_WAT,PROC_IN_ARRAY))
        THICKLAYERWC(:,:,:)=0.0_rsh
        
        ALLOCATE( THICKLAYERWW(NB_LAYER_WAT,PROC_IN_ARRAY))
        THICKLAYERWW(:,:,:)=0.0_rsh

     !=====================================================================
     !  Temperature and salinity variables
     !=====================================================================

      ALLOCATE( TEMP_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
      TEMP_BIOLink(:,:,:)=0.0_rsh

      ALLOCATE( SAL_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
      SAL_BIOLink(:,:,:)=0.0_rsh

     !=====================================================================
     !  Bottom current variables
     !=====================================================================

#  if defined BLOOM && defined key_benthos
      ALLOCATE(BOTTCURRENTBIO(PROC_IN_ARRAY)) 
#  endif /* BLOOM && key_benthos */


     !*************************************************************************!
     !*************************************************************************!
     !*************** Tables for the model helping routines ****************!
     !********************* Or sediment model *********************************!
     !*************************************************************************!
     !*************************************************************************!

     !=====================================================================
     !  Sinking velocity of tracers tables
     !=====================================================================

      ALLOCATE( WS_BIOLink(NB_LAYER_WAT,nv_adv,PROC_IN_ARRAY))
      WS_BIOLink(:,:,:,:)=0.0_rsh

     !=====================================================================
     !  Radiation tables
     !=====================================================================

     ! Photosynthetic available radiations

#  if defined  BIOLink_PAR_eval
      ALLOCATE( PAR_top_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )
      PAR_top_layer(:,:,:)=0.0_rsh

#    if defined PEPTIC
      ALLOCATE( PAR_avg_layer_phyto(1,NB_LAYER_WAT,PROC_IN_ARRAY) )
      PAR_top_layer_day(:,:,:)=0.0_rsh
      ALLOCATE( PAR_top_layer_day(NB_LAYER_WAT,PROC_IN_ARRAY))
      PAR_avg_layer_phyto(1,:,:,:)=0.0_rsh

#    elif defined METeOR
      ALLOCATE( Flimrad_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )

#    endif /* PEPTIC/METEOR */
#  endif /* BIOLink_PAR_eval */

     ! Extinction and attenuation

#  if defined  BIOLink_PAR_eval
      ALLOCATE( EXTINCTION_RAD(NB_LAYER_WAT,PROC_IN_ARRAY) )
      EXTINCTION_RAD(:,:,:)=0.0_rsh
#  endif /* BIOLink_PAR_eval */

#  if defined BLOOM
      ALLOCATE(extinction_ave4d(NB_LAYER_WAT,PROC_IN_ARRAY))
      extinction_ave4d(:,:,:)=0.0_rsh
      ALLOCATE(extinction_aveh(NB_LAYER_WAT,PROC_IN_ARRAY))
      extinction_aveh(:,:,:)=0.0_rsh
      ALLOCATE(extinction_tab(4,24,NB_LAYER_WAT,PROC_IN_ARRAY))
      extinction_tab(:,:,:,:,:)=0.0_rsh
      t_cum_extinctionh=0.0_rsh
      ihour_previous=0

#  endif /* BLOOM */

     ! Variables that affect the extinction
 
      ALLOCATE( SPMTOT_MGL(NB_LAYER_WAT,PROC_IN_ARRAY))
      SPMTOT_MGL(:,:,:)=0.0_rsh

#  if defined  BIOLink_PAR_eval
      ALLOCATE( BIOLink_chloro(NB_LAYER_WAT,PROC_IN_ARRAY) )  
      BIOLink_chloro(:,:,:)=0.0_rsh
#  endif /* BIOLink_PAR_eval */

#  if defined BLOOM

      ALLOCATE( effetlumiere_day_diat(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_diat(:,:,:)=0.0_rsh
      ALLOCATE( effetlumiere_day_dino(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_dino(:,:,:)=0.0_rsh
      ALLOCATE( effetlumiere_day_nano(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_nano(:,:,:)=0.0_rsh

#    ifdef key_karenia
      ALLOCATE( effetlumiere_day_karenia(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_karenia(:,:,:)=0.0_rsh
#    endif /* key_karenia */

#    ifdef key_psnz
      ALLOCATE( effetlumiere_day_psnz(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_psnz(:,:,:)=0.0_rsh
#    endif /* key_psnz */

#    ifdef key_phaeocystis
      ALLOCATE( effetlumiere_day_phaeocystis(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_phaeocystis(:,:,:)=0.0_rsh
#    endif /* key_phaeocystis */

#  endif /* BLOOM */

     !=====================================================================
     !  Wind variables
     !=====================================================================

#  if ! defined BULK_FLUX
      ALLOCATE( WIND_SPEED(PROC_IN_ARRAY))
      WIND_SPEED(:,:)=0.0_rsh
#  endif /* BULK_FLUX */

     !=====================================================================
     !  Estimation of the suspended matter via satellite
     !=====================================================================

#  ifdef key_messat
      ALLOCATE(messat(PROC_IN_ARRAY,51))
      messat(:,:,:)=0.0_rsh
#  endif /* key_messat */

     !=====================================================================
     !  Oyster related variables
     !=====================================================================

#  if defined BLOOM
#    ifdef key_oyster_SFG
      ALLOCATE(tpostpontecoq(PROC_IN_ARRAY))
      ALLOCATE(tpostpontegam(PROC_IN_ARRAY))
#    endif /* key_oyster_SFG */

#    if defined key_oyster_benthos || defined key_oyster_DEB || defined key_oyster_SFG
      ALLOCATE(nbhuitre(PROC_IN_ARRAY))
      hautable(:,:)=0.0_rsh
      ALLOCATE(hautable(PROC_IN_ARRAY))
      nbhuitre(:,:)=0.0_rsh
#    endif /* key_oyster_benthos */
#  endif /* BLOOM */


END SUBROUTINE  BIOLink_alloc

  SUBROUTINE BIOLink_update(ifirst,ilast,jfirst,jlast   &
#if defined key_MARS && (defined key_oyster_SFG || defined key_oyster_DEB)
          , CELL_SURF                                                &
#endif
         )
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_update  ***
  !&E
  !&E ** Purpose : Update of the BIOLink concentration, sources and sinks
  !&E              and helping variables (settling velocities, number of 
  !&E              oyster, etc...) at each time step.
  !&E
  !&E ** Description : 
  !&E
  !&E ** Called by : init
  !&E
  !&E ** External calls : peptic_sksc_wat, peptic_SPMtot_Chla from peptic
  !&E                     bloom_sksc_wat, bloom_eval_diag2d, 
  !&E                     bloom_SPMtot_Chla and bloom_extinction_avg
  !&E                     from bloom 
  !&E                     bloom_wavefile_MANGAbio from bloom
  !&E                     COUPLEUR_PHY_BIO from ECO3M
  !&E                     meteor_sksc_wat, meteor_reac_equi from METeOR
  !&E                     meteor_sksc_wat, meteor_reac_equi from comvars2d
  !&E
  !&E ** Reference : There is first an evaluation to know if the time step 
  !&E                the hydrodynamic model has evolved enough. If yes,
  !&E                there is an update of the biological time, concentration
  !&E                , height of the water column and helping variables. 
  !&E                And finally sometimes we reinitialize the oysters.
  !&E
  !&E ** History :
  !&E       !  2019-08  (B.Thouvenin)  creation 
  !&E       !  2022-02  (G. Koenig) Commenting
  !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined PEPTIC
  USE peptic,            ONLY : peptic_sksc_wat, peptic_SPMtot_Chla
                                ! Computation of source and sink terms in 
                                ! the water column and the total concentration
                                ! Of suspended matter and chlorophyll A
#elif defined BLOOM
  USE bloom,            ONLY : bloom_sksc_wat,bloom_eval_diag2d, bloom_SPMtot_Chla, &
                               bloom_extinction_avg
                               ! Computation of source and sink terms, diagnostic 2D
                               ! Concentration of total suspended particulate matter 
                               ! and average extinction of light

#  if defined key_MANGAbio && defined key_MANGAbiovague
  USE bloom,            ONLY : bloom_wavefile_MANGAbio
                               ! Reading of the orbital velocity of waves in an external
                               ! file
#  endif /* key_MANGAbio && key_MANGAbiovague */

#elif defined ECO3M
  USE COUPLEUR_PHY_BIO   ! Internal coupleur of ECO3M

#elif defined METeOR
  USE meteor,          ONLY : meteor_sksc_wat,meteor_reac_equi
                              ! Sources and sink terms and
                              ! Equilibrium reactions 
#endif /* PEPTIC/BLOOM/ECO3M/METeOR */

#if defined key_MARS 
#  if defined key_turbclim && defined key_daily_climato_kpar
    USE comvars2d,     ONLY : mes_sat ! I have no idea
#  endif /* key_turbclim && key_daily */  
#endif /* key_MARS */     

     !====================================================================
     ! External arguments
     !====================================================================

   INTEGER, INTENT(IN) :: ifirst,ilast,jfirst,jlast
#if defined key_MARS && (defined key_oyster_SFG || defined key_oyster_DEB)
   REAL(KIND=rsh),DIMENSION(ARRAY_CELL_SURF),INTENT(IN)          :: CELL_SURF
#endif /* key_MARS && (key_oyster_SFG || key_oyster_DEB ) */
   
     !====================================================================
     ! Local declarations of variables
     !====================================================================

   REAL(KIND=rlg)         :: dt_bio_cor ! Timestep for biology models
   INTEGER                :: i,j,k,iv ! Counter variables
   INTEGER                :: kmaxmod,itend ! Index of termination for
                                                       ! vertical and time loops
   INTEGER                :: tool_julien ! Function to determine the julian day
   CHARACTER(LEN=19)      :: tool_sectodat,cdate ! Function to change the second to date
   REAL(KIND=rsh), DIMENSION(PROC_IN_ARRAY)   :: forcSPM ! Something related to 
                                                         ! suspended particulate matter
#if defined key_MANGAbio && defined key_MANGAbiovague
   REAL(KIND=rsh), DIMENSION(NB_LAYER_WAT,PROC_IN_ARRAY)  :: forcSPMk ! Something related
                                                                      ! to suspended particulate
                                                                      ! Matter with a vertical
                                                                      ! Dimension
#endif 

     !====================================================================
     ! Execution of the function
     !====================================================================


     !************** Determination of the time steps *********************!

      itend=0 ! Initialization of the end time counter
      IF(CURRENT_TIME .GE. t_bio) THEN ! Updating of variables
      
        dt_bio_cor=MAX(ECO_TIME_STEP,TRANSPORT_TIME_STEP) ! Updating of
        BIO_TIME_STEP=dt_bio_cor+(CURRENT_TIME-t_bio) ! time steps

        cdate = tool_sectodat(CURRENT_TIME) ! Conversion of time to date 
        CALL tool_decompdate(cdate,ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK, &
                             iheure_BIOLINK,iminu_BIOLINK,isec_BIOLINK) 
        
        jjulien_BIOLINK=tool_julien(ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK)&
                                   - tool_julien(1,1,ian_BIOLINK)+1

        !*********** Computation of physical variables ********************!

        CALL BIOLink_water_column(ifirst,ilast,jfirst,jlast) ! Height of the
                                                             ! water column
      
        CALL BIOLink_convarray(ifirst,ilast,jfirst,jlast     & ! Conversion
                                                               ! of arrays 
                                                               ! in BIOLink
                                                               ! order

#if defined key_nosubstmodule

                  ,WAT_SETTL                                     &

#endif /* key_nosubstmodule */

#ifdef key_MARS

                  ,TEMPERATURE_MOD,SALINITY_MOD                  &

#endif /* key_MARS */
                  )     
      CALL BIOLink_sinking_rate(ifirst,ilast,jfirst,jlast) ! Evaluation and
                                                           ! limitation of 
                                                           ! sinking rates

      !*********************** Blue Öyster Cult ***************************!

      ! Reinitialization of the number of oyster for the 1rst january of each year
#if defined BLOOM && defined key_oyster_DEB
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
#endif /* BLOOM && key_oyster_DEB */

      !*********************** Estimation of PAR ***************************!

#if defined BIOLink_PAR_eval
!$OMP SINGLE
      forcSPM(:,:)=0.0_rsh ! Initialization of the SPM forcing

#  if defined key_messat
      CALL BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,1,forcSPM=forcSPM)
      ! SPM read from satellite input file
#  else
      ! SPM read from climato file
#    if defined key_MARS

#      if defined key_turbclim && defined key_daily_climato_kpar

      forcSPM(:,:)=mes_sat(:,:) ! Reading of the MES by the hydro code

#      endif /* key_turbclim && key_daily_climato_kpar */

#    endif /* key_MARS */

#  endif /* key_messat */

#  if defined key_MANGAbio && defined key_MANGAbiovague
      
      CALL bloom_wavefile_MANGAbio(ifirst,ilast,jfirst,          & ! Reading
                                   jlast,1,forcSPMk=forcSPMk) ! of waves fro
                                                              ! m an externa
                                                              ! file
#  endif /* key_MANGAbio && key_MANGAbiovague */
!$OMP END SINGLE

#  if defined PEPTIC    

      CALL peptic_SPMtot_Chla(ifirst,ilast,jfirst,jlast) ! Impact of 
                                                         ! SPM on Chla
                                                         ! from the PEPTIC
                                                         ! model

#  elif defined BLOOM

      CALL bloom_SPMtot_Chla(ifirst,ilast,jfirst,jlast          & ! Impact 
                                                         ! of SPM on Chla
                                                         ! from the BLOOM
                                                         ! model

#     if defined key_MANGAbio && defined key_MANGAbiovague

                             ,forcSPMk                          &

#     endif /* key_MANGAbio && key_MANGAbiovague */

#     if defined key_messat

                              ,forcSPM                          &

#     endif /* key_messat */

                         )

#  endif /* BLOOM */  

      CALL BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate) ! Computation 
                                                             ! of solar ra
                                                             ! -diation and
                                                             ! extinction

#  if defined BLOOM
      
      CALL bloom_extinction_avg(ifirst,ilast,jfirst,jlast) ! Computation of
                                                           ! 4 days extinct
                                                           ! -ion average in
                                                           ! BLOOM

#  endif /* BLOOM */

#  if defined PEPTIC && defined key_rand_extinc
   
      CALL peptic_abondance(ifirst,ilast,jfirst,jlast) ! I do not know

#  endif /* PEPTIC && key_rand_extinc */

#endif /* BIOLink_PAR_eval */

      !********************* Diagnostic variables **************************!

#if defined BLOOM

      CALL  bloom_eval_diag2d(ifirst,ilast,jfirst,jlast) ! Bloom function
                                                         ! for diagnostic
                                                         ! evaluations

#endif /* BLOOM */

      !**************** Wind speed for some BLOOM keys *********************!

#if defined key_oxygen || defined key_zostera || defined METeOR

#  if ! defined key_MARS && ! defined BULK_FLUX

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

#  endif /* key_MARS && BULK_FLUX */

#endif /* defined key_oxygen || defined key_zostera || defined METeOR */

 
      !******************* Tracers sources and sinks ***********************!
    
#if defined PEPTIC

     CALL peptic_sksc_wat(ifirst,ilast,jfirst,jlast) ! Source and sink 
                                                     ! terms for PEPTIC

#elif defined BLOOM

     CALL bloom_sksc_wat(ifirst,ilast,jfirst,jlast               & ! Source 
                                                                   ! and sin
                                                                   ! k terms
                                                                   ! for 
                                                                   ! BLOOM

#  if defined GLS_MIXING

                          ,CIN_TURBULENT_ENERGY                                    &

#  endif /* GLS_MIXING */

#  if defined key_zostera || defined key_oxygen

                          ,WIND_SPEED                                              &

#  endif /* key_zostera || key_oxygen */

#  if defined key_MARS && (defined key_oyster_SFG || defined key_oyster_DEB)

                       ,CELL_SURF                                            &

#  endif /* key_MARS && ( key_oyster_SFG || key_oyster_DEB) */

                                 ) 
#elif defined METeOR

    IF ( l_treat_re_dc) THEN
        CALL meteor_sksc_wat(ifirst,ilast,jfirst,jlast,WIND_SPEED)
    ENDIF                      

#endif /* BLOOM/METeOR */
   
      !**************** Conversion of the order of arrays  *****************!

#if ! defined key_MARS 

      CALL BIOLink2hydro(ifirst,ilast,jfirst,jlast)     

#endif /* key_MARS */

      !*********************** Final time stepping *************************!
    
!$OMP SINGLE
     t_bio=t_bio+BIO_TIME_STEP
!$OMP END SINGLE
     itend=1
 
   ENDIF  ! t>t_bio

     !****************** Update of tracer concentration ********************!
     !****************** if not computed by hydro model ********************!

#if defined BIOLink_UPDATE_CONCBIO 

#  if defined METeOR && ! defined BLOOM && ! defined PEPTIC && ! defined ECO3M

    IF ( l_treat_re_dc) THEN

#  endif /* METeOR && BLOOM && PEPTIC && ECO3M */   

      CALL BIOLink_updateconc_BIO(ifirst,ilast,jfirst,jlast)

#  if defined METeOR && ! defined BLOOM && ! defined PEPTIC && ! defined ECO3M

    ENDIF 

#  endif /* METeOR && BLOOM && PEPTIC && ECO3M */

#endif /* BIOLink_UPDATE_CONCBIO */


#if defined METeOR

       ! Update of the concentration at each time step for 
       ! Rapid chemical reactions

       IF ( nreacmax > 0 .AND. l_treat_re_eq) THEN

        IF (itend == 0)THEN

         cdate = tool_sectodat(CURRENT_TIME)

         CALL BIOLink_water_column(ifirst,ilast,jfirst,jlast)
         CALL BIOLink_convarray(ifirst,ilast,jfirst,jlast)

#  if defined BIOLink_PAR_eval
!$OMP SINGLE
         forcSPM(:,:)=0.0_rsh
!$OMP END SINGLE
         CALL BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate)
#  endif /* BIOLink_PAR_eval */
        ENDIF
       
        CALL meteor_reac_equi(ifirst,ilast,jfirst,jlast,WIND_SPEED)
 
       ENDIF

#endif /* METeOR */

      !***************** Evolution of fixed variables **********************!

     IF (nv_fix > 0 ) THEN
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)
       DO j=jfirst,jlast
         DO i=ifirst,ilast

            kmaxmod=NB_LAYER_WAT ! The fixed variables are only computed 
                                 ! where MUSTANG is not used

            DO k=1,kmaxmod

              FIXED_VAR_CONC(FIXED_VAR_INDEXkij)=FIXED_VAR_CONC(FIXED_VAR_INDEXkij) &
                                       +TRANSPORT_TIME_STEP*BIO_SKSC_FIX(FIXED_SKSC_INDEXkij)

            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
     ENDIF

  END SUBROUTINE BIOLink_update

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
  !&E ** Description : Everything is said in the purpose
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : None
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : Creation by B. Thouvenin ( date unknown)
  !&E              Commenting and editing by G. Koenig (February 2022)
  !&E
  !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined key_MARS
   USE toolgeom,     ONLY : f_dzu,f_dzw
   USE comvars2d,    ONLY : ig,id,jb,jh,hm
#endif /* key_MARS */
 
     !====================================================================
     ! External arguments
     !====================================================================

   INTEGER, INTENT(IN) :: ifirst,ilast,jfirst,jlast

     !====================================================================
     ! Local declarations of variables
     !====================================================================

   INTEGER :: i,j,k,kmaxmod,iv

     !====================================================================
     ! Execution of the function
     !====================================================================


#if defined PEPTIC || defined BLOOM 
! PEPTIC and BLOOM are in the order : iv (tracer), k (vertical),
! i (zonal) and j (meridional)
! but CROCO is in the order : iv (tracer), i (zonal), j (meridional)
! and k (vertical). 

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,iv)

     DO j=jfirst,jlast
     
#  if defined key_MARS

       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)

#  else

       DO i=ifirst,ilast

#  endif /* key_MARS */

         kmaxmod=NB_LAYER_WAT

         WATCONCPOS(:,:,i,j)=0.0_rsh
         FIXCONCPOS(:,:,i,j)=0.0_rsh

#  if defined key_benthic

         BENTCONCPOS(:,i,j)=0.0_rsh

#  endif /* key_benthic */
 
         DO k=1,kmaxmod

#  if ! defined key_MARS
            
           TEMP_BIOLink(k,i,j)=TEMPHYDRO_ijk ! Inverting of the place of
           SAL_BIOLink(k,i,j)=SALHYDRO_ijk   ! the vertical index

#  endif /* key_MARS */

#  if defined BLOOM && defined key_benthos
   
           bottom_current(i,j)=BOTTOM_CURRENT_ij

#  endif /* BLOOM && key_benthos */

           DO iv=1,nv_adv ! Advected variables
           
             WATCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCADV_ivkij+ABS(WAT_CONCADV_ivkij))
            
           END DO ! nv_adv
           
           DO iv=1,nv_fix ! Fixed variables
             
             FIXCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCFIX_ifixkij+ABS(WAT_CONCFIX_ifixkij))
           
           END DO ! nv_fix
         END DO ! kmax_mod

#  if defined key_benthic

         DO iv=1,nspb

           BENTCONCPOS(iv,i,j)=0.5_rsh*(WAT_CONCBENT_ivij+ABS(WAT_CONCBENT_ivij))
           
         END DO ! nspb
#  endif /* key_benthic */

#  if defined key_oyster_DEB
         ! Reinitialization of oysters for the 1rst of each year
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
#  endif /* key_oyster_DEB */

       END DO ! i
     END DO ! j
!$OMP END DO

#  if ! defined key_MARS
!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast   

       DO i=ifirst,ilast

        DO iv=1,nv_adv

          DO k=1,NB_LAYER_WAT

              WS_BIOLink(k,iv,i,j)=WAT_SETTL_ivkij ! Conversion of the 
                                                   ! index order 
            ENDDO ! k

          ENDDO ! iv

        END DO ! i

      END DO ! j
!$OMP END DO
#  endif /* key_MARS */

#elif defined METeOR && ! defined ECO3M

! METeOR is in the order : iv (tracer), k (vertical),
! i (zonal) and j (meridional)
! but CROCO is in the order : iv (tracer), i (zonal), j (meridional)
! and k (vertical). 

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,iv)
     DO j=jfirst,jlast

       DO i=ifirst,ilast

          kmaxmod=NB_LAYER_WAT

          WATCONCPOS(:,:,i,j)=0.0_rsh
          FIXCONCPOS(:,:,i,j)=0.0_rsh

          DO k=1,kmaxmod

            TEMP_BIOLink(k,i,j)=TEMPHYDRO_ijk
            SAL_BIOLink(k,i,j)=SALHYDRO_ijk

            DO iv=1,nv_adv

              WATCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCADV_ivkij+ABS(WAT_CONCADV_ivkij))

            END DO ! iv

            DO iv=1,nv_fix

              FIXCONCPOS(iv,k,i,j)=0.5_rsh*(WAT_CONCFIX_ifixkij+ABS(WAT_CONCFIX_ifixkij))

            END DO ! iv

          END DO ! k

        END DO ! i

      END DO ! j
!$OMP END DO

#elif defined ECO3M

#endif /* BLOOM/PEPTIC/ECO3M */

   END SUBROUTINE  BIOLink_convarray



    !!============================================================================== 
 


#if ! defined key_MARS 

  SUBROUTINE BIOLink2hydro(ifirst,ilast,jfirst,jlast  &

#  if defined key_nosubstmodule

                               ,WAT_SETTL  &

#  endif
                                 )     

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink2hydro ***
  !&E
  !&E ** Purpose : Conversion of setling velocities array from BIOLink 
  !&E              to hydro model 
  !&E
  !&E ** Description : The array in BIOLink are first indexed with the depth,
  !&E                  while the ones of hydro models are first indexed by
  !&E                  horizontal positions. Here we convert the array of 
  !&E                  BIOLink to the hydro model format
  !&E
  !&E ** Called by :  BIOLink_update
  !&E
  !&E ** External calls : None
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : ! Created by B. Thouvenin 
  !&E              ! Commented by G. Koenig ( february 2022)
  !&E
  !&E---------------------------------------------------------------------


     !====================================================================
     ! External arguments
     !====================================================================
 
   INTEGER, INTENT(IN)    :: ifirst,ilast,jfirst,jlast ! Limits of the MPI
                                                       ! subdomain

#  if defined key_nosubstmodule

   REAL(KIND=rsh),DIMENSION(ARRAY_WAT_SETTL), INTENT(INOUT)   :: WAT_SETTL
                                             ! Array for storing the sett-
                                             ! ling velocities
#  endif

     !====================================================================
     ! Local declarations of variables
     !====================================================================

    INTEGER                  :: i,j,k,iv ! Spatial and tracer counters
    INTEGER                  :: i1,i2,i3,i4 ! Internal BIOLink counters
    REAL(KIND=rsh), DIMENSION(ARRAY_WATER_CONC0) :: xnegtr ! Array for
                                                           ! negative
                                                           ! tracers
                                                           ! concentrations
    REAL(KIND=rsh)           :: ztra ! I do not know
    
     !====================================================================
     ! Execution of the function
     !====================================================================

 !****************** Conversion for the index order of *********************!
 !***************** BIOLink to the one of the hydro model ******************!


!$OMP DO SCHEDULE(RUNTIME)

     DO j=jfirst,jlast

       DO i=ifirst,ilast

          DO k=1,NB_LAYER_WAT

            DO iv=1,nv_adv

              WAT_SETTL_ivkij=WS_BIOLink(k,iv,i,j)

            END DO

          END DO

        END DO

      END DO

!$OMP END DO

END SUBROUTINE  BIOLink2hydro

#endif /* key_MARS */

#if defined BIOLink_UPDATE_CONCBIO 

  SUBROUTINE BIOLink_updateconc_BIO(ifirst,ilast,jfirst,jlast)     

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_updateconc_BIO ***
  !&E
  !&E ** Purpose : Updating tracers concentration
  !&E
  !&E ** Description : First there is a loop to determine if the conce
  !&E                  ntrations are going to become negative. If is the 
  !&E                  case we limit the sources and sink terms. And then
  !&E                  we compute the advective flux of tracers in the cell
  !&E                  with the hydro time steps.
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : ! Created by B. Thouvenin
  !&E              ! Commented by G. Koenig ( february 2022 )
  !&E
  !&E---------------------------------------------------------------------

     !====================================================================
     ! External arguments
     !====================================================================
 
   INTEGER, INTENT(IN)  :: ifirst,ilast,jfirst,jlast ! Limits of the MPI
                                                     ! subdomain

     !====================================================================
     ! Local declarations of variables
     !====================================================================

    INTEGER                  :: i,j,k,iv ! Spatial and tracer counters
    INTEGER                  :: i1,i2,i3,i4 ! Internal BIOLink counters
    REAL(KIND=rsh), DIMENSION(ARRAY_WATER_CONC0) :: xnegtr ! Variable to 
                                                           ! limit the 
                                                           ! flux out
                                                           ! or in the 
                                                           ! cell
    REAL(KIND=rsh)           :: ztra ! Ratio between the concentration
                                     ! and the flux
    REAL(KIND=rsh)           :: sinsksourcesdt ! Sink and source terms 
                                               ! multiplied by timestep                                               
    
     !====================================================================
     ! Execution of the function
     !====================================================================

      !***************** We test if the tracer concentration ***************!
      !************************* becomes negative **************************!

         xnegtr(:,:,:) = 1.e+0 ! 

!$OMP DO SCHEDULE(RUNTIME)

         DO i1=IRANGE1

           DO i2=IRANGE2

             DO i3=IRANGE3

               DO i4=IRANGE4

                  sinsksourcesdt=BIO_SKSC_ADV(BIOSKSC_INDEX_EQ)  & ! Total
                                 *TRANSPORT_TIME_STEP ! flux of tracer in 
                                                      ! the cell in one
                                                      ! time step
                  
                  ! Test if the flux makes the concentration go negative
                  IF( (WATER_CONCENTRATION(WATCONC_INDEX_EQ)+ sinsksourcesdt ) < 0.0_rsh ) THEN

                     ztra = ABS( ( WATER_CONCENTRATION(WATCONC_INDEX_EQ) &
                                   - epsilon_BIOLink ) &
                                / ( sinsksourcesdt + epsilon_BIOLink ) )
                            ! Computation of the absolute value of the
                            ! ratio between the concentration and the 
                            ! flux with a precision epsilon_BIOLink

                     xnegtr(i4,i3,i2) = MIN( xnegtr(i4,i3,i2), ztra )
                            ! Limiting test. If the ratio ztra is
                            ! above 1 we reduce it so that it does not
                            ! make the concentration goes to 0.
                            ! This being said it also filters the 
                            ! addition to the cell.

                  ENDIF

               END DO

             END DO

           END DO

         END DO
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)

         DO i1=IRANGE1

           DO i2=IRANGE2

             DO i3=IRANGE3

               DO i4=IRANGE4  

                 WATER_CONCENTRATION(WATCONC_INDEX_EQ) = WATER_CONCENTRATION(WATCONC_INDEX_EQ)    &
                                                         + xnegtr(i4,i3,i2) &
                                                         * BIO_SKSC_ADV(BIOSKSC_INDEX_EQ) &
                                                         * TRANSPORT_TIME_STEP    
                 ! Computation of the flux from the advection terms 
                 ! with the transport time step
               END DO

             END DO

           END DO

         END DO

!$OMP END DO
     
   
END SUBROUTINE  BIOLink_updateconc_BIO
#endif /* BIOLink_update_CONCBIO */


#if defined key_MARS
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
   !&E       !  2022-02  (G. Koenig) Commenting
   !&E
   !&E--------------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

   USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
   USE parameters,   ONLY : liminm1,limaxp2,ljminm1,ljmaxp2

     !====================================================================
     ! Execution of the function
     !====================================================================

      !******************* Calling of a function in hybrid *****************!
      !*************************** OMP/MPI mode ****************************!


OMPMPI barrier
OMPMPI master
     CALL_MPI ex_i_rsh(-1,2,nv_state*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,cv_wat(1:nv_state,:,liminm1:limaxp2,ljminm1:ljmaxp2))
     CALL_MPI ex_j_rsh(-1,2,nv_state*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,cv_wat(1:nv_state,:,liminm1:limaxp2,ljminm1:ljmaxp2))
OMPMPI end master
OMPMPI barrier

OMPMPI flush(cv_wat)


  END SUBROUTINE BIOLink_exchgMPI_cvwat
#endif /* key_MARS */
   
#endif /* BIOLink */

END MODULE coupleur_BIOLink

