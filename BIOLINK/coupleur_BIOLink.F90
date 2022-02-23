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
   !&E     subroutine BIOLink_update           ! evaluation of sinks and sources 
   !&E                                          terms  - routine called by step
   !&E                
   !&E     subroutine BIOLink_water_column     ! evaluation of total water height and 
   !&E                                          vertical meshes thickness - 
   !&E                                          called by BIOLink_update
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
                                                                           
    PUBLIC BIOLink_update         ! Update of BIOLink and biology model variables
                                  ! - called by main

#ifdef key_MARS
    PUBLIC BIOLink_exchgMPI_cvwat ! I do not know
#endif

#if defined key_nosubstmodule
    PUBLIC BIOLink_substance      ! Declaration of the characteristics of the tracers
                                  ! If the module substance is not declared
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

  SUBROUTINE BIOLink_water_column(ifirst,ilast,jfirst,jlast)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_water_column ***
  !&E
  !&E ** Purpose : Computation of total water thickness and thickness of 
  !&E              each vertical layer
  !&E
  !&E ** Description : We read the depth and water elevation from the 
  !&E                  hydro model and use it to compute the total water
  !&E                  thickness. Then we use the difference of height
  !&E                  of each cell to compute the thickness of each cell
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : limitation of the subgrid given by MPI ifirst,
  !&E                     ilast, jfirst and jlast
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : Creation by B. Thouvenin ( date unknown)
  !&E              Commenting by G. Koenig ( february 2022)
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

   INTEGER, INTENT(IN)      :: ifirst,ilast,jfirst,jlast
 
     !====================================================================
     ! Local declarations of variables
     !====================================================================

   INTEGER                  :: i,j,k,kmaxmod

    
     !====================================================================
     ! Execution of the function
     !====================================================================

      !******************* Water height computation ************************!

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)

   DO j=jfirst,jlast ! The loop is on the entire meridional direction

#if defined key_MARS

      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1) ! I do not know

        IF(j.GE.jb(i)+1 .AND. j .LE. jh(i)-1) THEN
#else

      DO i=ifirst,ilast ! And on the entire zonal direction

#endif /* key_MARS */

#if ! defined MUSTANG

#  if defined key_MARS

          TOTAL_WATER_HEIGHT(i,j)=BATHY_H0(i,j)+                 & ! The to
                                  WATER_ELEVATION(i,j) ! tal water thickness 
                                                       ! is computed from 
                                                       ! the depth
          
#  else

           TOTAL_WATER_HEIGHT(i,j)=z_w(i,j,N)+h(i,j) ! The total water 
                                                     ! thickness is computed
                                                     ! from the top of the 
                                                     ! last vertical layer


#  endif /* key_MARS */

#endif /* MUSTANG */

#if defined key_MARS

            IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN

              kmaxmod=1

            ELSE
 
              kmaxmod=NB_LAYER_WAT

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

          DO k=1,NB_LAYER_WAT-1          

                THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1) ! The thickness
                THICKLAYERWW(k,i,j)=z_r(i,j,k+1)-z_r(i,j,k) ! of each cell
                                                            ! is computed by
                                                            ! the height 
                                                            ! difference of 
                                                            ! the cells,
                                                            ! both from top
                                                            ! and center

          ENDDO

          k=NB_LAYER_WAT ! Last layer
          THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1)
          THICKLAYERWW(k,i,j)=0._rsh

#endif /* key_MARS */

      ENDDO
   ENDDO
!$OMP END DO

 
  END SUBROUTINE  BIOLink_water_column

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
 
  SUBROUTINE BIOLink_read_vardiag

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_read_vardiag  ***
   !&E
   !&E ** Purpose : Reading of the file describing the diagnostic variables
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : bloom_init_id,bloom_create_var_diagtracer 
   !&E                     from bloom_initdefine
   !&E                     nb_var_tracerN,nb_var_tracerP from parameters
   !&E                    
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-08 (B. Thouvenin) issued from sub_read_vardiag 
   !&E       !  2022-02 (G. Koenig) commenting 
   !&E
   !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined BLOOM

   USE bloom_initdefine, ONLY : bloom_init_id

#  if defined key_N_tracer

#    if defined key_MARS

   USE parameters, ONLY : nb_var_tracerN 

#    endif /* key_MARS */

   USE bloom_initdefine, ONLY : bloom_create_vardiagtracer

#  endif /* key_N_tracer */

#  if defined key_P_tracer

#    if defined key_MARS 

   USE parameters, ONLY : nb_var_tracerP 

#    endif /* key_MARS */

   USE bloom_initdefine, ONLY : bloom_create_vardiagtracer

#  endif /* key_P_tracer */

#endif /* BLOOM */

     !====================================================================
     ! Local declarations of variables
     !====================================================================

   LOGICAL               :: ex,l_diag_wat,l_diag_sed
   INTEGER               :: eof,isubs,isubs_r,dimvar,it,ind_white,IERR_MPI
   CHARACTER(LEN=lchain) :: namvar_r,long_name_var_r,standard_name_var_r,unitvar_r
   CHARACTER(LEN=5)      :: comment

#if defined key_N_tracer || defined key_P_tracer

   INTEGER               :: is,id,ivtra

#endif /* key_N_tracer || key_P_tracer */

     !====================================================================
     ! Execution of the function
     !====================================================================

      !********************** Saving into simu.log *************************!

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
      l_out_subs_diag=.false.  ! no saving of diagnoses variable by default

#if defined MUSTANG

      l_out_subs_diag_sed=.false.  ! no saving of diagnoses variable by default

#endif /* MUSTANG */

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

#if defined BLOOM

#  if defined key_N_tracer

       ndiag_tracerN=0

       DO is=1,nb_source_tracerN

            ndiag_tot=ndiag_tot+nb_var_tracerN+1
            ndiag_tracerN=ndiag_tracerN+nb_var_tracerN+1

#    if defined key_age_tracer

            ndiag_tot=ndiag_tot+nb_var_tracerN+1
            ndiag_tracerN=ndiag_tracerN+nb_var_tracerN+1

#    endif /* key_age_tracer */

       END DO ! is

#  endif /* key_N_tracer */

#  if defined key_P_tracer

       ndiag_tracerP=0

       DO is=1,nb_source_tracerP

            ndiag_tot=ndiag_tot+nb_var_tracerP+1
            ndiag_tracerP=ndiag_tracerP+nb_var_tracerP+1

#    if defined key_age_tracer

            ndiag_tot=ndiag_tot+nb_var_tracerP+1
            ndiag_tracerP=ndiag_tracerP+nb_var_tracerP+1

#    endif /* key_age_tracer */

       END DO ! is

#  endif / key_P_tracer */

#endif /* BLOOM */

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

#endif /* MUSTANG */

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

   END IF

#if defined BLOOM && (defined key_N_tracer || defined key_P_tracer)

       call bloom_create_vardiagtracer

#endif /* BLOOM && ( key_N_tracer || key_P_tracer ) */

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
#endif /* MUSTANG */

   ! allocate diagnostic variables
   ALLOCATE( diag_1d(1:ndiag_1d) )
   diag_1d(:)=0.0_rsh

   ALLOCATE( diag_2d(ndiag_1d+1:ndiag_1d+ndiag_2d,PROC_IN_ARRAY) )
   diag_2d(:,:,:)=0.0_rsh

   ALLOCATE( diag_3d_wat(ndiag_2d+1:ndiag_2d+ndiag_3d_wat,NB_LAYER_WAT,PROC_IN_ARRAY) )
   diag_3d_wat(:,:,:,:)=0.0_rsh

#if defined MUSTANG && defined key_BLOOM_insed

   ALLOCATE( diag_3D_sed(ndiag_tot-ndiag_3d_sed+1:ndiag_tot,ksdmin:ksdmax,PROC_IN_ARRAY) ) 
   diag_3d_sed(:,:,:,:)=0.0_rsh

   ALLOCATE( diag_2D_sed(ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d,PROC_IN_ARRAY) )
   diag_2D_sed(:,:,:)=0.0_rsh

#endif /* MUSTANG && key_BLOOM_insed */

   ! Storage of diagnostic variables within reading order
#if defined BLOOM

   it = 0
   ! dimvar=1 : diag1D, =2 diag2D
   ! dimvar=3 : diag3D in wat only
   ! dimvar=4 : diag 3D in wat and in sed
   ! dimvar=5 : diag3D in sed only
   ! dimvar=6 : diag2D in sed only

#  if defined key_BLOOM_insed
   ! dimvar=1 : diag1D, =2 diag2D

   DO dimvar = 1,2

     DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
!         irk_diag(it)=isubs
         irk_diag(isubs)=it

#    if defined key_N_tracer

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else

         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif /* key_N_tracer */

       END IF

     END DO

   END DO

   ! dimvar=6 : diag2D in sed only
   dimvar=6

   DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
         irk_diag(isubs)=it

#    if defined key_N_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else
         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif /* key_N_P_tracer && BLOOM */

       END IF

   END DO

   ! dimvar=3 : diag3D in wat only
   ! dimvar=4 : diag 3D in wat and in sed
   ! dimvar=5 : diag3D in sed only

   DO dimvar = 3,5

     DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1

         irk_diag(isubs)=it

#    if defined key_N_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else

         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif

       END IF

     END DO

   END DO

#  else

   DO dimvar = 1,6

     DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
         irk_diag(isubs)=it

#    if defined key_N_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else

         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif /* key_N_P_tracer && BLOOM */

       END IF

     END DO

   END DO
#  endif /* BLOOM */

#endif /* key_BLOOM_insed */

#if defined PEPTIC

   it = 0

   DO dimvar = 1,4

    DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
         irk_diag(isubs)=it

       END IF

    END DO

   END DO

#endif /* PEPTIC */

   IF_MPI (MASTER) THEN

     DO isubs = 1,ndiag_tot

        MPI_master_only WRITE(iscreenlog,*) isubs,TRIM(name_vardiag(isubs)),idimv_r(isubs),irk_diag(isubs)

     END DO

   ENDIF_MPI

  END SUBROUTINE BIOLink_read_vardiag



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

