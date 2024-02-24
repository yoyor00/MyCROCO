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

#include "coupleur_define_BIOLink.h" 
                                     ! Equivalence of names between the hydro and the different
                                     ! Tracer and biological model. Also contains a function
                                     ! For the allocation of the main tables

   USE module_BIOLink ! script that groups all the files of BIOLink together and allows the 
                      ! access to their functions/subroutines/variables
   USE comsubstance   ! Access to the functions/variables of SUBSTANCE ( from the MUSTANG
                      ! sediment model)

   USE comBIOLink     ! Common variables of the BIOLink coupler

   USE comBIOLink_helping ! Variables for the helping functions of BIOLink

   USE comBIOLink_physics ! Variables for the physics of BIOLink

   USE coupleur_BIOLink_helping ! For the helping functions of BIOLink

   USE coupleur_BIOLink_physics ! For the physics related functions of 
                                ! BIOLink
#if defined MUSTANG

   USE comMUSTANG ,  ONLY : htot ! Height of the water column
   USE comMUSTANG ,  ONLY : var3D_diagsed, var2D_diagsed

#endif /*MUSTANG*/

#if defined ECO3M   

   USE mod_eco3m_irrad, ONLY : irrad
   USE mod_eco3m, ONLY : nx_min, nx_max, ny_min, ny_max, nz_max, VAR
   USE mod_eco3m, ONLY : ECO_TIME_STEP, WATCONCPOS_tab, TEMP_BIOLink,&
                         SAL_BIOLink,BIO_SKSC_ADV, THICKLAYERWC 
   USE mod_eco3m, ONLY : TBIO 

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
 
    PUBLIC BIOLink_init_main      ! Initialization of some of BIOLink tables
                                  ! - called by main
                                                                           
    PUBLIC BIOLink_alloc          ! Allocation of different tables
    PUBLIC BIOLink_update         ! Update of BIOLink and biology model variables
                                  ! - called by main


!*****************************************************************************************!
!                                                                                         !
!                         Routines and functions                                          !
!                                                                                         !
!*****************************************************************************************!

 CONTAINS

  !!======================================================================

  SUBROUTINE BIOLink_initialization(icall, tile)
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
  !&E       !  2022-06  (M. Baklouti et C. Mazoyer) : Eco3M
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

#endif /* PEPTIC/BLOOM/METeOR */


     !====================================================================
     ! External arguments
     !====================================================================

  INTEGER, INTENT(IN) :: icall ! To differentiate the two calls to the function
                               ! The first call reads the parameter files
                               ! and the second initialize the tracer/biological
                               ! model 
  INTEGER, INTENT(IN) :: tile



     !====================================================================
     ! Local declarations of variables
     !====================================================================
  
  INTEGER               :: isubs
#if defined ECO3M
  INTEGER               :: i,j,k, kmaxmod, iv
  INTEGER               :: ifirst, ilast, jfirst, jlast
#endif
     !====================================================================
     ! Execution of the function
     !====================================================================
# include "compute_tile_bounds.h"

   IF (icall==0) THEN ! First call of the function, we read the parameter
                      ! files

#if defined PEPTIC

    CALL peptic_param('r')

#elif defined BLOOM

    CALL bloom_param('r')

#elif defined METeOR

    CALL meteor_param('r')

#endif /* PEPTIC/BLOOM/METeOR */

#if defined ECO3M
   ELSE IF (icall==1) THEN 
#else
   ELSE ! Second call of the function, now we initialize the variable of the 
        ! Biology/Tracer models
#endif

     TIME_BEGIN=CURRENT_TIME ! Time begin is the beginning time for BIOLink and 
                             ! biological/tracer models, while current time
                             ! is the time from the hydro model

     CALL BIOLink_alloc() ! Allocation of the variables used by BLOOM :
                          ! Source and sink terms,
                          ! Height of the water column
                          ! and light attenuations related variables
       

! Initialization of the tracer/fixed variables in biological/tracer models
! with their names and characteristics

#if defined BLOOM 

     DO isubs=1,nv_adv ! Counter for tracer variables

       CALL bloom_init_iv(isubs,standard_name_var(isubs),1)

     END DO ! isubs

     DO isubs=nv_adv+1,nv_adv+nv_fix !Counter for fixed variables

       CALL bloom_init_iv(isubs,standard_name_var_fix(isubs-nv_adv),2)

     END DO ! isubs

#  ifdef key_benthic

     DO isubs=nv_adv+nv_fix+nv_bent ! Counter for benthic variables

       CALL bloom_init_iv(isubs,standard_name_var_bent(isubs-nv_adv-nv_fix),3)

     END DO ! isubs

#  endif /* key_benthic */

#elif defined PEPTIC

     CALL peptic_alloc_var ! For PEPTIC all the variables 
                           ! are declared in the same function

#elif defined ECO3M 
     nx_min = Istr 
     nx_max = Iend
     ny_min = Jstr
     ny_max = Jend
     nz_max = NB_LAYER_WAT
     CALL eco3m_init_config 

#endif /* BLOOM/PEPTIC/ECO3M */
  
  ! Writing of the biological parameters used in an external file

#if defined PEPTIC

     CALL peptic_param('w')

#elif defined BLOOM

     CALL bloom_param('w')

#elif defined METeOR

     CALL meteor_param('w')

#endif /* PEPTIC/BLOOM/METeOR */

  ! Allocation of diagnostic variables for BLOOM

#if defined BLOOM && defined DIAGNOSTICS_BIO

     CALL BIOLink_read_vardiag

#endif /* BLOOM && DIAGNOSTICS_BIO*/

#if defined ECO3M 
   ELSE IF (icall==2) THEN 
        ! Third call of the function, now we initialize the variable of the 
        ! Biology/Tracer models with 1DV text files
        ! only for DYFAMED testcase

#if defined DYFAMED
     ! .dat containing values of DYFAMED is read through eco3m_init_config 
     ifirst=Istr
     ilast=Iend
     jfirst=Jstr
     jlast=Jend
     kmaxmod=NB_LAYER_WAT 
     DO iv=1,nv_adv
       DO k=1,kmaxmod
         DO j=jfirst,jlast
           DO i=ifirst,ilast
              ! variables d etat transportees
                ! duplication of water column (i=1, j=1) for all i,j
                WAT_CONCADV_ivkij_tini=WATCONCPOS_ivi1j1k
           END DO
         END DO
       END DO    
     END DO
#endif /* DYFAMED */
#endif /* ECO3M */


   ENDIF ! /* icall parameter */

   MPI_master_only PRINT*, 'END of BIOLink_initialization'


  END SUBROUTINE BIOLink_initialization

!!======================================================================

  SUBROUTINE BIOLink_init_main (tile)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_init_main  ***
  !&E
  !&E ** Purpose : passage de tile à Istr,Iend,Jstr,Jend
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
  !&E       !  2022-09 Solène
  !&E
  !&E---------------------------------------------------------------------

    integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"
    CALL BIOLink_init(Istr,Iend,Jstr,Jend)
  END SUBROUTINE BIOLink_init_main

!!========================================================================

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
#endif /* BLOOM */

#if defined METeOR

   USE meteor_initdefine, ONLY : meteor_read_react      ! I do not know

#endif /* METeOR */

  IMPLICIT NONE

     !====================================================================
     ! External arguments
     !====================================================================


!  INTEGER, INTENT(IN)  :: tile
     !====================================================================
     ! Local declarations of variables
     !====================================================================
                                                              
  INTEGER               :: i,j,k,iv ! Counter variables
  INTEGER               :: ifirst, ilast, jfirst, jlast
  CHARACTER(LEN=lchain) :: logfilename ! Name of the file to write the log

     !====================================================================
     ! Execution of the function
     !====================================================================
!# include "compute_tile_bounds.h"

     !************** Determination of the time steps *********************!
     
  t_bio=CURRENT_TIME+TRANSPORT_TIME_STEP ! The biology time steps is updated 
                                          ! From the hydrodynamical model 
                                          ! With a transport time step

# ifdef ECO3M
     ! transfer of current time to the biologial code
     TBIO=t_bio
# endif


     !********************* Sinking velocity *****************************!

!     ifirst=Istr
!     ilast=Iend
!     jfirst=Jstr
!     jlast=Jend
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

     !*********** Satellital measurement of suspended matter ****************!

#if defined key_messat
   ! initialization of the suspended matter measured by satellite
   messat(:,:,:)=0.0_rsh

   DO j=jfirst,jlast

     DO i=ifirst,ilast

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

#if ! defined ECO3M  

      ALLOCATE(WATCONCPOS(nv_adv,NB_LAYER_WAT,PROC_IN_ARRAY))
      WATCONCPOS(:,:,:,:)=0.0_rsh

#endif /* ECO3M */

     !===================================================================
     ! Table of positive concentrations for fixed variables
     !===================================================================

      ALLOCATE(FIXCONCPOS(nv_fix,NB_LAYER_WAT,PROC_IN_ARRAY))
      FIXCONCPOS(:,:,:,:)=0.0_rsh

     !===================================================================
     ! Table of positive concentrations for benthic variables
     !===================================================================

#ifdef key_benthic
      ALLOCATE(BENTCONCPOS(nv_bent,PROC_IN_ARRAY))
      BENTCONCPOS(:,:,:,:)=0.0_rsh

#endif /* key_benthic */ 

     !===================================================================
     ! Table of source and sink terms for tracer variables
     !===================================================================

#if ! defined ECO3M  
      ALLOCATE( BIO_SINKSOURCES(ARRAY_SINKSOURCES))
      BIO_SINKSOURCES(:,:,:,:)=0.0_rsh
#endif /* ECO3M */

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

#if ! defined MUSTANG

        ALLOCATE( TOTAL_WATER_HEIGHT(PROC_IN_ARRAY_m2p2))

#endif /* MUSTANG */
        
#if ! defined ECO3M  
        ALLOCATE( THICKLAYERWC(NB_LAYER_WAT,PROC_IN_ARRAY))
        THICKLAYERWC(:,:,:)=0.0_rsh
#endif        
        ALLOCATE( THICKLAYERWW(NB_LAYER_WAT,PROC_IN_ARRAY))
        THICKLAYERWW(:,:,:)=0.0_rsh

     !=====================================================================
     !  Temperature and salinity variables
     !=====================================================================
#  if !defined ECO3M
      ALLOCATE( TEMP_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
      TEMP_BIOLink(:,:,:)=0.0_rsh

      ALLOCATE( SAL_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
      SAL_BIOLink(:,:,:)=0.0_rsh
#  endif /* ECO3M */
     !=====================================================================
     !  Bottom current variables
     !=====================================================================

#if defined BLOOM && defined key_benthos

      ALLOCATE(BOTTCURRENTBIO(PROC_IN_ARRAY)) 

#endif /* BLOOM && key_benthos */


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

#if defined  BIOLink_PAR_eval

      ALLOCATE( PAR_top_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )
      PAR_top_layer(:,:,:)=0.0_rsh
      ALLOCATE( PAR(PROC_IN_ARRAY,NB_LAYER_WAT) ) ! It is smaller than 
      PAR(:,:,:)=0.0_rsh                            ! PAR_top_layer
                                                    ! because there seems
                                                    ! to be one MUSTANG
                                                    ! cell in PAR_top_layer

#  if defined PEPTIC

      ALLOCATE( PAR_avg_layer_phyto(1,NB_LAYER_WAT,PROC_IN_ARRAY) )
      PAR_top_layer_day(:,:,:)=0.0_rsh
      ALLOCATE( PAR_top_layer_day(NB_LAYER_WAT,PROC_IN_ARRAY))
      PAR_avg_layer_phyto(1,:,:,:)=0.0_rsh

#  elif defined METeOR

      ALLOCATE( Flimrad_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )

#  endif /* PEPTIC/METEOR */

#endif /* BIOLink_PAR_eval */

     ! Extinction and attenuation

#if defined  BIOLink_PAR_eval

      ALLOCATE( EXTINCTION_RAD(NB_LAYER_WAT,PROC_IN_ARRAY) )
      EXTINCTION_RAD(:,:,:)=0.0_rsh

#endif /* BIOLink_PAR_eval */

#if defined BLOOM

      ALLOCATE(extinction_ave4d(NB_LAYER_WAT,PROC_IN_ARRAY))
      extinction_ave4d(:,:,:)=0.0_rsh
      ALLOCATE(extinction_aveh(NB_LAYER_WAT,PROC_IN_ARRAY))
      extinction_aveh(:,:,:)=0.0_rsh
      ALLOCATE(extinction_tab(4,24,NB_LAYER_WAT,PROC_IN_ARRAY))
      extinction_tab(:,:,:,:,:)=0.0_rsh
      t_cum_extinctionh=0.0_rsh
      ihour_previous=0

#endif /* BLOOM */

     ! Variables that affect the extinction
 
      ALLOCATE( SPMTOT_MGL(NB_LAYER_WAT,PROC_IN_ARRAY))
      SPMTOT_MGL(:,:,:)=0.0_rsh

#if defined  BIOLink_PAR_eval

      ALLOCATE( BIOLink_chloro(NB_LAYER_WAT,PROC_IN_ARRAY) )  
      BIOLink_chloro(:,:,:)=0.0_rsh

#endif /* BIOLink_PAR_eval */

#if defined BLOOM

      ALLOCATE( effetlumiere_day_diat(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_diat(:,:,:)=0.0_rsh
      ALLOCATE( effetlumiere_day_dino(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_dino(:,:,:)=0.0_rsh
      ALLOCATE( effetlumiere_day_nano(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_nano(:,:,:)=0.0_rsh

#  if defined key_karenia

      ALLOCATE( effetlumiere_day_karenia(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_karenia(:,:,:)=0.0_rsh

#  endif /* key_karenia */

#  if defined key_psnz

      ALLOCATE( effetlumiere_day_psnz(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_psnz(:,:,:)=0.0_rsh

#  endif /* key_psnz */

#  if defined key_phaeocystis

      ALLOCATE( effetlumiere_day_phaeocystis(NB_LAYER_WAT,PROC_IN_ARRAY) )
      effetlumiere_day_phaeocystis(:,:,:)=0.0_rsh

#  endif /* key_phaeocystis */

#endif /* BLOOM */

     !=====================================================================
     !  Wind variables
     !=====================================================================

#if ! defined BULK_FLUX

      ALLOCATE( WIND_SPEED(PROC_IN_ARRAY))
      WIND_SPEED(:,:)=0.0_rsh

#endif /* BULK_FLUX */

     !=====================================================================
     !  Estimation of the suspended matter via satellite
     !=====================================================================

#ifdef key_messat

      ALLOCATE(messat(PROC_IN_ARRAY,51))
      messat(:,:,:)=0.0_rsh

#endif /* key_messat */

     !=====================================================================
     !  Oyster related variables
     !=====================================================================

#if defined BLOOM

#  if defined key_oyster_SFG

      ALLOCATE(tpostpontecoq(PROC_IN_ARRAY))
      ALLOCATE(tpostpontegam(PROC_IN_ARRAY))

#  endif /* key_oyster_SFG */

#  if defined key_oyster_benthos || defined key_oyster_DEB || defined key_oyster_SFG

      ALLOCATE(nbhuitre(PROC_IN_ARRAY))
      hautable(:,:)=0.0_rsh
      ALLOCATE(hautable(PROC_IN_ARRAY))
      nbhuitre(:,:)=0.0_rsh

#  endif /* key_oyster_benthos */

#endif /* BLOOM */


END SUBROUTINE  BIOLink_alloc

  SUBROUTINE BIOLink_update(ifirst,ilast,jfirst,jlast   &
#if defined key_oyster_SFG || defined key_oyster_DEB
          , CELL_SURF                                                &
#endif
         )
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_update  ***
  !&E
  !&E ** Purpose : Update of the BIOLink concentration, sources and sinks
  !&E              and helping variables (settling velocities, number of 
  !&E             oyster, etc...) at each time step.
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



#elif defined METeOR
  USE meteor,          ONLY : meteor_sksc_wat,meteor_reac_equi
                              ! Sources and sink terms and
                              ! Equilibrium reactions 
#endif /* PEPTIC/BLOOM/METeOR */

     !====================================================================
     ! External arguments
     !====================================================================

   INTEGER, INTENT(IN) :: ifirst,ilast,jfirst,jlast
   
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

#if defined PHYBIO_2ways
   integer                :: irgb
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

# ifdef ECO3M
        ! transfer of current time to the biologial code
        TBIO=t_bio+BIO_TIME_STEP
# endif


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

#  endif /* key_messat */

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

#     if defined key_messat

                              ,forcSPM                          &

#     endif /* key_messat */

                         )

#  endif /* BLOOM */  

      CALL BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate) ! Computation 
                                                             ! of solar ra
                                                             ! -diation and
                                                             ! extinction


      PAR = BIOLink2hydro_3D(ifirst,ilast,jfirst,jlast,1,NB_LAYER_WAT,     &
                                PAR_top_layer,0,NB_LAYER_WAT) ! We do not take the first
                                               ! layer because it is appa
                                               ! rently related to the 
                                               ! interface with MUSTANG
      
!        MPI_master_only PRINT*,"ifirst,ilast,jfirst,jlast",ifirst,ilast,jfirst,jlast
!        MPI_master_only PRINT*,"shape PAR_top_layer",shape(PAR_top_layer)
!        MPI_master_only PRINT*,"Shape PAR",shape(PAR)

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


      !!=================================================================================
      ! to compute the fraction of the solar shortwave flux "swdk" 
      ! penetrating to grid level depth (at vertical w-points)
      ! for the physical model (2 ways coupling)
      !!=================================================================================
#if defined PHYBIO_2ways 
      ALLOCATE(swdk_bio(ifirst:ilast,jfirst:jlast,NB_LAYER_WAT) )
      swdk_bio=0.0
      do k=1,NB_LAYER_WAT
          do j=jfirst, jlast
              do i=ifirst, ilast
                  do irgb=1,3
                      swdk_bio(i,j,k)=swdk_bio(i,j,k)+ &
    &                          (1.d0/3.d0)*E_PARZ_RGB(irgb,i,j,k)/E_PAR(i,j)
                  enddo
              enddo
          enddo
      enddo
#endif
      !!=================================================================================

#if defined ECO3M
      irrad(ifirst:ilast,jfirst:jlast)=SOLAR_RAD(ifirst:ilast,jfirst:jlast)/(RAD_SRFSCALE) ! net short wave solar flux (W/m2)
#endif



      !********************* Diagnostic variables **************************!

#if defined BLOOM

      CALL  bloom_eval_diag2d(ifirst,ilast,jfirst,jlast) ! Bloom function
                                                         ! for diagnostic
                                                         ! evaluations
      ! Here I convert the shape of the 3D diagnostics so that CROCO can use them
!      do i=1,ndiag_3d

!        diag_3d_CROCO(i,:,:,:) = BIOLink2hydro_3D(ifirst,ilast,jfirst,jlast,1,NB_LAYER_WAT, &
!                                diag_3d_wat(i,:,:,:),1,NB_LAYER_WAT)
  
!      end do



#endif /* BLOOM */

      !**************** Wind speed for some BLOOM keys *********************!

#if defined key_oxygen || defined key_zostera || defined METeOR

#  if ! defined BULK_FLUX

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)
     DO j=jfirst,jlast
     DO i=ifirst,ilast

        ! wind speed (m.s-1) evaluated from the surface stress values
        ! sustr and svstr  (m2.s-2)
        ! with rho_air=1.3  ; CD=0.0014  
        WIND_SPEED(i,j) = sqrt(sqrt( (0.5*(sustr(i,j)+sustr(i+1,j)))**2  &
                            +(0.5*(svstr(i,j)+svstr(i,j+1)))**2)        &
                            *RHOREF/(1.3*0.0014))
!     &                       *rho0/(rho_air*CD))
      ENDDO
      ENDDO
!$OMP END DO

#  endif /* BULK_FLUX */

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

                                 ) 
#elif defined METeOR

    IF ( l_treat_re_dc) THEN
        CALL meteor_sksc_wat(ifirst,ilast,jfirst,jlast,WIND_SPEED)
    ENDIF                      
#elif defined ECO3M
     CALL eco3m_main

#endif /* BLOOM/METeOR/ECO3M */
   
      !**************** Conversion of the order of arrays  *****************!

      CALL BIOLink2hydro(ifirst,ilast,jfirst,jlast)     

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

#  endif /* METeOR && ! BLOOM && ! PEPTIC && ! ECO3M */   

      CALL BIOLink_updateconc_BIO(ifirst,ilast,jfirst,jlast)

#  if defined METeOR && ! defined BLOOM && ! defined PEPTIC && ! defined ECO3M

    ENDIF 

#  endif /* METeOR && ! BLOOM && ! PEPTIC && ! ECO3M */

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

         ! I convert the array to send it in CROCO
         PAR = BIOLink2hydro_3D(ifirst,ilast,jfirst,jlast,1,NB_LAYER_WAT, &
                                PAR_top_layer,0,NB_LAYER_WAT) ! We do not take the first
                                               ! layer because it is appa
                                               ! rently related to the 
                                               ! interface with MUSTANG
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

            END DO

          END DO

        END DO

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modules BIO ECO3M : index ordre : (iv)%i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ECO3M 
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,iv)
! no fixed variables ocean model (i,j,k) to bio (i,j,k)
     kmaxmod=NB_LAYER_WAT
     DO k=1,kmaxmod
       DO j=jfirst,jlast
         DO i=ifirst,ilast
            TEMP_BIOLink(i,j,k)=TEMPHYDRO_ijk
            SAL_BIOLink(i,j,k)=SALHYDRO_ijk
         END DO
       END DO    
     END DO
     DO iv=1,nv_adv
       DO k=1,kmaxmod
         DO j=jfirst,jlast
           DO i=ifirst,ilast
              ! variables d etat transportees
                WATCONCPOS_ivijk=WAT_CONCADV_ivkij
           END DO
         END DO
       END DO    
     END DO
!$OMP END DO
#endif

#if defined PEPTIC || defined BLOOM 
! PEPTIC and BLOOM are in the order : iv (tracer), k (vertical),
! i (zonal) and j (meridional)
! but CROCO is in the order : iv (tracer), i (zonal), j (meridional)
! and k (vertical). 

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,iv)

     DO j=jfirst,jlast
     
       DO i=ifirst,ilast

         kmaxmod=NB_LAYER_WAT

         WATCONCPOS(:,:,i,j)=0.0_rsh
         FIXCONCPOS(:,:,i,j)=0.0_rsh

#  if defined key_benthic

         BENTCONCPOS(:,i,j)=0.0_rsh

#  endif /* key_benthic */
 
         DO k=1,kmaxmod

           TEMP_BIOLink(k,i,j)=TEMPHYDRO_ijk ! Inverting of the place of
           SAL_BIOLink(k,i,j)=SALHYDRO_ijk   ! the vertical index

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


#endif /* BLOOM/PEPTIC/ECO3M */

   END SUBROUTINE  BIOLink_convarray



    !!============================================================================== 
 
  FUNCTION BIOLink2hydro_3D(ifirst,ilast,jfirst,jlast,kfirst,klast,VAR,zstart_var,zend_var)     

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink2hydro ***
  !&E
  !&E ** Purpose : Conversion of a 3 variables from the BIOLink order (z,x,y)
  !&E              to the CROCO order (x,y,z)
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
  !&E ** History : ! Created by G. Koenig (march 2022)
  !&E
  !&E---------------------------------------------------------------------


     !====================================================================
     ! External arguments
     !====================================================================
 
   INTEGER, INTENT(IN)    :: ifirst,ilast,jfirst,jlast ! Limits of the MPI
                                                       ! subdomain
   INTEGER, INTENT(IN)    :: kfirst,klast ! Limits of the vertical subdomain

   INTEGER, INTENT(IN)    :: zstart_var,zend_var ! ( Vertical index for the 
                                                 ! z dimension )
   REAL(KIND=rsh),DIMENSION(zstart_var:zend_var,PROC_IN_ARRAY), INTENT(IN)   :: VAR
                                             ! Input array to be modified

     !====================================================================
     ! Local declarations of variables
     !====================================================================

    INTEGER                  :: i,j,k ! Spatial and tracer counters
    INTEGER                  :: i1,i2,i3 ! Internal BIOLink counters
    REAL(KIND=rsh), DIMENSION(PROC_IN_ARRAY,kfirst:klast) :: BIOLink2hydro_3D ! Array for
                                                           ! returning the 
                                                           ! transformed
                                                           ! variable
    
     !====================================================================
     ! Execution of the function
     !====================================================================

 !******************Allocation of the array for storing the variables*******!

     BIOLink2hydro_3D(:,:,:)=0.0_rsh


 !****************** Conversion for the index order of *********************!
 !***************** BIOLink to the one of the hydro model ******************!


!$OMP DO SCHEDULE(RUNTIME)
     DO i = ifirst,ilast

       DO j = jfirst,jlast

          DO k=kfirst,klast
              
              BIOLink2hydro_3D(i,j,k) = VAR(k,i,j)

          END DO

        END DO

      END DO


!$OMP END DO

END FUNCTION  BIOLink2hydro_3D



  SUBROUTINE BIOLink2hydro(ifirst,ilast,jfirst,jlast  &

#if defined key_nosubstmodule

                               ,WAT_SETTL  &

#endif
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

    INTEGER                  :: i,j,ks,iv ! Spatial and tracer counters
    INTEGER                  :: i1,i2,i3,i4 ! Internal BIOLink counters
    INTEGER                  :: isubs,ind_diag2d,ind_diag3d_sed,ind_diag3d_wat
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
# include "diagnostics.h"

      !***************** We test if the tracer concentration ***************!
      !************************* becomes negative **************************!
#if ! defined ECO3M
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
#endif

!$OMP DO SCHEDULE(RUNTIME)
         DO i1=IRANGE1

           DO i2=IRANGE2

             DO i3=IRANGE3

               DO i4=IRANGE4  
#if !defined ECO3M
                 WATER_CONCENTRATION(WATCONC_INDEX_EQ) = WATER_CONCENTRATION(WATCONC_INDEX_EQ)    &
                                                         + xnegtr(i4,i3,i2) &
                                                         * BIO_SKSC_ADV(BIOSKSC_INDEX_EQ) &
                                                         * TRANSPORT_TIME_STEP    
#else
                 WATER_CONCENTRATION(WATCONC_INDEX_EQ) = WATER_CONCENTRATION(WATCONC_INDEX_EQ)    &
                                                         + BIO_SKSC_ADV(BIOSKSC_INDEX_EQ) 
#endif
                 ! Computation of the flux from the advection terms 
                 ! with the transport time step
               END DO

             END DO

           END DO

         END DO

!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Mise à jourdes variables Diagnostiques Bio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (mplus jan 2024)
#ifdef DIAGNOSTICS_BIO

!$OMP DO SCHEDULE(RUNTIME)
        ind_diag2d = 0
        ind_diag3d_wat = 0
        ind_diag3d_sed = 0
        DO isubs=1,ndiag_tot
          IF(idimv_r(isubs) == 2 .OR. idimv_r(isubs) == 6) THEN
            ind_diag2d = ind_diag2d + 1 
            DO i3=IRANGE3
              DO i4=IRANGE4
                biodiag2d(i4,i3,ind_diag2d) = diag_2d(irk_diag(isubs),i4,i3)
              END DO
            END DO
          ELSEIF(idimv_r(isubs) == 3) THEN
            ind_diag3d_wat = ind_diag3d_wat + 1
            DO i2=IRANGE2
              DO i3=IRANGE3
                DO i4=IRANGE4
                  biodiag3d_wat(i4,i3,i2,ind_diag3d_wat) = diag_3d_wat(irk_diag(isubs),i2,i4,i3)
                END DO
              END DO
            END DO
#ifdef key_BLOOM_insed
          ELSEIF(idimv_r(isubs) == 5) THEN
            ind_diag3d_sed = ind_diag3d_sed + 1
            DO ks=ksdmin,ksdmax
              DO i3=IRANGE3
                DO i4=IRANGE4
                  biodiag3d_sed(i4,i3,ks,ind_diag3d_sed) = diag_3d_sed(irk_diag(isubs),ks,i4,i3)
                END DO
              END DO
            END DO
#endif
          ENDIF
        END DO

!$OMP END DO
#endif
     
END SUBROUTINE  BIOLink_updateconc_BIO
#endif /* BIOLink_update_CONCBIO */

#endif /* BIOLink */

END MODULE coupleur_BIOLink

