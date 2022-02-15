!---------------------------------------------------------------------------
!
                     MODULE comBIOLink
!
!---------------------------------------------------------------------------
#include "cppdefs.h"

#if defined SUBSTANCE && defined BIOLink

  !!======================================================================================!!
  !!                   ***  MODULE  comBIOLink  ***                                       !!
  !! Ocean dynamics Bio :  declare and initialize all common variables                    !!
  !!                       related to BIOLink module                                      !!
  !!                       + routines related to radiation, extinction, etc               !!
  !!                       + routine related to water column                              !!
  !!                                                                                      !!
  !!   History :                                                                          !!
  !!    !  2019-08 (B. Thouvenin) issued from modules BLOOM (old ECOMARS/ key_biolo)      !!
  !!       and peptic for portability adaptation                                          !!
  !!      !            (V. Garnier, M. Sourisseau, A. Menesguen, P. Cugier, A. Arancio)   !!
  !!                                                                                      !!
  !!======================================================================================!!

#include "coupleur_define_BIOLink.h"

      !! * Modules used
      USE comsubstance

      IMPLICIT NONE

     !*************************************************************************!
     !*************************************************************************!
     !*************** BIOLink tracer and time variables ***********************!
     !*********************** and subroutines *********************************!
     !*************************************************************************!
     
     !====================================================================
     ! Public subroutines
     !====================================================================

      PUBLIC BIOLink_alloc ! Subroutines for allocating the space of variables in BIOLink

     !====================================================================
     ! Local variables of BIOLink
     !====================================================================

      REAL(kind=rsh)  ,PARAMETER                :: epsilon_BIOLink=1.e-10
      ! Variable for the comparison of float internal to BIOLink
      CHARACTER(LEN=lchain)                      :: suffix_fileres
      ! File for ?
 
     !=====================================================================
     !  Time-related variables: Timesteps and dates
     !=====================================================================

      INTEGER,PUBLIC        ::  ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK,  &
                                iheure_BIOLINK,iminu_BIOLINK,isec_BIOLINK, &
                                jjulien_BIOLINK ! Variables used in the BLOOM module that 
                                                ! Relies on julian days

      REAL(KIND=rlg),PUBLIC                     :: t_bio,BIO_TIME_STEP,DT_CONSERV_BIOLINK
      !Time and time steps variables used in the biological models 
      REAL(KIND=rlg)                             :: ECO_TIME_STEP,TIME_START_ECO 
      ! Timestep of the hydrodynamical model
      REAL(KIND=rlg)                             :: TIME_BEGIN  
      ! Date of the start of the model run
      
     !===================================================================
     ! Positive concentrations for tracer variables
     !===================================================================

#  if !defined ECO3M 
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: WATCONCPOS 
      ! Concentration in the water column, the order is the one of BIOLink for the storage 
      ! which means (species_index,vertical_direction,zonal_direction,meridional_direction)
#  endif /* ECO3M */

     !===================================================================
     ! Positive concentrations for fixed variables
     !===================================================================

      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: FIXCONCPOS  
      ! Concentration of fixed variables, the order is the one of BIOLink for the storage

     !===================================================================
     ! Positive concentrations for benthic variables
     !===================================================================

      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BENTCONCPOS 
      ! Concentration of benthic variables, the order is the one of BIOLink for the storage

     !===================================================================
     ! Sources and sink terms for tracer variables
     !===================================================================
      
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: BIO_SINKSOURCES 
      ! For the variables in the water column

     !===================================================================
     ! Sources and sink terms for fixed variables
     !===================================================================

      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: BIO_SKSC_FIX 
      ! For the fixed variables

     !=====================================================================
     !  Booleans for the activation of helping functions
     !=====================================================================   
      
      LOGICAL, PUBLIC                           :: l_bioretro_extinct
      ! For the light attenuation (PAR) routines
      LOGICAL, PUBLIC                           :: l_waterdensity_known
      ! For the computation of water density if it is not provided
      ! by the hydrodynamic model
      LOGICAL,ALLOCATABLE,DIMENSION(:),PUBLIC :: l_diagBIOLink_out
      ! Booleans for the writing in external file of diagnostic variables of BIOLink
     
     !=====================================================================
     !  Variables related to the verification of the conservation routine
     !=====================================================================     

      INTEGER, PUBLIC                           :: IVERIF_BIOLINK,JVERIF_BIOLINK
      ! Counter variables for the verification loops

#  if defined key_BIOLink_verif_conserv 
      INTEGER, PUBLIC , PARAMETER               :: iscreenlog_conserv=66
      ! Printing variable for the verification of conservativity
#  endif /* key_BIOLink_verif_conserv */
 
     
     !*************************************************************************!
     !*************************************************************************!
     !****************** Variables for the biological or **********************!
     !*********************** chemical model **********************************!
     !*************************************************************************!
     !*************************************************************************!
     
     !=====================================================================
     !  Mystery variables
     !=====================================================================  

      INTEGER,ALLOCATABLE,DIMENSION(:),PUBLIC :: idimv_r
      ! I am clueless about this variable

     !=====================================================================
     !  Internal variables of the biological models
     !=====================================================================  

#  if defined BLOOM
#include "combloom.h"

#  elif defined PEPTIC
#include "compeptic.h"

#  elif defined METeOR 
#include "commeteor.h"

#  endif /* Biological models */

     !=====================================================================
     !  General diagnostic variables 
     !=====================================================================  

      CHARACTER(LEN=lchain)                      :: filevardiag
      ! File for the diagnostic variables
      LOGICAL               :: l_out_subs_diag 
      ! Boolean for the use of diagnostic variables
#  ifdef MUSTANG
      LOGICAL               :: l_out_subs_diag_sed
      ! Boolean for the use of diagnostics in the sediment
#  endif MUSTANG

      INTEGER,ALLOCATABLE,DIMENSION(:),PUBLIC :: irk_diag ! (reading order)
      ! Table for the reading order of diagnostic variables. This will help when 
      ! interacting with the hydro model which may have a different order for 
      ! Saving them as tracers.

      CHARACTER(LEN=lchain),ALLOCATABLE,DIMENSION(:),PUBLIC :: name_vardiag,         &
                                                               long_name_vardiag,    &
                                                               standard_name_vardiag,&
                                                               unit_vardiag
      ! names, attributes and units of diagnostic variables
 
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:),PUBLIC :: valid_min_vardiag, valid_max_vardiag
      ! Numerical attributes of diagnostic variables

     !=====================================================================
     !  Number of diagnostic variables and storage tables
     !=====================================================================

      INTEGER, PUBLIC :: ndiag_1d, ndiag_2d, ndiag_2d_sed, ndiag_3d, ndiag_3d_wat, &
                         ndiag_3d_sed, ndiag_tot
                         ! Number of diagnostics with different types of dimensions. 
                         ! Following the order of variables : 1d (time), 2d(horizontal),
                         ! 3d (without further precision),
                         ! 3d (vertical dimension, horizontal dimensions),
                         ! 3d (horizontal dimensions, vertical dimension), 
                         ! total number of diagnostics

      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:)      ,PUBLIC :: diag_1d     
      !Table for storing diagnostics variables of dimensions (index of var, time)
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: diag_2d     
      !Table for storing diagnostics variables of dimensions (index of var, horizontal dimensions)
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_wat
      !Table for storing diagnostics variables of dimensions 
      !(index of var, vertical dimension, horizontal dimensions)

#  ifdef MUSTANG
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_sed
      !Table for storing diagnostics variables of dimensions 
      ! (index of var, horizontal dimensions, vertical dimension)
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:),PUBLIC   :: diag_2d_sed
      !Table for storing diagnostics variables of dimensions 
      ! (index of var, horizontal dimensions)
#  endif /* MUSTANG */

     !=====================================================================
     !  Variables that apparently are declared in MUSTANG,
     !  I really should investigate on those
     !=====================================================================

#  if ! defined MUSTANG
      LOGICAL :: l_testcase=.FALSE.
      REAL(kind=rsh)  ,PARAMETER :: valmanq=999.0
      REAL(kind=riosh),PARAMETER :: rg_valmanq_io=999.0
      REAL(kind=rlg)  ,PARAMETER :: epsilon=0.00000001
#  endif /* MUSTANG */

 
     !*************************************************************************!
     !*************************************************************************!
     !******************** Variables from the hydro model *********************!
     !*************************************************************************!
     !*************************************************************************!

     !=====================================================================
     !  Height of the water column variables
     !=====================================================================

      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE  :: THICKLAYERWC,THICKLAYERWW 
      ! Thickness of the water column (?) and of the wave layer (?), I am not sure
#  if ! defined MUSTANG
      REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: TOTAL_WATER_HEIGHT 
      ! Total water height of the column
#  endif /* MUSTANG */

     !=====================================================================
     !  Bottom current variables
     !=====================================================================

#  if defined BLOOM && defined key_benthos
      REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: BOTTCURRENTBIO 
      ! Current in the bottom layer
#  endif /* BLOOM && key_benthos */

     !=====================================================================
     !  Temperature and salinity variables
     !=====================================================================

      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:)               :: SAL_BIOLink,TEMP_BIOLink
      ! Salinity and temperature in the water column

     !*************************************************************************!
     !*************************************************************************!
     !*************** Variables for the model helping routines ****************!
     !********************* Or sediment model *********************************!
     !*************************************************************************!
     !*************************************************************************!

     !=====================================================================
     !  Sinking velocity of tracers
     !=====================================================================

      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:,:)             :: WS_BIOLink 
      ! Settling velocities, it depends on particles density and vertical currents

     !=====================================================================
     !  Radiation variables
     !=====================================================================

     ! Photosynthetic available radiations

#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer
      ! Photosynthetic Available Radiation in the top layer
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: PAR_avg_layer_phyto 
      ! Photosynthetic Available Radiation in a lyaer concerned with phytoplankton,
      ! I am not sure
#    if defined PEPTIC
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer_day 
      ! Photosynthetic Available Radiation in the top layer per day 
!      REAL(KIND=rsh),DIMENSION(:,:)  ,ALLOCATABLE :: light_ave_daily,light_integ
      ! Light averaged per day and total light ( per day ? I am not sure)
#    endif /* PEPTIC */
#  endif /* BIOLink_PAR_eval */
     
     ! Extinction and attenuation

#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EXTINCTION_RAD
      ! Coefficient of light extinction of the water column, I am not sure
#  endif /* BIOLink_PAR_eval */

#  ifdef BLOOM
      REAL(KIND=rsh),DIMENSION(:,:,:,:,:)  ,ALLOCATABLE :: extinction_tab
      ! Table of light extinction, but I do not entirely understand the dimensions
      REAL(KIND=rsh),DIMENSION(:,:,:)  ,ALLOCATABLE :: extinction_ave4d,extinction_aveh
      ! Extinction averaged on four days and per hour
      REAL(KIND=rsh)                                :: t_cum_extinctionh
      ! Extinction cumulated on time
#  endif /* BLOOM */

     ! Variables that affect the extinction

      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: SPMTOT_MGL
      ! Concentration of suspended matter

#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BIOLink_chloro
      ! Chlorophyll in the water column, I am not sure
#  endif /* BIOLink_PAR_eval */

#  ifdef BLOOM
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_diat,effetlumiere_day_dino,effetlumiere_day_nano
      ! Integration over one day of the effect of light on the phytoplankton growth ( from midnight
      !  to midnight) of diatoms, dinoflagellates and nanoplankton
#    ifdef key_phaeocystis
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_phaeocystis
      ! Effect of light on the growth of phaeocystis
#    endif /* key_phaeocystis */
#    ifdef key_karenia
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_karenia
      ! Effect of light on the growth of karenia
#    endif /* key_karenia */
#    ifdef key_psnz
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_psnz
      ! Effect of light on the growth of psnz
#    endif /* key_psnz */
#  endif /* BLOOM */

     !=====================================================================
     !  Estimation of the suspended matter via satellite
     !=====================================================================

#  ifdef key_messat
      CHARACTER(LEN=lchain)                       :: filemessat_clim,filemessat_obs
      ! Files of climatology and observation of suspended matter 
      LOGICAL                                     :: l_messat_clim,l_messat_obs 
      ! Boolean for the activation of climatology of suspended mattter climatology and observations
      INTEGER, PUBLIC                             :: idateinf_messat,idatesup_messat
      ! Date of beginning of climatology and satellite observations
      INTEGER, PUBLIC                             :: date_s_annee_messat,idimt_messat,idimt_messat_obs
      ! Temporal variables for the satellite observations, I am note entirely sure
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat
      ! Concentration of suspended matter observed by satellite
      REAL(KIND=rlg),DIMENSION(1:51)              :: t_clim_messat
      ! Time of the satellite observations, it is one per week I guess
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat_obs
      ! Observation of suspended matter via satellite ?
      REAL(KIND=rlg),DIMENSION(:),ALLOCATABLE      :: t_obs_messat
      ! Time of observation of the satellite measurements ?
      CHARACTER(LEN=lchain)                       :: file_mes_obs_path
      ! Path to the file of observations measurements ?
#  endif /* key_messat */

     !=====================================================================
     !  Wind variables
     !=====================================================================

#  if ! defined BULK_FLUX
      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:)                :: WIND_SPEED 
      ! Wind velocity
#  endif /* BULK_FLUX */

     !*************************************************************************!
     !*************************************************************************!
     !************************** Subroutines **********************************!
     !*************************************************************************!
     !*************************************************************************!


CONTAINS

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

      ALLOCATE( BIO_SKSC_FIX (ARRAY_FIXED_SKSC))
      BIO_SKSC_FIX (:,:,:,:)=0.0_rsh

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

      ALLOCATE( BIO_SKSC_FIX (ARRAY_FIXED_SKSC))
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

  !!======================================================================



#endif /* SUBSTANCE && BIOLink */

END MODULE
