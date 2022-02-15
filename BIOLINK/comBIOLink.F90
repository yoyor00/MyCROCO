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


      !---------------------------------------------------------------
      !   Declarations of public subroutines and variables from the 
      !   Program itself.
      !---------------------------------------------------------------

      ! Public subroutines
      PUBLIC BIOLink_alloc ! Subroutines for allocating the space of variables in BIOLink
  
      ! Public variables
      INTEGER,PUBLIC        ::  ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK,  &
                                iheure_BIOLINK,iminu_BIOLINK,isec_BIOLINK, &
                                jjulien_BIOLINK ! Variables used in the BLOOM module that 
                                                ! Relies on julian days

      
      REAL(KIND=rlg),PUBLIC                     :: t_bio,BIO_TIME_STEP,DT_CONSERV_BIOLINK
      !Time and time steps variables used in the biological models 
      LOGICAL, PUBLIC                           :: l_bioretro_extinct,l_waterdensity_known
      ! Boolean variable for the activation of some subroutines of BIOLink
      INTEGER, PUBLIC                           :: IVERIF_BIOLINK,JVERIF_BIOLINK
      ! Counter variables for the verification loops
      
#  if defined key_BIOLink_verif_conserv 
      INTEGER, PUBLIC , PARAMETER               :: iscreenlog_conserv=66
      ! Printing variable for the verification of conservativity
#  endif /* key_BIOLink_verif_conserv */
 

   ! ----------------------------------------------------------------------------
   ! Variables used for the diagnostic of the biological models
   ! ----------------------------------------------------------------------------

      INTEGER, PUBLIC :: ndiag_1d, ndiag_2d, ndiag_2d_sed, ndiag_3d, ndiag_3d_wat, &
                         ndiag_3d_sed, ndiag_tot
                         ! Number of diagnostics with different types of dimensions. 
                         ! Following the order of variables : 1d (time), 2d(horizontal),
                         ! 3d (without further precision),
                         ! 3d (vertical dimension, horizontal dimensions),
                         ! 3d (horizontal dimensions, vertical dimension), 
                         ! total number of diagnostics

      LOGICAL               :: l_out_subs_diag 
      ! Boolean for the use of diagnostic variables

      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:)      ,PUBLIC :: diag_1d     
      !Table for storing diagnostics variables of dimensions (index of var, time)
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: diag_2d     
      !Table for storing diagnostics variables of dimensions (index of var, horizontal dimensions)
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_wat
      !Table for storing diagnostics variables of dimensions 
      !(index of var, vertical dimension, horizontal dimensions)

#  ifdef MUSTANG
      LOGICAL               :: l_out_subs_diag_sed
      ! Boolean for the use of diagnostics in the sediment
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_sed
      !Table for storing diagnostics variables of dimensions 
      ! (index of var, horizontal dimensions, vertical dimension)
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:),PUBLIC   :: diag_2d_sed
      !Table for storing diagnostics variables of dimensions 
      ! (index of var, horizontal dimensions)
#  endif /* MUSTANG */


      INTEGER,ALLOCATABLE,DIMENSION(:),PUBLIC :: idimv_r
      ! I am clueless about this variable
      LOGICAL,ALLOCATABLE,DIMENSION(:),PUBLIC :: l_diagBIOLink_out
      ! Booleans for the writing in external file of diagnostic variables of BIOLink
      
      CHARACTER(LEN=lchain),ALLOCATABLE,DIMENSION(:),PUBLIC :: name_vardiag,         &
                                                               long_name_vardiag,    &
                                                               standard_name_vardiag,&
                                                               unit_vardiag
      ! names, attributes and units of diagnostic variables
           
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:),PUBLIC :: valid_min_vardiag, valid_max_vardiag
      ! Numerical attributes of diagnostic variables

      INTEGER,ALLOCATABLE,DIMENSION(:),PUBLIC :: irk_diag ! (reading order)
      ! Table for the reading order of diagnostic variables. This will help when 
      ! interacting with the hydro model which may have a different order for 
      ! Saving them as tracers.
      



      !---------------------------------------------------------------
      !   Declarations of variables in case of missing values or 
      !   for the comparison of floats (epsilon). Those are normally
      !   Declared in MUSTANG. 
      !---------------------------------------------------------------

#  if ! defined MUSTANG

      LOGICAL :: l_testcase=.FALSE.
      REAL(kind=rsh)  ,PARAMETER :: valmanq=999.0
      REAL(kind=riosh),PARAMETER :: rg_valmanq_io=999.0
      REAL(kind=rlg)  ,PARAMETER :: epsilon=0.00000001
#  endif /* MUSTANG */

  
      !---------------------------------------------------------------
      !   Private part of the program, declaration of local variables     
      !---------------------------------------------------------------
      
      REAL(kind=rsh)  ,PARAMETER                :: epsilon_BIOLink=1.e-10
      ! Variable for the comparison of float internal to BIOLink
      REAL(KIND=rlg)                             :: ECO_TIME_STEP,TIME_START_ECO 
      ! Timestep of the hydrodynamical model
      REAL(KIND=rlg)                             :: TIME_BEGIN  
      ! Date of the start of the model run
      CHARACTER(LEN=lchain)                      :: suffix_fileres,filevardiag
      ! File names for the diagnostic variables
      
      !---------------------------------------------------------------
      !   Private part of the program, declaration of local tables
      !   before their allocations     
      !---------------------------------------------------------------
      
      !***********************************************************************
      ! Those are the most important tables of BIOLink, the ones that make   !
      ! The link between the hydro model and the biological model through    !
      ! Sources and sink terms that are computed in the bio model and        !
      ! Transferred to the hydro model                                       !
      !**********************************************************************!
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: BIO_SINKSOURCES 
      ! For the variables in the water column
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: BIO_SKSC_FIX 
      ! For the fixed variables

      !---------------------------------------------------------------
      !  Variables related to light availability     
      !---------------------------------------------------------------
                                      
#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer
      ! Photosynthetic Available Radiation in the top layer
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: PAR_avg_layer_phyto 
      ! Photosynthetic Available Radiation in a lyaer concerned with phytoplankton,
      ! I am not sure
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EXTINCTION_RAD
      ! Coefficient of light extinction of the water column, I am not sure
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BIOLink_chloro
      ! Chlorophyll in the water column, I am not sure
#    if defined PEPTIC
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer_day 
      ! Photosynthetic Available Radiation in the top layer per day 
!      REAL(KIND=rsh),DIMENSION(:,:)  ,ALLOCATABLE :: light_ave_daily,light_integ
      ! Light averaged per day and total light ( per day ? I am not sure)
#    endif /* PEPTIC */
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

   ! ----------------------------------------------------------------------------
   ! Variables relative to light penetration
   ! ----------------------------------------------------------------------------

      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: SPMTOT_MGL
      ! Concentration of suspended matter

#  ifdef BLOOM
      REAL(KIND=rsh),DIMENSION(:,:,:,:,:)  ,ALLOCATABLE :: extinction_tab
      ! Table of light extinction, but I do not entirely understand the dimensions
      REAL(KIND=rsh),DIMENSION(:,:,:)  ,ALLOCATABLE :: extinction_ave4d,extinction_aveh
      ! Extinction averaged on four days and per hour
      REAL(KIND=rsh)                                :: t_cum_extinctionh
      ! Extinction cumulated on time
#  endif /* BLOOM */

   ! ----------------------------------------------------------------------------
   ! Subsection dedicated to the estimation of the suspended matter via satellite
   ! ----------------------------------------------------------------------------

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
      ! Time of the satellite observations, why 51 though ?
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat_obs
      ! Observation of suspended matter via satellite ?
      REAL(KIND=rlg),DIMENSION(:),ALLOCATABLE      :: t_obs_messat
      ! Time of observation of the satellite measurements ?
      CHARACTER(LEN=lchain)                       :: file_mes_obs_path
      ! Path to the file of observations measurements ?
#  endif /* key_messat */



      !---------------------------------------------------------------
      !  Variables related to the water thickness    
      !---------------------------------------------------------------

      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE  :: THICKLAYERWC,THICKLAYERWW 
      ! Thickness of the water column (?) and of the wave layer (?), I am not sure
#  if ! defined MUSTANG
      REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: TOTAL_WATER_HEIGHT 
      ! Total water height of the column
#  endif /* MUSTANG */


      !---------------------------------------------------------------
      ! Variables related to the concentrations of chemical/biological
      ! Species in the water column     
      !---------------------------------------------------------------

#  if !defined ECO3M 
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: WATCONCPOS 
      ! Concentration in the water column, the order is the one of BIOLink for the storage 
      ! which means (species_index,vertical_direction,zonal_direction,meridional_direction)
#  endif /* ECO3M */

      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: FIXCONCPOS  
      ! Concentration of fixed variables, the order is the one of BIOLink for the storage
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BENTCONCPOS 
      ! Concentration of benthic variables, the order is the one of BIOLink for the storage

      !---------------------------------------------------------------
      ! Variables related to physico-chemical conditions (temperature
      ! ,salinity,etc...)  
      !---------------------------------------------------------------


#  if defined BLOOM && defined key_benthos
      REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: BOTTCURRENTBIO 
      ! Current in the bottom layer
#  endif /* BLOOM && key_benthos */
  
      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:)               :: SAL_BIOLink,TEMP_BIOLink
      ! Salinity and temperature in the water column
      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:,:)             :: WS_BIOLink 
      ! Settling velocities, it depends on particles density and vertical currents
#  if ! defined BULK_FLUX
      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:)                :: WIND_SPEED 
      ! Wind velocity
#  endif /* BULK_FLUX */


   ! ----------------------------------------------------------------------------
   ! Internal variables of the biological model used, those are stored in 
   ! external files.
   ! ----------------------------------------------------------------------------

#  if defined BLOOM
#include "combloom.h"

#  elif defined PEPTIC
#include "compeptic.h"

#  elif defined METeOR 
#include "commeteor.h"

#  endif /* Biological models */

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


   !---------------------------------------------------------------------
   ! Here we are going to group the allocation and initialization 
   ! As it will make the code clearer. But it could be even easier
   ! To produce two subroutines, one for the allocation and one for the 
   ! Initialization
   !---------------------------------------------------------------------
  

   !---------------------------------------------------------------------
   ! Allocation of tables related to the height of the water column
   !---------------------------------------------------------------------

#  if ! defined MUSTANG
        ALLOCATE( TOTAL_WATER_HEIGHT(PROC_IN_ARRAY_m2p2))
#  endif /* MUSTANG */
        
        ALLOCATE( THICKLAYERWC(NB_LAYER_WAT,PROC_IN_ARRAY))
        ALLOCATE( THICKLAYERWW(NB_LAYER_WAT,PROC_IN_ARRAY))
        ALLOCATE( SPMTOT_MGL(NB_LAYER_WAT,PROC_IN_ARRAY))
        
   !---------------------------------------------------------------------
   ! Initialization of the table of water height to dummy values
   !---------------------------------------------------------------------
        THICKLAYERWC(:,:,:)=0.0_rsh
        THICKLAYERWW(:,:,:)=0.0_rsh
        SPMTOT_MGL(:,:,:)=0.0_rsh

   !---------------------------------------------------------------------
   ! Allocation of tables related to Photosynthetic Available Radiation
   ! And light penetration/absorption
   !---------------------------------------------------------------------

#  if defined  BIOLink_PAR_eval

      ALLOCATE( EXTINCTION_RAD(NB_LAYER_WAT,PROC_IN_ARRAY) )
      ALLOCATE( PAR_top_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )
      ALLOCATE( BIOLink_chloro(NB_LAYER_WAT,PROC_IN_ARRAY) )  

#    if defined PEPTIC
      ALLOCATE( PAR_avg_layer_phyto(1,NB_LAYER_WAT,PROC_IN_ARRAY) )
      ALLOCATE( PAR_top_layer_day(NB_LAYER_WAT,PROC_IN_ARRAY))

#    elif defined METeOR
      ALLOCATE( Flimrad_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )

#    endif /* PEPTIC/METEOR */
#  endif /* BIOLink_PAR_eval */

#  if defined BLOOM
      ALLOCATE(extinction_ave4d(NB_LAYER_WAT,PROC_IN_ARRAY))
      ALLOCATE(extinction_aveh(NB_LAYER_WAT,PROC_IN_ARRAY))
      ALLOCATE(extinction_tab(4,24,NB_LAYER_WAT,PROC_IN_ARRAY))
      
      ALLOCATE( effetlumiere_day_diat(NB_LAYER_WAT,PROC_IN_ARRAY) )
      ALLOCATE( effetlumiere_day_dino(NB_LAYER_WAT,PROC_IN_ARRAY) )
      ALLOCATE( effetlumiere_day_nano(NB_LAYER_WAT,PROC_IN_ARRAY) )
      
#    ifdef key_karenia
      ALLOCATE( effetlumiere_day_karenia(NB_LAYER_WAT,PROC_IN_ARRAY) )
#    endif /* key_karenia */

#    ifdef key_psnz
      ALLOCATE( effetlumiere_day_psnz(NB_LAYER_WAT,PROC_IN_ARRAY) )
#    endif /* key_psnz */

#    ifdef key_phaeocystis
      ALLOCATE( effetlumiere_day_phaeocystis(NB_LAYER_WAT,PROC_IN_ARRAY) )
#    endif /* key_phaeocystis */
#  endif /* BLOOM */


   !---------------------------------------------------------------------
   ! Initialization of tables related to the height of the water column
   ! to dummy values. Some tables do are allocated but not initialized
   ! because they are initialized directly in the code they belong to
   !---------------------------------------------------------------------
#  if defined BIOLink_PAR_eval
      EXTINCTION_RAD(:,:,:)=0.0_rsh
      PAR_top_layer(:,:,:)=0.0_rsh
      BIOLink_chloro(:,:,:)=0.0_rsh

#    if defined PEPTIC
      PAR_top_layer_day(:,:,:)=0.0_rsh
      PAR_avg_layer_phyto(1,:,:,:)=0.0_rsh
#    endif /* PEPTIC */

#  endif /* BIOLink_PAR_EVAL */

#  ifdef BLOOM
      extinction_ave4d(:,:,:)=0.0_rsh
      extinction_aveh(:,:,:)=0.0_rsh
      extinction_tab(:,:,:,:,:)=0.0_rsh
      t_cum_extinctionh=0.0_rsh
      ihour_previous=0

      !ALLOCATE(phytomoy(PROC_IN_ARRAY))
      !ALLOCATE(phytocarre(PROC_IN_ARRAY))
      !phytomoy(:,:)=0.0_rsh
      
      effetlumiere_day_diat(:,:,:)=0.0_rsh
      effetlumiere_day_dino(:,:,:)=0.0_rsh
      effetlumiere_day_nano(:,:,:)=0.0_rsh

#    ifdef key_karenia
      effetlumiere_day_karenia(:,:,:)=0.0_rsh
#    endif /* key_karenia */

#    ifdef key_psnz
      effetlumiere_day_psnz(:,:,:)=0.0_rsh
#    endif /* key_psnz */

#    ifdef key_phaeocystis
      effetlumiere_day_phaeocystis(:,:,:)=0.0_rsh
#    endif /* key_phaeocystis */
#  endif /* BLOOM */


   !---------------------------------------------------------------------
   ! Allocation of tables related to physico-chemical conditions
   !---------------------------------------------------------------------

      ALLOCATE( TEMP_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
      ALLOCATE( SAL_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
      ALLOCATE( WS_BIOLink(NB_LAYER_WAT,nv_adv,PROC_IN_ARRAY))
      ALLOCATE( BIO_SKSC_FIX (ARRAY_FIXED_SKSC))

#  if defined BLOOM && defined key_benthos

      ALLOCATE(BOTTCURRENTBIO(PROC_IN_ARRAY)) 

#  endif /* BLOOM && key_benthos */


#  if ! defined BULK_FLUX
      ALLOCATE( WIND_SPEED(PROC_IN_ARRAY))
#  endif /* BULK_FLUX */


   !---------------------------------------------------------------------
   ! Initialization of tables related to physico-chemical conditions
   ! To dummy values
   !---------------------------------------------------------------------

      TEMP_BIOLink(:,:,:)=0.0_rsh
      SAL_BIOLink(:,:,:)=0.0_rsh
      WS_BIOLink(:,:,:,:)=0.0_rsh
      BIO_SKSC_FIX (:,:,:,:)=0.0_rsh

#  if ! defined BULK_FLUX
      WIND_SPEED(:,:)=0.0_rsh
#  endif /* BULK_FLUX */

   !---------------------------------------------------------------------
   ! Allocation of tables related to the concentration of tracers
   ! and the sink/source terms. The exception for ECO3M are because the 
   ! Allocation and initialization are internal for this model.
   !---------------------------------------------------------------------

      ALLOCATE( BIO_SINKSOURCES(ARRAY_SINKSOURCES) )
      ALLOCATE(FIXCONCPOS(nv_fix,NB_LAYER_WAT,PROC_IN_ARRAY))

#  if ! defined ECO3M  
      ALLOCATE(WATCONCPOS(nv_adv,NB_LAYER_WAT,PROC_IN_ARRAY))
#  endif /* ECO3M */

#  ifdef key_benthic
      ALLOCATE(BENTCONCPOS(nv_bent,PROC_IN_ARRAY))
#  endif /* key_benthic */

   !---------------------------------------------------------------------
   ! Initialization of tables related to the concentration of tracers
   ! and the sink/source terms to dummy values
   !---------------------------------------------------------------------
    
      BIO_SINKSOURCES(:,:,:,:)=0.0_rsh
      FIXCONCPOS(:,:,:,:)=0.0_rsh

#  if ! defined ECO3M  
      WATCONCPOS(:,:,:,:)=0.0_rsh
#  endif /* ECO3M */

#  ifdef key_benthic
      BENTCONCPOS(:,:,:,:)=0.0_rsh
#  endif /* key_benthic */

   !---------------------------------------------------------------------
   ! Allocation of the table related to suspended matter measured by
   ! satellite
   !---------------------------------------------------------------------
  
#  ifdef key_messat
      ALLOCATE(messat(PROC_IN_ARRAY,51))
#  endif /* key_messat */

   !---------------------------------------------------------------------
   ! Initialization of the table related to suspended matter measured by
   ! satellite to dummy value
   !---------------------------------------------------------------------
  
#  ifdef key_messat
      messat(:,:,:)=0.0_rsh
#  endif /* key_messat */

   !---------------------------------------------------------------------
   ! Allocation of the tables related to oysters
   !---------------------------------------------------------------------


#  if defined BLOOM
#    ifdef key_oyster_SFG
      ALLOCATE(tpostpontecoq(PROC_IN_ARRAY))
      ALLOCATE(tpostpontegam(PROC_IN_ARRAY))
#    endif /* key_oyster_SFG */

#    if defined key_oyster_benthos || defined key_oyster_DEB || defined key_oyster_SFG
      ALLOCATE(nbhuitre(PROC_IN_ARRAY))
      ALLOCATE(hautable(PROC_IN_ARRAY))
#    endif /* key_oyster_benthos */
#  endif /* BLOOM */

   !---------------------------------------------------------------------
   ! Initialization of the tables related to oysters to dummy values
   !---------------------------------------------------------------------


#  if defined BLOOM
#    if defined key_oyster_benthos || defined key_oyster_DEB || defined key_oyster_SFG
      hautable(:,:)=0.0_rsh
      nbhuitre(:,:)=0.0_rsh
#    endif /* key_oyster_benthos */
#  endif /* BLOOM */

END SUBROUTINE  BIOLink_alloc

  !!======================================================================



#endif /* SUBSTANCE && BIOLink */

END MODULE
