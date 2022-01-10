!---------------------------------------------------------------------------
!
                     MODULE comBIOLink
!
!---------------------------------------------------------------------------
#include "cppdefs.h"

#if defined SUBSTANCE && defined BIOLink

  !!======================================================================
  !!                   ***  MODULE  comBIOLink  ***
  !! Ocean dynamics Bio :  declare and initialize all common variables related to BIOLink module
  !!                            + routines related to radiation, extinction, PAR etc...
  !!                            + routine related to water column
  !!
  !!   History :
  !!    !  2019-08 (B. Thouvenin) issued from modules BLOOM (old ECOMARS/ key_biolo) and peptic for portability adaptation
  !!      !            (V. Garnier, M. Sourisseau, A. Menesguen, P. Cugier, A. Arancio)

  !!======================================================================

#ifdef key_MARS
#include "toolcpp.h"
#include "coupleur_dimhydro_BIOLink.h"
#endif
#include "coupleur_define_BIOLink.h"

  !! * Modules used
#if ! defined key_MARS
   USE comsubstance
#endif

  IMPLICIT NONE

  !! * Accessibility
  PUBLIC BIOLink_alloc
  
  !! * Shared or public module variables

 INTEGER,PUBLIC        ::  ijour_BIOLINK,imois_BIOLINK,ian_BIOLINK,  &
                           iheure_BIOLINK,iminu_BIOLINK,isec_BIOLINK, &
                           jjulien_BIOLINK

#if ! defined key_MARS && ! defined MUSTANG
 !-------------------------------------------------------------------
 !  definition of rsh, rlg...
 !-------------------------------------------------------------------
   !INTEGER,PARAMETER :: riosh=4,riolg=8,rlg=8,rsh=8
 !   + fixed data used by MUSTANG but specific to MARS and therefore not necessarily known for another model
   LOGICAL :: l_testcase=.FALSE.
   REAL(kind=rsh)  ,PARAMETER :: valmanq=999.0
   REAL(kind=riosh),PARAMETER :: rg_valmanq_io=999.0
   REAL(kind=rlg)  ,PARAMETER :: epsilon=0.00000001
!!!!
#endif
  
  REAL(kind=rsh)  ,PARAMETER                :: epsilon_BIOLink=1.e-10

  REAL(KIND=rlg),PUBLIC                     :: t_bio,           &   ! time of next bio dynamic step in water column 
                                               BIO_TIME_STEP,           &   ! effective time step for evaluation bio sources and sinkks
                                               DT_CONSERV_BIOLINK           ! time step for conservativity verification in one mesh

  LOGICAL, PUBLIC                           :: l_bioretro_extinct,l_waterdensity_known
  INTEGER, PUBLIC                           :: IVERIF_BIOLINK,JVERIF_BIOLINK

#if defined key_BIOLink_verif_conserv 
  INTEGER, PUBLIC , PARAMETER               :: iscreenlog_conserv=66
#endif

  REAL(KIND=rlg)                             :: ECO_TIME_STEP,TIME_START_ECO  ! for Hydro
#if ! defined key_MARS
  REAL(KIND=rlg)                             :: TIME_BEGIN  ! declared in BIOLink for CROCO (=tdeb in MARS)
#endif  
  CHARACTER(LEN=lchain)                      :: suffix_fileres,filevardiag

!!! avec MARS : ces termes sont declares dans comsubstance et utilises dans le code hydro en les incluant dans la resolution des equations de transport
! attention avec d autres codes hydro, peuvent etre deja declares et alloues dans le code hydro et donc supprimer la ligne
!                                      ou ne pas exister dans le code hyrdo car la resolution se fait dans BIOLink - dc/dt= sinksources (pas fractionnaire)
#if ! defined key_MARS
  REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: BIO_SINKSOURCES ! variation temporelle des concentrations dc/dt due aux transformations bio
                                                                   ! = dcdt_adv si interne a BIOLink
                                                                   ! mais attention si externe a BIOLink, peut etre deja declare dans modele hote
                                                                   ! dans MARS, declaration dans comsubstance
#endif
!!!!!!!!!!!!!
  ! si calcul de PAR via BIOLink (sinon, PAR calcule par module BIO a partir de la radiation solaire incidente)
#if defined  BIOLink_PAR_eval
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer
  REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: PAR_avg_layer_phyto !
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EXTINCTION_RAD !
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BIOLink_chloro  ! chlorophylle (attention unity)
#if defined PEPTIC
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer_day 
!  REAL(KIND=rsh),DIMENSION(:,:)  ,ALLOCATABLE :: light_ave_daily,light_integ
#endif
#endif

  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE  :: THICKLAYERWC,THICKLAYERWW 
#if ! defined MUSTANG
  REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: TOTAL_WATER_HEIGHT 
#endif
# if !defined ECO3M 
  REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: WATCONCPOS  ! positive advected conc. with index order of BIOLink module
#endif
  REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: FIXCONCPOS  ! positive fixed conc. with index order of BIOLink module
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BENTCONCPOS  ! positive benthic conc. with index order of BIOLink module

#if defined BLOOM && defined key_benthos
  REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: BOTTCURRENTBIO 
#endif
#if ! defined key_MARS
  REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:)               :: SAL_BIOLink,TEMP_BIOLink
  REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:,:)             :: WS_BIOLink,BIO_SKSC_FIX 
!! pour CROCO if not BULK_FLUX : wind_speed not known and evaluated in BIOLink
#if ! defined BULK_FLUX
  REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:)                :: WIND_SPEED 
#endif
#endif

#ifdef BLOOM
   ! integration sur 1 jour (minuit_minuit) de l effet lumiere sur la croissance du phytoplancton
   REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_diat,effetlumiere_day_dino,effetlumiere_day_nano
#ifdef key_phaeocystis
   REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_phaeocystis
#endif
#ifdef key_karenia
   REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_karenia
#endif
#ifdef key_psnz
   REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_psnz
#endif
#endif
   ! ----------------------------------------------------------------------------
   ! PARAMETRES RELATIFS A PENETRATION DE LA LUMIERE
   ! ----------------------------------------------------------------------------

   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: SPMTOT_MGL
#ifdef BLOOM
   REAL(KIND=rsh),DIMENSION(:,:,:,:,:)  ,ALLOCATABLE :: extinction_tab
   REAL(KIND=rsh),DIMENSION(:,:,:)  ,ALLOCATABLE :: extinction_ave4d,extinction_aveh
   REAL(KIND=rsh)                                :: t_cum_extinctionh
#endif

#ifdef key_messat
   CHARACTER(LEN=lchain)                       :: filemessat_clim,filemessat_obs
   LOGICAL                                     :: l_messat_clim,l_messat_obs
   INTEGER, PUBLIC                             :: idateinf_messat,idatesup_messat
   INTEGER, PUBLIC                             :: date_s_annee_messat,idimt_messat,idimt_messat_obs
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat
   REAL(KIND=rlg),DIMENSION(1:51)              :: t_clim_messat
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat_obs
   REAL(KIND=rlg),DIMENSION(:),ALLOCATABLE      :: t_obs_messat
   CHARACTER(LEN=lchain)                       :: file_mes_obs_path
#endif

   ! ----------------------------------------------------------------------------
   ! PARAMETRES RELATIFS AU MODULE BIO CHOISI
   ! ----------------------------------------------------------------------------

#if defined BLOOM

!!!!!! module option BLOOM (ECOMARS)   !!!!!!!!!!!!!!!!!!!!!!

#include "combloom.h"

#elif defined PEPTIC

!!!!!! module option peptic    !!!!!!!!!!!!!!!!!!!!!!

#include "compeptic.h"

#elif defined METeOR 

!!!!!! module option meteor (contaminants : metals, Organic, batcerian, radioactiv    !!!!!!!!!!!!!!!!!!!!!!

#include "commeteor.h"

#endif

   ! ----------------------------------------------------------------------------
   ! PARAMETRES RELATIFS AUX VARIABLES DIAGNOSTIQUES
   ! ----------------------------------------------------------------------------

  ! number of diagnostic variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, PUBLIC :: ndiag_1d,     &  ! 1 dimension (time)
                     ndiag_2d,     &  ! 2 dimensions (i,j)
                     ndiag_2d_sed, &  ! 3 d (i,j)
                     ndiag_3d,     &  ! 3 d
                     ndiag_3d_wat, &  ! 3 d (k,i,j)
                     ndiag_3d_sed, &  ! 3 d (i,j,k)
                     ndiag_tot
#if ! defined key_MARS
 LOGICAL               :: l_out_subs_diag   ! dans MARS, declaration dans comvars2d
#endif
  ! arrays of diagnostic variables
  REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:)      ,PUBLIC :: diag_1d     !(index of var, time)
  REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: diag_2d     ! (index of var, i, j)
  REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_wat !(index of var,k,i,j)
#ifdef MUSTANG
  LOGICAL               :: l_out_subs_diag_sed
  REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_sed !(index of var,i,j,k)
  REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:),PUBLIC   :: diag_2d_sed !(index of var,i,j,k)
#endif


  INTEGER,ALLOCATABLE,DIMENSION(:),PUBLIC :: idimv_r
  LOGICAL,ALLOCATABLE,DIMENSION(:),PUBLIC :: l_diagBIOLink_out

  ! names and attributes of diagnostic variables
  CHARACTER(LEN=lchain),ALLOCATABLE,DIMENSION(:),PUBLIC :: name_vardiag,         &
                                                            long_name_vardiag,    &
                                                            standard_name_vardiag

  ! units of diagnostic variables
  CHARACTER(LEN=lchain),ALLOCATABLE,DIMENSION(:),PUBLIC :: unit_vardiag

  ! rank of diagnostic variables
  INTEGER,ALLOCATABLE,DIMENSION(:),PUBLIC :: irk_diag ! (reading order)

  ! values of attributes of diagnostic variables
  REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:),PUBLIC :: valid_min_vardiag, valid_max_vardiag



  CONTAINS

  !!======================================================================
  SUBROUTINE BIOLink_alloc(nv_state)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_alloc ***
  !&E
  !&E ** Purpose : allocation common variables 
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
     INTEGER ,INTENT(IN)                 :: nv_state


  !! * Local declarations

  !!----------------------------------------------------------------------
  !! * Executable part
  
#if ! defined MUSTANG
   ALLOCATE( TOTAL_WATER_HEIGHT(PROC_IN_ARRAY_m2p2))
#endif
   ALLOCATE( THICKLAYERWC(NB_LAYER_WAT,PROC_IN_ARRAY))   !centered around C,T,Sal
   ALLOCATE( THICKLAYERWW(NB_LAYER_WAT,PROC_IN_ARRAY))   !centered around W, interface
   ALLOCATE( SPMTOT_MGL(NB_LAYER_WAT,PROC_IN_ARRAY))
   THICKLAYERWC(:,:,:)=0.0_rsh
   THICKLAYERWW(:,:,:)=0.0_rsh
   SPMTOT_MGL(:,:,:)=0.0_rsh

#if defined  BIOLink_PAR_eval
   ALLOCATE( EXTINCTION_RAD(NB_LAYER_WAT,PROC_IN_ARRAY) )         ! extinction
   ALLOCATE( PAR_top_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )         ! alumplafond
   ALLOCATE( BIOLink_chloro(NB_LAYER_WAT,PROC_IN_ARRAY) )     ! chlorophylle concentration for radiation attenuation 
                                                      ! and PAR estimation (unity depends on bio module)
   EXTINCTION_RAD(:,:,:)=0.0_rsh
   PAR_top_layer(:,:,:)=0.0_rsh
   BIOLink_chloro(:,:,:)=0.0_rsh

#if defined PEPTIC
   ALLOCATE( PAR_avg_layer_phyto(1,NB_LAYER_WAT,PROC_IN_ARRAY) )     !PAR total : indice 1
   ALLOCATE( PAR_top_layer_day(NB_LAYER_WAT,PROC_IN_ARRAY))
   PAR_top_layer_day(:,:,:)=0.0_rsh
   PAR_avg_layer_phyto(1,:,:,:)=0.0_rsh
#elif defined BLOOM
   ! no need PAR_avg_layer_phyto fro tis module. Average light effect is estimated into the module for each phyto species
#endif

#if defined METeOR
   ALLOCATE( Flimrad_layer(0:NB_LAYER_WAT,PROC_IN_ARRAY) )         ! alumplafond with RAD=1
#endif
#endif

#if defined BLOOM && defined key_benthos
  ALLOCATE(BOTTCURRENTBIO(PROC_IN_ARRAY)) 
#endif

#if ! defined key_MARS
   ! salinity and temperature arrays converted for BIOLink (if index different from in hydro model)
   ALLOCATE( TEMP_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
   ALLOCATE( SAL_BIOLink(NB_LAYER_WAT,PROC_IN_ARRAY))
   ALLOCATE( WS_BIOLink(NB_LAYER_WAT,nv_adv,PROC_IN_ARRAY))
   ALLOCATE( BIO_SKSC_FIX (ARRAY_FIXED_SKSC))
   TEMP_BIOLink(:,:,:)=0.0_rsh
   SAL_BIOLink(:,:,:)=0.0_rsh
   WS_BIOLink(:,:,:,:)=0.0_rsh
   BIO_SKSC_FIX (:,:,:,:)=0.0_rsh
!! pour CROCO if not BULK_FLUX : wind_speed not known and evaluated in BIOLink
#if ! defined BULK_FLUX
   ALLOCATE( WIND_SPEED(PROC_IN_ARRAY))
   WIND_SPEED(:,:)=0.0_rsh
#endif

#endif


!!! allocation uniquement si pas alloue dans le code hydro ou dans comsubstance
#if ! defined key_MARS
   ALLOCATE( BIO_SINKSOURCES(ARRAY_SINKSOURCES) )
      BIO_SINKSOURCES(:,:,:,:)=0.0_rsh
#endif
#if ! defined ECO3M  
!    /* internal allocate in ECO3M */
   ALLOCATE(WATCONCPOS(nv_adv,NB_LAYER_WAT,PROC_IN_ARRAY))
   WATCONCPOS(:,:,:,:)=0.0_rsh
#endif
   ALLOCATE(FIXCONCPOS(nv_fix,NB_LAYER_WAT,PROC_IN_ARRAY))
   FIXCONCPOS(:,:,:,:)=0.0_rsh
#ifdef key_benthic
   ALLOCATE(BENTCONCPOS(nv_bent,PROC_IN_ARRAY))
   BENTCONCPOS(:,:,:,:)=0.0_rsh
#endif
     
#ifdef BLOOM
   ALLOCATE(extinction_ave4d(NB_LAYER_WAT,PROC_IN_ARRAY))
   ALLOCATE(extinction_aveh(NB_LAYER_WAT,PROC_IN_ARRAY))
   ALLOCATE(extinction_tab(4,24,NB_LAYER_WAT,PROC_IN_ARRAY))
   extinction_ave4d(:,:,:)=0.0_rsh
   extinction_aveh(:,:,:)=0.0_rsh
   extinction_tab(:,:,:,:,:)=0.0_rsh
   t_cum_extinctionh=0.0_rsh
   ihour_previous=0
   !ALLOCATE(phytomoy(PROC_IN_ARRAY))
   !ALLOCATE(phytocarre(PROC_IN_ARRAY))
   !phytomoy(:,:)=0.0_rsh
#endif


#ifdef key_messat
   ALLOCATE(messat(PROC_IN_ARRAY,51))
   messat(:,:,:)=0.0_rsh
#endif


#ifdef BLOOM
   !allocate effet lumiere integre sur la journee
   ALLOCATE( effetlumiere_day_diat(NB_LAYER_WAT,PROC_IN_ARRAY) )
   ALLOCATE( effetlumiere_day_dino(NB_LAYER_WAT,PROC_IN_ARRAY) )
   ALLOCATE( effetlumiere_day_nano(NB_LAYER_WAT,PROC_IN_ARRAY) )
   effetlumiere_day_diat(:,:,:)=0.0_rsh
   effetlumiere_day_dino(:,:,:)=0.0_rsh
   effetlumiere_day_nano(:,:,:)=0.0_rsh
#ifdef key_karenia
   ALLOCATE( effetlumiere_day_karenia(NB_LAYER_WAT,PROC_IN_ARRAY) )
   effetlumiere_day_karenia(:,:,:)=0.0_rsh
#endif
#ifdef key_psnz
   ALLOCATE( effetlumiere_day_psnz(NB_LAYER_WAT,PROC_IN_ARRAY) )
   effetlumiere_day_psnz(:,:,:)=0.0_rsh
#endif

#ifdef key_phaeocystis
   ALLOCATE( effetlumiere_day_phaeocystis(NB_LAYER_WAT,PROC_IN_ARRAY) )
   effetlumiere_day_phaeocystis(:,:,:)=0.0_rsh
#endif

#ifdef key_oyster_SFG
   ALLOCATE(tpostpontecoq(PROC_IN_ARRAY))
   ALLOCATE(tpostpontegam(PROC_IN_ARRAY))
#endif

#if defined key_oyster_benthos || defined key_oyster_DEB || defined key_oyster_SFG
    ALLOCATE(nbhuitre(PROC_IN_ARRAY))
!    ALLOCATE(nbhuitre2(PROC_IN_ARRAY))
!    ALLOCATE(nbhuitre3(PROC_IN_ARRAY))
    ALLOCATE(hautable(PROC_IN_ARRAY))
    hautable(:,:)=0.0_rsh
    nbhuitre(:,:)=0.0_rsh
!    nbhuitre2(:,:)=0.0_rsh
!    nbhuitre3(:,:)=0.0_rsh
#endif

#endif


END SUBROUTINE  BIOLink_alloc

  !!======================================================================



#endif

END MODULE
