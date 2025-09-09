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

      REAL(KIND=rlg),PUBLIC                     :: t_bio,DT_CONSERV_BIOLINK
      REAL(KIND=rlg),PUBLIC                     :: BIO_TIME_STEP          ! effective time step for evaluation bio sources and sinkks
      !Time and time steps variables used in the biological models 
#if !defined ECO3M                                       
      REAL(KIND=rlg)                             :: ECO_TIME_STEP
#endif
      REAL(KIND=rlg)                             :: TIME_START_ECO 
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
      
#if ! defined ECO3M
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: BIO_SINKSOURCES 
      ! For the variables in the water column
      ! ECO3M has it own array because it is a code that can run by itself
#endif /* ECO3M */

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
#  if defined MUSTANG
      LOGICAL               :: l_out_subs_diag_sed
      ! Boolean for the use of diagnostics in the sediment
#  endif /* MUSTANG */

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
                         ndiag_3d_sed, ndiag_tot, ndiag_2d_wat
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
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC :: diag_3d_CROCO
      ! Table for storing diagnostics 3D variables and exporting them in CROCO
      

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

 

CONTAINS


  !!======================================================================



#endif /* SUBSTANCE && BIOLink */

END MODULE
