! $Id: comsubstance.F90  2018-08- bthouvenin 
#include "cppdefs.h"

 MODULE comsubstance


   !!======================================================================
   !!                   ***  MODULE comsubstance   ***
   !! 
   !!  ** History :
   !!        ! 2018-08 (B. Thouvenin) adapted from MARS for key SUBSTANCE   !! 
   !!======================================================================


   !! * Modules used
   USE module_substance

#if defined SUBSTANCE
   !-------------------------------------------------------------------
   !  definition of rsh, rlg...
   !-------------------------------------------------------------------
   !INTEGER,PARAMETER :: riosh=4,riolg=8,rlg=8,rsh=8
   !INTEGER,PARAMETER :: lchain=200
             
   ! ----------------------------------------------------------------------------
   ! GESTION DU TYPE DE VARIABLE
   ! ----------------------------------------------------------------------------
     ! USE module_substance,  ONLY : lchain

      INTEGER           :: nv_grav,nv_sand,nv_mud,nspb
      INTEGER           :: igrav1,igrav2,isand1,isand2,imud1,imud2
      INTEGER           :: nv_ncp,nv_dis,nv_sorb,nv_fix,nvp,nvpc
      INTEGER           :: nv_adv,nv_state,nv_tot,nv_bent
          ! number of different maximum settling velocities
       !INTEGER      :: nb_ws_max  ! a group per different value of the maximum velocities
      !INTEGER,DIMENSION(:),ALLOCATABLE :: itypv

        ! indices des particulaires constitutives associees aux particulaires adsorbees, 
      INTEGER,DIMENSION(:),ALLOCATABLE           :: irkm_var_assoc

      ! ----------------------------------------------------------------------------
      ! GESTION CORRESPONDANCE ORDRE DES VARIABLES LUES ET VARIABLES DANS MODELE
      ! ----------------------------------------------------------------------------

      ! irk_mod : rang (position) dans le modele de la variable dans tableaux gerant substances
      ! irk_fil : rang (position) dans le fichier de donnees de la variable
  
      !INTEGER       ,DIMENSION(:),ALLOCATABLE    :: irk_mod,irk_fil
   

      !  gestion des variables forcantes variables with time
      !LOGICAL                                  :: l_driv_cst
      !REAL(KIND=rlg)                           :: dt_driv,t_driv

      ! ----------------------------------------------------------------------------
      ! MATRIX
      ! ----------------------------------------------------------------------------
   
      REAL(KIND=rsh),DIMENSION(:)  ,ALLOCATABLE :: cini_wat,cini_air,cobc_wat
      REAL(KIND=rsh),DIMENSION(:)  ,ALLOCATABLE :: typart,typdiss
      REAL(KIND=rsh),DIMENSION(:)  ,ALLOCATABLE :: sub_flx_atm,cv_rain
      CHARACTER(LEN=lchain),DIMENSION(:) ,ALLOCATABLE :: obc_cv_name
      CHARACTER(LEN=lchain),DIMENSION(:) ,ALLOCATABLE :: init_cv_name
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE           :: ws_part
      REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE       :: ws_free_min, ws_free_max
   
      ! REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: phicon
      ! REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: phicon_drycell
      ! REAL,DIMENSION(:,:,:),ALLOCATABLE   :: phicon_atm

#ifdef MUSTANG
      REAL(KIND=rsh),DIMENSION(:,:)    ,ALLOCATABLE :: ws_free_para, ws_hind_para
      INTEGER,DIMENSION(:)   ,ALLOCATABLE :: ws_free_opt,ws_hind_opt
      REAL(KIND=rsh),DIMENSION(:)      ,ALLOCATABLE :: tocd      
      REAL(KIND=rsh), DIMENSION(:)     ,ALLOCATABLE :: cini_sed_r
      REAL(KIND=rsh),DIMENSION(:)      ,ALLOCATABLE :: diam_sed,ros

      LOGICAL, PUBLIC,DIMENSION(:),ALLOCATABLE                 :: l_subs2D
      INTEGER,PUBLIC,DIMENSION(:),ALLOCATABLE                  :: irk_fil
      REAL(KIND=rsh),DIMENSION(:), PUBLIC  ,ALLOCATABLE        :: unit_modif_mudbio_N2dw
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
   INTEGER                          :: ibedload1,ibedload2
#endif
#endif

   ! Management of group of substances with settling velocity
   ! a group corresponds to a value of the maximum settling velocities
   !INTEGER       ,DIMENSION(:,:),ALLOCATABLE     :: iv_ws_max ! index of variable with settling velocity
   !INTEGER       ,DIMENSION(:),ALLOCATABLE       :: ndt_part  ! sub-time step relative to a max settling velocity
   !INTEGER       ,DIMENSION(:),ALLOCATABLE       :: nv_ws_max ! number of variables with the same max settling velocity
   !INTEGER       ,DIMENSION(:),ALLOCATABLE       :: ap_ws_gr  ! 
   !REAL,DIMENSION(:),ALLOCATABLE       :: sv_part

   ! ----------------------------------------------------------------------------
   ! PARAMETRIZATION of SUBSTANCE OUTFLOWS
   ! ----------------------------------------------------------------------------

   !LOGICAL              ,DIMENSION(:)      ,ALLOCATABLE :: l_rejflux
   !INTEGER                                              :: nbrej,nbobc_cv
   !INTEGER              ,DIMENSION(:)      ,ALLOCATABLE :: krej,numriv
   !INTEGER              ,DIMENSION(:)      ,ALLOCATABLE :: nucon
   !REAL       ,DIMENSION(:)      ,ALLOCATABLE :: xrej,yrej
   !REAL       ,DIMENSION(:,:)    ,ALLOCATABLE :: valrej

   ! ----------------------------------------------------------------------------
   ! fixed substances in water column
   ! ----------------------------------------------------------------------------
      REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE                  :: cini_wat_fix
      LOGICAL,DIMENSION(:),ALLOCATABLE               :: l_out_subs_fix
      CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE :: init_cv_name_fix
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE	         :: cv_watfix

#ifdef key_benthic
   ! ---------------------------------------------------------------------------
   ! VARIABLES BENTHIQUES
   ! ---------------------------------------------------------------------------
      CHARACTER(LEN=lchain):: filespcbenthic 
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE      :: cv_bent
   !REAL,DIMENSION(:,:,:),ALLOCATABLE      :: c_bent2D
   !REAL,DIMENSION(:,:,:,:),ALLOCATABLE    :: c_bent3D
      CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE :: name_var_bent,long_name_var_bent, &
                                                      standard_name_var_bent,unit_var_bent
      REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE :: cbentmax,hbentmin,hbentmax,cini_bent
      LOGICAL,DIMENSION(:),ALLOCATABLE :: l_bent_sedvert, l_bent_drive,l_bent_out,l_out_subs_bent
#endif

     ! ----------------------------------------------------------------------------
     ! atmospheric fluxes for substance
     ! ----------------------------------------------------------------------------
     ! file_flxatm_subs ! file name defini if atmospheric flux for substances
      CHARACTER(LEN=lchain):: file_flxatm_subs
      LOGICAL         :: l_subflxatm,l_subflxatm_xyt
      LOGICAL         :: l_cvrain_readfile,l_subflxatm_readfile
      REAL(KIND=rsh)            :: sflx_sub_atm_depth

   

#endif /* SUBSTANCE */
# undef  F90CODE

   !!==============================================================================

 END MODULE comsubstance
