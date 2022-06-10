#include "cppdefs.h"

MODULE comMUSTANG

#ifdef MUSTANG

!!============================================================================
!! ***  MODULE  comMUSTANG  ***
!! Purpose : declare all common variables related to sediment dynamics
!!============================================================================

!! * Modules used
USE comsubstance ! for lchain, rsh, rlg, riosh

implicit none

! default
public

#include "coupler_define_MUSTANG.h"
  
    !! * Shared or public variables for MUSTANG 

    ! parameters
    REAL(kind=rlg), PARAMETER :: epsi30_MUSTANG = 1.e-30
    REAL(kind=rlg), PARAMETER :: epsilon_MUSTANG = 1.e-09 
    REAL(kind=rsh), PARAMETER :: valmanq = 999.0
    REAL(kind=riosh), PARAMETER :: rg_valmanq_io = 999.0

    ! namelists

    ! namsedim_init
    CHARACTER(len=19) :: date_start_dyninsed ! starting date for dynamic processes in sediment format '01/01/0000 00:00:00'
    CHARACTER(len=19) :: date_start_morpho ! starting date for morphodynamic format '01/01/0000 00:00:00'
    LOGICAL  :: l_repsed
    LOGICAL  :: l_initsed_vardiss
    LOGICAL  :: l_unised
    LOGICAL  :: l_init_hsed
    CHARACTER(len=lchain) :: filrepsed 
    CHARACTER(len=lchain) :: fileinised
    REAL(KIND=rsh) :: cseduni
    REAL(KIND=rsh) :: hseduni     
    REAL(KIND=rsh) :: csed_mud_ini     
    INTEGER        :: ksmiuni
    INTEGER        :: ksmauni  
    REAL(KIND=rsh) :: sini_sed
    REAL(KIND=rsh) :: tini_sed
    REAL(KIND=rsh) :: poro_mud_ini !only if key_MUSTANG_V2


    ! namsedim_layer
    LOGICAL  :: l_dzsmaxuni
    LOGICAL  :: l_dzsminuni !only if key_MUSTANG_V2
    REAL(KIND=rsh) :: dzsminuni !only if key_MUSTANG_V2
    REAL(KIND=rsh) :: dzsmin
    REAL(KIND=rsh) :: dzsmaxuni
    REAL(KIND=rsh) :: dzsmax_bottom
    REAL(KIND=rsh) :: k1HW97 !only if key_MUSTANG_V2
    REAL(KIND=rsh) :: k2HW97 !only if key_MUSTANG_V2
    REAL(KIND=rsh) :: fusion_para_activlayer !only if key_MUSTANG_V2
    INTEGER :: nlayer_surf_sed


    ! namsedim_bottomstress
    LOGICAL  :: l_z0seduni
    REAL(KIND=rsh) :: z0seduni
    REAL(KIND=rsh) :: z0sedmud
    REAL(KIND=rsh) :: z0sedbedrock
    LOGICAL :: l_fricwave
    LOGICAL :: l_z0hydro_coupl_init
    LOGICAL :: l_z0hydro_coupl
    REAL(KIND=rsh) :: fricwav
    REAL(KIND=rsh) :: coef_z0_coupl
    REAL(KIND=rsh) :: z0_hydro_mud
    REAL(KIND=rsh) :: z0_hydro_bed


    ! namsedim_deposition
    REAL(KIND=rsh) :: cfreshmud
    REAL(KIND=rsh) :: csedmin
    REAL(KIND=rsh) :: cmudcr    
    REAL(KIND=rsh) :: aref_sand  ! parameter used in sandconcextrap
    REAL(KIND=rsh) :: cvolmaxsort
    REAL(KIND=rsh) :: cvolmaxmel
    REAL(KIND=rsh) :: slopefac


    ! namsedim_erosion
    REAL(KIND=rsh) :: activlayer
    REAL(KIND=rsh) :: frmudcr2
    REAL(KIND=rsh) :: coef_frmudcr1
    REAL(KIND=rsh) :: x1toce_mud
    REAL(KIND=rsh) :: x2toce_mud
    REAL(KIND=rsh) :: E0_sand_para
    REAL(KIND=rsh) :: n_eros_sand
    REAL(KIND=rsh) :: E0_mud
    REAL(KIND=rsh) :: n_eros_mud
    INTEGER        :: ero_option
    INTEGER        :: E0_sand_option
    REAL(KIND=rsh) :: xexp_ero
    REAL(KIND=rsh) :: E0_sand_Cst
    REAL(KIND=rsh) :: E0_mud_para_indep !only if key_MUSTANG_V2
    LOGICAL        :: l_peph_suspension !only if key_MUSTANG_V2
    LOGICAL        :: l_eroindep_noncoh !only if key_MUSTANG_V2
    LOGICAL        :: l_eroindep_mud !only if key_MUSTANG_V2
    LOGICAL        :: l_xexp_ero_cst !only if key_MUSTANG_V2
    INTEGER        :: tau_cri_option !only if key_MUSTANG_V2
    INTEGER        :: tau_cri_mud_option_eroindep !only if key_MUSTANG_V2


#ifdef key_MUSTANG_V2
    ! namsedim_poro 
    INTEGER :: poro_option
    REAL(KIND=rsh) :: Awooster
    REAL(KIND=rsh) :: Bwooster
    REAL(KIND=rsh) :: Bmax_wu
    REAL(KIND=rsh) :: poro_min
#endif


#ifdef key_MUSTANG_V2 && key_MUSTANG_bedload
    ! namsedim_bedload 
    LOGICAL :: l_peph_bedload
    LOGICAL :: l_slope_effect_bedload
    LOGICAL :: l_fsusp
    REAL(KIND=rsh) :: alphabs
    REAL(KIND=rsh) :: alphabn
    REAL(KIND=rsh) :: hmin_bedload  
#endif


!**TODO** put under cpp key #ifdef key_MUSTANG_lateralerosion
    ! namsedim_lateral_erosion 
    REAL(KIND=rsh) :: htncrit_eros
    REAL(KIND=rsh) :: coef_erolat
    REAL(KIND=rsh) :: coef_tauskin_lat
    LOGICAL        :: l_erolat_wet_cell
!**TODO** put under cpp key #endif


    ! namsedim_consolidation 
    LOGICAL        :: l_consolid
    REAL(KIND=rsh) :: dt_consolid
    REAL(KIND=rlg) :: subdt_consol ! sub time step for consolidation and 
                                   ! particulate bioturbation  in sediment
    REAL(KIND=rsh) :: csegreg
    REAL(KIND=rsh) :: csandseg
    REAL(KIND=rsh) :: xperm1
    REAL(KIND=rsh) :: xperm2
    REAL(KIND=rsh) :: xsigma1
    REAL(KIND=rsh) :: xsigma2


    ! namsedim_diffusion     
    LOGICAL        :: l_diffused
    REAL(KIND=rsh) :: dt_diffused
    INTEGER        :: choice_flxdiss_diffsed
    REAL(KIND=rsh) :: xdifs1
    REAL(KIND=rsh) :: xdifs2
    REAL(KIND=rsh) :: xdifsi1
    REAL(KIND=rsh) :: xdifsi2
    REAL(KIND=rsh) :: epdifi
    REAL(KIND=rsh) :: fexcs


    ! namsedim_bioturb  
    LOGICAL        :: l_bioturb
    LOGICAL        :: l_biodiffs 
    REAL(KIND=rsh) :: dt_bioturb
    REAL(KIND=rsh) :: subdt_bioturb
    REAL(KIND=rsh) :: xbioturbmax_part
    REAL(KIND=rsh) :: xbioturbk_part
    REAL(KIND=rsh) :: dbiotu0_part
    REAL(KIND=rsh) :: dbiotum_part
    REAL(KIND=rsh) :: xbioturbmax_diss
    REAL(KIND=rsh) :: xbioturbk_diss
    REAL(KIND=rsh) :: dbiotu0_diss
    REAL(KIND=rsh) :: dbiotum_diss
    REAL(KIND=rsh) :: frmud_db_max
    REAL(KIND=rsh) :: frmud_db_min


    ! namsedim_morpho   
    LOGICAL :: l_morphocoupl
    LOGICAL :: l_MF_dhsed
    LOGICAL :: l_bathy_actu
    REAL(KIND=rsh) :: MF
    REAL(KIND=rlg) :: dt_morpho ! time step for morphodynamic 

#if  ! defined key_noTSdiss_insed
    ! namtempsed  
    REAL(KIND=rsh) :: mu_tempsed1
    REAL(KIND=rsh) :: mu_tempsed2
    REAL(KIND=rsh) :: mu_tempsed3
    REAL(KIND=rsh) :: epsedmin_tempsed
    REAL(KIND=rsh) :: epsedmax_tempsed
#endif


    ! namsedoutput  
    LOGICAL :: l_outsed_flx_WS_all
    LOGICAL :: l_outsed_flx_WS_int
    LOGICAL :: l_outsed_saltemp
    INTEGER :: nk_nivsed_out
    INTEGER :: choice_nivsed_out
    REAL(KIND=rsh), DIMENSION(5) :: ep_nivsed_out
    REAL(KIND=rsh) :: epmax_nivsed_out ! Max thickness from the sediment surface


#if defined key_MUSTANG_debug && defined key_MUSTANG_V2
    ! namsedim_debug 
    LOGICAL        :: l_debug_effdep
    LOGICAL        :: l_debug_erosion
    REAL(kind=rlg)   :: lon_debug
    REAL(kind=rlg)   :: lat_debug
    INTEGER          :: i_MUSTANG_debug
    INTEGER          :: j_MUSTANG_debug
    CHARACTER(len=19):: date_start_debug
#endif


#ifdef key_MUSTANG_flocmod
    ! namflocmod  
    LOGICAL :: l_ASH
    LOGICAL :: l_ADS
    LOGICAL :: l_COLLFRAG
    INTEGER :: f_ero_iv
    REAL(KIND=rsh) :: f_ater
    REAL(KIND=rsh) :: f_dmin_frag
    REAL(KIND=rsh) :: f_ero_frac
    REAL(KIND=rsh) :: f_ero_nbfrag
    REAL(KIND=rsh) :: f_mneg_param
    REAL(KIND=rsh) :: f_collfragparam
    REAL(KIND=rsh) :: f_cfcst
    REAL(KIND=rsh) :: f_fp
    REAL(KIND=rsh) :: f_fy
    REAL(KIND=rsh) :: f_dp0
    REAL(KIND=rsh) :: f_alpha
    REAL(KIND=rsh) :: f_beta
    REAL(KIND=rsh) :: f_nb_frag
    REAL(KIND=rsh) :: f_nf
#endif

! end namelist variables


    REAL(KIND=rsh) :: h0fond  ! RESIDUAL_THICKNESS_WAT

    ! fwet =1 if not used  
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: fwet

    ! Initialization 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: cini_sed
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: cv_sedini
    REAL(KIND=rsh) :: hsed_new

    ! Fluxes at the interface water-sediment
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_sum
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: EROS_FLUX_s2w
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: SETTL_FLUX_w2s
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: SETTL_FLUXSUM_w2s

    ! Sediment parameters
    REAL(KIND=rsh)            :: ros_sand_homogen 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: typart
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: diamstar
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: ws_sand
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: rosmrowsros
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: stresscri0
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: tetacri0


    INTEGER :: nv_use
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ksmi
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ksma
    REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE   :: cv_sed
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: c_sedtot
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: poro
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: dzs       
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: dzsmax
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: gradvit       

    ! Sediment height
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed0
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed_previous

    ! Bottom stress variables
    REAL(KIND=rsh) :: fws2  ! fricwav/2   
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE     :: z0sed
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin ! max bottom stress due to the combinaison current/wave (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_c ! bottom stress due to current (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_w ! bottom stress due to wave (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_x ! bottom stress - component on x axis
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_y ! bottom stress - component on y axis
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: ustarbot
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: raphbx
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: raphby

    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w
    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w_consol
    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w_drycell

    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: htot
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: alt_cw1

    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: sal_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: temp_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: epn_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: cw_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: ws3_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: roswat_bot

!**TODO** put under cppkey MUSTANG_CORFLUX
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: corflux
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: corfluy

#ifdef key_sand2D
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: rouse2D ! Rouse2D number
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: sum_tmp ! SUM(dzcche*((htot-hzed)/hzed)**rouse) 
#endif


    ! Dynamic in sediment (consolidation/diffusion/bioturbation)
    REAL(KIND=rlg)   :: tstart_dyninsed ! time beginning dynamic in sediment
    REAL(KIND=rlg)   :: t_dyninsed      ! time of next dynamic in sediment step
    REAL(KIND=rlg)   :: dt_dyninsed     ! time step for dynamic in sediment (min of dt for each process)
    LOGICAL :: l_dyn_insed ! true if (l_consolid .OR. l_bioturb .OR. l_diffused .OR. l_biodiffs)
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fludif
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fluconsol
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fluconsol_drycell
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flu_dyninsed
   
    ! Diffusion
    REAL(KIND=rsh) :: cexcs

    ! Morphodynamic
    REAL(KIND=rlg) :: tstart_morpho   ! time beginning morphodynamic
    REAL(KIND=rlg) :: t_morpho        ! time of next morphodynamic step
    REAL(KIND=rsh) :: MF_dhsed
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: morpho0
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: h0_bedrock

#ifdef key_MUSTANG_V2
    REAL(KIND=rsh) :: coeff_dzsmin
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: l_isitcohesive
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: psi_sed
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: poro_mud
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: crel_mud
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: sigmapsg
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: stateconsol
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: permeab
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: E0_sand
#ifdef  key_MUSTANG_bedload
        REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE  :: flx_bx
        REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE  :: flx_by
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: slope_dhdx
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: slope_dhdy
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: sedimask_h0plusxe
#if defined MORPHODYN_MUSTANG_byHYDRO
            INTEGER :: it_morphoYes
#endif
#endif
#ifdef key_MUSTANG_debug
        REAL(KIND=rlg)   :: t_start_debug
#endif
#endif
   

    ! Sedim output
    INTEGER :: nv_out3Dk_specif
    INTEGER :: nv_out3Dnv_specif
    INTEGER :: nv_out2D_specif
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: ep_nivsed_outp1
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: nivsed_out
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_cvsed
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_dzs
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_TEMP
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_SAL
#if defined key_BLOOM_insed
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_diagsed
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_diagsed
#endif
#ifdef key_MUSTANG_specif_outputs
    REAL(KIND=rsh), DIMENSION(:,:,:,:), ALLOCATABLE    :: varspecif3Dk_save
    REAL(KIND=rsh), DIMENSION(:,:,:,:), ALLOCATABLE    :: varspecif3Dnv_save
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: varspecif3Dnv_out
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE      :: varspecif2D_save
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE    :: varspecif2D_out
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_specifout
#endif


!**TODO** put under cpp key #if defined key_MUSTANG_lateralerosion
!  used in erosion only but exchange and dimensions could depend on grid model 
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corim1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corip1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corjm1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corjp1
!**TODO** put under cpp key #if ! defined key_nofluxwat_IWS
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corim1
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corip1
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corjm1
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corjp1
!**TODO** put under cpp key #endif
!**TODO** put under cpp key #endif

! slipdeposit : **TODO** put under cpp key key_MUSTANG_slipdeposit
   !  used in accretion (settling) only bud exchange and dimensions could depend on grid model 
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corin
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corjp1


#if ! defined key_noTSdiss_insed
    ! Temperature in sediment 
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phitemp_s
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ivdiss
    INTEGER       , DIMENSION(:), ALLOCATABLE :: D0_funcT_opt
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: D0_m0
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: D0_m1
#endif

#if ! defined key_nofluxwat_IWS && ! defined key_noTSdiss_insed
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: WATER_FLUX_INPUTS ! not operationnal, stil to code **TODO**
#endif


#ifdef key_MUSTANG_flocmod
    ! Explicit FLOCULATION 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: f_diam
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: f_ws
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: f_vol
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: f_rho
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: f_mass
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: f_l3
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: f_coll_prob_sh
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: f_coll_prob_ds
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: f_l1_sh
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: f_l1_ds
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: f_g3
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_g1_sh
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_g1_ds
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_d50
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_d90
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_d10
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_davg
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: f_dtmin
#endif


#ifdef key_BLOOM_insed
    LOGICAL :: l_out_subs_diag_sed
#endif


    CONTAINS
 
#endif /* ifdef MUSTANG */

END MODULE comMUSTANG
