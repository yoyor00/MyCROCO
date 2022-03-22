#include "cppdefs.h"

  MODULE comMUSTANG

#ifdef MUSTANG

   !&E==========================================================================
   !&E                   ***  MODULE  comMUSTANG  ***
   !&E
   !&E ** Purpose : declare all common variables related to sediment dynamics
   !&E==========================================================================

   !! * Modules used
   USE comsubstance ! for lchain, rsh...

   implicit none

   ! default
   public

#include "coupler_define_MUSTANG.h"
  
   !! * Shared or public variables for MUSTANG (common at all threads) but spatialized 

   REAL(kind=rlg)  ,PARAMETER :: epsi30_MUSTANG=1.e-30 , epsilon_MUSTANG=1.e-09 

   REAL(KIND=rlg)          :: tstart_dyninsed,  &   ! time beginning consolidation/diffusion/bioturbation/morphodynamic
                              tstart_morpho,    &   ! time beginning consolidation/diffusion/bioturbation/morphodynamic
                              t_dyninsed,       &   ! time of next dynamic in sediment step
                              subdt_consol,     &   ! sub time step for consolidation and particulate bioturbation  in sediment
                              dt_dyninsed,      &   ! time step for dynamic in sediment
                              t_morpho,         &   ! time of next morphodynamic  step
                              dt_morpho             ! time step for morphodynamic

   LOGICAL            :: l_fricwave,l_diffused,l_consolid,l_bioturb,l_biodiffs,l_dyn_insed
   LOGICAL            :: l_morphocoupl,l_bathy_smoothing,l_erolat_wet_cell,l_morphomesh
   LOGICAL            :: l_transfer2hydro_dhsed,l_z0hydro_coupl_init,l_z0hydro_coupl
   LOGICAL            :: l_dredging,l_MF_dhsed,l_bathy_actu,l_out_subs_diag_sed
   INTEGER            :: choice_flxdiss_diffsed,ero_option,E0_sand_option,nv_use,nlayer_surf_sed

! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   INTEGER            :: tau_cri_option,tau_cri_mud_option_eroindep
#ifdef key_MUSTANG_V2
   INTEGER            :: poro_option
#endif
   CHARACTER(len=19)  :: date_start_dyninsed,date_start_morpho
   
   REAL(KIND=rsh)                             :: aref_sand  ! parametre used in sandconcextrap

#ifdef key_MUSTANG_flocmod
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE   :: f_ws
#endif

   REAL(KIND=rsh)            :: fricwav,fws2,dzsmin,cfreshmud,csedmin,cmudcr,ros_sand_homogen, &
                                cvolmaxsort,cvolmaxmel,xperm1,xperm2,xsigma1,xsigma1sg,        &
                                xsigma2, xdifs1,xdifs2,xdifsi1,xdifsi2,                        &
                                epdifi,fexcs,cexcs,activlayer,frmudcr2,coef_frmudcr1,          &
                                x1toce_mud,x2toce_mud,E0_sand_para,n_eros_sand,E0_mud,         &
                                n_eros_mud, MF, MF_dhsed, htncrit_eros,       &
                                slopefac,csegreg,csandseg,xexp_ero,E0_sand_Cst,                &
                                coef_z0_coupl,z0_hydro_mud,z0_hydro_bed,                       &
                                xbioturbmax_part,xbioturbk_part,dbiotu0_part,dbiotum_part,     &
                                xbioturbmax_diss,xbioturbk_diss,dbiotu0_diss,dbiotum_diss,     &
                                xbioturbmax,xbioturbk,dbiotu0,dbiotum,frmud_db_max,frmud_db_min, &
                                dt_consolid,dt_diffused,dt_bioturb,subdt_bioturb
   REAL(KIND=rsh)                                :: hsed_new,coef_erolat,coef_tenfon_lat
   REAL(KIND=rsh)                                :: sed_difint,sed_difsed 
   
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE       :: alp1,alp2,alp3,alp4,alp5,typart
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE       :: diamstar,ws_sand,rosmrowsros,       &
                                                           stresscri0,tetacri0,xnielsen
   INTEGER       ,DIMENSION(:)      ,ALLOCATABLE        :: D0_funcT_opt
   REAL(KIND=rsh),DIMENSION(:)      ,ALLOCATABLE        :: D0_m0,D0_m1

! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   REAL(KIND=rsh)     :: k1HW97,k2HW97,E0_mud_para_indep,fusion_para_activlayer
#ifdef key_MUSTANG_V2
   REAL(KIND=rsh)     :: Awooster,Bwooster,Bmax_wu,     &
                                poro_min,alphabs,alphabn
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE   :: sigmapsg,stateconsol,permeab,E0_sand
#endif
                                                      
   INTEGER,DIMENSION(:,:),ALLOCATABLE              :: ksmi,ksma
   INTEGER, DIMENSION(:,:), ALLOCATABLE            :: dry_cell
#if  ! defined key_noTSdiss_insed
   INTEGER,DIMENSION(:), ALLOCATABLE               :: ivdiss
#endif

   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE   :: cv_sed

   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: c_sedtot,poro,dzs       
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: corflux,corfluy,gradvit       
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: fludif,fluconsol,fluconsol_drycell,flu_dyninsed
   
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed,z0sed,hsed0,dzsmax,hsed_previous
#if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN_MUSTANG_byHYDRO
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: dhsed_save
#endif
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: tenfon,tenfonc,tenfonw
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: raphbx,raphby,frofonx,frofony
#if defined key_tenfon_upwind
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: tenfonx,tenfony
#endif
   REAL(KIND=rlg),DIMENSION(:,:),ALLOCATABLE       :: phieau_s2w,phieau_s2w_consol,phieau_s2w_drycell
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: ustarbot,htot,alt_cw1

   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: morpho0,h0_bedrock

 
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: sal_bottom_MUSTANG,temp_bottom_MUSTANG,epn_bottom_MUSTANG
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: cw_bottom_MUSTANG,ws3_bottom_MUSTANG
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: roswat_bot

#ifdef key_sand2D
   REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE     :: rouse2D,sum_tmp  ! definition nombre de Rouse2D et SUM(dzcche*((htot-hzed)/hzed)**rouse) 
#endif

   ! ---------------------------------------------------------------------------
   ! VARIABLES FLUXES AT THE INTERFACE WATER-SEDIMENT
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_s2w,flx_w2s,flx_w2s_sum

   ! ---------------------------------------------------------------------------
   ! Initialization
   ! ---------------------------------------------------------------------------
   LOGICAL          :: l_repsed,l_unised,l_z0seduni,l_initsed_vardiss,l_dzsmaxuni,l_init_hsed
! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   LOGICAL          :: l_peph_suspension,l_dzsminuni
   LOGICAL          :: l_eroindep_noncoh,l_eroindep_mud,l_xexp_ero_cst
#ifdef key_MUSTANG_V2
   LOGICAL          :: l_peph_bedload
   LOGICAL          :: l_slope_effect_bedload,l_fsusp
   REAL(KIND=rsh)   :: coeff_dzsmin
#ifdef key_MUSTANG_debug
   LOGICAL        :: l_debug_effdep,l_debug_erosion
   REAL(kind=rlg)   :: lon_debug,lat_debug
   INTEGER          :: i_MUSTANG_debug,j_MUSTANG_debug
   CHARACTER(len=19):: date_start_debug
   REAL(KIND=rlg)   :: t_start_debug
#endif
#endif
   INTEGER          :: ksmiuni,ksmauni
   REAL(KIND=rsh)   :: cseduni,hseduni,sini_sed,tini_sed,dzsmaxuni,csed_mud_ini
   REAL(KIND=rsh)   :: z0seduni,z0sedmud,z0sedbedrock
   REAL(KIND=rsh)   :: dzsmax_bottom
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE    :: cini_sed,cv_sedini

! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   REAL(KIND=rsh)                                    :: dzsminuni,poro_mud_ini
#ifdef key_MUSTANG_V2
   LOGICAL,DIMENSION(:,:),ALLOCATABLE :: l_isitcohesive
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE       :: psi_sed
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE   :: poro_mud
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE   :: crel_mud
   REAL(KIND=rsh)                                :: hmin_bedload ! out of key_sedim_bedload because in paraMUSTANGV2.txt
#endif

   ! ---------------------------------------------------------------------------
   ! Sedim output
   ! ---------------------------------------------------------------------------
   LOGICAL                               :: l_outsed_saltemp,l_outsed_poro
   CHARACTER(len=lchain)   :: name_out_hsed,name_out_nbniv,name_out_dzs,name_out_tenfon,&
                                     name_out_tenfonc,name_out_tenfonw
   REAL(kind=riosh)             :: riog_valid_min_hsed,riog_valid_max_hsed,riog_valid_min_nbniv, &
                                   riog_valid_max_nbniv,riog_valid_min_dzs,riog_valid_max_dzs,     &
                                   riog_valid_min_tenfon,riog_valid_max_tenfon
   INTEGER                                    :: nk_nivsed_out,choice_nivsed_out
   INTEGER                                    :: nv_out3Dk_specif,nv_out3Dnv_specif,nv_out2D_specif

   REAL(KIND=rsh),DIMENSION(5)                :: ep_nivsed_out
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE    :: ep_nivsed_outp1,nivsed_out
   REAL(KIND=rsh)                             :: epmax_nivsed_out ! Max thickness from the sediment surface  
   REAL(KIND=riosh),  DIMENSION(:,:,:,:),ALLOCATABLE  :: var3D_cvsed
   REAL(KIND=riosh),  DIMENSION(:,:,:),ALLOCATABLE    :: var3D_dzs,var3D_TEMP,var3D_SAL
#if defined key_BLOOM_insed
   REAL(KIND=riosh),  DIMENSION(:,:,:,:),ALLOCATABLE  :: var3D_diagsed
   REAL(KIND=riosh),  DIMENSION(:,:,:),ALLOCATABLE  :: var2D_diagsed
#endif
#ifdef key_MUSTANG_specif_outputs
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE    :: varspecif3Dk_save,varspecif3Dnv_save
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE      :: varspecif2D_save
   REAL(KIND=riosh),DIMENSION(:,:,:),ALLOCATABLE    :: varspecif2D_out
   REAL(KIND=riosh),DIMENSION(:,:,:,:),ALLOCATABLE  :: var3D_specifout,varspecif3Dnv_out
#endif

   CHARACTER(len=lchain)   :: filrepsed,filoutsed,fileinised

   LOGICAL          :: l_outsed_flx_WS_all,l_outsed_toce,l_outsed_frmudsup
   LOGICAL          :: l_outsed_dzs_ksmax
! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   LOGICAL          :: l_outsed_flx_Bload_all,l_outsed_bil_Bload_all,l_outsed_fsusp
   LOGICAL          :: l_outsed_activlayer,l_outsed_surf,l_outsed_peph,l_outsed_eroiter
   LOGICAL          :: l_outsed_z0sed,l_outsed_flx_WS_int,l_outsed_flx_Bload_int
!
   LOGICAL          :: l_outsed_bil_Bload_int
#ifdef key_MUSTANG_V2
#ifdef key_MUSTANG_bedload
!  bedload 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE  :: flx_bx,flx_by
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE    :: slope_dhdx,slope_dhdy
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE    :: sedimask_h0plusxe
#if defined MORPHODYN_MUSTANG_byHYDRO
   INTEGER  :: it_morphoYes
#endif
#endif
#endif
   ! ---------------------------------------------------------------------------

!  used in erosion only but exchange and dimensions could depend on grid model 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_s2w_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_s2w_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_s2w_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_s2w_corjp1
#if ! defined key_nofluxwat_IWS
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: phieau_s2w_corim1
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: phieau_s2w_corip1
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: phieau_s2w_corjm1
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: phieau_s2w_corjp1
#endif
   !  used in accretion (settling) only bud exchange and dimensions could depend on grid model 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_w2s_corin
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_w2s_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_w2s_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_w2s_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_w2s_corjp1

#ifdef key_MUSTANG_flocmod
  !--------------------------
  !   explicit FLOCULATION 
  !--------------------------
  LOGICAL                :: l_ASH,l_ADS,l_COLLFRAG,l_out_MUDtot,l_out_f_dtmin,l_out_G, &
                                   l_out_f_d90,l_out_f_d10

  INTEGER               :: f_ero_iv
  REAL(KIND=rsh)         :: f_nf, f_dp0,f_alpha,f_beta,f_nb_frag,f_fter,f_ater,f_dmin_frag, &
                                   f_ero_frac,f_ero_nbfrag,f_mneg_param,f_collfragparam,f_cfcst,f_fp,f_fy
  REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE     :: f_diam,f_vol,f_rho,f_mass,f_cv,f_l3

  REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE   :: f_coll_prob_sh,f_coll_prob_ds,f_l1_sh,f_l1_ds,f_g3
  
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE   :: f_g1_sh,f_g1_ds
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE   :: f_d50,f_d90,f_d10,f_davg,f_dtmin
#endif

   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: emissivity_s
   REAL(KIND=rsh)                            :: alb,emissivity_sed
 
#if ! defined key_noTSdiss_insed
 ! Temperature in sediment 
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: phitemp_s,phitemp_sout
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE :: cp_s,mu_tempsedsurf,poro_sedsurf
   REAL(KIND=rsh)                            :: mu_tempsed1,mu_tempsed2,mu_tempsed3,  &
                                                        epsedmin_tempsed,epsedmax_tempsed,cp_suni   !!!FG(29/06/2018)
#endif

 !-------------------------------------------------------------------
 !  declaration of variables used both in MUSTANG and hydro modele
 !  and used only with module MUSTANG
 !  and which must have specific names and dimensions, different from those in MUSTANG
 !-------------------------------------------------------------------
   
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EROS_FLUX_s2w , SETTL_FLUX_w2s ,SETTL_FLUXSUM_w2s
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: CORFLUX_SAND,CORFLUY_SAND 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: WATER_FLUX_INPUTS 
   ! if temperature and salinity not ranged in the same table as substances concentrations
   !REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EROS_FLUX_s2w_TEMP,EROS_FLUX_s2w_SAL
   
  !   + fixed data used by MUSTANG but specific to MARS and therefore not necessarily known for another model
   REAL(kind=rsh)  ,PARAMETER :: valmanq=999.0
   REAL(kind=riosh),PARAMETER :: rg_valmanq_io=999.0
    ! fwet =1 if not used  
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE   :: fwet

   CONTAINS
 
#endif /* ifdef MUSTANG */

END MODULE comMUSTANG
