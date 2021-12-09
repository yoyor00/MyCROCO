!***************************************************************************
!***************************************************************************
!Copyright or (c) or Copr. : IFREMER
!contributor(s) : IFREMER/DYNECO/DHYSED
!
!contact Ifremer : mustang@ifremer.fr
!
!This software (MUSTANG, MUd and Sand TrAnsport modelliNG) is a Fortran F90 
!computer program whose purpose is to perform sediment transport process 
!modelling coupled to hydrodynamic models.
!Full details can be obtained on https://wwz.ifremer.fr/dyneco/MUSTANG
!
!This software is governed by the CeCILL-C license under French law and
!abiding by the rules of distribution of free software. You can use, 
!modify and/ or redistribute the software under the terms of the CeCILL-C
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************

#include "cppdefs.h"
!---------------------------------------------------------------------------
!
                     MODULE comMUSTANG
!
!---------------------------------------------------------------------------
  

#ifdef MUSTANG

   !&E==========================================================================
   !&E                   ***  MODULE  comMUSTANG  ***
   !&E
   !&E
   !&E ** Purpose : declare all common variables related to sediment dynamics
   !&E              
   !&E 
   !&E ** Description :
   !&E     subroutine MUSTANG_alloc          ! allocates variables in sediment 
   !&E
   !&E ** History :
   !&E     ! 2015-12  (B.Thouvenin )    : creation from sedim.F90 -reorganization of module SEDIMARS
   !&E     ! 2018-11  (B.Thouvenin )    : reorganization for module MUSTANG
   !&E     ! 2019-06  (B.Thouvenin, P. Le Hir, B. Mengual ) :  MUSTANG V2 with bedload and new porosity evaluation
   !&E
   !&E==========================================================================


   !! * Modules used
   USE module_MUSTANG

   IMPLICIT NONE
   
!!#include "coupleur_dimhydro_MUSTANG.h"
#include "coupler_define_MUSTANG.h"

  !! * Accessibility

   ! functions & routines of this module, called outside :
   PUBLIC MUSTANG_alloc
   
   !! * Shared or public variables for MUSTANG (common at all threads) but spatialized 

   REAL(kind=rlg)  ,PARAMETER :: epsi30_MUSTANG=1.e-30 , epsdep_MUSTANG = 1.e-14, epsilon_MUSTANG=1.e-09 

   REAL(KIND=rlg),PUBLIC   :: tstart_dyninsed,  &   ! time beginning consolidation/diffusion/bioturbation/morphodynamic
                              tstart_morpho,    &   ! time beginning consolidation/diffusion/bioturbation/morphodynamic
                              t_dyninsed,       &   ! time of next dynamic in sediment step
                              subdt_consol,     &   ! sub time step for consolidation and particulate bioturbation  in sediment
                              dt_dyninsed,      &   ! time step for dynamic in sediment
                              t_morpho,         &   ! time of next morphodynamic  step
                              dt_morpho             ! time step for morphodynamic

   LOGICAL,PUBLIC            :: l_fricwave,l_diffused,l_consolid,l_bioturb,l_biodiffs,l_dyn_insed
   LOGICAL,PUBLIC            :: l_morphocoupl,l_bathy_smoothing,l_erolat_wet_cell,l_morphomesh
   LOGICAL,PUBLIC            :: l_transfer2hydro_dhsed,l_z0hydro_coupl_init,l_z0hydro_coupl
   LOGICAL,PUBLIC            :: l_dredging,l_MF_dhsed,l_bathy_actu,l_out_subs_diag_sed
   INTEGER,PUBLIC            :: choice_flxdiss_diffsed,ero_option,E0_sand_option,nv_use,nlayer_surf_sed

! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   INTEGER,PUBLIC            :: tau_cri_option,tau_cri_mud_option_eroindep
#ifdef key_MUSTANG_V2
   INTEGER,PUBLIC            :: poro_option
#endif
   CHARACTER(len=19),PUBLIC  :: date_start_dyninsed,date_start_morpho
   
#ifdef key_MARS
   ! for a hydrodynamic model other than MARS
   ! these 2 variables must be declared in comsubstance if module substance is installed
   !    (diam_sed and ros initalized in substance.F90)
   ! but in MARS, they are declared in comMUSTANG and allocate/initialized in subreaddat
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC    :: diam_sed,ros
#endif
   REAL(KIND=rsh),PUBLIC                             :: aref_sand  ! parametre used in sandconcextrap

#ifdef key_MUSTANG_flocmod
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC   :: f_ws
#endif

   REAL(KIND=rsh),PUBLIC     :: fricwav,fws2,dzsmin,cfreshmud,csedmin,cmudcr,ros_sand_homogen, &
                                cvolmaxsort,cvolmaxmel,xperm1,xperm2,xsigma1,xsigma1sg,        &
                                xsigma2, xdifs1,xdifs2,xdifsi1,xdifsi2,                        &
                                epdifi,fexcs,cexcs,activlayer,frmudcr2,coef_frmudcr1,          &
                                x1toce_mud,x2toce_mud,E0_sand_para,n_eros_sand,E0_mud,         &
                                n_eros_mud,corfluer1,corfluer2,MF,MF_dhsed,htncrit_eros,       &
                                slopefac,csegreg,csandseg,xexp_ero,E0_sand_Cst,                &
                                coef_z0_coupl,z0_hydro_mud,z0_hydro_bed,                       &
                                xbioturbmax_part,xbioturbk_part,dbiotu0_part,dbiotum_part,     &
                                xbioturbmax_diss,xbioturbk_diss,dbiotu0_diss,dbiotum_diss,     &
                                xbioturbmax,xbioturbk,dbiotu0,dbiotum,frmud_db_max,frmud_db_min, &
                                dt_consolid,dt_diffused,dt_bioturb,subdt_bioturb
   REAL(KIND=rsh),PUBLIC                                :: hsed_new,coef_erolat,coef_tenfon_lat
   REAL(KIND=rsh),PUBLIC                                :: sed_difint,sed_difsed 
   
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC       :: alp1,alp2,alp3,alp4,alp5,typart
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC       :: diamstar,ws_sand,rosmrowsros,       &
                                                           stresscri0,tetacri0,xnielsen
   INTEGER       ,DIMENSION(:)      ,ALLOCATABLE        :: D0_funcT_opt
   REAL(KIND=rsh),DIMENSION(:)      ,ALLOCATABLE        :: D0_m0,D0_m1

! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   REAL(KIND=rsh),PUBLIC     :: k1HW97,k2HW97,E0_mud_para_indep,fusion_para_activlayer
#ifdef key_MUSTANG_V2
   REAL(KIND=rsh),PUBLIC     :: Awooster,Bwooster,Bmax_wu,     &
                                poro_min,alphabs,alphabn
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC   :: sigmapsg,stateconsol,permeab,E0_sand
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
#if defined key_MARS
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: morphox,morphoy,hxi,hyi,h0_bedrock
#else
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: morpho0,h0_bedrock
#endif
 
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: sal_bottom_MUSTANG,temp_bottom_MUSTANG,epn_bottom_MUSTANG
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: cw_bottom_MUSTANG,ws3_bottom_MUSTANG
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: roswat_bot

#ifdef key_sand2D
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE    :: rouse2D,sum_tmp  ! definition nombre de Rouse2D et SUM(dzcche*((htot-hzed)/hzed)**rouse) 
#endif

!#ifdef key_dredging
!   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE        :: dredgmassloc
!   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE        :: dredgmassloc_cum
!   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE          :: dredgmass
!   INTEGER,DIMENSION(imin:imax,jmin:jmax)           :: dredgflag
!   REAL(KIND=rsh),DIMENSION(9)                      :: dredglev
!   CHARACTER(len=lchain),PUBLIC   :: filsauvdredg,filrepdredg,filoutdredg
!#endif
#ifdef key_MARS
#if defined key_castest_2DVSN
   REAL(KIND=rsh),PUBLIC     :: tenfoncshtx_jmin
#ifdef key_wave_crossshore
   REAL(KIND=rsh),PUBLIC     :: hwaveoff,phas_fracan_wave,break_factor_cas2DV,fact_fws2_cas2DV
   CHARACTER(LEN=lchain),PUBLIC  :: name_wave_hs,name_out_wave_hs
   REAL(KIND=riosh), PUBLIC      :: riog_valid_min_wave_hs,riog_valid_max_wave_hs
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC      :: hwave,uwave,htw,sqrthtw
   REAL(KIND=rsh), PUBLIC, DIMENSION(:,:), ALLOCATABLE ::  wave_hs
#endif
#endif
#endif

   ! ---------------------------------------------------------------------------
   ! VARIABLES FLUXES AT THE INTERFACE WATER-SEDIMENT
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: flx_s2w,flx_w2s,flx_w2s_sum

   ! ---------------------------------------------------------------------------
    ! relatif to initialization
   ! ---------------------------------------------------------------------------
   LOGICAL,PUBLIC          :: l_repsed,l_unised,l_z0seduni,l_initsed_vardiss,l_dzsmaxuni,l_init_hsed
! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   LOGICAL,PUBLIC          :: l_peph_suspension,l_dzsminuni
   LOGICAL,PUBLIC          :: l_eroindep_noncoh,l_eroindep_mud,l_xexp_ero_cst
#ifdef key_MUSTANG_V2
   LOGICAL,PUBLIC          :: l_peph_bedload
   LOGICAL,PUBLIC          :: l_slope_effect_bedload,l_fsusp
   REAL(KIND=rsh),PUBLIC   :: coeff_dzsmin
#ifdef key_MUSTANG_debug
   LOGICAL, PUBLIC         :: l_debug_effdep,l_debug_erosion
   REAL(kind=rlg),PUBLIC   :: lon_debug,lat_debug
   INTEGER,PUBLIC          :: i_MUSTANG_debug,j_MUSTANG_debug
   CHARACTER(len=19),PUBLIC:: date_start_debug
   REAL(KIND=rlg),PUBLIC   :: t_start_debug
#endif
#endif
   INTEGER,PUBLIC          :: ksmiuni,ksmauni
   REAL(KIND=rsh),PUBLIC   :: cseduni,hseduni,sini_sed,tini_sed,dzsmaxuni,csed_mud_ini
   REAL(KIND=rsh),PUBLIC   :: z0seduni,z0sedmud,z0sedbedrock
   REAL(KIND=rsh),PUBLIC   :: dzsmax_bottom
   REAL(KIND=rsh),PUBLIC,DIMENSION(:),ALLOCATABLE    :: cini_sed,cv_sedini

! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   REAL(KIND=rsh),PUBLIC                                    :: dzsminuni,poro_mud_ini
#ifdef key_MUSTANG_V2
   LOGICAL,DIMENSION(:,:),ALLOCATABLE,PUBLIC :: l_isitcohesive
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC           :: psi_sed
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,    PUBLIC   :: poro_mud
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,    PUBLIC   :: crel_mud
   REAL(KIND=rsh),PUBLIC                                    :: hmin_bedload ! out of key_sedim_bedload because in paraMUSTANGV2.txt
#endif


   ! ---------------------------------------------------------------------------
   ! relatif to sedim output
   ! ---------------------------------------------------------------------------
   LOGICAL,PUBLIC                               :: l_outsed_saltemp,l_outsed_poro
   CHARACTER(len=lchain),PUBLIC   :: name_out_hsed,name_out_nbniv,name_out_dzs,name_out_tenfon,&
                                     name_out_tenfonc,name_out_tenfonw
   REAL(kind=riosh),PUBLIC      :: riog_valid_min_hsed,riog_valid_max_hsed,riog_valid_min_nbniv, &
                                   riog_valid_max_nbniv,riog_valid_min_dzs,riog_valid_max_dzs,     &
                                   riog_valid_min_tenfon,riog_valid_max_tenfon
   INTEGER,PUBLIC                                    :: nk_nivsed_out,choice_nivsed_out
   INTEGER,PUBLIC                                    :: nv_out3Dk_specif,nv_out3Dnv_specif,nv_out2D_specif

   REAL(KIND=rsh),DIMENSION(5),PUBLIC                :: ep_nivsed_out
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC    :: ep_nivsed_outp1,nivsed_out
   REAL(KIND=rsh),PUBLIC                             :: epmax_nivsed_out ! Max thickness from the sediment surface  
   REAL(KIND=riosh),  DIMENSION(:,:,:,:),ALLOCATABLE  :: var3D_cvsed
   REAL(KIND=riosh),  DIMENSION(:,:,:),ALLOCATABLE    :: var3D_dzs,var3D_TEMP,var3D_SAL
#if defined key_BLOOM_insed
   REAL(KIND=riosh),  DIMENSION(:,:,:,:),ALLOCATABLE  :: var3D_diagsed
   REAL(KIND=riosh),  DIMENSION(:,:,:),ALLOCATABLE  :: var2D_diagsed
#endif
#ifdef key_MUSTANG_specif_outputs
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,    PUBLIC   :: varspecif3Dk_save,varspecif3Dnv_save
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,      PUBLIC   :: varspecif2D_save
   REAL(KIND=riosh),DIMENSION(:,:,:),ALLOCATABLE,      PUBLIC   :: varspecif2D_out
   REAL(KIND=riosh),DIMENSION(:,:,:,:),ALLOCATABLE,    PUBLIC :: var3D_specifout,varspecif3Dnv_out
#endif

   CHARACTER(len=lchain),PUBLIC   :: filrepsed,filoutsed,fileinised
#ifdef key_sand2D
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC                :: l_outsandrouse
#endif

   LOGICAL,PUBLIC          :: l_outsed_flx_WS_all,l_outsed_toce,l_outsed_frmudsup
   LOGICAL,PUBLIC          :: l_outsed_dzs_ksmax
! used only if key_MUSTANG_V2 but need to be declared  for MUSTANG input file
   LOGICAL,PUBLIC          :: l_outsed_flx_Bload_all,l_outsed_bil_Bload_all,l_outsed_fsusp
   LOGICAL,PUBLIC          :: l_outsed_activlayer,l_outsed_surf,l_outsed_peph,l_outsed_eroiter
   LOGICAL,PUBLIC          :: l_outsed_z0sed,l_outsed_flx_WS_int,l_outsed_flx_Bload_int
!
   LOGICAL,PUBLIC          :: l_outsed_bil_Bload_int
#ifdef key_MUSTANG_V2
#ifdef key_MUSTANG_bedload
!  bedload (charriage)
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC  :: flx_bx,flx_by
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC    :: slope_dhdx,slope_dhdy
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC    :: sedimask_h0plusxe
#if defined MORPHODYN_MUSTANG_byHYDRO
   INTEGER  :: it_morphoYes
#endif
#endif
#endif
   ! ---------------------------------------------------------------------------

!  used in erosion only but exchange and dimensions could depend on grid model 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_s2w_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_s2w_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_s2w_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_s2w_corjp1
#if ! defined key_nofluxwat_IWS
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: phieau_s2w_corim1
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: phieau_s2w_corip1
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: phieau_s2w_corjm1
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: phieau_s2w_corjp1
#endif
   !  used in accretion (settling) only bud exchange and dimensions could depend on grid model 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_w2s_corin
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_w2s_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_w2s_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_w2s_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: flx_w2s_corjp1

#ifdef key_MUSTANG_flocmod
  !--------------------------
  !   explicit FLOCULATION 
  !--------------------------
  LOGICAL,PUBLIC                :: l_ASH,l_ADS,l_COLLFRAG,l_out_MUDtot,l_out_f_dtmin,l_out_G, &
                                   l_out_f_d90,l_out_f_d10

  INTEGER, PUBLIC               :: f_ero_iv
  REAL(KIND=rsh),PUBLIC         :: f_nf, f_dp0,f_alpha,f_beta,f_nb_frag,f_fter,f_ater,f_dmin_frag, &
                                   f_ero_frac,f_ero_nbfrag,f_mneg_param,f_collfragparam,f_cfcst,f_fp,f_fy
  REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC     :: f_diam,f_vol,f_rho,f_mass,f_cv,f_l3

  REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC   :: f_coll_prob_sh,f_coll_prob_ds,f_l1_sh,f_l1_ds,f_g3
  
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC   :: f_g1_sh,f_g1_ds
  REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC   :: f_d50,f_d90,f_d10,f_davg,f_dtmin
#endif

   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: emissivity_s
   REAL(KIND=rsh),PUBLIC                            :: alb,emissivity_sed
 
#if ! defined key_noTSdiss_insed
 ! relatif to Temperature in sediment 
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: phitemp_s,phitemp_sout
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: cp_s,mu_tempsedsurf,poro_sedsurf
   REAL(KIND=rsh),PUBLIC                            :: mu_tempsed1,mu_tempsed2,mu_tempsed3,  &
                                                        epsedmin_tempsed,epsedmax_tempsed,cp_suni   !!!FG(29/06/2018)
#endif


#if ! defined key_MARS
 !-------------------------------------------------------------------
 !  declaration of variables used both in MUSTANG and hydro modele
 !  and used only with module MUSTANG
 !  and which must have specific names and dimensions, different from those in MUSTANG
 !  not used in MARS because same variable (name and dimensions) in MUSTANG and  in hydro modele
 !  (ksdmin and ksdmax defined in parameters in MARS)
 !-------------------------------------------------------------------
   ! a eliminer pour CROCO
   ! INTEGER, PARAMETER,PUBLIC :: ksdmin=1,ksdmax=100
   
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: EROS_FLUX_s2w , SETTL_FLUX_w2s ,SETTL_FLUXSUM_w2s
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: CORFLUX_SAND,CORFLUY_SAND 
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: WATER_FLUX_INPUTS 
   ! if temperature and salinity not ranged in the same table as substances concentrations
   !REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EROS_FLUX_s2w_TEMP,EROS_FLUX_s2w_SAL
   
  !   + fixed data used by MUSTANG but specific to MARS and therefore not necessarily known for another model
   LOGICAL :: l_testcase=.FALSE.
   REAL(kind=rsh)  ,PARAMETER :: valmanq=999.0
   REAL(kind=riosh),PARAMETER :: rg_valmanq_io=999.0
    ! fwet et fwetp =1 if not used  
   REAL(KIND=rsh),PUBLIC,DIMENSION(:,:),ALLOCATABLE   :: fwet,fwetp

#endif


 CONTAINS
 
   !!===========================================================================
 
  SUBROUTINE MUSTANG_alloc(l_filesubs)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_alloc  ***
   !&E
   !&E ** Purpose : allocation of arrays relative to sediment
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : MUSTANG_initialization 
   !&E                (and in MARS by cas_init_subs_disspart for test cases
   !&E                        and  by Agrif_User for Agrif runs)
   !&E
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#include "coupler_define_MUSTANG.h"

#ifdef key_MARS
   USE parameters,   ONLY : l_testcase
   USE comvars2d,    ONLY : iscreenlog
   USE comsubstance, ONLY : irk_fil,diam_r,nvpc,   &
                            cini_sed_r, irkm_var_assoc,ros_r,unit_modif_mudbio_N2dw
#ifdef key_MUSTANG_flocmod
   USE comsubstance, ONLY : nv_mud
   USE parameters,   ONLY : num_testcase

#endif

#else
   ! if module substance is installed in hydro modele
   !USE comsubstance, ONLY : nv_tot,nv_adv,nv_state,nvp,nv_mud,           &
   !                              irk_fil,irkm_var_assoc,unit_modif_mudbio_N2dw
   USE comsubstance
#endif

   !! * Arguments
   LOGICAL, INTENT(IN), OPTIONAL   :: l_filesubs

   !! * Local declarations
   INTEGER               :: iv
   CHARACTER(len=lchain) :: filepc
   LOGICAL               :: l_filesubsr

   !!--------------------------------------------------------------------------
   !! * Executable part
   IF(PRESENT(l_filesubs)) THEN
     l_filesubsr=l_filesubs
   ELSE
     l_filesubsr=.FALSE.
   ENDIF

!  allocation of 1DV variables 
   !ALLOCATE(ros(nvp))
   !ALLOCATE(diam_sed(nvp))
   
   !ros(1:nvp)=0.0_rsh
   !diam_sed(1:nvp)=0.0_rsh
   
   ALLOCATE(ws_sand(nvp))        
   ALLOCATE(diamstar(nvp))
   ALLOCATE(rosmrowsros(nvp))
   ALLOCATE(stresscri0(nvp))
   ALLOCATE(tetacri0(nvp))
   ALLOCATE(xnielsen(nvp))
#if ! defined key_MARS
   ALLOCATE (typart(-1:nv_adv))
   typart(-1:0)=0.0_rsh
   typart(1:nvpc)=1.0_rsh
   typart(nvpc+1:nv_adv)=0.0_rsh
#endif
#ifdef key_MUSTANG_V2
   ALLOCATE(E0_sand(nvp))
   E0_sand(1:nvp)=0.0_rsh
#endif   
#if  ! defined key_noTSdiss_insed
   ALLOCATE(ivdiss(-1:nv_adv-nvp))
#endif

      
#ifdef key_MUSTANG_flocmod
   ALLOCATE(f_ws(1:nv_mud))
   f_ws(1:nv_mud)=0.0_rsh
#endif   

#ifdef key_MUSTANG_V2
   ALLOCATE(psi_sed(nvp))
   psi_sed(1:nvp)=0.0_rsh
#endif   

!  allocation of spatial variables  
!  dimensions defined dans coupler_define_MUSTANG.h
!  dimensions in MARS :  PROC_IN_ARRAY       = limin:limax,ljmin:ljmax
!                        PROC_IN_ARRAY_m1p2  = liminm1:limaxp2,ljminm1:ljmaxp2
!                        PROC_IN_ARRAY_m1p1  = liminm1:limaxp1,ljminm1:ljmaxp1
!                        PROC_IN_ARRAY_0p1   = limin:limaxp1,ljmin:ljmaxp1

   ALLOCATE(ksmi(PROC_IN_ARRAY))
   ALLOCATE(ksma(PROC_IN_ARRAY))
   ALLOCATE(hsed(PROC_IN_ARRAY))  
   ALLOCATE(z0sed(PROC_IN_ARRAY))
   ALLOCATE(tenfon(PROC_IN_ARRAY))
   ALLOCATE(tenfonc(PROC_IN_ARRAY))
   ALLOCATE(tenfonw(PROC_IN_ARRAY))
   ALLOCATE(ustarbot(PROC_IN_ARRAY))
   ALLOCATE(dzsmax(PROC_IN_ARRAY))
#if ! defined key_noTSdiss_insed
   ALLOCATE(phitemp_s(PROC_IN_ARRAY))
   ALLOCATE(phitemp_sout(PROC_IN_ARRAY))
   ALLOCATE(cp_s(PROC_IN_ARRAY))
   !ALLOCATE(mu_tempsedsurf(PROC_IN_ARRAY))
   ALLOCATE(poro_sedsurf(PROC_IN_ARRAY))
#endif
   ALLOCATE(emissivity_s(PROC_IN_ARRAY))
   ALLOCATE(htot(PROC_IN_ARRAY_m2p2))
   ALLOCATE(alt_cw1(PROC_IN_ARRAY))
   ALLOCATE(epn_bottom_MUSTANG(PROC_IN_ARRAY_m1p2))  
   ALLOCATE(sal_bottom_MUSTANG(PROC_IN_ARRAY_m1p2))  
   ALLOCATE(temp_bottom_MUSTANG(PROC_IN_ARRAY_m1p2))
   ALLOCATE(cw_bottom_MUSTANG(nv_tot,PROC_IN_ARRAY_m1p2))
   ALLOCATE(ws3_bottom_MUSTANG(nvp,PROC_IN_ARRAY_m1p2))  
   ALLOCATE(roswat_bot(PROC_IN_ARRAY))  
#ifdef key_MUSTANG_V2
   ALLOCATE(sigmapsg(ksdmin:ksdmax))
   ALLOCATE(stateconsol(ksdmin:ksdmax))
   ALLOCATE(permeab(ksdmin:ksdmax))
#endif
#ifdef key_sand2D
   ALLOCATE(rouse2D(nv_adv,PROC_IN_ARRAY))
   ALLOCATE(sum_tmp(nv_adv,PROC_IN_ARRAY))
   rouse2D(1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   sum_tmp(1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
#endif

   ksmi(PROC_IN_ARRAY)=0
   ksma(PROC_IN_ARRAY)=0
   hsed(PROC_IN_ARRAY)=0.0_rsh
   z0sed(PROC_IN_ARRAY)=0.0_rsh
   tenfon(PROC_IN_ARRAY)=0.0_rsh
   tenfonc(PROC_IN_ARRAY)=0.0_rsh
   tenfonw(PROC_IN_ARRAY)=0.0_rsh
   ustarbot(PROC_IN_ARRAY)=0.0_rsh
   
#if defined key_MARS && defined key_castest_2DVSN
    tenfoncshtx_jmin=0.0_rsh
#ifdef key_wave_crossshore
    ALLOCATE(hwave(jmin:jmax))
     hwave(:)=0.0_rsh
    ALLOCATE(uwave(jmin:jmax))
    ALLOCATE(wave_hs(PROC_IN_ARRAY))
    riog_valid_min_wave_hs=-0.1
    riog_valid_max_wave_hs=20.0
    name_out_wave_hs="hs"
#endif
#endif

   dzsmax(PROC_IN_ARRAY)=0.0_rsh
   htot(PROC_IN_ARRAY_m2p2)=0.0_rsh
   alt_cw1(PROC_IN_ARRAY)=0.0_rsh
   emissivity_s(PROC_IN_ARRAY)=0.0_rsh
#if ! defined key_noTSdiss_insed
   phitemp_s(PROC_IN_ARRAY)=0.0_rsh
   phitemp_sout(PROC_IN_ARRAY)=0.0_rsh
   cp_s(PROC_IN_ARRAY)=0.0_rsh
   !mu_tempsedsurf(PROC_IN_ARRAY)=0.0_rsh
   poro_sedsurf(PROC_IN_ARRAY)=0.0_rsh
#endif
   epn_bottom_MUSTANG(PROC_IN_ARRAY_m1p2)=0.0_rsh
   sal_bottom_MUSTANG(PROC_IN_ARRAY_m1p2)=0.0_rsh
   temp_bottom_MUSTANG(PROC_IN_ARRAY_m1p2)=0.0_rsh
   cw_bottom_MUSTANG(nv_tot,PROC_IN_ARRAY_m1p2)=0.0_rsh
   ws3_bottom_MUSTANG(nvp,PROC_IN_ARRAY_m1p2)=0.0_rsh
   roswat_bot(PROC_IN_ARRAY)=0.0_rsh

   ALLOCATE(cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY))
   ALLOCATE(c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY))
   ALLOCATE(poro(ksdmin:ksdmax,PROC_IN_ARRAY))
   ALLOCATE(dzs(ksdmin:ksdmax,PROC_IN_ARRAY))
   ALLOCATE(flx_s2w(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(flx_w2s(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(flx_w2s_sum(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(corflux(nv_adv,PROC_IN_ARRAY_m1p1  ))
   ALLOCATE(corfluy(nv_adv,PROC_IN_ARRAY_m1p1))
   ALLOCATE(fludif(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(fluconsol(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(fluconsol_drycell(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(flu_dyninsed(-1:nv_adv,PROC_IN_ARRAY))
   ALLOCATE(gradvit(NB_LAYER_WAT,PROC_IN_ARRAY))

   cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
   c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
   poro(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
   dzs(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
   flx_s2w(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   flx_w2s(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   flx_w2s_sum(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   corflux(1:nv_adv,PROC_IN_ARRAY_m1p1  )=1.0_rsh
   corfluy(1:nv_adv,PROC_IN_ARRAY_m1p1)=1.0_rsh
   fludif(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   fluconsol(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   fluconsol_drycell(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   flu_dyninsed(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
   gradvit(NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh

#ifdef key_MUSTANG_specif_outputs
  ! outputs
  ! variables 3D /k
   nv_out3Dk_specif=1
      ! 1 : poro_save  
#if defined key_MUSTANG_add_consol_outputs && defined key_MUSTANG_V2
   nv_out3Dk_specif=nv_out3Dk_specif+8
      ! 2 : loadograv_save
      ! 3 : permeab_save
      ! 4 : sigmapsg_save
      ! 5 : dtsdzs_save
      ! 6 : hinder_save
      ! 7 : sed_rate_save
      ! 8 : sigmadjge_save
      ! 9 : stateconsol_save
#endif
   ALLOCATE(varspecif3Dk_save(nv_out3Dk_specif,ksdmin:ksdmax,PROC_IN_ARRAY))
   varspecif3Dk_save(1:nv_out3Dk_specif,ksdmin:ksdmax,PROC_IN_ARRAY)= 0.0_rsh

    nv_out3Dnv_specif=3
      ! 1 : toce_save
      ! 2 : flx_s2w_save
      ! 3 : flx_w2s_save
#ifdef key_MUSTANG_V2
    nv_out3Dnv_specif=nv_out3Dnv_specif+1
      ! 4 : pephm_fcor_save  
#ifdef key_MUSTANG_bedload
    nv_out3Dnv_specif=nv_out3Dnv_specif+4
      ! 5 : flx_bx
      ! 6 : flx_by
      ! 7 : bil_bedload
      ! 8 : fsusp
#endif
#endif
   ALLOCATE(varspecif3Dnv_save(nv_out3Dnv_specif,nvpc,PROC_IN_ARRAY))
   ALLOCATE(varspecif3Dnv_out(nv_out3Dnv_specif,nvpc,PROC_IN_ARRAY))
   varspecif3Dnv_save(1:nv_out3Dnv_specif,1:nvpc,PROC_IN_ARRAY)= 0.0_rsh
   
    nv_out2D_specif=2
      ! 1 : frmudsup 
      ! 2 : dzs_ksmax 
#ifdef key_MUSTANG_V2
       nv_out2D_specif=nv_out2D_specif+13
      ! 3 : dzs_aclay_comp_save
      ! 4 : dzs_aclay_kept_save
      ! 5 : tero_noncoh (cumulated time (in hours) elapsed in non cohesive regime)
      ! 6 : tero_coh (cumulated time (in hours) elapsed in cohesive regime)
      ! 7 : pct_iter_noncoh
      ! 8 : pct_iter_coh
      ! 9 : niter_ero
      ! 10: z0sed
      ! 11 : flx_s2w_noncoh
      ! 12 : flx_w2s_noncoh
      ! 13 : flx_s2w_coh
      ! 14 : flx_w2s_coh
      ! 15: z0hydro (if l_z0hydro_coupl)
#ifdef key_MUSTANG_bedload
    nv_out2D_specif=nv_out2D_specif+3
      ! 16 : flx_bx_int
      ! 17 : flx_by_int
      ! 18 : bil_bedload_int
#endif
!   end version MUSTANG V2
#endif

    ALLOCATE(varspecif2D_save(nv_out2D_specif,PROC_IN_ARRAY))
    ALLOCATE(varspecif2D_out(nv_out2D_specif,PROC_IN_ARRAY))
    varspecif2D_save(1:nv_out2D_specif,PROC_IN_ARRAY)= 0.0_rsh
    
!   end specifs output
#endif

#ifdef key_MUSTANG_V2
   ALLOCATE(poro_mud(ksdmin:ksdmax,PROC_IN_ARRAY))
   ALLOCATE(crel_mud(ksdmin:ksdmax,PROC_IN_ARRAY))
   ALLOCATE(l_isitcohesive(PROC_IN_ARRAY))
   poro_mud(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
   crel_mud(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
   l_isitcohesive(PROC_IN_ARRAY)=.FALSE.
 
#ifdef key_MUSTANG_bedload
   ALLOCATE( flx_bx(1:nvp,PROC_IN_ARRAY_m1p1)) ! attention /Baptiste : m1p1 sur les 2 indices au lieu d1
   ALLOCATE( flx_by(1:nvp,PROC_IN_ARRAY_m1p1)) ! attention /Baptiste : m1p1 sur les 2 indices au lieu d1
   ALLOCATE( slope_dhdx(PROC_IN_ARRAY))
   ALLOCATE( slope_dhdy(PROC_IN_ARRAY))
   ALLOCATE( sedimask_h0plusxe(PROC_IN_ARRAY_m1p1)) ! attention /Baptiste : m1p1 sur les 2 indices au lieu d1

   flx_bx(1:nvp,PROC_IN_ARRAY_m1p1)=0.0_rsh 
   flx_by(1:nvp,PROC_IN_ARRAY_m1p1)=0.0_rsh
   slope_dhdx(PROC_IN_ARRAY)=0.0_rsh
   slope_dhdy(PROC_IN_ARRAY)=0.0_rsh
   sedimask_h0plusxe(PROC_IN_ARRAY_m1p1)=0.0_rsh
#endif
#endif

   ALLOCATE(phieau_s2w(PROC_IN_ARRAY))
   ALLOCATE(phieau_s2w_drycell(PROC_IN_ARRAY))
   ALLOCATE(phieau_s2w_consol(PROC_IN_ARRAY))
   phieau_s2w(PROC_IN_ARRAY)=0.0_rlg
   phieau_s2w_drycell(PROC_IN_ARRAY)=0.0_rlg
   phieau_s2w_consol(PROC_IN_ARRAY)=0.0_rlg

   ALLOCATE( raphbx(PROC_IN_ARRAY_m1p1),raphby(PROC_IN_ARRAY_m1p1) )
   ALLOCATE( frofonx(PROC_IN_ARRAY_m1p1),frofony(PROC_IN_ARRAY_m1p1) )
#ifdef key_tenfon_upwind
   ALLOCATE( tenfonx(PROC_IN_ARRAY),tenfony(PROC_IN_ARRAY) )
#endif
   ALLOCATE( dry_cell(PROC_IN_ARRAY))
   raphbx(PROC_IN_ARRAY_m1p1)=0.0_rsh
   raphby(PROC_IN_ARRAY_m1p1)=0.0_rsh
   dry_cell(PROC_IN_ARRAY)=0

   ALLOCATE(flx_s2w_corim1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_s2w_corip1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_s2w_corjm1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_s2w_corjp1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_w2s_corin(nvp,PROC_IN_ARRAY))
   ALLOCATE(flx_w2s_corim1(nvp,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_w2s_corip1(nvp,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_w2s_corjm1(nvp,PROC_IN_ARRAY_m1p1))
   ALLOCATE(flx_w2s_corjp1(nvp,PROC_IN_ARRAY_m1p1))
#if ! defined key_nofluxwat_IWS
   ALLOCATE(phieau_s2w_corim1(PROC_IN_ARRAY_m1p1))
   ALLOCATE(phieau_s2w_corip1(PROC_IN_ARRAY_m1p1))
   ALLOCATE(phieau_s2w_corjm1(PROC_IN_ARRAY_m1p1))
   ALLOCATE(phieau_s2w_corjp1(PROC_IN_ARRAY_m1p1))
#endif
   
! relatif to initialization
   ALLOCATE(cini_sed(nv_state))
   
#ifdef key_MARS
   IF ((.NOT.l_testcase) .OR. (l_testcase .AND. l_filesubsr)) THEN
      DO iv=1,nvp
        ros(iv)=ros_r(irk_fil(iv))
        diam_sed(iv)=diam_r(irk_fil(iv))
      ENDDO
      ! def of inital concentrations in sediment and for particulate variable MUDB 
      ! conversion of initial value in sediment expressed in mmole/kg (as NoCP variable)
      ! in kg/kg of sediment
      DO iv=1,nv_state
        cini_sed(iv)=cini_sed_r(irk_fil(iv))*unit_modif_mudbio_N2dw(irk_fil(iv))
      ENDDO
   ENDIF
#else
      DO iv=1,nv_adv
        cini_sed(iv)=cini_sed_r(irk_fil(iv))*unit_modif_mudbio_N2dw(irk_fil(iv))
      ENDDO
#endif    

#ifdef key_Pconstitonly_insed 
      nv_use=nvpc
#else
      nv_use=nvp
#endif  

#if ! defined key_noTSdiss_insed
! counting of dissolved variables for diffusion in the sediment
   ivdiss(:)=0
   ivdiss(-1)=-1
   ivdiss(0)=0
#if ! defined key_Pconstitonly_insed
   DO iv=1,nv_adv-nvp
      ivdiss(iv)=iv+nvp
   ENDDO  
#endif
#endif

!  option morpho
!!!!!!!!!!!!!!!!!!!!
   IF(l_morphocoupl) THEN
!! ATTENTION : no hx, hy in other model than MARS 
#ifdef key_MARS
       ALLOCATE(hsed0(PROC_IN_ARRAY))
       ALLOCATE(hsed_previous(PROC_IN_ARRAY))
#if defined key_oasis && defined key_oasis_mars_ww3  
       ALLOCATE(dhsed_save(PROC_IN_ARRAY))  
#endif
       ALLOCATE(morphox(ARRAY_morphox))
       ALLOCATE(morphoy(ARRAY_morphoy))
       ALLOCATE(hxi(ARRAY_hxi))
       ALLOCATE(hyi(ARRAY_hyi))
       ALLOCATE(h0_bedrock(ARRAY_h0_bedrock))
       hxi(ARRAY_hxi)=0.0_rsh
       hyi(ARRAY_hyi)=0.0_rsh
       hsed0(PROC_IN_ARRAY)=0.0_rsh
       hsed_previous(PROC_IN_ARRAY)=0.0_rsh
#if defined key_oasis && defined key_oasis_mars_ww3  
       dhsed_save(PROC_IN_ARRAY)=0.0_rsh
#endif
#else
       ALLOCATE(morpho0(ARRAY_morpho))
       ALLOCATE(h0_bedrock(ARRAY_h0_bedrock))
       ALLOCATE(hsed0(PROC_IN_ARRAY))
       ALLOCATE(hsed_previous(PROC_IN_ARRAY))
#if defined MORPHODYN_MUSTANG_byHYDRO
       ALLOCATE(dhsed_save(PROC_IN_ARRAY))
#endif
       hsed0(PROC_IN_ARRAY)=0.0_rsh
       hsed_previous(PROC_IN_ARRAY)=0.0_rsh
#endif
   ENDIF

#if ! defined key_MARS
  !! declaration of the MUSTANG variables needed in the hydro model 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ALLOCATE(EROS_FLUX_s2w(ARRAY_EROS_FLUX_s2w))
   ALLOCATE(SETTL_FLUX_w2s(ARRAY_SETTL_FLUX_w2s))
   ALLOCATE(SETTL_FLUXSUM_w2s(ARRAY_SETTL_FLUXSUM_w2s))
   ALLOCATE(CORFLUX_SAND(ARRAY_CORFLUX_SAND))
   ALLOCATE(CORFLUY_SAND(ARRAY_CORFLUY_SAND))
   ALLOCATE(WATER_FLUX_INPUTS(ARRAY_WATER_FLUX_INPUTS))
   EROS_FLUX_s2w(ARRAY_EROS_FLUX_s2w)=0.0_rsh
   SETTL_FLUX_w2s(ARRAY_SETTL_FLUX_w2s)=0.0_rsh
   SETTL_FLUXSUM_w2s(ARRAY_SETTL_FLUXSUM_w2s)=0.0_rsh
   CORFLUX_SAND(ARRAY_CORFLUX_SAND )=1.0_rsh
   CORFLUY_SAND(ARRAY_CORFLUY_SAND)=1.0_rsh
   WATER_FLUX_INPUTS(ARRAY_WATER_FLUX_INPUTS)=0.0_rsh
   
   ALLOCATE(fwet(PROC_IN_ARRAY))
   ALLOCATE(fwetp(PROC_IN_ARRAY))
   fwet(:,:)=1.0_rsh
   fwetp(:,:)=1.0_rsh

   ! if temperature and salinity not ranged in the same table as substances concentrations
   !ALLOCATE(EROS_FLUX_s2w_TEMP(ARRAY_EROS_FLUX_s2w_TEMPSAL))
   !ALLOCATE(EROS_FLUX_s2w_SAL(ARRAY_EROS_FLUX_s2w_TEMPSAL))
   !EROS_FLUX_s2w_TEMP(ARRAY_EROS_FLUX_s2w_TEMPSAL)=0.0_rsh
   !EROS_FLUX_s2w_SAL(ARRAY_EROS_FLUX_s2w_TEMPSAL)=0.0_rsh

#endif
   
#ifdef key_MUSTANG_flocmod
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!! module FLOCULATION  !!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  ! floc characteristics
  ALLOCATE(f_diam(1:nv_mud))     ! floc diameter
  ALLOCATE(f_vol(1:nv_mud))      ! floc volume
  ALLOCATE(f_rho(1:nv_mud))      ! floc density
  ALLOCATE(f_mass(0:nv_mud+1))     ! floc mass
  
  ! mass concentration 
  !ALLOCATE(f_cv(1:nv_mud))       ! extracted from cv_wat(:,k,j,j) mass concentration for every mud variables

  ! agregation kernels
  ALLOCATE(f_coll_prob_sh(1:nv_mud,1:nv_mud)) !  shear agregation collision probability
  ALLOCATE(f_coll_prob_ds(1:nv_mud,1:nv_mud)) ! differential settling collision probability
  
  ALLOCATE(f_g1_sh(1:nv_mud,1:nv_mud,1:nv_mud)) ! shear agregation gain term
  ALLOCATE(f_g1_ds(1:nv_mud,1:nv_mud,1:nv_mud)) ! differential settling agregation gain term
  ALLOCATE(f_l1_sh(1:nv_mud,1:nv_mud)) ! shear agregation loss term
  ALLOCATE(f_l1_ds(1:nv_mud,1:nv_mud)) ! differential settling agregation loss term  
  ALLOCATE(f_g3(1:nv_mud,1:nv_mud)) ! fragmentation gain term     
  ALLOCATE(f_l3(1:nv_mud)) ! fragmentation loss term
    
  ALLOCATE(f_davg(1:NB_LAYER_WAT,PROC_IN_ARRAY))
  ALLOCATE(f_d50(1:NB_LAYER_WAT,PROC_IN_ARRAY))
  ALLOCATE(f_d90(1:NB_LAYER_WAT,PROC_IN_ARRAY))
  ALLOCATE(f_d10(1:NB_LAYER_WAT,PROC_IN_ARRAY))
  ALLOCATE(f_dtmin(1:NB_LAYER_WAT,PROC_IN_ARRAY))
  
  f_diam(1:nv_mud)=0.0_rsh
  f_vol(1:nv_mud)=0.0_rsh
  f_rho(1:nv_mud)=0.0_rsh
  f_mass(0:nv_mud+1)=0.0_rsh
  
  !f_cv(1:nv_mud)=0.0_rsh
  
  f_coll_prob_sh(1:nv_mud,1:nv_mud)=0.0_rsh
  f_coll_prob_ds(1:nv_mud,1:nv_mud)=0.0_rsh
  
  f_g1_sh(1:nv_mud,1:nv_mud,1:nv_mud)=0.0_rsh
  f_g1_ds(1:nv_mud,1:nv_mud,1:nv_mud)=0.0_rsh
  f_l1_sh(1:nv_mud,1:nv_mud)=0.0_rsh
  f_l1_ds(1:nv_mud,1:nv_mud)=0.0_rsh
  f_g3(1:nv_mud,1:nv_mud)=0.0_rsh
  f_l3(1:nv_mud)=0.0_rsh

  f_davg(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
  f_d50(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
  f_d10(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
  f_d90(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
  f_dtmin(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh  

#endif
   
   PRINT_DBG*, 'END MUSTANG_ALLOC'
   
   END SUBROUTINE MUSTANG_alloc

  !!===========================================================================

#endif

END MODULE comMUSTANG

