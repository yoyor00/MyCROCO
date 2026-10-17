!=======================================================================
!  croco_namelist_agrif.h: Multi-grid storage definition for AGRIF
!=======================================================================
   type t_croco_nml_grid
      character(len=200) :: fname_nml
      character(len=64) :: testcase_name
      logical :: lts_hadv
      logical :: lts_rot_or_diag
      logical :: lts_per
      logical :: lekman_apg
      logical :: loscnonrot
      logical :: lnonrot
      real :: UG_forced
      real :: wave_ramp_duration
      character(len=80) :: title
      character(len=180) :: logname
      real :: dt
      integer :: ntimes
      integer :: ndtfast
      integer :: ninfo
      real :: m2filter_alpha
      logical :: ldefhis
      integer :: nwrt
      integer :: nrpfhis
      character(len=180) :: hisname
      integer :: nrrec
      character(len=180) :: ininame
      integer :: nrst
      integer :: nrpfrst
      character(len=180) :: rstname
#if defined NBQ || defined K3FAST   /* H1-fix 2026-08-26: directive invalide (#ifdef NBQ || ...) */
      integer :: ndtnbq_in
      real :: csound_nbq
      real :: visc2read_nbq
#endif
#ifdef SOLVE3D
      real :: theta_s
      real :: theta_b
      real :: hc
      real :: Tcline
#endif
#ifdef USE_CALENDAR
      character(len=19) :: start_date
      character(len=19) :: end_date
      real :: dt_his
      real :: dt_avg
      real :: dt_rst
#endif
#ifndef ANA_GRID
      character(len=180) :: grdname
#endif
      character(len=180) :: frcname
      logical :: use_frcname
      character(len=180) :: bulkname
#if (  defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || (defined ZCLIMATOLOGY&& !defined ANA_SSH) || (defined M2CLIMATOLOGY&& !defined ANA_M2CLIMA) || (defined M3CLIMATOLOGY&& !defined ANA_M3CLIMA)
      character(len=180) :: clmname
#endif
#if !defined ANA_BRY && defined FRC_BRY
      character(len=180) :: bry_file
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      character(len=180) :: brywkb_file
#endif
#ifdef WKB_WWAVE
      real :: wkb_amp
      real :: wkb_ang
      real :: wkb_prd
      real :: wkb_tide
      real :: wkb_btg
      real :: wkb_gam
#  ifdef WAVE_ROLLER
      real :: wkb_rsb
      real :: wkb_roller
#  endif
#endif
#ifdef WAVE_MAKER
      real :: wmaker_amp
      real :: wmaker_prd
      real :: wmaker_dir
      real :: wmaker_dsp
      real :: wmaker_fsp
      real :: wmaker_wf1
      real :: wmaker_wf2
      real :: wmaker_wa1
      real :: wmaker_wa2
#endif
#ifdef AVERAGES
      integer :: ntsavg
      integer :: navg
      integer :: nrpfavg
      character(len=180) :: avgname
#endif
#if defined OUTPUTS_SURFACE && !defined XIOS
      logical :: ldefsurf
      integer :: nwrtsurf
      integer :: nrpfsurf
      character(len=180) :: surfname
      logical :: out_his_surf
#  ifdef AVERAGES
      logical :: ldefsurf_avg
      integer :: ntssurf_avg
      integer :: nwrtsurf_avg
      integer :: nrpfsurf_avg
      character(len=180) :: surfname_avg
      logical :: out_avg_surf
#  endif
#endif
#if defined ABL1D && !defined XIOS
      logical :: ldefablhis
      integer :: nwrtablhis
      integer :: nrpfablhis
      character(len=180) :: ablname
#  ifdef AVERAGES
      logical :: ldefablavg
      integer :: ntsablavg
      integer :: nwrtablavg
      integer :: nrpfablavg
      character(len=180) :: ablname_avg
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      character(len=180) :: wave_file
#endif
#if defined BIOLOGY && defined PISCES
      character(len=180) :: bioname
#endif
#ifdef BODYFORCE
      integer :: levsfrc
      integer :: levbfrc
#endif
      real :: R0
      real :: T0
      real :: S0
      real :: Tcoef
      real :: Scoef
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
      real :: ltra_min
      real :: ltra_max
#endif
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
      real :: ldyn_min
      real :: ldyn_max
#endif
#ifdef SEDIMENT
      character(len=180) :: sedname
#endif
#ifdef MUSTANG
      character(len=180) :: sedname_must
#endif
#ifdef SUBSTANCE
      character(len=180) :: subsfilename
#endif
#ifdef OBSTRUCTION
      character(len=180) :: obstname
#endif
#ifdef XIOS
      character(len=80) :: xios_origin_date
#endif
#ifdef ASSIMILATION
      character(len=180) :: aparnam
      character(len=180) :: assname
#endif
      real :: rho0
      real :: rdrg
      real :: rdrg2
      real :: Zobt
      real :: Cdb_min
      real :: Cdb_max
      real :: gamma2
      real :: visc2
      real :: visc4
#ifdef SOLVE3D
#  ifdef TRACERS
      real, allocatable :: tnu2(:)
      real, allocatable :: tnu4(:)
#  endif
#  if !defined LMD_MIXING
      real :: Akv_bak
#    ifdef TRACERS
      real, allocatable :: Akt_bak(:)
#    endif
#  endif
#endif
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
      integer :: psource_Nsrc
      integer, allocatable :: psource_Isrc(:)
      integer, allocatable :: psource_Jsrc(:)
      integer, allocatable :: psource_Dsrc(:)
      character(len=180) :: psource_qbarname
      real, allocatable :: psource_qbardir(:)
      real, allocatable :: psource_Qbar(:)
# if defined TRACERS
      logical, allocatable :: psource_Lsrc(:, :)
      real, allocatable :: psource_Tsrc0(:, :)
# endif
#endif
#if defined DIAGNOSTICS_TS
      logical :: ldefdia
      integer :: nwrtdia
      integer :: nrpfdia
      character(len=180) :: dianame
#  ifdef AVERAGES
      logical :: ldefdia_avg
      integer :: ntsdia_avg
      integer :: nwrtdia_avg
      integer :: nrpfdia_avg
      character(len=180) :: dianame_avg
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
      real :: mld_crit_D
      real :: mld_crit_T
#  endif
#endif
#if defined DIAGNOSTICS_UV
      logical :: ldefdiaM
      integer :: nwrtdiaM
      integer :: nrpfdiaM
      character(len=180) :: dianameM
#  ifdef AVERAGES
      logical :: ldefdiaM_avg
      integer :: ntsdiaM_avg
      integer :: nwrtdiaM_avg
      integer :: nrpfdiaM_avg
      character(len=180) :: dianameM_avg
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      logical :: ldefdiags_vrt
      integer :: nwrtdiags_vrt
      integer :: nrpfdiags_vrt
      character(len=180) :: diags_vrtname
#  ifdef AVERAGES
      logical :: ldefdiags_vrt_avg
      integer :: ntsdiags_vrt_avg
      integer :: nwrtdiags_vrt_avg
      integer :: nrpfdiags_vrt_avg
      character(len=180) :: diags_vrtname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_KE
      logical :: ldefdiags_ek
      integer :: nwrtdiags_ek
      integer :: nrpfdiags_ek
      character(len=180) :: diags_ekname
#  ifdef AVERAGES
      logical :: ldefdiags_ek_avg
      integer :: ntsdiags_ek_avg
      integer :: nwrtdiags_ek_avg
      integer :: nrpfdiags_ek_avg
      character(len=180) :: diags_ekname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      logical :: ldefdiags_pv
      integer :: nwrtdiags_pv
      integer :: nrpfdiags_pv
      character(len=180) :: diags_pvname
#  ifdef AVERAGES
      logical :: ldefdiags_pv_avg
      integer :: ntsdiags_pv_avg
      integer :: nwrtdiags_pv_avg
      integer :: nrpfdiags_pv_avg
      character(len=180) :: diags_pvname_avg
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
#  ifdef AVERAGES
      logical :: ldefdiags_eddy_avg
      integer :: ntsdiags_eddy_avg
      integer :: nwrtdiags_eddy_avg
      integer :: nrpfdiags_eddy_avg
      character(len=180) :: diags_eddyname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      logical :: ldefdiabio
      integer :: nwrtdiabio
      integer :: nrpfdiabio
      character(len=180) :: dianamebio
#  ifdef AVERAGES
      logical :: ldefdiabio_avg
      integer :: ntsdiabio_avg
      integer :: nwrtdiabio_avg
      integer :: nrpfdiabio_avg
      character(len=180) :: dianamebio_avg
#  endif
#endif
#ifdef STATIONS
      logical :: ldefsta
      integer :: nsta
      integer :: nrpfsta
      character(len=180) :: staposname
      character(len=180) :: staname
      logical :: sta_grd
      logical :: sta_temp
      logical :: sta_salt
      logical :: sta_rho
      logical :: sta_vel
#endif
#if (defined SPONGE && !defined SPONGE_GRID)
      real :: x_sponge
#endif
#ifdef WET_DRY
      real :: D_wetdry
#endif
#if defined MRL_WCI && defined WAVE_DRY
      real :: D_wavedry
#endif
#if defined T_FRC_BRY     || defined M2_FRC_BRY    || defined M3_FRC_BRY||defined Z_FRC_BRY|| defined W_FRC_BRY||defined NBQ_FRC_BRY|| defined TCLIMATOLOGY||defined M2CLIMATOLOGY|| defined M3CLIMATOLOGY||defined ZCLIMATOLOGY|| defined WCLIMATOLOGY||defined NBQCLIMATOLOGY|| defined TNUDGING||defined ZNUDGING|| defined M2NUDGING||defined M3NUDGING|| defined NBQ_NUDGING||defined SPONGE
      real :: tauT_in
      real :: tauT_out
      real :: tauM_in
      real :: tauM_out
#endif
#if defined BHFLUX || (defined BWFLUX && defined SALINITY)
      character(len=180) :: btfname
#endif
      logical :: out_his_zeta
      logical :: out_his_ubar
      logical :: out_his_vbar
#ifdef SOLVE3D
      logical :: out_his_u
      logical :: out_his_v
# ifdef TRACERS
      logical, allocatable :: out_his_tracer(:)
# endif
#endif
#ifdef AVERAGES
      logical :: out_avg_zeta
      logical :: out_avg_ubar
      logical :: out_avg_vbar
# ifdef SOLVE3D
      logical :: out_avg_u
      logical :: out_avg_v
#  ifdef TRACERS
      logical, allocatable :: out_avg_tracer(:)
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_TS
# ifdef TRACERS
      logical, allocatable :: out_his_dia3D_tracer(:)
#  ifdef DIAGNOSTICS_TS_MLD
      logical, allocatable :: out_his_dia2D_tracer(:)
#  endif
# endif
# if defined AVERAGES && defined TRACERS
      logical, allocatable :: out_avg_dia3D_tracer(:)
#  ifdef DIAGNOSTICS_TS_MLD
      logical, allocatable :: out_avg_dia2D_tracer(:)
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_UV
      logical :: out_his_diagM_u
      logical :: out_his_diagM_v
# ifdef AVERAGES
      logical :: out_avg_diagM_u
      logical :: out_avg_diagM_v
# endif
#endif
#ifdef DIAGNOSTICS_VRT
      logical :: out_his_diags_vrt
# ifdef AVERAGES
      logical :: out_avg_diags_vrt
# endif
#endif
#ifdef DIAGNOSTICS_KE
      logical :: out_his_diags_ek
# ifdef AVERAGES
      logical :: out_avg_diags_ek
# endif
#endif
#if defined DIAGNOSTICS_PV && defined TRACERS
      logical, allocatable :: out_his_diags_pv_tracer(:)
# ifdef AVERAGES
      logical, allocatable :: out_avg_diags_pv_tracer(:)
# endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
# ifdef AVERAGES
      logical :: out_avg_diags_eddy
# endif
#endif
#ifdef DIAGNOSTICS_BIO
      logical, allocatable :: out_his_diagbioFlux(:)
      logical, allocatable :: out_his_diagbioVSink(:)
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      logical, allocatable :: out_his_diagbioGasExc(:)
# endif
# ifdef AVERAGES
      logical, allocatable :: out_avg_diagbioFlux(:)
      logical, allocatable :: out_avg_diagbioVSink(:)
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      logical, allocatable :: out_avg_diagbioGasExc(:)
#  endif
# endif
#endif
#ifdef ONLINE
      integer :: yearnum
      integer :: monthnum
      integer :: recordsperday
      integer :: yearend
      integer :: monthend
      character(len=250) :: pathbulk
#endif
#ifdef STOGEN
      logical :: out_his_xi2d
      logical :: out_his_xi3d
#endif
#if defined SOLVE3D && defined GLS_MIXING
      logical :: out_his_tke
      logical :: out_his_gls
      logical :: out_his_lscale
# ifdef AVERAGES
      logical :: out_avg_tke
      logical :: out_avg_gls
      logical :: out_avg_lscale
# endif
#endif
#if defined SOLVE3D && defined SEDIMENT
      logical :: out_his_sed_athk
      logical :: out_his_sed_bthk
      logical :: out_his_sed_bpor
      logical, allocatable :: out_his_sed_bfra(:)
# ifdef SUSPLOAD
      logical, allocatable :: out_his_sed_dflx(:)
      logical, allocatable :: out_his_sed_eflx(:)
# endif
# ifdef BEDLOAD
      logical, allocatable :: out_his_sed_bdlu(:)
      logical, allocatable :: out_his_sed_bdlv(:)
# endif
# if defined MIXED_BED || defined COHESIVE_BED
      logical :: out_his_sed_btcr
# endif
#endif
#ifdef BBL
      logical :: out_his_abed
      logical :: out_his_hripple
      logical :: out_his_lripple
      logical :: out_his_zbnot
      logical :: out_his_zbapp
      logical :: out_his_bostrw
#endif
#ifdef MRL_WCI
      logical :: out_his_sup
      logical :: out_his_ust2d
      logical :: out_his_vst2d
# ifdef SOLVE3D
      logical :: out_his_ust
      logical :: out_his_vst
      logical :: out_his_wst
      logical :: out_his_akb
      logical :: out_his_akw
      logical :: out_his_kvf
      logical :: out_his_calp
      logical :: out_his_kaps
# endif
# ifdef AVERAGES
      logical :: out_avg_sup
      logical :: out_avg_ust2d
      logical :: out_avg_vst2d
#  ifdef SOLVE3D
      logical :: out_avg_ust
      logical :: out_avg_vst
      logical :: out_avg_wst
      logical :: out_avg_akb
      logical :: out_avg_akw
      logical :: out_avg_kvf
      logical :: out_avg_calp
      logical :: out_avg_kaps
#  endif
# endif
#endif
#if defined MRL_WCI || defined OW_COUPLING
      logical :: out_his_hrm
      logical :: out_his_frq
      logical :: out_his_action
      logical :: out_his_k_xi
      logical :: out_his_k_eta
      logical :: out_his_eps_b
      logical :: out_his_eps_d
      logical :: out_his_erol
      logical :: out_his_eps_r
# ifdef AVERAGES
      logical :: out_avg_hrm
      logical :: out_avg_frq
      logical :: out_avg_action
      logical :: out_avg_k_xi
      logical :: out_avg_k_eta
      logical :: out_avg_eps_b
      logical :: out_avg_eps_d
      logical :: out_avg_erol
      logical :: out_avg_eps_r
# endif
#endif
#if defined ABL1D && !defined XIOS
      logical :: out_his_abl_pu_dta
      logical :: out_his_abl_pv_dta
      logical :: out_his_abl_pt_dta
      logical :: out_his_abl_pq_dta
      logical :: out_his_abl_pgu_dta
      logical :: out_his_abl_pgv_dta
      logical :: out_his_abl_u_abl
      logical :: out_his_abl_v_abl
      logical :: out_his_abl_t_abl
      logical :: out_his_abl_q_abl
      logical :: out_his_abl_tke_abl
      logical :: out_his_abl_mxlm_abl
      logical :: out_his_abl_mxld_abl
      logical :: out_his_abl_avm_abl
      logical :: out_his_abl_avt_abl
      logical :: out_his_abl_ablh_abl
      logical :: out_his_abl_zr_abl
      logical :: out_his_abl_zw_abl
      logical :: out_his_abl_Hzr_abl
      logical :: out_his_abl_Hzw_abl
# ifdef AVERAGES
      logical :: out_avg_abl_pu_dta
      logical :: out_avg_abl_pv_dta
      logical :: out_avg_abl_pt_dta
      logical :: out_avg_abl_pq_dta
      logical :: out_avg_abl_pgu_dta
      logical :: out_avg_abl_pgv_dta
      logical :: out_avg_abl_u_abl
      logical :: out_avg_abl_v_abl
      logical :: out_avg_abl_t_abl
      logical :: out_avg_abl_q_abl
      logical :: out_avg_abl_tke_abl
      logical :: out_avg_abl_mxlm_abl
      logical :: out_avg_abl_mxld_abl
      logical :: out_avg_abl_avm_abl
      logical :: out_avg_abl_avt_abl
      logical :: out_avg_abl_ablh_abl
      logical :: out_avg_abl_zr_abl
      logical :: out_avg_abl_zw_abl
      logical :: out_avg_abl_Hzr_abl
      logical :: out_avg_abl_Hzw_abl
# endif
#endif
# if defined SOLVE3D || defined RIP
      logical :: out_his_rho
      logical :: out_his_omega
      logical :: out_his_w
      logical :: out_his_akv
      logical :: out_his_bostr
      logical :: out_his_bustr
      logical :: out_his_bvstr
      logical :: out_his_wstr
      logical :: out_his_ustr
      logical :: out_his_vstr
# endif
# if defined SOLVE3D && defined TEMPERATURE
      logical :: out_his_akt
      logical :: out_his_shflx
      logical :: out_his_shflx_rsw
# endif
# if defined SOLVE3D && defined SALINITY
      logical :: out_his_aks
      logical :: out_his_swflx
# endif
# if defined SOLVE3D && defined BULK_FLUX
      logical :: out_his_shflx_rlw
      logical :: out_his_shflx_lat
      logical :: out_his_shflx_sen
# endif
# if defined SOLVE3D && defined BHFLUX
      logical :: out_his_bhflx
# endif
# if defined SOLVE3D && defined BWFLUX && defined SALINITY
      logical :: out_his_bwflx
# endif
# if defined SOLVE3D && (defined ANA_VMIX||defined LMD_MIXING|| defined LMD_SKPP||defined LMD_BKPP|| defined GLS_MIXING)
      logical :: out_his_bvf
# endif
# if defined SOLVE3D && (defined LMD_SKPP || defined GLS_MIXING)
      logical :: out_his_hbl
# endif
# if defined SOLVE3D && defined LMD_BKPP
      logical :: out_his_hbbl
# endif
# if defined SOLVE3D && defined VIS_COEF_3D
      logical :: out_his_visc3d
# endif
# if defined SOLVE3D && defined DIF_COEF_3D
      logical :: out_his_diff3d
# endif
# if defined SOLVE3D && defined BIOLOGY && !defined PISCES
      logical :: out_his_hel
# endif
# if defined SOLVE3D && defined BIO_NChlPZD
      logical :: out_his_chc
      logical :: out_his_u10
      logical :: out_his_kvo2
      logical :: out_his_o2sat
# endif
# if defined SOLVE3D && defined BIO_BioEBUS
      logical :: out_his_aou
      logical :: out_his_wind10
# endif
# if defined SOLVE3D && defined MORPHODYN
      logical :: out_his_hm
# endif
# if defined AVERAGES && defined SOLVE3D
      logical :: out_avg_rho
      logical :: out_avg_omega
      logical :: out_avg_w
      logical :: out_avg_akv
      logical :: out_avg_bostr
      logical :: out_avg_bustr
      logical :: out_avg_bvstr
      logical :: out_avg_wstr
      logical :: out_avg_ustr
      logical :: out_avg_vstr
# endif
# if defined AVERAGES && defined SOLVE3D && defined TEMPERATURE
      logical :: out_avg_akt
      logical :: out_avg_shflx
      logical :: out_avg_shflx_rsw
# endif
# if defined AVERAGES && defined SOLVE3D && defined SALINITY
      logical :: out_avg_aks
      logical :: out_avg_swflx
# endif
# if defined AVERAGES && defined SOLVE3D && defined BULK_FLUX
      logical :: out_avg_shflx_rlw
      logical :: out_avg_shflx_lat
      logical :: out_avg_shflx_sen
# endif
# if defined AVERAGES && defined SOLVE3D && defined BHFLUX
      logical :: out_avg_bhflx
# endif
# if defined AVERAGES && defined SOLVE3D && defined BWFLUX && defined SALINITY
      logical :: out_avg_bwflx
# endif
# if defined AVERAGES && defined SOLVE3D && (defined ANA_VMIX||defined LMD_MIXING|| defined LMD_SKPP||defined LMD_BKPP|| defined GLS_MIXING)
      logical :: out_avg_bvf
# endif
# if defined AVERAGES && defined SOLVE3D && (defined LMD_SKPP||defined GLS_MIXING)
      logical :: out_avg_hbl
# endif
# if defined AVERAGES && defined SOLVE3D && defined LMD_BKPP
      logical :: out_avg_hbbl
# endif
# if defined AVERAGES && defined SOLVE3D && defined VIS_COEF_3D
      logical :: out_avg_visc3d
# endif
# if defined AVERAGES && defined SOLVE3D && defined DIF_COEF_3D
      logical :: out_avg_diff3d
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIOLOGY && !defined PISCES
      logical :: out_avg_hel
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_NChlPZD
      logical :: out_avg_chc
      logical :: out_avg_u10
      logical :: out_avg_kvo2
      logical :: out_avg_o2sat
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_BioEBUS
      logical :: out_avg_aou
      logical :: out_avg_wind10
# endif
# if defined AVERAGES && defined SOLVE3D && defined MORPHODYN
      logical :: out_avg_hm
# endif
   end type t_croco_nml_grid

   integer, parameter :: MAX_AGRIF_GRIDS = 10
   type(t_croco_nml_grid), save :: agrif_nml_storage(0:MAX_AGRIF_GRIDS)
   logical, save :: agrif_nml_stored(0:MAX_AGRIF_GRIDS) = .false.
   integer, save :: agrif_nml_current_grid = -1

   public :: croco_namelist_save_grid
   public :: croco_namelist_restore_grid
   public :: croco_namelist_switch_grid

