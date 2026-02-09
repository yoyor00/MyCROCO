/* This is include file "das_covar.h". 
  --------------------------------------------
*/
! z-coord
      real bz_das(GLOBAL_2D_ARRAY)
      real bu_das(GLOBAL_2D_ARRAY,ndas)
      real bv_das(GLOBAL_2D_ARRAY,ndas)
      real bpsi_das(GLOBAL_2D_ARRAY,ndas)
      real bchi_das(GLOBAL_2D_ARRAY,ndas)
      real bt_das(GLOBAL_2D_ARRAY,ndas,NT)
      common 
     &     /ocean_bu_das/bu_das /ocean_bv_das/bv_das
     &     /ocean_bpsi_das/bpsi_das /ocean_bchi_das/bchi_das
     &     /ocean_bt_das/bt_das /ocean_bz_das/bz_das
!
! zeta
!
      real czE_das(GLOBAL_1D_ETA,GLOBAL_1D_ETA)
      real czX_das(GLOBAL_1D_XI,GLOBAL_1D_XI)
!
! u & v
!
      real cuE_das(GLOBAL_1D_ETA,GLOBAL_1D_ETA)
      real cuX_das(GLOBAL_1D_XI,GLOBAL_1D_XI)
      real cuZ_das(NDAS,NDAS)
      real cvE_das(GLOBAL_1D_ETA,GLOBAL_1D_ETA)
      real cvX_das(GLOBAL_1D_XI,GLOBAL_1D_XI)
      real cvZ_das(NDAS,NDAS)
!
! psi & chi
!
      real cpsiE_das(GLOBAL_1D_ETA,GLOBAL_1D_ETA)
      real cpsiX_das(GLOBAL_1D_XI,GLOBAL_1D_XI)
      real cpsiZ_das(NDAS,NDAS)
      real cchiE_das(GLOBAL_1D_ETA,GLOBAL_1D_ETA)
      real cchiX_das(GLOBAL_1D_XI,GLOBAL_1D_XI)
      real cchiZ_das(NDAS,NDAS)
!
! temp and salt
!
      real ctE_das(GLOBAL_1D_ETA,GLOBAL_1D_ETA,NT)
      real ctX_das(GLOBAL_1D_XI,GLOBAL_1D_XI,NT)
      real ctZ_das(ndas,ndas,NT)
!
      common /ocean_cze_das/czE_das /ocean_czX_das/czX_das
     &       /ocean_cue_das/cuE_das /ocean_cuX_das/cuX_das
     &       /ocean_cuz_das/cuZ_das
     &       /ocean_cve_das/cvE_das /ocean_cvX_das/cvX_das
     &       /ocean_cvz_das/cvZ_das
     &       /ocean_cpsie_das/cpsiE_das
     &       /ocean_cpsix_das/cpsiX_das
     &       /ocean_cpsiz_das/cpsiZ_das
     &       /ocean_cchie_das/cchiE_das
     &       /ocean_cchix_das/cchiX_das
     &       /ocean_cchiz_das/cchiZ_das
     &       /ocean_cte_das/ctE_das /ocean_ctX_das/ctX_das
     &       /ocean_ctz_das/ctZ_das
!
! cross T/S
!
      real corr_ts(GLOBAL_2D_ARRAY,ndas)
      common /ocean_cross_ts/corr_ts

! variance and correlation files
      character*250 file_corr
      common /corr_file/file_corr
      character*250 file_corr_vert
      common /corr_file_vert/file_corr_vert
      character*250 file_bvar
      common /bvar_file_vert/file_bvar
