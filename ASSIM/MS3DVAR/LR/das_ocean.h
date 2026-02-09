/* This is include file "das_ocean.h". 
  --------------------------------------------
*/
! z-coord
      real zeta_das(GLOBAL_2D_ARRAY)
      real u_das(GLOBAL_2D_ARRAY,ndas)
      real v_das(GLOBAL_2D_ARRAY,ndas)
      real psi_das(GLOBAL_2D_ARRAY,ndas)
      real chi_das(GLOBAL_2D_ARRAY,ndas)
      real t_das(GLOBAL_2D_ARRAY,ndas,NT)
      common /ocean_u_das/u_das /ocean_v_das/v_das
     &    /ocean_psi_das/psi_das /ocean_chi_das/chi_das
     &    /ocean_t_das/t_das /ocean_zeta_das/zeta_das
!
! intermediately saved arrays
!
      real zeta_s(GLOBAL_2D_ARRAY)
      real zeta_h(GLOBAL_2D_ARRAY)
      real u_s(GLOBAL_2D_ARRAY,ndas)
      real v_s(GLOBAL_2D_ARRAY,ndas)
      real psi_s(GLOBAL_2D_ARRAY,ndas)
      real chi_s(GLOBAL_2D_ARRAY,ndas)
      real rho_s(GLOBAL_2D_ARRAY,ndas)
      real p_s(GLOBAL_2D_ARRAY,ndas)
      real t_s(GLOBAL_2D_ARRAY,ndas,NT)
      real t_w(GLOBAL_2D_ARRAY,ndas,NT)
      common /ocean_u_s/u_s /ocean_v_s/v_s
     &      /ocean_psi_s/psi_s /ocean_chi_s/chi_s
     &      /ocean_rho_s/rho_s /ocean_p_s/p_s
     &     /ocean_t_s/t_s /ocean_zeta_s/zeta_s
     &     /ocean_zeta_h/zeta_h /ocean_t_w/t_w
!
! adjoint variables
!
      real zeta_adj(GLOBAL_2D_ARRAY)
      real psi_adj(GLOBAL_2D_ARRAY,ndas)
      real chi_adj(GLOBAL_2D_ARRAY,ndas)
      real t_adj(GLOBAL_2D_ARRAY,ndas,NT)
      common
     &     /ocean_psi_adj/psi_adj /ocean_chi_adj/chi_adj
     &     /ocean_t_adj/t_adj /ocean_zeta_adj/zeta_adj
!
      real zeta_s_adj(GLOBAL_2D_ARRAY)
      real zeta_h_adj(GLOBAL_2D_ARRAY)
      real u_s_adj(GLOBAL_2D_ARRAY,ndas)
      real v_s_adj(GLOBAL_2D_ARRAY,ndas)
      real psi_s_adj(GLOBAL_2D_ARRAY,ndas)
      real chi_s_adj(GLOBAL_2D_ARRAY,ndas)
      real rho_s_adj(GLOBAL_2D_ARRAY,ndas)
      real p_s_adj(GLOBAL_2D_ARRAY,ndas)
      real t_s_adj(GLOBAL_2D_ARRAY,ndas,NT)
      real t_w_adj(GLOBAL_2D_ARRAY,ndas,NT)
      common /ocean_u_s_adj/u_s_adj /ocean_v_s_adj/v_s_adj
     &   /ocean_psi_s_adj/psi_s_adj /ocean_chi_s_adj/chi_s_adj
     &     /ocean_rho_s_adj/rho_s_adj /ocean_p_s_adj/p_s_adj
     &     /ocean_t_s_adj/t_s_adj /ocean_zeta_s_adj/zeta_s_adj
     &     /ocean_zeta_h_adj/zeta_h_adj /ocean_t_w_adj/t_w_adj

#ifdef MASKING
      real rmask_das(GLOBAL_2D_ARRAY,ndas)
      real pmask_das(GLOBAL_2D_ARRAY,ndas)
      real umask_das(GLOBAL_2D_ARRAY,ndas)
      real vmask_das(GLOBAL_2D_ARRAY,ndas)
      common /mask_r_das/rmask_das /mask_p_das/pmask_das
     &       /mask_u_das/umask_das /mask_v_das/vmask_das
#endif

!
! for interpolation from s-coord to z-coord
!
! z_das(k) is in water when k >= nz_das  
!
      integer nzr_das(GLOBAL_2D_ARRAY)
      integer nzu_das(GLOBAL_2D_ARRAY)
      integer nzv_das(GLOBAL_2D_ARRAY)
      common /ocean_nzr_das/nzr_das
     &       /ocean_nzu_das/nzu_das /ocean_nzv_das/nzv_das
!
      integer nsr_das(GLOBAL_2D_ARRAY)
      integer nsu_das(GLOBAL_2D_ARRAY)
      integer nsv_das(GLOBAL_2D_ARRAY)
      common /ocean_nsr_das/nsr_das
     &       /ocean_nsu_das/nsu_das /ocean_nsv_das/nsv_das
!
! Geo ratio
!     
      real georatio(GLOBAL_2D_ARRAY)
      common /cgeoratio/georatio

/**** 1d *****/
      real z_das(ndas),tc1d_das(ndas),sc1d_das(ndas)
      common /grid_z_das/z_das,tc1d_das,sc1d_das  

/**** some work array*****/
!      real s1d(N),z1d(NDAS),xs1d(N),xz1d(NDAS)
!      common /work_das/s1d,z1d,xs1d,xz1d
