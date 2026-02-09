/* This is include file "das_ocean9.h". 
  --------------------------------------------
*/
      real zeta_sm(GLOBAL_2D_ARRAY)
      common /ocean_z_smooth/zeta_sm
      real zeta_sm_adj(GLOBAL_2D_ARRAY)
      common /adj_ocean_z_smooth/zeta_sm_adj
#if defined DAS_PWTHZETATOT
      real zeta_h_sm(GLOBAL_2D_ARRAY)
      common /ocean_zh_smooth/zeta_h_sm
      real zeta_h_sm_adj(GLOBAL_2D_ARRAY)
      common /adj_ocean_zh_smooth/zeta_h_sm_adj
#endif
/**** 3d *****/
      real p_sm(GLOBAL_2D_ARRAY,ndas)
      real p_sm_adj(GLOBAL_2D_ARRAY,ndas)
      real psi_sm(GLOBAL_2D_ARRAY,ndas)
      real chi_sm(GLOBAL_2D_ARRAY,ndas)
      common /ocean_p_smooth/p_sm /ocean_padj_smooth/p_sm_adj
     &       /ocean_psi_smooth/psi_sm /ocean_chi_smooth/chi_sm
      real t_sm(GLOBAL_2D_ARRAY,ndas,NT)
      common /ocean_t_smooth/t_sm 
      real t_sm_adj(GLOBAL_2D_ARRAY,ndas,NT)
      common /adj_ocean_t_smooth/t_sm_adj
      real u_sm(GLOBAL_2D_ARRAY,ndas)
      common /ocean_u_smooth/u_sm
      real v_sm(GLOBAL_2D_ARRAY,ndas)
      common /ocean_v_smooth/v_sm
