/* This is include file "das_ocean9.h". 
  --------------------------------------------
*/
/**** 3d *****/
      real zeta_das9(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE zeta_das9(BLOCK_PATTERN) BLOCK_CLAUSE
!      real u_das9(GLOBAL_2D_ARRAY,ndas)
!CSDISTRIBUTE_RESHAPE u_ads9(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      real v_das9(GLOBAL_2D_ARRAY,ndas)
!CSDISTRIBUTE_RESHAPE v_das9(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real t_das9(GLOBAL_2D_ARRAY,ndas,NT)
CSDISTRIBUTE_RESHAPE t_das9(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real rho_das9(GLOBAL_2D_ARRAY,ndas)
CSDISTRIBUTE_RESHAPE rho_das9(BLOCK_PATTERN,*) BLOCK_CLAUSE
#if defined DAS_BKGZETA
CSDISTRIBUTE_RESHAPE zetab_das9(BLOCK_PATTERN) BLOCK_CLAUSE
#endif
      common /ocean_t_das9/t_das9 /ocean_zeta_das9/zeta_das9
     &       /ocean_rho_das9/rho_das9 
!     &       /ocean_u_das9/u_das9 /ocean_v_das9/v_das9
      real zetab_das9(GLOBAL_2D_ARRAY)
#if defined DAS_BKGZETA
      common /ocean_zetab_das9/zetab_das9
#endif
