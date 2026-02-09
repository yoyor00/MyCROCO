/* This is include file "das_ocean.h". 
  --------------------------------------------
*/
! z-coord
      real zeta_test(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE zeta_test(BLOCK_PATTERN) BLOCK_CLAUSE
      real u_test(GLOBAL_2D_ARRAY,ndas)
CSDISTRIBUTE_RESHAPE u_ads(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real v_test(GLOBAL_2D_ARRAY,ndas)
CSDISTRIBUTE_RESHAPE v_test(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real t_test(GLOBAL_2D_ARRAY,ndas,NT)
CSDISTRIBUTE_RESHAPE t_test(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      common /ocean_u_test/u_test /ocean_v_test/v_test
     &     /ocean_t_test/t_test /ocean_zeta_test/zeta_test
