/* This is include file "das_ocean_diag.h". 
  --------------------------------------------
*/
! intermediately saved arrays
!

      real u_ap(GLOBAL_2D_ARRAY,ndas)
      real v_ap(GLOBAL_2D_ARRAY,ndas)
      real u_ac(GLOBAL_2D_ARRAY,ndas)
      real v_ac(GLOBAL_2D_ARRAY,ndas)
      common /ocean_u_a/u_ap,u_ac  
     &       /ocean_v_a/v_ap,v_ac
