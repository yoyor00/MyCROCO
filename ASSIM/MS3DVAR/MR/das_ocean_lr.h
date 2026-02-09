/* This is include file "das_ocean_lr.h". 
--------------------------------------------
*/
! z-coord
      real zeta_lr(GLOBAL_2D_ARRAY_LR)
      real u_lr(GLOBAL_2D_ARRAY_LR,ndas)
      real v_lr(GLOBAL_2D_ARRAY_LR,ndas)
      real t_lr(GLOBAL_2D_ARRAY_LR,ndas,NT)
      real lat_lr(GLOBAL_2D_ARRAY_LR)
      real lon_lr(GLOBAL_2D_ARRAY_LR)

      common /ocean_u_lr/u_lr /ocean_v_lr/v_lr
     &     /ocean_t_lr/t_lr /ocean_zeta_lr/zeta_lr
     &     /ocean_lat_lr/lat_lr /ocean_lon_lr/lon_lr

      real zeta_lr_s(GLOBAL_2D_ARRAY)
      real u_lr_s(GLOBAL_2D_ARRAY,ndas)
      real v_lr_s(GLOBAL_2D_ARRAY,ndas)
      real t_lr_s(GLOBAL_2D_ARRAY,ndas,NT)
      common /ocean_u_lr_s/u_lr_s /ocean_v_lr_s/v_lr_s
     &     /ocean_t_lr_s/t_lr_s /ocean_zeta_lr_s/zeta_lr_s

      integer flag_lr_da
      common /ocean_flag_lr_da/flag_lr_da
