C     -*- fortran -*-
      integer ad_array_size
      parameter (ad_array_size=1)

C     size of the assimilation window (number of steps)
      integer ad_nt
      parameter (ad_nt = 7200)

C     start of assimilation
      integer ad_ast
      parameter (ad_ast = 7200)

      double precision ad_obs(GLOBAL_2D_ARRAY, 20000)
      
      double precision ad_x(ad_array_size)
      double precision ad_xd(ad_array_size)

      integer sim_iicroot
      
      common /obs_data/ ad_obs

      common /state_info/ sim_iicroot
