C     -*- fortran -*-
      integer ad_array_size
      parameter (ad_array_size=1 + (3+Lm+3+padd_X)*(3+Mm+3+padd_E))

      integer NT
      parameter (NT = 8641)

      real*4 obs(GLOBAL_2D_ARRAY, NT)
      
      double precision x(ad_array_size)
      double precision x0(ad_array_size)

      integer sim_iicroot
      
      common /obs_data/ obs

      common /state_info/ sim_iicroot
