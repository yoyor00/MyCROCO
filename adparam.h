C     -*- fortran -*-
      integer ad_array_size
      parameter (ad_array_size=1)

      integer NT
      parameter (NT = 8641)

      double precision obs(GLOBAL_2D_ARRAY, NT)
      
      double precision x(ad_array_size)
      double precision xd(ad_array_size)

      integer sim_iicroot
      
      common /obs_data/ obs

      common /state_info/ sim_iicroot
