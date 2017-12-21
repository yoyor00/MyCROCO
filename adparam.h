C     -*- fortran -*-
      integer ad_array_size
      parameter (ad_array_size=1 + (3+Lm+3+padd_X)*(3+Mm+3+padd_E))

      integer NT, NX, NY
      parameter (NT = 8641, NX = 17, NY = 17)

      real*4 obsz(Lm, Mm, NT)
      real*4 robsz(NX, NY, NT)
      real*4 obsu(NX, NY, NT)
      real*4 robsu(NX, NY, NT)
      real*4 obsv(NX, NY, NT)
      real*4 robsv(NX, NY, NT)
      
      double precision x(ad_array_size)

      integer sim_iicroot
      
      common /obs_data/ obsz, robsz, obsu, robsu, obsv, robsv

      common /state_info/ sim_iicroot
