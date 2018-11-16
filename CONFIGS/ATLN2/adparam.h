C     -*- fortran -*-

C     size of the problem (number of control variables)
      integer ad_array_size
      parameter (ad_array_size=(3+Lm+3+padd_X)*(3+Mm+3+padd_E))
c      parameter (ad_array_size=(lm+1+padd_x)*(mm+1+padd_e))

C     size of the assimilation window (number of steps)
C     on change check BINOMIAL-CKP param in cost_fun
      integer ad_nt
      parameter (ad_nt = 120)

C     start of assimilation in the obs file
      integer ad_ast
      parameter (ad_ast = 14401)

C     number of time steps in the main file before assimilation
      integer ad_main_st
      parameter (ad_main_st = 1)

C     observations
      double precision ad_obs(GLOBAL_2D_ARRAY, ad_nt+2)

C     control vector
      double precision ad_x(ad_array_size)

C     tangential vector
      double precision ad_xd(ad_array_size)

C     time step of the main simulation
      integer sim_iicroot

C     commons
      common /obs_data/ ad_obs
      common /state_info/ sim_iicroot
