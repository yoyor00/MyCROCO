C     -*- fortran -*-

C     size of the problem (number of control variables)
      integer ad_array_size
      parameter (ad_array_size=1)
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

c     coordinates of control points (i.e collocation points)
      integer ad_i(ad_array_size)
      integer ad_j(ad_array_size)

c     weighted coefficients
      double precision W(GLOBAL_2D_ARRAY,ad_array_size)
      double precision SkW(GLOBAL_2D_ARRAY)

c     direct affectation of collocation points
      integer ad_colloc(GLOBAL_2D_ARRAY)
      
C     tangential vector
      double precision ad_xd(ad_array_size)

C     time step of the main simulation
      integer sim_iicroot

c     tidal period (M2)
      double precision TM2
      parameter(TM2 = 12.4206012)
      
C     commons
      common /colloc_id/ ad_colloc
      common /collocation_coords/ ad_i,ad_j
      common /weighted_coefs/ W,SkW
      common /obs_data/ ad_obs
      common /state_info/ sim_iicroot
