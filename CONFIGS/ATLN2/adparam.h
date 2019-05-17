C     -*- fortran -*-

C     size of the optimization problem
      integer ad_array_size
      parameter (ad_array_size=100*nnodes)
c      parameter (ad_array_size=(lm+1+padd_x)*(mm+1+padd_e))


C     number of steps between cost function computation
      integer ad_ns
      parameter (ad_ns = 180)

C     size of the assimilation window
C     on change check BINOMIAL-CKP param in cost_fun
      integer ad_nt
      parameter (ad_nt = 48)

C     start of assimilation in the obs file
      integer ad_ast
      parameter (ad_ast = 140)

C     number of time steps in the main file before assimilation
      integer ad_main_st
      parameter (ad_main_st = 1)

c     observations
      double precision ad_obs(GLOBAL_2D_ARRAY, ad_nt+2)

c     control vector / process
      double precision ad_x(ad_array_size)

c     gradient vector / process
      double precision ad_g(ad_array_size)
      
c     full control vector
      double precision ad_x_f(ad_array_size)

c     full gradient vector
      double precision ad_g_f(ad_array_size)

c     sum of all full gradient vectors
      double precision ad_sg_f(ad_array_size)

c     number of measure points (<= ad_array_size/nnodes)
      integer ncpoints

c     coordinates of control points
      integer ad_i(ad_array_size)
      integer ad_j(ad_array_size)

c     latitudes/longitudes of control points on whole grid
      double precision ad_latr(ad_array_size)
      double precision ad_lonr(ad_array_size)

      double precision ad_latr_f(ad_array_size)
      double precision ad_lonr_f(ad_array_size)
c     depth of control points
      double precision ad_h(ad_array_size)
      double precision ad_h_f(ad_array_size)

c     weighted coefficients
      double precision W(GLOBAL_2D_ARRAY,ad_array_size)
      double precision SkW(GLOBAL_2D_ARRAY)

c     direct affectation of collocation points
      integer ad_colloc(GLOBAL_2D_ARRAY)

C     tangential vector
      double precision ad_xd(ad_array_size)

C     time step of the main simulation
      integer sim_iicroot

C     general iteration counter
      integer ad_counter

C     cost function counter
      integer ad_cost_counter
      
c     tidal period (M2)
      double precision TM2
      parameter(TM2 = 12.4206012)

c     backup
      real ubar_bck(GLOBAL_2D_ARRAY,4)
      real vbar_bck(GLOBAL_2D_ARRAY,4)
      real zeta_bck(GLOBAL_2D_ARRAY,4)
      real zob_bck(GLOBAL_2D_ARRAY)
      real Zobt_bck
      real ad_cost
      
      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

C     commons
      common /backup/ ubar_bck, vbar_bck, zeta_bck, zob_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, Zobt_bck

      common /colloc_id/ ad_colloc
      common /collocation_coords/ ncpoints,ad_i,ad_j,ad_latr_f,ad_lonr_f
     &     ,ad_h_f
      common /weighted_coefs/ W,SkW
      common /obs_data/ ad_obs
      common /state_info/ sim_iicroot,ad_counter,ad_cost_counter,ad_cost
