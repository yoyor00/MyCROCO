C     -*- fortran -*-

C     size of the optimization problem
      integer ad_array_size
      parameter (ad_array_size=(lm+1+padd_x)*(mm+1+padd_e)*nnodes)

c     real size of the problem per node (<= ad_array_size/nnodes)
      integer ad_array_node_size

C     number of steps between cost function computations
      integer ad_ns
#if defined AD_DL_Z0B_CTRL
      parameter (ad_ns = 10)
#elif defined AD_ATLN
      parameter (ad_ns = 180)
#endif

C     number of cost function computations
      integer ad_nt
c     parameter (ad_nt = 2400)
#if defined AD_DL_Z0B_CTRL
      parameter (ad_nt = 240)
#elif defined AD_ATLN
      parameter (ad_nt = 48)
#endif

C     start of assimilation in the obs file
      integer ad_ast
#if defined AD_DL_Z0B_CTRL
      parameter (ad_ast = 1200)
#elif defined AD_ATLN
      parameter (ad_ast = 388)
#endif

C     number of time steps in the main file before assimilation
      integer ad_main_st
      parameter (ad_main_st = 1)

c     observations
      double precision ad_obs(GLOBAL_2D_ARRAY, ad_nt*ad_ns+3)
      double precision ad_obs_time(ad_nt*ad_ns+3)

c     state vector / process
      double precision ad_x(ad_array_size/nnodes)
      double precision ad_dz(ad_array_size/nnodes)

c     gradient vector / process
      double precision ad_g(ad_array_size/nnodes)

c     full control vector
      double precision ad_x_f(ad_array_size)

c     full gradient vector
      double precision ad_g_f(ad_array_size)

c     sum of all full gradient vectors
      double precision ad_sg_f(ad_array_size)

C     tangential vector
      double precision ad_xd(ad_array_size)

C     time step of the main simulation
      integer ad_sim_iicroot

C     general iteration counter
      integer ad_counter

C     cost function counter
      integer ad_cost_counter

c     timings
      double precision ad_dir_time
      double precision ad_adj_time

c     tidal period (M2)
      double precision TM2
      parameter(TM2 = 12.4206012*3600)

c     backup
      real ad_ubar_bck(GLOBAL_2D_ARRAY,4)
      real ad_vbar_bck(GLOBAL_2D_ARRAY,4)
      real ad_zeta_bck(GLOBAL_2D_ARRAY,4)
      real zob_bck(GLOBAL_2D_ARRAY)
      real Zobt_bck
      real ad_cost

c     rms
      real ad_rms,ad_irms,ad_irms_f
      integer ad_array_node_size_f
      integer ad_ta

      real ad_spval
      parameter(ad_spval = -999)

      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

C     commons
      common /ad_backup/ ad_ubar_bck, ad_vbar_bck, ad_zeta_bck,
     &     zob_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, Zobt_bck

      common /ad/ ad_array_node_size, ad_dz

      common /ad_timings/ ad_dir_time,ad_adj_time

      common /ad_obs_data/ ad_obs, ad_obs_time
      common /ad_state_info/ ad_sim_iicroot,ad_counter,ad_cost_counter,
     &     ad_cost,ad_rms,ad_irms,ad_irms_f,ad_ta
