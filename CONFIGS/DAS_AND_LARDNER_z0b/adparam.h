C     -*- fortran -*-
      
C     size of the optimization problem
      integer ad_array_size
c      parameter (ad_array_size=900)
      parameter (ad_array_size=(lm+1+padd_x)*(mm+1+padd_e)*nnodes)

c     real size of the problem per node (<= ad_array_size/nnodes)
      integer ad_array_node_size

c     number of proposed measure points
      integer npcpoints
      parameter (npcpoints=10)

C     number of steps between cost function computations
      integer ad_ns
      parameter (ad_ns = 10)

C     number of cost function computations
      integer ad_nt
c     parameter (ad_nt = 2400)
      parameter (ad_nt = 240)

C     start of assimilation in the obs file
      integer ad_ast
      parameter (ad_ast = 1200)

C     number of time steps in the main file before assimilation
      integer ad_main_st
      parameter (ad_main_st = 1)

c     observations
      double precision ad_obs(GLOBAL_2D_ARRAY, ad_nt*ad_ns+3)

c     state vector / process
      double precision ad_x(ad_array_size/nnodes)

c     gradient vector / process
      double precision ad_g(ad_array_size/nnodes)

c     full control vector
      double precision ad_x_f(ad_array_size)

c     full gradient vector
      double precision ad_g_f(ad_array_size)

c     sum of all full gradient vectors
      double precision ad_sg_f(ad_array_size)

c     number of validated measure points (<= ad_array_size/nnodes)
      integer ncpoints

c     MPI node of validated measure points (-1 for unvalidated points)
      integer ad_cpoint_node(npcpoints)

c     id of validated measure points
      integer ad_cpoint_id(npcpoints)

c     numbers validated measure points per nodes
      integer ncpoints_f(nnodes)

c     coordinates of measure points
      integer ad_i(npcpoints)
      integer ad_j(npcpoints)

c     latitudes/longitudes of control points on whole grid
      double precision ad_latr_f(npcpoints)
      double precision ad_lonr_f(npcpoints)
c     depth of control points
      double precision ad_h(npcpoints)
      double precision ad_h_f(npcpoints)

c     weighted coefficients
      double precision W(GLOBAL_2D_ARRAY,npcpoints)
      double precision SkW(GLOBAL_2D_ARRAY)

c     direct affectation of collocation points
      integer ad_node_colloc(GLOBAL_2D_ARRAY)
      integer ad_colloc(GLOBAL_2D_ARRAY)

C     tangential vector
      double precision ad_xd(ad_array_size)

C     time step of the main simulation
      integer sim_iicroot

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
      real ubar_bck(GLOBAL_2D_ARRAY,4)
      real vbar_bck(GLOBAL_2D_ARRAY,4)
      real zeta_bck(GLOBAL_2D_ARRAY,4)
      real zob_bck(GLOBAL_2D_ARRAY)
      real Zobt_bck
      real ad_cost

c     rms
      real ad_rms
      integer ad_ta

      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

C     commons
      common /ad_backup/ ubar_bck, vbar_bck, zeta_bck, zob_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, Zobt_bck

      common /ad/ ad_array_node_size

      common /ad_timings/ ad_dir_time,ad_adj_time

      common /ad_colloc_id/ ad_colloc,ad_node_colloc
      common /ad_collocation_coords/ ncpoints,ncpoints_f,ad_i,ad_j,
     &     ad_latr_f,ad_lonr_f,ad_h_f,ad_cpoint_node,ad_cpoint_id
      common /ad_weighted_coefs/ W,SkW
      common /ad_obs_data/ ad_obs
      common /ad_state_info/ sim_iicroot,ad_counter,ad_cost_counter,
     &     ad_cost,ad_rms,ad_ta
