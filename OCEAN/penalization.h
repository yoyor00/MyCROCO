      real pena_u(GLOBAL_2D_ARRAY,N)
      real pena_v(GLOBAL_2D_ARRAY,N)
      real pena_u2(GLOBAL_2D_ARRAY,N)
      real pena_v2(GLOBAL_2D_ARRAY,N)
      real pena_r(GLOBAL_2D_ARRAY,N)
      real h_pena_target(GLOBAL_2D_ARRAY)
      real Hzf(GLOBAL_2D_ARRAY,N)
      real rdrg_target
      integer last_index(GLOBAL_2D_ARRAY)
      integer last_index_w(GLOBAL_2D_ARRAY)
#if defined POROSITY
      real z_min_poro(GLOBAL_2D_ARRAY)
      real z_max_poro(GLOBAL_2D_ARRAY)
      real z_poro_discrete(GLOBAL_2D_ARRAY,0:N_poro_discrete)
      real z_int_poro(GLOBAL_2D_ARRAY,0:N_poro_discrete)
      real poro_discrete(GLOBAL_2D_ARRAY,1:N_poro_discrete)
#endif
      real pororatio(GLOBAL_2D_ARRAY,N)
      real epsilon_pena
#if defined POROSITY
      common/pena/pena_u, pena_v, pena_r,
     &   h_pena_target,Hzf,epsilon_pena,
     &  z_poro_discrete, poro_discrete,
     &   rdrg_target,z_min_poro,z_max_poro,z_int_poro,
     &   pororatio,
     & last_index, last_index_w
#else
      common/pena/pena_u,pena_v,pena_r,
     &   h_pena_target,Hzf,epsilon_pena,     
     &   rdrg_target,
     &   last_index, last_index_w
#endif
