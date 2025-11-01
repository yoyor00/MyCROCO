!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
! Include file "penalization.h".
! ==============================
!
      real pena_u(GLOBAL_2D_ARRAY,N)
      real pena_v(GLOBAL_2D_ARRAY,N)
      real pena_r(GLOBAL_2D_ARRAY,N)
      real h_pena_target(GLOBAL_2D_ARRAY)
      real Hzf(GLOBAL_2D_ARRAY,N)
      real z_w_tmp(GLOBAL_2D_ARRAY,0:N)
      real z_lastindex(GLOBAL_2D_ARRAY)
      integer last_index(GLOBAL_2D_ARRAY)
      integer last_index_w(GLOBAL_2D_ARRAY)
#ifdef POROSITY
      real z_min_poro(GLOBAL_2D_ARRAY)
      real z_max_poro(GLOBAL_2D_ARRAY)
      real z_poro_discrete(GLOBAL_2D_ARRAY,0:N_poro_discrete)
      real z_int_poro(GLOBAL_2D_ARRAY,0:N_poro_discrete)
      real poro_discrete(GLOBAL_2D_ARRAY,1:N_poro_discrete)
#endif
      real pororatio(GLOBAL_2D_ARRAY,N)
      real epsilon_pena
#ifdef POROSITY
      common/pena/pena_u, pena_v, pena_r,
     &            h_pena_target, Hzf, z_w_tmp, epsilon_pena,
     &            z_poro_discrete, poro_discrete, pororatio,
     &            z_min_poro, z_max_poro, z_int_poro,
     &            z_lastindex, last_index, last_index_w
#else
      common/pena/pena_u, pena_v, pena_r,
     &            h_pena_target, Hzf, z_w_tmp, epsilon_pena,     
     &            z_lastindex, last_index, last_index_w
#endif
