! $Id:$
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!! Needed for common2device.py
# define NSLP1N -N_sl+1:N
# define NSLN -N_sl:N

# ifdef K3FAST
! !------------------------------------------------------------------
! ! nbq.h: Logical
! !------------------------------------------------------------------
      logical M2bc_nbq_flag
      common /nbq_M2bc/ M2bc_nbq_flag
! !------------------------------------------------------------------
! ! nbq.h: Integer scalars
! !------------------------------------------------------------------
      integer iteration_nbq      
      common /nbq_var2/ iteration_nbq      
      integer ifl_nbq     
      common /nbq_var3/ ifl_nbq     
      integer slip_nbq   
      common /nbq_var4/ slip_nbq    
#  ifdef KNHINT_3M
      integer nsdtnbq 
      common /nbq_nsdtnbq/ nsdtnbq
#  endif
      integer ndtnbq
      common /time_nbq1/ ndtnbq
      integer inc_faststep
      common/nbq_inc_faststep/inc_faststep
      integer nb_faststep
      common/nbq_nb_faststep/nb_faststep
! !------------------------------------------------------------------
! ! nbq.h: Real scalars
! !------------------------------------------------------------------
#  ifndef K3FAST_CSVISC2K
      real soundspeed_nbq
      common /nbq_param1/ soundspeed_nbq         
      real soundspeed2_nbq
      common /nbq_param2/ soundspeed2_nbq
#  else
      real soundspeed_nbq(GLOBAL_2D_ARRAY,NSLP1N) 
      common /nbq_param1/ soundspeed_nbq         
      real soundspeed2_nbq(GLOBAL_2D_ARRAY,NSLP1N) 
      common /nbq_param2/ soundspeed2_nbq
!$acc declare create( soundspeed_nbq, soundspeed2_nbq )      
#  endif
      double precision time_nbq               
      common /nbq_param3/ time_nbq
      double precision csvisc1_nbq                
      common /nbq_param4/ csvisc1_nbq 
      double precision csvisc2_nbq            
      common /nbq_param5/ csvisc2_nbq 
      double precision cw_int_nbq       
      common /nbq_param6/ cw_int_nbq
      double precision ifl_imp_nbq
      common /nbq_param7/ ifl_imp_nbq
      real dtnbq
      common /time_NBQ2/ dtnbq 
      real csound_nbq
      common /nbq_csound/ csound_nbq  
#  ifdef K3FAST_SACOUS
      real  period_exp  
      common/ACOUS1/period_exp
      real  for_a_exp   
      common/ACOUS2/for_a_exp
      real  dg_exp     
      common/ACOUS3/dg_exp 
      real  hmax_exp    
      common/ACOUS4/hmax_exp
      real  amp_exp
      common/ACOUS5/amp_exp
      real  rx_exp,ry_exp,rz_exp
      common/ACOUS6/rx_exp,ry_exp,rz_exp
#  endif
#  ifdef K3FAST_SEDLAYERS
      real rho_sdl
      common /SDL_RHO/ rho_sdl
#  endif
	real CFL_nbq
	common/nbq_CFL/CFL_nbq
! !------------------------------------------------------------------
! ! nbq.h: nbq.h: Real scalars or arrays
! !------------------------------------------------------------------
#  ifndef K3FAST_CSVISC2K
      real visc2_nbq
      common /nbq_visc2/ visc2_nbq
      real visc2v_nbq
      common /test_visc2v/ visc2v_nbq
#  else
      real visc2_nbq (GLOBAL_2D_ARRAY,NSLP1N)
      common /nbq_visc2/ visc2_nbq
      real visc2v_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common /test_visc2v/ visc2v_nbq
!$acc declare create( visc2_nbq , visc2v_nbq )      
#  endif
#  ifdef NBQ_SPONGE
      real visc2_nbq_sponge (GLOBAL_2D_ARRAY)
      common /nbq_visc2_sponge/ visc2_nbq_sponge
#  endif
      real visc2read_nbq
      common /nbq_visc2read/ visc2read_nbq
      real dtgrid_nbq
      common /nbq_dtgrid/ dtgrid_nbq
#  ifdef KNHINT_CORR
      real alphaw_nbq 
      common/KNHINT_1/alphaw_nbq
      integer alphaNw_nbq  
      common/KNHINT_2/alphaNw_nbq
#  endif
! !-------------------------------------------------------------------
! ! nbq.h: Real arrays 2D
! !-------------------------------------------------------------------
#  if defined NBQ_FREESLIP && ! defined K3FAST_SEDLAYERS
      real qdmw0_nbq(GLOBAL_2D_ARRAY)
      common /nbq_qdmw0_nbq/ qdmw0_nbq
#  endif
#  if defined K3FAST_ZETAW || defined PSOURCE || defined WET_DRY
      real DU_nbq(GLOBAL_2D_ARRAY)
      common /nbq_DU_nbq/ DU_nbq
      real DV_nbq(GLOBAL_2D_ARRAY)
      common /nbq_DV_nbq/ DV_nbq
#  endif
#  if defined K3FAST_ZETAW || defined PSOURCE || defined WET_DRY || defined KNHINT_ZETAW
      real dzeta_nbq(GLOBAL_2D_ARRAY)
      common /nbq_dzeta/dzeta_nbq
      real wsurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_wsurf/ wsurf_nbq
      real usurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_usurf/ usurf_nbq
      real vsurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_vsurf/ vsurf_nbq 
#  endif
#  ifdef NBQ_NUDGING
      real NBQnudgcof(GLOBAL_2D_ARRAY)
      common /nbq_nudg/ NBQnudgcof
#  endif
! !-------------------------------------------------------------------
! ! nbq.h: Real arrays 2D/3D
! !-------------------------------------------------------------------
#  ifdef K3FAST_UV
#   ifdef NBQ_GRID_SLOW
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY,NSLN,2)
      common /nbq_nods3/ dthetadiv_nbqdz
#   else
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY)
      common /nbq_nods3/ dthetadiv_nbqdz
#   endif /* NBQ_GRID_SLOW */
#  endif
! !-------------------------------------------------------------------
! ! nbq.h: Real arrays 3D+
! !-------------------------------------------------------------------
      real qdmu_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common /nbq_qdmu_nbq/ qdmu_nbq
      real qdmv_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common /nbq_qdmv_nbq/ qdmv_nbq
      real qdmw_nbq(GLOBAL_2D_ARRAY,NSLN)
      common /nbq_qdmw_nbq/ qdmw_nbq
#     ifdef NBQ_THETAIMP
       real qdmwold_nbq(GLOBAL_2D_ARRAY,NSLN)
       common /nbq_qdmwold_nbq/ qdmwold_nbq
#     endif
#  ifdef K3FAST_UV
!$acc declare create( dthetadiv_nbqdz )      
      real dZdxq_w(GLOBAL_2D_ARRAY,NSLN+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,NSLN+1)
      common /nbq_nods7/ dZdyq_w
!$acc declare create( dZdxq_w, dZdyq_w )
#  endif 
      real thetadiv_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common /nbq_thetadiv_nbq/ thetadiv_nbq
#  if defined NBQ_HZ_PROGNOSTIC 
      real thetadiv2_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common /nbq_thetadiv2_nbq/ thetadiv2_nbq
#  endif  
#  if defined K3FAST_C3D_UVSF &&  defined K3FAST_COUPLING3D
      real ru_int2d_nbq_bak(GLOBAL_2D_ARRAY,2)
      common /nbq_ru2d/ru_int2d_nbq_bak
      real rv_int2d_nbq_bak(GLOBAL_2D_ARRAY,2)
      common /nbq_rv2d/rv_int2d_nbq_bak
!$acc declare create( ru_int2d_nbq_bak, rv_int2d_nbq_bak)      
#  endif
      real ru_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ruint/ ru_int_nbq
      real rv_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rvint/ rv_int_nbq
#  ifdef K3FAST_COUPLING_SCH1
      real ru_intt_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ruintt/ ru_intt_nbq
      real rv_intt_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rvintt/ rv_intt_nbq
#  elif defined K3FAST_COUPLING_SCH2
      real ru_intt_nbq(GLOBAL_2D_ARRAY,N,2)
      common /nbq_ruintt/ ru_intt_nbq
      real rv_intt_nbq(GLOBAL_2D_ARRAY,N,2)
      common /nbq_rvintt/ rv_intt_nbq
#  endif
!$acc declare create( ru_int_nbq, rv_int_nbq )     
      real ru_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ru/ ru_nbq
      real rv_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rv/ rv_nbq
      real ru_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_runbq/ ru_nbq_avg2
      real rv_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_rvnbq/ rv_nbq_avg2
      real Hzw_nbq(GLOBAL_2D_ARRAY,NSLN)
      common /grid_Hzw_nbq/ Hzw_nbq
      real Hzu_nbq_inv(GLOBAL_2D_ARRAY,NSLN)
      common /grid_Hzu_nbq/ Hzu_nbq_inv
      real Hzv_nbq_inv(GLOBAL_2D_ARRAY,NSLN)
      common /grid_Hzv_nbq/ Hzv_nbq_inv
!$acc declare create(  Hzw_nbq, Hzu_nbq_inv, Hzv_nbq_inv )      
      real rw_int_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rwint/ rw_int_nbq   
#  ifdef K3FAST_COUPLINGW_SCH1
      real rw_intt_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rwintt/ rw_intt_nbq   
#  elif defined K3FAST_COUPLINGW_SCH2
      real rw_intt_nbq(GLOBAL_2D_ARRAY,0:N,2)
      common /nbq_rwintt/ rw_intt_nbq   
#  endif
      real rw_nbq(GLOBAL_2D_ARRAY,NSLN)
      common /nbq_rw/ rw_nbq
!$acc declare create( rw_int_nbq )      
      real rw_nbq_avg2(GLOBAL_2D_ARRAY,0:N)
      common /avg2_rwnbq/ rw_nbq_avg2
      real rho_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common/nbq_rho_nbq/rho_nbq
#  ifdef K3FAST_DIAGACOUS
      real p_nbq(GLOBAL_2D_ARRAY,NSLP1N)
      common/p_TL_nbq/p_nbq
      real p_nbq_max(GLOBAL_2D_ARRAY,NSLP1N)
      common/p_TL_nbq/p_nbq_max
#  endif
#  ifdef ONLINE_ANALYSIS
      real mvoa1(GLOBAL_2D_ARRAY,1:N)
      common/p_TL_nbq/mvoa1
      real mvoa2(GLOBAL_2D_ARRAY,1:N)
      common/p_TL_nbq/mvoa2
#  endif
#  ifdef NBQ_GRAV
      real rho_nh(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_nh/rho_nh
!$acc declare create( rho_nh )
#   ifdef CONVECT
      real rhoi_nh(GLOBAL_2D_ARRAY,N)
      common/nbq_rhoi_nh/rhoi_nh
!$acc declare create( rhoi_nh )      
#   endif
#  endif
      real rho_grd(GLOBAL_2D_ARRAY,NSLP1N)
      common/nbq_rho_grd/rho_grd
#  ifdef NBQ_MASS
      real rho_nbq_avg1(GLOBAL_2D_ARRAY,0:N)
      common /avg1_rhonbq/ rho_nbq_avg1
      real rhobar_nbq(GLOBAL_2D_ARRAY,4)
      common /nbq_rhobar/ rhobar_nbq
      real rhobar_nbq_avg1(GLOBAL_2D_ARRAY)
      common /nbq_rhobar_AVG1/ rhobar_nbq_avg1
#  endif
      real zw_nbq(GLOBAL_2D_ARRAY,NSLN)
      common /nbq_zw/ zw_nbq
# ifdef NBQ_HZCORRECT
       real Hz_correct(GLOBAL_2D_ARRAY,NSLP1N)
       common /grid_Hz_correct/ Hz_correct
!$acc declare create( Hz_correct )       
#  ifdef NBQ_HZCORR_DEBUG
      real  Hz_corr(GLOBAL_2D_ARRAY,N) 
      common/corr_Hz/Hz_corr
#  endif
# endif
#  ifdef NBQ_HZ_PROGNOSTIC
#   ifndef K3FAST_SEDLAYERS
      real Hz_bak2(GLOBAL_2D_ARRAY,1:N)
#   else
      real Hz_bak2(GLOBAL_2D_ARRAY,NSLP1N)
#   endif
      common /nbq_H_bak2/ Hz_bak2
#  endif
      real FC3D(GLOBAL_2D_ARRAY,NSLN+1)
      common /dum_FC3D/FC3D
      real DC3D(GLOBAL_2D_ARRAY,NSLN)
      common /dum_DC3D/DC3D
      real CF3D(GLOBAL_2D_ARRAY,NSLN)
      common /dum_CF3D/CF3D
!$acc declare create( FC3D, DC3D, CF3D )
#  ifdef ANA_MVB
      real rhoi_nbq(GLOBAL_2D_ARRAY,N)
      common/nbq_rhoi/rhoi_nbq
#  endif
#  ifdef KNHINT_WH
      real wzh_nbq(GLOBAL_2D_ARRAY,NSLN)
      common /nbq_wh/ wzh_nbq
#  endif

# ifdef XIOS
      real worknbq(GLOBAL_2D_ARRAY,0:N+N_sl)
      common /work3dnbq/ worknbq
# endif

! !-------------------------------------------------------------------
! ! nbq.h: Real BC
! !-------------------------------------------------------------------
#  if defined OBC_NBQ 
#   ifdef OBC_COM_WEST
      real qdmu_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_unbq_west/ qdmu_nbq_west
      real qdmv_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_vnbq_west/ qdmv_nbq_west
!$acc declare create( qdmu_nbq_west, qdmv_nbq_west )
#    ifdef K3FAST_W
      real qdmw_nbq_west(GLOBAL_1D_ARRAYETA,0:N,2)
      common /bry_wnbq_west/ qdmw_nbq_west
      real  rho_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_rnbq_west/ rho_nbq_west
!$acc declare create( qdmw_nbq_west, rho_nbq_west )
#    endif
#   endif
#   ifdef OBC_COM_EAST
      real qdmu_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_unbq_east/ qdmu_nbq_east
      real qdmv_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_vnbq_east/ qdmv_nbq_east
!$acc declare create( qdmu_nbq_east, qdmv_nbq_east )
#    ifdef K3FAST_W
      real qdmw_nbq_east(GLOBAL_1D_ARRAYETA,0:N,2)
      common /bry_wnbq_east/ qdmw_nbq_east
      real  rho_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_rnbq_east/ rho_nbq_east
!$acc declare create( qdmw_nbq_east, rho_nbq_east )
#    endif
#   endif
#   ifdef OBC_COM_SOUTH
      real qdmu_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_unbq_south/ qdmu_nbq_south
      real qdmv_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_vnbq_south/ qdmv_nbq_south
!$acc declare create( qdmu_nbq_south, qdmv_nbq_south )
#    ifdef K3FAST_W
      real qdmw_nbq_south(GLOBAL_1D_ARRAYXI,0:N,2)
      common /bry_wnbq_south/ qdmw_nbq_south
      real  rho_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_rnbq_south/ rho_nbq_south
!$acc declare create( qdmw_nbq_south, rho_nbq_south )
#    endif
#   endif
#   ifdef OBC_COM_NORTH
      real qdmu_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_unbq_north/ qdmu_nbq_north
      real qdmv_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_vnbq_north/ qdmv_nbq_north
!$acc declare create( qdmu_nbq_north, qdmv_nbq_north )
#    ifdef K3FAST_W
      real qdmw_nbq_north(GLOBAL_1D_ARRAYXI,0:N,2)
      common /bry_wnbq_north/ qdmw_nbq_north
      real  rho_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_rnbq_north/ rho_nbq_north
!$acc declare create( qdmw_nbq_north, rho_nbq_north )
#    endif
#   endif
#  endif    
# endif /* K3FAST */
