
#ifdef NBQ
                            
      
!**********************************************************************
      logical WEST_INTER_NBQ         
      common /nbq_log1/ WEST_INTER_NBQ 
      logical EAST_INTER_NBQ         
      common /nbq_log2/ EAST_INTER_NBQ                                            
      logical SOUTH_INTER_NBQ       
      common /nbq_log3/ SOUTH_INTER_NBQ                                               
      logical NORTH_INTER_NBQ       
      common /nbq_log4/ NORTH_INTER_NBQ        

!**********************************************************************
      integer   istr_nh    
      common /nbq_int1/ istr_nh                                              
      integer   jstr_nh    
      common /nbq_int2/ jstr_nh                                                     
      integer   iend_nh   
      common /nbq_int3/ iend_nh                                                       
      integer   jend_nh    
      common /nbq_int4/ jend_nh                                                       
      integer   istru_nh     
      common /nbq_int5/  istru_nh                                                   
      integer   jstru_nh    
      common /nbq_int6/   jstru_nh                                                   
      integer   istrv_nh   
      common /nbq_int7/ istrv_nh                                                      
      integer   jstrv_nh     
      common /nbq_int8/ jstrv_nh                                                    
      integer   iendu_nh      
      common /nbq_int9/ iendu_nh                                                    
      integer   jendu_nh      
      common /nbq_int10/ jendu_nh                                                    
      integer   iendv_nh    
      common /nbq_int11/ iendv_nh                                                       
      integer   jendv_nh      
      common /nbq_int12/ jendv_nh                                                      
      integer   istrq_nh    
      common /nbq_int13/ istrq_nh                                                      
      integer   iendq_nh     
      common /nbq_int14/iendq_nh                                                      
      integer   jstrq_nh     
      common /nbq_int15/ jstrq_nh                                                     
      integer   jendq_nh    
      common /nbq_int16/ jendq_nh  
                                   
!**********************************************************************              
      integer iteration_nbq_max  
      common /nbq_var1/ iteration_nbq_max             
      integer iteration_nbq      
      common /nbq_var2/ iteration_nbq      
      integer ifl_nbq     
      common /nbq_var3/ ifl_nbq     
      integer slip_nbq   
      common /nbq_var4/ slip_nbq    

!**********************************************************************
      real soundspeed_nbq(GLOBAL_2D_ARRAY) 
      common /nbq_param1/ soundspeed_nbq                              
      real soundspeed2_nbq(GLOBAL_2D_ARRAY) 
      common /nbq_param2/ soundspeed2_nbq                                                
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

!**********************************************************************
      integer ndtnbq
      common /time_nbq1/ ndtnbq
      real dtnbq
      common /time_nbq2/ dtnbq 
      real csound_nbq
      common /nbq_csound/ csound_nbq
      real visc2_nbq
      common /nbq_visc2/ visc2_nbq

      real dtgrid_nbq
      common /nbq_dtgrid/ dtgrid_nbq

!**********************************************************************
      real qdmu_nbq(GLOBAL_2D_ARRAY,N)
      common/nbq_qdmu_nbq/qdmu_nbq
      real qdmv_nbq(GLOBAL_2D_ARRAY,N)
      common/nbq_qdmv_nbq/qdmv_nbq
      real qdmw_nbq(GLOBAL_2D_ARRAY,0:N)
      common/nbq_qdmw_nbq/qdmw_nbq
       
      real thetadiv_nbq(GLOBAL_2D_ARRAY,0:N)
      common/nbq_thetadiv_nbq/thetadiv_nbq
      real thetadiv2_nbq(GLOBAL_2D_ARRAY,0:N)
      common/nbq_thetadiv2_nbq/thetadiv2_nbq
      real thetadiv3_nbq(GLOBAL_2D_ARRAY,0:N)
      common/nbq_thetadiv3_nbq/thetadiv3_nbq

      real ru_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ruint/ ru_int_nbq
      real rv_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rvint/ rv_int_nbq
      real rw_int_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rwint/ rw_int_nbq   

      real rho_nbq(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_nbq/rho_nbq

      real ru_nbq_ext(GLOBAL_2D_ARRAY,N)
      common /nbq_ru_ext/ ru_nbq_ext
      real rv_nbq_ext(GLOBAL_2D_ARRAY,N)
      common /nbq_rv_ext/ rv_nbq_ext
      real rw_nbq_ext(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rw_ext/ rw_nbq_ext

      real ru_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_runbq/ ru_nbq_avg2
      real rv_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_rvnbq/ rv_nbq_avg2
      real rw_nbq_avg2(GLOBAL_2D_ARRAY,0:N)
      common /avg2_rwnbq/ rw_nbq_avg2

       
!**********************************************************************
      real e(GLOBAL_2D_ARRAY)
      common /nbq_e/ e
      real eomn(GLOBAL_2D_ARRAY)
      common /nbq_eomn/ eomn
      real cosa(GLOBAL_2D_ARRAY)
      common /nbq_cosa/ cosa
      real sina(GLOBAL_2D_ARRAY)
      common /nbq_sina/ sina

!**********************************************************************
       real qdm_u_ext (GLOBAL_2D_ARRAY)  
       common /nbq_qdmu/qdm_u_ext
       real qdm_v_ext (GLOBAL_2D_ARRAY)  
       common /nbq_qdmv/qdm_v_ext 
      real rubar_nbq(GLOBAL_2D_ARRAY)
      common /nbq_rubar/ rubar_nbq
      real rvbar_nbq(GLOBAL_2D_ARRAY)
      common /nbq_rvbar/ rvbar_nbq

!**********************************************************************

!     real qdm_u2 (GLOBAL_2D_ARRAY,0:N+1)
!     common /nbq_u2/ qdm_u2  
!     real qdm_v2 (GLOBAL_2D_ARRAY,0:N+1)
!     common /nbq_v2/ qdm_v2  
!     real qdm_w2 (GLOBAL_2D_ARRAY,0:N+1) 
!     common /nbq_w2/ qdm_w2
      real Hzw_half_nbq(GLOBAL_2D_ARRAY,0:N)
      common /grid_Hzw_half_nbq/ Hzw_half_nbq
# ifndef NBQ_ZETAW
      real zr_half_nbq(GLOBAL_2D_ARRAY,N)
      common /grid_zr_half_nbq/ zr_half_nbq
      real zw_half_nbq(GLOBAL_2D_ARRAY,0:N)
      common /grid_zw_half_nbq/ zw_half_nbq     
# endif    

!**********************************************************************
# ifndef NBQ_IJK
      integer   ifl_nh    
      common /nbq_int16/ ifl_nh                                                         
      integer l_nh                                                      
      integer l1_nh                                                           
      integer l2_nh        
      integer l_nbq    
      integer l1_nbq  
      integer l2_nbq   
# endif
  
!**********************************************************************
# ifdef NBQ_MASS
       real Hzr_half_nbq(GLOBAL_2D_ARRAY,N)
       common /grid_Hzr_half_nbq/ Hzr_half_nbq
# else
#  define Hzr_half_nbq Hz_half
# endif

# ifdef NBQ_TRACERS
       real Hz_tra(GLOBAL_2D_ARRAY,N)
       common /grid_Hz_tra/ Hz_tra
#endif

  
!**********************************************************************
      real DU_nbq  (GLOBAL_2D_ARRAY)
      common /nbq_DU_nbq/ DU_nbq
      real DV_nbq  (GLOBAL_2D_ARRAY)
      common /nbq_DV_nbq/ DV_nbq
# ifdef NBQ_MASS
      real rho_nbq_ext (GLOBAL_2D_ARRAY,N)
      real rho_nbq_avg1(GLOBAL_2D_ARRAY,0:N)
      real rhos_nbq_int(GLOBAL_2D_ARRAY)
      common /nbq_rho_ext/ rho_nbq_ext
      common /avg1_rhonbq/ rho_nbq_avg1
      common /int_rhonbq/ rhos_nbq_int
      real rhobar_nbq(GLOBAL_2D_ARRAY,4)
      common /nbq_rhobar/ rhobar_nbq
      real rhobar_nbq_avg1(GLOBAL_2D_ARRAY)
      common /nbq_rhobar_AVG1/ rhobar_nbq_avg1
      real rhobar_nbq_int(GLOBAL_2D_ARRAY)
      common /nbq_rhobar_int/ rhobar_nbq_int
# endif
  
!**********************************************************************
# ifdef NBQ_NODS
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY,0:N,2)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods7/ dZdyq_w
# endif

!**********************************************************************
# ifdef NBQ_ZETAW
      real wmean_nbq(GLOBAL_2D_ARRAY,4)
      common /nbq_wmean/wmean_nbq
      real umean_nbq(GLOBAL_2D_ARRAY)
      common /nbq_umean/umean_nbq
      real vmean_nbq(GLOBAL_2D_ARRAY)
      common /nbq_vmean/vmean_nbq
      integer :: knew2,kstp2,kbak2,kold2
      common /gridext1/knew2,kstp2,kbak2,kold2   
      real cff4
      common /com_cff4/ cff4
      real cff5
      common /com_cff5/ cff5
      real cff6
      common /com_cff6/ cff6
      real cff7
      common /com_cff7/ cff7
      real cff8
      common /com_cff8/ cff8
      real cff9
      common /com_cff9/ cff9
      real cff10
      common /com_cff10/ cff10
      integer flag_grid
      common /grid_flag/ flag_grid
      integer IstrU2,JstrV2
      integer IstrR2,IendR2
      integer JstrR2,JendR2

#  if !defined NBQ_NODS
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods7/ dZdyq_w
#   endif

#  endif
#  ifdef NBQ_NUDGING
      real nudg_coef_nbq(GLOBAL_2D_ARRAY)
      common /nbq_nudg1/nudg_coef_nbq
#  endif

!**********************************************************************
# ifdef NBQ_DTDRHO
      real hrho_nbq(GLOBAL_2D_ARRAY,1:N,4)
      common /nbq_hrho/hrho_nbq
      real z_nbq(GLOBAL_2D_ARRAY,0:N,4)
      common /nbq_z/z_nbq
# endif

# ifdef NBQ_DTDRHO2
      real zr_nbq(GLOBAL_2D_ARRAY,N,4)
      common /nbq_zr/zr_nbq
      real z_nbq(GLOBAL_2D_ARRAY,0:N,4)
      common /nbq_z/z_nbq

#  if defined NBQ_ZETAW 
      real rho_bak(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_bak/rho_bak
#  endif

# endif

!# ifdef NBQ_DTDRHO2B
      real Hz_bak2(GLOBAL_2D_ARRAY,1:N)
      common /nbq_H_bak2/ Hz_bak2
!# endif

!**********************************************************************
# if defined ACOUSTIC && defined NBQ_IJK

      real  period_exp  
      common/ACOUS1/period_exp
      real  for_a_exp   
      common/ACOUS2/for_a_exp
      real  dg_exp     
      common/ACOUS3/dg_exp 
      real  hmax_exp    
      common/ACOUS4/hmax_exp
      real  amp_exp
      common/ACOUS4/amp_exp

#endif

      real  Hz_corr(GLOBAL_2D_ARRAY,N) 
      common/corr_Hz/Hz_corr

!      real Hz_t(GLOBAL_2D_ARRAY,1:N,4)
!      common /nbq_H_t/ Hz_t

#endif /* NBQ */

  
