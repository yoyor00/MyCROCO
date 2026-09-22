! !
! !===================================================================
! !
! ! Momentum equations: time stepping to time n+1/2
! !
! !===================================================================
! !
# undef DC
# undef FC
# undef CF
! !
! !-------------------------------------------------------------------
! ! U momentum
! !-------------------------------------------------------------------
! !
        do j=Jstr,Jend
        do i=IstrU,Iend
          WORK(i,j)=pm_u(i,j)*pn_u(i,j)
        enddo
        enddo 
! !        
! !*******************************************
! !************* U: 1st timestep ************* 
! !*******************************************
! !       
        if (FIRST_TIME_STEP) then
          do k=1,N
           do j=Jstr,Jend
            do i=IstrU,Iend    
! !
! !******** Predictor @ t+dt/2 ***************      
! !
# ifdef VADV_ADAPT_IMP
              cff = 1.
# else
              cff = 2./(Hz_half(i,j,k)+Hz_half(i-1,j,k))
# endif           
              cff0=WORK(i,j)*ru(i,j,k)
              u(i,j,k,nnew)=( u(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                       +cdt*cff0 )*cff
! !
! !******** RHS Slow to Fast *****************      
! !
# if defined M3FAST || (defined K3FAST_COUPLING_SCH0 && defined K3FAST_C3D_UVSF )
              ru_int_nbq(i,j,k) =cff0
# elif defined K3FAST_C3D_UVSF 
              ru_int_nbq(i,j,k) =cff0
#  ifndef K3FAST_COUPLING_SCH2
              ru_intt_nbq(i,j,k,nstp)=cff0
#  else
              ru_intt_nbq(i,j,k,nstp)=cff0
#  endif
# endif
! !
! !******** Depth-average ********************      
! !
# if defined M3FAST || (defined K3FAST_COUPLING2D && defined K3FAST_C3D_UVSF)
              ru_int_nbq_2d(i,j)= ru_int_nbq_2d(i,j) 
     &                           +ru_int_nbq(i,j,k) 
# endif
! !
! !******** u(t) x dz ************************      
! !
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5*( Hz(i  ,j,k)
     &                                         +Hz(i-1,j,k))
            enddo
          enddo
         enddo
! !*******************************************
! !************* U: timesteps > 1 ************ 
! !*******************************************
        else
          do k=1,N
            do j=Jstr,Jend
             do i=IstrU,Iend
! !
! !******** Predictor @ t+dt/2 ***************      
! !
# ifdef VADV_ADAPT_IMP
              cff = 1.
# else
              cff = 2./(Hz_half(i,j,k)+Hz_half(i-1,j,k))
# endif 
              qdm_nstp = u(i,j,k,nstp)*0.5*( Hz(i  ,j,k)
     &                                      +Hz(i-1,j,k))
              qdm_indx = u(i,j,k,indx)*0.5*( Hz_bak(i  ,j,k)
     &                                      +Hz_bak(i-1,j,k))
              u(i,j,k,nnew)=( cff14*qdm_nstp
     &                       +cff15*qdm_indx
     &                       +cdt*WORK(i,j)*(ru(i,j,k)
# if defined M3FAST || defined K3FAST_C3D_UVFS
     &                                     +ru_nbq_avg1(i,j,k)
# endif     
     &                                      ))*cff 
! !
! !******** RHS Slow to Fast *****************      
! !
# if defined M3FAST  || defined K3FAST_COUPLING_SCH0
              ru_int_nbq(i,j,k)=((qdm_nstp-qdm_indx) /dt
     &               -WORK(i,j)*ru_nbq_avg2(i,j,k))
#  ifdef MASKING
     &                          *umask(i,j)
#  endif   
# elif defined K3FAST_C3D_UVSF
#  ifndef K3FAST_COUPLING_SCH2
              cff0=((qdm_nstp-qdm_indx) /dt
     &               -WORK(i,j)*ru_nbq_avg2(i,j,k))
              ru_int_nbq(i,j,k)=(cff6*cff0
     &                          +cff7*ru_intt_nbq(i,j,k,3-nstp)
     &                          +cff8*ru_intt_nbq(i,j,k,nstp))
#   ifdef MASKING
     &                          *umask(i,j)
#   endif        
              ru_intt_nbq(i,j,k,nstp)=cff0
#  else
           cff0=WORK(i,j)*ru(i,j,k)! ! *(1.-gamma)   
           ru_int_nbq(i,j,k)=(cff1 *cff0
     &                      + cff2 *ru_intt_nbq(i,j,k,3-nstp)
     &                      + cff13*ru_intt_nbq(i,j,k,nstp))  
#   ifdef MASKING
     &                         *umask(i,j)
#   endif        
           ru_intt_nbq(i,j,k,nstp)=cff0
#  endif /* ! !K3FAST_COUPLING_SCH2 */
! !
! !******** Test dissipation *****************      
! !                  
#  ifdef K3FAST_DISSCOUPLING
              cff10=((qdm_nstp-qdm_indx) /dt
     &              -WORK(i,j)*ru_nbq_avg2(i,j,k))
              ru_int_nbq(i,j,k) = (1.-alpha_uv) * ru_int_nbq(i,j,k) ! !ICI
     &                            +alpha_uv * cff10
#  endif     
# endif /* M3FAST || K3FAST_C3D_UVSF */
! !
! !******** Depth-average ********************      
! !
# if defined M3FAST || (defined K3FAST_COUPLING2D && defined K3FAST_C3D_UVSF)
              ru_int_nbq_2d(i,j)= ru_int_nbq_2d(i,j) 
     &                           +ru_int_nbq(i,j,k) 
# endif
! !
! !******** u(t) x dz ************************      
! !
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            enddo
          enddo  
         enddo
        endif ! ! .not. FIRST_TIME_STEP
! !
! !******** Substract depth average **********
! !    
# if defined M3FAST || (defined K3FAST_COUPLING2D && defined K3FAST_C3D_UVSF)
          do k=1,N
            do j=Jstr,Jend
             do i=IstrU,Iend
               ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)
     &            - ru_int_nbq_2d(i,j) *(Hz(i,j,k)+Hz(i-1,j,k))
     &                               /( Zt_avg1(i  ,j)+h(i  ,j)
     &                                 +Zt_avg1(i-1,j)+h(i-1,j))
            enddo
          enddo
        enddo
# endif
! !
! !-------------------------------------------------------------------
! ! V momentum
! !-------------------------------------------------------------------
! !
        do j=JstrV,Jend
          do i=Istr,Iend
            WORK(i,j)=pm_v(i,j)*pn_v(i,j)
          enddo
        enddo
! !
! !*******************************************
! !************* V: 1st timestep ************* 
! !*******************************************
! !
          if (FIRST_TIME_STEP) then
            do k=1,N
             do j=JstrV,Jend
              do i=Istr,Iend
! !
! !******** Predictor @ t+dt/2 ***************      
! !
# ifdef VADV_ADAPT_IMP
                cff = 1.
# else
                cff = 2./(Hz_half(i,j,k)+Hz_half(i,j-1,k))
# endif 
                cff0=WORK(i,j)*rv(i,j,k)
                v(i,j,k,nnew)=( v(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                         +cdt*cff0)*cff
! !
! !******** RHS Slow to Fast *****************      
! !
# if defined M3FAST || (defined K3FAST_COUPLING_SCH0 && defined K3FAST_C3D_UVSF)
                rv_int_nbq(i,j,k) =cff0
# elif defined K3FAST_C3D_UVSF
                rv_int_nbq(i,j,k) =cff0
#  ifndef K3FAST_COUPLING_SCH2
                rv_intt_nbq(i,j,k,nstp)=cff0
#  else
                rv_intt_nbq(i,j,k,nstp)=cff0
#  endif
# endif
! !
! !******** Depth-average ********************      
! !
# if defined M3FAST || (defined K3FAST_COUPLING2D && defined K3FAST_C3D_UVSF)
                rv_int_nbq_2d(i,j)= rv_int_nbq_2d(i,j)
     &                             +rv_int_nbq(i,j,k) 
# endif
! !
! !******** v(t) x dz ************************      
! !
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5*( Hz(i,j  ,k)
     &                                           +Hz(i,j-1,k))
              enddo
            enddo
          enddo
! !
! !*******************************************
! !************* V: timesteps > 1 ************ 
! !*******************************************
! !
          else
            do k=1,N
              do j=JstrV,Jend
               do i=Istr,Iend     
! !
! !******** Predictor @ t+dt/2 ***************      
! !
# ifdef VADV_ADAPT_IMP              
                cff = 1.
# else
                cff = 2./(Hz_half(i,j,k)+Hz_half(i,j-1,k))
# endif              
                qdm_nstp = v(i,j,k,nstp)*0.5*( Hz(i,j  ,k)
     &                                        +Hz(i,j-1,k))
                qdm_indx = v(i,j,k,indx)*0.5*( Hz_bak(i  ,j,k)
     &                                        +Hz_bak(i,j-1,k))
                v(i,j,k,nnew)=( cff14*qdm_nstp
     &                         +cff15*qdm_indx
     &                         +cdt*WORK(i,j)*(rv(i,j,k)
# if defined M3FAST || defined K3FAST_C3D_UVFS
     &                                  +rv_nbq_avg1(i,j,k)
# endif     
     &                                  ))*cff
! !
! !******** RHS Slow to Fast *****************      
! !
# if defined M3FAST || defined K3FAST_COUPLING_SCH0
                rv_int_nbq(i,j,k)=((qdm_nstp-qdm_indx) /dt
     &               -WORK(i,j)*rv_nbq_avg2(i,j,k))
#   ifdef MASKING
     &                           *vmask(i,j)   
#   endif   
# elif defined K3FAST_C3D_UVSF
#  ifndef K3FAST_COUPLING_SCH2
                cff0=((qdm_nstp-qdm_indx) /dt
     &               -WORK(i,j)*rv_nbq_avg2(i,j,k))
                rv_int_nbq(i,j,k)=(cff6*cff0
     &                            +cff7*rv_intt_nbq(i,j,k,3-nstp)
     &                            +cff8*rv_intt_nbq(i,j,k,nstp))
#   ifdef MASKING
     &                            *vmask(i,j)   
#   endif     
                rv_intt_nbq(i,j,k,nstp)=cff0
#  else
                cff0=WORK(i,j)*rv(i,j,k)
                rv_int_nbq(i,j,k)= (cff1*cff0
     &                            + cff2 *rv_intt_nbq(i,j,k,3-nstp)
     &                            + cff13*rv_intt_nbq(i,j,k,nstp))
#   ifdef MASKING
     &                           *vmask(i,j)   
#   endif     
                rv_intt_nbq(i,j,k,nstp)=cff0
#  endif /* ! !K3FAST_COUPLING_SCH2 */
! !
! !******** Test dissipation *****************      
! !    
#  ifdef K3FAST_DISSCOUPLING
                cff10=((qdm_nstp-qdm_indx) /dt
     &                -WORK(i,j)*rv_nbq_avg2(i,j,k))
                rv_int_nbq(i,j,k) = (1.-alpha_uv) * rv_int_nbq(i,j,k)  
     &                             +alpha_uv * cff10
#  endif   
# endif /* M3FAST || K3FAST_C3D_UVSF */ 
! !
! !******** Depth-average ********************      
! !  
# if defined M3FAST || (defined K3FAST_COUPLING2D && defined K3FAST_C3D_UVSF)
              rv_int_nbq_2d(i,j)= rv_int_nbq_2d(i,j)
     &                           +rv_int_nbq(i,j,k) 
# endif
! !
! !******** v(t) x dz ************************      
! !
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              enddo
            enddo
           enddo    
          endif        ! !   .not. FIRST_TIME_STEP 
! !
! !******** Substract depth average **********
! !    
# if defined M3FAST || (defined K3FAST_COUPLING2D && defined K3FAST_C3D_UVSF)
            do k=1,N
             do j=JstrV,Jend
              do i=Istr,Iend
               rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)
     &            - rv_int_nbq_2d(i,j) *(Hz(i,j,k)+Hz(i,j-1,k))
     &                               /( Zt_avg1(i,j  )+h(i,j  )
     &                                 +Zt_avg1(i,j-1)+h(i,j-1))
            enddo
          enddo
        enddo
# endif   
! !
! !-------------------------------------------------------------------
! ! Solve tridiag system for U & V (if required)
! !-------------------------------------------------------------------
! !
# ifdef VADV_ADAPT_IMP
!$acc loop private(DC)
      do j=Jstr,Jend
        do i=IstrU,Iend
          DC(i,0)=pm_u(i,j)*pn_u(i,j)
        enddo
#  undef  TRIDIAG_TRA
#  define TRIDIAG_U
#  undef  TRIDIAG_V
#  include "tridiag_pred.h"
        if (j.ge.JstrV) then
          do i=Istr,Iend
            DC(i,0)=pm_v(i,j)*pn_v(i,j)
          enddo
#  undef  TRIDIAG_TRA
#  undef  TRIDIAG_U
#  define TRIDIAG_V
#  include "tridiag_pred.h"
        endif
        
      enddo     ! !<-- j
# endif
!$acc end kernels      
! !
! !-------------------------------------------------------------------
! ! W momentum
! !-------------------------------------------------------------------
! !
# if defined NBQ || defined K3SLOW_W
! !
! !----------------------------------------------------------
! !<-- at this point
! !<-- rw            contains internal 3D advection + Coriolis
! !<-- rw_nbq_avg1   contains NBQ pressure gradient (+ gravity)
! !                                        + second viscosity
! !----------------------------------------------------------
! !
!$acc kernels if(compute_on_device) default(present)
      do j=Jstr,Jend
       do i=Istr,Iend
        WORK(i,j)=pm(i,j)*pn(i,j)
       enddo
! !
! !*******************************************
! !************* W: 1st timestep ************* 
! !*******************************************
! !
       if (FIRST_TIME_STEP) then
        do k=1,N-1
          do i=Istr,Iend
! !
! !******** Predictor @ t+dt/2 ***************      
! !
            cff0=WORK(i,j)*rw(i,j,k)
            wz(i,j,k,nnew)=( wz(i,j,k,nstp)*(Hz(i,j,k)+Hz(i,j,k+1))
     &                              +2.*dt*cff0) 
     &                     /  (Hz_half(i,j,k)+Hz_half(i,j,k+1)) 
#  ifdef MASKING
     &                           *rmask(i,j)
#  endif    
! !
! !******** RHS Slow to Fast *****************      
! !
#  if defined NBQ || (defined K3FAST_COUPLINGW_SCH0 && defined K3FAST_C3D_WSF)
            rw_int_nbq(i,j,k) =cff0
#  elif defined K3FAST_C3D_WSF
            rw_int_nbq(i,j,k) =cff0
#   ifndef K3FAST_COUPLINGW_SCH2            
            rw_intt_nbq(i,j,k,nstp)=cff0
#   else
            rw_intt_nbq(i,j,k,nstp)=cff0
#   endif
#  endif
! !
! !******** w(t) x dz ************************      
! !
            wz(i,j,k,indx)=wz(i,j,k,nstp)*0.5*(Hz(i,j,k  )+
     &                                         Hz(i,j,k+1))
          enddo
        enddo    
! !
! !== Special treatment for k=N because the control 
! !== volume at the top is Hz(N)/2
! !
        do i=Istr,Iend
! !
! !******** Predictor @ t+dt/2 ***************      
! !
            cff0=WORK(i,j)*rw(i,j,N)
            wz(i,j,N,nnew)=(wz(i,j,N,nstp)*Hz(i,j,N)+2.*dt*cff0)
     &                     /Hz_half(i,j,N) ! !<-- wz has units m.s-1
#  ifdef MASKING
     &                     *rmask(i,j)
#  endif     
! !
! !******** RHS Slow to Fast *****************      
! !
#  if defined NBQ || (defined K3FAST_COUPLINGW_SCH0 && defined K3FAST_C3D_WSF)
            rw_int_nbq(i,j,N) =cff0
#  elif defined K3FAST_C3D_WSF
            rw_int_nbq(i,j,N) =cff0
#   ifndef K3FAST_COUPLINGW_SCH2            
            rw_intt_nbq(i,j,N,nstp)=cff0
#   else
            rw_intt_nbq(i,j,N,nstp)=cff0
#   endif
#  endif     
! !
! !******** v(t) x dz ************************      
! !       
            wz(i,j,N,indx)=wz(i,j,N,nstp)*0.5*Hz(i,j,N) ! !<-- wz(indx) has units kg.m-1.s-1
! !
! !== Special treatment for k=0 because the control 
! !== volume at the top is Hz(1)/2
! !
! !
! !******** Predictor @ t+dt/2 ***************      
! !
#  ifdef NBQ_FREESLIP
            wz(i,j,0,nnew)=  0.5*
     &                       ( u(i,j,1,nnew) * pm_u(i,j) 
     &                                    * ( z_w(i  ,j,0)-z_w(i-1,j,0))
     &                        +u(i+1,j,1,nnew) * pm_u(i+1,j) 
     &                                    * ( z_w(i+1,j,0)-z_w(i  ,j,0))
     &                        +v(i,j,1,nnew) * pm_v(i,j) 
     &                                    * ( z_w(i,j  ,0)-z_w(i,j-1,0)) 
     &                        +v(i,j+1,1,nnew) * pm_v(i,j+1) 
     &                                    * ( z_w(i,j+1,0)-z_w(i,j,0)))
! !
! !******** RHS Slow to Fast *****************      
! !
            rw_int_nbq (i,j,0)=0.
! !
! !******** w(t) x dz ************************      
! !
            wz(i,j,0,indx)=wz(i,j,0,nnew)*0.5*Hz(i,j,1) ! !<-- wz(indx) has units kg.m-1.s-1
#  else
            wz(i,j,0,nnew)=0.
            wz(i,j,0,indx)=0.
#  endif
        enddo  
! !
! !*******************************************
! !************* W: timesteps > 1 ************ 
! !*******************************************
! !
      else 
        do k=1,N-1
          do i=Istr,Iend
! !
! !******** Predictor @ t+dt/2 ***************      
! !
            qdm_nstp = wz(i,j,k,nstp)*(Hz(i,j,k  )+             
     &                                             Hz(i,j,k+1))
            qdm_indx = wz(i,j,k,indx)*(Hz_bak(i,j,k  )+     
     &                                             Hz_bak(i,j,k+1))
            wz(i,j,k,nnew)=( (cff14*qdm_nstp
     &                       +cff15*qdm_indx)  
     &                      +2.*cdt*WORK(i,j)       
     &                             *(rw(i,j,k)
#  if defined NBQ || defined K3FAST_C3D_WFS
     &                               +rw_nbq_avg1(i,j,k)  
#  endif     
     &                     )) / (Hz_half(i,j,k)+Hz_half(i,j,k+1))  ! !<-- wz has units m.s-1
! !
! !******** RHS Slow to Fast *****************      
! !
#  if defined NBQ || defined K3FAST_COUPLINGW_SCH0
           rw_int_nbq(i,j,k)=((qdm_nstp-qdm_indx)*0.5/dt
     &                 - WORK(i,j)*rw_nbq_avg2(i,j,k))
#    ifdef MASKING
     &                              *rmask(i,j)
#    endif     
#  elif defined K3FAST_C3D_WSF
#   ifndef K3FAST_COUPLINGW_SCH2  
           cff0 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                 - WORK(i,j)*rw_nbq_avg2(i,j,k))
           rw_int_nbq(i,j,k)=(cff6*cff0
     &                       +cff7*rw_intt_nbq(i,j,k,3-nstp)
     &                       +cff8*rw_intt_nbq(i,j,k,nstp))
#    ifdef MASKING
     &                              *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,k,nstp)=cff0
#   else
            cff0=WORK(i,j)*rw(i,j,k) ! ! wzdiff here if necessary
            rw_int_nbq(i,j,k)= (cff1*cff0
     &                        + cff2*rw_intt_nbq(i,j,k,3-nstp)
     &                        + cff13*rw_intt_nbq(i,j,k,nstp))
#    ifdef MASKING
     &                              *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,k,nstp)=cff0
#   endif /* K3FAST_COUPLINGW_SCH2 */
! !
! !******** Test dissipation *****************      
! !    
#   ifdef K3FAST_DISSCOUPLING
           cff10 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                 - WORK(i,j)*rw_nbq_avg2(i,j,k))
           rw_int_nbq(i,j,k)=(1.-alpha_w)*rw_int_nbq(i,j,k)   
     &                      +alpha_w * cff10
#   endif
#  endif /* NBQ ||K3FAST_C3D_WSF */
! !
! !******** w(t) x dz ************************      
! !
            wz(i,j,k,indx)=wz(i,j,k,nstp)*0.5*(Hz(i,j,k  )+
     &                                         Hz(i,j,k+1)) 
          enddo
        enddo
        do i=Istr,Iend
! !
! !== Special treatment for k=N because 
! !== the control volume at the top is Hz(N)/2  
! !
! !
! !******** Predictor @ t+dt/2 ***************      
! !
           qdm_nstp = wz(i,j,N,nstp)*Hz    (i,j,N)
           qdm_indx = wz(i,j,N,indx)*Hz_bak(i,j,N)
           wz(i,j,N,nnew)=((cff14*qdm_nstp
     &                     +cff15*qdm_indx)  
     &                     +2.*cdt*WORK(i,j)
     &                        *(rw(i,j,N)
#  if defined NBQ || defined K3FAST_C3D_WFS
     &                         +rw_nbq_avg1(i,j,N)
#  endif     
     &                   ))/ Hz_half(i,j,N)          ! !<-- wz has units m.s-1   
! !
! !******** RHS Slow to Fast *****************      
! !
#  if defined NBQ || defined K3FAST_COUPLINGW_SCH0
           rw_int_nbq(i,j,N)= ((qdm_nstp-qdm_indx)*0.5/dt
     &                      - WORK(i,j)*rw_nbq_avg2(i,j,N))
#    ifdef MASKING
     &                           *rmask(i,j)
#    endif     
#  elif defined K3FAST_C3D_WSF
#   ifndef K3FAST_COUPLINGW_SCH2  
           cff0 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                      - WORK(i,j)*rw_nbq_avg2(i,j,N))
           rw_int_nbq(i,j,N)=(cff6*cff0
     &                       +cff7*rw_intt_nbq(i,j,N,3-nstp)
     &                       +cff8*rw_intt_nbq(i,j,N,nstp))
#    ifdef MASKING
     &                           *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,N,nstp)=cff0
#   else
           cff0=WORK(i,j)*rw(i,j,N) 
           rw_int_nbq(i,j,N)= (cff1*cff0
     &                        +cff2 *rw_intt_nbq(i,j,N,3-nstp)
     &                        +cff13*rw_intt_nbq(i,j,N,nstp))
#    ifdef MASKING
     &                           *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,N,nstp)=cff0
#   endif /* K3FAST_COUPLINGW_SCH2 */
! !
! !******** Test dissipation *****************      
! !    
#   ifdef K3FAST_DISSCOUPLING
           cff10 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                      - WORK(i,j)*rw_nbq_avg2(i,j,N))
           rw_int_nbq(i,j,N)=  (1.-alpha_w)*rw_int_nbq(i,j,N)     ! ! ICI
     &                         +alpha_w*cff10
#   endif
#  endif /* NBQ || K3FAST_C3D_WSF */
! !
! !******** w(t) x dz ************************      
! !
           wz(i,j,N,indx)=wz(i,j,N,nstp)*0.5*Hz(i,j,N)   ! !<-- wz(indx) has units kg.m-1.s-1
! !
! !== Special treatment for k=0 because 
! ! the control volume at the bottom is Hz(1)/2  
! !
#  ifdef NBQ_FREESLIP
! !
! !******** Predictor @ t+dt/2 ***************      
! !
            wz(i,j,0,nnew)=  0.5*
     &                       ( u(i,j,1,nnew) * pm_u(i,j) 
     &                                    * ( z_w(i  ,j,0)-z_w(i-1,j,0))
     &                        +u(i+1,j,1,nnew) * pm_u(i+1,j) 
     &                                    * ( z_w(i+1,j,0)-z_w(i  ,j,0))
     &                        +v(i,j,1,nnew) * pm_v(i,j) 
     &                                    * ( z_w(i,j  ,0)-z_w(i,j-1,0)) 
     &                        +v(i,j+1,1,nnew) * pm_v(i,j+1) 
     &                                    * ( z_w(i,j+1,0)-z_w(i,j,0)))
! !
! !******** w(t) x dz ************************      
! !
          wz(i,j,0,indx)=wz(i,j,0,nnew)*0.5*Hz(i,j,1) ! !<-- wz(indx) has units kg.m-1.s-1
#  else
            wz(i,j,0,nnew)=0.
            wz(i,j,0,indx)=0.
#  endif

        enddo          
      endif
  
      enddo  
!$acc end kernels
! !==  
# endif /* NBQ || K3SLOW_W */   
  
# if defined M3FAST || defined K3FAST_C3D_UVFS
#  undef ru_nbq_avg1 
#  undef rv_nbq_avg1
# endif
# if defined NBQ || defined K3FAST_C3D_WFS
#  undef rw_nbq_avg1
# endif
! !
! !---------------------- Test ?---------------------------------
! !
# if defined TS_HADV_TEST && ( defined DIAGONAL_ADV || defined SOLID_BODY_ROT )
!$acc kernels if(compute_on_device) default(present)
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            u(i,j,k,nnew) = u(i,j,k,nstp)
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            v(i,j,k,nnew) = v(i,j,k,nstp)
          enddo
        enddo
      enddo
!$acc end kernels
# endif   
! !
! !-------------------------------------------------------------------
! ! Solid body ?
! !-------------------------------------------------------------------
! !
# if defined TS_HADV_TEST && defined SOLID_BODY_PER
      cff = cos(2.*pi*(time+dt/2.)/(float(ntimes)*dt))
!$acc kernels if(compute_on_device) default(present)
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            u(i,j,k,nnew) = u(i,j,k,nstp) * cff
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            v(i,j,k,nnew) = v(i,j,k,nstp) * cff
          enddo
        enddo
      enddo
!$acc end kernels
# endif
! !
! !===================================================================
! ! Set PHYSICAL lateral boundary conditions for tracers.
! !===================================================================
! !
# if defined TRACERS
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#  ifdef THREE_GHOST_POINTS_TS
        call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                               t(START_2D_ARRAY,1,nnew,itrc))
#   else
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
#  endif
# endif
      enddo
# endif

      call u3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
      call v3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
# ifdef K3SLOW_W
      call w3dbc_tile (Istr,Iend,Jstr,Jend, WORK)      
# endif
! !
! !===================================================================
! ! Coupling, include ghost points associated with PHYSICAL
! ! boundaries ONLY. Do not touch periodic ghost points or
! ! internal computational margins (MPI code).
! !===================================================================
! !
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
! !
! !===================================================================
! ! Set PHYSICAL lateral boundary conditions for momentum.
! !===================================================================
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#  ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v(START_2D_ARRAY,1,nnew))
#   if defined NBQ || defined K3SLOW_W
      call exchange_w3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             wz(START_2D_ARRAY,0,nnew))
#   endif
#   if defined M3FAST || defined K3FAST_C3D_UVSF
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &            ru_int_nbq(START_2D_ARRAY,1))   
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &            rv_int_nbq(START_2D_ARRAY,1))   
#   endif
#   if defined NBQ || defined K3FAST_C3D_WSF
      call exchange_w3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &            rw_int_nbq(START_2D_ARRAY,0))   
#   endif
#  else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,nnew))
#   if defined NBQ || defined K3SLOW_W   
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        wz(START_2D_ARRAY,0,nnew))   
#   endif
#   if defined M3FAST || defined K3FAST_C3D_UVSF
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &            ru_int_nbq(START_2D_ARRAY,1))   
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &            rv_int_nbq(START_2D_ARRAY,1))   
#   endif
#   if defined NBQ || defined K3FAST_C3D_WSF
       call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &            rw_int_nbq(START_2D_ARRAY,0))   
#   endif
#  endif
# endif

