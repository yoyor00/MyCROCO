
#  ifndef NBQ_THETAIMP
! !------------------------------------------------------------------
! !  Store qdmw_nbq into working array
! !------------------------------------------------------------------
! !
      if (LAST_FAST_STEP) then
      
! ! KERNEL_18  rw_nbq <=  ( qdmw_nbq )

!$acc kernels default( present )     
        do k=0,N 
          do j=Jstr,Jend             
            do i=Istr,Iend
               rw_nbq(i,j,k)=qdmw_nbq(i,j,k)
            enddo
          enddo
        enddo
#   ifdef M3FAST_SEDLAYERS
        do k=-N_SL,-1
          do j=Jstr,Jend             
            do i=Istr,Iend
               rw_nbq(i,j,k)=qdmw_nbq(i,j,k)
            enddo
          enddo
        enddo
#   endif
!$acc end kernels
!$acc update host( rw_nbq )     !!iif=last
      endif
#  endif  /* NBQ_THETAIMP */
#  ifdef NBQ_IMP

! ! KERNEL_23  FC <= ( soundspeed2_nbq, rho_nbq, thetadiv_nbq )
! ! KERNEL_23  FC <= ( FC, Hzr_nbq_inv )
! ! KERNEL_23  dum_s <= ( FC )
! ! KERNEL_23  qdmw_nbq <= ( qdmw_nbq, dm_s,rw_int_nbq )
! ! KERNEL_23  qdmw_nbq <= ( qdmw_nbq, rmask)

!$acc kernels default( present )
      do j=Jstr,Jend
        do i=Istr,Iend
            FC3D(i,j,N+1)=0
        enddo
      enddo
#  ifndef M3FAST_SEDLAYERS        
      do k=1,N
#  else   
      do k=-N_sl+1,N
#  endif
      !do j=Jstr,Jend ! <-- j loop
! !
! !....................................................................
! ! Vertical fluxes
! !....................................................................
! !
!         do k=1,N
          do j=Jstr,Jend ! <-- j loop
          do i=Istr,Iend
#  ifndef M3FAST_CSVISC2K
            FC3D(i,j,k)= soundspeed2_nbq*rho_nbq(i,j,k)
     &           -( thetaimp_nbq*soundspeed2_nbq*dtfast
     &             +visc2v_nbq)
     &           *thetadiv_nbq(i,j,k)   
#  else
            FC3D(i,j,k)= soundspeed2_nbq(i,j,k)*rho_nbq(i,j,k)
     &           -( thetaimp_nbq*soundspeed2_nbq(i,j,k)*dtfast
     &             +visc2v_nbq(i,j,k))
     &           *thetadiv_nbq(i,j,k)   
#  endif 
#  ifdef NBQ_THETAIMP 
#  ifdef NBQ_MASS
	stop "CAUTION when THETAIMP and NBQ_MASS."
#  endif
            FC3D(i,j,k)=FC3D(i,j,k)
     &       -thetaimp_nbq*(1.-thetaimp_nbq)*dtfast
     &                    *( Hzw_nbq_inv(i,j,k  )*qdmw_nbq(i,j,k)      ! CAUTION when used with NBQ_MASS
     &                      -Hzw_nbq_inv(i,j,k-1)*qdmw_nbq(i,j,k-1) )
#  ifndef M3FAST_CSVISC2K
     &                    *soundspeed2_nbq
#  else
     &                    *soundspeed2_nbq(i,j,k)
#  endif
#  endif
            FC3D(i,j,k)=FC3D(i,j,k)*Hzr_nbq_inv(i,j,k) 
          enddo
        enddo
      enddo 
      
!$acc end kernels

#  endif /* NBQ_IMP */
! !
! !--------------------------------------------------------------------
! ! Inner layers (ocean)
! !--------------------------------------------------------------------
! !
!$acc kernels default( present )      
!#  ifndef NBQ_IMP
#  if !defined  M3FAST_SPDUP && defined M3FAST_BOTH 
      if (IstrU.le.Iend) then
        do j=Jstr,Jend
          do i=IstrU-1,Iend+1     
            usurf_nbq(i,j)=(qdmu_nbq(i,j,N)*Hzu_nbq_inv(i,j,N)
#   ifdef MRL_WCI
     &                     +ust(i,j,N)
#    ifdef WET_DRY
     &                      *umask_wet(i,j)
#    endif
#   endif
     &                     )
#   ifdef MASKING
     &                     *umask(i,j) 
#   endif
          enddo 
        enddo 
      endif

      if (JstrV.le.Jend) then
        do j=JstrV-1,Jend+1
          do i=Istr,Iend     
            vsurf_nbq(i,j)=(qdmv_nbq(i,j,N)*Hzv_nbq_inv(i,j,N)
#   ifdef MRL_WCI
     &                     +vst(i,j,N)
#    ifdef WET_DRY
     &                      *vmask_wet(i,j)
#    endif
#   endif
     &                      )
#   ifdef MASKING
     &                     *vmask(i,j) 
#   endif
          enddo
        enddo 
      endif

#   if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      if (IstrU.le.Iend) then
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          usurf_nbq(START_2D_ARRAY))
      endif
      if (JstrV.le.Jend) then
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          vsurf_nbq(START_2D_ARRAY))
      endif
#   endif
#  endif  /* M3FAST_SPDUP && M3FAST_BOTH */
!
      do j=Jstr,Jend   ! Large j loop starting here
        do k=1,N-1
          do i=Istr,Iend  
! !
! ! Compute vertical gradient
! !
#   ifndef NBQ_IMP
            dum_s = thetadiv_nbq(i,j,k) - thetadiv_nbq(i,j,k+1) 
#   else
            dum_s = FC3D(i,j,k) - FC3D(i,j,k+1)
#   endif
! !
! ! Compute Non-traditionnal component of Coriolis
! !
#   ifdef UV_COR_NT
            dum_s = dum_s + ntcorw(i,j,k)
#   endif
! !
! ! Compute weight
! !
#   ifdef NBQ_GRAV
            cff=sign(1.,qdmw_nbq(i,j,k))
            dum_s = dum_s  
     &            -0.5*g*Hzw_nbq(i,j,k)*(
     &      ((cff+abs(cff))*rho_nh(i,j,k  )
     &      -(cff-abs(cff))*rho_nh(i,j,k+1))/rho0
     &      +(cff+abs(cff))*rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k)
     &      -(cff-abs(cff))*rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1)
     &                                  )
#   endif
! ! 
! !  Update qdmw_nh (fast and slow components)
! !
            qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)   
     &                      + dtfast * ( dum_s 
#   ifdef M3FAST_C3D_WSF
     &                      + rw_int_nbq(i,j,k)
#   endif     
     &                                )
! ! 
! !  Masking
! !
#   ifdef MASKING
          qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k) * rmask(i,j)
#   endif  
#  if defined NBQ_NUDGING && defined NBQCLIMATOLOGY
           qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)*(1.-NBQnudgcof(i,j))
     &                        +wz(i,j,k,nrhs)
     &                         *0.5*(Hzr(i,j,k)+Hzr(i,j,k+1))*pm(i,j)
     &                         *NBQnudgcof(i,j)
#   ifdef MASKING
     &                         *rmask(i,j)
#   endif  
#  endif
          enddo  ! i loop           
        enddo    ! k loop
!
#   ifdef M3FAST_SEDLAYERS
! !
! !--------------------------------------------------------------------
! ! Sediment Layers 
! !--------------------------------------------------------------------
! !
        do k=-N_sl+1,-1
          do i=Istr,Iend  
! !
! ! Compute vertical gradient
! !
#   ifndef NBQ_IMP
            dum_s = thetadiv_nbq(i,j,k) - thetadiv_nbq(i,j,k+1) 
#   else
            dum_s = FC3D(i,j,k) - FC3D(i,j,k+1)
#   endif
            qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k) + dtfast * dum_s 
          enddo
        enddo
! !
! !--------------------------------------------------------------------
! ! Interface layer:
! !     1/2 ocean and 1/2 SdL
! !--------------------------------------------------------------------
! !
        k = 0
        do i=Istr,Iend           
!
#   ifndef NBQ_IMP
            dum_s = thetadiv_nbq(i,j,k) - thetadiv_nbq(i,j,k+1) 
#   else
            dum_s = FC3D(i,j,k) - FC3D(i,j,k+1)
#   endif
!
!#   ifdef UV_COR_NT
!            dum_s = dum_s + ntcorw(i,j,k)
!#   endif
!#   ifdef NBQ_THETAIMP
!            rw_nbq(i,j,k)=qdmw_nbq(i,j,k)
!#   endif 
#   ifdef NBQ_GRAV
!            cff=sign(1.,qdmw_nbq(i,j,k))     ! TBF CAUTION
!            dum_s = dum_s  
!     &            -0.25*g*Hzw_nbq(i,j,k)*(    ! /2.
!     &      ((cff+abs(cff))*rho_nh(i,j,k  )
!     &      -(cff-abs(cff))*rho_nh(i,j,k+1))/rho0
!     &       (cff+abs(cff))*rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k)
!     &      -(cff-abs(cff))*rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1)
!     &                                  )
#   endif
            qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)   
     &                      + dtfast * ( dum_s 
#   ifdef M3FAST_C3D_WSF
!     &                      + rw_int_nbq(i,j,k)  ! Slow mode is not integrated here!
#   endif     
     &                                )
#   ifdef MASKING
          qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k) * rmask(i,j)
#   endif  
        enddo ! i loop
#   endif
! !
! !--------------------------------------------------------------------
! ! Bottom BC on w
! !--------------------------------------------------------------------
! !
       k = -N_sl
#  ifdef NBQ_FREESLIP
          do i=Istr,Iend    
                qdmw_nbq(i,j,k)=
#    if defined MVB && ! defined M3FAST_SEDLAYERS
     &       -0.5*(dh_mvb(i,j,kbak2)-dh_mvb(i,j,kstp2))/dtfast
#    endif
     &       +(0.5*( qdmu_nbq(i  ,j,k+1)*Hzu_nbq_inv(i  ,j,k+1)
     &             *(z_w(i,j,k)-z_w(i-1,j,k))
     &             +qdmu_nbq(i+1,j,k+1)*Hzu_nbq_inv(i+1,j,k+1)
     &             *(z_w(i+1,j,k)-z_w(i,j,k))
     &           )*pm_u(i,j)
     &       +0.5*( qdmv_nbq(i,j  ,k+1)*Hzv_nbq_inv(i,j  ,k+1)
     &             *(z_w(i,j,k)-z_w(i,j-1,k))
     &             +qdmv_nbq(i,j+1,k+1)*Hzv_nbq_inv(i,j+1,k+1)
     &             *(z_w(i,j+1,k)-z_w(i,j,k))
     &           )*pm_v(i,j))
     &           *Hzw_nbq(i,j,k)
#    ifdef MASKING
           qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)*rmask(i,j)
#    endif 
        enddo
#  else
          do i=Istr,Iend   
#    ifdef MVB
          qdmw_nbq(i,j,k) = 0.5*w_mvb(i,j,knew2)*Hz(i,j,k+1)   ! CAUTION HERE
#    else
          qdmw_nbq(i,j,k) = 0.
#    endif
         enddo 
#  endif
! !
! !--------------------------------------------------------------------
! ! Surface boundary condition
! !--------------------------------------------------------------------
! !
        k=N
        do i=Istr,Iend
! ! 
! !  Vertical gradient
! !
#    ifndef NBQ_IMP
          dum_s =   thetadiv_nbq(i,j,N)
#    else
          dum_s =   FC3D(i,j,k) 
#    endif
! ! 
! !  Non-traditional component of Coriolis
! !
#    ifdef UV_COR_NT
          dum_s = dum_s + ntcorw(i,j,k)
#    endif
! !
! !  Compute weight
! !
#    ifdef NBQ_GRAV
          dum_s = dum_s 
     &         -g*Hzw_nbq(i,j,N)*(
     &            rho_nh (i,j,N)/rho0
     &           +rho_nbq(i,j,N)*Hzr_nbq_inv(i,j,N) 
     &                           )
#    endif    
! ! 
! !  Update qdmw(N): fast and slow components
! !
          qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)   
     &                   + dtfast * ( dum_s 
#    ifdef M3FAST_C3D_WSF
     &                   + rw_int_nbq(i,j,N)
#    endif
     &                              )
! !
! !  Masking
! !
#   ifdef MASKING
          qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N) * rmask(i,j)
#   endif
#   if defined NBQ_NUDGING && defined NBQCLIMATOLOGY
           qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)*(1.-NBQnudgcof(i,j))
     &                        +wz(i,j,N,nrhs)
     &                         *Hzr(i,j,N)*pm(i,j)
     &                         *NBQnudgcof(i,j)
#    ifdef MASKING
     &                         *rmask(i,j)
#    endif
#   endif  
     
        enddo
! !
! !  End of Large j loop
! !
      enddo !<-- j loop
!$acc end kernels
! !
! !--------------------------------------------------------------------
! ! Vertical momentum open boundary conditions
! !--------------------------------------------------------------------
! !
!      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)

!#   ifdef UV_COR_NT
!#    if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!#      ifndef M3FAST_SEDLAYERS
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmw_nbq(START_2D_ARRAY,0))   
!#      else
!      call exchange_w3d_sedlayb_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmw_nbq(START_2D_ARRAY,-N_sl))   
!#      endif
!#    endif
!#   endif
#   ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
#    ifdef M3FAST_SEDLAYER      
      call check_tab3d_sedlay(qdmw_nbq,'qdmw_nbq','wint',-N_sl,N)
      call check_tab3d_sedlay(Hz,'Hz','r',-N_sl+1,N)
#    endif      
C$OMP END MASTER
#   endif
!#  endif  /* ! NBQ_IMP */
