! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_divh.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !
! !********************************
! !********************************
! !        Divergence 
! ! Precompute dZdx*qdmu terms 
! !       & Bottom BC
! !********************************
! !********************************
! !
! !********************************
! ! Initializations
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)

      if (IstrU.gt.Iend) then
        do j=Jstr,Jend
          do i=Istr,Iend+1
            FX(i,j)=0.
          enddo
        enddo
      endif
      if (JstrV.gt.Jend) then
        do j=Jstr,Jend+1
          do i=Istr,Iend 
            FY(i,j)=0.
          enddo
        enddo
      endif
!$acc end kernels
! !
! !********************************
! !  X-component dZdxq_u*qdmu for 
! !    horizontal divergence
! !    Algo: 
! !    - if grid-fast or grid-slow 
! !      at update time: d/dstermes are 
! !      computed
! !    - if grid-slow with no update:
! !      only BBC
! !    - if grid-fast: BBC
! !    - d/ds terms are added to
! !      horizontal divergence at m+1
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)

#   ifdef NBQ_GRID_SLOW
      if (NSTEP_DS) then
#   endif
! !
! !--------------------------------
! ! Pre-commpute 
! !     dZdxq_u = \Delta z * qdmu
! !--------------------------------
! !
      do k=-N_sl,N-1
        kp1 = k + 1
        do j=Jstr,Jend
!            do i=IstrU,Iend+1
            do i=Istr,Iend+1
                dZdxq_u(i,j,k)= (z_r(i,j,kp1)-z_r(i-1,j,kp1)) 
     &                         *qdmu_nbq(i,j,kp1)  
            enddo
        enddo
      enddo
! !
! !--------------------------------
! ! Inner layers (-Nsl<k<N )  :  
! !         \Sigma k
! !--------------------------------
! !
      do k=-N_sl+1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend+1
                dZdxq_w(i,j,k )= Hzw_nbq_inv_u(i,j,k)
     &                          *(dZdxq_u(i,j,k-1)+dZdxq_u(i,j,k)) 
#    ifdef MASKING
     &                          *umask(i,j)
#    endif
          enddo 
         enddo
       enddo
! !
! !--------------------------------
! ! k = N   :   \Sigma k
! !             special scheme 
! !             to obtain Surface 
! !             Kinematic Relation
! !--------------------------------
! ! 
        do j=Jstr,Jend
            do i=Istr,Iend+1
                dZdxq_w(i,j,N)= (z_w(i,j,N)-z_w(i-1,j,N))
     &                           *qdmu_nbq(i,j,N)
     &                           /(Hz(i,j,N)+Hz(i-1,j,N)) 
#    ifdef MASKING
     &                           *umask(i,j)
#    endif
             enddo
        enddo  
#   ifdef NBQ_GRID_SLOW
      endif
#   endif
! !
! !--------------------------------
! ! k = -N_sl  : \Sigma k
! !--------------------------------
! !
      k = -N_sl
!
#   ifdef NBQ_FREESLIP
            do j=Jstr,Jend
              do i=Istr,Iend+1 
!                 dZdxq_w(i,j,k )= (z_w(i,j,k)-z_w(i-1,j,k))
!     &           * (
#    ifdef MVB
!     &              -2.*u_mvb(i,j)
#    endif
!     &              +qdmu_nbq(i,j,k+1)/(Hzr(i,j,k+1)+Hzr(i-1,j,k+1)))
#    ifdef MASKING
!     &              *umask(i,j)
#    endif
                 dZdxq_w(i,j,k) = 
     &           0.5*(z_w(i,j,k)-z_w(i-1,j,k))*(
     &           qdmu_nbq(i,j,k+1)
     &           *Hzu_nbq_inv(i,j,k+1)   
#    ifdef MVB
!    &           +u_mvb(i,j,knew2)   ! MVBFA to be finished
#    endif
     &           )
              enddo
            enddo
#   else /* ! NBQ_FREESLIP */
            do j=Jstr,Jend          ! No-slip condition
              do i=Istr,Iend+1
#    ifdef MVB
                dZdxq_w(i,j,k)=0.5*(z_w(i,j,k)-z_w(i-1,j,k))
     &                            *u_mvb(i,j,knew2) 
#    ifdef MASKING
     &              *umask(i,j)
#    endif
#    else   /* ! MVB */
                dZdxq_w(i,j,k)=0.  
#    endif  /* MVB */    
              enddo
            enddo
            ! Bottom BC for w in the following.
#   endif /* NBQ_FREESLIP */
! !
! !--------------------------------
! !       \Sigma i
! !--------------------------------
! ! FX ou FC3D ou private( FX )  à regarder
! !
! !
! ! Ocean layer (inner domain)
! !
          do k=1,N-1
            do j=Jstr,Jend
              do i=Istr,Iend+1
                FC3D(i,j,k)=-pm_u(i,j)*(dZdxq_w(i,j,k )-dZdxq_w(i,j,k-1))
#   ifdef MASKING
     &                      *umask(i,j)
#   endif                
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=FC3D(i,j,k)+FC3D(i+1,j,k)
              enddo
            enddo
          enddo ! k=1,N-1
! !
! ! Ocean layer (surface)
! !
            do j=Jstr,Jend
              do i=Istr,Iend+1
                FC3D(i,j,N)=-pm_u(i,j)*(-dZdxq_w(i,j,N-1))
#   ifdef MASKING
     &                      *umask(i,j)
#   endif                
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,N)=FC3D(i,j,N)+FC3D(i+1,j,N)
#   ifdef NBQ_MASS
     &         -(1.+rho_grd(i,j,N))*
#   else
     &         -
#   endif     
     &         (pm_u(i,j)*dZdxq_w(i,j,N)
#   ifdef MASKING
     &                      *umask(i,j)
#   endif                
     &          +pm_u(i+1,j)*dZdxq_w(i+1,j,N)
#   ifdef MASKING
     &                      *umask(i+1,j)
#   endif                
     &          )
              enddo
            enddo
          
#  ifdef M3FAST_SEDLAYERS 
! !
! ! Sediment layer (inner domain)
! !
          do k=-N_sl+1,0
            do j=Jstr,Jend
              do i=Istr,Iend+1
                FC3D(i,j,k)=-pm_u(i,j)*(dZdxq_w(i,j,k )-dZdxq_w(i,j,k-1))
#   ifdef MASKING
     &             *umask(i,j)
#   endif                
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=FC3D(i,j,k)+FC3D(i+1,j,k)
              enddo
            enddo
          enddo ! k=-N_sl+1,-1
          
#   ifdef NBQ_FREESLIP
! !
! ! Interface in SdL
! !
          do j=Jstr,Jend
            do i=Istr,Iend+1
              k = 1
              FC3D(i,j,k)=-pm_u(i,j)*(dZdxq_w(i,j,k )
     &         -(z_w(i,j,0)-z_w(i-1,j,0))
     &         *qdmu_nbq(i,j,1)/(Hzr(i,j,1)+Hzr(i-1,j,1))
     &                           )
#    ifdef MASKING
     &              *umask(i,j)
#    endif
              k = 0
              FC3D(i,j,k)=-pm_u(i,j)*(
     &              (z_w(i,j,0)-z_w(i-1,j,0))
     &         *qdmu_nbq(i,j,0)/(Hzr(i,j,0)+Hzr(i-1,j,0)) 
     &              -dZdxq_w(i,j,k-1) 
     &                           )
#    ifdef MASKING
     &              *umask(i,j)
#    endif
            enddo
          enddo
          do j=Jstr,Jend
             do i=Istr,Iend
                k = 1
                thetadiv_nbq(i,j,k)=FC3D(i,j,k)+FC3D(i+1,j,k)
                k = 0
                thetadiv_nbq(i,j,k)=FC3D(i,j,k)+FC3D(i+1,j,k)
             enddo
          enddo
#   endif /* NBQ_FREESLIP */
#  endif   /* M3FAST_SEDLAYERS */
!$acc end kernels
! !
! !********************************
! !  Y-component dZdyq_v*qdmv 
! !    for horizontal divergence
! !    Algo:
! !    - if grid-fast or grid-slow 
! !       at update time: d/dstermes are 
! !      computed
! !    - if grid-slow with no update: only BBC
! !    - if grid-fast: BBC
! !    - d/ds terms are added to horizontal 
! !      divergence at m+1
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)

#   ifdef NBQ_GRID_SLOW
      if (NSTEP_DS) then
#   endif
! !
! !--------------------------------
! ! Pre-commpute dZdxq_v
! !--------------------------------
! !
      do k=-N_sl,N-1
        kp1 = k + 1
            do j=Jstr,Jend+1
              do i=Istr,Iend
                dZdyq_v(i,j,k)=(z_r(i,j,kp1)-z_r(i,j-1,kp1)) 
     &                           *qdmv_nbq(i,j,kp1)  ! (dZdy * (rho v))_v
              enddo
            enddo
      enddo
! !
! !--------------------------------
! ! Inner layers (-Nsl<k<N )
! !--------------------------------
! !
      do k=-N_sl+1,N-1  
        do j=Jstr,Jend+1
          do i=Istr,Iend
                dZdyq_w(i,j,k )=Hzw_nbq_inv_v(i,j,k)
     &                *(dZdyq_v(i,j,k-1)+dZdyq_v(i,j,k)) 
#    ifdef MASKING
     &              *vmask(i,j)
#    endif
          enddo 
         enddo
       enddo
! !
! !--------------------------------
! ! k = N   :   \Sigma k
! !             special scheme to obtain 
! !             Surface Kinematic Relation
! !--------------------------------
! ! 
            do j=Jstr,Jend+1
              do i=Istr,Iend
                dZdyq_w(i,j,N )= (z_w(i,j,N)-z_w(i,j-1,N))
     &                          *qdmv_nbq(i,j,N)
     &                          /(Hz(i,j,N)+Hz(i,j-1,N))
#    ifdef MASKING
     &              *vmask(i,j)
#    endif
              enddo  
            enddo
#   ifdef NBQ_GRID_SLOW
      endif
#   endif 
! !
! !--------------------------------
! ! k = -N_sl  : \Sigma k
! !--------------------------------
! !
            k=-N_sl
!
#   ifdef NBQ_FREESLIP
            do j=Jstr,Jend+1
              do i=Istr,Iend
!                 dZdyq_w(i,j,k )= (z_w(i,j,k)-z_w(i,j-1,k))
!     &           * (
#    ifdef MVB
!     &              -2.*v_mvb(i,j)
#    endif
!     &            +qdmv_nbq(i,j,k+1)/(Hzr(i,j,k+1)+Hzr(i,j-1,k+1)))
#    ifdef MASKING
!     &              *vmask(i,j)
#    endif
                 dZdyq_w(i,j,k) = 
     &           0.5*(z_w(i,j,k)-z_w(i,j-1,k))
     &           *qdmv_nbq(i,j,k+1)
     &           *Hzv_nbq_inv(i,j,k+1)       
              enddo
            enddo
#   else /* ! NBQ_FREESLIP */
            do j=Jstr,Jend+1          ! No-slip condition
              do i=Istr,Iend
#    ifdef MVB
                dZdyq_w(i,j,k)=0.5*(z_w(i,j,k)-z_w(i,j-1,k))
     &                            *v_mvb(i,j,knew2) 
#    ifdef MASKING
     &              *vmask(i,j)
#    endif
#    else   /* ! MVB */
                dZdyq_w(i,j,k)=0.  
#    endif  /* MVB */    
              enddo
            enddo
            ! Bottom BC for w in the following.
#   endif /* NBQ_FREESLIP */
!! <-- k=0
! !
! !--------------------------------
! !       \Sigma i
! !--------------------------------
! !
! !
! ! Ocean layer (inner domain)
! !
          do k=1,N-1
!          if (JstrV.le.Jend) then
            do j=Jstr,Jend+1
              do i=Istr,Iend
                FC3D(i,j,k)=-pn_v(i,j)*(dZdyq_w(i,j,k )-dZdyq_w(i,j,k-1))
#   ifdef MASKING
     &                      *vmask(i,j)
#   endif                
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                             +FC3D(i,j,k)+FC3D(i,j+1,k)
              enddo
            enddo
!	   endif 
          enddo ! k=1,N-1
      
! !
! ! Ocean layer (surface)
! !
            do j=Jstr,Jend+1
              do i=Istr,Iend
                FC3D(i,j,N)=-pn_v(i,j)*(-dZdyq_w(i,j,N-1))
#   ifdef MASKING
     &                      *vmask(i,j)
#   endif                
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,N)=FC3D(i,j,N)+FC3D(i,j+1,N)
#   ifdef NBQ_MASS
     &         -(1.+rho_grd(i,j,N))*
#   else
     &         -
#   endif     
     &         (pn_v(i,j)*dZdyq_w(i,j,N)
#   ifdef MASKING
     &                      *vmask(i,j)
#   endif                
     &          +pn_v(i+1,j)*dZdyq_w(i,j+1,N)
#   ifdef MASKING
     &                      *vmask(i,j+1)
#   endif                
     &          )
              enddo
            enddo    
            
#   ifdef M3FAST_SEDLAYERS
! !
! ! Sediment layer (inner domain)
! !
          do k=-N_sl+1,0
            do j=Jstr,Jend+1
              do i=Istr,Iend
                FC3D(i,j,k)=-pn_v(i,j)*(dZdyq_w(i,j,k )-dZdyq_w(i,j,k-1))
#   ifdef MASKING
     &             *vmask(i,j)
#   endif                
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                             +FC3D(i,j,k)+FC3D(i,j+1,k)
              enddo
            enddo
          enddo ! k=-N_sl+1,-1
          
#   ifdef NBQ_FREESLIP
! !
! ! Interface in SdL
! !
          do j=Jstr,Jend+1
            do i=Istr,Iend
              k = 1
              FC3D(i,j,k)=-pm_v(i,j)*(dZdyq_w(i,j,k )
     &         -(z_w(i,j,0)-z_w(i,j-1,0))
     &         *qdmv_nbq(i,j,1)/(Hzr(i,j,1)+Hzr(i,j-1,1))
     &                           )
#    ifdef MASKING
     &              *vmask(i,j)
#    endif
              k = 0
              FC3D(i,j,k)=-pm_v(i,j)*(
     &              (z_w(i,j,0)-z_w(i,j-1,0))
     &         *qdmv_nbq(i,j,0)/(Hzr(i,j,0)+Hzr(i,j-1,0)) 
     &              -dZdyq_w(i,j,k-1) 
     &                           )
#    ifdef MASKING
     &              *vmask(i,j)
#    endif
            enddo
          enddo
          do j=Jstr,Jend
             do i=Istr,Iend
                k = 1
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                              +FC3D(i,j,k)+FC3D(i,j+1,k)
                k = 0
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                              +FC3D(i,j,k)+FC3D(i,j+1,k)
             enddo
          enddo
#   endif /* NBQ_FREESLIP */
#  endif   /* M3FAST_SEDLAYERS */
!$acc end kernels
! !
! !********************************
! !  Horizontal Divergence (qdmH(m+1)): 
! !      add d/dx and D/dy terms
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)

#  ifndef M3FAST_SEDLAYERS
        do k=1,N!<-- k loop
#  else
        do k=-N_sl+1,N!<-- k loop
#  endif
        if (IstrU.le.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend+1
              FC3D(i,j,k)=on_u(i,j)*qdmu_nbq(i,j,k)
            enddo
          enddo
        endif
        if (JstrV.le.Jend) then
          do j=Jstr,Jend+1  
            do i=Istr,Iend
              DC3D(i,j,k)=om_v(i,j)*qdmv_nbq(i,j,k)
            enddo
          enddo
        endif
        if (IstrU.gt.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend  
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                           +pm(i,j)*pn(i,j)
     &                           *(DC3D(i,j+1,k)-DC3D(i,j,k))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        elseif (JstrV.gt.Jend) then
          do j=Jstr,Jend  
            do i=Istr,Iend   
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)
     &                            *(FC3D(i+1,j,k)-FC3D(i,j,k))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)
     &                            *(FC3D(i+1,j,k)-FC3D(i,j,k)
     &                            + DC3D(i,j+1,k)-DC3D(i,j,k))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        endif  

      enddo ! <-- k=1,N
!$acc end kernels
! ! 
! !********************************
! !  thetadiv2_nbq: complet time-corrective term  (dh/dt included) 
! !  thetadiv_nbq: reduced time-corrective term  (no dh/dt)
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)

#  ifndef HCOMP 
#   ifdef NBQ_GRID_SLOW
       if (FIRST_FAST_STEP_3M) then
#   endif
!$acc loop independent private ( FC, CF )
        do j=Jstr,Jend !<-- j loop
          do i=Istr,Iend

#   ifndef M3FAST_SEDLAYERS          
#    ifdef NBQ_HZ_PROGNOSTIC
            FC3D(i,j,0)=0.    ! Bottom boundary condition
#    endif
            CF3D(i,j,0)=0.
#   else
#    ifdef NBQ_HZ_PROGNOSTIC
            FC3D(i,j,-N_sl)=0.    ! Bottom boundary condition
#    endif
            CF3D(i,j,-N_sl)=0.
#   endif  /* M3FAST_SEDLAYERS */
          enddo
        enddo    
        do k=-N_sl+1,N-1
          do j=Jstr,Jend !<-- j loop
            do i=Istr,Iend
#   ifdef NBQ_HZ_PROGNOSTIC
              FC3D(i,k)=   
     &          -(z_w(i,j,k)-zw_nbq(i,j,k))/dtgrid_nbq
     &           *0.5*( (1.+rho_grd(i,j,k  ))
     &                  +rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k  )  
     &                  +(1.+rho_grd(i,j,k+1))
     &                  +rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1))
#   endif
              CF3D(i,j,k)=   
     &          -(z_w(i,j,k)-zw_nbq(i,j,k))/dtgrid_nbq     ! ATTENTION (Francis): these terms are under M3FAST_W
     &          *0.5*(    rho_grd(i,j,k  ) 
     &                   +rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k  ) 
     &                   +rho_grd(i,j,k+1)
     &                   +rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1) )

            enddo
          enddo
        enddo
        
#   ifndef M3FAST_SEDLAYERS        
        do k=1,N-1
#   else   
        do k=2,N-1
#   endif
          do j=Jstr,Jend !<-- j loop
            do i=Istr,Iend
              thetadiv_nbq(i,j,k) =thetadiv_nbq(i,j,k)+CF3D(i,j,k)-CF3D(i,j,k-1)
#   ifdef NBQ_HZ_PROGNOSTIC
              thetadiv2_nbq(i,j,k)=thetadiv_nbq(i,j,k)+FC3D(i,j,k)-FC3D(i,j,k-1) 
#   endif
              zw_nbq(i,j,k)=z_w(i,j,k)
            enddo
          enddo
        enddo
        
#   ifdef M3FAST_SEDLAYERS
        do k=-N_sl+1,-1
          do j=Jstr,Jend !<-- j loop
            do i=Istr,Iend
!             thetadiv_nbq(i,j,k) =thetadiv_nbq(i,j,k)+CF3D(i,j,k)-CF3D(i,j,k-1)
#    ifdef NBQ_HZ_PROGNOSTIC
!             thetadiv2_nbq(i,j,k)=thetadiv_nbq(i,j,k)+FC3D(i,j,k)-FC3D(i,j,k-1) 
#    endif
              zw_nbq(i,j,k)=z_w(i,j,k)
            enddo
          enddo
        enddo
#   endif /* M3FAST_SEDLAYERS */  

#   ifdef M3FAST_SEDLAYERS
          do j=Jstr,Jend !<-- j loop
            do i=Istr,Iend
              k=0
              thetadiv_nbq(i,j,k) =thetadiv_nbq(i,j,k)
     &          -CF3D(i,j,k-1) *0
     &          -(z_w(i,j,k)-zw_nbq(i,j,k))/dtgrid_nbq    
     &          *(  rho_grd(i,j,k  ) 
     &                 +rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k  )  )
#    ifdef NBQ_HZ_PROGNOSTIC
              thetadiv2_nbq(i,j,k)=thetadiv_nbq(i,j,k)+FC3D(i,j,k)-FC3D(i,j,k-1) 
#    endif
              zw_nbq(i,j,k)=z_w(i,j,k)
              k=1
              thetadiv_nbq(i,j,k) =thetadiv_nbq(i,j,k)+CF3D(i,j,k) 
     &          +(z_w(i,j,k-1)-zw_nbq(i,j,k-1))/dtgrid_nbq  
     &          *(    rho_grd(i,j,k  ) 
     &                   +rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k  ) )
#    ifdef NBQ_HZ_PROGNOSTIC
              thetadiv2_nbq(i,j,k)=thetadiv_nbq(i,j,k)+FC3D(i,j,k)-FC3D(i,j,k-1) 
              stop 'To be done'
#    endif
              zw_nbq(i,j,k)=z_w(i,j,k)
            enddo
          enddo
#   endif /* M3FAST_SEDLAYERS */       
         
          do j=Jstr,Jend !<-- j loop
            do i=Istr,Iend
#   ifdef NBQ_HZ_PROGNOSTIC
             FC3D(i,j,N)=   
     &        -(z_w(i,j,N)-zw_nbq(i,j,N))/dtgrid_nbq
     &         *( 1.+rho_grd(i,j,N)
     &          )
#   endif
            CF3D(i,j,N)=  
     &        -(z_w(i,j,N)-zw_nbq(i,j,N))/dtgrid_nbq
     &         *( rho_grd(i,j,N)
     &          )
            thetadiv_nbq(i,j,N)=thetadiv_nbq(i,j,N)+CF3D(i,j,N)-CF3D(i,j,N-1) 
#   ifdef NBQ_HZ_PROGNOSTIC
            thetadiv2_nbq(i,j,N)=thetadiv_nbq(i,j,N)+FC3D(i,j,N)-FC3D(i,j,N-1) 
#   endif
              zw_nbq(i,j,N)=z_w(i,j,N)
            enddo
          enddo 
          
#   ifdef NBQ_GRID_SLOW
       endif !<-- FIRST_FAST_STEP
#   endif
#  endif /* HCOMP */
!$acc end kernels
! !
! !********************************
! ! Bottom BC (and interface BC in SdL) on w 
! ! both with explicit and implicit schemes
! !********************************
! !
#  ifdef NBQ_FREESLIP
#   ifndef M3FAST_SEDLAYERS
       k = 0
#   else
       k = -N_sl
#   endif
            do j=Jstr,Jend  
              do i=Istr,Iend    
                qdmw_nbq(i,j,k)=
#   if defined MVB && ! defined M3FAST_SEDLAYERS
     &       -0.5*(dh_mvb(i,j,kbak2)-dh_mvb(i,j,kstp2))/dtfast
#   endif
     &                         +0.5*( dZdxq_w(i  ,j,k)*pm_u(i  ,j)
     &                               +dZdxq_w(i+1,j,k)*pm_u(i+1,j) )
     &                             * Hzr(i,j,k+1)
#   ifdef MASKING
     &                             *rmask(i,j)
#   endif 
              enddo
            enddo 
#  endif
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_divh.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
     
