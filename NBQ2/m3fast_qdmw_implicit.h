! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_qdmw_implicit.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! !********************************
! ! Tridiag (implicit) system RHS:
! ! Gaussian Elimination
! !********************************
! !
!$acc kernels default( present )
! !
! !********************************
! !..........Bottom BC
! !********************************
! !
#   ifndef M3FAST_SEDLAYERS
        k=1
#   else
        k=-N_sl+1
#   endif
        do j=Jstr,Jend  !<-- j loop
          do i=Istr,Iend
            !-------------------------
            ! \beta (normalization)
            !-------------------------      
#   ifndef M3FAST_CSVISC2K
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtfast+visc2v_nbq)) 
#   else
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k)
     &              *dtfast+visc2v_nbq(i,j,k))) 
            ! \beta(k+1)/\beta(k) 
            cff2=dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k+1)
     &                  *dtfast+visc2v_nbq(i,j,k+1))
     &                *cff1
#   endif
            !-------------------------
            ! b*(k)=b(k)-a(k)*c(k-1)
            !-------------------------        
            cff=cff1+Hzw_nbq_inv(i,j,k)       ! b(k)
     &              *(
     &                 Hzr_nbq_inv(i,j,k)
     &                +Hzr_nbq_inv(i,j,k+1)  
#   ifdef M3FAST_CSVISC2K
     &                *cff2  
#   endif
     &              )
            !-------------------------
            ! c*(k)=c(k)/b*(k)
            !-------------------------
            FC3D(i,j,k)=-Hzw_nbq_inv(i,j,k+1)
     &              *Hzr_nbq_inv(i,j,k+1)/cff
#   ifdef M3FAST_CSVISC2K
     &                *cff2  
#   endif
            !-------------------------
            ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
            !-------------------------
            DC3D(i,j,k)=qdmw_nbq(i,j,k)*cff1/cff  ! RHS=f(k)
     &              +Hzw_nbq_inv(i,j,k-1)         ! -ak(k)
     &              *Hzr_nbq_inv(i,j,k)
     &              *qdmw_nbq(i,j,k-1)/cff
              
          enddo
        enddo
! !
! !********************************
! !..........Inner layers
! !********************************
! !
! !--------------------------------
! ! Sediment Layer (implicit)
! !--------------------------------
! !
#   ifdef M3FAST_SEDLAYERS
        do k=-N_sl+2,0
! ! !$acc loop independent private( cff1, cff )        
!$acc loop independent collapse(2)
          do j=Jstr,Jend  !<-- j loop
          do i=Istr,Iend
            !-------------------------
            ! \beta(k) (normalization)
            !-------------------------  
#   ifndef M3FAST_CSVISC2K
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtfast+visc2v_nbq)) 
#   else
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k)
     &              *dtfast+visc2v_nbq(i,j,k))) 
            ! \beta(k+1)/\beta(k) 
            cff2=dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k+1)
     &                  *dtfast+visc2v_nbq(i,j,k+1))
     &                *cff1
#   endif
            !-------------------------
            ! b*(k)=b(k)-a(k)*c*(k-1)
            !-------------------------        
            cff=cff1+Hzw_nbq_inv(i,j,k)       ! b(k)
     &              *(
     &                Hzr_nbq_inv(i,j,k)
     &                +Hzr_nbq_inv(i,j,k+1)
#   ifdef M3FAST_CSVISC2K
     &                *cff2       
#   endif
     &              )
     &              + Hzw_nbq_inv(i,j,k-1)    ! -a(k)
     &               *Hzr_nbq_inv(i,j,k)
     &               *FC3D(i,j,k-1)           !c*(k-1)
            !-------------------------
            ! c*(k)=c(k)/b*(k)
            !-------------------------
            FC3D(i,j,k)=-Hzw_nbq_inv(i,j,k+1)
     &              *Hzr_nbq_inv(i,j,k+1)/cff
#   ifdef M3FAST_CSVISC2K
     &              *cff2
#   endif
            !-------------------------
            ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
            !-------------------------
            DC3D(i,j,k)=(qdmw_nbq(i,j,k)*cff1        ! RHS=f(k)
     &              +Hzw_nbq_inv(i,j,k-1)            ! -ak(k)
     &              *Hzr_nbq_inv(i,j,k)
!    &              *DC(i,k-1)) /cff           
     &              *DC3D(i,j,k-1)) /cff           
          enddo
          enddo
        enddo
#   endif
! !
! !--------------------------------
! ! Ocean inner Layer (implicit)
! !--------------------------------
! ! 
#   ifndef M3FAST_SEDLAYERS
        do k=2,N-1
#   else
        do k=1,N-1
#   endif
! ! !$acc loop independent private( cff1, cff )        
!$acc loop independent collapse(2)
          do j=Jstr,Jend  !<-- j loop
          do i=Istr,Iend
            !-------------------------
            ! \beta(k) (normalization)
            !-------------------------  
#   ifndef M3FAST_CSVISC2K
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtfast+visc2v_nbq)) 
#   else
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k)
     &              *dtfast+visc2v_nbq(i,j,k))) 
            ! \beta(k+1)/\beta(k) 
            cff2=dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k+1)
     &                  *dtfast+visc2v_nbq(i,j,k+1))
     &                *cff1
#   endif
            !-------------------------
            ! b*(k)=b(k)-a(k)*c(k-1)
            !-------------------------        
            cff=cff1+Hzw_nbq_inv(i,j,k)       ! b(k)
     &              *(
     &                Hzr_nbq_inv(i,j,k)
     &                +Hzr_nbq_inv(i,j,k+1)
#   ifdef M3FAST_CSVISC2K
     &                  *cff2       
#   endif
     &              )
     &              + Hzw_nbq_inv(i,j,k-1)    ! -a(k)
     &               *Hzr_nbq_inv(i,j,k)
     &               *FC3D(i,j,k-1)           !c(k-1)
            !-------------------------
            ! c*(k)=c(k)/b*(k)
            !-------------------------
            FC3D(i,j,k)=-Hzw_nbq_inv(i,j,k+1)
     &              *Hzr_nbq_inv(i,j,k+1)/cff
#   ifdef M3FAST_CSVISC2K
     &              *cff2
#   endif
            !-------------------------
            ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
            !-------------------------
            DC3D(i,j,k)=(qdmw_nbq(i,j,k)*cff1        ! RHS=f(k)
     &              +Hzw_nbq_inv(i,j,k-1)            ! -ak(k)
     &              *Hzr_nbq_inv(i,j,k)
!    &              *DC(i,k-1)) /cff           
     &              *DC3D(i,j,k-1)) /cff           
          enddo
          enddo
        enddo
! !
! !--------------------------------
! !..........Surface BC
! !--------------------------------
! !
        k=N
        do j=Jstr,Jend  !<-- j loop
        do i=Istr,Iend
            !-------------------------
            ! \beta (normalization)
            !-------------------------     
#   ifndef M3FAST_CSVISC2K
          cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &                    *dtfast+visc2v_nbq)) 
#   else
          cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,N)
     &                    *dtfast+visc2v_nbq(i,j,k))) 
#   endif
            !-------------------------
            ! b*(k)=b(k)-a(k)*c(k-1)
            !-------------------------        
          cff=cff1+Hzw_nbq_inv(i,j,N)
     &            *Hzr_nbq_inv(i,j,N) 
     &       +Hzw_nbq_inv(i,j,N-1)
     &             *Hzr_nbq_inv(i,j,N)
     &             *FC3D(i,j,N-1) 
            !-------------------------
            ! c*(k)=c(k)/b*(k)
            !-------------------------
          FC3D(i,j,N)=0. 
            !-------------------------
            ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
            !-------------------------
            ! WSURF
          DC3D(i,j,N)=(qdmw_nbq(i,j,N)*cff1
     &             +Hzw_nbq_inv(i,j,N-1)
     &             *Hzr_nbq_inv(i,j,N)  
     &             *DC3D(i,j,N-1))/cff
          qdmw_nbq(i,j,N)=DC3D(i,j,N) 
        enddo
        enddo
! !
! !********************************
! !.....Solves tri-diag system
! ! (starts above with qdmw_nbq(N))
! !********************************
! !
!         k=N-1
!         do i=Istr,Iend
!             qdmw_nbq(i,j,k)=DC(i,k)/(1.+CF(i,k))
!             qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N-1)
!          enddo
! 
!           do i=Istr,Iend
!             qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
!           enddo            
#   ifndef M3FAST_SEDLAYERS
        do k=N-1,1,-1
#   else
        do k=N-1,-N_sl+1,-1
#   endif
          do j=Jstr,Jend 
          do i=Istr,Iend
            qdmw_nbq(i,j,k)=DC3D(i,j,k)-FC3D(i,j,k)*qdmw_nbq(i,j,k+1)
          enddo            
          enddo                        
        enddo 
!
!$acc end kernels
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_qdmw_implicit.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
