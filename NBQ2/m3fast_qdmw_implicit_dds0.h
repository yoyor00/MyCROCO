#  ifdef NBQ_IMP
!---------------------
! Gaussian Elimination
!---------------------

!.......Compute coef.

      do j=Jstr,Jend  !<-- j loop

!..........Bottom BC
        k=1
        do i=Istr,Iend
          cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq
     &            *dtnbq+visc2_nbq)) 
          cff=(cff1+Hzw_nbq_inv(i,j,1)*(Hzr_nbq_inv(i,j,1)
     &             +Hzr_nbq_inv(i,j,2)))
          CF(i,1)=(-Hzw_nbq_inv(i,j,2)*Hzr_nbq_inv(i,j,2))
     &                                                          /cff
          DC(i,1)=qdmw_nbq(i,j,1)*cff1/cff   
     &           +qdmw_nbq(i,j,0)/cff*Hzw_nbq_inv(i,j,0)
     &                               *Hzr_nbq_inv(i,j,1)
        enddo

!..........Inner layers
        do k=2,N-1
          do i=Istr,Iend
            cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtnbq+visc2_nbq)) 
            cff=(cff1+
     &               Hzw_nbq_inv(i,j,k)*(Hzr_nbq_inv(i,j,k)
     &              +Hzr_nbq_inv(i,j,k+1))   
     &              +Hzw_nbq_inv(i,j,k-1)*Hzr_nbq_inv(i,j,k)
     &              *CF(i,k-1))
            CF(i,k)=(-Hzw_nbq_inv(i,j,k+1)
     &               *Hzr_nbq_inv(i,j,k+1))/cff
            DC(i,k)=(qdmw_nbq(i,j,k)*cff1+Hzw_nbq_inv(i,j,k-1)
     &               *Hzr_nbq_inv(i,j,k)*DC(i,k-1)) /cff           
          enddo            
        enddo

!..........Surface BC
        k=N
        do i=Istr,Iend
          cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq
     &            *dtnbq+visc2_nbq)) 
          cff=(cff1+Hzw_nbq_inv(i,j,N)*Hzr_nbq_inv(i,j,N) 
     &             +Hzw_nbq_inv(i,j,N-1)
     &             *Hzr_nbq_inv(i,j,N)*CF(i,N-1))  
          CF(i,N)=0. 
          DC(i,k)=(qdmw_nbq(i,j,N)*cff1+Hzw_nbq_inv(i,j,N-1)
     &             *Hzr_nbq_inv(i,j,N)*DC(i,N-1))/cff
        enddo 

!..........Solves tri-diag system
        do i=Istr,Iend
          qdmw_nbq(i,j,N)=DC(i,k)   
#   ifdef NBQ_GRAV
     &                    -rho_nbq(i,j,N)*0.5*g*dtnbq
#   endif     
        enddo

        do k=N-1,1,-1
          do i=Istr,Iend
            qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
#   ifdef NBQ_GRAV
     &            -0.25*(rho_nbq(i,j,k)*Hzr_nbq_inv(i,j,k)
     &                 +rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1))
     &                *(Hzr(i,j,k)+Hzr(i,j,k+1))
     &                 *g*dtnbq
#   endif
          enddo            
        enddo                        

      enddo   !<-- j loop 
!
!-------------------------------------------------------------------
!  W-momentum open boundary conditions
!-------------------------------------------------------------------
!
#   ifdef OBC_NBQ
      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
#   endif

#   ifdef UV_COR_NT
#    if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmw_nbq(START_2D_ARRAY,0))
#    endif
#   endif

#  endif /* NBQ_IMP */
