!
!-------------------------------------------------------------------
!  Solve implicit W-momentum equation
!-------------------------------------------------------------------
!
#  ifdef NBQ_IMP
      do k=0,N
        do j=Jstr,Jend
          do i=Istr,Iend
            qdmw_nbq_old(i,j,k)=qdmw_nbq(i,j,k)
          enddo
        enddo
      enddo
  
      do j=Jstr,Jend ! <-- j loop

        do k=1,N
          do i=Istr,Iend
            FC(i,k)= soundspeed2_nbq*rho_nbq(i,j,k)
     &              -(thetaimp_nbq*soundspeed2_nbq*dtnbq+visc2_nbq)
     &                       *(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k))
    
#   ifdef NBQ_THETAIMP 
            FC(i,k)=FC(i,k)
     &       -thetaimp_nbq*(1.-thetaimp_nbq)*soundspeed2_nbq*dtnbq*
     &                   ( Hzw_nbq_inv(i,j,k  )*qdmw_nbq(i,j,k)
     &                    -Hzw_nbq_inv(i,j,k-1)*qdmw_nbq(i,j,k-1) )
#   endif

            FC(i,k)=FC(i,k)*Hzr_nbq_inv(i,j,k) 
          enddo
        enddo
     
!.........Inner layers

        do k=1,N-1
          do i=Istr,Iend   
            dum_s = FC(i,k) - FC(i,k+1)
#   ifdef UV_COR_NT
            dum_s = dum_s + ntcorw(i,j,k)
#   endif
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &                        + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
#   ifdef MASKING
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#   endif               
          enddo             
        enddo

!.........Surface BC

        k=N
        do i=Istr,Iend
          dum_s =   FC(i,k)
#   ifdef UV_COR_NT
          dum_s = dum_s + ntcorw(i,j,k)
#   endif
          qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &                      + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
#   ifdef MASKING
          qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#   endif              
        enddo   

      enddo !<-- j loop
# endif 
