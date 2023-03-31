! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_divv.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
!$acc kernels default( present )    
! ! 
      do j=Jstr,Jend  !<-- j loop
! ! *******************************
! ! Vertical flux in bottom layer
! ! *******************************
      k = -N_sl
#  ifdef NBQ_FREESLIP
        do i=Istr,Iend
#   ifdef NBQ_THETAIMP
          FC3D(i,j,k)=Hzw_nbq_inv(i,j,k)*
     &            (thetaimp_nbq*qdmw_nbq(i,j,k)
     &           +(1.-thetaimp_nbq)*rw_nbq(i,j,k)) 
#   else
          FC3D(i,j,k)=Hzw_nbq_inv(i,j,k)*qdmw_nbq(i,j,k)       
#   endif
        enddo
#  else /*  ! NBQ_FREESLIP */
        do i=Istr,Iend
          FC3D(i,j,k)=0.   ! Bottom BC
        enddo
#  endif /* NBQ_FREESLIP */
! ! *******************************
! ! Vertical flux in the water column
! ! Sum divergence
! ! *******************************
        do k=-N_sl+1,N
          do i=Istr,Iend
#  ifdef NBQ_THETAIMP
            FC3D(i,j,k)=Hzw_nbq_inv(i,j,k)*
     &              (thetaimp_nbq*qdmw_nbq(i,j,k)
     &              +(1.-thetaimp_nbq)*rw_nbq(i,j,k))
#  else
            FC3D(i,j,k)=Hzw_nbq_inv(i,j,k)*qdmw_nbq(i,j,k)
#  endif   
            thetadiv_nbq (i,j,k)= thetadiv_nbq (i,j,k)
     &                           +FC3D(i,j,k)-FC3D(i,j,k-1) 
#  ifdef NBQ_HZ_PROGNOSTIC
            thetadiv2_nbq(i,j,k)= thetadiv2_nbq(i,j,k)
     &                           +FC3D(i,j,k)-FC3D(i,j,k-1)  
#  endif            
          enddo
        enddo
      enddo  !<-- j loop
!$acc end kernels
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_divv.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !           
