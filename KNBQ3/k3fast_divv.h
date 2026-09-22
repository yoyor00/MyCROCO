! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_divv.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
# if defined CVTK_DEBUG_ADV1 && defined KNBQ
!$acc wait
      call check_tab3d(thetadiv_nbq,'3d_fast bef divv thetadiv_nbq',
     &  'r',ondevice=.TRUE.)
      call check_tab3d(Hzw_nbq_inv,'3d_fast bef divv Hzw_nbq_inv',
     &  'r',ondevice=.TRUE.)
      call check_tab3d(qdmw_nbq,'3d_fast bef divv qdmw_nbq',
     &  'r',ondevice=.TRUE.)
# endif    
!$acc kernels if(compute_on_device) default(present) async(1)    
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
     &           +(1.-thetaimp_nbq)*qdmwold_nbq(i,j,k)) 
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
     &              +(1.-thetaimp_nbq)*qdmwold_nbq(i,j,k))
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
# if defined CVTK_DEBUG_ADV1 && defined KNBQ
!$acc wait
      call check_tab3d(thetadiv_nbq,'3d_fast aft divv thetadiv_nbq',
     &  'r',ondevice=.TRUE.)
# endif    
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_divv.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !           
