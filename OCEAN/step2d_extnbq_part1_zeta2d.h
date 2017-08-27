 
!
!-----------------------------------------------------------------------
!
! Prepares Coupling with NBQ
!
!   ru_int_nbq        : RHS (3D) ( *mask & 2D correction)
!
!   ruext_nbq_2d     : RHS (2D)
!   ruext_nbq_2d_old : RHS (2D) at previous time-step
!   ruext_nbq_2d_sum : to compute time-averaged RHS (2D)
!   
!
!-----------------------------------------------------------------------
!
# if defined M2FILTER_NONE
!
      if (FIRST_2D_STEP) then
!
!------------------------------
! First 2D time-step only
!------------------------------
!
        do j=Jstr,Jend
          do i=IstrU,Iend
            ruext_nbq_2d_sum(i,j)=0.
            ruext_nbq_2d_old(i,j)=0.
!            Drhs_half_u(i,j) = 0.
          enddo
        enddo

        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
!              Drhs_half_u(i,j)=Drhs_half_u(i,j)
!     &           *(Hz_half(i-1,j,k)+Hz_half(i,j,k))
#  ifdef MASKING
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)*umask(i,j)
#  endif
            enddo
          enddo
        enddo
 
        do j=JstrV,Jend
          do i=Istr,Iend
            rvext_nbq_2d_sum(i,j)=0.
            rvext_nbq_2d_old(i,j)=0.
!            Drhs_half_v(i,j) = 0.
          enddo
        enddo 

        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
!              Drhs_half_v(i,j)=Drhs_half_v(i,j)
!     &           *(Hz_half(i,j-1,k)+Hz_half(i,j,k))
#  ifdef MASKING
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)*vmask(i,j)
#  endif
            enddo
          enddo
        enddo

      endif
# endif
!
!------------------------------
! All 2D time-steps 
!------------------------------
!
# define ruext_nbq_2d UFx

      do j=Jstr,Jend
        do i=IstrU,Iend
          ruext_nbq_2d(i,j)=(rufrc(i,j)+rubar(i,j))*pm_u(i,j)*pn_u(i,j)
     &              /(Drhs(i,j)+Drhs(i-1,j))
# ifdef MASKING
          ruext_nbq_2d(i,j)=ruext_nbq_2d(i,j)*umask(i,j)
# endif
# if defined M2FILTER_NONE     
          ruext_nbq_2d_sum(i,j)=ruext_nbq_2d_sum(i,j)+ruext_nbq_2d(i,j)
# endif          
        enddo
      enddo
       
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)  +
     &        (ruext_nbq_2d(i,j)-ruext_nbq_2d_old(i,j))
     &        *(Hz_half(i-1,j,k)+Hz_half(i,j,k))
          enddo
        enddo
      enddo      

      do j=Jstr,Jend
        do i=IstrU,Iend
          ruext_nbq_2d_old(i,j)=ruext_nbq_2d(i,j)
        enddo
      enddo
      
# undef ruext_nbq_2d

# define rvext_nbq_2d UFx
      do j=JstrV,Jend
        do i=Istr,Iend
          rvext_nbq_2d(i,j)=(rvfrc(i,j)+rvbar(i,j))*pm_v(i,j)*pn_v(i,j)
     &              /(Drhs(i,j)+Drhs(i,j-1))
# ifdef MASKING
          rvext_nbq_2d(i,j)=rvext_nbq_2d(i,j)*vmask(i,j)
# endif
# if defined M2FILTER_NONE     
          rvext_nbq_2d_sum(i,j)=rvext_nbq_2d_sum(i,j)+rvext_nbq_2d(i,j)
# endif          
        enddo
      enddo
        
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)+
     &        (rvext_nbq_2d(i,j)-rvext_nbq_2d_old(i,j))
     &         *(Hz_half(i,j-1,k)+Hz_half(i,j,k))
          enddo
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          rvext_nbq_2d_old(i,j)=rvext_nbq_2d(i,j)
        enddo
      enddo
# undef rvext_nbq_2d 

!#ifdef NBQ
!# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!     call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  
!    &                                 ru_int_nbq(START_2D_ARRAY,1))
!     call exchange_u3d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  
!    &                                 rv_int_nbq(START_2D_ARRAY,1))
!     call exchange_u3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,  
!    &                                 rw_int_nbq(START_2D_ARRAY,0))
!# endif     
!#endif     
