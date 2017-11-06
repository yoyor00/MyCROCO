 
# if defined NBQ 
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
          enddo
        enddo

#   ifdef MASKING
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)*umask(i,j)
            enddo
          enddo
        enddo
#   endif
 
        do j=JstrV,Jend
          do i=Istr,Iend
            rvext_nbq_2d_sum(i,j)=0.
            rvext_nbq_2d_old(i,j)=0.
          enddo
        enddo 

#   ifdef MASKING
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)*vmask(i,j)
            enddo
          enddo
        enddo
#   endif

      endif
#  endif
!
!------------------------------
! All 2D time-steps 
!------------------------------
!

!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!     call exchange_u3d_tile (Istr,Iend,Jstr,Jend,         ! TO be removed
!    &                   ru_int_nbq(START_2D_ARRAY,1))
!     call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
!    &                   rv_int_nbq(START_2D_ARRAY,1))
!#  endif  

#  define ruext_nbq_2d UFx
      do j=Jstr,Jend
        do i=IstrU,Iend
          ruext_nbq_2d(i,j)=(rufrc(i,j)+rubar(i,j))*pm_u(i,j)*pn_u(i,j)
     &
     &              /(Drhs(i,j)+Drhs(i-1,j))
#  ifdef MASKING
          ruext_nbq_2d(i,j)=ruext_nbq_2d(i,j)*umask(i,j)
#  endif
          ruext_nbq_2d_old(i,j)=ruext_nbq_2d(i,j)-ruext_nbq_2d_old(i,j)
#  if defined M2FILTER_NONE     
          ruext_nbq_2d_sum(i,j)=ruext_nbq_2d_sum(i,j)+ruext_nbq_2d(i,j)
#  endif          
        enddo
      enddo
      
!#  ifdef RVTK_DEBUG
!      call check_tab3d(ru_int_nbq,'ru_int_nbq','u')
!#  endif  

      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)= ru_int_nbq(i,j,k)+
     &        ruext_nbq_2d_old(i,j)
     &        *(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
     
      do j=Jstr,Jend
        do i=IstrU,Iend
          ruext_nbq_2d_old(i,j)=ruext_nbq_2d(i,j)
        enddo
      enddo
    
#  undef ruext_nbq_2d

#  define rvext_nbq_2d UFx
      do j=JstrV,Jend
        do i=Istr,Iend
          rvext_nbq_2d(i,j)=(rvfrc(i,j)+rvbar(i,j))*pm_v(i,j)*pn_v(i,j)
     &              /(Drhs(i,j)+Drhs(i,j-1))
#  ifdef MASKING
          rvext_nbq_2d(i,j)=rvext_nbq_2d(i,j)*vmask(i,j)
#  endif
          rvext_nbq_2d_old(i,j)=rvext_nbq_2d(i,j)-rvext_nbq_2d_old(i,j)
#  if defined M2FILTER_NONE     
          rvext_nbq_2d_sum(i,j)=rvext_nbq_2d_sum(i,j)+rvext_nbq_2d(i,j)
#  endif          
        enddo
      enddo
        
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)+
     &          rvext_nbq_2d_old(i,j)
     &         *(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          rvext_nbq_2d_old(i,j)=rvext_nbq_2d(i,j)
        enddo
      enddo
#  undef rvext_nbq_2d 
# endif

# ifdef NBQ
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
c LAURENT: Is this useful ? the RHS do not have to be valid in ghost cells (TO BE TESTED with RVTK)
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                                 ru_int_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                                 rv_int_nbq(START_2D_ARRAY,1))
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
!     &                                 rw_int_nbq(START_2D_ARRAY,0))
#  endif     
# endif 

# ifdef RVTK_DEBUG
      call check_tab3d(ru_int_nbq,'ru_int_nbq (A)','u')
      call check_tab3d(rv_int_nbq,'rv_int_nbq (A)','v')
      call check_tab3d(rw_int_nbq,'rw_int_nbq (A)','wint')
# endif  
