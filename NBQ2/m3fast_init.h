! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_init.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! 
! !
! !********************************
! !  Store qdm(u,v,w).nbq
! !********************************
! !
# ifdef M3FAST_C3D_UVFS
      if (FIRST_FAST_STEP) then
!$acc kernels default(present)
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_nbq_avg2(i,j,k)=qdmu_nbq(i,j,k)
            enddo
          enddo 
        enddo
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_nbq_avg2(i,j,k)=qdmv_nbq(i,j,k)
            enddo
          enddo 
        enddo
!$acc end kernels
      endif
# endif      
# ifdef M3FAST_C3D_WFS
      if (FIRST_FAST_STEP) then
!$acc update device( qdmw_nbq )   !! iif=1  
!$acc kernels default(present)
          do k=0,N
           do j=Jstr,Jend
            do i=Istr,Iend
#  ifdef NHINT_WH
              wzh_nbq(i,j,k)=0.
#  endif
              rw_nbq_avg2(i,j,k)=qdmw_nbq(i,j,k)
            enddo
          enddo 
        enddo
!$acc end kernels
      endif    ! FIRST_FAST_STEP
# endif
! !
! !********************************
! !  zw_nbq backups
! !********************************
! !
# ifdef M3SLOW_W
      if (FIRST_TIME_STEP .and. FIRST_FAST_STEP) then
!$acc update device( z_w )   !! iif=1 ic=1
!$acc kernels default(present)
        do k=-N_sl,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              zw_nbq(i,j,k)=z_w(i,j,k)
            enddo
          enddo
        enddo
!$acc end kernels
      endif
# endif  /* M3SLOW_W */
! !
! !********************************
! !  Implicit part: system setup
! !********************************
! !    
# ifdef M3SLOW_W
!$acc kernels default( present )
      do j=Jstr,Jend
        do i=Istr,Iend
          work(i,j)=pm(i,j)*pn(i,j)
#  ifdef M3FAST_ZETAW
          DU_nbq(i,j)=0.
          DV_nbq(i,j)=0.
#  endif
        enddo
      enddo
!$acc end kernels
# endif /* M3SLOW_W */
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_init.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! 
