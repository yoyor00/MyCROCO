! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_init.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! 
! !
! !********************************
! !  Store qdm(u,v,w).nbq
! !********************************
! !
# ifdef K3FAST_C3D_UVFS
      if (FIRST_FAST_STEP) then
!$acc kernels if(compute_on_device) default(present)

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
# ifdef K3FAST_C3D_WFS
      if (FIRST_FAST_STEP) then

!$acc kernels if(compute_on_device) default(present)

          do k=0,N
           do j=Jstr,Jend
            do i=Istr,Iend
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
# ifdef K3SLOW_W
      if (FIRST_TIME_STEP .and. FIRST_FAST_STEP) then

!$acc kernels if(compute_on_device) default(present)
        do k=-N_sl,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              zw_nbq(i,j,k)=z_w(i,j,k)
            enddo
          enddo
        enddo
!$acc end kernels
      endif
# endif  /* K3SLOW_W */
! !
! !********************************
! ! Update ubar and vbar
! !     and their BC (AGRIF)
! ! BC for DU(V)_nbq treated here
! !********************************
! !
!!$acc kernels if(compute_on_device) default(present)
!      do j=Jstr,Jend
!        do i=IstrU,Iend
!          ubar(i,j,knew)=urhs(i,j)
!        enddo
!      enddo
!      do j=JstrV,Jend
!        do i=Istr,Iend
!          vbar(i,j,knew)=vrhs(i,j)
!        enddo
!      enddo
!!$acc end kernels
      M2bc_nbq_flag=.true. ! apply boundary wet/dry conditions
                           ! and boundaries for DU(V)_nbq
      call u2dbc_tile   (Istr,Iend,Jstr,Jend, work)
      call v2dbc_tile   (Istr,Iend,Jstr,Jend, work)
!$acc wait	
! !
! !********************************
! !  Implicit part: system setup
! !********************************
! !    
# ifdef K3SLOW_W
!$acc kernels if(compute_on_device) default(present)
      do j=Jstr,Jend
        do i=Istr,Iend
          work(i,j)=pm(i,j)*pn(i,j)
#  if defined K3FAST_DUVNBQ || defined K3FAST_DUVNBQ2
#   ifdef KNHINT_3M	
       if ((iic==1.and.iif==1).or.mod(iif,nsdtnbq).eq.0) then
#   endif
            DU_nbq(i,j)=0.
            DV_nbq(i,j)=0.
#   ifdef KNHINT_3M	
       endif
#   endif
#  endif
        enddo
      enddo
!$acc end kernels
# endif /* K3SLOW_W */
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_init.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! 
