!
        if (FIRST_TIME_STEP) then
          cff3=0.                        ! This version is designed
          cff2=0.                        ! for coupling during 3D
          cff1=1.                        ! predictor sub-step: here
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
          cff3= 0.                       ! computed as instantaneous
          cff2=-0.5                      ! value at 3D time step
          cff1= 1.5                      ! "nstp" first, and then
        else                             ! extrapolated half-step
          cff3= 0.281105                 ! forward using  AM3-like
          cff2=-0.5-2.*cff3              ! weights optimized for
          cff1= 1.5+cff3                 ! maximum stability (with
        endif                            ! special care for startup)
        

! ! KERNEL_7  cff <= ( rubar )
! ! KERNEL_7  ru_int2d_nbq <= ( cff , ru_int2d_nbq_bak )
! ! KERNEL_7  ru_int2d_nbq_bak <= ( cff )
! ! KERNEL_7  cff <= ( rvbar )
! ! KERNEL_7  rv_int2d_nbq <= ( cff , rv_int2d_nbq_bak )
! ! KERNEL_7  rv_int2d_nbq_bak <= ( cff )

!$acc kernels default(present)
!# ifndef M3FAST_ZETAW
# if defined M3FAST_COUPLING2D  || defined NHINT
        do j=Jstr,Jend
         do i=IstrU,Iend
#  ifdef M3FAST_BOTH
           cff=rufrc(i,j)-rubar(i,j)-rubarh(i,j)
#  else           
           cff=rufrc(i,j)-rubar(i,j)
#  endif           
           rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                         + cff3*rufrc_bak(i,j,nstp)
           rufrc_bak(i,j,nstp)=cff
         enddo
        enddo
# endif
# if defined M3FAST_C3D_UVSF && defined M3FAST_COUPLING3D
        do j=Jstr,Jend
         do i=IstrU,Iend
           cff=rubar(i,j)! -rufrc(i,j)
# ifndef BODYFORCE
!#  ifdef BSTRESS_FAST
!     &        - sustr(i,j)*om_u(i,j)*on_u(i,j)
!#  else
!     &        -(sustr(i,j)-bustr(i,j))*om_u(i,j)*on_u(i,j)
!#  endif
# endif
           ru_int2d_nbq(i,j)=cff1*cff+cff2*ru_int2d_nbq_bak(i,j,3-nstp)
     &                               +cff3*ru_int2d_nbq_bak(i,j,nstp)
           ru_int2d_nbq_bak(i,j,nstp)=cff
         enddo
        enddo
# endif      
!# ifndef M3FAST_ZETAW
# if defined M3FAST_COUPLING2D  || defined NHINT
        do j=JstrV,Jend
         do i=Istr,Iend
#  ifdef M3FAST_BOTH  
           cff=rvfrc(i,j)-rvbar(i,j)-rvbarh(i,j)
#  else           
           cff=rvfrc(i,j)-rvbar(i,j)
#  endif           
           rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                         + cff3*rvfrc_bak(i,j,nstp)
           rvfrc_bak(i,j,nstp)=cff
         enddo
        enddo
# endif
# if defined M3FAST_C3D_UVSF &&  defined M3FAST_COUPLING3D
        do j=JstrV,Jend
         do i=Istr,Iend
           cff=rvbar(i,j)!-rvfrc(i,j)
# ifndef BODYFORCE
!#  ifdef BSTRESS_FAST
!     &        - svstr(i,j)*om_v(i,j)*on_v(i,j)
!#  else
!     &        -(svstr(i,j)-bvstr(i,j))*om_v(i,j)*on_v(i,j)
!#  endif
# endif
           rv_int2d_nbq(i,j)=cff1*cff+cff2*rv_int2d_nbq_bak(i,j,3-nstp)
     &                               +cff3*rv_int2d_nbq_bak(i,j,nstp)
           rv_int2d_nbq_bak(i,j,nstp)=cff
         enddo
        enddo
# endif        
#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(rufrc_bak(:,:,nstp),'rufrc_bak_st_fast_d','u')
         call check_tab2d(rufrc_bak(:,:,3-nstp),'rufrc_bak_d','u')
         call check_tab2d(rufrc(:,:),'rufrc st_fast_d','u')
         call check_tab2d(rvfrc(:,:),'rvfrc st_fast_d','v')
#endif
