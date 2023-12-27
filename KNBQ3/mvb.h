! !
! !***********************************************************************
! !
! !              Construction of an analytical forcing for
! !                      a moving bathy. test-case
! !
! !         These lines are called by ana_initial & analytical
! !         subroutines, new variables are to be declared in these
! !                           subroutines.
! !                               
! !***********************************************************************
! !
! ! 
! !====================================================================
! !                      Increment MVB time-step
! !====================================================================
! !
       i=kbak2
       kbak2=kstp2
       kstp2=knew2
       knew2=i
! ! 
! !====================================================================
! !                      Displacement field
! !====================================================================
! !
       cff=2.*pi/10.05
       do j=JR_RANGE
        do i=IR_RANGE
         !   CAUTION: displ. has to be in advance (m+2) -> kbak
!           if (abs(xr(i,j)).le.0.03686*2) then
            x_mvb(i,j,kbak2)=
     &                    min(1., max(0.,tanh((time_mvb+dtfast)/10.)))*(
     &                    +1.e-3*sin(cff*(time_mvb+dtfast))  ! oscillation
!    &                    +1.*time_mvb                       ! translation
     &                                                    )
!           else
!           x_mvb(i,j,kbak2)=0.
!           endif
!!          y_mvb(i,j,kbak2)=0.
        enddo
       enddo
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
!     &                        x_mvb(START_2D_ARRAY,kbak2))
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
!     &                        y_mvb(START_2D_ARRAY,kbak2))
#  endif
! ! 
! !====================================================================
! !                      Time-dependant bottom topography
! !====================================================================
! !
      cff=1./(0.03686)**2
      do j=JstrR,JendR
        do i=IstrR,IendR
          dh_mvb(i,j,kbak2)=
     &         0.394 - 0.1*exp(-cff*(xr(i,j)-x_mvb(i,j,kbak2))**2)
        enddo
      enddo
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        dh_mvb(START_2D_ARRAY,kbak2))
#  endif
#  ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(dh_mvb(:,:,knew2),'dh_mvb(knew2) mvb.h','r')
        call check_tab2d(dh_mvb(:,:,kbak2),'dh_mvb(kbak2) mvb.h','r')
C$OMP END MASTER
#  endif   
! ! 
! !====================================================================
! !                           Horizontal velocity
! !====================================================================
! !
      do j=JR_RANGE
        do i=IU_RANGE
          ! u at m+1
! LF scheme:
!           u_mvb(i,j,knew2)=  
!     &                      0.25*( (x_mvb(i  ,j,kbak2)-x_mvb(i  ,j,kstp2))    !FA
!     &                           +(x_mvb(i-1,j,kbak2)-x_mvb(i-1,j,kstp2))
!     &                          ) /dtfast
! FW scheme:
           u_mvb(i,j,knew2)=  
     &                      0.5*( (x_mvb(i  ,j,knew2)-x_mvb(i  ,j,kstp2))    !FA
     &                           +(x_mvb(i-1,j,knew2)-x_mvb(i-1,j,kstp2))
     &                          ) /dtfast
        enddo
      enddo
!      do j=JV_RANGE
!        do i=IR_RANGE
!          v_mvb(i,j,knew2)=0.
!        enddo
!      enddo
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        u_mvb(START_2D_ARRAY,knew2))
!      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
!     &                        v_mvb(START_2D_ARRAY,knew2))
#  endif
! ! 
! !====================================================================
! !                           Vertical velocity
! !====================================================================
! !
      do j=JstrR,Jend
        jm1=max(j-1,0)
        do i=IstrR,Iend
           im1=max(i-1,0)
! LF scheme:
!           w_mvb(i,j,knew2)=   
!     &                    -0.5*(dh_mvb(i,j,kbak2)-dh_mvb(i,j,kstp2))
!     &                      /dtfast
!     &                    -0.5*u_mvb(i  ,j,knew2)*pm_u(i  ,j)
!     &                    *(dh_mvb(i  ,j,knew2)-dh_mvb(im1,j,knew2))
!     &                    -0.5*u_mvb(i+1,j,knew2)*pm_u(i+1,j)
!     &                    *(dh_mvb(i+1,j,knew2)-dh_mvb(i  ,j,knew2))
!!     &                    -0.5*v_mvb(i,j  ,knew2)*pn_v(i,j  )
!!     &                    *(dh_mvb(i,j,knew2)-dh_mvb(i,jm1,knew2))
!!     &                    -0.5*v_mvb(i,j+1,knew2)*pn_v(i,j+1)
!!     &                    *(dh_mvb(i,j+1,knew2)-dh_mvb(i,j  ,knew2))
! FW scheme:
           w_mvb(i,j,knew2)=   
     &                    -(dh_mvb(i,j,knew2)-dh_mvb(i,j,kstp2))
     &                      /dtfast
     &                    -0.5*u_mvb(i  ,j,knew2)*pm_u(i  ,j)
     &                    *(dh_mvb(i  ,j,knew2)-dh_mvb(im1,j,knew2))
     &                    -0.5*u_mvb(i+1,j,knew2)*pm_u(i+1,j)
     &                    *(dh_mvb(i+1,j,knew2)-dh_mvb(i  ,j,knew2))
!!     &                    -0.5*v_mvb(i,j  ,knew2)*pn_v(i,j  )
!!     &                    *(dh_mvb(i,j,knew2)-dh_mvb(i,jm1,knew2))
!!     &                    -0.5*v_mvb(i,j+1,knew2)*pn_v(i,j+1)
!!     &                    *(dh_mvb(i,j+1,knew2)-dh_mvb(i,j  ,knew2))
        enddo
      enddo
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        w_mvb(START_2D_ARRAY,knew2))
#  endif
         
