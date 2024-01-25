! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_sacous.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
#    if defined RVTK_DEBUG && defined KNBQ
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(rho_nbq,'step3d_fastrho_nbq(a)','r')
      call check_tab3d(z_w,'z_w (a)','r')
      call check_tab2d(xr,'xr (a)','r')
      call check_tab2d(yr,'yr (a)','r')
C$OMP END MASTER
#    endif   
!$acc kernels if(compute_on_device) default(present)
        time_nbq = time_nbq + dtfast
        do k=-N_sl+1,N
        do j=Jstr,Jend
        do i=Istr,Iend
              dist_d=sqrt((xr(i,j)-xl/rx_exp)**2+(yr(i,j)-el/ry_exp)**2
     &                             +(abs(z_w(i,j,k))-hmax_exp/rz_exp)**2)

             rho_nbq(i,j,k)=rho_nbq(i,j,k)
     &               +amp_exp*sin(2*acos(-1.)*time_nbq/period_exp)
     &                       *exp(-dist_d**2/for_a_exp**2)
     &                       *Hzr(i,j,k)
          enddo
         enddo
        enddo
!$acc end kernels     
! !
! !--------------------------------------------------------------------
! !         Exchange rho_nbq
! !--------------------------------------------------------------------
! !     
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq(-2,-2,1))


! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_sacous.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
