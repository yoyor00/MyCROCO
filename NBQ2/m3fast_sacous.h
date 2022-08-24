#    if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(rho_nbq,'step3d_fastrho_nbq(a)','r')
      call check_tab3d(z_w,'z_w (a)','r')
      call check_tab2d(xr,'xr (a)','r')
      call check_tab2d(yr,'yr (a)','r')
C$OMP END MASTER
#    endif   
!$acc kernels default( present )
        if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do j=jvar1,jvar2
      ! do i=IstrU-2,Iend+1 ! Patch temporaire pour sorties (Francis)
        do i=IstrU-1,Iend
!       cff=0.
!       cff1=0.

        do k=-N_sl+1,N
              dist_d=sqrt((xr(i,j)-xl/2.)**2+(0.*(yr(i,j)-el/2.))**2
     &                             +(abs(z_w(i,j,k))-hmax_exp/2.)**2)

             rho_nbq(i,j,k)=rho_nbq(i,j,k)
     &               +amp_exp*sin(2*acos(-1.)*time_nbq/period_exp)
     &                       *exp(-dist_d**2/for_a_exp**2)
     &                       *Hzr(i,j,k)
!            cffz(k)=+amp_exp
!    &                       *sin(2*acos(-1.)*time_nbq/period_exp)
!    &                       *exp(-dist_d**2/for_a_exp**2)
!    &                       *Hzr(i,j,k)
!           cff=cff+cffz(k)*Hzr(i,j,k)
!           cff1=cff1+Hzr(i,j,k)
        enddo
!       do k=1,N
!            rho_nbq(i,j,k)=rho_nbq(i,j,k)+cffz(k)-cff/cff1
!       enddo
        enddo
        enddo
!$acc end kernels     