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

      if (sacous_nbq.lt.0) then  
! !
! !********************************
! !        User (ana_initial)
! !********************************
! !
        if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do j=jvar1,jvar2
        do i=IstrU-1,Iend
        do k=-N_sl+1,N
             dist_d=sqrt((xr(i,j)-xl/rx_exp)**2+(yr(i,j)-el/ry_exp)**2
     &                             +(abs(z_w(i,j,k))-hmax_exp/rz_exp)**2)

             rho_nbq(i,j,k)=rho_nbq(i,j,k)
     &               +amp_exp*sin(2*acos(-1.)*time_nbq/period_exp)
     &                       *exp(-dist_d**2/for_a_exp**2)
     &                       *Hzr(i,j,k)
          enddo
         enddo
        enddo
      elseif (sacous_nbq==0) then  
! !
! !********************************
! !         Ray
! !********************************
! !
         if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do j=jvar1,jvar2
        do i=IstrU-1,Iend
        do k=-N_sl+1,N
            dist_d=sqrt((xr(i,j)-xl*2./100.)**2+(0.*(yr(i,j)-el/2.))**2
     &                             +(abs(z_w(i,j,k))-hmax_exp*1000./5000.)**2)
            if (time_nbq.ge.0.and.time_nbq.le.1000) then
             rho_nbq(i,j,k)=rho_nbq(i,j,k)
     &               +amp_exp*sin(2*acos(-1.)*time_nbq/period_exp)
     &                       *exp(-dist_d**2/for_a_exp**2)
     &                       *Hzr(i,j,k)
	    endif
        enddo
        enddo
        enddo  
        
       elseif (sacous_nbq==1) then       
! !
! !********************************
! !          Mode 1
! !********************************
! !
        if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do k=1,N  ! Changer ici pour source dans Sdl
          do j=jvar1,jvar2
            do i=IstrU-2,Iend+1
              dist_d=sqrt((xr(i,j)-xl*2./10.)**2)
              if (time_nbq.le.1.) then
               if (dist_d.le.2.) then
              rho_nbq(i,j,k)=rho_nbq(i,j,k)+amp_exp*
     &                         (-9.81*.14429**2/(70*3.14)**2*0.02595**2*
     &                         *cos(0.02595*z_w(i,j,k))+1/0.02595*
     &                         sin(0.02595*z_w(i,j,k)))
     &                         *sin(2*acos(-1.)*time_nbq/period_exp)
     &                         *exp(-dist_d**2/2.**2)
               endif
              endif
            enddo
          enddo
        enddo
       elseif (sacous_nbq==2) then
! !
! !********************************
! !          Mode 2
! !********************************
! !
        if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do k=1,N  ! Changer ici pour source dans Sdl
          do j=jvar1,jvar2
            do i=IstrU-2,Iend+1
              dist_d=sqrt((xr(i,j)-xl*2./10.)**2)
              if (time_nbq.le.1.) then
               if (dist_d.le.2.) then
              rho_nbq(i,j,k)=rho_nbq(i,j,k)+amp_exp*
     &                         (-9.81*.13675**2/((70*3.14)**2*0.05286**2)
     &                         *cos(0.05286*z_w(i,j,k))+1/0.05286*
     &                         sin(0.05286*z_w(i,j,k)))
     &                         *sin(2*acos(-1.)*time_nbq/period_exp)
     &                         *exp(-dist_d**2/2.**2)
               endif
              endif
            enddo
          enddo
        enddo
       elseif (sacous_nbq==3) then 
! !
! !********************************
! !          Mode 3
! !********************************
! !
         if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do k=1,N  ! Changer ici pour source dans Sdl
          do j=jvar1,jvar2
            do i=IstrU-2,Iend+1
              dist_d=sqrt((xr(i,j)-xl*2./10.)**2)
              if (time_nbq.le.1.) then
               if (dist_d.le.2.) then
              rho_nbq(i,j,k)=rho_nbq(i,j,k)+amp_exp*
     &                         (-9.81*.12316**2/((70*3.14)**2*0.07953**2)
     &                         *cos(0.07953*z_w(i,j,k))+1/0.07953*
     &                         sin(0.07953*z_w(i,j,k)))
     &                         *sin(2*acos(-1.)*time_nbq/period_exp)
     &                         *exp(-dist_d**2/2.**2)
               endif
              endif
            enddo
          enddo
        enddo     
       elseif (sacous_nbq==10) then
! !
! !********************************
! !       Wedge 3D type source
! !********************************
! !         
        if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
        time_nbq = time_nbq + dtfast
        do k=-N_sl+1,N
        do j=jvar1,jvar2
        do i=IstrU-2,Iend+1
            dist_d=sqrt((xr(i,j)-1000.)**2+((yr(i,j)-1000.))**2+
     &              (abs(z_w(i,j,k))-40)**2)
            if (time_nbq.le.15) then
             rho_nbq(i,j,k)=rho_nbq(i,j,k)
     &               +amp_exp*0.5*(1-cos(2*acos(-1.)*time_nbq/period_exp/4))
     &                       *sin(2*acos(-1.)*time_nbq/period_exp)
     &                       *exp(-dist_d**2/for_a_exp**2)
     &                       *Hzr(i,j,k)
	    endif
        enddo
        enddo
        enddo
       endif    
!$acc end kernels     
