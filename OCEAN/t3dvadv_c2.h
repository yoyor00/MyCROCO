!
!----------------------------------------------------------
! Compute vertical advective fluxes
! using 2th-order centered scheme
!----------------------------------------------------------
!

         do k=1,N-1
           do i=Istr,Iend
             FC(i,k)=0.5*We(i,j,k)*(t(i,j,k  ,nadv,itrc)
     &                           +  t(i,j,k+1,nadv,itrc))
           enddo
         enddo
         do i=Istr,Iend
            FC(i,0)=0.
           FC(i,N )=0.
           CF(i,0 )=dt*pm(i,j)*pn(i,j)
         enddo