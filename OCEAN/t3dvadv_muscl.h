!
!----------------------------------------------------------
! Compute vertical advective fluxes
! using 2th-order centered scheme
!----------------------------------------------------------
!
# if defined PREDICTOR
          rdtstp = (1.0-gamma)*dt
# else
          rdtstp = dt
# endif

         DO i=Istr,Iend
	    zdzt_kp1(i) = ( t(i,j,N,nadv,itrc) - t(i,j,N-1,nadv,itrc) )
            zslpz(i) = 0.0
         END DO

         DO k = N, 3, -1
            do i=Istr,Iend
            !  !-- Slopes of tracer
            !  ! masked vertical gradient at jk+2
            zdzt_kp2  = ( t(i,j,k-1,nadv,itrc) - t(i,j,k-2,nadv,itrc) )
            !               ! vertical slope at jk+1
            zslpz_kp1 = ( zdzt_kp1(i) + zdzt_kp2 )
     &        * (  0.25 + SIGN( 0.25, zdzt_kp1(i) * zdzt_kp2 )  )
            !  ! slope limitation
            zslpz_kp1 = SIGN( 1.0, zslpz_kp1 ) 
     &                 * MIN(    ABS( zslpz_kp1 ),
     &                 2.*ABS( zdzt_kp2 ),   
     &                 2.*ABS( zdzt_kp1(i) )  )
            !  !-- vertical advective flux at k+1
            z0w = SIGN( 0.5, We(i,j,k-1) )
            zalphaw = 0.5 + z0w
            zw  = z0w - 0.5 * We(i,j,k-1) * rdtstp 
     &           / ( on_r(i,j) * om_r(i,j) * 0.5*(Hz(i,j,k)+Hz(i,j,k-1)) )
            zzwzx = t(i,j,k-1,nadv,itrc) + zw * zslpz_kp1
            zzwzy = t(i,j,k,nadv,itrc) + zw * zslpz(i)
            FC(i,k-1) = We(i,j,k-1) * ( zalphaw * zzwzx + (1.-zalphaw) * zzwzy )
            !    !-- vertical advective trend at jk
            zdzt_kp1(i) = zdzt_kp2
            zslpz   (i) = zslpz_kp1
	    end do
         END DO  ! end of jk loop

	 do i=Istr,Iend
	    zdzt_kp2 = 0.0
            zslpz_kp1 = ( zdzt_kp1(i) + zdzt_kp2 )
     &        * (  0.25 + SIGN( 0.25, zdzt_kp1(i) * zdzt_kp2 )  )
            !  ! slope limitation
            zslpz_kp1 = SIGN( 1.0, zslpz_kp1 )
     &                 * MIN(    ABS( zslpz_kp1 ),
     &                 2.*ABS( zdzt_kp2 ),
     &                 2.*ABS( zdzt_kp1(i) )  )
            !  !-- vertical advective flux at k+1
            z0w = SIGN( 0.5, We(i,j,1) )
            zalphaw = 0.5 + z0w
            zw  = z0w - 0.5 * We(i,j,1) * rdtstp
     &           / ( on_r(i,j) * om_r(i,j) * 0.5*(Hz(i,j,1)+Hz(i,j,2)) )
            zzwzx = t(i,j,1,nadv,itrc) + zw * zslpz_kp1
            zzwzy = t(i,j,2,nadv,itrc) + zw * zslpz(i)
            FC(i,1) = We(i,j,1) * ( zalphaw * zzwzx + (1.-zalphaw) * zzwzy )
         end do

         do i=Istr,Iend
           FC(i,0 )=0.
           FC(i,N )=0.
           CF(i,0 )=dt*pm(i,j)*pn(i,j)
         enddo
