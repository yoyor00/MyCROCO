!----------------------------------------------------------
!
!----------------------------------------------------------
! Second order advection scheme
!----------------------------------------------------------
# if defined PREDICTOR
          rdtstp = (1.0-gamma)*dt
# else
	  rdtstp = dt
# endif

          do j=Jstr-2,Jend+1
            do i=Istr-2,Iend+1
              zdxt(i,j)=(t(i+1,j,k,nadv,itrc)-t(i,j,k,nadv,itrc))
#  ifdef MASKING
     &                                             *umask(i+1,j)
#  endif
              zdyt(i,j)=(t(i,j+1,k,nadv,itrc)-t(i,j,k,nadv,itrc))
#  ifdef MASKING
     &                                             *vmask(i,j+1)
#  endif
            enddo
          enddo

          do j=Jstr-1,Jend+1
            do i=Istr-1,Iend+1
              zzslpx = ( zdxt(i,j) + zdxt(i-1,j) )
     &          * ( 0.25 + SIGN( 0.25, zdxt(i,j) * zdxt (i-1,j) ) )
              zzslpy = ( zdyt(i,j) + zdyt(i,j-1) )
     &          * ( 0.25 + SIGN( 0.25, zdyt(i,j) * zdyt (i,j-1) ) )
              zslpx(i,j) = SIGN( 1.0, zzslpx) * MIN ( ABS( zzslpx),
     &                     2.0 * ABS(zdxt(i-1,j)), 2*ABS(zdxt(i,j)) )
              zslpy(i,j) = SIGN( 1.0, zzslpy) * MIN ( ABS( zzslpy),
     &                     2.0 * ABS(zdyt(i,j-1)), 2*ABS(zdyt(i,j)) )
             ENDDO
          ENDDO

          DO j=Jstr,Jend
             DO i=Istr-1,Iend
               z0u = SIGN( 0.5,Huon(i+1,j,k) )
               zalpha = 0.5 - z0u
               zu  = z0u-0.5*Huon(i+1,j,k)*rdtstp*pn_u(i+1,j)
     &               *pm_u(i+1,j)/(0.5*(Hz(i,j,k)+Hz(i+1,j,k)))
               zzwx = t(i+1,j,k,nadv,itrc) + zu * zslpx(i+1,j)
               zzwy = t(i,j,k,nadv,itrc) + zu * zslpx(i,j)
               FX(i+1,j) = Huon(i+1,j,k)*(zalpha*zzwx+(1.-zalpha)*zzwy)
             END DO
	  END DO

	  DO j=Jstr-1,Jend
             DO i=Istr,Iend
               !
               !
               z0v = SIGN( 0.5, Hvom(i,j+1,k) )
               zalpha = 0.5 - z0v
               zv  = z0v-0.5*Hvom(i,j+1,k)*rdtstp*pn_v(i,j+1)
     &               *pm_v(i,j+1)/(0.5*(Hz(i,j,k)+Hz(i,j+1,k)))
               zzwx = t(i,j+1,k,nadv,itrc) + zv * zslpy(i+1,j)
               zzwy = t(i,j,k,nadv,itrc) + zv * zslpy(i,j)
               FE(i,j+1) = Hvom(i,j+1,k)*(zalpha*zzwx+(1.-zalpha)* zzwy)
            enddo
        enddo
