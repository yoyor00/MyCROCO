!
!===============================================================
!
! Compute 7th order horizontal advection (u-component)
!
!===============================================================
!
#  ifdef NS_PERIODIC
          jmin=1
          jmax=LOCALMM+1
#  else
#   ifdef MPI
          if (SOUTH_INTER) then
            jmin=1
          else
            jmin=5
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi+1
          else
            jmax=Mmmpi-3
          endif
#   else
          jmin=5
          jmax=Mm-3
#   endif
#  endif

#  ifdef EW_PERIODIC
          imin=0
          imax=LOCALLM
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=0
          else
            imin=5
          endif
          if (EAST_INTER) then
            imax=Lmmpi
          else
            imax=Lmmpi-4
          endif
#   else
          imin=5
          imax=Lm-4
#   endif
#  endif
!
!----------------------------------------------------------------------
!  UFe interior (j loop)
!----------------------------------------------------------------------
!
!$acc loop independent
          DO j = max(Jstr,jmin),min(Jend+1,jmax)
            DO i = IstrU,Iend

              if ( i.ge.imin+1 .and. i.le.imax ) then
                vel = flux10(Hvom(i-5,j,k),Hvom(i-4,j,k),
     &                       Hvom(i-3,j,k),Hvom(i-2,j,k),
     &                       Hvom(i-1,j,k),Hvom(i  ,j,k),
     &                       Hvom(i+1,j,k),Hvom(i+2,j,k),
     &                       Hvom(i+3,j,k),Hvom(i+4,j,k),1.)
              else
                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
              endif
              flx9 = vel*FLUX9(
     &             u(i,j-5,k,nrhs), u(i,j-4,k,nrhs),
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs),
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs),
     &             u(i,j+3,k,nrhs), u(i,j+4,k,nrhs), vel )
#  ifdef MASKING
              mask2=umask(i,j-2)*umask(i,j+1)
              mask3=umask(i,j-3)*umask(i,j+2)*mask2
              mask4=umask(i,j-4)*umask(i,j+3)*mask3
              mask5=umask(i,j-5)*umask(i,j+4)*mask4
#   ifdef UP9_MASKING
              IF (vel.gt.0) THEN
                mask5=umask(i,j-5)*mask4
                mask4=umask(i,j-4)*mask3
                mask3=umask(i,j-3)*mask2
                mask2=umask(i,j-2)
              ELSE
                mask5=umask(i,j+4)*mask4
                mask4=umask(i,j+3)*mask3
                mask3=umask(i,j+2)*mask2
                mask2=umask(i,j+1)
              ENDIF
#   endif
              flx7 = vel*FLUX7(
     &             u(i,j-4,k,nrhs), u(i,j-3,k,nrhs),
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs),
     &             u(i,j+2,k,nrhs), u(i,j+3,k,nrhs), vel )
              flx5 = vel*FLUX5(
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs),
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs), vel )
              flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)

              UFe(i,j)=flx9*mask5
     &               +flx7*(1-mask5)*mask4
     &               +flx5*(1-mask5)*(1-mask4)*mask3
     &               +flx3*(1-mask5)*(1-mask4)*(1-mask3)*mask2
     &               +flx2*(1-mask5)*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              UFe(i,j)=flx9
#  endif
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
!  UFe degradation
!----------------------------------------------------------------------
!
!$acc loop independent
          DO j = Jstr,Jend+1
            DO i = IstrU,Iend

              if ( i.ge.imin+1 .and. i.le.imax ) then
                vel = flux10(Hvom(i-5,j,k),Hvom(i-4,j,k),
     &                       Hvom(i-3,j,k),Hvom(i-2,j,k),
     &                       Hvom(i-1,j,k),Hvom(i  ,j,k),
     &                       Hvom(i+1,j,k),Hvom(i+2,j,k),
     &                       Hvom(i+3,j,k),Hvom(i+4,j,k),1.)
              else
                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
              endif

              IF ( j.eq.jmin-4 .or. j.eq.jmax+4 ) THEN
!
! ---- 2nd order ----
!
                UFe(i,j) = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)

              ELSE IF ( j.eq.jmin-3 .or. j.eq.jmax+3 ) THEN
!
! ---- 3rd order with masking ----
!
                flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), vel)
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
                mask2=umask(i,j-2)*umask(i,j+1)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask2=umask(i,j-2)
                ELSE
                  mask2=umask(i,j+1)
                ENDIF
#   endif
                UFe(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                UFe(i,j)=flx3
#  endif

              ELSE IF ( j.eq.jmin-2 .or. j.eq.jmax+2 ) THEN
!
! ---- 5th order with masking ----
!
                flx5 = vel*FLUX5(
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs),
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs), vel)
#  ifdef MASKING
                flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j,k,nrhs),   u(i,j+1,k,nrhs), vel)
                flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
                mask2=umask(i,j-2)*umask(i,j+1)
                mask3=umask(i,j-3)*umask(i,j+2)*mask2
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask3=umask(i,j-3)*mask2
                  mask2=umask(i,j-2)
                ELSE
                  mask3=umask(i,j+2)*mask2
                  mask2=umask(i,j+1)
                ENDIF
#   endif
                UFe(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                              (1-mask3)*(1-mask2)*flx2
#  else
                UFe(i,j)=flx5
#  endif
              ELSE IF ( j.eq.jmin-1 .or. j.eq.jmax+1 ) THEN
!
! ---- 7th order with masking ----
!
                flx7 = vel*FLUX7(
     &             u(i,j-4,k,nrhs), u(i,j-3,k,nrhs),
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs),
     &             u(i,j+2,k,nrhs), u(i,j+3,k,nrhs), vel )
#  ifdef MASKING
                flx5 = vel*FLUX5(
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs),
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs), vel)
                flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j,k,nrhs),   u(i,j+1,k,nrhs), vel)
                flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
                mask2=umask(i,j-2)*umask(i,j+1)
                mask3=umask(i,j-3)*umask(i,j+2)*mask2
                mask4=umask(i,j-4)*umask(i,j+3)*mask3
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask4=umask(i,j-4)*mask3
                  mask3=umask(i,j-3)*mask2
                  mask2=umask(i,j-2)
                ELSE
                  mask4=umask(i,j+3)*mask3
                  mask3=umask(i,j+2)*mask2
                  mask2=umask(i,j+1)
                ENDIF
#   endif
                UFe(i,j)=flx7*mask4
     &                  +flx5*(1-mask4)*mask3
     &                  +flx3*(1-mask4)*(1-mask3)*mask2
     &                  +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
                UFe(i,j)=flx7
#  endif

              ENDIF
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
!  UFx interior (i loop)
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend
!$acc loop independent
            DO i = max(IstrU-1,imin),min(Iend,imax)

              vel = flux10(Huon(i-4,j,k),Huon(i-3,j,k),
     &                     Huon(i-2,j,k),Huon(i-1,j,k),
     &                     Huon(i  ,j,k),Huon(i+1,j,k),
     &                     Huon(i+2,j,k),Huon(i+3,j,k),
     &                     Huon(i+4,j,k),Huon(i+5,j,k),1.)
              flx9 = vel*FLUX9(
     &             u(i-4,j,k,nrhs), u(i-3,j,k,nrhs),
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs),
     &             u(i+4,j,k,nrhs), u(i+5,j,k,nrhs), vel )

#  ifdef MASKING
              mask2=umask(i-1,j)*umask(i+2,j)
              mask3=umask(i-2,j)*umask(i+3,j)*mask2
              mask4=umask(i-3,j)*umask(i+4,j)*mask3
              mask5=umask(i-4,j)*umask(i+5,j)*mask4
#   ifdef UP9_MASKING
              IF (vel.gt.0) THEN
                mask5=umask(i-4,j)*mask4
                mask4=umask(i-3,j)*mask3
                mask3=umask(i-2,j)*mask2
                mask2=umask(i-1,j)
              ELSE
                mask5=umask(i+5,j)*mask4
                mask4=umask(i+4,j)*mask3
                mask3=umask(i+3,j)*mask2
                mask2=umask(i+2,j)
              ENDIF
#   endif
              IF (mask4.eq.1.) THEN
                vel = flux8(Huon(i-3,j,k),Huon(i-2,j,k),
     &                      Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),
     &                      Huon(i+3,j,k),Huon(i+4,j,k),1.)
                flx7 = vel*FLUX7(
     &             u(i-3,j,k,nrhs), u(i-2,j,k,nrhs),
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs),
     &             u(i+3,j,k,nrhs), u(i+4,j,k,nrhs), vel )
              ELSE IF (mask3.eq.1.) THEN
                vel = flux6(Huon(i-2,j,k),Huon(i-1,j,k),
     &                      Huon(i  ,j,k),Huon(i+1,j,k),
     &                      Huon(i+2,j,k),Huon(i+3,j,k),1.)
                flx5 = vel*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), vel )
              ELSE IF (mask2.eq.1.) THEN
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
              ELSE
                vel  = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
              ENDIF
              UFx(i,j)=flx9*mask5
     &                +flx7*(1-mask5)*mask4
     &                +flx5*(1-mask5)*(1-mask4)*mask3
     &                +flx3*(1-mask5)*(1-mask4)*(1-mask3)*mask2
     &                +flx2*(1-mask5)*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              UFx(i,j)=flx9
#  endif
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
!  UFx degradation
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend
!$acc loop independent
            DO i = IstrU-1,Iend

              vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))

              IF ( i.eq.imin-4 .or. i.eq.imax+4 ) THEN
!
! ---- 2nd order ----
!
                UFx(i,j) = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)

              ELSE IF ( i.eq.imin-3 .or. i.eq.imax+3 ) THEN
!
! ---- 3rd order with masking ----
!
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
#  ifdef MASKING
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
                mask2=umask(i-1,j)*umask(i+2,j)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask2=umask(i-1,j)
                ELSE
                  mask2=umask(i+2,j)
                ENDIF
#   endif
                UFx(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                UFx(i,j)=flx3
#  endif

              ELSE IF ( i.eq.imin-2 .or. i.eq.imax+2 ) THEN
!
! ---- 5th order with masking ----
!
                vel = flux6(Huon(i-2,j,k),Huon(i-1,j,k),
     &                      Huon(i  ,j,k),Huon(i+1,j,k),
     &                      Huon(i+2,j,k),Huon(i+3,j,k),1.)
                flx5 = vel*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), vel )
#  ifdef MASKING
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
                mask2=umask(i-1,j)*umask(i+2,j)
                mask3=umask(i-2,j)*umask(i+3,j)*mask2
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask3=umask(i-2,j)*mask2
                  mask2=umask(i-1,j)
                ELSE
                  mask3=umask(i+3,j)*mask2
                  mask2=umask(i+2,j)
                ENDIF
#   endif
                UFx(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                              (1-mask3)*(1-mask2)*flx2
#  else
                UFx(i,j)=flx5
#  endif
              ELSE IF ( i.eq.imin-1 .or. i.eq.imax+1 ) THEN
!
! ---- 7th order with masking ----
!
                vel = flux8(Huon(i-3,j,k),Huon(i-2,j,k),
     &                      Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),
     &                      Huon(i+3,j,k),Huon(i+4,j,k),1.)
                flx7 = vel*FLUX7(
     &             u(i-3,j,k,nrhs), u(i-2,j,k,nrhs),
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs),
     &             u(i+3,j,k,nrhs), u(i+4,j,k,nrhs), vel )
#  ifdef MASKING
                vel = flux6(Huon(i-2,j,k),Huon(i-1,j,k),
     &                      Huon(i  ,j,k),Huon(i+1,j,k),
     &                      Huon(i+2,j,k),Huon(i+3,j,k),1.)
                flx5 = vel*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), vel )
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
                mask2=umask(i-1,j)*umask(i+2,j)
                mask3=umask(i-2,j)*umask(i+3,j)*mask2
                mask4=umask(i-3,j)*umask(i+4,j)*mask3
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask4=umask(i-3,j)*mask3
                  mask3=umask(i-2,j)*mask2
                  mask2=umask(i-1,j)
                ELSE
                  mask4=umask(i+4,j)*mask3
                  mask3=umask(i+3,j)*mask2
                  mask2=umask(i+2,j)
                ENDIF
#   endif
                UFx(i,j)=flx7*mask4
     &                  +flx5*(1-mask4)*mask3
     &                  +flx3*(1-mask4)*(1-mask3)*mask2
     &                  +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
                UFx(i,j)=flx7
#  endif
              ENDIF
            ENDDO
          ENDDO

