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
            jmin=4
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi+1
          else
            jmax=Mmmpi-2
          endif
#   else
          jmin=4
          jmax=Mm-2
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
            imin=4
          endif
          if (EAST_INTER) then
            imax=Lmmpi
          else
            imax=Lmmpi-3
          endif
#   else
          imin=4
          imax=Lm-3
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

              if ( i.ge.imin+2 .and. i.le.imax-1 ) then
                vel = flux8(Hvom(i-4,j,k),Hvom(i-3,j,k),
     &                      Hvom(i-2,j,k),Hvom(i-1,j,k),
     &                      Hvom(i  ,j,k),Hvom(i+1,j,k),
     &                      Hvom(i+2,j,k),Hvom(i+3,j,k),1.)
              else
                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
              endif

              flx7 = vel*FLUX7(
     &             u(i,j-4,k,nrhs), u(i,j-3,k,nrhs),
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs),
     &             u(i,j+2,k,nrhs), u(i,j+3,k,nrhs), vel )

#  ifdef MASKING
              vel  = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
              flx5 = vel*FLUX5(
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs),
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs), vel )
              flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
              mask2=umask(i,j-2)*umask(i,j+1)
              mask3=umask(i,j-3)*umask(i,j+2)*mask2
#   ifdef UP7_MASKING
              IF (vel.gt.0) THEN
                mask4=umask(i,j-4)*mask3
                mask3=umask(i,j-3)*mask2
                mask2=umask(i,j-2)
              ELSE
                mask4=umask(i,j+3)*mask3
                mask3=umask(i,j+2)*mask2
                mask2=umask(i,j+1)
              ENDIF
#   else
              mask4=umask(i,j-4)*umask(i,j+3)*mask3
#   endif
              UFe(i,j)=flx7*mask4
     &                +flx5*(1-mask4)*mask3
     &                +flx3*(1-mask4)*(1-mask3)*mask2
     &                +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              UFe(i,j)=flx7
#  endif
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
!  UFe degradation
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend+1
!$acc loop independent
            DO i = IstrU,Iend

              vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))

              IF ( j.eq.jmin-3 .or. j.eq.jmax+3 ) THEN
!
! ---- 2nd order ----
!
                UFe(i,j) = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)

              ELSE IF ( j.eq.jmin-2 .or. j.eq.jmax+2 ) THEN
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
                UFe(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                UFe(i,j)=flx3
#  endif

              ELSE IF ( j.eq.jmin-1 .or. j.eq.jmax+1 ) THEN
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
#   ifdef UP7_MASKING
                mask2=umask(i,j-2)*umask(i,j+1)
                IF (vel.gt.0) THEN
                  mask3=umask(i,j-3)*mask2
                  mask2=umask(i,j-2)
                ELSE
                  mask3=umask(i,j+2)*mask2
                  mask2=umask(i,j+1)
                ENDIF
#   else
                mask2=umask(i,j-2)*umask(i,j+1)
                mask3=umask(i,j-3)*umask(i,j+2)*mask2
#   endif
                UFe(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                              (1-mask3)*(1-mask2)*flx2
#  else
                UFe(i,j)=flx5
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

              vel = flux8(Huon(i-3,j,k),Huon(i-2,j,k),
     &                    Huon(i-1,j,k),Huon(i  ,j,k),
     &                    Huon(i+1,j,k),Huon(i+2,j,k),
     &                    Huon(i+3,j,k),Huon(i+4,j,k),1.)

              flx7 = vel*FLUX7(
     &             u(i-3,j,k,nrhs), u(i-2,j,k,nrhs),
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs),
     &             u(i+3,j,k,nrhs), u(i+4,j,k,nrhs), vel )

#  ifdef MASKING
              vel  = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
              flx5 = vel*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), vel )
              flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
              mask2=umask(i-1,j)*umask(i+2,j)
              mask3=umask(i-2,j)*umask(i+3,j)*mask2
#   ifdef UP7_MASKING
              IF (vel.gt.0) THEN
                mask4=umask(i-3,j)*mask3
                mask3=umask(i-2,j)*mask2
                mask2=umask(i-1,j)
              ELSE
                mask4=umask(i+4,j)*mask3
                mask3=umask(i+3,j)*mask2
                mask2=umask(i+2,j)
              ENDIF
#   else
              mask4=umask(i-3,j)*umask(i+4,j)*mask3
#   endif
              UFx(i,j)=flx7*mask4
     &                +flx5*(1-mask4)*mask3
     &                +flx3*(1-mask4)*(1-mask3)*mask2
     &                +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              UFx(i,j)=flx7
#  endif
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
!  UFx egradation
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend
!$acc loop independent
            DO i = IstrU-1,Iend

              vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))

              IF ( i.eq.imin-3 .or. i.eq.imax+3 ) THEN
!
! ---- 2nd order ----
!
                UFx(i,j) = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)

              ELSE IF ( i.eq.imin-2 .or. i.eq.imax+2 ) THEN
!
! ---- 3rd order with masking ----
!
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
                mask2 = umask(i-1,j)*umask(i+2,j)
                UFx(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                UFx(i,j)=flx3
#  endif

              ELSE IF ( i.eq.imin-1 .or. i.eq.imax+1 ) THEN
!
! ---- 5th order with masking ----
!
              flx5 = vel*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), vel )
#  ifdef MASKING
              flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
#   ifdef UP7_MASKING
                mask2=umask(i-1,j)*umask(i+2,j)
                IF (vel.gt.0) THEN
                  mask3=umask(i-2,j)*mask2
                  mask2=umask(i-1,j)
                ELSE
                  mask3=umask(i+3,j)*mask2
                  mask2=umask(i+2,j)
                ENDIF
#   else
                mask2=umask(i-1,j)*umask(i+2,j)
                mask3=umask(i-2,j)*umask(i+3,j)*mask2
#   endif
                UFx(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                              (1-mask3)*(1-mask2)*flx2
#  else
                UFx(i,j)=flx5
#  endif
              ENDIF
            ENDDO
          ENDDO


