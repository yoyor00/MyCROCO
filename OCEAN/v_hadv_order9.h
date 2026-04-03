!
!===============================================================
!
! Compute 7th order horizontal advection (v-momentum)
!
!===============================================================
!
#  ifdef NS_PERIODIC
          jmin=0
          jmax=LOCALMM
#  else
#   ifdef MPI
          if (SOUTH_INTER) then
            jmin=0
          else
            jmin=5
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi
          else
            jmax=Mmmpi-4
          endif
#   else
          jmin=5
          jmax=Mm-4
#   endif
#  endif

#  ifdef EW_PERIODIC
          imin=1
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=1
          else
            imin=5
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-3
          endif
#   else
          imin=5
          imax=Lm-3
#   endif
#  endif
!
!----------------------------------------------------------------------
!  VFe interior (j loop)
!----------------------------------------------------------------------
!
!$acc loop independent
          DO j = max(JstrV-1,jmin),min(Jend,jmax)
            DO i = Istr,Iend

              vel = flux10(Hvom(i,j-4,k),Hvom(i,j-3,k),
     &                     Hvom(i,j-2,k),Hvom(i,j-1,k),
     &                     Hvom(i,j  ,k),Hvom(i,j+1,k),
     &                     Hvom(i,j+2,k),Hvom(i,j+3,k),
     &                     Hvom(i,j+4,k),Hvom(i,j+5,k),1.)
              flx9 = vel*FLUX9(
     &             v(i,j-4,k,nrhs), v(i,j-3,k,nrhs),
     &             v(i,j-2,k,nrhs), v(i,j-1,k,nrhs),
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs),
     &             v(i,j+2,k,nrhs), v(i,j+3,k,nrhs),
     &             v(i,j+4,k,nrhs), v(i,j+5,k,nrhs), vel )

#  ifdef MASKING
              vel  = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
              flx7 = vel*FLUX7(
     &             v(i,j-3,k,nrhs), v(i,j-2,k,nrhs),
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs),
     &             v(i,j+3,k,nrhs), v(i,j+4,k,nrhs), vel )
              flx5 = vel*FLUX5(
     &             v(i,j-2,k,nrhs), v(i,j-1,k,nrhs),
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs),
     &             v(i,j+2,k,nrhs), v(i,j+3,k,nrhs), vel )
              flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
              mask2=vmask(i,j-1)*vmask(i,j+2)
              mask3=vmask(i,j-2)*vmask(i,j+3)*mask2
              mask4=vmask(i,j-3)*vmask(i,j+4)*mask3
#   ifdef UP9_MASKING
              IF (vel.gt.0) THEN
                mask5=vmask(i,j-4)*mask4
                mask4=vmask(i,j-3)*mask3
                mask3=vmask(i,j-2)*mask2
                mask2=vmask(i,j-1)
              ELSE
                mask5=vmask(i,j+5)*mask4
                mask4=vmask(i,j+4)*mask3
                mask3=vmask(i,j+3)*mask2
                mask2=vmask(i,j+2)
              ENDIF
#   else
              mask5=vmask(i,j-4)*vmask(i,j+5)*mask4
#   endif
              VFe(i,j)=flx9*mask5
     &                +flx7*(1-mask5)*mask4
     &                +flx5*(1-mask5)*(1-mask4)*mask3
     &                +flx3*(1-mask5)*(1-mask4)*(1-mask3)*mask2
     &                +flx2*(1-mask5)*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              VFe(i,j)=flx9
#  endif
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
! VFe boundary degradation
!----------------------------------------------------------------------
!
!$acc loop independent
          DO j = JstrV-1,Jend
            DO i = Istr,Iend

              vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))

              IF ( j.eq.jmin-4 .or. j.eq.jmax+4 ) THEN
!
! ---- 2nd order ----
!
                VFe(i,j) = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)

              ELSE IF ( j.eq.jmin-3 .or. j.eq.jmax+3 ) THEN
!
! ---- 3rd order with masking ----
!
                flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel)

#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask2=vmask(i,j-1)
                ELSE
                  mask2=vmask(i,j+2)
                ENDIF
#   else
                mask2=vmask(i,j-1)*vmask(i,j+2)
#   endif
                VFe(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                VFe(i,j)=flx3
#  endif

              ELSE IF ( j.eq.jmin-2 .or. j.eq.jmax+2 ) THEN
!
! ---- 5th order with masking ----
!
                flx5 = vel*FLUX5(
     &             v(i,j-2,k,nrhs), v(i,j-1,k,nrhs),
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs),
     &             v(i,j+2,k,nrhs), v(i,j+3,k,nrhs), vel)
#  ifdef MASKING
                flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel)
                flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                mask2=vmask(i,j-1)*vmask(i,j+2)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask3=vmask(i,j-2)*mask2
                  mask2=vmask(i,j-1)
                ELSE
                  mask3=vmask(i,j+3)*mask2
                  mask2=vmask(i,j+2)
                ENDIF
#   else
                mask3=vmask(i,j-2)*vmask(i,j+3)*mask2
#   endif
                VFe(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                              (1-mask3)*(1-mask2)*flx2
#  else
                VFe(i,j)=flx5
#  endif
              ELSE IF ( j.eq.jmin-1 .or. j.eq.jmax+1 ) THEN
!
! ---- 7th order with masking ----
!
                flx7 = vel*FLUX7(
     &             v(i,j-3,k,nrhs), v(i,j-2,k,nrhs),
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs),
     &             v(i,j+3,k,nrhs), v(i,j+4,k,nrhs), vel )
#  ifdef MASKING
                flx5 = vel*FLUX5(
     &             v(i,j-2,k,nrhs), v(i,j-1,k,nrhs),
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs),
     &             v(i,j+2,k,nrhs), v(i,j+3,k,nrhs), vel)
                flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel)
                flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                mask2=vmask(i,j-1)*vmask(i,j+2)
                mask3=vmask(i,j-2)*vmask(i,j+3)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask4=vmask(i,j-3)*mask3
                  mask3=vmask(i,j-2)*mask2
                  mask2=vmask(i,j-1)
                ELSE
                  mask4=vmask(i,j+4)*mask3
                  mask3=vmask(i,j+3)*mask2
                  mask2=vmask(i,j+2)
                ENDIF
#   else
                mask4=vmask(i,j-3)*vmask(i,j+4)*mask3
#   endif
                VFe(i,j)=flx7*mask4
     &                  +flx5*(1-mask4)*mask3
     &                  +flx3*(1-mask4)*(1-mask3)*mask2
     &                  +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
                VFe(i,j)=flx7
#  endif
              ENDIF
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
!  VFx interior (i loop)
!----------------------------------------------------------------------
!
          DO j = JstrV,Jend
!$acc loop independent
            DO i = max(Istr,imin),min(Iend+1,imax)

              if ( j.ge.jmin+1 .and. j.le.jmax ) then
                vel = flux10(Huon(i,j-5,k),Huon(i,j-4,k),
     &                       Huon(i,j-3,k),Huon(i,j-2,k),
     &                       Huon(i,j-1,k),Huon(i,j  ,k),
     &                       Huon(i,j+1,k),Huon(i,j+2,k),
     &                       Huon(i,j+3,k),Huon(i,j+4,k),1.)
              else
                vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
              endif

              flx9 = vel*FLUX9(
     &             v(i-5,j,k,nrhs), v(i-4,j,k,nrhs),
     &             v(i-3,j,k,nrhs), v(i-2,j,k,nrhs),
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs),
     &             v(i+1,j,k,nrhs), v(i+2,j,k,nrhs),
     &             v(i+3,j,k,nrhs), v(i+4,j,k,nrhs), vel )

#  ifdef MASKING
              vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
              flx7 = vel*FLUX7(
     &             v(i-4,j,k,nrhs), v(i-3,j,k,nrhs),
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs),
     &             v(i+2,j,k,nrhs), v(i+3,j,k,nrhs), vel )
              flx5 = vel*FLUX5(
     &             v(i-3,j,k,nrhs), v(i-2,j,k,nrhs),
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs),
     &             v(i+1,j,k,nrhs), v(i+2,j,k,nrhs), vel )
              flx3 = vel*FLUX3(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
              mask2=vmask(i-2,j)*vmask(i+1,j)
              mask3=vmask(i-3,j)*vmask(i+2,j)*mask2
              mask4=vmask(i-4,j)*vmask(i+3,j)*mask3
#   ifdef UP9_MASKING
              IF (vel.gt.0) THEN
                mask5=vmask(i-5,j)*mask4
                mask4=vmask(i-4,j)*mask3
                mask3=vmask(i-3,j)*mask2
                mask2=vmask(i-2,j)
              ELSE
                mask5=vmask(i+4,j)*mask4
                mask4=vmask(i+3,j)*mask3
                mask3=vmask(i+2,j)*mask2
                mask2=vmask(i+1,j)
              ENDIF
#   else
              mask5=vmask(i-5,j)*vmask(i+4,j)*mask4
#   endif
              VFx(i,j)=flx9*mask5
     &                +flx7*(1-mask5)*mask4
     &                +flx5*(1-mask5)*(1-mask4)*mask3
     &                +flx3*(1-mask5)*(1-mask4)*(1-mask3)*mask2
     &                +flx2*(1-mask5)*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              VFx(i,j)=flx9
#  endif
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------
! VFx boundary degradation
!----------------------------------------------------------------------
!
          DO j = JstrV,Jend
!$acc loop independent
            DO i = Istr,Iend+1

              vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))

              IF ( i.eq.imin-4 .or. i.eq.imax+4 ) THEN
!
! ---- 2nd order ----
!
                VFx(i,j) = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)

              ELSE IF ( i.eq.imin-3 .or. i.eq.imax+3 ) THEN
!
! ---- 3rd order with masking ----
!
                flx3 = vel*FLUX3(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel)
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask2=vmask(i-2,j)
                ELSE
                  mask2=vmask(i+1,j)
                ENDIF
#   else
                mask2=vmask(i-2,j)*vmask(i+1,j)
#   endif
                VFx(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                VFx(i,j)=flx3
#  endif

              ELSE IF ( i.eq.imin-2 .or. i.eq.imax+2 ) THEN
!
! ---- 5th order with masking ----
!
                flx5 = vel*FLUX5(
     &             v(i-3,j,k,nrhs), v(i-2,j,k,nrhs),
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs),
     &             v(i+1,j,k,nrhs), v(i+2,j,k,nrhs), vel)
#  ifdef MASKING
                flx3 = vel*FLUX3(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel)
                flx2 = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
                mask2=vmask(i-2,j)*vmask(i+1,j)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask3=vmask(i-3,j)*mask2
                  mask2=vmask(i-2,j)
                ELSE
                  mask3=vmask(i+2,j)*mask2
                  mask2=vmask(i+1,j)
                ENDIF
#   else
                mask3=vmask(i-3,j)*vmask(i+2,j)*mask2
#   endif
                VFx(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                              (1-mask3)*(1-mask2)*flx2
#  else
                VFx(i,j)=flx5
#  endif
              ELSE IF ( i.eq.imin-2 .or. i.eq.imax+2 ) THEN
!
! ---- 7th order with masking ----
!
                flx7 = vel*FLUX7(
     &             v(i-4,j,k,nrhs), v(i-3,j,k,nrhs),
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs),
     &             v(i+2,j,k,nrhs), v(i+3,j,k,nrhs), vel )
#  ifdef MASKING
                flx5 = vel*FLUX5(
     &             v(i-3,j,k,nrhs), v(i-2,j,k,nrhs),
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs),
     &             v(i+1,j,k,nrhs), v(i+2,j,k,nrhs), vel)
                flx3 = vel*FLUX3(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel)
                flx2 = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
                mask2=vmask(i-2,j)*vmask(i+1,j)
                mask3=vmask(i-3,j)*vmask(i+2,j)
#   ifdef UP9_MASKING
                IF (vel.gt.0) THEN
                  mask4=vmask(i-4,j)*mask3
                  mask3=vmask(i-3,j)*mask2
                  mask2=vmask(i-2,j)
                ELSE
                  mask4=vmask(i+3,j)*mask3
                  mask3=vmask(i+2,j)*mask2
                  mask2=vmask(i+1,j)
                ENDIF
#   else
                mask4=vmask(i-4,j)*vmask(i+3,j)*mask3
#   endif
                VFx(i,j)=flx7*mask4
     &                  +flx5*(1-mask4)*mask3
     &                  +flx3*(1-mask4)*(1-mask3)*mask2
     &                  +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
                VFx(i,j)=flx7
#  endif
              ENDIF
            ENDDO
          ENDDO

