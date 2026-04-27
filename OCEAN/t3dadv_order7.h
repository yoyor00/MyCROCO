!
!===============================================================
!
! Compute 7th order horizontal advection for tracers
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
          imin=1
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=1
          else
            imin=4
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-2
          endif
#   else
          imin=4
          imax=Lm-2
#   endif
#  endif
!
!----------------------------------------------------------------------
!  FE interior (j loop)
!----------------------------------------------------------------------
!
#ifdef OPENACC
!$acc parallel loop if(compute_on_device) default(present) private(vel)
!$acc&	independent	
          DOEXTEND(k,1,N,FX,FE,WORK)
#endif
!$acc loop independent

          DO j = max(Jstr,jmin),min(Jend+1,jmax)

!$acc loop independent private(vel)
            DO i = Istr,Iend
              vel = Hvom(i,j,k)
              flx7 = vel*FLUX7(
     &             t(i,j-4,k,nrhs,itrc), t(i,j-3,k,nrhs,itrc),
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc),
     &             t(i,j+2,k,nrhs,itrc), t(i,j+3,k,nrhs,itrc), vel )
#  ifdef MASKING
              flx5 = vel*FLUX5(
     &             t(i,j-3,k,nrhs,itrc), t(i,j-2,k,nrhs,itrc),
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc),
     &             t(i,j+1,k,nrhs,itrc), t(i,j+2,k,nrhs,itrc), vel )
              flx3 = vel*FLUX3(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc), vel )
              flx2 = vel*FLUX2(
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc), vel, cdif)
              mask2=rmask(i,j-2)*rmask(i,j+1)
              mask3=rmask(i,j-3)*rmask(i,j+2)*mask2
#   ifdef UP7_MASKING
              IF (vel.gt.0) THEN
                mask4=rmask(i,j-4)*mask3
                mask3=rmask(i,j-3)*mask2
                mask2=rmask(i,j-2)
              ELSE
                mask4=rmask(i,j+3)*mask3
                mask3=rmask(i,j+2)*mask2
                mask2=rmask(i,j+1)
              ENDIF
#   else
              mask4=rmask(i,j-4)*rmask(i,j+3)*mask3
#   endif
              FE(i,j)=flx7*mask4
     &               +flx5*(1-mask4)*mask3
     &               +flx3*(1-mask4)*(1-mask3)*mask2
     &               +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              FE(i,j)=flx7
#  endif
            ENDDO
          ENDDO
#ifdef OPENACC
          ENDDOEXTEND
#endif
!
!----------------------------------------------------------------------
!  FE degradation
!----------------------------------------------------------------------
!
#ifdef OPENACC
!$acc parallel loop if(compute_on_device) default(present) private(vel)
!$acc& independent	
          DOEXTEND(k,1,N,FX,FE,WORK)
#endif          
!$acc loop independent

          DO j = Jstr,Jend+1

!$acc loop independent private(vel)
            DO i = Istr,Iend

              vel = Hvom(i,j,k)

              IF ( j.eq.jmin-3 .or. j.eq.jmax+3 ) THEN
!
! ---- 2nd order ----
!
                FE(i,j)=vel*FLUX2(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)

              ELSE IF ( j.eq.jmin-2 .or. j.eq.jmax+2 ) THEN
!
! ---- 3rd order with masking ----
!
                flx3 = vel*FLUX3(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc), vel)
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc), vel, cdif)
#   ifdef UP7_MASKING
                IF (vel.gt.0) THEN
                  mask2=rmask(i,j-2)
                ELSE
                  mask2=rmask(i,j+1)
                ENDIF
#   else
                mask2 = rmask(i,j-2)*rmask(i,j+1)
#   endif
                FE(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                FE(i,j)=flx3
#  endif

              ELSE IF ( j.eq.jmin-1 .or. j.eq.jmax+1 ) THEN
!
! ---- 5th order with masking ----
!
                flx5 = vel*FLUX5(
     &             t(i,j-3,k,nrhs,itrc), t(i,j-2,k,nrhs,itrc),
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc),
     &             t(i,j+1,k,nrhs,itrc), t(i,j+2,k,nrhs,itrc), vel)
#  ifdef MASKING
                flx3 = vel*FLUX3(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc), vel)
                flx2 = vel*FLUX2(
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc), vel, cdif)
                mask2=rmask(i,j-2)*rmask(i,j+1)
#   ifdef UP7_MASKING
                IF (vel.gt.0) THEN
                  mask3=rmask(i,j-3)*mask2
                  mask2=rmask(i,j-2)
                ELSE
                  mask3=rmask(i,j+2)*mask2
                  mask2=rmask(i,j+1)
                ENDIF
#   else
                mask3=rmask(i,j-3)*rmask(i,j+2)*mask2
#   endif
                FE(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                             (1-mask3)*(1-mask2)*flx2
#  else
                FE(i,j)=flx5
#  endif
              ENDIF
            ENDDO
          ENDDO
#ifdef OPENACC
          ENDDOEXTEND
#endif
!
!----------------------------------------------------------------------
!  FX interior (i loop)
!----------------------------------------------------------------------
!
#ifdef OPENACC
!$acc parallel loop if(compute_on_device) default(present) private(vel)
!$acc& independent	
          DOEXTEND(k,1,N,FX,FE,WORK)
#endif          
!$acc loop independent

          DO i = max(Istr,imin),min(Iend+1,imax)

!$acc loop independent  private(vel)
            DO j = Jstr,Jend
              vel = Huon(i,j,k)
              flx7 = vel*FLUX7(
     &             t(i-4,j,k,nrhs,itrc), t(i-3,j,k,nrhs,itrc),
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc),
     &             t(i+2,j,k,nrhs,itrc), t(i+3,j,k,nrhs,itrc), vel )
#  ifdef MASKING
              flx5 = vel*FLUX5(
     &             t(i-3,j,k,nrhs,itrc), t(i-2,j,k,nrhs,itrc),
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc),
     &             t(i+1,j,k,nrhs,itrc), t(i+2,j,k,nrhs,itrc), vel )
              flx3 = vel*FLUX3(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc), vel )
              flx2 = vel*FLUX2(
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc), vel, cdif)
              mask2=rmask(i-2,j)*rmask(i+1,j)
              mask3=rmask(i-3,j)*rmask(i+2,j)*mask2
#   ifdef UP7_MASKING
              IF (vel.gt.0) THEN
                mask4=rmask(i-4,j)*mask3
                mask3=rmask(i-3,j)*mask2
                mask2=rmask(i-2,j)
              ELSE
                mask4=rmask(i+3,j)*mask3
                mask3=rmask(i+2,j)*mask2
                mask2=rmask(i+1,j)
              ENDIF
#   else
              mask4=rmask(i-4,j)*rmask(i+3,j)*mask3
#   endif
              FX(i,j)=flx7*mask4
     &               +flx5*(1-mask4)*mask3
     &               +flx3*(1-mask4)*(1-mask3)*mask2
     &               +flx2*(1-mask4)*(1-mask3)*(1-mask2)
#  else
              FX(i,j)=flx7
#  endif
            ENDDO
          ENDDO
#ifdef OPENACC
          ENDDOEXTEND
#endif
!
!----------------------------------------------------------------------
!  FX degradation
!----------------------------------------------------------------------
!
#ifdef OPENACC
!$acc parallel loop if(compute_on_device) default(present) private(vel)
!$acc& independent	
          DOEXTEND(k,1,N,FX,FE,WORK)
#endif          
!$acc loop independent

          DO i = Istr,Iend+1
            DO j = Jstr,Jend

!$acc loop independent  private(vel)
              vel = Huon(i,j,k)

              IF ( i.eq.imin-3 .or. i.eq.imax+3 ) THEN
!
! ---- 2nd order ----
!
                FX(i,j) = vel*FLUX2(
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc), vel, cdif)

              ELSE IF ( i.eq.imin-2 .or. i.eq.imax+2 ) THEN
!
! ---- 3rd order with masking ----
!
                flx3 = vel*FLUX3(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc), vel)
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc), vel, cdif)
#   ifdef UP7_MASKING
                IF (vel.gt.0) THEN
                  mask2=rmask(i-2,j)
                ELSE
                  mask2=rmask(i+1,j)
                ENDIF
#   else
                mask2=rmask(i-2,j)*rmask(i+1,j)
#   endif
                FX(i,j)=mask2*flx3+(1-mask2)*flx2
#  else
                FX(i,j)=flx3
#  endif

              ELSE IF ( i.eq.imin-1 .or. i.eq.imax+1 ) THEN
!
! ---- 5th order with masking ----
!
                flx5 = vel*FLUX5(
     &             t(i-3,j,k,nrhs,itrc), t(i-2,j,k,nrhs,itrc),
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc),
     &             t(i+1,j,k,nrhs,itrc), t(i+2,j,k,nrhs,itrc), vel)
#  ifdef MASKING
                flx3 = vel*FLUX3(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc), vel)
                flx2 = vel*FLUX2(
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc), vel, cdif)
                mask2=rmask(i-2,j)*rmask(i+1,j)
#   ifdef UP7_MASKING
                IF (vel.gt.0) THEN
                  mask3=rmask(i-3,j)*mask2
                  mask2=rmask(i-2,j)
                ELSE
                  mask3=rmask(i+2,j)*mask2
                  mask2=rmask(i+1,j)
                ENDIF
#   else
                mask3=rmask(i-3,j)*rmask(i+2,j)*mask2
#   endif
                FX(i,j)=mask3*flx5+(1-mask3)*mask2*flx3+
     &                             (1-mask3)*(1-mask2)*flx2
#  else
                FX(i,j)=flx5
#  endif
              ENDIF
            ENDDO
          ENDDO
#ifdef OPENACC
          ENDDOEXTEND
#endif



