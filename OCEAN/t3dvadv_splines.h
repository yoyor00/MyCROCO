!
!----------------------------------------------------------
! Compute vertical advective fluxes using parabolic splines:
! FC=W*[spline-interpolated tracer]
!----------------------------------------------------------
!
                              ! Construct parabolic splines: here
          do i=Istr,Iend      ! CF is the set of vertical derivatives
            FC(i,0)=0.        ! of the tracer field t(:,:,:,nadv,:),
            CF(i,0)=0.        ! FC is an auxiliary scratch variable.
          enddo
          do k=1,N-1,+1
            do i=Istr,Iend
              cff    = 1./(2.*HZR(i,j,k+1)+HZR(i,j,k)*(2.-FC(i,k-1)))
              FC(i,k)= cff*HZR(i,j,k+1)
              CF(i,k)= cff*( 6.*( t(i,j,k+1,nadv,itrc)
     &                           -t(i,j,k  ,nadv,itrc) )-HZR(i,j,k)*CF(i,k-1)
     &                                                                      )
            enddo
          enddo
          do i=Istr,Iend
            CF(i,N)=0.
          enddo
          do k=N-1,1,-1       !<-- irreversible
            do i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            enddo
          enddo               !--> discard FC, keep CF

          cff=1./3.           ! Compute vertical advective fluxes
          do k=1,N-1          ! FC=W*[spline-interpolated tracer]
            do i=Istr,Iend
              FC(i,k)=We(i,j,k)*( t(i,j,k,nadv,itrc)+cff*HZR(i,j,k)
     &                                  *(CF(i,k)+0.5*CF(i,k-1))
     &                                                            )
            enddo
          enddo               !--> discard CF
          do i=Istr,Iend
            FC(i,N)=0.
            FC(i,0)=0.
            CF(i,0)=dt*pm(i,j)*pn(i,j)
          enddo