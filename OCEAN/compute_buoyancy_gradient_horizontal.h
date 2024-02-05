!
!----------------------------------------------------------
! Compute horizontal adiabatic buoyancy gradients
!
! used to diagnose diapycnal fluxes and diapycnal velocity 
!----------------------------------------------------------
!


!----------------------------
! Compute vertical component
!----------------------------
!
      cff=g/rho0
      do j=JstrV-1,Jend
        do k=N,1,-1
          do i=IstrU-1,Iend
            dbdz(i,j,k)=0.5*(dbdz(i,j,k)+dbdz(i,j,k-1))
          enddo ! i
        enddo ! k
       enddo ! j


!-----------------------------------------
! Horizontal buoyancy gradient (adiabatic)
!-----------------------------------------

#   ifndef EW_PERIODIC
      if (WESTERN_EDGE) then     ! Restrict extended ranges one
        imin=IstrU               ! point inward near the physical
      else                       ! boundary. Note that this version
        imin=IstrU-1             ! of code works in MPI configuration
      endif                      ! too, while a more straightforward
      if (EASTERN_EDGE) then     ! loop range setting
        imax=Iend                !
      else                       !   i=max(2,IstrU-1),min(Iend+1,Lm)
        imax=Iend+1              !
      endif                      ! does not.
#   else
      imin=Istr-1
      imax=Iend+1
#   endif

#   ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=JstrV
        else
          jmin=JstrV-1
        endif
        if (NORTHERN_EDGE) then
          jmax=Jend
        else
          jmax=Jend+1
        endif
#   else
        jmin=Jstr-1
        jmax=Jend+1
#   endif

!
! Compute XI-component
!-------- ------------ 
!
      do k=N,1,-1    !<-- k loop

        do j=Jstr,Jend
          do i=imin,imax
            dZx(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
     &                * 0.5 * (pm(i,j)+pm(i-1,j))
#   ifdef MASKING
     &                              *umask(i,j)
#   endif
          enddo
        enddo

#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then         ! Extrapolate elementary
          do j=Jstr,Jend               ! differences near physical
            dZx(imin-1,j)=dZx(imin,j)  ! boundaries to compencate.
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=Jstr,Jend
            dZx(imax+1,j)=dZx(imax,j)
          enddo
        endif
#   endif

        do j=Jstr,Jend
          do i=IstrU-1,Iend

            dZx(i,j)=0.5*(dZx(i,j)+dZx(i+1,j))
            rx(i,j)=0.5*(dbdx(i,j,k)+dbdx(i+1,j,k))

          enddo ! i
        enddo ! j

        do j=Jstr,Jend
          do i=IstrU-1,Iend

            dbdx(i,j,k) = ( rx(i,j)
     &                    - dZx(i,j) * dbdz(i,j,k) )

          enddo ! i
        enddo ! j

!
! Compute ETA-component
!-------- ------------ 
!


        do j=jmin,jmax
          do i=Istr,Iend
            dZx(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
     &                * 0.5 * (pn(i,j)+pn(i,j-1))
#   ifdef MASKING
     &                              *vmask(i,j)
#   endif
          enddo
        enddo

        if (SOUTHERN_EDGE) then
          do i=Istr,Iend
            dZx(i,jmin-1)=dZx(i,jmin)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=Istr,Iend
            dZx(i,jmax+1)=dZx(i,jmax)
          enddo
        endif

        do j=JstrV-1,Jend
          do i=Istr,Iend

            dZx(i,j)=0.5*(dZx(i,j)+dZx(i,j+1))
            rx(i,j)=0.5*(dbdy(i,j,k)+dbdy(i,j+1,k))

          enddo ! i
        enddo ! j

        do j=JstrV-1,Jend
          do i=Istr,Iend

            dbdy(i,j,k) = rx(i,j)
     &                  - dZx(i,j) * dbdz(i,j,k)

          enddo ! i
        enddo ! j

      enddo   ! k loop







