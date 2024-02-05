!
!----------------------------------------------------------
! Compute adiabatic buoyancy gradients
! used to diagnose diapycnal fluxes and diapycnal velocity 
!----------------------------------------------------------
!

!----------------------------
! Compute vertical component
!----------------------------
!
      cff=g/rho0
      do j=JstrV-1,Jend
        do k=1,N-1
          do i=IstrU-1,Iend
#  ifdef SPLIT_EOS
            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            
            cff2=rho1(i,j,k+1)-rho1(i,j,k)               ! Elementary
     &              +(qp1(i,j,k+1)-qp1(i,j,k))           ! adiabatic
     &                     *dpth*(1.-qp2*dpth)           ! difference

#  else
            cff2=rho(i,j,k+1)-rho(i,j,k)
#  endif
            dbdz(i,j,k)=-cff*cff2 / (z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo
        do i=IstrU-1,Iend
          dbdz(i,j,N)=dbdz(i,j,N-1)
          dbdz(i,j,0)=dbdz(i,j,1)
        enddo
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


! Compute XI-component
!-------- ------------ 
!
      do k=N,1,-1    !<-- k loop

        do j=Jstr,Jend
          do i=imin,imax
#  ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))

            cff2=( rho1(i,j,k)-rho1(i-1,j,k)          ! Elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )           ! difference
#   else
            cff2=(rho(i,j,k)-rho(i-1,j,k))
#   endif
            dbdx(i,j,k)=-cff*cff2
     &              * 0.5 * (pm(i,j)+pm(i-1,j))

#   ifdef MASKING
     &                              *umask(i,j)
#   endif



          enddo
        enddo

#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then         ! Extrapolate elementary
          do j=Jstr,Jend               ! differences near physical
            dbdx(imin-1,j,k)=dbdx(imin,j,k)    ! for reduced loop ranges.
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=Jstr,Jend
            dbdx(imax+1,j,k)=dbdx(imax,j,k)
          enddo
        endif
#   endif


!
! Compute ETA-component
!-------- ------------ 
!


        do j=jmin,jmax
          do i=Istr,Iend
#  ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_r(i,j,k)-z_r(i,j-1,k))

            cff2=( rho1(i,j,k)-rho1(i,j-1,k)             ! Elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
#   else
            cff2=(rho(i,j,k)-rho(i,j-1,k))
#   endif
            dbdy(i,j,k)=-cff*cff2
     &                * 0.5 * (pn(i,j)+pn(i,j-1))
#   ifdef MASKING
     &                              *vmask(i,j)
#   endif
          enddo
        enddo

        if (SOUTHERN_EDGE) then
          do i=Istr,Iend
            dbdy(i,jmin-1,k)=dbdy(i,jmin,k)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=Istr,Iend
            dbdy(i,jmax+1,k)=dbdy(i,jmax,k)
          enddo
        endif

      enddo   ! k loop







