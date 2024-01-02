!
!----------------------------------------------------------
! Compute adiabatic buoyancy gradients
! used to diagnose diapycnal fluxes and diapycnal velocity 
!----------------------------------------------------------
!
      do itrc=1,NTA
!----------------------------
! Compute vertical component
!----------------------------
!
      do j=JstrV-1,Jend
        do k=N,1,-1
          do i=IstrU-1,Iend
            cff=0.5*pm(i,j)*pn(i,j)
            TF_zHmix_rot(i,j,k,itrc)= cff *
     &                     ( TF_zHmix(i,j,k  ,itrc) * dbdz(i,j,k)
     &                     + TF_zHmix(i,j,k-1,itrc) * dbdz(i,j,k-1) )
            TF_zVmix_rot(i,j,k,itrc)=cff *
     &                     ( TF_zVmix(i,j,k  ,itrc) * dbdz(i,j,k)
     &                     + TF_zVmix(i,j,k-1,itrc) * dbdz(i,j,k-1) )
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
            cff = dbdx(i,j,k) * 2./(Hz(i,j,k)+Hz(i-1,j,k))
            rx(i,j) = TF_xHmix(i,j,k,itrc) * cff
#   ifdef MASKING
     &                              *umask(i,j)
#   endif
          enddo
        enddo

#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then         ! Extrapolate elementary
          do j=Jstr,Jend               ! differences near physical
            rx(imin-1,j)=rx(imin,j)    ! for reduced loop ranges.
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=Jstr,Jend
            rx(imax+1,j)=rx(imax,j)
          enddo
        endif
#   endif

        do j=Jstr,Jend
          do i=IstrU-1,Iend
            TF_xHmix_rot(i,j,k,itrc)=0.5*(rx(i,j)+rx(i+1,j)) * pn(i,j)
          enddo ! i
        enddo ! j

!
! Compute ETA-component
!-------- ------------ 
!


        do j=jmin,jmax
          do i=Istr,Iend
            cff = dbdy(i,j,k) * 2./(Hz(i,j,k)+Hz(i,j-1,k))
            rx(i,j)=TF_yHmix(i,j,k,itrc) * cff
#   ifdef MASKING
     &                              *vmask(i,j)
#   endif
          enddo
        enddo

        if (SOUTHERN_EDGE) then
          do i=Istr,Iend
            rx(i,jmin-1)=rx(i,jmin)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=Istr,Iend
            rx(i,jmax+1)=rx(i,jmax)
          enddo
        endif

        do j=JstrV-1,Jend
          do i=Istr,Iend
            TF_yHmix_rot(i,j,k,itrc)= 0.5*(rx(i,j)+rx(i,j+1)) * pm(i,j)
          enddo ! i
        enddo ! j

      enddo   ! k loop

      enddo   ! itrc (for T and S)





