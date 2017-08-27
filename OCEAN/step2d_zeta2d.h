
      if (FIRST_2D_STEP) then
#ifdef SOLVE3D
        cff0=0.                !---> Compute pressure-gradient
        cff1=1.                !  terms using just zeta(:,:,kstp)
#else
        cff0=1.
        cff1=0.
#endif
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        cff0= 1.0833333333333    ! Logically AB2-AM3 forward-backward
        cff1=-0.1666666666666    ! scheme with coefficients chosen for
        cff2= 0.0833333333333    ! maximum stability, while maintaining
        cff3= 0.                 ! third-accuracy; alpha_max=1.73
      else
#ifdef M2FILTER_NONE
        cff0=0.5+2.*myepsilon+mygamma+2.*myalpha
        cff1=1.-cff0-mygamma-myepsilon
        cff2=mygamma
        cff3=myepsilon
#else
        cff0=0.614
        cff1=0.285
        cff2=0.088
        cff3=0.013
#endif
      endif
 
#define zwrk UFx
#define rzeta  UFe
#define rzeta2  VFe
#define rzetaSA VFx
!
!-----------------------------------------------------------------------
! Computes zeta(n+1)
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zeta(n+1) from depth-integral 
! mass (or volume) conservation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1))
# ifdef ZONAL_NUDGING
          zeta(i,j,knew)=zeta_new(i,j)
# endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
! Add nudging terms
!-----------------------------------------------------------------------
!
#ifdef ZNUDGING
# ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew),zetazon)
        endif
      endif
      if (iic.eq.ntstart) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ssh(START_2D_ARRAY),sshzon)
        endif
      endif
# endif  /* ZONAL_NUDGING */
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j) + dtfast*Znudgcof(i,j)
# ifdef ZONAL_NUDGING
     &                                 *(sshzon(j)-zetazon(j))
# else
     &                                 *(ssh(i,j)-zeta_new(i,j))
# endif /* ZONAL_NUDGING */
        enddo
      enddo
#endif /* ZNUDGING */
 
