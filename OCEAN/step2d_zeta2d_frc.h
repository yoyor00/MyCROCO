 
# ifdef PRED_COUPLED_MODE
        if (FIRST_TIME_STEP) then
          cff3=0.                        ! This version is designed
          cff2=0.                        ! for coupling during 3D
          cff1=1.                        ! predictor sub-step: here
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
          cff3=0.                        ! computed as instantaneous
          cff2=-0.5                      ! value at 3D time step
          cff1=1.5                       ! "nstp" first, and then
        else                             ! extrapolated half-step
          cff3=0.281105                  ! forward using  AM3-like
          cff2=-0.5-2.*cff3              ! weights optimized for
          cff1=1.5+cff3                  ! maximum stability (with
        endif                            ! special care for startup)
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=rufrc(i,j)-rubar(i,j)
            rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                          + cff3*rufrc_bak(i,j,nstp)
            rufrc_bak(i,j,nstp)=cff
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=rvfrc(i,j)-rvbar(i,j)
            rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                          + cff3*rvfrc_bak(i,j,nstp)
            rvfrc_bak(i,j,nstp)=cff
          enddo
        enddo
# else
        do j=Jstr,Jend                       ! This version is
          do i=Istr,Iend                     ! designed for coupling
            rufrc(i,j)=rufrc(i,j)-rubar(i,j) ! during 3D corrector
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j) ! sub-step: no forward
          enddo                              ! extrapolation needs
        enddo                                ! to be performed.
# endif
# ifndef NBQ_ZETAW
        if (FIRST_TIME_STEP) then
          cff2=0.                        ! 
          cff1=1.                        !
        else                             !
          cff2=-0.5
          cff1=1.5
        endif
# endif
!
!-----------------------------------------------------------------------
! Since coupling requires that pressure gradient term is computed
! using zeta(:,:,kstp) instead of zeta_new(:,:) needed to achieve
! numerical stability, apply compensation to shift pressure gradient
! terms from "kstp" to "knew": in essense, convert the fist 2D step
! from Forward Euler to Forward-Backward].
!-----------------------------------------------------------------------
!
# define zwrk UFx
# define rzeta  UFe
# define rzeta2  VFe
# define rzetaSA VFx

!# ifdef NBQ_MASS
!        do j=JstrV-1,Jend
!          do i=IstrU-1,Iend
!            zwrk(i,j)=zeta_new(i,j) /rhobar_nbq(i,j,knew)
!     &               -zeta(i,j,kstp)/rhobar_nbq(i,j,kstp)
!#  if defined VAR_RHO_2D && defined SOLVE3D
!            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
!            rzeta2(i,j)=rzeta(i,j)*
!     &          (( zeta_new(i,j) /rhobar_nbq(i,j,knew)
!     &            +zeta(i,j,kstp)/rhobar_nbq(i,j,kstp) )
!     &            -2.*h(i,j))
!            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
!#  else
!            rzeta(i,j)=zwrk(i,j)
!            rzeta2(i,j)=zwrk(i,j)*
!     &           ((zeta_new(i,j) /rhobar_nbq(i,j,knew)
!     &            +zeta(i,j,kstp)/rhobar_nbq(i,j,kstp))
!     &            -2.*h(i,j))
!#  endif
!          enddo
!        enddo
!# else
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            zwrk(i,j)=zeta_new(i,j)-zeta(i,j,kstp)
#  if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
#  else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
#  endif
          enddo
        enddo
!# endif /* NBQ_MASS */

        cff=0.5*g
        do j=Jstr,Jend
          do i=Istr,Iend
            rubar(i,j)=rubar(i,j) +cff*on_u(i,j)*( (h(i-1,j)+h(i,j))
     &          *(rzeta(i-1,j)-rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
# endif
     &                                                              )
!
            rvbar(i,j)=rvbar(i,j) +cff*om_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(rzeta(i,j-1)-rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

#  if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(zwrk(i,j-1)-zwrk(i,j)) )
#  endif
     &                                                              )
          enddo
        enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

# undef rzetaSA
# undef rzeta2
# undef rzeta
# undef zwrk
