
# ifdef M3FAST_AM4
! !
! !--------------------------------------------------------------------
! ! Surface pressure gradient
! !--------------------------------------------------------------------
! !
! ! Interpolate zeta fields half-step backward (AM4) for the subsequent 
! ! computation of barotropic pressure-gradient
! !
# ifdef TANK
         myalpha   = 0.01
# else
         myalpha   = 0.25
# endif
      myepsilon = 0.00976186 - 0.13451357*myalpha
      mygamma   = 0.08344500 - 0.51358400*myalpha
!
      if (FIRST_FAST_STEP) then
        cff0=0.                  !---> Compute pressure-gradient
        cff1=1.                  !     terms using just zeta(:,:,kstp)
        cff2=0.
        cff3=0.
      elseif (FIRST_FAST_STEP+1) then
        cff0= 1.0833333333333    ! AM3 backward scheme
        cff1=-0.1666666666666    ! with coefficients chosen for
        cff2= 0.0833333333333    ! maximum stability, while maintaining
        cff3= 0.                 ! third-accuracy; alpha_max=1.73
      else
        cff0=0.5+2.*myepsilon+mygamma+2.*myalpha  ! AM4 backward scheme
        cff1=1.-cff0-mygamma-myepsilon            ! with implicit diffusion
        cff2=mygamma                              ! given by myalpha
        cff3=myepsilon
      endif
# endif  /* M3FAST_AM4 */
!
# define zwrk    UFx
# define rzeta   UFe
# define rzeta2  VFe
# define rzetaSA VFx

! ! KERNEL_5  zwrk <= ( zeta )
! ! KERNEL_5  rzeta <= ( rhoS, zwrk )
! ! KERNEL_5  rzeta2 <= ( rzeta, zwrk )
! ! KERNEL_5  rzetaSA <= ( zwrk, rhoS, rhoA )

!$acc kernels default(present)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
# ifndef M3FAST_AM4
          zwrk(i,j)=zeta(i,j,kstp)
# else
          zwrk(i,j)=cff0*zeta(i,j,knew)+cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold)
# endif 
# ifdef NBQ_MASS
     &                                *rhobar_nbq(i,j,kstp) 
# endif
# ifdef VAR_RHO_2D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
        enddo
      enddo
!$acc end kernels

! #ifdef RVTK_DEBUG_ADVANCED
! C$OMP BARRIER
! C$OMP MASTER
!        call check_tab2d(rubar,'rubar st_fast_h','uint')
! #endif

!$acc kernels default(present)
!
! Compute surface pressure gradient
!
      cff=0.5*g
      do j=Jstr,Jend
        do i=Istr,Iend
# ifdef M3FAST_BOTH
          rubar(i,j)=0.
          rubarh(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)    ! ATTENTION: h here
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
#  else
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)    ! ATTENTION: h here
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)   
#  endif     
# ifdef VAR_RHO_2D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(zwrk(i-1,j)-zwrk(i,j)))
# endif
# ifdef MRL_WCI
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i-1,j) )
# endif
     &                                                              )
 
# ifdef M3FAST_BOTH
          rvbar(i,j)=0.
          rvbarh(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)
# else
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)
# endif     
# ifdef VAR_RHO_2D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(zwrk(i,j-1)-zwrk(i,j)))
# endif
# ifdef MRL_WCI
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i,j-1) )
# endif
     &                                                              )
        enddo
      enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA
!$acc end kernels
# undef rzetaSA
# undef rzeta2
# undef rzeta
# undef zwrk
!
#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
!        call check_tab2d(rubar,'rubar st_fast_b','uint')
!        call check_tab2d(rvbar,'rvbar st_fast_b','vint')
       call check_tab2d(zeta(:,:,kstp),'zeta st_fast_b','r')
       call check_tab2d(on_u,'on_u st_fast_b','u')
       call check_tab2d(h,'h st_fast_b','r')
       call check_tab2d(rhoS,'rhoS st_fast_b','r')
       call check_tab2d(rhoA,'rhoA st_fast_b','r')
!        call check_tab2d(rubarh,'rubarh st_fast_b','uint')
!        call check_tab2d(rvbarh,'rvbarh st_fast_b','vint')
#endif

! ! KERNEL_6  UFx <= ( DUon, urhs )
! ! KERNEL_6  VFx <= ( DUon, vrhs, pmask )
! ! KERNEL_6  VFe <= ( DVom, vrhs )
! ! KERNEL_6  UFe <= ( DVom, urhs, pmask )
! ! KERNEL_6  rubar <= ( rubar, UFx, UFe )
! ! KERNEL_6  rvbar <= ( rvbar, VFx, VFe )
! ! KERNEL_6  cff <= ( Drhs, fomn, dndx, vrhs, dmde, urhs)
! ! KERNEL_6  UFx <= ( cff, vrhs )
! ! KERNEL_6  VFe <= ( cff, urhs )
! ! KERNEL_6  rubar <= ( rubar, UFx )
! ! KERNEL_6  rvbar <= ( rvbar, VFe )
! ! KERNEL_6  cff <= ( vbar )
! ! KERNEL_6  rubar <= ( rubar, ubar, cff, om_u, on_u )
! ! KERNEL_6  cff <= ( ubar )
! ! KERNEL_6  rvbar <= ( rvbar, vabr, cff, om_v, on_v )
! ! KERNEL_6  rubar <= ( rubar, om_u, on_u )
! ! KERNEL_6  rvbar <= ( rvbar, om_v, on_v )

!$acc kernels default(present)
# ifdef UV_ADV
! !--------------------------------------------------------------------
! ! Compute horizontal advection terms for momentum equations (2D only)
! !-------- ---------- --------- ----- --- -------- --------- --- -----
! !
! ! Centered second order advection scheme
! !
! ! Numerical diffusion of momentum is implicitely added through 3D
! ! forcing of advection in rufrc and rvfrc (i.e., diffusion is
! ! at slow time scale)
! !
! ! NOTE: mathematically necessary (minimal) index ranges for momentum-
! ! flux components are 
! !
! !      UFx(IstrU-1:Iend,Jstr:Jend)   VFx(Istr:Iend+1,JstrV:Jend)
! !      UFe(IstrU:Iend,Jstr:Jend+1)   VFe(Istr,Iend,JstrV-1,Jend)
! !
! ! however, for computational efficiency, these ranges are
! ! unified by suppressing U,V-suffices in order to allow fusion of the
! ! consecutive loops. This leads to slight increase of the redundant
! ! computations near western and southern boundaries in non-periodic
! ! directions.
! !--------------------------------------------------------------------
! !

      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))

          VFx(i+1,j)=0.25*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i+1,j)
#  endif
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          VFe(i,j)=0.25*(DVom(i,j)+DVom(i,j+1))
     &                 *(vrhs(i,j)+vrhs(i,j+1))

          UFe(i,j+1)=0.25*(DVom(i,j+1)+DVom(i-1,j+1))
     &                   *(urhs(i,j+1)+urhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i,j+1)
#  endif
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)

          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo    !--> discard UFx,VFe,UFe,VFx, DUon,DVom
# endif /* UV_ADV */
! !
! !--------------------------------------------------------------------
! ! Compute Coriolis (2D and 3D) term and advective curvilinear metric
! ! terms (2D only).
! !--------------------------------------------------------------------
! !
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=Drhs(i,j)*(
#  ifdef UV_COR
     &                   fomn(i,j)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &          +0.5*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &                -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
#  endif
     &                   )
#  ifdef MRL_WCI
#   if defined CURVGRID && defined UV_ADV
          cff1 = Drhs(i,j)*(
     &    0.5*( dndx(i,j)*(vst2d(i,j)+vst2d(i,j+1))
     &         -dmde(i,j)*(ust2d(i,j)+ust2d(i+1,j)) ))
#   else
          cff1 = 0.0
#   endif
          UFx(i,j)=(cff+cff1)*(vrhs(i,j)+vrhs(i,j+1))
     &                 +cff*(vst2d(i,j)+vst2d(i,j+1))
          VFe(i,j)=(cff+cff1)*(urhs(i,j)+urhs(i+1,j))
     &                 +cff*(ust2d(i,j)+ust2d(i+1,j))
#  else
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
#  endif 
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar(i,j)=rubar(i,j)+0.25*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j)-0.25*(VFe(i,j)+VFe(i,j-1))
        enddo 
      enddo
# endif /* UV_COR */
! !
! !--------------------------------------------------------------------
! ! Linear and/or quadratic bottom stress.
! !--------------------------------------------------------------------
! !
# ifndef BSTRESS_FAST
#  ifdef BBL
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - bustr(i,j)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - bvstr(i,j)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
#  else /* ! BBL */
      if (rdrg2.gt.0.) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
 
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
 
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
#  endif
# endif /* ! BSTRESS_FAST */
! !
! !--------------------------------------------------------------------
! ! Add 2D vortex-force terms combined with advection terms
! !--------------------------------------------------------------------
! !
# ifdef MRL_WCI
      do j=Jstr,Jend
        do i=IstrU,Iend
          vstu = 0.25*( vst2d(i  ,j)+vst2d(i  ,j+1)
     &                 +vst2d(i-1,j)+vst2d(i-1,j+1) )
          dudx = 0.50*( urhs(i+1,j)-urhs(i-1,j  ) )
          dvdx = 0.50*( vrhs(i  ,j)-vrhs(i-1,j  )
     &                 +vrhs(i,j+1)-vrhs(i-1,j+1) )
          rubar(i,j) = rubar(i,j) + 0.5*on_u(i,j)*
     &                      ( Drhs(i-1,j)+Drhs(i,j) )
     &               *( ust2d(i,j)*dudx + vstu*dvdx )
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          ustv = 0.25*( ust2d(i,j  )+ust2d(i+1,j  )
     &                 +ust2d(i,j-1)+ust2d(i+1,j-1) )
          dude = 0.50*( urhs(i,j  )-urhs(i  ,j-1)
     &                 +urhs(i+1,j)-urhs(i+1,j-1) )
          dvde = 0.50*( vrhs(i,j+1)-vrhs(i  ,j-1) )
          rvbar(i,j) = rvbar(i,j) + 0.5*om_v(i,j)*
     &                      ( Drhs(i,j-1)+Drhs(i,j) )
     &               *( ustv*dude + vst2d(i,j)*dvde )
        enddo
      enddo
# endif
!$acc end kernels
!#ifdef NBQ_TBT 
      if (FIRST_FAST_STEP) then 
! !
! !--------------------------------------------------------------------
! ! Since coupling requires that pressure gradient term is computed
! ! using zeta(:,:,kstp) instead of zeta_new(:,:) needed to achieve
! ! numerical stability, apply compensation to shift pressure gradient
! ! terms from "kstp" to "knew": in essense, convert the fist 2D step
! ! from Forward Euler to Forward-Backward].
! !--------------------------------------------------------------------
! !   
#  define zwrk UFx
#  define rzeta  UFe
#  define rzeta2  VFe
#  define rzetaSA VFx
!$acc kernels default(present)
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            zwrk(i,j)=zeta(i,j,knew)-zeta(i,j,kstp)
# if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*(zeta(i,j,knew)+zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*(zeta(i,j,knew)+zeta(i,j,kstp))
# endif
          enddo
        enddo
        cff=0.5*g
        do j=Jstr,Jend
          do i=Istr,Iend
# ifdef M3FAST_BOTH          
            rubarh(i,j)=rubarh(i,j) +cff*on_u(i,j)*( (h(i-1,j)+h(i,j))  
# else            
            rubar(i,j)=rubar(i,j) +cff*on_u(i,j)*( (h(i-1,j)+h(i,j))  
# endif            
     &          *(rzeta(i-1,j)-rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
# endif
     &                                                              )
!
# ifdef M3FAST_BOTH          
            rvbarh(i,j)=rvbarh(i,j) +cff*om_v(i,j)*( (h(i,j-1)+h(i,j))
# else            
            rvbar(i,j)=rvbar(i,j) +cff*om_v(i,j)*( (h(i,j-1)+h(i,j))
# endif            
     &          *(rzeta(i,j-1)-rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(zwrk(i,j-1)-zwrk(i,j)) )
# endif
     &                                                              )
          enddo
        enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA
!$acc end kernels        

#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
!        call check_tab2d(rubarh(:,:),'rubarh st_fast_f','uint')
#endif

# undef rzetaSA
# undef rzeta2
# undef rzeta
# undef zwrk
      endif   !<-- FIRST_FAST_STEP
!#endif  /* NBQ_TBT */
