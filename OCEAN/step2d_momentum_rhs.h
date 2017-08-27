!
!
!======================================================================
! AM4 Backward Step: Computes ubar,vbar at time n+1
!                    first computes rubar,rvbar
!======================================================================
!
!

! Compute pressure-gradient terms  NOTE that "rubar" and "rvbar"
!-------- -------- -------- -----  are computed within the same
! fused loop despite the fact that their normal index ranges are
! different. Fusing loops causes redundant computation of one
! column of "rubar" on the western physical boundary and one row
! of "rvbar" on the southern, but, at the same time it allows to
! share references to array elements (i,j) which results in an
! increase of computational density by almost a factor of 1.5
! resulting in overall more efficient code pipelined in 26 cycles
! (61% of peak speed) on R10000 vs. 16+16 cycles of separate loop
! version for the case when both CPP switches below are defined.
!
      cff=0.5*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(zwrk(i-1,j)-zwrk(i,j)))
#endif
#ifdef MRL_WCI
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i-1,j) )
#endif
#if defined POT_TIDES && !defined SOLVE3D
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                   *( Ptide(i,j)-Ptide(i-1,j) )
#endif
     &                                                              )    
! 
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)
#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(zwrk(i,j-1)-zwrk(i,j)))
#endif
#ifdef MRL_WCI
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i,j-1) )
#endif
#if defined POT_TIDES && !defined SOLVE3D
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                   *( Ptide(i,j)-Ptide(i,j-1) )
#endif
     &                                                              )
        enddo
      enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

#undef rzetaSA
#undef rzeta2
#undef rzeta
#undef zwrk
!
!========================================================================
! Compute horizontal advection terms for momentum equations (2D only)
!-------- ---------- --------- ----- --- -------- --------- --- -----
! NOTE: mathematically necessary (minimal) index ranges for momentum-
! flux components are 
!
!      UFx(IstrU-1:Iend,Jstr:Jend)   VFx(Istr:Iend+1,JstrV:Jend)
!      UFe(IstrU:Iend,Jstr:Jend+1)   VFe(Istr,Iend,JstrV-1,Jend)
!
! however, for the purpose computational efficiency, these ranges are
! unified by suppressing U,V-suffices in order to allow fusion of the
! consecutive loops. This leads to slight increase of the redundant
! computations near western and southern boundaries in non-periodic
! directions. 
!========================================================================
!
#ifdef UV_ADV
!
# if defined SOLVE3D || !defined M2_HADV_UP3
!-----------------------------------------------------------------------
! Centered Second order advection scheme
!
! Numerical diffusion of momentum is implicitely added through 3D
! forcing of advection in rufrc and rvfrc (i.e., diffusion is
! at slow time scale)
!-----------------------------------------------------------------------
      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))

          VFx(i+1,j)=0.25*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i+1,j)
#  endif
#  ifdef WET_DRY0
     &                                 *pmask_wet(i+1,j)
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
#  ifdef WET_DRY0
     &                                 *pmask_wet(i,j+1)
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
      enddo !--> discard UFx,VFe,UFe,VFx, DUon,DVom

# else   
!-----------------------------------------------------------------------
!  UP3 (QUICK) spatial advection scheme as in 3D part
!
!  in this case, there is no numerical requirement for explicit viscosity
!-----------------------------------------------------------------------
#  define uxx wrk1
#  define Huxx wrk2
!
#  ifdef EW_PERIODIC
#   define IU_EXT_RANGE IstrU-1,Iend+1
#  else
#   ifdef MPI
        if (WEST_INTER) then
          imin=IstrU-1
        else
          imin=max(IstrU-1,2)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
#    define IU_EXT_RANGE imin,imax
#   else
#    define IU_EXT_RANGE max(IstrU-1,2),min(Iend+1,Lm)
#   endif
#  endif

        do j=Jstr,Jend
          do i=IU_EXT_RANGE
            uxx(i,j)=urhs(i-1,j)-2.*urhs(i,j)+urhs(i+1,j)
            Huxx(i,j)=Duon(i-1,j)-2.*Duon(i,j)+Duon(i+1,j)
          enddo
        enddo
#   undef IU_EXT_RANGE
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=Jstr,Jend
            uxx(1,j)=uxx(2,j)
            Huxx(1,j)=Huxx(2,j)
          enddo
        endif
        if (EASTERN_EDGE) then
#   ifdef MPI        
          do j=Jstr,Jend
            uxx(Lmmpi+1,j)=uxx(Lmmpi,j)
            Huxx(Lmmpi+1,j)=Huxx(Lmmpi,j)
          enddo
#   else
          do j=Jstr,Jend
            uxx(Lm+1,j)=uxx(Lm,j)
            Huxx(Lm+1,j)=Huxx(Lm,j)
          enddo
#   endif          
        endif
#  endif

        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cffX=urhs(i,j)+urhs(i+1,j)
            if (cffX.gt.0.) then
              curvX=uxx(i,j)
            else
              curvX=uxx(i+1,j)
            endif
            UFx(i,j)=0.25*(cffX+gamma*curvX)*( Duon(i,j)+
     &               Duon(i+1,j)-0.125*(Huxx(i,j)+Huxx(i+1,j)))
          enddo
        enddo
#  undef Huxx
#  undef uxx
#  define vee wrk1
#  define Hvee wrk2
      
#  ifdef NS_PERIODIC
#   define JV_EXT_RANGE JstrV-1,Jend+1
#  else
#   ifdef MPI
        if (SOUTH_INTER) then
          jmin=JstrV-1
        else
          jmin=max(JstrV-1,2)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
#    define JV_EXT_RANGE jmin,jmax
#   else
#    define JV_EXT_RANGE max(JstrV-1,2),min(Jend+1,Mm)
#   endif
#  endif

        do j=JV_EXT_RANGE
          do i=Istr,Iend
            vee(i,j) =vrhs(i,j-1)-2.*vrhs(i,j)+vrhs(i,j+1)
            Hvee(i,j)=Dvom(i,j-1)-2.*Dvom(i,j)+Dvom(i,j+1)
          enddo
        enddo
#   undef JV_EXT_RANGE
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=Istr,Iend
            vee(i,1)=vee(i,2)
            Hvee(i,1)=Hvee(i,2)
          enddo
        endif
        if (NORTHERN_EDGE) then
#   ifdef MPI
          do i=Istr,Iend
            vee(i,Mmmpi+1)=vee(i,Mmmpi)
            Hvee(i,Mmmpi+1)=Hvee(i,Mmmpi)
          enddo
#   else        
          do i=Istr,Iend
            vee(i,Mm+1)=vee(i,Mm)
            Hvee(i,Mm+1)=Hvee(i,Mm)
          enddo
#   endif          
        endif
#  endif
        do j=JstrV-1,Jend
          do i=Istr,Iend
            cffE=vrhs(i,j)+vrhs(i,j+1)
            if (cffE.gt.0.) then
              curvE=vee(i,j)
            else
              curvE=vee(i,j+1)
            endif
            VFe(i,j)=0.25*(cffE+gamma*curvE)*( Dvom(i,j)+
     &               Dvom(i,j+1)-0.125*(Hvee(i,j)+Hvee(i,j+1)))
          enddo
        enddo
#  undef Hvee
#  undef vee
#  define uee wrk1
#  define Hvxx wrk2

#  ifdef NS_PERIODIC
#   define JU_EXT_RANGE Jstr-1,Jend+1
#  else
#   ifdef MPI
        if (SOUTH_INTER) then
          jmin=Jstr-1
        else
          jmin=max(Jstr-1,1)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
#    define JU_EXT_RANGE jmin,jmax
#   else
#    define JU_EXT_RANGE max(Jstr-1,1),min(Jend+1,Mm)
#   endif
#  endif

        do j=JU_EXT_RANGE
          do i=IstrU,Iend
            uee(i,j)=urhs(i,j-1)-2.*urhs(i,j)+urhs(i,j+1)
          enddo
        enddo
#   undef JU_EXT_RANGE
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=IstrU,Iend
            uee(i,0)=uee(i,1)
          enddo
        endif
        if (NORTHERN_EDGE) then
#   ifdef MPI
          do i=IstrU,Iend
            uee(i,Mmmpi+1)=uee(i,Mmmpi)
          enddo
#   else        
          do i=IstrU,Iend
            uee(i,Mm+1)=uee(i,Mm)
          enddo
#   endif          
        endif
#  endif
        do j=Jstr,Jend+1
          do i=IstrU-1,Iend
           Hvxx(i,j)=Dvom(i-1,j)-2.*Dvom(i,j)+Dvom(i+1,j)
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=IstrU,Iend
            cffX=urhs(i,j)+urhs(i,j-1)
            cffE=Dvom(i,j)+Dvom(i-1,j)
            if (cffE.gt.0.) then
              curvX=uee(i,j-1)
            else
              curvX=uee(i,j)
            endif
            UFe(i,j)=0.25*(cffX+gamma*curvX)*(cffE-0.125*(
     &                             Hvxx(i,j)+Hvxx(i-1,j) ))
          enddo
        enddo
#  undef Hvxx
#  undef uee
#  define vxx wrk1
#  define Huee wrk2

#  ifdef EW_PERIODIC 
#   define IV_EXT_RANGE Istr-1,Iend+1
#  else
#   ifdef MPI
        if (WEST_INTER) then
          imin=Istr-1
        else
          imin=max(Istr-1,1)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
#    define IV_EXT_RANGE imin,imax
#   else
#    define IV_EXT_RANGE max(Istr-1,1),min(Iend+1,Lm)
#   endif
#  endif

        do j=JstrV,Jend
          do i=IV_EXT_RANGE
            vxx(i,j)=vrhs(i-1,j)-2.*vrhs(i,j)+vrhs(i+1,j)
          enddo
        enddo
#   undef IV_EXT_RANGE
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=JstrV,Jend
            vxx(0,j)=vxx(1,j)
          enddo
        endif
        if (EASTERN_EDGE) then
#   ifdef MPI
          do j=JstrV,Jend
            vxx(Lmmpi+1,j)=vxx(Lmmpi,j)
          enddo
#   else        
          do j=JstrV,Jend
            vxx(Lm+1,j)=vxx(Lm,j)
          enddo
#   endif          
        endif
        
#  endif
        do j=JstrV-1,Jend
          do i=Istr,Iend+1
           Huee(i,j)=Duon(i,j-1)-2.*Duon(i,j)+Duon(i,j+1)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend+1
            cffE=vrhs(i,j)+vrhs(i-1,j)
            cffX=Duon(i,j)+Duon(i,j-1)
            if (cffX.gt.0.) then
              curvE=vxx(i-1,j)
            else
              curvE=vxx(i,j)
            endif
            VFx(i,j)=0.25*(cffE+gamma*curvE)*(cffX-0.125*(
     &                             Huee(i,j)+Huee(i,j-1) ))
          enddo
        enddo
#undef Huee
#undef vxx
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j)-UFx(i,j  )+UFx(i-1,j)
     &                           -UFe(i,j+1)+UFe(i  ,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j  )
     &                           -VFe(i  ,j)+VFe(i,j-1)
          enddo
        enddo
# endif /* !SOLVE3D && M2_HADV_UP3 */
#endif /* UV_ADV */
!
!-----------------------------------------------------------------------
! Compute Coriolis (2D and 3D) term and advective curvilinear metric
! terms (2D only).
!-----------------------------------------------------------------------
!
#if defined UV_COR || (defined CURVGRID && defined UV_ADV)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=Drhs(i,j)*(
# ifdef UV_COR
     &                   fomn(i,j)
# endif
# if (defined CURVGRID && defined UV_ADV)
     &          +0.5*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &                -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
# endif
     &                   )
# ifdef MRL_WCI
#  if defined CURVGRID && defined UV_ADV
          cff1 = Drhs(i,j)*(
     &    0.5*( dndx(i,j)*(vst2d(i,j)+vst2d(i,j+1))
     &          -dmde(i,j)*(ust2d(i,j)+ust2d(i+1,j)) ))
#  else
          cff1 = 0.0
#  endif
          UFx(i,j)=(cff+cff1)*(vrhs(i,j)+vrhs(i,j+1))
     &                 +cff*(vst2d(i,j)+vst2d(i,j+1))
          VFe(i,j)=(cff+cff1)*(urhs(i,j)+urhs(i+1,j))
     &                 +cff*(ust2d(i,j)+ust2d(i+1,j))
# else
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
# endif 
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar(i,j)=rubar(i,j)+0.25*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j) -0.25*(VFe(i,j)+VFe(i,j-1))
        enddo 
      enddo
#endif
!
!-----------------------------------------------------------------------
! Compute horizontal viscous stress terms (2D only).
!-----------------------------------------------------------------------
!
#ifdef SOLVE3D
# undef UV_VIS2
#endif
#ifdef UV_VIS2
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          cff=2.*Drhs(i,j)*visc2_r(i,j)
          UFx(i,j)=cff*(ubar(i+1,j,kstp)-ubar(i,j,kstp))
     &                                 *pm(i,j)*on_r(i,j)
          VFe(i,j)=cff*(vbar(i,j+1,kstp)-vbar(i,j,kstp))
     &                                 *pn(i,j)*om_r(i,j)

          cff1=0.0625*visc2_p(i+1,j+1)*( Drhs(i,j)
     &       +Drhs(i+1,j)+Drhs(i,j+1)+Drhs(i+1,j+1) )*(
     &          (pn(i+1,j+1)+pn(i,j+1)+pn(i+1,j)+pn(i,j))
     &             *(ubar(i+1,j+1,kstp)-ubar(i+1,j,kstp))
     &         +(pm(i+1,j+1)+pm(i,j+1)+pm(i+1,j)+pm(i,j))
     &             *(vbar(i+1,j+1,kstp)-vbar(i,j+1,kstp))
     &                                                  )
# ifdef MASKING
     &                     *pmask(i+1,j+1)
# endif
          UFe(i+1,j+1)=cff1*om_p(i+1,j+1)
          VFx(i+1,j+1)=cff1*on_p(i+1,j+1)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)+UFx(i,j)-UFx(i-1,j)
     &                         +UFe(i,j+1)-UFe(i,j)

          rvbar(i,j)=rvbar(i,j)+VFx(i+1,j)-VFx(i,j)
     &                         +VFe(i,j)-VFe(i,j-1)
        enddo
      enddo
#endif /* UV_VIS2 */
!
!-----------------------------------------------------------------------
! Linear and/or quadratic bottom stress.
!-----------------------------------------------------------------------
!
# ifdef BBL
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
# else
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
# endif
!
!-----------------------------------------------------------------------
! Add 2D vortex-force terms combined with advection terms
!-----------------------------------------------------------------------
!
#ifdef MRL_WCI
      do j=Jstr,Jend
        do i=IstrU,Iend
          vstu = 0.25*( vst2d(i,j) + vst2d(i,j+1)
     &                 +vst2d(i-1,j)+vst2d(i-1,j+1) )
          dudx = 0.5*( urhs(i+1,j)-urhs(i-1,j) )
          dvdx = 0.5*( vrhs(i,j) - vrhs(i-1,j)
     &                 +vrhs(i,j+1) - vrhs(i-1,j+1) )
          rubar(i,j) = rubar(i,j) + 0.5*om_u(i,j)*
     &                      ( Drhs(i-1,j)+Drhs(i,j) )
     &               *( ust2d(i,j)*dudx + vstu*dvdx )
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          ustv = 0.25*( ust2d(i,j) + ust2d(i+1,j)
     &                 +ust2d(i,j-1)+ust2d(i+1,j-1) )
          dude = 0.5*( urhs(i,j) - urhs(i,j-1)
     &                 +urhs(i+1,j) - urhs(i+1,j-1) )
ccc          dvde = 0.5*( vrhs(i,j-1)-vrhs(i,j+1) )
          dvde = 0.5*( vrhs(i,j+1)-vrhs(i,j-1) )
          rvbar(i,j) = rvbar(i,j) + 0.5*on_v(i,j)*
     &                      ( Drhs(i,j-1)+Drhs(i,j) )
     &               *( ustv*dude + vst2d(i,j)*dvde )
        enddo
      enddo
#endif
