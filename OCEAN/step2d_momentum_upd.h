!
!-----------------------------------------------------------------------
! Perform time step for the 2D momentum equations. 

! Also compute fast-time averaged barotropic mass fluxes. 
! Doing so on the fly yields a more computationally dense code and 
! eliminates repeated multiplication by Dnew (since mass fluxes are 
! actually available as volatile variables DUnew,DVnew at this moment. 
! However, as a result of this arrangement, a special code is needed 
! to compute fast-time averages along the physical boundaries, which is 
! done below.
!-----------------------------------------------------------------------
!

#define Dstp DUon
 
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
#ifdef NBQ_MASS
          Dstp(i,j)=zeta(i,j,kstp)
#else
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
#endif
        enddo
      enddo

      cff=0.5*dtfast
#ifdef SOLVE3D
      cff1=0.5*weight(1,iif)
#else
      cff2=2.*dtfast
#endif

      do j=Jstr,Jend
        do i=IstrU,Iend
!          DUnew=( (Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)
!     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!#ifdef SOLVE3D
!# ifdef NBQ
!     &              *(rubar(i,j)+rufrc(i,j)+rubar_nbq(i,j))
!# else
!     &                             *(rubar(i,j)+rufrc(i,j))
!# endif
!#else
!     &                          *rubar(i,j)+cff2*sustr(i,j)
!# ifdef MRL_WCI
!     &                                    +cff2*brk2dx(i,j)
!#  ifdef WAVE_STREAMING
!     &                                    +cff2*frc2dx(i,j)
!#  endif
!# endif
!#endif
!     &                                                    )
!#ifdef MASKING
!     &                                         *umask(i,j)
!#endif

          DUnew=DU_nbq(i,j) *2.
#ifdef MASKING
     &                                         *umask(i,j)
#endif

#ifdef WET_DRY
          cff1_WD=ABS(ABS(umask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DUnew)*umask_wet(i,j)
          umask_wet(i,j)=0.5*umask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.-cff1_WD)
          DUnew=DUnew*umask_wet(i,j)
# ifdef MRL_WCI
          ust2d(i,j)=ust2d(i,j)*umask_wet(i,j)
# endif
#endif
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
#ifdef SOLVE3D
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*on_u(i,j)*( DUnew
# ifdef MRL_WCI
     &                 +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
# endif
     &                                                   )
#endif
        enddo
      enddo   
  
      do j=JstrV,Jend
        do i=Istr,Iend
!          DVnew=( (Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)
!     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!#ifdef SOLVE3D
!# ifdef NBQ
!     &              *(rvbar(i,j)+rvfrc(i,j)+rvbar_nbq(i,j))
!# else
!     &                             *(rvbar(i,j)+rvfrc(i,j))
!# endif
!#else
!     &                          *rvbar(i,j)+cff2*svstr(i,j)
!# ifdef MRL_WCI
!     &                                    +cff2*brk2de(i,j)
!#  ifdef WAVE_STREAMING
!     &                                    +cff2*frc2de(i,j)
!#  endif
!# endif
!#endif
!     &                                                    )
!#ifdef MASKING
!     &                                         *vmask(i,j)
!#endif

          DVnew=DV_nbq(i,j) *2.
#ifdef MASKING
     &                                         *vmask(i,j)
#endif

#ifdef WET_DRY
          cff1_WD=ABS(ABS(vmask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DVnew)*vmask_wet(i,j)
          vmask_wet(i,j)=0.5*vmask_wet(i,j)*cff1_WD
     &                        +cff2_WD*(1.-cff1_WD)
          DVnew=DVnew*vmask_wet(i,j)
# ifdef MRL_WCI
          vst2d(i,j)=vst2d(i,j)*vmask_wet(i,j)
# endif
#endif
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
#ifdef SOLVE3D
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*om_v(i,j)*(DVnew
# ifdef MRL_WCI
     &                 +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
# endif
     &                                                   )
#endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
!  Set 2D Momemtum nudging
!-----------------------------------------------------------------------
!
#if defined M2NUDGING && defined M2CLIMATOLOGY

# ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew),ubzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew),vbzon)
        endif
      endif
      if (iic.eq.ntstart) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubclm(START_2D_ARRAY),ubclmzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbclm(START_2D_ARRAY),vbclmzon)
        endif
      endif
# endif /* ZONAL_NUDGING */

      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef ZONAL_NUDGING        
          DUnew = dtfast*M2nudgcof(i,j)*(ubclmzon(j)-ubzon(j))
# else          
          DUnew = dtfast*M2nudgcof(i,j)*(ubclm(i,j)-ubar(i,j,knew))
# endif
# ifdef MASKING
     &                 * umask(i,j)
# endif 
# ifdef WET_DRY
     &                 * umask_wet(i,j)
# endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
# ifdef SOLVE3D
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# endif
        enddo
      enddo
      
      do j=JstrV,Jend
        do i=Istr,Iend
# if defined ZONAL_NUDGING  
          DVnew = dtfast*M2nudgcof(i,j)*(vbclmzon(j)-vbzon(j)) 
# else
          DVnew = dtfast*M2nudgcof(i,j)*(vbclm(i,j)-vbar(i,j,knew))    
# endif      
# ifdef MASKING
     &                 * vmask(i,j)
# endif 
# ifdef WET_DRY
     &                 * vmask_wet(i,j)
# endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
# ifdef SOLVE3D
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# endif

        enddo
      enddo
#endif
!
!-----------------------------------------------------------------------
!  Body force for the Internal Tide test case
!-----------------------------------------------------------------------
!
#if defined INTERNAL || defined BODYTIDE
      omega=2.*pi/(12.4*3600.)
      U0=0.02
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew = dtfast*omega*U0*cos(omega*time)
# ifdef MASKING
     &                 * umask(i,j)
# endif
# ifdef WET_DRY
     &                 * umask_wet(i,j)
# endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
#ifdef SOLVE3D
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
#endif
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          DVnew = dtfast*0.5*(f(i,j)+f(i,j-1))*
     &                   U0*sin(omega*time)
# ifdef MASKING
     &                 * vmask(i,j)
# endif
# ifdef WET_DRY
     &                 * vmask_wet(i,j)
# endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
#ifdef SOLVE3D
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
#endif

        enddo
      enddo
#endif  

!
!-----------------------------------------------------------------------
! Set boundary conditions and compute integral mass flux accross
! all open boundaries, if any.
!-----------------------------------------------------------------------
!
      call u2dbc_tile (Istr,Iend,Jstr,Jend, UFx) 
      call v2dbc_tile (Istr,Iend,Jstr,Jend, UFx)

          
#ifdef OBC_VOLCONS
      call obc_flux_tile (Istr,Iend,Jstr,Jend)
#endif
!
!-----------------------------------------------------------------------
! Compute fast-time averaged barotropic mass fluxes along physical
! boundaries.
!-----------------------------------------------------------------------
!
#ifdef SOLVE3D
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=Jstr-1,JendR
#  ifndef NBQ_MASS
          Dnew(Istr-1,j)=h(Istr-1,j)+zeta(Istr-1,j,knew)
#  else
          Dnew(Istr-1,j)=zeta(Istr-1,j,knew)
#  endif
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr-1,JendR
#  ifndef NBQ_MASS
          Dnew(Iend+1,j)=h(Iend+1,j)+zeta(Iend+1,j,knew)
#  else
          Dnew(Iend+1,j)=zeta(Iend+1,j,knew)
#  endif
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=Istr-1,IendR
#  ifndef NBQ_MASS
          Dnew(i,Jstr-1)=h(i,Jstr-1)+zeta(i,Jstr-1,knew)
#  else
          Dnew(i,Jstr-1)=zeta(i,Jstr-1,knew)
#  endif
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=Istr-1,IendR
#  ifndef NBQ_MASS
          Dnew(i,Jend+1)=h(i,Jend+1)+zeta(i,Jend+1,knew)
#  else
          Dnew(i,Jend+1)=zeta(i,Jend+1,knew)
#  endif
        enddo
      endif
# endif
      cff1=0.5*weight(1,iif)
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)
     &         +cff1*(Dnew(IstrU-1,j)
     &         +Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
# ifdef MRL_WCI
     &                                             +ust2d(IstrU-1,j)
# endif
     &                                             )*on_u(IstrU-1,j)
        enddo
        do j=JstrV,Jend
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)
     &       +cff1*(Dnew(Istr-1,j)
     &       +Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
# ifdef MRL_WCI
     &                                              +vst2d(Istr-1,j)
# endif
     &                                              )*om_v(Istr-1,j)
        enddo
      endif
        
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)
     &            +cff1*( Dnew(Iend+1,j)
     &            +Dnew(Iend,j) )*(ubar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                              +ust2d(Iend+1,j)
# endif
     &                                              )*on_u(Iend+1,j)
        enddo
        do j=JstrV,Jend
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)
     &        +cff1*( Dnew(Iend+1,j)
     &        +Dnew(Iend+1,j-1) )*(vbar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                              +vst2d(Iend+1,j)
# endif
     &                                              )*om_v(Iend+1,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=IstrU,Iend
          DU_avg1(i,Jstr-1,nnew)=DU_avg1(i,Jstr-1,nnew)
     &        +cff1*( Dnew(i,Jstr-1)
     &        +Dnew(i-1,Jstr-1) )*(ubar(i,Jstr-1,knew)
# ifdef MRL_WCI
     &                                              +ust2d(i,Jstr-1)
# endif
     &                                              )*on_u(i,Jstr-1)
        enddo
        do i=IstrR,IendR
          DV_avg1(i,JstrV-1,nnew)=DV_avg1(i,JstrV-1,nnew)
     &         +cff1*(Dnew(i,JstrV-1)
     &         +Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
# ifdef MRL_WCI
     &                                              +vst2d(i,JstrV-1)
# endif
     &                                              )*om_v(i,JstrV-1)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=IstrU,Iend
          DU_avg1(i,Jend+1,nnew)=DU_avg1(i,Jend+1,nnew)
     &        +cff1*( Dnew(i,Jend+1)
     &        +Dnew(i-1,Jend+1) )*(ubar(i,Jend+1,knew)
# ifdef MRL_WCI
     &                                               +ust2d(i,Jend+1)
# endif
     &                                               )*on_u(i,Jend+1)
        enddo
        do i=IstrR,IendR
          DV_avg1(i,Jend+1,nnew)=DV_avg1(i,Jend+1,nnew)
     &            +cff1*( Dnew(i,Jend+1)
     &            +Dnew(i,Jend) )*(vbar(i,Jend+1,knew)
# ifdef MRL_WCI
     &                                               +vst2d(i,Jend+1)
# endif
     &                                               )*om_v(i,Jend+1)
        enddo
      endif
# endif 
#endif
!
!-----------------------------------------------------------------------
! Apply point sources for river runoff simulations
!-----------------------------------------------------------------------
!
#ifdef PSOURCE
      do is=1,Nsrc 
# ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
# else
        i=Isrc(is)
        j=Jsrc(is)
# endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            ubar(i,j,knew)=2.*Qbar(is)/( on_u(i,j)
     &                   *(Dnew(i-1,j)+Dnew(i,j)) )
# ifdef SOLVE3D
            DU_avg1(i,j,nnew)=Qbar(is)
# endif
          else
            vbar(i,j,knew)=2.*Qbar(is)/( om_v(i,j)
     &                   *(Dnew(i,j-1)+Dnew(i,j)) )
# ifdef SOLVE3D
            DV_avg1(i,j,nnew)=Qbar(is)
# endif
          endif
        endif
      enddo
#endif
!
!-----------------------------------------------------------------------
!  Diagnostics
!-----------------------------------------------------------------------
!
#ifndef SOLVE3D
      call diag_tile (Istr,Iend,Jstr,Jend, UFx,UFe)
#endif
     
