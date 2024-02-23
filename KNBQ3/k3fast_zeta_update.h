! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_zeta_update.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
#  if defined K3FAST_ZETAW || defined KNHINT_ZETAW
! !
! !********************************
! ! Computes surface velocities
! !********************************
! ! 
!$acc kernels if(compute_on_device) default(present)
#   ifndef K3FAST_SPDUP
      if (IstrU.le.Iend) then
        do j=Jstr,Jend
          do i=IstrU-1,Iend+1     
            usurf_nbq(i,j)=(qdmu_nbq(i,j,N)
     &                     *2./(Hz(i,j,N)+Hz(i-1,j,N))

#   ifdef MRL_WCI
     &                     +ust(i,j,N)
#    ifdef WET_DRY
     &                      *umask_wet(i,j)
#    endif
#   endif
     &                     )
#   ifdef MASKING
     &                     *umask(i,j) 
#   endif
          enddo 
        enddo 
      endif

      if (JstrV.le.Jend) then
        do j=JstrV-1,Jend+1
          do i=Istr,Iend     
            vsurf_nbq(i,j)=(qdmv_nbq(i,j,N)
     &                     *2./(Hz(i,j,N)+Hz(i,j-1,N))
#   ifdef MRL_WCI
     &                     +vst(i,j,N)
#    ifdef WET_DRY
     &                      *vmask_wet(i,j)
#    endif
#   endif
     &                      )
#   ifdef MASKING
     &                     *vmask(i,j) 
#   endif
          enddo
        enddo 
      endif
#  endif  /* K3FAST_SPDUP  */

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          wsurf_nbq(i,j)=(qdmw_nbq(i,j,N)         
#   ifdef NBQ_MASS     
     &                       /(1.+rho_grd(i,j,N))
#   endif
     &                       *Hzw_nbq_inv(i,j,N)
!         wsurf_nbq(i,j)=(We(i,j,N)*pm(i,j)*pn(i,j)   ! CAUTION
#  ifdef MRL_WCI
     &                     +wst(i,j,N)*rmask_wet(i,j)
#  endif
     &                    )
#  ifdef MASKING
     &                       *rmask(i,j) 
#  endif
        enddo
      enddo 
!$acc end kernels   
! !
! !********************************
! ! ATTENTION replaced by zeta 3pts  
! !        FRANCIS
! !********************************
! !
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      if (IstrU.le.Iend) then
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          usurf_nbq(START_2D_ARRAY))
      endif
      if (JstrV.le.Jend) then
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          vsurf_nbq(START_2D_ARRAY))
      endif
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        wsurf_nbq(START_2D_ARRAY))
#  endif
! !
! !********************************
! ! Apply point sources for 
! !   river runoff simulations
! !********************************
! !
#  ifdef PSOURCE
      do is=1,Nsrc 
#   ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#   else
        i=Isrc(is)
        j=Jsrc(is)
#   endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            urhs(i,j)=2.*Qbar(is)/( on_u(i,j)
     &                             *(Drhs(i-1,j)+Drhs(i,j)) )
            DUon(i,j)=Qbar(is)
            usurf_nbq(i,j)=2.*Qsrc(is,N)
     &                     /( on_u(i,j)*(Hz(i,j,N)+Hz(i-1,j,N)) )
          else
            vrhs(i,j)=2.*Qbar(is)/( om_v(i,j)
     &                             *(Drhs(i,j-1)+Drhs(i,j)) )
            DVom(i,j)=Qbar(is)
            vsurf_nbq(i,j)=2.*Qsrc(is,N)
     &                     /( om_v(i,j)*(Hz(i,j,N)+Hz(i,j-1,N)) )
          endif
        endif 
      enddo
#  endif  /* PSOURCE */
! !
! !********************************
! ! Advance zeta at m+1 with 
! !  kinematic condition
! !********************************
! ! 
!$acc kernels if(compute_on_device) default(present)
#   ifndef K3FAST_SPDUP
#    define zab3 UFx
      if (FIRST_FAST_STEP) then
        cff1 = 1.0
        cff2 = 0.0
        cff3 = 0.0
      elseif (FIRST_FAST_STEP+1) then  ! AB2
        cff1 = 1.5
        cff2 =-0.5
        cff3 = 0.0
      else                             ! AB3
        cff1 = 1.5+mybeta
        cff2 =-2.0*mybeta-0.5
        cff3 = mybeta
      endif

      do j=JstrV-2,Jend+1
        do i=IstrU-2,Iend+1
          zab3(i,j) =   cff1 * zeta(i,j,kstp)
     &                + cff2 * zeta(i,j,kbak)
     &                + cff3 * zeta(i,j,kold)
        enddo
      enddo
#   endif /* !K3FAST_SPDUP */
#   ifdef KNHINT_ZETAW
      if (iic==1.and.iif==1) dzeta_nbq=0.
#   endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
#   ifdef K3FAST_ZETAW
	  zeta(i,j,knew)=(zeta(i,j,kstp) + dtfast*( wsurf_nbq(i,j)
#   elif defined KNHINT_ZETAW
	  zeta(i,j,knew)=wsurf_nbq(i,j)
#   endif
#   ifndef K3FAST_SPDUP
     &                     -0.5*(usurf_nbq(i  ,j)
     &                                     *(zab3(i  ,j)
     &                                      -zab3(i-1,j))*pm_u(i,j)
     &                          +usurf_nbq(i+1,j)
     &                                     *(zab3(i+1,j)
     &                                      -zab3(i  ,j))*pm_u(i+1,j)
     &                          )
     &                     -0.5*(vsurf_nbq(i  ,j)
     &                                     *(zab3(i,j  )
     &                                      -zab3(i,j-1))*pn_v(i,j)
     &                          +vsurf_nbq(i,j+1)
     &                                     *(zab3(i,j+1)
     &                                      -zab3(i,j  ))*pn_v(i,j+1)
#   ifdef K3FAST_ZETAW
     &                          )))
#   elif defined KNHINT_ZETAW
     &                          )-dzeta_nbq(i,j)
#   endif
#    endif    /* K3FAST_SPDUP */
#    ifdef MASKING
     &                                                  *rmask(i,j)
#    endif
        enddo
      enddo
!$acc end kernels
# endif /* ! K3FAST_ZETAW && ! KNHINT_ZETAW*/

# ifndef K3FAST_ZETAW 
! !
! !********************************
! !  Advance zeta at m+1 from continuity 
! !   for the hydrostatic case 
! !********************************
! !
! !  First, apply point sources
! !
#  ifdef PSOURCE
      do is=1,Nsrc 
#   ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#   else
        i=Isrc(is)
        j=Jsrc(is)
#   endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            urhs(i,j)=2.*Qbar(is)/( on_u(i,j)
     &                             *(Drhs(i-1,j)+Drhs(i,j)) )
            DUon(i,j)=Qbar(is)
          else
            vrhs(i,j)=2.*Qbar(is)/( om_v(i,j)
     &                             *(Drhs(i,j-1)+Drhs(i,j)) )
            DVom(i,j)=Qbar(is)
          endif
        endif
      enddo
#  endif
!$acc kernels if(compute_on_device) default(present)
#  ifdef KNHINT_ZETAW
!     if (iic==1.and.iif==1) zeta=0.
#  endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
#  ifndef KNHINT_ZETAW
          zeta(i,j,knew)=(zeta(i,j,kstp)
#  else
          cff=(zeta(i,j,kstp)
#  endif
#  ifdef MVB
     &                    -(dh_mvb(i,j,knew2)-dh_mvb(i,j,kstp2))
#  endif
     &                                   + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1)))
#  ifdef MASKING
     &                                                   *rmask(i,j)
#  endif
#  ifdef K3FAST_ZETADISS
     &              + diss_zta*( zeta(i+1,j,kstp)-2.*zeta(i,j,kstp)
     &                       +zeta(i-1,j,kstp)) /2.*pm(i,j)**2
#  endif
#  ifdef KNHINT_ZETAW
           dzeta_nbq(i,j)=(cff-zeta(i,j,kstp))/dtfast
           zeta(i,j,knew)=cff -dtfast*zeta(i,j,knew)
#  endif

        enddo
      enddo
!$acc end kernels      
# endif /* !K3FAST_ZETAW */
! !
! !********************************
! !  Add nudging terms
! !********************************
! !
# ifdef ZNUDGING
#  ifdef ZONAL_NUDGING
      if (FIRST_TIME_STEP .or. mod(iic,10).eq.0) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew),zetazon)
        endif
      endif
      if (FIRST_TIME_STEP) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ssh(START_2D_ARRAY),sshzon)
        endif
      endif
#  endif  /* ZONAL_NUDGING */
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=zeta(i,j,knew) + dtfast*Znudgcof(i,j)
#  ifdef ZONAL_NUDGING
     &                                 *(sshzon(j)-zetazon(j))
#  else
     &                              *(ssh(i,j)-zeta(i,j,knew))
#  endif /* ZONAL_NUDGING */
#  ifdef MASKING
     &                                             *rmask(i,j)
#  endif
        enddo
      enddo
# endif /* ZNUDGING */
! !
! !********************************
! ! Compute wet/dry masks
! !********************************
! !
! First: modify new free-surface to ensure that depth 
!        is > Dcrit in masked cells.
!
# if defined WET_DRY && defined MASKING
!$acc kernels if(compute_on_device) default(present)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=0.5+SIGN(0.5,Dcrit(i,j)-h(i,j))
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   cff*(Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
        enddo
      enddo 
!$acc end kernels      
# endif
!
! Then compute wet/dry masks
!
# ifdef WET_DRY
      call wetdry_tile (Istr,Iend,Jstr,Jend)
# endif
! !
! !********************************
! !  Set boundary conditions 
! !   for the free-surface
! !********************************
! !
      call zetabc_tile (Istr,Iend,Jstr,Jend)
# if !defined NBQ_HZCORRECT
      if (LAST_FAST_STEP) then
!$acc update host( zeta )        
      endif
# endif      
! !
! !********************************
! !  Perform exchanges
! !   CAUTION: over 1 point only
! !********************************
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
# endif
# if !defined NBQ_HZCORRECT
      if (LAST_FAST_STEP) then
!$acc update host( zeta )        
      endif
# endif   
! !
! !********************************
! !  Compute time averaged fields 
! !   over all short timesteps.
! !
! ! Reset/initialise arrays for averaged fields during the first
! ! barotropic time step; Accumulate averages after that. Include
! ! physical boundary points, but not periodic ghost points or
! ! computation  MPI computational margins.
! !********************************
! !
#ifdef K3FAST_AVG_CLASSIC
# ifdef SOLVE3D
      cff1=weight(1,iif)
      cff2=weight(2,iif)
!$acc kernels if(compute_on_device) default(present) async(1)
      if (FIRST_2D_STEP) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
            DU_avg2(i,j)=cff2*DUon(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
          enddo
        enddo
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
!$acc end kernels      
# else      
!$acc kernels if(compute_on_device) default(present) async(1)
# endif
# endif /* K3FAST_AVG_CLASSIC */

# ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
      call check_tab2d(zeta(:,:,kstp),'zeta step3d_fast #1','r'
     &  ,ondevice=.TRUE.)
C$OMP END MASTER
# endif 
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_zeta_update.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
