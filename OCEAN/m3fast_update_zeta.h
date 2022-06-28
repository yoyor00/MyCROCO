#  ifdef M3FAST_ZETAW

! ! KERNEL_3  usurf_nbq <= ( qdmu_nbq, Hz, umask )
! ! KERNEL_3  vsurf_nbq <= ( qdmv_nbq, Hz, vmask )
! ! KERNEL_3  wsurf_nbq <= ( qdmw_nbq, Hzw_nbq_inv, rmask )

!$acc kernels default(present)
#   ifndef M3FAST_SPDUP
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
#  endif  /* M3FAST_SPDUP  */

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

!------------------------------------------------
!   ATTENTION replaced by zeta 3pts  FRANCIS
!------------------------------------------------
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
! !-------------------------------------------------------------------
! ! Apply point sources for river runoff simulations
! !-------------------------------------------------------------------
!
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
! !-------------------------------------------------------------------
! ! Advance zeta at m+1 with kinematic condition
! !-------------------------------------------------------------------
! !
! ! KERNEL_4  zab3 <= ( zeta )
! ! KERNEL_4  zeta <= ( zeta, wsurf_nbq, usurf_nbq, vsurf_nbq, zab3, 
! ! KERNEL_4  - - -     pm_u, pn_v,  rmask ) 

!$acc kernels default(present)
#   ifndef M3FAST_SPDUP
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
#   endif /* !M3FAST_SPDUP */

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
	  zeta(i,j,knew)=(zeta(i,j,kstp) + dtfast*( wsurf_nbq(i,j)
#   ifndef M3FAST_SPDUP
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
     &                          )))
#    endif    /* M3FAST_SPDUP */
#    ifdef MASKING
     &                                                  *rmask(i,j)
#    endif

#   ifdef M3FAST_ZETADISS
     &              + diss_zta*( zeta(i+1,j,kstp)-2.*zeta(i,j,kstp)
     &                       +zeta(i-1,j,kstp)) /2.
#   endif
        enddo
      enddo
!$acc end kernels
# else /* ! M3FAST_ZETAW */
! !
! !******************************************************************
! !
! !  Advance zeta at m+1 from continuity for the hydrostatic case ***
! !
! !******************************************************************
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
!$acc kernels default( present )
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=(zeta(i,j,kstp)
#  ifdef MVB
     &                    -(dh_mvb(i,j,knew2)-dh_mvb(i,j,kstp2))
#  endif
     &                                   + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1)))
#  ifdef MASKING
     &                                                   *rmask(i,j)
#  endif
#  ifdef M3FAST_ZETADISS
     &              + diss_zta*( zeta(i+1,j,kstp)-2.*zeta(i,j,kstp)
     &                       +zeta(i-1,j,kstp)) /2.*pm(i,j)**2
#  endif

        enddo
      enddo
!$acc end kernels      
# endif /* M3FAST_ZETAW */
!
!-----------------------------------------------------------------------
!  Add nudging terms
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!  Compute wet/dry masks
!-----------------------------------------------------------------------
!
! First: modify new free-surface to ensure that depth 
!        is > Dcrit in masked cells.
!
# if defined WET_DRY && defined MASKING
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=0.5+SIGN(0.5,Dcrit(i,j)-h(i,j))
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   cff*(Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
        enddo
      enddo 
# endif
!
! Then compute wet/dry masks
!
# ifdef WET_DRY
      call wetdry_tile (Istr,Iend,Jstr,Jend)
# endif
! !
! !--------------------------------------------------------------------
! !  Set boundary conditions for the free-surface
! !--------------------------------------------------------------------
! !
      call zetabc_tile (Istr,Iend,Jstr,Jend)
# if !defined NBQ_HZCORRECT
      if (LAST_FAST_STEP) then
!$acc update host( zeta )        
      endif
# endif      
! !
! !--------------------------------------------------------------------
! !  Perform exchanges
! !   CAUTION: over 1 point only!
! !--------------------------------------------------------------------
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!#  ifdef M3FAST_ZETAW
!       call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,
!      &                        zeta(START_2D_ARRAY,knew))
!#  else     
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
!#  endif     
# endif
# if !defined NBQ_HZCORRECT
      if (LAST_FAST_STEP) then
!$acc update host( zeta )        
      endif
# endif      


