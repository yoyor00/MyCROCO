# ifdef NBQ_HZCORRECT_ZETA
!-----------------------------------------------------------------------
!  Update zeta(m+1)
!-----------------------------------------------------------------------
!
! ! KERNEL_33  Dnew <= ( zeta, h, pm, pn, Dnew, ubar, on_u, vbar, om_v )
! ! KERNEL_33  zeta <= ( zeta, rmask )
!$acc kernels default( present )
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
#  ifdef NBQ_MASS
     &                                          *rhobar_nbq(i,j,knew)
#  endif
          enddo
        enddo

        do j=Jstr,Jend
          do i=Istr,Iend
            zeta(i,j,knew)=( 
#  ifndef MVB             
     &                  (h(i,j)+zeta(i,j,kstp))
#  else
     &                  (h0_mvb(i,j)+dh_mvb(i,j,kstp2)+zeta(i,j,kstp))
#  endif
#  ifdef NBQ_MASS
     &                                          *rhobar_nbq(i,j,kstp)
#  endif
     &                                + (dtfast*pm(i,j)*pn(i,j)*0.5*(
     &      (Dnew(i  ,j)+Dnew(i-1,j))*(ubar(i  ,j,knew)
#  ifdef MRL_WCI
     &                               +ust2d(i  ,j)
#  endif  
     &                                                 )*on_u(i  ,j)
     &     -(Dnew(i+1,j)+Dnew(i  ,j))*(ubar(i+1,j,knew)
#  ifdef MRL_WCI
     &                               +ust2d(i+1,j)
#  endif 
     &                                                 )*on_u(i+1,j)
     &     +(Dnew(i,j  )+Dnew(i,j-1))*(vbar(i,j  ,knew)
#  ifdef MRL_WCI
     &                               +vst2d(i,j  )
#  endif
     &                                                 )*om_v(i,j  )
     &     -(Dnew(i,j+1)+Dnew(i,j  ))*(vbar(i,j+1,knew)
#  ifdef MRL_WCI
     &                               +vst2d(i,j+1)
#  endif
     &                                                 )*om_v(i,j+1))) )
#  ifdef NBQ_MASS
     &                                             /rhobar_nbq(i,j,knew)
#  endif
     &                                                          - h(i,j)
          enddo
        enddo
!$acc end kernels
!$acc kernels default( present )
! !
! !--------------------------------------------------------------------
! !  Set masking for zeta, including wet/dry conditions
! !--------------------------------------------------------------------
! !
#  ifdef MASKING
! ! KERNEL_34  zeta <= ( zeta, rmask )
      do j=Jstr,Jend
        do i=Istr,Iend
          zeta(i,j,knew)=zeta(i,j,knew)*rmask(i,j)
#   ifdef WET_DRY
!    modify new free-surface to ensure that depth 
!    is > Dcrit in masked cells.
          cff=0.5+SIGN(0.5,Dcrit(i,j)-h(i,j))
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   cff*(Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
#   endif
        enddo
      enddo 
#  endif /* MASKING */
!$acc end kernels
! !
! !--------------------------------------------------------------------
! !  Set boundary conditions for the free-surface
! !  --> ensure closed boundaries
! !--------------------------------------------------------------------
! !
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifndef OBC_NBQ
       call zetabc_tile (Istr,Iend,Jstr,Jend)
#  endif
# endif /* NBQ_HZCORRECT_ZETA */
! !
! !--------------------------------------------------------------------
! !  Update Zt_avg1 at last fast step
! !--------------------------------------------------------------------
! !
      if (LAST_FAST_STEP) then
!$acc kernels default(present)
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,knew)
          enddo
        enddo
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
!     &                          Zt_avg1(START_2D_ARRAY))
!#  endif
!$acc end kernels
! !$acc update host( Zt_avg1,zeta(:,:,knew)  )  !! iif=last
!$acc update host( Zt_avg1,zeta  )  !! iif=last
      endif
! !
! !--------------------------------------------------------------------
! ! Update grid parameters at m+1: Hz, z_r, z_w
! ! in prognostic or diagnostic way
! !--------------------------------------------------------------------
! !
#  ifdef NBQ_GRID_SLOW
      if (LAST_FAST_STEP) then
!$acc update host( Hz )      
#  endif

#  ifdef NBQ_HZ_PROGNOSTIC
!
!  Prognostic evaluation using momentum divergence
!
        do k=1,N
          do j=JstrV-1,Jend
            do i=IstrU-1,Iend
              Hz(i,j,k)=Hz_bak2(i,j,k) - dtfast*thetadiv2_nbq(i,j,k)
              Hzr(i,j,k)=(Hz(i,j,k)-rho_nbq(i,j,k))/(1.+rho(i,j,k)/rho0)
              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
              z_r(i,j,k)=0.5*(z_w(i,j,k)+z_w(i,j,k-1))
            enddo
          enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ATTENTION FRANCIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#     ifndef M3FAST_SEDLAYERS
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hzr(START_2D_ARRAY,1))     
        call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                          z_w(START_2D_ARRAY,0))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          z_r(START_2D_ARRAY,1))
#      else
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,-N_sl+1))
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          Hzr(START_2D_ARRAY,-N_sl+1))     
        call exchange_w3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          z_w(START_2D_ARRAY,-N_sl))
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          z_r(START_2D_ARRAY,-N_sl+1))
#      endif
#   endif
#  else   /* ! NBQ_HZ_PROGNOSTIC */
! !
! !--------------------------------------------------------------------
! !  Diagnostic evaluation from zeta(m+1)
! !--------------------------------------------------------------------
! !
# ifdef OPENACC
#  undef exchange_r2d_tile 
#  undef exchange_u2d_tile 
#  undef exchange_v2d_tile 
#  undef exchange_r3d_tile 
#  undef exchange_u3d_tile 
#  undef exchange_v3d_tile 
#  undef exchange_w3d_tile 
# endif
        call set_depth_tile(Istr,Iend,Jstr,Jend)
#  endif  /* NBQ_HZ_PROGNOSTIC */ 
! !
! !--------------------------------------------------------------------
! ! Compute derived grid parameters if fast update
! !--------------------------------------------------------------------
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ATTENTION FRANCIS: this call should not be needed
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
#  ifndef NBQ_GRID_SLOW
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u  , Hzw_nbq_inv_v)
#  endif
#  ifdef NBQ_GRID_SLOW
      endif !<-- LAST_FAST_STEP
#  endif
! !
! !--------------------------------------------------------------------
! !  Exchange boundary information.
! !   FRANCIS EXCHANGE OUT to be tested
! !--------------------------------------------------------------------
! !
# ifdef NBQ_HZCORRECT 
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      if (LAST_FAST_STEP) then  
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))    
      else
#   ifdef OPENACC      
      call exchange_r2d_tile_device (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
#   else
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
#   endif
      endif
#  endif
# endif
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! FRANCIS EXCHANGE OUT to be tested
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (LAST_FAST_STEP) then  
!     call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
!    &                        DU_avg1(START_2D_ARRAY,nnew))
!     call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
!    &                        DV_avg1(START_2D_ARRAY,nnew))
!     call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
!    &                        DU_avg2(START_2D_ARRAY))
!     call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
!    &                        DV_avg2(START_2D_ARRAY))
#  if defined MRL_WCI && defined WET_DRY
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ust2d(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vst2d(START_2D_ARRAY))
#  endif
      endif
