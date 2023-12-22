# ifdef NBQ_HZCORRECT_ZETA
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_zeta_correct.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! !  Update zeta(m+1)
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)
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
!$acc kernels if(compute_on_device) default(present)
! !
! !********************************
! !  Set masking for zeta, 
! !  including wet/dry conditions
! !********************************
! !
#  ifdef MASKING
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
! !********************************
! !  Set boundary conditions 
! !   for the free-surface
! !  --> ensure closed boundaries
! !********************************
! !
#  ifndef OBC_NBQ
       call zetabc_tile (Istr,Iend,Jstr,Jend)
#  endif
# endif /* NBQ_HZCORRECT_ZETA */
! !
! !********************************
! ! Update Zt_avg1 at last fast step
! !********************************
! !
      if (LAST_FAST_STEP) then
!$acc kernels if(compute_on_device) default(present)
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,knew)
          enddo
        enddo
!$acc end kernels
! !$acc update host( Zt_avg1,zeta(:,:,knew)  )  !! iif=last
!$acc update host( Zt_avg1,zeta  )  !! iif=last
      endif
! !
! !********************************
! ! Update grid parameters at m+1: 
! !     Hz, z_r, z_w
! ! in prognostic or diagnostic way
! !********************************
! !
# ifdef NBQ_GRID_SLOW
      if (LAST_FAST_STEP) then
!$acc update host( Hz )      
# endif

# ifdef NBQ_HZ_PROGNOSTIC
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
! ! 
! !********************************
! !  Exchange:  ATTENTION FRANCIS
! !********************************
! !
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#    ifndef K3FAST_SEDLAYERS
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hzr(START_2D_ARRAY,1))     
        call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                          z_w(START_2D_ARRAY,0))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          z_r(START_2D_ARRAY,1))
#    else
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,-N_sl+1))
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          Hzr(START_2D_ARRAY,-N_sl+1))     
        call exchange_w3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          z_w(START_2D_ARRAY,-N_sl))
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          z_r(START_2D_ARRAY,-N_sl+1))
#    endif
#   endif
#  else   /* ! NBQ_HZ_PROGNOSTIC */
! !
! !********************************
! ! Diagnostic evaluation from zeta(m+1)
! !********************************
! !
#  ifdef OPENACC
#   undef exchange_r2d_tile 
#   undef exchange_u2d_tile 
#   undef exchange_v2d_tile 
#   undef exchange_r3d_tile 
#   undef exchange_u3d_tile 
#   undef exchange_v3d_tile 
#   undef exchange_w3d_tile 
#  endif
        call set_depth_tile(Istr,Iend,Jstr,Jend)
# endif  /* NBQ_HZ_PROGNOSTIC */ 
! !
! !********************************
! ! Compute derived grid parameters 
! !     if fast update
! !********************************
! ! ATTENTION FRANCIS: this call 
! !   should not be needed
! !********************************
! !
# ifndef NBQ_GRID_SLOW
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u  , Hzw_nbq_inv_v)
# endif
# ifdef NBQ_GRID_SLOW
      endif !<-- LAST_FAST_STEP
# endif
! !
! !********************************
! ! Exchange ZETA
! !********************************
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      if (LAST_FAST_STEP) then  
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))    
      else
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
      endif
# endif
! !
! !********************************
! ! Copy density for extrapolation
! !********************************
! !
# ifdef NBQ_MASS
      if (LAST_FAST_STEP) then
        do k=1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              rho_grd(i,j,k)=rho(i,j,k)
            enddo
           enddo
         enddo
      endif
# endif
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_zeta_correct.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
