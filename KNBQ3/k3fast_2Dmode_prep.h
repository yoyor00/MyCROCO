! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_2Dmode_prep.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! ! Store rho (from slow mode) 
! !  at first slow 
! !  and fast time-step
! !********************************
! !
#ifdef OPENACC
       if (FIRST_FAST_STEP) then
!!!!$acc wait( sync_rho_rufrc_z_w )
!!! !$acc wait( sync_ruv_int_nbq )
!!!!$acc update device(ru_int_nbq,rv_int_nbq,rw_int_nbq) !!iff=1
       endif
#endif
# if defined K3FAST_UV || defined K3FAST_W || defined K3FAST_RHO
      if (FIRST_FAST_STEP) then
       if (FIRST_TIME_STEP) then
!$acc kernels if(compute_on_device) default(present) async(1)
         do j=JstrR,JendR
         do k=1,N
           do i=IstrR,IendR
               rho_grd(i,j,k)=rho(i,j,k) /rho0
#  ifdef K3FAST_NOBPG
               rho_bpg(i,j,k)=rho(i,j,k)/rho0
#   ifdef NBQ_GRAV
     &                       -rho_nh(i,j,k)/rho0
#   endif
#  endif
             enddo
            enddo
          enddo
#  ifdef K3FAST_SEDLAYERS
         do k=-N_sl+1,0
           do j=JstrR,JendR
           do i=IstrR,IendR
               rho_grd(i,j,k) = rho_sdl
             enddo
            enddo
          enddo
#  endif
!$acc end kernels
       endif
! !      
! !********************************
! ! Extrapolation in time:
! !********************************
! !
!$acc kernels if(compute_on_device) default(present) async(1)
         do j=JstrR,JendR
           do k=1,N
             do i=IstrR,IendR
               rho_grd(i,j,k)=rho(i,j,k)*1.5/rho0-rho_grd(i,j,k)*0.5
#  ifdef K3FAST_NOBPG
      rho_bpg(i,j,k)=1.5*rho(i,j,k)/rho0-rho_bpg(i,j,k)*0.5
#   ifdef NBQ_GRAV
     &               -rho_nh(i,j,k)*1.5/rho0
#   endif
#  endif
             enddo
           enddo
         enddo
!$acc end kernels
#  ifdef MPI
             call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_grd(-2,-2,1))
#   ifdef K3FAST_NOBPG
             call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_bpg(-2,-2,1))
#   endif
#endif

      endif
# endif
! !
! !********************************
! ! AB3 Forward Step: compute total depth of water column and 
! !                   vertically integrated mass fluxes which
! ! --- ------- ----   are needed to compute 
! ! rhs terms of the barotropic momentum equations (rubar,rvbar).
! !********************************
! !
! !--------------------------------
! !  Set indices to extrapolate 
! !  (D,ubar,vbar) at m+1/2 (AB3)
! !--------------------------------
! !
# ifdef K3FAST_2DCONT
       if (FIRST_FAST_STEP.and.FIRST_TIME_STEP) then   
# else
       if (FIRST_FAST_STEP) then     
# endif
                                      ! Meaning of temporal indices
        kbak=kstp                     ! ------- -- -------- -------
        kold=kstp                     ! m-2     m-1      m      m+1
        cff1= 1.0                     ! kold    kbak     kstp   knew
        cff2= 0.0
        cff3= 0.0
# ifdef K3FAST_2DCONT
       elseif (FIRST_FAST_STEP+1.and.FIRST_TIME_STEP) then  
# else
       elseif (FIRST_FAST_STEP+1) then  
# endif
        kbak=kstp-1                   ! AB2 forward scheme
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1= 1.5
        cff2=-0.5
        cff3= 0.0
      else                             ! AB3 forward scheme
        kbak=kstp-1 
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1= 1.5+mybeta
        cff2=-2.0*mybeta-0.5
        cff3= mybeta
      endif
      if (FIRST_FAST_STEP.and.FIRST_TIME_STEP) then
# ifdef K3FAST_AB3
       kab3_1=2
       kab3_2=1
# endif
# if defined K3FAST_AM4 || defined K3FAST_AM4b || defined K3FAST_AM4c
       kam4_1=3
       kam4_2=2
       kam4_3=1
# endif
      else
# ifdef K3FAST_AB3
       kab3_1=kab3_1+1
       if (kab3_1.ge.3) kab3_1=1
       kab3_2=kab3_2+1
       if (kab3_2.ge.3) kab3_2=1
# endif
# if defined K3FAST_AM4 || defined K3FAST_AM4b || defined K3FAST_AM4c
       kam4_1=kam4_1+1
       if (kam4_1.ge.4) kam4_1=1
       kam4_2=kam4_2+1
       if (kam4_2.ge.4) kam4_2=1
       kam4_3=kam4_3+1
       if (kam4_3.ge.4) kam4_3=1
# endif
      endif

! !
! !--------------------------------
! ! Extrapolate (D,ubar,vbar) at m+1/2
! !--------------------------------
! !
! !--------------------------------
! ! Total depth/mass at m+1/2
! !--------------------------------
! !
!      if (FIRST_FAST_STEP) then
!!$acc update device( h, zeta, ubar, vbar ) !iif=1
!      endif
! !
!$acc kernels if(compute_on_device) default(present) async(1)
      do j=JstrV-2,Jend+1
       do i=IstrU-2,Iend+1
# ifdef NBQ_MASS
#  ifndef MVB
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+h(i,j))*rhobar_nbq(i,j,kstp)
     &             +cff2*(zeta(i,j,kbak)+h(i,j))*rhobar_nbq(i,j,kbak)
     &             +cff3*(zeta(i,j,kold)+h(i,j))*rhobar_nbq(i,j,kold)
#  else
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+dh_mvb(i,j,knew2))
     &                                  *rhobar_nbq(i,j,kstp)
     &             +cff2*(zeta(i,j,kbak)+dh_mvb(i,j,kstp2))
     &                                  *rhobar_nbq(i,j,kbak)
     &             +cff3*(zeta(i,j,kold)+dh_mvb(i,j,kbak2))
     &                                  *rhobar_nbq(i,j,kold)
#  endif
# else       
#  ifndef MVB
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+h(i,j))
     &             +cff2*(zeta(i,j,kbak)+h(i,j))
     &             +cff3*(zeta(i,j,kold)+h(i,j))
#  else     
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+dh_mvb(i,j,knew2))   
     &             +cff2*(zeta(i,j,kbak)+dh_mvb(i,j,kstp2))
     &             +cff3*(zeta(i,j,kold)+dh_mvb(i,j,kbak2))
#  endif     
# endif /* NBQ_MASS */
       enddo
      enddo
! !
! !--------------------------------
! ! Depth-average ubar velocity at m+1/2
! !--------------------------------
! !
      do j=Jstr-1,Jend+1
       do i=IstrU-1,Iend+1
          urhs(i,j)=cff1*ubar(i,j,kstp) 
     &             +cff2*ubar(i,j,kbak)
     &             +cff3*ubar(i,j,kold)
# if defined MRL_WCI && defined MASKING
          urhs(i,j)=urhs(i,j)*umask(i,j)+ust2d(i,j)*(umask(i,j)-1.0)
# endif
          DUon(i,j)=0.5*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
# ifdef MRL_WCI
     &                                                   + ust2d(i,j)
# endif
     &                                                              )
        enddo
      enddo
! !
! !--------------------------------
! ! Depth-average vbar velocity at m+1/2
! !--------------------------------
! !
      do j=JstrV-1,Jend+1
        do i=Istr-1,Iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp)
     &             +cff2*vbar(i,j,kbak)
     &             +cff3*vbar(i,j,kold)
# if defined MRL_WCI && defined MASKING
          vrhs(i,j)=vrhs(i,j)*vmask(i,j)+vst2d(i,j)*(vmask(i,j)-1.0)
# endif
          DVom(i,j)=0.5*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
# ifdef MRL_WCI
     &                                                   + vst2d(i,j)
# endif
     &                                                              )
        enddo
      enddo
!$acc end kernels
# ifdef OBC_VOLCONS
      call set_DUV_bc_tile (Istr,Iend,Jstr,Jend, Drhs, DUon,DVom)
# endif
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_2Dmode_prep.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
