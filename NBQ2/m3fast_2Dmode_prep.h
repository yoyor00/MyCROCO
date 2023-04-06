! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_2Dmode_prep.h (begin)
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
!$acc wait( sync_rho_rufrc_z_w )
! !$acc wait( sync_ruv_int_nbq )
!$acc update device(ru_int_nbq,rv_int_nbq,rw_int_nbq) !!iff=1
       endif
#endif
# if defined M3FAST_UV || defined M3FAST_W || defined M3FAST_RHO
!$acc kernels default(present)
      if (FIRST_FAST_STEP) then
       if (FIRST_TIME_STEP) then
         do k=1,N
           do j=JstrV-2,Jend+1
             do i=IstrU-2,Iend+1
               rho_grd(i,j,k)=rho(i,j,k) 
             enddo
            enddo
          enddo
#  ifdef M3FAST_SEDLAYERS
         do k=-N_sl+1,0
           do j=JstrV-2,Jend+1
             do i=IstrU-2,Iend+1
               rho_grd(i,j,k) = rho_sdl
             enddo
            enddo
          enddo
#  endif
       endif
! !      
! !********************************
! ! Extrapolation in time:
! !********************************
! !
         do j=JstrR,JendR
           do k=1,N
             do i=IstrR,IendR
               rho_grd(i,j,k)=(rho(i,j,k)*1.5-rho_grd(i,j,k)*0.5)/rho0
             enddo
           enddo
         enddo
      endif
!$acc end kernels
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
      mybeta=0.281105 ! parameter for AB3 extrapolation

       if (FIRST_FAST_STEP) then     
                                      ! Meaning of temporal indices
        kbak=kstp                     ! ------- -- -------- -------
        kold=kstp                     ! m-2     m-1      m      m+1
        cff1= 1.0                     ! kold    kbak     kstp   knew
        cff2= 0.0
        cff3= 0.0
       elseif (FIRST_FAST_STEP+1) then  
        kbak=kstp-1                   ! AB2 forward scheme
        if (kbak.lt.1) kbak=4
        kold=kbak
# ifdef M3FAST_ZETAW
        cff1= 1.5
        cff2=-0.5
# else
        cff1= 1.   ! Just to agree with step2d
        cff2= 0. 
# endif  
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
! !
! !--------------------------------
! ! Extrapolate (D,ubar,vbar) at m+1/2
! !--------------------------------
! !
! !--------------------------------
! ! Total depth/mass at m+1/2
! !--------------------------------
! !
      if (FIRST_FAST_STEP) then
!$acc update device( h, zeta, ubar, vbar ) !iif=1
      endif
! !
!$acc kernels default(present)
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
! ! m3fast_2Dmode_prep.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
