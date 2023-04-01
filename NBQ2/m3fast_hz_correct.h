! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_hz_correct.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
      if (LAST_FAST_STEP) then
! !
! !--------------------------------------------------------------------
! ! Store HZ before correction
! !--------------------------------------------------------------------
! !
       do k=-N_sl+1,N
        do j=JstrV-2,Jend+1
         do i=IstrU-2,Iend+1
              Hz_correct(i,j,k)=Hz(i,j,k)
          enddo
        enddo
       enddo
!$acc update device( Hz_correct )      
! !
! !--------------------------------------------------------------------
! ! Exchanges (1)
! !--------------------------------------------------------------------
! !
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg2(START_2D_ARRAY))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg2(START_2D_ARRAY))
#  endif
! ! 
! !--------------------------------------------------------------------
! ! Adjust depth average
! !--------------------------------------------------------------------
! ! 
        do j=Jstr,Jend
          do i=Istr,Iend
            dum_s=0.
            do k=1,N
              dum_s=dum_s+Hz(i,j,k)-Hz_bak(i,j,k)
            enddo
            do k=1,N
              Hz(i,j,k)=Hz(i,j,k)
     &                -(dum_s
     &             + (DU_avg2(i+1,j)-DU_avg2(i,j)
     &               +DV_avg2(i,j+1)-DV_avg2(i,j)
     &               )*pm(i,j)*pn(i,j)
     &               *dt)
     &               /(z_w(i,j,N)-z_w(i,j,0))
     &               *(z_w(i,j,k)-z_w(i,j,k-1))
#  ifdef MASKING
     &               *rmask(i,j)
#  endif
#  ifdef NBQ_HZCORR_DEBUG
              Hz_corr(i,j,k)=Hz(i,j,k)-Hz_correct(i,j,k)
#  endif
            enddo
          enddo
        enddo
!# undef Hz_correct  
! !
! !--------------------------------------------------------------------
! ! Exchanges (2)
! !--------------------------------------------------------------------
! !
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#   ifndef M3FAST_SEDLAYERS
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,1))
#   else
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,-N_sl+1))
#   endif
#  endif
!$acc update device( Hz )      
! !
      endif ! iif.eq.nfast
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_zeta_correct.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
