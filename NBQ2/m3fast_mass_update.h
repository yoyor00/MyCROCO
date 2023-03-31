! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_mass_update.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! ! Solve rho_nbq from mass conservation 
! ! using a Forward-Backward scheme:
! ! rho_nbq(m+1) = rho_nbq(m) - dtfast*DIV(m+1)
! !********************************
! !
# ifdef M3FAST_RHO
!$acc kernels default( present )
      do k=-N_sl+1,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            rho_nbq(i,j,k) = rho_nbq(i,j,k)  
     &                       - dtfast*thetadiv_nbq(i,j,k)
          enddo
        enddo
      enddo
#  ifdef NHINT_CORR
! !
! !********************************
! ! NHINT numerical mode control: 
! ! remove potential surface component
! !********************************
! !
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
             
              k=N
              cff2=qdmw_nbq(i,j,N)*(alphaw_nbq-1.)
              qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)+cff2
              
              do k=N-1,N-2*alphaNw_nbq,-1
               cff3=cff2
               cff=(1.-alphaw_nbq)
     &               *exp(-(z_w(i,j,k)            -z_w(i,j,N))**2
     &                    /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2)
               cff2= cff*(
     &    -      2.*qdmw_nbq(i,j,k)
     &    +         qdmw_nbq(i,j,N)
     &         *(z_w(i,j,k)+H(i,j))/(z_w(i,j,N)+H(i,j))
     &         *(Hz(i,j,k)+Hz(i,j,k+1))/Hz(i,j,N))

               qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)+cff2

!               rho_nbq(i,j,k+1)=rho_nbq(i,j,k+1)
!    &           -dtfast*( cff3*Hzw_nbq_inv(i,j,k+1)
!    &                    -cff2*Hzw_nbq_inv(i,j,k  ))
             enddo
           enddo
         enddo
#    if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq(START_2D_ARRAY,1))   
#    endif
#   endif /* NHINT_CORR  */
!$acc end kernels
! !
! !********************************
! !  BC on rho_nbq
! !********************************
! !
!       call rnbq_bc_tile(Istr,Iend,Jstr,Jend, work) !! leave commented
! !
! !********************************
! !  Acoustic wave emission
! !********************************
! !
#  ifdef M3FAST_SACOUS
#   include "m3fast_sacous.h"
#  endif
! !
! !********************************
! !  rhobar_nbq: depth-mean density 
! !              (/rho0)
! !********************************
! !
#  ifdef NBQ_MASS
!
! Compute rhobar_nbq(m+1) used in zeta diagnostic from 
! depth-integrated continuity equation
!
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          rhobar_nbq(i,j,knew)=0.
        enddo  
      enddo
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            rhobar_nbq(i,j,knew)= rhobar_nbq(i,j,knew)
     &                             +rho_nbq(i,j,k)
     &                             +rho_grd(i,j,k)*Hzr(i,j,k)
          enddo  
        enddo  
      enddo
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          rhobar_nbq(i,j,knew) = 1.+(rhobar_nbq(i,j,knew)) 
     &                            / (z_w(i,j,N)-z_w(i,j,0))
        enddo
      enddo
! !
! !********************************
! !   ATTENTION exchange !!!!!!! FRANCIS
! !********************************
! !
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend, 
     &                           rhobar_nbq(START_2D_ARRAY,knew))
#   endif

#   ifdef RVTK_DEBUG
!       call check_tab2d(rhobar_nbq(:,:,knew),'rhobar_nbq','r')
#   endif    
#  endif /* NBQ_MASS */

# endif  /* M3FAST_RHO */
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_mass_update.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
