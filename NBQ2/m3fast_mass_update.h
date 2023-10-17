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
!$acc kernels if(compute_on_device) default(present)
      do k=-N_sl+1,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            rho_nbq(i,j,k) = rho_nbq(i,j,k)  
     &                       - dtfast*thetadiv_nbq(i,j,k)
          enddo
        enddo
      enddo
! !
! !********************************
! ! NHINT numerical mode control: 
! ! remove potential surface component
! !********************************
! !
#  ifdef NHINT_CORR
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
       
              do k=N-1,Max(1,N-2*alphaNw_nbq),-1
               cff=(alphaw_nbq-1.)
!    &               *exp(-(z_w(i,j,k)            -z_w(i,j,N))**2
!    &                    /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2)
               cff2= cff*(
     &    +         qdmw_nbq(i,j,N)
     &         *(z_w(i,j,k)+H(i,j))/(z_w(i,j,N)+H(i,j))
     &         *(Hz(i,j,k)+Hz(i,j,k+1))/Hz(i,j,N))

               qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)+cff2
             enddo
             k=N
             cff2=qdmw_nbq(i,j,N)*(alphaw_nbq-1.)
             qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)+cff2
           enddo
         enddo
#   endif /* NHINT_CORR  */
!$acc end kernels
! !
! !********************************
! !  BC on rho_nbq (AGRIF only)
! !********************************
! !
#   ifdef AGRIF
!       call rnbq_bc_tile(Istr,Iend,Jstr,Jend, work) 
#   endif
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

#  ifdef M3FAST_DIAGACOUS
! !
! !********************************
! ! Diag. Acoustic: pressure field
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)
        do j=Jstr,Jend
          do i=Istr,Iend
            do k=-N_sl+1,N
              p_nbq(i,j,k)=soundspeed2_nbq(i,j,k)*rho_nbq(i,j,k)
              p_nbq_max(i,j,k)=max(p_nbq_max(i,j,k),
     &                      abs(soundspeed2_nbq(i,j,k)*rho_nbq(i,j,k)))
            enddo
          enddo
        enddo
!$acc end kernels        
#  endif /* M3FAST_DIAGACOUS */
# endif  /* M3FAST_RHO */
! !
! !********************************
! ! Mass/volume conservation (diag)
! !********************************
! !
# ifdef  NBQ_DIAGMASS
      masstot =0.
      masstot2=0.
      do j=Jstr,Jend
        do i=Istr,Iend
          do k=1,N
#  ifdef NBQ_MASS
          masstot =masstot +Hz (i,j,k)
          masstot2=masstot2+Hzr(i,j,k)
#  else
          masstot =masstot 
     &                             +rho_nbq(i,j,k)
     &                             +(1.+rho_grd(i,j,k))*Hzr(i,j,k)
          masstot2=masstot2+Hz (i,j,k)
#  endif
          enddo
        enddo
      enddo
#  ifdef MPI
      call MPI_ALLGATHER(masstot,1,MPI_DOUBLE_PRECISION,
     &                allmasstot,1,MPI_DOUBLE_PRECISION,
     &                          MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(masstot2,1,MPI_DOUBLE_PRECISION,
     &                allmasstot2,1,MPI_DOUBLE_PRECISION,
     &                          MPI_COMM_WORLD,ierr)
      masstot=QuadZero
      masstot2=QuadZero
      do i=1,NNODES
        masstot =masstot +allmasstot (1,i)
        masstot2=masstot2+allmasstot2(1,i)
      enddo
#  endif
      if (mynode==0) then
       if (FIRST_TIME_STEP.and.FIRST_FAST_STEP) then
        masstot0=masstot
        masstot02=masstot2
       endif
       open(unit=10,file='diag_mass.dat',access='append')
       write(10,*) 1.d0-masstot/masstot0,1.d0-masstot2/masstot02
       close(10)
      endif
# endif
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_mass_update.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
