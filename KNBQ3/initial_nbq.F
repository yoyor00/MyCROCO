! $Id:$
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#ifdef K3FAST
!
      subroutine initial_nbq (tile)
!
      implicit none
      integer tile, icall, trd
# include "param.h"
# include "private_scratch.h"
!$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      trd=0
!$    trd=omp_get_thread_num()

      call initial_nbq_tile (Istr,Iend,Jstr,Jend
# ifdef K3FAST
#  ifdef NONLIN_EOS
     &                     ,A2d(1,1,trd), A2d(1,2,trd)
#  endif
# endif
     &                      )
      end subroutine initial_nbq
!
      subroutine initial_nbq_tile (Istr,Iend,Jstr,Jend
# ifdef K3FAST
#  ifdef NONLIN_EOS
     &                           ,K_up, K_dw
#  endif
# endif
     &                         )
!
!======================================================================
!
!                       NBQ initialization
!
!======================================================================
!
      implicit none
      integer  Istr, Iend, Jstr, Jend
# ifdef MPI      
      include 'mpif.h'
# endif      
# include "param.h"
# include "scalars.h"
# include "private_scratch.h"
# if defined K3FAST_UV || defined K3FAST_W || defined K3FAST_RHO
#  include "nbq.h"
# endif
# include "work.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# ifdef K3FAST
#  ifdef NONLIN_EOS      
      real K_up(PRIVATE_1D_SCRATCH_ARRAY,0:N)  ! work arrays for call
      real K_dw(PRIVATE_1D_SCRATCH_ARRAY,0:N)  ! to nonlinear EOS
#  endif
# endif
# ifdef PEKERIS_ACOUS3
      real zbar
# endif
      integer i,j,k

# include "compute_extended_bounds.h"

# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
!
# if defined K3FAST_UV || defined K3FAST_W

#  ifdef NBQ_IMP
        ifl_imp_nbq = 1
        MPI_master_only write(6,*)
        MPI_master_only write(6,*) '--------------------------------'
        MPI_master_only write(6,*) ' NBQ: semi-implicit integration '
        MPI_master_only write(6,*) '--------------------------------'
        MPI_master_only write(6,*)
#  else
        ifl_imp_nbq = 0
        MPI_master_only write(6,*)
        MPI_master_only write(6,*) '---------------------------'
        MPI_master_only write(6,*) ' NBQ: explicit integration '
        MPI_master_only write(6,*) '---------------------------'
        MPI_master_only write(6,*)
#  endif
! !
! !--------------------------------------------------------------------
! !  Initialize parameters
! !--------------------------------------------------------------------
! !
        ifl_nbq  = 1
        slip_nbq = 0
! !
! !  Pseudoacoustic speed
! !     
        if (csound_nbq.lt.0.) then
# ifdef KNBQ
           csound_nbq=sqrt(g*hmax)*2.
# endif
# ifdef KNHINT
           csound_nbq=0.5/CFL_nbq
#  ifdef K3FAST_PG2
!           csound_nbq=csound_nbq*1.7
#  endif
# endif
        endif
        if (visc2read_nbq.lt.0.) then
!           visc2read_nbq=0.95*min(grdmin
!    &      ,hmin/real(N))**2/dtnbq
!          visc2read_nbq=0.95*dtnbq/CFL_nbq*csound_nbq
           visc2read_nbq=0.1*csound_nbq*min(grdmin
     &      ,hmin/real(N))
#  ifdef K3FAST_PG2
!          visc2read_nbq=visc2read_nbq*1.7
#  endif
        endif
        if (mynode==0) then
         write(stdout,1) "[K3FAST] Cs, Visc2 recomputed: "
     & ,csound_nbq," m/s, ",visc2read_nbq," m2/s"
   1          format(A,F6.2,A,F6.2,A)
         write(stdout,*)
        endif

        visc2_nbq=visc2read_nbq
        visc2v_nbq=visc2read_nbq   
        soundspeed_nbq  =csound_nbq
        soundspeed2_nbq =csound_nbq**2
        
#ifdef K3FAST_CSVISC2K        
!$acc update device( visc2_nbq,  visc2v_nbq
!$acc&          , soundspeed_nbq, soundspeed2_nbq
!$acc&           )       
#endif        
!$acc kernels if(compute_on_device) default(present)
#  if defined K3FAST_SEDLAYERS && defined K3FAST_CSVISC2K
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
               soundspeed_nbq (i,j,k)=csound_nbq
               soundspeed2_nbq(i,j,k)=csound_nbq**2
            enddo
          enddo
        enddo
#  endif
! !
! ! Depth KNHINT
! !
#  ifdef KNHINT_CORR
        alphaw_nbq  = 1. 
        alphaNw_nbq = 5  
#  endif
!
! Bulk viscosity (horizontal & vertical components)
!     
#  if defined NBQ_SPONGE || defined KNHINT
#   ifdef K3FAST_CSVISC2K
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
               visc2_nbq(i,j,k) =visc2_nbq(i,j,k)
     &          *(1.
#   ifdef NBQ_SPONGE 
     &              +visc2_nbq_sponge(i,j)
#   endif
#   ifdef KNHINT_CORR_visc2
     &           +exp(-(z_w(i,j,k   )-         z_w(i,j,N))**2
     &                /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2)
#   endif
     &           )
               visc2v_nbq(i,j,k) =visc2v_nbq(i,j,k)
     &          *(1.
#   ifdef NBQ_SPONGE 
     &              +visc2_nbq_sponge(i,j) 
#   endif
#   ifdef KNHINT_CORR_visc2
     &           +exp(-(z_w(i,j,k   )-         z_w(i,j,N))**2
     &                /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2)
#   endif
     &           )
            enddo
          enddo
        enddo
#   endif
#  endif
#  if defined K3FAST_SEDLAYERS && defined K3FAST_CSVISC2K
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
               visc2_nbq(i,j,k) =visc2read_nbq
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
#  endif
!$acc end kernels
c c c !$acc update host( visc2_nbq, visc2v_nbq )

!$acc kernels if(compute_on_device) default(present)
!
! Density of sediments
!
#  ifdef K3FAST_SEDLAYERS
!     rho_sdl = 1. / rho0 
      rho_sdl = 800. / rho0 
#  endif
!
!  Grid update time-step
!
#  ifdef NBQ_GRID_SLOW
        dtgrid_nbq = dt
#  else
        dtgrid_nbq = dtfast
#  endif
!
!----------------------------------------------------------------------
!   Compressible density initializations
!----------------------------------------------------------------------
!
#  ifdef K3FAST_RHO
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              rho_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
#  endif
#  ifdef NBQ_MASS
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              rho_nbq_avg1(i,j,k)=(rho0+rho(i,j,k))/rho0
            enddo
          enddo
        enddo 
        do j=JR_RANGE
          do i=IR_RANGE
            work2d(i,j)      =0.
            rhobar_nbq(i,j,:)=0.
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              work2d(i,j)      =work2d(i,j)+Hzr(i,j,k)
              rhobar_nbq(i,j,:)=rhobar_nbq(i,j,:)+
     &                          rho(i,j,k)*Hzr(i,j,k)/rho0
            enddo
          enddo
        enddo
        do j=JR_RANGE       ! Add rho0 for added precision
          do i=IR_RANGE
            rhobar_nbq(i,j,:)   =rhobar_nbq(i,j,:)/work2d(i,j) + 1.
            rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,1) 
          enddo
        enddo
#  endif
!$acc end kernels
# endif /* K3FAST_RHO */
!
!----------------------------------------------------------------------
!  NBQ Momentum initialization
!----------------------------------------------------------------------
!
!$acc kernels if(compute_on_device) default(present)
        do k=1,N
          do j=Jstr-1,Jend+1
            do i=Istr,Iend+1
# if defined K3FAST_UV && defined K3FAST_C3D_UVSF
              qdmu_nbq(i,j,k)=0.5*u(i,j,k,nrhs)*(Hz(i,j,k)+Hz(i-1,j,k))
# elif defined K3FAST_UV            
              qdmu_nbq(i,j,k)=0.
# endif              
          enddo
         enddo
        enddo

        do k=1,N
          do j=Jstr,Jend+1
            do i=Istr-1,Iend+1
# if defined K3FAST_UV && defined K3FAST_C3D_UVSF
              qdmv_nbq(i,j,k)=0.5*v(i,j,k,nrhs)*(Hz(i,j,k)+Hz(i,j-1,k))
# elif defined K3FAST_UV                      
              qdmv_nbq(i,j,k)=0.
# endif       
            enddo
          enddo
        enddo
        do k=1,N-1
          do j=Jstr-1,Jend+1
            do i=Istr-1,Iend+1
# if defined K3FAST_W && defined K3FAST_C3D_WSF
              qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*(Hz(i,j,k)+Hz(i,j,k+1))
# elif defined K3FAST_W          
              qdmw_nbq(i,j,k)=0.
# endif              
            enddo
          enddo
        enddo
        k=0 
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
# if defined K3FAST_W && defined K3FAST_C3D_WSF
            qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*Hz(i,j,k+1)
# elif defined K3FAST_W              
            qdmw_nbq(i,j,k)=0.
# endif              
          enddo
        enddo
        k=N 
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
# if defined K3FAST_W && defined K3FAST_C3D_WSF
            qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*Hz(i,j,k)
# elif defined K3FAST_W               
            qdmw_nbq(i,j,k)=0.
# endif              
          enddo
        enddo
# if defined K3FAST_RHO
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              thetadiv_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
# endif
!$acc end kernels
! !
! ! Acoustic configurations
! !
!$acc kernels if(compute_on_device) default(present)
# ifdef K3FAST_DIAGACOUS
        do k=-N_sl+1,N
          do j=JR_RANGE
            do i=IR_RANGE
               p_nbq(i,j,k)=0.
               p_nbq_max(i,j,k)=0.
            enddo
          enddo
        enddo
# endif
# ifdef PEKERIS_ACOUS1
        sacous_nbq=0
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1500. ! isocelere
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1500.
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1000.) then
               visc2_nbq(i,j,k) =0.+1.e2*cos(xr(i,j)/1000.*pi/2.)
               elseif (xr(i,j).ge.9000.) then
               visc2_nbq(i,j,k) =0.+1.e2*sin((xr(i,j)-9000)/1000.*
     &                           pi/2.)
               else
               visc2_nbq(i,j,k)=0.
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
#  ifdef M3FAST_SEDLAYERS
        rho_sdl = 0. / rho0 
#  endif
# elif defined PEKERIS_ACOUS2
        sacous_nbq=1
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1500. ! isocelere
!                 soundspeed_nbq(i,j,k)=1450.-z_r(i,j,k)/1000.*32. ! profil lineaire
!                 if (-z_r(i,j,k).le.200.) then
!                    soundspeed_nbq(i,j,k)=1500.+20*z_r(i,j,k)/200. ! sofar
!                 else
!                    soundspeed_nbq(i,j,k)=1480-40*(z_r(i,j,k)-200.)/(1000.-200.) ! sofar
!                 endif
!                 if (-z_r(i,j,k).le.100.) then
!                    soundspeed_nbq(i,j,k)=1480.-20*z_r(i,j,k)/100. ! chenalsurface
!                 else
!                    soundspeed_nbq(i,j,k)=1500-20*(z_r(i,j,k)-100.)/(1000.-100.) ! chenalsurface
!                 endif
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1800.
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1000.) then
               visc2_nbq(i,j,k) =3.5e1+1.e2*cos(xr(i,j)/1000.*pi/2.)
               elseif (xr(i,j).ge.9000.) then
               visc2_nbq(i,j,k) =3.5e1+1.e2*sin((xr(i,j)-9000)/1000.*
     &                           pi/2.)
               else
               visc2_nbq(i,j,k)=3.5e1!3.5e0!7e1!3.5e1
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1000.) then
               visc2_nbq(i,j,k) =0.+1.e2*cos(xr(i,j)/1000.*pi/2.)
               elseif (xr(i,j).ge.9000.) then
               visc2_nbq(i,j,k) =0.+1.e2*sin((xr(i,j)-9000)/1000.*
     &                           pi/2.)
               else
               visc2_nbq(i,j,k)=0.
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
# elif defined PEKERIS_ACOUS3
        sacous_nbq=0
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
                 zbar = 2*(-z_r(i,j,k)-1300.0)/1300.0
                 soundspeed_nbq(i,j,k)=1500.0*(1.0+0.00737*
     &             (zbar-1.0+exp(-zbar)))
!                 soundspeed_nbq(i,j,k)=1500.
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1600.
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
#  ifdef K3FAST_DIAGACOUS
        do k=-N_sl+1,N
          do j=JR_RANGE
            do i=IR_RANGE
               p_nbq(i,j,k)=0.
               p_nbq_max(i,j,k)=0.
            enddo
          enddo
        enddo
#  endif
                do k=-N_sl+1,-N_sl+250
          do j=JR_RANGE
            do i=IR_RANGE
               visc2_nbq(i,j,k)=5e2!3.5e0!7e1!3.5e1
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=-N_sl+250,0
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1000.) then
               visc2_nbq(i,j,k) =3.5e1+1.e2*cos(xr(i,j)/1000.*pi/2.)
               elseif (xr(i,j).ge.99000.) then
               visc2_nbq(i,j,k) =3.5e1+1.e2*sin((xr(i,j)-99000)/1000.*
     &                           pi/2.)
               else
               visc2_nbq(i,j,k)=0.!3.5e1!3.5e0!7e1!3.5e1
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1000.) then
               visc2_nbq(i,j,k) =0.+1.e2*cos(xr(i,j)/1000.*pi/2.)
               elseif (xr(i,j).ge.99000.) then
               visc2_nbq(i,j,k) =0.+1.e2*sin((xr(i,j)-99000)/1000.*
     &                           pi/2.)
               else
               visc2_nbq(i,j,k)=0.
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
!
! Density of sediments
!
#  ifdef M3FAST_SEDLAYERS
      rho_sdl = 0. / rho0 
!      rho_sdl = 1. / rho0 
#  endif
# elif defined WEDGE3D
        sacous_nbq=10
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1500. ! isocelere
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo
        do k=-N_sl+1,0
          do j=JR_RANGE
            do i=IR_RANGE
                 soundspeed_nbq(i,j,k)=1700.
                 soundspeed2_nbq(i,j,k)=soundspeed_nbq(i,j,k)**2
            enddo
          enddo
        enddo   
        
        do k=-N_sl+1,-N_sl+40
          do j=JR_RANGE
            do i=IR_RANGE
               visc2_nbq(i,j,k)=5e3!3.5e0!7e1!3.5e1
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=-N_sl+40,0
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1100.) then
               visc2_nbq(i,j,k) =5e3
               elseif (xr(i,j).ge.5500.) then
               visc2_nbq(i,j,k) =5e3
               elseif (yr(i,j).le.1100.) then
               visc2_nbq(i,j,k) =5e3
               elseif (yr(i,j).ge.26000) then
               visc2_nbq(i,j,k) =5e3
               else
               visc2_nbq(i,j,k)=1.25e3!1e3
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
               if (xr(i,j).le.1000.) then
               visc2_nbq(i,j,k) =5e3
               elseif (xr(i,j).ge.5500.) then
               visc2_nbq(i,j,k) =5e3
               elseif (yr(i,j).le.500.) then
               visc2_nbq(i,j,k) =5e3
               elseif (yr(i,j).ge.26000) then
               visc2_nbq(i,j,k) =5e3
               else
               visc2_nbq(i,j,k)=1e-6
               endif
               visc2v_nbq(i,j,k)=visc2_nbq(i,j,k)
            enddo
          enddo
        enddo
!
! Density of sediments
!
#  ifdef K3FAST_SEDLAYERS
      rho_sdl = 500. / rho0 
#  endif
# endif
!$acc end kernels      
!
!----------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!----------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifdef NBQ_MASS
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rho_nbq_avg1(START_2D_ARRAY,1))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,1))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,2))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,3))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,4))
#  endif
#  if defined OBC_NBQ && defined OBC_NBQORLANSKI && defined K3FAST_UV 
       call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                         qdmu_nbq(START_2D_ARRAY,1))
       call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                         qdmv_nbq(START_2D_ARRAY,1))
#  endif
# endif

      return
      end subroutine initial_nbq_tile
#else
      subroutine initial_nbq_empty
      return
      end
#endif
