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
! !   should be around 5 times external phase speed sqrt(g*h)
!
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
        alphaw_nbq  = 0.  
        alphaNw_nbq = N 
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
