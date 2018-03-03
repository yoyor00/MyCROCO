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
#ifdef NBQ
!
      subroutine initial_nh (tile, icall)
!
      use module_param
      implicit none
      integer tile, icall, trd
#include "private_scratch.h"
!$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
      trd=0
!$    trd=omp_get_thread_num()

      call initial_nh_tile (icall, Istr,Iend,Jstr,Jend       &
                         , A3d(1,1,trd), A3d(1,2,trd)        &
                         , A3d(1,3,trd), A3d(1,4,trd)        &
#ifdef NBQ_ZETAW
                         , A3d(1,5,trd), A3d(1,6,trd)        &
#endif
#ifdef NONLIN_EOS
                         , A2d(1,1,trd), A2d(1,2,trd)        &
#endif
                         )
  
      end subroutine initial_nh
  
      subroutine initial_nh_tile (icall, Istr,Iend,Jstr,Jend           &
                               ,Hzw_half_nbq_inv, Hzr_half_nbq_inv     &
                               ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v  &
#ifdef NBQ_ZETAW
                               ,Hzu_half_qdmu,Hzv_half_qdmv            &
#endif
#ifdef NONLIN_EOS
                               ,K_up, K_dw                             &
#endif
                               )
!
!======================================================================
!
!                 NH and NBQ initialization
!
!======================================================================
!
# ifdef TRACETXT
      use module_tracetxt_out
# endif
      implicit none
      integer :: Istr, Iend, Jstr, Jend

# ifdef MPI      
      include 'mpif.h'
# endif      
# include "param_F90.h"
# include "scalars_F90.h"
# include "private_scratch.h"
# include "nbq.h"
# include "work.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# ifdef NBQ_ZETAW
       real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# endif
       real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
       real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# ifndef NBQ_ZETAW	   	   
       real work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,N,5)
# endif
# ifdef NONLIN_EOS      
      real K_up(PRIVATE_1D_SCRATCH_ARRAY,0:N)  ! work arrays for call
      real K_dw(PRIVATE_1D_SCRATCH_ARRAY,0:N)  ! to nonlinear EOS
# endif  

      integer :: i,j,k,ierr
      integer :: icall

# include "compute_auxiliary_bounds.h"
!
# if defined EW_PERIODIC && !defined MPI
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
# endif
# if defined NS_PERIODIC && !defined MPI
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif

      if (icall.eq.1) then
!**********************************************************************
!
!                 Initializations (PART I)
!
!**********************************************************************

!----------------------------------------------------------------------
!........Semi-implicit scheme (0/1):
!----------------------------------------------------------------------
# ifdef NBQ_IMP
        ifl_imp_nbq = 1
        MPI_master_only write(6,*)
        MPI_master_only write(6,*) '--------------------------------'
        MPI_master_only write(6,*) ' NBQ: semi-implicit integration !'
        MPI_master_only write(6,*) '--------------------------------'
        MPI_master_only write(6,*)
# else
        ifl_imp_nbq = 0
        MPI_master_only write(6,*)
        MPI_master_only write(6,*) '---------------------------'
        MPI_master_only write(6,*) ' NBQ: explicit integration !'
        MPI_master_only write(6,*) '---------------------------'
        MPI_master_only write(6,*)
# endif
!
!----------------------------------------------------------------------
!  Initialize density perturbation and momentum arrays
!----------------------------------------------------------------------
!
        rho_nbq  = 0.0
        qdmu_nbq = 0.0
        qdmv_nbq = 0.0
        qdmw_nbq = 0.0
# ifdef NBQ_MASS
        rhobar_nbq  = 1.
# endif
!
!----------------------------------------------------------------------
!  Initialize parameters: should be done in a NH-namelist
!----------------------------------------------------------------------
!
        ifl_nbq  = 1
        slip_nbq = 0

        iteration_nbq_max=ndtnbq
        soundspeed_nbq =csound_nbq !!! pseudoacoustic speed for tank
        soundspeed2_nbq=csound_nbq**2

!       do j=lbound(soundspeed_nbq,2),ubound(soundspeed_nbq,2)
!       do i=lbound(soundspeed_nbq,1),ubound(soundspeed_nbq,1)
!         soundspeed_nbq(i,j)=5.*sqrt(g*h(i,j))
!         if (soundspeed_nbq(i,j)>=300) then
!          print *,i,j, soundspeed_nbq(i,j),csound_nbq,h(i,j)
!         endif
!         soundspeed2_nbq(i,j)=soundspeed_nbq(i,j)**2
!       enddo
!       enddo

!       cw_int_nbq=soundspeed_nbq !!! ~ 2-10 sqrt(gH)_max
!       cw_int_nbq=sqrt(9.81*4000.) !soundspeed_nbq !!! ~ 2-10 sqrt(gH)_max


!       MPI_master_only write(stdout,'3(A,I4/)')
!                       'NBQ: INITIALIZING ifl_nbq      =',ifl_nbq
!                       '     INITIALIZING slip_nbq     =',slip_nbq
!                       '     INITIALIZING ifl_imp_nbq  =',ifl_imp_nbq
!
!----------------------------------------------------------------------
!  Pre-numbering of grid points and matrices:
!----------------------------------------------------------------------
!
!....NH-Grid definition:
        call grid_def_nh
!
!----------------------------------------------------------------------
!... Initializes Grid-coef
!----------------------------------------------------------------------
!
        Hzw_half_nbq_inv  =1.e-30
        Hzr_half_nbq_inv  =1.e-30
        Hzw_half_nbq_inv_u=1.e-30
        Hzw_half_nbq_inv_v=1.e-30           
# ifdef NBQ_ZETAW  
        Hzu_half_qdmu     =0.
        Hzv_half_qdmv     =0.
        call grid_coef_nh( Istr,Iend,Jstr,Jend                      &
                          ,Hzw_half_nbq_inv,Hzr_half_nbq_inv        &
                          ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v    &
                          ,Hzu_half_qdmu, Hzv_half_qdmv             &
                         )
# else
        call grid_coef_nh( Istr,Iend,Jstr,Jend                      &
                          ,Hzw_half_nbq_inv,Hzr_half_nbq_inv        &
                          ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v    &
                          ,work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,1) &
                          ,work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,2) &
                         )
# endif
!
!----------------------------------------------------------------------
!... Set second viscosity coefficient:
!----------------------------------------------------------------------
!
!       csvisc1_nbq  = dtnbq * soundspeed2_nbq + visc2_nbq
!       csvisc2_nbq  = dtnbq * soundspeed2_nbq / csvisc1_nbq 
!
      endif     ! icall == 1

      if (icall.eq.2) then

!**********************************************************************
!
!                 rhobar_nbq initializations (PART II)
!
!**********************************************************************
!.........Grid update:
# ifdef NBQ_PRECISE
        dtgrid_nbq = dtfast
# else
        dtgrid_nbq = dt
# endif

!.........EOS to compute rho (if not already done):
# ifdef NONLIN_EOS
        call rho_eos_tile(Istr,Iend,Jstr,Jend,K_up,K_dw)
# else
        call rho_eos_tile(Istr,Iend,Jstr,Jend)
# endif

# ifdef NBQ_MASS
!.........Initialize NBQ density field:
        do k=1,N
          do j=jstrq_nh-1,jendq_nh+1
            do i=istrq_nh-1,iendq_nh+1
              rho_nbq_avg1(i,j,k)=(rho0+rho(i,j,k))/rho0
            enddo
          enddo
        enddo 

        rhobar_nbq     (:,:,:)=1.
        rhobar_nbq_avg1(:,:  )=1.

        do j=jstrq_nh-1,jendq_nh+1
          do i=istrq_nh-1,iendq_nh+1
            work2d(i,j)         = 0.
            rhobar_nbq(i,j,:)   = 0.
          enddo
        enddo
        do j=jstrq_nh-1,jendq_nh+1
          do i=istrq_nh-1,iendq_nh+1
            do k=1,N
              work2d(i,j)         = work2d(i,j)+Hzr(i,j,k)
              rhobar_nbq(i,j,:)   = rhobar_nbq(i,j,:)+          &
                                    rho(i,j,k)*Hzr(i,j,k)/rho0
            enddo
          enddo
        enddo
 
!....... add rho0 for added precision:
        do j=jstrq_nh-1,jendq_nh+1
          do i=istrq_nh-1,iendq_nh+1
            rhobar_nbq(i,j,:)   = rhobar_nbq(i,j,:)/work2d(i,j) + 1.
            rhobar_nbq_avg1(i,j)= rhobar_nbq(i,j,1) 
          enddo
        enddo
# else
!      rhobar_nbq     (:,:,:)=1.
!      rhobar_nbq_avg1(:,:  )=1.
# endif
        
!.......Some remaining initializations:
        thetadiv_nbq        = 0.d0

# ifdef NBQ_MASS
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh &
                               ,rhobar_nbq(START_2D_ARRAY,1))   
        call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh &
                               ,rhobar_nbq(START_2D_ARRAY,2))   
        call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh &
                               ,rhobar_nbq(START_2D_ARRAY,3))   
        call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh &
                               ,rhobar_nbq(START_2D_ARRAY,4))   
#  endif
# endif

      endif 

      if (icall == 3) then 
!**********************************************************************
!
!               NBQ initializations (PART III)
!
!**********************************************************************
        qdmu_nbq=0.
        qdmv_nbq=0.
        qdmw_nbq=0.
        rho_nbq=0.
# ifdef NBQ_ZETAW
        usurf_nbq=0.
        vsurf_nbq=0.
        wsurf_nbq=0.
# else
        zr_half_nbq=0.
        zw_half_nbq=0.
# endif
# ifdef KH_INST
        if (iic.le.1.and.iif.le.1) then
          do k=1,N
            do j=JR_RANGE
              do i=IU_RANGE
# ifdef NBQ_MASS
                qdmu_nbq(i,j,k)=(1.+0.5*(rho(i,j,k)+rho(i-1,j,k))/rho0) &
                             *0.5*u(i,j,k,nrhs)*(hz(i,j,k)+hz(i-1,j,k))
# else
                qdmu_nbq(i,j,k)=0.5*u(i,j,k,nrhs)*(hz(i,j,k)+hz(i-1,j,k))
# endif
              enddo
            enddo
          enddo

          do k=1,N-1
            do j=JR_RANGE
              do i=IR_RANGE
# ifdef NBQ_MASS
                qdmw_nbq(i,j,k)=(1.+0.5*(rho(i,j,k)+rho(i,j,k+1))/rho0) &
                            *0.5*wz(i,j,k,nrhs)*(hz(i,j,k)+hz(i,j,k+1))
# else
                qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*(hz(i,j,k)+hz(i,j,k+1))
# endif
              enddo
            enddo
          enddo
          k=0 
          do j=JR_RANGE 
            do i=IR_RANGE
# ifdef NBQ_MASS
              qdmw_nbq(i,j,k)=(1.+rho(i,j,k+1)/rho0) &
                              *0.5*wz(i,j,k,nrhs)*hz(i,j,k+1)
# else
              qdmw_nbq(i,j,k)= 0.5*wz(i,j,k,nrhs)*hz(i,j,k+1)
# endif
            enddo
          enddo
          k=N 
          do j=JR_RANGE
            do i=IR_RANGE
# ifdef NBQ_MASS
              qdmw_nbq(i,j,k)=(1.+rho(i,j,k)/rho0) &
                              *0.5*wz(i,j,k,nrhs)*hz(i,j,k)
# else
              qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*hz(i,j,k)
# endif
            enddo
          enddo

          call parallele_nbq(51)
          call parallele_nbq(151)

        endif ! iic
# endif

      endif  ! icall
      

      return
      end subroutine initial_nh_tile
#else
      subroutine initial_nh_empty
      return
      end
#endif
