! $Id: step3d_fast.F 1615 2014-12-17 13:27:07Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! ! and Rutgers University (Arango et al) are under MIT/X style license
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

#ifdef M3FAST

      subroutine step3d_fast (tile)
      implicit none
      integer tile, trd
# include "param.h"
# include "private_scratch.h"
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"

      trd=0

C$    trd=omp_get_thread_num()
!
      call step3d_fast_tile ( Istr,Iend,Jstr,Jend
     &                   , A2d(1, 1,trd), A2d(1, 2,trd), A2d(1, 3,trd)
     &                   , A2d(1, 4,trd), A2d(1, 5,trd), A2d(1, 6,trd)
     &                   , A2d(1, 7,trd), A2d(1, 8,trd)
     &                   , A2d(1, 9,trd), A2d(1,10,trd)
     &                   , A2d(1,11,trd), A2d(1,12,trd)
     &                   , A2d(1,13,trd)
#  ifdef M3FAST_C3D_UVSF
     &                   , A2d(1,14,trd), A2d(1,15,trd)
     &                   , A2d(1,16,trd), A2d(1,17,trd)
#   ifdef M3FAST_COUPLING3D
     &                   , A2d(1,18,trd), A2d(1,19,trd)
#   endif
#  endif
#  if defined M3FAST_UV && ! defined M3FAST_SEDLAYERS
     &                   , A3d(1, 1,trd), A3d(1, 2,trd)
     &                   , A3d(1, 3,trd), A3d(1, 4,trd)   ! CAUTION SDL
#  endif
     &                   )

      return
      end
! ! 
! !====================================================================
! !                      step3d_fast_tile
! !====================================================================
! !
       subroutine step3d_fast_tile (Istr,Iend,Jstr,Jend
     &                          ,Dnew,rubar,rvbar
     &                          ,Drhs, UFx,UFe
     &                          ,VFx,VFe
     &                          ,urhs,vrhs 
     &                          ,DUon,DVom ,work
#  ifdef M3FAST_C3D_UVSF
     &                          ,ru_ext_nbq_sum, rv_ext_nbq_sum
     &                          ,ru_ext_nbq_old, rv_ext_nbq_old
#   ifdef M3FAST_COUPLING3D
     &                          ,ru_int2d_nbq, rv_int2d_nbq
#   endif
#  endif
#  if defined M3FAST_UV && ! defined M3FAST_SEDLAYERS
     &                          ,Hzw_nbq_inv,  Hzr_nbq_inv 
     &                          ,Hzw_nbq_inv_u,  Hzw_nbq_inv_v  ! CAUTION SDL
#  endif
     &                          )
     
!
!
!***********************************************************************
!
!                    SOLVE FAST MODE 3D EQUATIONS
!
!***********************************************************************
!
!  This routines:
!  1- Computes non-NBQ RHS forcing terms of momentum equations. First
!     computes the barotropic (external) RHS forcing term (rubar,rvbar)
!     then adds it to the internal RHS forcing (computed in pre_step3d). 
!  2- Solves the 3D momentum conservation equations for fast-mode 
!     components (qdmu_nbq, qdmv_nbq, qdmw_nbq) by time integration of 
!     all forces:
!       Compressible pressure force + second viscosity + gravity
!       + NT Coriolis force + restoring force + non-NBQ RHS forces
!  3- Solves mass conservation equation, i.e., computes compressible 
!     density rho_nbq by time integratkdion of momentum divergence
!
!  In this version, a first guest of zeta is derived from the surface 
!  vertical velocity (surface characteristic relation) instead of the 
!  depth-averaged conservation of mass. This satisfies dynamical coupling 
!  with the surface layer. After solving the 3D momentum equations, a 
!  final zeta field is diagnozed from mass conservation (then Hz is also
!  corrected for the internal time step).
!
!  W-momentum equation is solved with explicit or implicit methods:
!  - Explicit scheme: w-momentum is updated right after (and the same 
!                     way as) u- and v-momentum.
!  - Implicit scheme: horizontal component of divergence is first 
!                     precomputed (as required by fast-mode mass 
!                     conservation) before tridiagonal Gauss Elimination 
!                     is carried out for qdmw_nbq(m).
!
!  For all components, a Forward-Backward scheme is implemented:
!  - Explicit scheme: Forward:  zeta, qdmu_nbq, qdmw_nbq. 
!                     Backward:       rho_nbq.
!  - Implicit scheme: Forward:  zeta, qdmu_nbq. 
!                     Backward:       qdmw_nbq, rho_nbq.
!
!  In the NBQ_PERF option, the vertical grid is not evolving at fast
! !  time step to gain computational time.
!
!***********************************************************************
!
      use debug
      implicit none
# include "param.h"

      integer Istr,Iend,Jstr,Jend, i,j,k, kbak,kold,
     &        imin,imax,jmin,jmax,
     &        kp1
# ifdef M3FAST_DDS0
       integer k1, k2
# endif
# ifdef MPI
     &       ,err
# endif
# ifdef PSOURCE
     &       ,is
# endif
      real    mybeta,myalpha,myepsilon,mygamma,
     &        VMAX,VMAXL, cff,cff0,cff1,cff2,cff3,cff4, 
     &        DUnew,DVnew, dum_s, dum2_s

# ifdef NBQ_THETAIMP
      real,parameter :: thetaimp_nbq = 0.5 
# else
      real,parameter :: thetaimp_nbq = 1.
# endif

# if defined SWASH || defined MOVING_BATHY || defined MILES\
  || defined KH_INST || defined REGIONAL
      real,parameter :: gammau=0.45   ! NBQ pressure gradient scheme
# else
      real,parameter :: gammau=1.  ! --> second-order
# endif
      real,parameter :: gammau_2=(1./3.)*(1.-gammau)
      real
     &                Dnew(PRIVATE_2D_SCRATCH_ARRAY),
     &               rubar(PRIVATE_2D_SCRATCH_ARRAY),
     &               rvbar(PRIVATE_2D_SCRATCH_ARRAY),
# ifdef M3FAST_BOTH
     &               rubarh(PRIVATE_2D_SCRATCH_ARRAY),
     &               rvbarh(PRIVATE_2D_SCRATCH_ARRAY),
# endif     
     &                Drhs(PRIVATE_2D_SCRATCH_ARRAY),
     &                 UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &                 UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &                 VFx(PRIVATE_2D_SCRATCH_ARRAY),
     &                 VFe(PRIVATE_2D_SCRATCH_ARRAY),
     &                urhs(PRIVATE_2D_SCRATCH_ARRAY),  
     &                vrhs(PRIVATE_2D_SCRATCH_ARRAY),
     &                DUon(PRIVATE_2D_SCRATCH_ARRAY),
     &                DVom(PRIVATE_2D_SCRATCH_ARRAY)
# ifdef M3FAST_BOTH
#  ifdef OPENACC
!$acc declare create( rubarh, rvbarh)
#  endif
# endif     
      real
     &                 work(PRIVATE_2D_SCRATCH_ARRAY)
# ifdef M3FAST_C3D_UVSF
#  ifdef M3FAST_COUPLING2D
      real 
     &       ru_ext_nbq(PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_ext_nbq(PRIVATE_2D_SCRATCH_ARRAY)
!$acc declare create ( ru_ext_nbq, rv_ext_nbq )
#  endif
      real 
     &       ru_ext_nbq_sum(PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_ext_nbq_sum(PRIVATE_2D_SCRATCH_ARRAY),
     &       ru_ext_nbq_old(PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_ext_nbq_old(PRIVATE_2D_SCRATCH_ARRAY)
#  ifdef M3FAST_COUPLING3D
      real
     &       ru_int2d_nbq  (PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_int2d_nbq  (PRIVATE_2D_SCRATCH_ARRAY)
#  endif
# endif
        real 
     &   Hzw_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl  :N),
     &   Hzr_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl+1:N),
!    &   Hzu_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl+1:N),
!    &   Hzv_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl+1:N),
     &   Hzw_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &   Hzw_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N)
!$acc declare create( Hzw_nbq_inv, Hzr_nbq_inv )
!$acc declare create( Hzw_nbq_inv_u, Hzw_nbq_inv_v )        
      real
     &    dthetadiv_nbqdz_u(PRIVATE_2D_SCRATCH_ARRAY,2),
     &    dthetadiv_nbqdz_v(PRIVATE_2D_SCRATCH_ARRAY,2),
     &              dZdxq_u(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &              dZdyq_v(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &                   FX(PRIVATE_2D_SCRATCH_ARRAY),
     &                   FY(PRIVATE_2D_SCRATCH_ARRAY),
     &                   FC(PRIVATE_1D_SCRATCH_ARRAY,-N_sl:N),
     &                   DC(PRIVATE_1D_SCRATCH_ARRAY,-N_sl:N),
     &                   CF(PRIVATE_1D_SCRATCH_ARRAY,-N_sl:N)
#  ifdef M3FAST_DDS0
     &         ,qdmw_nbq_old(PRIVATE_2D_SCRATCH_ARRAY,0:N)
#  endif
!$acc declare create( FX, FY )     
!$acc declare create( FC, CF )     
!$acc declare create( dZdxq_u, dZdyq_v )

# ifdef UV_COR_NT
      real
     &               ntcoru(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &               ntcorv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &               ntcorw(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N)
# endif 

#  ifndef NBQ_GRID_SLOW
#   define keepu1 FC3D
#   define keepv1 FC3D
#   define keepu2 FC3D
#   define keepv2 FC3D
#  endif
      
# ifdef M3FAST_ZETADISS
      real diss_zta
# endif
    
   !   real mytime(100)
# include "set_obc_definitions.h" 
! needed by nbq_bry_store, and Xnbq_bc
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
       real :: sum_nbq
# ifdef MRL_WCI
      real vstu,ustv,dudx,dvdx,dude,dvde
# endif
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
#ifdef OPENACC
# include "openacc.h" 
!!!! OpenACC boundary a enlever
# include "boundary.h"
#endif
# include "sources.h"
# ifdef AGRIF
#  include "zoom.h"
        integer irhox, irhoy, irhot
# endif
# if defined INTERNAL || defined BODYTIDE
      real,parameter :: U0 = 0.02
      real,parameter :: omega = 2.*pi/(12.4*3600.)
# endif
# ifdef M3FAST_SACOUS
      real dist_d  
      integer jvar1,jvar2
      real cffz(N)
# endif
# ifdef WET_DRY
      real cff1_WD,cff2_WD
# endif

# include "nbq.h"
# if defined MASS || defined MOVING_BATHY
      real*QUAD masstot
#  ifdef MPI
      real*QUAD allmasstot(1,NNODES)
#  endif
      integer ierr
# endif      
      integer ksl,kbot
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
# ifdef MPI
#  include "mpi_cpl.h"
      include 'mpif.h'
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif   

# include "compute_auxiliary_bounds.h"

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
!
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
# ifdef OPENACC
#  define exchange_r2d_tile exchange_r2d_tile_device
#  define exchange_u2d_tile exchange_u2d_tile_device
#  define exchange_v2d_tile exchange_v2d_tile_device
#  define exchange_r3d_tile exchange_r3d_tile_device
#  define exchange_u3d_tile exchange_u3d_tile_device
#  define exchange_v3d_tile exchange_v3d_tile_device
#  define exchange_w3d_tile exchange_w3d_tile_device
#  ifdef M3FAST_SEDLAYERS
#   define exchange_r3d_sedlay_tile exchange_r3d_sedlay_tile_device
#   define exchange_u3d_sedlay_tile exchange_u3d_sedlay_tile_device
#   define exchange_v3d_sedlay_tile exchange_v3d_sedlay_tile_device
#   define exchange_w3d_sedlay_tile exchange_w3d_sedlay_tile_device7
#   define exchange_r3d_sedlayb_tile exchange_r3d_sedlayb_tile_device
#   define exchange_u3d_sedlayb_tile exchange_u3d_sedlayb_tile_device
#   define exchange_v3d_sedlayb_tile exchange_v3d_sedlayb_tile_device
#   define exchange_w3d_sedlayb_tile exchange_w3d_sedlayb_tile_device
#
#  endif
#  define check_tab3d check_tab3d_device
#  define check_tab2d check_tab2d_device
# endif

# ifndef NHINT_3M
#  define FIRST_FAST_STEP_3M FIRST_FAST_STEP
# else
#  define FIRST_FAST_STEP_3M iif.eq.nsdtnbq 
# endif
!
!# define Hz_correct Hzw_nbq
!
!*********************************************************************
! Set several approximations to speed-up computation
!*********************************************************************
!
!# ifdef NBQ
!----------------------------------------------------------------------
! d./ds terms update frequency
!----------------------------------------------------------------------
!
#  ifdef NBQ_GRID_SLOW
!#   define NSTEP_DS mod(iic,1).eq.1 .and. iif.eq.1
!#   define NSTEP_DS .false. 
#   define NSTEP_DS FIRST_FAST_STEP_3M 
!#   define NSTEP_DS mod(iic,1).eq.0
#  else
#   define NSTEP_DS .true. 
#  endif
!
!----------------------------------------------------------------------
! Grid update frequency (Hz ...)
!----------------------------------------------------------------------
!
#  ifdef NBQ_GRID_SLOW
#   define NSTEP_GRID .or. iif.eq.nfast
#  else
#   define NSTEP_GRID .or. mod(iif,1).eq.0
#  endif
!# endif /* NBQ */
!
!----------------------------------------------------------------------
! Zeta: filter 
!----------------------------------------------------------------------
!
#  ifdef M3FAST_ZETADISS
       diss_zta = 0.e-5
#  endif  
!
!*********************************************************************
!
!                 EXTERNAL (2D) FAST MODE PROCESSING
!
!  Compute external (2D) rhs terms rubar,rvbar of 2D fast-mode equations
!  using Generalized Forward-Backward AB3-AM4 algorithm, and update 
!  internal and external forcing terms for 3D fast mode equations: 
!  ru_int, ru_ext~(rufrc+rubar).
!
!  1- AB3 forward step for D,ubar,vbar
!  2- Advance zeta(m+1)
!  3- AM4 backward step for rubar,rvbar 
!  4- update rufrc,rvfrc
!  5- update ru_int,rv_int
!  6- make some backups
!
! Reference for Generalized FB scheme:
! ------------------------------------
! Shchepetkin, A.F., and J.C. McWilliams, 2009: Computational kernel 
! algorithms for fine-scale, multiprocess, longtime oceanic simulations. 
! Pp. 119–182 in Handbook of Numerical Analysis: Computational Methods 
! for the Atmosphere and Oceans. R.M. Teman and J.J. Tribbia, eds, 
! Elsevier Science.
!
!*********************************************************************
!           call system_clock(cput1, count_rate, count_max)
!
!=====================================================================
! Store rho.h at first slow and fast time-step
!=====================================================================
!
#ifdef OPENACC
       if (FIRST_FAST_STEP) then
!$acc wait( sync_rho_rufrc_z_w )
! !$acc wait( sync_ruv_int_nbq )
!$acc update device(ru_int_nbq,rv_int_nbq,rw_int_nbq) !!iff=1
        endif
#endif

# if defined M3FAST_UV || defined M3FAST_W || defined M3FAST_RHO

! ! KERNEL_1  rho_grd <= ( rho )

!$acc kernels default(present)
!---------------------------------------------------------------------
! Initialisations
!---------------------------------------------------------------------
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
       
!---------------------------------------------------------------------
! Extrapolation in time:
!---------------------------------------------------------------------
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
! !====================================================================
! ! AB3 Forward Step: compute total depth of water column and 
! !                   vertically integrated mass fluxes which
! ! --- ------- ----   are needed to compute 
! ! rhs terms of the barotropic momentum equations (rubar,rvbar).
! !====================================================================
! !
! !--------------------------------------------------------------------
! !  Set indices to extrapolate (D,ubar,vbar) at m+1/2 (AB3)
! !--------------------------------------------------------------------
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
! 
! !
! !--------------------------------------------------------------------
! ! Extrapolate (D,ubar,vbar) at m+1/2
! !--------------------------------------------------------------------
! 
! ! Total depth/mass at m+1/2
!
! ! KERNEL_2  Drhs <= ( zeta, ubar, h  )
! ! KERNEL_2  urhs <= ( ubar   )
! ! KERNEL_2  DUon <= ( Drhs, on_u, urhs  )
! ! KERNEL_2  vrhs <= ( vbar   )
! ! KERNEL_2  DVom <= ( Drhs, om_v, vrhs  )

       if (FIRST_FAST_STEP) then
!$acc update device( h, zeta, ubar, vbar ) !iif=1
       endif
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
!
! Depth-average ubar velocity at m+1/2
!
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
!
! Depth-average vbar velocity at m+1/2
!
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
!
!----------------------------------------------------------------------
! Compute time averaged fields over all short timesteps.
!
! Reset/initialise arrays for averaged fields during the first
! barotropic time step; Accumulate averages after that. Include
! physical boundary points, but not periodic ghost points or
! computation  MPI computational margins.
!----------------------------------------------------------------------
!
#ifdef M3FAST_AVG_CLASSIC
# ifdef SOLVE3D
      cff1=weight(1,iif)
      cff2=weight(2,iif)
!$acc kernels if(compute_on_device) default(present) async(1)
      if (FIRST_2D_STEP) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
            DU_avg2(i,j)=cff2*DUon(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
          enddo
        enddo
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
# else      
!$acc kernels if(compute_on_device) default(present) async(1)
# endif
# endif /* M3FAST_AVG_CLASSIC */

# ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
      call check_tab2d(zeta(:,:,kstp),'zeta step3d_fast #1','r'
     &  ,ondevice=.TRUE.)
C$OMP END MASTER
# endif 
!
!-----------------------------------------------------------------------
! Load RHS values into additional AGRIF shared array for nesting
!-----------------------------------------------------------------------
! 
# ifdef AGRIF
      irhox = Agrif_Irhox()
      irhoy = Agrif_Irhoy()
      irhot = Agrif_Irhot()

      if (FIRST_FAST_STEP) then
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,0)=zeta(i,j,kstp)       
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
            DU_avg3(i,j,0)=DUon(i,j)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
            DV_avg3(i,j,0)=DVom(i,j)
          enddo
        enddo 
      endif

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Zt_avg3(START_2D_ARRAY,0))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg3(START_2D_ARRAY,0))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg3(START_2D_ARRAY,0))
#  endif

#  ifdef RVTK_DEBUG_ADVANCED
      if (.not.agrif_Root()) then
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(Zt_avg3(:,:,0),'Zt_avg3 (0) step3d_fast','r'
      &    ,ondevice=.TRUE.)
        call check_tab2d(DU_avg3(:,:,0),'DU_avg3 (0) step3d_fast','u'
      &    ,ondevice=.TRUE.)  
        call check_tab2d(DV_avg3(:,:,0),'DV_avg3 (0) step3d_fast','v'
      &    ,ondevice=.TRUE.)  
C$OMP END MASTER  
      endif
#  endif  
# endif /* AGRIF */
! !
! !====================================================================
! !
! ! Compute zeta(m+1) -- first guest -- and update grid
! ! 
! !====================================================================
! !
! !	
! !------------------------------------------------------------------
! ! Recover Hz at first fast step (if final Hz correction needed)
! !------------------------------------------------------------------
! !
# ifdef NBQ_HZCORRECT
!$acc kernels default(present)
! A TESTER .......
      if (iic.gt.ntstart .and. FIRST_FAST_STEP) then
        do k=-N_sl+1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              Hz(i,j,k) = Hz_correct(i,j,k)
            enddo
          enddo
        enddo
      endif
!! A TESTER .......
!$acc end kernels
# endif
! !
! !--------------------------------------------------------------------
! ! Get derived z grid variables at first fast step
! !--------------------------------------------------------------------
! !
# if ! defined M3FAST_COUPLING2D || defined M3FAST_UV || defined M3FAST_W
       if ( FIRST_FAST_STEP) then
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u  , Hzw_nbq_inv_v)
       endif
# endif
! !
! !-----------------------------------------------------------------
! ! zeta(m+1) is derived here from the surface vertical velocity 
! ! (surface characteristic relation) instead of depth-averaged 
! !  conservation of mass. This satisfies dynamical coupling with 
! !  the surface layer.
! !
! ! The surface kinematic relation is:
! !
! ! zeta(m+1) = zeta(m) + dtfast * [ wsurf_nbq(m)
! !                                 - usurf_nbq(m)*dzeta/dx(m+0.5)
! !                                 - vsurf_nbq(m)*dzeta/dy(m+0.5) ]
! !
! ! with zeta(m+0.5)=(1.5+beta)*zeta(m)-(0.5+2*beta)*zeta(m-1)
! !                                            +beta*zeta(m-2) --> AB3
! !-------------------------------------------------------------------
! ! Computes surface velocities
! !-------------------------------------------------------------------
! !
#include "m3fast_update_zeta.h"
! !
! !--------------------------------------------------------------------
! ! Compute time averaged fields over all short timesteps.
! !
! ! Reset/initialise arrays for averaged fields during the first
! ! barotropic time step; Accumulate averages after that. Include
! ! physical boundary points, but not periodic ghost points or
! ! computation  MPI computational margins.
! !--------------------------------------------------------------------
! !
!#if defined SOLVE3D  && defined HCOMP 
#if defined SOLVE3D   && defined toto
      cff1=weight(1,iif)
      cff2=weight(2,iif)
!$acc kernels default(present)
      if (FIRST_2D_STEP) then
! To be tested
        do j=JstrR,JendR
          do i=IstrR,IendR
!# ifndef M3SLOW_W         
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
!# endif        
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
            DU_avg2(i,j)=cff2*DUon(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
          enddo
        enddo
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
!# ifdef M3FAST_COUPLING2D          
!            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
!# endif            
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
!$acc end kernels      
#endif
! !
! !--------------------------------------------------------------------
! !  Update vertical grid
! !
! !  As soon as zeta(m+1) is known, the grid can be updated at fast
! !  step m+1. In PERF option, the grid is updated at lower frequency 
! !  and only during re-evaluation of zeta(m+1).
! !
! !  Caution: Hz_bak must be assigned only once in set_depth. The
! !  following code must thus be consistant with set_depth routine.
! !--------------------------------------------------------------------
! !
!# if defined NBQ && !defined NBQ_GRID_SLOW
# if !defined NBQ_GRID_SLOW
! !-------------------------------------------------------------------
! !Update moving bathymetry analytically
! !-------------------------------------------------------------------
! !
# ifdef MVB
         call ana_mvb_tile(Istr,Iend,Jstr,Jend)
# endif /* MVB */

        if ( (FIRST_TIME_STEP .and. FIRST_FAST_STEP) NSTEP_GRID ) then
!
! Update main grid parameters
!
! ! z_w, z_r, Hz, Hz_bak <= h, zeta, hinv, Hz 
        call set_depth_tile(Istr,Iend,Jstr,Jend)
!
! Update derived grid variables
!
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u, Hzw_nbq_inv_v)
     
      endif
# endif /* !NBQ_GRID_SLOW */
! !
! ! 
! !====================================================================
! ! AM4 backward step: compute depth-averaged RHS of 2D fast-mode 
! ! --- -------- ----  momentum equations
! !
! ! Compute here:
! ! - Surface pressure gradient
! ! - Horizontal advection terms for momentum equations (2D only)
! ! - Compute Coriolis (2D and 3D) term and advective curvilinear metric
! !    terms (2D only).
! ! - Linear and/or quadratic bottom stress.
! ! - Add 2D vortex-force terms combined with advection terms
! !====================================================================
! !
# include "m3fast_2Dmode.h"
! !
! !=====================================!==============================
! ! Compute rufrc & rvfrc: internal mode forcing for barotropic fast 
! ! mode
! !
! ! During the first fast time step convert rufrc & fvfrc into forcing
! ! terms by subtracting the fast-time "rubar" and "rvbar" from them;
! ! These forcing terms are then extrapolated forward in time using
! ! optimized Adams-Bashforth weights, so that the resultant rufrc
! ! and rvfrc are centered effectively at time n+1/2. From now on,
! ! these newly computed forcing terms will remain constant during
! ! the fast time stepping and will be added to "rubar" and "rvbar"
! ! during all subsequent fast time steps.
! !====================================================================
! !
      if (FIRST_FAST_STEP) then

! !    
! !  Pre-step2D: extrapolate rufrc
! !    
# include "m3fast_pre_step2d.h"
! !     
! !  Still first fast step
! !
#ifdef M3FAST_DRAFT 
         do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)
     &                        /(Drhs(i,j)+Drhs(i-1,j))
             enddo
           enddo
         enddo

         do k=1,N
           do j=JstrV,Jend
             do i=Istr,Iend
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)
     &                        /(Drhs(i,j)+Drhs(i,j-1))
             enddo
           enddo
         enddo
#endif
!#endif
!$acc end kernels
      endif   !<-- FIRST_FAST_STEP
! !
! !====================================================================
! !   Update internal and external forcing terms for NBQ mode
! !
! !   Compute external forcing terms ru_ext_nbq and updated internal 
! !   forcing terms ru_int_nbq for NBQ equations
! !
! !   ru_int_nbq     : RHS (3D) ( *mask & 2D correction)
! !   ru_ext_nbq     : RHS (2D)
! !   ru_ext_nbq_old : RHS (2D) at previous time-step
! !   ru_ext_nbq_sum : time-integrated RHS (2D)
! !====================================================================
! !
! !--------------------------------------------------------------------
! !  First fast time step only
! !--------------------------------------------------------------------
! !
# ifdef M3FAST_C3D_UVSF
       if (FIRST_FAST_STEP) then
      
! ! KERNEL_8  ru_ext_nbq_sum <= O
! ! KERNEL_8  ru_ext_nbq_old <= 0
      
!$acc kernels default(present)
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_ext_nbq_sum(i,j)=0.
            ru_ext_nbq_old(i,j)=0.
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_ext_nbq_sum(i,j)=0.
            rv_ext_nbq_old(i,j)=0.
          enddo
        enddo 
# if defined MASKING && defined M3FAST_2D
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)*umask(i,j)
            enddo
          enddo
        enddo
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)*vmask(i,j)
            enddo
          enddo
        enddo
# endif
!$acc end kernels
       endif ! FIRST_FAST_STEP
# endif
! !
! !--------------------------------------------------------------------
! !  All fast time steps
! !--------------------------------------------------------------------
! !

# ifndef M3FAST_COUPLING2D
#  define ru_ext_nbq UFx
# endif

! ! KERNEL_9  UFx <= ( rubar, ru_int2d_nbq, pm_u, Drhs, umask )
! ! KERNEL_9  ru_ext_nbq_old <= ( UFx, ru_ext_nbq_old )
! ! KERNEL_9  ru_ext_nbq_sum <= (ru_ext_nbq_sum, UFx ) 
! ! KERNEL_9  ru_int_nbq <= ( ru_int_nbq, ru_ext_nbq_old, Hz )
! ! KERNEL_9  ru_ext_nbq_old <= ( UFx )
! ! KERNEL_9  UFx <= ( rvbar, rv_int2d_nbq, pm_v, pn_v, Drhs, vmask )
! ! KERNEL_9  rv_ext_nbq_old <= ( UFx, rv_ext_nbq_old )
! ! KERNEL_9  rv_ext_nbq_sum <= ( rv_ext_nbq_sum, UFx )
! ! KERNEL_9  rv_int_nbq <= (rv_int_nbq, rv_ext_nbq_old,  Hz )
! ! KERNEL_9  rv_ext_nbq_old <= ( UFx )
! !       if ( FIRST_FAST_STEP) then
! ! !$acc update device( ru_int_nbq, rv_int_nbq )
! !       endif
#  if defined RVTK_DEBUG 
C$OMP BARRIER
C$OMP MASTER
!        call Check_tab3d(ru_int_nbq,'ru_int_nbq (A0)','uint')
!        call Check_tab3d(rv_int_nbq,'rv_int_nbq (A0)','vint')
       call check_tab3d(ru_int_nbq,'ru_int_nbq (A1)','uint')
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A1)','vint')
C$OMP END MASTER     
#  endif  
# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
!       call check_tab3d_sedlay(Hz,'Hz','r',N_sl+1,N)
!       call check_tab2d(ru_ext_nbq_old,'ru_ext_nbq_old (A1)','u')
!       call check_tab2d(rubar,'rubar (A1)','uint')
C$OMP END MASTER     
#  endif  
  
!$acc kernels default(present)
#ifdef M3FAST_C3D_UVSF
      do j=Jstr,Jend
        do i=IstrU,Iend
          ru_ext_nbq(i,j)=
# ifndef M3FAST_COUPLING3D
     &                     (rufrc(i,j)+rubar(i,j))     
# else
     &              (rubar(i,j)-ru_int2d_nbq(i,j)) ! 2* delta with slow 
# endif     
     &                                *pm_u(i,j)*pn_u(i,j)
     &                                /(Drhs(i,j)+Drhs(i-1,j))
#  ifdef MASKING
     &                                *umask(i,j)
#  endif
   !  &                    (rubar(i,j)+rufrc(i,j))
          ru_ext_nbq_old(i,j)=ru_ext_nbq(i,j)-ru_ext_nbq_old(i,j)  ! 2* delta (m+1 - m)
          ru_ext_nbq_sum(i,j)=ru_ext_nbq_sum(i,j)+ru_ext_nbq(i,j)
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)
     &                        +ru_ext_nbq_old(i,j)       ! 2* delta (m+1 - m)
     &                        *(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
       do j=Jstr,Jend
         do i=IstrU,Iend
           ru_ext_nbq_old(i,j)=ru_ext_nbq(i,j)   ! delta with slow
         enddo
       enddo
# endif /* M3FAST_C3D_UVSF */
!
# ifndef M3FAST_COUPLING2D
#  undef  ru_ext_nbq      
#  define rv_ext_nbq UFx
# endif
!
# ifdef M3FAST_C3D_UVSF
      do j=JstrV,Jend
        do i=Istr,Iend
          rv_ext_nbq(i,j)=
# ifndef M3FAST_COUPLING3D
     &                     (rvfrc(i,j)+rvbar(i,j))
# else
     &                    (rvbar(i,j)-rv_int2d_nbq(i,j))
# endif     
     &                                *pm_v(i,j)*pn_v(i,j)
     &                                /(Drhs(i,j)+Drhs(i,j-1))
#  ifdef MASKING
     &                                *vmask(i,j)
#  endif
          rv_ext_nbq_old(i,j)=rv_ext_nbq(i,j)-rv_ext_nbq_old(i,j)
          rv_ext_nbq_sum(i,j)=rv_ext_nbq_sum(i,j)+rv_ext_nbq(i,j)
        enddo
      enddo    
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)
     &                       +rv_ext_nbq_old(i,j)
     &                         *(Hz(i,j,k)+Hz(i,j-1,k))
          enddo
        enddo
      enddo
!#  if defined RVTK_DEBUG 
!C$OMP BARRIER
!C$OMP MASTER
!       call check_tab2d(ru_ext_nbq_old,'ru_ext_nbq_old (A3)','uint')
!       call check_tab3d(ru_int_nbq(:,:,1),'ru_int_nbq (A3)','uint')
!       call check_tab3d(rv_int_nbq,'rv_int_nbq (A3)','vint')
!C$OMP END MASTER     
!#  endif  
      do j=JstrV,Jend
        do i=Istr,Iend
          rv_ext_nbq_old(i,j)=rv_ext_nbq(i,j)
        enddo
      enddo
# endif /* M3FAST_C3D_UVSF */
!
# ifndef M3FAST_COUPLING2D
#  undef rv_ext_nbq  
# endif
!
!$acc end kernels
# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
       call check_tab3d(ru_int_nbq(:,:,1),'ru_int_nbq (A2)','uint'
     &    ,ondevice=.TRUE.)
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A2)','vint'
     &    ,ondevice=.TRUE.)
C$OMP END MASTER     
#  endif  

#  ifdef M3FAST_C3D_UVSF
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                        ru_int_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                        rv_int_nbq(START_2D_ARRAY,1))
#   endif 
#  endif 

# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
       call check_tab3d(ru_int_nbq,'ru_int_nbq (A)','uint'
     &    ,ondevice=.TRUE.)
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A)','vint'
     &    ,ondevice=.TRUE.)
C$OMP END MASTER     
#  endif  
! !
! !====================================================================
! !  Initializations and backups for 3D NBQ equations
! !====================================================================
! !
! !--------------------------------------------------------------------
! !  Initializations of ru_nbq_avg2 (as qdmu_nbq) at first fast step
! !  It is used at last fast step for computation of new ru_nbq_avg2
! !--------------------------------------------------------------------
! !

# ifdef M3FAST_C3D_UVFS
      if (FIRST_FAST_STEP) then

! ! KERNEL_10  ru_nbq_avg2 <= ( qdmu_nbq )
! ! KERNEL_10  rv_nbq_avg2 <= ( qdmv_nbq )

!$acc kernels default(present)
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_nbq_avg2(i,j,k)=qdmu_nbq(i,j,k)
            enddo
          enddo 
        enddo
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_nbq_avg2(i,j,k)=qdmv_nbq(i,j,k)
            enddo
          enddo 
        enddo
!$acc end kernels
      endif
# endif      
# ifdef M3FAST_C3D_WFS
      if (FIRST_FAST_STEP) then
!$acc update device( qdmw_nbq )   !! iif=1   

! ! KERNEL_11  rw_nbq_avg2 <= ( qdmw_nbq )

!$acc kernels default(present)
          do k=0,N
           do j=Jstr,Jend
            do i=Istr,Iend
              rw_nbq_avg2(i,j,k)=qdmw_nbq(i,j,k)
#  ifdef M3FAST_BOTH              
     &                          +qdmwh_nbq(i,j,k)
#  endif     
            enddo
#  ifdef M3FAST_BOTH               
         !    qdmwh2_nbq(i,j)=0.  ! DANGER change k
#  endif              
          enddo 
        enddo
!$acc end kernels
      endif    ! FIRST_FAST_STEP
# endif
! !
! !--------------------------------------------------------------------
! !  zw_nbq backups
! !--------------------------------------------------------------------
! !
# ifdef M3SLOW_W
      if (FIRST_TIME_STEP .and. FIRST_FAST_STEP) then
!$acc update device( z_w )   !! iif=1 ic=1

! ! KERNEL_12  zw_nbq <= ( z_w )

!$acc kernels default(present)
        do k=-N_sl,N
          do j=JstrR,JendR
            do i=IstrR,IendR
#  ifndef M3FAST_DDS0
              zw_nbq(i,j,k)=z_w(i,j,k)
#  else
              zw_nbq(i,j,k,:)=z_w(i,j,k)
#  endif
            enddo
          enddo
        enddo
!$acc end kernels
      endif
# endif  /* M3SLOW_W */

! !
! !--------------------------------------------------------------------
! !  Fast bottom friction
! !
! !  Set bottom stress using logarithmic or linear 
! !  and/or quadratic formulation. 
! !--------------------------------------------------------------------
! !
# if defined BSTRESS_FAST && !defined BBL
#  include "m3fast_fbf.h"
# endif
! !
! !--------------------------------------------------------------------
! !  Implicit part: system setup
!----------------------------------------------------------------------
! !    

# ifdef M3SLOW_W

! ! KERNEL_13  work <= ( pm, pn )
! ! KERNEL_13  DU_nbq <= 0
! ! KERNEL_13  DV_nbq <= 0

!$acc kernels default( present )
      do j=Jstr,Jend
        do i=Istr,Iend
          work(i,j)=pm(i,j)*pn(i,j)
          DU_nbq(i,j)=0.
          DV_nbq(i,j)=0.
        enddo
      enddo
!$acc end kernels
# endif /* M3SLOW_W */
!	
!------------------------------------------------------------------
!  Store qdmw_nbq into working array
!------------------------------------------------------------------
!
# ifdef M3FAST_DRAFT 
# ifdef NBQ
      if (LAST_FAST_STEP) then
        do k=0,N 
          do j=Jstr,Jend             
            do i=Istr,Iend
               rw_nbq(i,j,k)=qdmw_nbq(i,j,k) ! CAUTION qdmwh_nbq
            enddo
          enddo
        enddo
      endif
# endif
# endif
! !	
! !--------------------------------------------------------------------
! !  Recover Hz at first fast step (if final Hz correction needed)
! !--------------------------------------------------------------------
! !
!# ifdef NBQ_HZCORRECT
!!$acc kernels default( present )       
!        if (iic.gt.ntstart .and. iif.eq.1) then
!        do k=1,N
!          do j=JstrV-2,Jend+1
!            do i=IstrU-2,Iend+1
!              Hz(i,j,k) = Hz_correct(i,j,k)
!            enddo
!          enddo
!        enddo
!        endif
!!$acc end kernels
!# endif
! !
! !--------------------------------------------------------------------
! !  Store boundary values of nbq variables at previous 
! !  time-step for use in radiation boundary conditions
! !--------------------------------------------------------------------
! !
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
      call nbq_bry_store_tile (Istr,Iend,Jstr,Jend)
# endif

# ifdef NHINT_3M
      if (mod(iif,nsdtnbq).eq.0) then  ! HEREFA

       dtfast=dtfast*float(nsdtnbq)
       dtgrid_nbq=dtgrid_nbq*float(nsdtnbq)
       ndtfast=ndtfast/nsdtnbq
# endif
           
! !
! !********************************************************************
! !
! !                 SOLVE FAST MODE 3D NBQ EQUATIONS
! !
! !********************************************************************
! !
! ! Fast mode 3D momentum and mass-conservation equations can now be 
! ! solved.
! !
! ! W-momentum equation is solved with explicit or implicit methods:
! ! - Explicit scheme: w-momentum is updated right after (and the same 
! !                    way as) u- and v-momentum.
! ! - Implicit scheme: horizontal component of divergence is first 
! !                    precomputed (as required by fast-mode mass 
! !                    conservation) before tridiagonal Gauss  
! !                    Elimination is carried out for qdmw_nbq(m).
! !
! !  Compressible pressure-force and second viscosity are calculated 
! ! in thetadiv_nbq. Caution: this variable contains theta in the 
! ! first part of the algorithm and momentum divergence in the 
! ! remaining.
! !
! !  A Forward-backward scheme is implemented:
! !  - Explicit scheme: Forward:  zeta, qdmu_nbq, qdmw_nbq. 
! !                    Backward:       rho_nbq.
! !  - Implicit scheme: Forward:  zeta, qdmu_nbq. 
! !                     Backward:       qdmw_nbq, rho_nbq.
! !
! !********************************************************************
! !
! !
! !====================================================================
! ! Compute Non-traditional Coriolis pseudo-force
! !
! !   ntcoru = -e W cos(a) * Hz      [m2/s2]
! !   ntcorv = +e W sin(a) * Hz
! !   ntcorw = +e (U cos(a) - V sin(a) ) * Hz
! ! 
! ! with e = 2 Omega cos(Phi)
! !      a = angle between North and meridional grid axis
! !
! !====================================================================
!
# ifdef UV_COR_NT
#  include "m3fast_ntcoriolis.h"
# endif
! !
! !====================================================================
! !  Compute "Pressure - Viscosity" component (theta)
! !====================================================================
! ! 
# ifdef M3FAST_RHO 
!      if (FIRST_TIME_STEP.and.FIRST_FAST_STEP) then
! !     
! ! KERNEL_13b  thetadiv_nbq <= ( soundspeed2_nbq, rho_nbq, 
! ! KERNEL_13b  ------------ -- - Hzr_nbq_inv )     
! !
!$acc kernels default( present )    
      if (FIRST_TIME_STEP.and.FIRST_FAST_STEP_3M) then 
       do k=-N_sl+1,N
         do j=JstrV-2,Jend+1
           do i=IstrU-2,Iend+1
#  ifndef M3FAST_CSVISC2K
             thetadiv_nbq(i,j,k)=(
     &                     +rho_nbq(i,j,k) 
     &                           *soundspeed2_nbq
     &                          )*Hzr_nbq_inv(i,j,k) 
#  else
             thetadiv_nbq(i,j,k)=(
     &                     +rho_nbq(i,j,k) 
     &                           *soundspeed2_nbq(i,j,k)
     &                          )*Hzr_nbq_inv(i,j,k)
#  endif 
           enddo
         enddo
       enddo   
      else
       do k=-N_sl+1,N
         do j=JstrV-2,Jend+1
           do i=IstrU-2,Iend+1
#  ifndef M3FAST_CSVISC2K
             thetadiv_nbq(i,j,k)=( -visc2_nbq*( thetadiv_nbq(i,j,k)
#   ifdef M3FAST_DDS0
     &                                         +thetadiv3_nbq(i,j,k)
#   endif
     &                                        )
     &                     +rho_nbq(i,j,k) 
     &                           *soundspeed2_nbq
     &                          )*Hzr_nbq_inv(i,j,k) 
#  else
             thetadiv_nbq(i,j,k)=( -visc2_nbq(i,j,k)
     &                                *( thetadiv_nbq(i,j,k)
#   ifdef M3FAST_DDS0
     &                                  +thetadiv3_nbq(i,j,k)
#   endif
     &                                        )
     &                     +rho_nbq(i,j,k) 
     &                           *soundspeed2_nbq(i,j,k)
     &                          )*Hzr_nbq_inv(i,j,k)
#  endif 
           enddo
         enddo
       enddo     
      endif
!$acc end kernels  
# endif
! !
! !====================================================================
! !====================================================================
! !  Integrate fast-mode momentum equations
! !
! !  dqdm/dt=dtfast*(Compressible pressure force + second viscosity
! !                  + gravity + non-NBQ forces + nudging)
! !
! !  If explicit scheme: all (x,y,z) mom components are done here
! !  If implicit scheme: z-mom component is done after
! !====================================================================
! !====================================================================
! ! 
! !====================================================================
! ! U and V fast-momentum
! !====================================================================
! !
#  ifdef M3FAST_DDS0
#   include "m3fast_update_qdmuv_dds0.h"
#  else
#   include "m3fast_update_qdmuv_dds1.h"
#  endif
! !
! !====================================================================
! !      Compute horizontal part of divergence operator
! !              if implicit scheme
! !  (if so, from now on, thetadiv_nbq array is used
! !               for divergence, not theta)
! !====================================================================
! !
# ifdef NBQ_IMP
#  ifdef M3FAST_DDS0
#   include "m3fast_divh_dds0.h"
#  else
#   include "m3fast_divh_dds1.h"
#  endif
# endif
! !
# ifdef M3FAST_W
! !
! !====================================================================
! !      Fast mode W-Momentum equation: qdmw_nbq 
! !        - Explicit integration: compute qdmw_nh
! !        - Implicit integration: compute RHS of tri-diag system 
! !====================================================================
! !
#  ifdef M3FAST_DDS0
#   include "m3fast_update_qdmw_dds0.h"
#  else
#   include "m3fast_update_qdmw_dds1.h"
#  endif
! !
! !====================================================================
! !   Implicit integration of W-momentum
! !====================================================================
! !
! !   Initialisation FX and FYX out here
! !
! ! KERNEL_20  dZdxq_u <= ( z_r, qdmu_nbq ) 
! ! KERNEL_20  dZdxq_w <= ( z_w, qdmu_nbq, Hz )
! ! KERNEL_20  dZdxq_w <= ( Hzw_nbq_inv_u, dZdxq_u, umask )
! ! KERNEL_20  FX <= (pm_u, dZdxq_w, umask )
! ! KERNEL_20  thetadiv_nbq <= ( FX )
! !
! !
! ! Compute and solve tridiag (implicit) system for W
! !
#  ifdef NBQ_IMP
#   ifdef M3FAST_DDS0
#    include "m3fast_qdmw_implicit_dds0.h"
#   else
#    include "m3fast_qdmw_implicit_dds1.h"
#   endif                                  
#  endif                                  
! !
! !====================================================================
! !  W-momentum open boundary conditions
! !====================================================================
! !
#  ifdef OBC_NBQ
      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
#  endif
! !
      if (LAST_FAST_STEP) then
!$acc update host( qdmw_nbq )   !! iif=last
      endif 
! !
! !====================================================================
! !  Exchange qdmw (if necessary)
! !====================================================================
! ! 
#  ifdef UV_COR_NT
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_w3d_sedlayb_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmw_nbq(START_2D_ARRAY,-N_sl))   
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmw_nbq(START_2D_ARRAY,0))   
#   endif
#  endif
#  ifdef M3FAST_BOTH
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmwh_nbq(START_2D_ARRAY,0))     ! CAUTION needed ?
#   endif
#  endif
! !
! !====================================================================
! !      Compute horizontal components of divergence operator
! !              if explicit scheme
! !  (if so, from now on, thetadiv_nbq array is used
! !               for divergence, not theta)
! !====================================================================
! !
#  ifndef NBQ_IMP
#  ifdef M3FAST_DDS0
#   include "m3fast_divh_dds0.h"
#  else
#   include "m3fast_divh_dds1.h"
#  endif
#  endif
! !
! !====================================================================
! !====================================================================
! !      Fast-mode conservation of mass
! !====================================================================
! !====================================================================
! !
! !
! !====================================================================
! !         Vertical part of divergence
! !      div(m+1) = div(m+1) + divZ(qdmw_nbq(m+1))
! !====================================================================
! !		
! !.......Computes fluxes:  
! !
! ! KERNEL_25  FC <= ( Hzw_nbq_inv, qdmw_nbq )
! ! KERNEL_25  thetadiv_nbq <= ( thetadiv_nbq, FC )
! !
#  include "m3fast_divv.h"
! !
! !====================================================================
! !         Exchange divergence
! !====================================================================
! !	
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   ifndef M3FAST_SEDLAYERS
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv_nbq(START_2D_ARRAY,1))
#     if defined NBQ_HZ_PROGNOSTIC  || defined M3FAST_DDS0 
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv2_nbq(START_2D_ARRAY,1))
#     endif
#     ifdef M3FAST_DDS0 
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv3_nbq(START_2D_ARRAY,1))
#     endif
#   else
      call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv_nbq(START_2D_ARRAY,-N_sl+1))
#     if defined NBQ_HZ_PROGNOSTIC  || defined M3FAST_DDS0 
       call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv2_nbq(START_2D_ARRAY,-N_sl+1))
#     endif
#   endif
#  endif

# endif  /* M3FAST_W */
! !
! !
! !=====================================================================
! !      Solve rho_nbq from mass conservation 
! !
! !      using a Forward-Backward scheme:
! !          rho_nbq(m+1) = rho_nbq(m) - dtfast*DIV(m+1)
! !====================================================================
! !
# ifdef M3FAST_RHO
!$acc kernels default( present )

! ! KERNEL_26 rho_nbq <= ( rho_nbq, thetadiv_nbq ) 
! ! KERNEL_26 thetadiv_nbq <= (  thetadiv_nbq, 
! ! KERNEL_26 ------------ -- - soundspeed2_nbq, rho_nbq, Hzr_nbq_inv )
      do k=-N_sl+1,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            rho_nbq(i,j,k) = rho_nbq(i,j,k)  
#   ifndef M3FAST_DDS0
     &                       - dtfast*thetadiv_nbq(i,j,k)
#   else
     &                       - dtfast*(thetadiv_nbq(i,j,k)+
     &                                 thetadiv3_nbq(i,j,k))
#   endif

!#   ifndef M3FAST_CSVISC2K
!#    ifndef M3FAST_DDS0
!             thetadiv_nbq(i,j,k)=(-visc2_nbq*thetadiv_nbq(i,j,k)
!     &                            +soundspeed2_nbq*rho_nbq(i,j,k)) 
!     &                           *Hzr_nbq_inv(i,j,k)  
!#    else
! !           thetadiv_nbq(i,j,k)=( -visc2_nbq*( thetadiv_nbq(i,j,k)
! !    &                                       +thetadiv3_nbq(i,j,k))
! !    &                         +soundspeed2_nbq*rho_nbq(i,j,k) )
! !    &                                     *Hzr_nbq_inv(i,j,k) 
!#    endif
!#   else
!#    ifndef M3FAST_DDS0
!             thetadiv_nbq(i,j,k)=(-visc2_nbq(i,j,k)*thetadiv_nbq(i,j,k)
!     &                            +soundspeed2_nbq(i,j,k)
!     &                             *rho_nbq(i,j,k))
!     &                           *Hzr_nbq_inv(i,j,k)  
!#    else
!  !          thetadiv_nbq(i,j,k)=( -visc2_nbq(i,j,k)*( thetadiv_nbq(i,j,k)
!  !   &                                       +thetadiv3_nbq(i,j,k))
!  !   &                         +soundspeed2_nbq(i,j,k)*rho_nbq(i,j,k) )
!  !   &                                     *Hzr_nbq_inv(i,j,k) 
!#    endif
!#   endif
          enddo
        enddo
      enddo
!$acc end kernels
! !
! !-------------------------------------------------------------------
! !  BC on rho_nbq
! !-------------------------------------------------------------------
! !
#   ifdef M3FAST_DDS0
        call rnbq_bc_tile(Istr,Iend,Jstr,Jend, work) !! leave commented
#   endif
! !
! !-------------------------------------------------------------------
! !  Acoustic wave emission
! !-------------------------------------------------------------------
!!
#  ifdef M3FAST_SACOUS
#   include "m3fast_sacous.h"
#  endif
! !
! !-------------------------------------------------------------------
! !  rhobar_nbq: depth-mean density (/rho0)
! !-------------------------------------------------------------------
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ATTENTION exchange !!!!!!! FRANCIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend, 
     &                           rhobar_nbq(START_2D_ARRAY,knew))
#   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   ifdef RVTK_DEBUG
!       call check_tab2d(rhobar_nbq(:,:,knew),'rhobar_nbq','r')
#   endif    
#  endif /* NBQ_MASS */

# endif  /* M3FAST_RHO */

# ifdef NHINT_3M
        dtfast=dtfast/float(nsdtnbq)
        dtgrid_nbq=dtgrid_nbq/float(nsdtnbq)
        ndtfast=ndtfast*nsdtnbq
       else  !HEREFA
! !
! !--------------------------------------------------------------------
! ! Update ubar     TBT   06/02/2023
! !--------------------------------------------------------------------
! !

!# ifndef M3FAST_COUPLING2D

! ! KERNEL_17  ubar <= ( urhs )
! ! KERNEL_17  vbar <= ( vrhs )

#ifdef NBQ
!$acc kernels default( present )
      do j=Jstr,Jend
        do i=IstrU,Iend
          ubar(i,j,knew)=urhs(i,j)
        enddo 
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          vbar(i,j,knew)=vrhs(i,j)
        enddo 
      enddo
!$acc end kernels
#endif

# if defined M2_FRC_BRY || defined Z_FRC_BRY
      if (FIRST_FAST_STEP) then
!$acc update device(                    
#   ifdef M2_FRC_BRY
!$acc&        ubarbry_west, ubarbry_east
#   endif
#  ifdef Z_FRC_BRY
!$acc&       ,zetabry_west, zetabry_east
#  endif
!$acc&              ) !! iif=1
      endif
# endif      
! !
! !--------------------------------------------------------------------
! !  U & V momentum open boundary conditions
! !--------------------------------------------------------------------
! !
!# ifdef NBQ 
!     M2bc_nbq_flag=.true. ! apply boundary wet/dry conditions
                           ! and compute DU_nbq
!     call u2dbc_tile   (Istr,Iend,Jstr,Jend, work)
!     call v2dbc_tile   (Istr,Iend,Jstr,Jend, work)

     ! call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
     ! call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!# endif /* ! M3FAST_COUPLING2D */
       endif !HEREFA
# endif
! !
! !
! !
! !********************************************************************
! !
! !
! !                     FAST-MODE POST-PROCESSING
! !
! !
! !********************************************************************
! !
! !====================================================================
! !  Update total mass of water volume Dnew 
! !====================================================================
! !
! ! KERNEL_27  Dnew <= ( zeta, h ) 

!$acc kernels default( present )
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
# ifdef NBQ_MASS
     &               *rhobar_nbq(i,j,knew)
# endif
        enddo
      enddo
!$acc end kernels 
# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
c       call check_tab3d(thetadiv_nbq,'step3d_fastthetadiv_nbq','rint')
# ifdef M3FAST_SEDLAYERS
      call check_tab3d_sedlay(rho_nbq,'step3d_fastrho_nbq_1',
     &  'r',N_sl+1,N,ondevice=.TRUE.)
# else
      call check_tab3d(rho_nbq,'step3d_fastrho_nbq_1','r',ondevice=.TRUE.)
# endif      
C$OMP END MASTER
# endif    

! !====================================================================
! !  Update W NBQ forcing for internal mode (rw_nbq ~ rw_nbq_avg1)
! !  Note: here rw_nbq contains qdmw_nbq(m)
! !====================================================================
! !
# ifdef M3FAST_W

#  ifdef NHINT 
! !
! !--------------------------------------------------------------------
! ! NHINT numerical mode control, v2:
! !--------------------------------------------------------------------
! !
            do j=Jstr,Jend
             do i=Istr,Iend
             
              k=N
              cff2=qdmw_nbq(i,j,k)*alphaw_nbq-qdmw_nbq(i,j,N)
              qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)+cff2
              
              do k=N-1,N-2*alphaNw_nbq,-1
               cff3=cff2
               cff=(1.-alphaw_nbq)
     &               *exp(-(z_w(i,j,k)            -z_w(i,j,N))**2
     &                    /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2)

               cff2=
     &    -      cff*qdmw_nbq(i,j,k)
     &    -      cff *(qdmw_nbq(i,j,k)-qdmw_nbq(i,j,N)
     &         *(z_w(i,j,k)+H(i,j))/(z_w(i,j,N)+H(i,j))
     &         *0.5*(Hz(i,j,k)+Hz(i,j,k+1))/Hz(i,j,N))

               qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)+cff2

                rho_nbq(i,j,k+1)=rho_nbq(i,j,k+1)
     &           -dtfast*( cff3*Hzw_nbq_inv(i,j,k+1)
     &                    -cff2*Hzw_nbq_inv(i,j,k  ))
             enddo
           enddo
         enddo
     
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq(START_2D_ARRAY,1))   
#   endif

#   endif /* NHINT  */
      if (LAST_FAST_STEP) then
      
! ! KERNEL_28  rw_nbq <= ( qdmw_nbq, rw_nbq, rw_int_nbq, pm, pn )

!$acc kernels default( present )
        do k=0,N 
          do j=Jstr,Jend              
            do i=Istr,Iend
              rw_nbq(i,j,k)=((qdmw_nbq(i,j,k)
#  ifdef M3FAST_BOTH               
     &                       +qdmwh_nbq(i,j,k)
#  endif              
     &                    -rw_nbq(i,j,k))
     &                       /dtfast
#  ifdef NHINT_3M
     &                       /float(nsdtnbq)
#  endif
#   ifdef M3FAST_C3D_WSF
     &                     -rw_int_nbq(i,j,k)
#   endif     
     &                                  )/(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo
!$acc end kernels 
!$acc update host( rw_nbq )         !! iif=last
      endif
!# if defined RVTK_DEBUG && defined NBQ
!      call check_tab3d(rw_nbq,'rw_nbq','rint')
!#  endif  
# endif
! !
! !====================================================================
! ! Get filtered RHS terms
! ! and multiply by dx*dy to get units of rho*Hz*dx*dy*ru
! !====================================================================
! !
      if (LAST_FAST_STEP) then
! !
! !--------------------------------------------------------------------
! ! Store average fields AVG1 of rho and rhobar
! !--------------------------------------------------------------------
! !     
# ifdef NBQ_MASS
        do j=Jstr,Jend
          do i=Istr,Iend
            rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,knew)
          enddo
        enddo 
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rho_nbq_avg1(i,j,k)=1.d0
     &                            + ( rho_nbq(i,j,k)/Hzr(i,j,k)
     &                               +rho_grd(i,j,k) )  
            enddo
          enddo 
        enddo
# endif /* NBQ_MASS */
! !
! !--------------------------------------------------------------------
! !  Compute average fields AVG2 of RHS NBQ forcing
! !  Note: here ru_nbq_avg2, ru_nbq_2d_old ... are working arrays
! !--------------------------------------------------------------------
! !
!$acc kernels default( present )
# ifdef M3FAST_UV 

! ! KERNEL_29  ru_int_nbq <= ( ru_int_nbq, ru_ext_nbq_old, Hz )
! ! KERNEL_29  ru_nbq_avg2 <= ( qdmu_nbq, ru_nbq_avg2, ru_int_nbq,
! ! KERNEL_29  ----------- -- - ru_ext_nbq_sum, Hz, on_u, om_u )
! ! KERNEL_29  rv_int_nbq <= ( rv_int_nbq, rv_ext_nbq_old, Hz )
! ! KERNEL_29  rv_nbq_avg2 <= ( qdmv_nbq, rv_nbq_avg2, rv_int_nbq,
! ! KERNEL_29  ----------- -- - rv_ext_nbq_sum, Hz, on_v, om_v )
! ! KERNEL_29  rw_nbq_avg2 <= ( qdmw_nbq, rw_nbq_avg2, 
! ! KERNEL_29  ----------- -- - rw_int_nbq, on_r, om_r)

       do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend 
!             ru_nbq_avg1(i,j,k)= ru_nbq(i,j,k)
#  ifdef M3FAST_C3D_UVSF
              ru_int_nbq(i,j,k) = ru_int_nbq(i,j,k)
#   ifndef M3FAST_COUPLING2D
     &                     -rubar(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))   ! CAUTION HERE
#   else
     &                     -ru_ext_nbq_old(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
#   endif     
#  endif     
#  ifdef M3FAST_C3D_UVFS
              ru_nbq_avg2(i,j,k)=
     &                       ((qdmu_nbq(i,j,k)-ru_nbq_avg2(i,j,k))/dt  ! CAUTION: use define
#   ifdef M3FAST_C3D_UVSF
     &                   -ru_int_nbq(i,j,k)-(ru_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i-1,j,k))
#   endif
     &                       )*on_u(i,j)*om_u(i,j)
#  endif     
            enddo
          enddo 
        enddo  
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend             
    !          rv_nbq_avg1(i,j,k)= rv_nbq(i,j,k) 
#  ifdef M3FAST_C3D_UVSF
              rv_int_nbq(i,j,k) = rv_int_nbq(i,j,k)
#   ifndef M3FAST_COUPLING2D
     &                     -rvbar(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
#   else
     &                     -rv_ext_nbq_old(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
#   endif
#  endif     
#  ifdef M3FAST_C3D_UVFS
              rv_nbq_avg2(i,j,k)=
     &                        ((qdmv_nbq(i,j,k)-rv_nbq_avg2(i,j,k))/dt ! CAUTION: use define
#   ifdef M3FAST_C3D_UVSF
     &                   -rv_int_nbq(i,j,k)-(rv_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i,j-1,k))
#   endif
     &                        )*on_v(i,j)*om_v(i,j)
#  endif    
            enddo
          enddo 
        enddo
# endif /* M3FAST_UV */
# ifdef M3FAST_W
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
!             rw_nbq_avg1(i,j,k)= rw_nbq(i,j,k)
              rw_nbq_avg2(i,j,k)=
     &                          ((qdmw_nbq(i,j,k)
#   ifdef M3FAST_BOTH
     &        +qdmwh_nbq(i,j,k)
#   endif
     &        -rw_nbq_avg2(i,j,k))/dt
#   ifdef M3FAST_C3D_WSF
     &                           -rw_int_nbq(i,j,k)
#   endif
     &                          )*on_r(i,j)*om_r(i,j)
            enddo
          enddo 
        enddo
# endif /* M3FAST_W */
!$acc end kernels
      endif ! LAST_FAST_STEP
! !
! !--------------------------------------------------------------------
! ! Dismiss coupling of NBQ, NBQ2EXT & NBQ2INT to debug
! !--------------------------------------------------------------------
! !
# ifdef NBQ_NOCOUPLING 
      ru_nbq      =0.   ! 3D
      rv_nbq      =0.
      ru_nbq_avg2 =0.
      rv_nbq_avg2 =0.
#  ifdef M3FAST 
      rw_nbq      =0.
      rw_nbq_avg2 =0.
#  endif
#  ifdef NBQ_MASS
      rhobar_nbq  =1.
      rho_nbq     =0.
      rho_nbq_avg1=1.
#  endif
# endif /* NBQ_NOCOUPLING */
! !
! !--------------------------------------------------------------------
! ! Exchange NBQ coupling
! !--------------------------------------------------------------------
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifdef NBQ_MASS
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        rhobar_nbq_avg1(START_2D_ARRAY))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq_avg1(START_2D_ARRAY,1))
#  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ATTENTION exchange !!!!!!! FRANCIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  ifdef M3FAST_C3D_UVFS
      if (LAST_FAST_STEP) then
       call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                         ru_nbq(START_2D_ARRAY,1))
       call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rv_nbq(START_2D_ARRAY,1))
       call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                         ru_nbq_avg2(START_2D_ARRAY,1))
       call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rv_nbq_avg2(START_2D_ARRAY,1))
      endif
#  endif
!#  ifdef NBQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ATTENTION exchange !!!!!!! FRANCIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  ifdef M3FAST_C3D_WFS
      if (LAST_FAST_STEP) then
       call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rw_nbq(START_2D_ARRAY,0))
       call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rw_nbq_avg2(START_2D_ARRAY,0))
      endif
#  endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# endif

# ifdef M3FAST_UV 
      if (LAST_FAST_STEP) then
!$acc update host( ru_int_nbq, ru_nbq_avg2      !! iif=last
!$acc&            ,rv_int_nbq, rv_nbq_avg2 
!$acc&            ,rw_nbq_avg2
!$acc&           ) !async(sync_ruv_nbq_avg2)
! ! ru_int_nbq => pre_step
! ! ru_nbq_avg2 => step3d_uv1
      endif
#endif

! !
! !====================================================================
! !  Depth-averaged velocity & forcing from fast mode
! !====================================================================
! !
! ! Output: (ubar,vbar), (DU_avg1,DV_avg1)
! !
! ! Cytil Attention sert a rien dans le cas NBQ
!# ifdef M3FAST_ZETAW
!$acc kernels default( present )
#  define Dstp DUon
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
#  ifndef MVB        
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
#  else
          Dstp(i,j)=zeta(i,j,kstp)+dh_mvb(i,j,kstp2)
#  endif          
        enddo
      enddo
!$acc end kernels
!# endif /* M3FAST_ZETAW */
      cff=0.5*dtfast
    !  cff=2.*dtfast
      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)
      
! #ifdef RVTK_DEBUG
! C$OMP BARRIER
! C$OMP MASTER
!        call check_tab2d(ubar(:,:,knew),'ubar st_fast_a','u')
!        call check_tab2d(ubar(:,:,kstp),'ubarkstp    _a','u')
!        call check_tab2d(rubar(:,:),'rubar st_fast_a','uint')
!        call check_tab2d(rvbar(:,:),'rvbar st_fast_a','vint')
! !       call check_tab2d(rufrc(:,:),'rufrc st_fast_a','u')
! !       call check_tab2d(rubarh(:,:),'rubarh st_fast_a','uint')
!        call check_tab2d(DUon(:,:),'DUon st_fast_a','r')
!        call check_tab2d(pn(:,:),'pn st_fast_a','r')
!        call check_tab2d(pm(:,:),'pm st_fast_a','r')
!        call check_tab2d(rvfrc(:,:),'rvfrc st_fast_a','v')
! #endif
         
!$acc kernels default( present )
# ifndef M3FAST_AVG_CLASSIC
      if (iif.eq.1) then
! ! KERNEL_30  DU_avg1 <= 0    
! ! KERNEL_30  DU_avg2 <= 0   

         DU_avg1=0.
         DV_avg1=0.
         DU_avg2=0.
         DV_avg2=0.

      endif
# endif


! ! KERNEL_31  DUnew <= ( DU_nbq )
! ! KERNEL_31  ubar <= ( DUnew, Dnew )
! ! KERNEL_31  DU_avg1 <= ( DU_avg1, on_u, DUnew )
! ! KERNEL_31  DU_avg2 <= ( DU_avg2, on_u, DUnew )
! ! KERNEL_31  DVnew <= ( DV_nbq )
! ! KERNEL_31  vbar <= ( DVnew, Dnew )
! ! KERNEL_31  DV_avg1 <= ( DV_avg1, om_v, DVnew )
! ! KERNEL_31  DV_avg2 <= ( DV_avg2, om_v, DVnew )
! 
      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef M3FAST_ZETAW
           DUnew=DU_nbq(i,j) *2. 
# else
          DUnew=( (Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)
     &     +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                                *( rubar(i,j)
     &                                  +rufrc(i,j)
     &                          !   +ru_int2d_nbq(i,j)
#  ifdef M3FAST_BOTH
     &                                  +rubarh(i,j)
#  endif     
     &                                 ))
#  ifdef MASKING
     &                                     *umask(i,j)
#  endif
# endif
# ifdef WET_DRY
     &                                     *umask_wet(i,j)
# endif
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) 
     &                             +cff1*on_u(i,j)*( DUnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#  ifdef WET_DRY
     &  *umask_wet(i,j)
#  endif
# endif
     &                                                   )
# ifndef M3FAST_AVG_CLASSIC
          DU_avg2(i,j)=DU_avg2(i,j)+cff2*on_u(i,j)*( DUnew 
#  ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#   ifdef WET_DRY
     & *umask_wet(i,j)
#   endif
#  endif
     &                                                   )
# endif /*   M3FAST_AVG_CLASSIC */
        enddo
      enddo 
      do j=JstrV,Jend
        do i=Istr,Iend
# ifdef M3FAST_ZETAW
          DVnew=DV_nbq(i,j) *2.
# else
          DVnew=( (Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)
     &     +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                                *( rvbar(i,j)
     &                                  +rvfrc(i,j)
     &                          !   +rv_int2d_nbq(i,j)
#  ifdef M3FAST_BOTH
     &                                  +rvbarh(i,j)
#  endif     
     &                                 ))
#  ifdef MASKING
     &                                     *vmask(i,j)
#  endif
# endif
# ifdef WET_DRY
     &                                     *vmask_wet(i,j)
# endif
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) 
     &                              +cff1*om_v(i,j)*(DVnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
#  ifdef WET_DRY
     & *vmask_wet(i,j)
#  endif
# endif
     &                                                   )
# ifndef M3FAST_AVG_CLASSIC
          DV_avg2(i,j)=DV_avg2(i,j)+cff2*om_v(i,j)*( DVnew
#  ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
#   ifdef WET_DRY
     & * vmask_wet(i,j)
#   endif
#  endif
     &                                                   )
# endif /*   M3FAST_AVG_CLASSIC */
        enddo
      enddo

!$acc end kernels

! !
! !--------------------------------------------------------------------
! ! Apply point sources for hydrostatic case
! !--------------------------------------------------------------------
! !
# if defined PSOURCE && !defined M3FAST 
      do is=1,Nsrc 
#  ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#  else
        i=Isrc(is)
        j=Jsrc(is)
#  endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            ubar(i,j,knew)=2.*Qbar(is)/( on_u(i,j)
     &                       *(Dnew(i-1,j)+Dnew(i,j)) )
            DU_avg1(i,j,nnew)=Qbar(is)
          else
            vbar(i,j,knew)=2.*Qbar(is)/( om_v(i,j)
     &                       *(Dnew(i,j-1)+Dnew(i,j)) )
            DV_avg1(i,j,nnew)=Qbar(is)
          endif
        endif
      enddo
# endif
! !
! !---------------------------------------------------------------------
! !  Set 2D Momemtum nudging
! !---------------------------------------------------------------------
! !
# if defined M2NUDGING && defined M2CLIMATOLOGY

#  ifdef ZONAL_NUDGING
      if (FIRST_TIME_STEP .or. mod(iic,10).eq.0) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew),ubzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew),vbzon)
        endif
      endif
      if (FIRST_TIME_STEP) then
        if (FIRST_FAST_STEP) then 
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubclm(START_2D_ARRAY),ubclmzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbclm(START_2D_ARRAY),vbclmzon)
        endif
      endif
#  endif /* ZONAL_NUDGING */

      do j=Jstr,Jend
        do i=IstrU,Iend
#  ifdef ZONAL_NUDGING        
          DUnew = dtfast*M2nudgcof(i,j)*(ubclmzon(j)-ubzon(j))
#  else          
          DUnew = dtfast*M2nudgcof(i,j)*(ubclm(i,j)-ubar(i,j,knew))
#  endif
#  ifdef MASKING
     &                 * umask(i,j)
#  endif 
#  ifdef WET_DRY
     &                 * umask_wet(i,j)
#  endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# ifndef M3FAST_AVG_CLASSIC
          DU_avg2(i,j)     =DU_avg2(i,j)      +cff2*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# endif /*   M3FAST_AVG_CLASSIC */    
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
#  if defined ZONAL_NUDGING  
          DVnew = dtfast*M2nudgcof(i,j)*(vbclmzon(j)-vbzon(j)) 
#  else
          DVnew = dtfast*M2nudgcof(i,j)*(vbclm(i,j)-vbar(i,j,knew))    
#  endif      
#  ifdef MASKING
     &                 * vmask(i,j)
#  endif 
#  ifdef WET_DRY
     &                 * vmask_wet(i,j)
#  endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# ifndef M3FAST_AVG_CLASSIC
          DV_avg2(i,j)     =DV_avg2(i,j)      +cff2*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# endif /*   M3FAST_AVG_CLASSIC */
        enddo
      enddo
# endif /* M2NUDGING */
! !
! !-------------------------------------------------------------------- 
! ! Set boundary conditions and compute integral mass flux accross
! ! all open boundaries, if any.
! !--------------------------------------------------------------------
! !
!# ifdef NBQ
      M2bc_nbq_flag=.false.  ! skip wet/dry conditions
  !                           ! and DU_nbq computation
!# endif
      call u2dbc_tile (Istr,Iend,Jstr,Jend, work) 
      call v2dbc_tile (Istr,Iend,Jstr,Jend, work)
# ifdef WET_DRY
#  ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
        DO j=Jstr,Jend
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask_wet(Istr,j)
#   ifdef MRL_WCI
          ust2d(Istr,j)=ust2d(Istr,j)*umask_wet(Istr,j)
#   endif
        END DO
        DO j=JstrV,Jend
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask_wet(Istr-1,j)
#   ifdef MRL_WCI
          vst2d(Istr-1,j)=vst2d(Istr-1,j)*vmask_wet(Istr-1,j)
#   endif
        END DO
      END IF
      if (EASTERN_EDGE) then
        DO j=Jstr,Jend
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask_wet(Iend+1,j)
#   ifdef MRL_WCI
          ust2d(Iend+1,j)=ust2d(Iend+1,j)*umask_wet(Iend+1,j)
#   endif
        END DO
        DO j=JstrV,Jend
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask_wet(Iend+1,j)
#   ifdef MRL_WCI
          vst2d(Iend+1,j)=vst2d(Iend+1,j)*vmask_wet(Iend+1,j)
#   endif
        END DO
      END IF
#  endif
#  ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
        DO i=IstrU,Iend
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask_wet(i,Jstr-1)
#   ifdef MRL_WCI
          ust2d(i,Jstr-1)=ust2d(i,Jstr-1)*umask_wet(i,Jstr-1)
#   endif
        END DO
        DO i=IstrU,Iend
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask_wet(i,Jstr)
#   ifdef MRL_WCI
          vst2d(i,Jstr)=vst2d(i,Jstr)*vmask_wet(i,Jstr)
#   endif
        END DO
      END IF
      if (NORTHERN_EDGE) then
        DO i=Istr,Iend
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask_wet(i,Jend+1)
#   ifdef MRL_WCI
          ust2d(i,Jend+1)=ust2d(i,Jend+1)*umask_wet(i,Jend+1)
#   endif
        END DO
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask_wet(i,Jend+1)
#   ifdef MRL_WCI
          vst2d(i,Jend+1)=vst2d(i,Jend+1)*vmask_wet(i,Jend+1)
#   endif
        END DO
      END IF
#  endif
# endif
!
! zeta vill be recomputed via depth-integrated continuity equation
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,knew))
# endif
#ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(ubar(:,:,knew),'ubar step3d_fast','u'
     &    ,ondevice=.TRUE.)
       call check_tab2d(vbar(:,:,knew),'vbar step3d_fast','v'
     &    ,ondevice=.TRUE.)
#endif
         
# ifdef OBC_VOLCONS
      call obc_flux_tile (Istr,Iend,Jstr,Jend)
# endif

!      call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!      call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
!      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmu_nbq(START_2D_ARRAY,1))
!      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmv_nbq(START_2D_ARRAY,1))
!#  endif
! !
! !--------------------------------------------------------------------
! ! Compute fast-time averaged barotropic mass fluxes along physical
! ! boundaries.
! !--------------------------------------------------------------------
! !
! ! KERNEL_32  Dnew <= ( h, zeta )
! ! KERNEL_32  cff <= ( Dnew, ubar, on_u, vbar, om_v )
! ! KERNEL_32  DU_avg1 <= ( DU_avg1, cff, cff1 )
! ! KERNEL_32  DU_avg2 <= ( DU_avg2, cff, cff2 )
! ! KERNEL_32  DV_avg1 <= ( DV_avg1, cff, cff1 )
! ! KERNEL_32  DV_avg2 <= ( DV_avg2, cff, cff2 )

!$acc kernels default( present )
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Istr-1,j)=(h(Istr-1,j)+zeta(Istr-1,j,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(Istr-1,j,knew)
#  endif
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Iend+1,j)=(h(Iend+1,j)+zeta(Iend+1,j,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(Iend+1,j,knew)
#  endif
        enddo
      endif
# endif /* !EW_PERIODIC */

# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jstr-1)=(h(i,Jstr-1)+zeta(i,Jstr-1,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(i,Jstr-1,knew)
#  endif
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jend+1)=(h(i,Jend+1)+zeta(i,Jend+1,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(i,Jend+1,knew)
#  endif
        enddo
      endif
# endif /* !NS_PERIODIC */
!$acc end kernels

      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)

!$acc kernels default( present )      
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          cff=(Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
# ifdef MRL_WCI
     &                                              +ust2d(IstrU-1,j)
# endif
     &                                              )*on_u(IstrU-1,j)
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DU_avg2(IstrU-1,j)=DU_avg2(IstrU-1,j)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */      
        enddo
        do j=JstrV,Jend
          cff=(Dnew(Istr-1,j)+Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
# ifdef MRL_WCI
     &                                               +vst2d(Istr-1,j)
# endif
     &                                               )*om_v(Istr-1,j)
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DV_avg2(Istr-1,j)=DV_avg2(Istr-1,j)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */      
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          cff=(Dnew(Iend+1,j)+Dnew(Iend,j))*(ubar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                          +ust2d(Iend+1,j)
# endif
     &                                          )*on_u(Iend+1,j)
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DU_avg2(Iend+1,j)=DU_avg2(Iend+1,j)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */
        enddo
        do j=JstrV,Jend
          cff=(Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*(vbar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                              +vst2d(Iend+1,j)
# endif
     &                                              )*om_v(Iend+1,j)
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DV_avg2(Iend+1,j)=DV_avg2(Iend+1,j)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */     
        enddo
      endif
# endif /* !EW_PERIODIC */
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=IstrU,Iend
          cff=(Dnew(i,Jstr-1)+Dnew(i-1,Jstr-1))*(ubar(i,Jstr-1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jstr-1)
#  endif
     &                                              )*on_u(i,Jstr-1)
          DU_avg1(i,Jstr-1,nnew)=DU_avg1(i,Jstr-1,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DU_avg2(i,Jstr-1)=DU_avg2(i,Jstr-1)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */       
        enddo
        do i=IstrR,IendR
          cff=(Dnew(i,JstrV-1)+Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
#  ifdef MRL_WCI
     &                                              +vst2d(i,JstrV-1)
#  endif
     &                                              )*om_v(i,JstrV-1)
          DV_avg1(i,JstrV-1,nnew)=DV_avg1(i,JstrV-1,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DV_avg2(i,JstrV-1)=DV_avg2(i,JstrV-1)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */     
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=IstrU,Iend
          cff=(Dnew(i,Jend+1)+Dnew(i-1,Jend+1))*(ubar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jend+1)
#  endif
     &                                              )*on_u(i,Jend+1)
          DU_avg1(i,Jend+1,nnew)=DU_avg1(i,Jend+1,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DU_avg2(i,Jend+1)=DU_avg2(i,Jend+1)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */       
        enddo
        do i=IstrR,IendR
          cff=(Dnew(i,Jend+1)+Dnew(i,Jend))*(vbar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                          +vst2d(i,Jend+1)
#  endif
     &                                          )*om_v(i,Jend+1)
          DV_avg1(i,Jend+1,nnew)=DV_avg1(i,Jend+1,nnew)+cff1*cff
# ifndef M3FAST_AVG_CLASSIC
          DV_avg2(i,Jend+1)=DV_avg2(i,Jend+1)+cff2*cff
# endif /*   M3FAST_AVG_CLASSIC */     
        enddo
      endif
# endif /* !NS_PERIODIC */
!$acc end kernels
      if (iif.eq.nfast) then
!$acc update host( DU_avg1, DU_avg2, DV_avg1, DV_avg2 )   !! iif=last
      endif
! !
! !====================================================================
! !
! !  Adjust ZETA using depth-integrated continuity equation 
! !  and update grid
! !
! !  Once rhobar_nbq and depth-averaged momentum is updated, surface 
! !  anomalies can be adjusted to satisfy the low-frequency mass 
! !  conservation equation. However, this adjustment do not satisfy 
! !  conservation at machine precision and a final correction is 
! !  needed.
! !
! !  This recomputation of zeta(m+1) using div(ubar) is only done for 
! !  PRECISE ! !  option. In PERF option, zeta(m+1) is computed only 
! ! once based on the surface characteristic relation. But the final 
! ! numerical correction is applied in all cases.
! !
! !====================================================================
! !
# ifdef NBQ_HZCORRECT_ZETA
!-----------------------------------------------------------------------
!  Update zeta(m+1)
!-----------------------------------------------------------------------
!
! ! KERNEL_33  Dnew <= ( zeta, h, pm, pn, Dnew, ubar, on_u, vbar, om_v )
! ! KERNEL_33  zeta <= ( zeta, rmask )

!     if (LAST_FAST_STEP) then
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
!     endif
!#  endif
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
!# ifdef NBQ
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
!
!  Diagnostic evaluation from zeta(m+1)
!
# ifdef OPENACC
#  undef exchange_r2d_tile 
#  undef exchange_u2d_tile 
#  undef exchange_v2d_tile 
#  undef exchange_r3d_tile 
#  undef exchange_u3d_tile 
#  undef exchange_v3d_tile 
#  undef exchange_w3d_tile 
# endif
!# ifdef M3FAST_COUPLING2D  
!                          CAUTION: necessary also when !defined (TBT)
        call set_depth_tile(Istr,Iend,Jstr,Jend)
!#   endif        
#  endif  /* NBQ_HZ_PROGNOSTIC */
!
! Compute derived grid parameters if fast update
!
!#  ifdef NBQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ATTENTION FRANCIS: this call should not be needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  ifndef NBQ_GRID_SLOW
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u  , Hzw_nbq_inv_v)
#  endif
!#  endif

#  ifdef NBQ_GRID_SLOW
      endif !<-- LAST_FAST_STEP
#  endif
!# endif /* NBQ */
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
! # ifdef RVTK_DEBUG
! C$OMP BARRIER
! C$OMP MASTER
!        call Check_tab2d(zeta(:,:,knew),'zeta step3d_fast','r')
!  # endif      
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
!       if (LAST_FAST_STEP) then 
! # ifdef RVTK_DEBUG
! C$OMP BARRIER
! C$OMP MASTER
!          call Check_tab2d(zeta(:,:,knew),'zeta step3d_fast last','rint')
! # endif      
!       endif
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !--------------------------------------------------------------------
! !  Apply conservation requirements for nesting
! !--------------------------------------------------------------------
! !
# ifdef AGRIF
      if (.NOT.Agrif_Root()) THEN
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,iif)=zeta(i,j,knew)
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
            DU_avg3(i,j,iif) = 0.5*(h(i  ,j)+zeta(i  ,j,knew)+
     &                              h(i-1,j)+zeta(i-1,j,knew)) 
     &                               *on_u(i,j)*ubar(i,j,knew)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
            DV_avg3(i,j,iif) = 0.5*(h(i,j  )+zeta(i,j  ,knew)+
     &                              h(i,j-1)+zeta(i,j-1,knew)) 
     &                               *om_v(i,j)*vbar(i,j,knew)
          enddo
        enddo
 
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          Zt_avg3(START_2D_ARRAY,iif))
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg3(START_2D_ARRAY,iif))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg3(START_2D_ARRAY,iif))
#  endif

#  ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(Zt_avg3(:,:,iif),'Zt_avg3 step3d_fast','r'
      &    ,ondevice=.TRUE.)
        call check_tab2d(DU_avg3(:,:,iif),'DU_avg3 step3d_fast','u'
      &    ,ondevice=.TRUE.)
        call check_tab2d(DV_avg3(:,:,iif),'DV_avg3 step3d_fast','v'
      &    ,ondevice=.TRUE.)
C$OMP END MASTER
#  endif   
      endif

#  ifdef AGRIF_CONSERV_VOL
      if (iif.eq.nfast) then
        if (agrif_root()) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              DU_avg1(i,j,5) = dt * DU_avg2(i,j)
              DV_avg1(i,j,5) = dt * DV_avg2(i,j)
            enddo
          enddo
        else
          do j=JstrR,JendR
            do i=IstrR,IendR
              DU_avg1(i,j,5) = dt * DU_avg2(i,j)
              DV_avg1(i,j,5) = dt * DV_avg2(i,j)
              DU_avg1(i,j,4) = DU_avg1(i,j,4) + DU_avg1(i,j,5)
              DV_avg1(i,j,4) = DV_avg1(i,j,4) + DV_avg1(i,j,5)
            enddo
          enddo       
        endif
      endif
#  endif
# endif /* AGRIF */
! !
! !--------------------------------------------------------------------
! !  Copy density for extrapolation: 
! !--------------------------------------------------------------------
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
! !--------------------------------------------------------------------
! !  
! !--------------------------------------------------------------------
! !
# ifdef NBQ_FREESLIP
       if (iif==1) then
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              qdmw0_nbq(i,j)=qdmw_nbq(i,j,0)
            enddo
           enddo
        elseif (iif.lt.ndtfast) then
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              qdmw0_nbq(i,j)=qdmw0_nbq(i,j)+qdmw_nbq(i,j,0)
            enddo
           enddo
        else
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              qdmw0_nbq(i,j)=(qdmw0_nbq(i,j)+qdmw_nbq(i,j,0))
     &               /float(ndtfast)
            enddo
           enddo
           
        endif
# endif
!
!-----------------------------------------------------------------------
!  Correct Hz(m+1) for internal mode
!  by inverting internal continuity equation
!-----------------------------------------------------------------------
!
# ifdef NBQ_HZCORRECT
      if (LAST_FAST_STEP) then
       do k=-N_sl+1,N
        do j=JstrV-2,Jend+1
         do i=IstrU-2,Iend+1
              Hz_correct(i,j,k)=Hz(i,j,k)
          enddo
        enddo
       enddo
!$acc update device( Hz_correct )      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !   ATTENTION  FRANCIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg2(START_2D_ARRAY))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg2(START_2D_ARRAY))
#  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do j=Jstr,Jend
          do i=Istr,Iend
            dum_s=0.
         !   dum2_s=0.
            do k=1,N
              dum_s=dum_s+Hz(i,j,k)-Hz_bak(i,j,k)
        !      dum2_s=dum2_s+Hz(i,j,k)
            enddo
            do k=1,N
!              Hz(i,j,k)=Hz(i,j,k)*(1.
!     &                -(dum_s
!     &             + (DU_avg2(i+1,j)-DU_avg2(i,j)
!     &               +DV_avg2(i,j+1)-DV_avg2(i,j)
!     &               )*pm(i,j)*pn(i,j)
!     &               *dt)
!     &               /dum2_s )
!!     &               /(z_w(i,j,N)-z_w(i,j,0))
!!     &               *(z_w(i,j,k)-z_w(i,j,k-1))
!#  ifdef MASKING
!     &               *rmask(i,j)
!#  endif

!#  ifdef NBQ_HZCORR_DEBUG
!              Hz_corr(i,j,k)=Hz(i,j,k)-Hz_correct(i,j,k)
!#  endif

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


!# ifdef NBQ_MASS
!              Hzr(i,j,k)=(Hz(i,j,k)-rho_nbq(i,j,k)) ! adjust Zw,Zr?
!     &                        /(1.+rho(i,j,k)/rho0)
!# endif
!              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
!              z_r(i,j,k)=0.5*(z_w(i,j,k)+z_w(i,j,k-1))

            enddo
          enddo
        enddo

!# undef Hz_correct

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#   ifndef M3FAST_SEDLAYERS
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,1))
#   else
        call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,-N_sl+1))
#   endif
!# ifdef NBQ_MASS
!         call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                           Hzr(START_2D_ARRAY,1))
!# endif
!        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                          z_w(START_2D_ARRAY,0))
!        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                          z_r(START_2D_ARRAY,1))
#  endif 
!$acc update device( Hz )      

      endif ! iif.eq.nfast

# endif  /* NBQ_HZCORRECT */
!
# undef Hz_correct
!-----------------------------------------------------------------------
!  Mass/volume conservation
!-----------------------------------------------------------------------
!
# if defined NBQ_MASS && defined MOVING_BATHY
!     if (LAST_FAST_STEP) then
      masstot=0.
      do j=Jstr,Jend
        do i=Istr,Iend
          do k=1,N
          masstot=masstot+Hz(i,j,k)
          enddo
        enddo
      enddo
#  ifdef MPI
      call MPI_ALLGATHER(masstot,1,MPI_DOUBLE_PRECISION,
     &                allmasstot,1,MPI_DOUBLE_PRECISION,
     &                          MPI_COMM_WORLD,ierr)
      masstot=QuadZero
      do i=1,NNODES
        masstot=masstot+allmasstot(1,i)
      enddo
#  endif
      if (mynode==0) then
       open(unit=10,file='diag_mass.dat',access='append')
       write(10,*) masstot
       close(10)
      endif
!     endif
# endif

# ifdef M3FAST_TEST
!  OUTPUT  FRANCIS ICI
      if (mynode==0) then
       i=istr+10
       j=jstr+10
       write(6,*) '==>',iic,iif
       write(6,*) '   zeta',zeta(i,j,knew),ubar(i,j,knew),vbar(i,j,knew)
       write(6,*) '   avg1',DU_avg1(i,j,nnew),DV_avg1(i,j,nnew)
       write(6,*) '   avg2',DU_avg2(i,j),DV_avg2(i,j)
      endif
# endif
     
! !
! !--------------------------------------------------------------------
! !  TEST FOR CFL VIOLATION. IF SO, PRINT AND STOP
! !--------------------------------------------------------------------
! !
      VMAXL=100.
      VMAX=0.
      do j=Jstr,Jend
        do i=Istr,Iend
          cff1=ubar(i,j,knew)
          cff2=vbar(i,j,knew)
          cff=max(abs(cff1),abs(cff2))
          IF (cff.GE.VMAX .or. cff1.ne.cff1 .or. cff2.ne.cff2) THEN
            IF (cff.GE.VMAX .and. cff1.eq.cff1 .and. cff2.eq.cff2) THEN 
              VMAX=cff
            ELSE
              VMAX=666.
            ENDIF
# ifdef MPI
            imax=i+iminmpi-1
            jmax=j+jminmpi-1
# else
            imax=i
            jmax=j
# endif
          ENDIF
        enddo
      enddo
      
      IF (VMAX.GT.VMAXL) THEN
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =  STEP3D_FAST:   ABNORMAL JOB END    = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '

# ifdef MPI
        write(stdout,'(A,I4)') '      mynode =',mynode
# endif

        if (VMAX.eq.666.) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)') 
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')  
     &                                       '  IMAX JMAX   =',imax,jmax
        write(stdout,'(A,2I6/)')
     &                                       '  IINT IFAST  =',iic,iif

!     GC quik fix Comment this because it blocks the correct mpi_abort
!     croco do not hang out  
!!      call wrt_his            ! get output during blow-up

        may_day_flag=1
# ifdef MPI
        call mpi_abort (MPI_COMM_WORLD, err)
# else
        stop                            !-->  EXIT
# endif
      ENDIF

# undef zwrk 
# undef rzeta  
# undef rzeta2  
# undef rzetaSA 

      return
      end
#else
      subroutine step3d_fast_empty
      end
#endif /* M3FAST */
