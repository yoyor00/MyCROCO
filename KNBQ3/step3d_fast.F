! $Id: step3d_fast.F 1615 2014-12-17 13:27:07Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
!  and Rutgers University (Arango et al) are under MIT/X style license
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

#ifdef K3FAST
! ! 
! !*********************************************************************
! !*********************************************************************
! !                      step3d_fast (interface)
! !*********************************************************************
! !*********************************************************************
! !
      subroutine step3d_fast (tile)
! !
! !*********************************************************************
! ! Variables & parameters (interface)
! !*********************************************************************
! !
      implicit none
      integer tile, trd
# include "param.h"
# include "private_scratch.h"
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      trd=0
C$    trd=omp_get_thread_num()
      call step3d_fast_tile ( Istr,Iend,Jstr,Jend
     &                   , A2d(1, 1,trd), A2d(1, 2,trd), A2d(1, 3,trd)
     &                   , A2d(1, 4,trd), A2d(1, 5,trd), A2d(1, 6,trd)
     &                   , A2d(1, 7,trd), A2d(1, 8,trd)
     &                   , A2d(1, 9,trd), A2d(1,10,trd)
     &                   , A2d(1,11,trd), A2d(1,12,trd)
     &                   , A2d(1,13,trd)
# ifdef K3FAST_C3D_UVSF
     &                   , A2d(1,14,trd), A2d(1,15,trd)
     &                   , A2d(1,16,trd), A2d(1,17,trd)
#  ifdef K3FAST_COUPLING3D
     &                   , A2d(1,18,trd), A2d(1,19,trd)
#  endif
# endif
# if defined K3FAST_UV
     &                   , A3dHz(1, 1,trd), A3dHZ(1, 2,trd)
     &                   , A3dHz(1, 3,trd), A3dHZ(1, 4,trd)   
# endif
     &                   )

      return
      end
! ! 
! !*********************************************************************
! !*********************************************************************
! !           Step3d_fast: tile routine
! !*********************************************************************
! !*********************************************************************
! !
       subroutine step3d_fast_tile (Istr,Iend,Jstr,Jend
     &                          ,Dnew,rubar,rvbar
     &                          ,Drhs, UFx,UFe
     &                          ,VFx,VFe
     &                          ,urhs,vrhs 
     &                          ,DUon,DVom ,work
# ifdef K3FAST_C3D_UVSF
     &                          ,ru_ext_nbq_sum, rv_ext_nbq_sum
     &                          ,ru_ext_nbq_old, rv_ext_nbq_old
#  ifdef K3FAST_COUPLING3D
     &                          ,ru_int2d_nbq, rv_int2d_nbq
#  endif
# endif
# if defined K3FAST_UV 
     &                          ,Hzw_nbq_inv,  Hzr_nbq_inv 
     &                          ,Hzw_nbq_inv_u,  Hzw_nbq_inv_v  
# endif
     &                          )
     
! !
! !
! !*********************************************************************
! !*********************************************************************
! !                    SOLVE FAST MODE 3D EQUATIONS
! !*********************************************************************
! !*********************************************************************
! !
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
!     density rho_nbq by time integration of continuity.
!
!  Treatment of the free-surface:
!  A first guest of zeta is derived from the surface kinematic relation
!  After solving the 3D momentum equations, zeta is updated to conserve
!  mass for each water column.
!
!  Forward-backward time-stepping:
!  W-momentum equation is solved with explicit or implicit methods:
!  - Explicit scheme: w-momentum is updated right after (and the same 
!                     way as) u- and v-momentum. Only rho_nbq is
!                     computed with backward scheme.
!  - Implicit scheme: rho_nbq(m+1) is substituted in w-momentum Eq.
!                     and a tridiagonal Gauss Elimination is carried
!                     out for qdmw_nbq(m).
!  For all components, a Forward-Backward scheme is implemented:
!  - Explicit scheme: Forward:  zeta, qdmu_nbq, qdmw_nbq. 
!                     Backward:       rho_nbq.
!  - Implicit scheme: Forward:  zeta, qdmu_nbq. 
!                     Backward:       qdmw_nbq, rho_nbq.
!***********************************************************************
! !  Overview:
! !  PART A: Declare variables, parameters...
! !  PART B: 2D-mode processing.
! !          B.0: Prepare integration
! !          B.1: (AGRIF) prepare 2D mode integration    
! !          B.2: Compute zeta(m+1) first guess.
! !          B.3: Update vertical grid.
! !          B.4: AM4 Backward step.
! !          B.5: Pre-step 2D.
! !  PART C: Solve fast-mode 3D compressible equations.
! !          C.0: Initialization & backups for 3D eq.
! !          C.1: Compute NT Coriolis pseudo-force.
! !          C.2: Compute "Pressure - Viscosity" component (theta).
! !          C.3: Compute fast-mode U & V momentum equations.
! !          C.4: Compute fast-mode W-momentum equation.
! !          C.5: Compute divergence operator.
! !          C.6: Fast-mode conservation of mass.
! !  PART D: Fast-mode post-processing.
! !          D.1: Compute time average & adjust depth-averaged variables
! !          D.2: Adjust zeta.
! !          D.3: (AGRIF) Apply conservation requirements for nesting
! !          D.4: Adjust Hz(m+1) for internal mode.
! !  PART E: Test for CFL violation.
! !*********************************************************************
! !      
! !*********************************************************************
! !*********************************************************************
! ! PART A: declare variables, parameters...
! !*********************************************************************
! !*********************************************************************
! !
      use debug
      implicit none
! ! 
! !*********************************************************************
! ! PART A.1: Parameters
! !*********************************************************************
! ! 
# include "param.h"
! !      
! !*********************************************************************
! ! PART A.2: Local variables
! !*********************************************************************
! ! 
      integer Istr,Iend,Jstr,Jend, i,j,k, kbak,kold,
     &        imin,imax,jmin,jmax,
     &        kp1
# ifdef MPI
     &       ,err
# endif
# ifdef PSOURCE
     &       ,is
# endif
      real    mybeta,myalpha,myepsilon,mygamma,
     &        VMAX,VMAXL, cff,cff0,cff1,cff2,cff3,cff4, 
     &        DUnew,DVnew, dum_s, dum2_s
      real  sum_nbq

# ifdef NBQ_THETAIMP
      real,parameter :: thetaimp_nbq = 0.6 
# else
      real,parameter :: thetaimp_nbq = 1.
# endif

# if defined K3FAST_PG2
      real,parameter :: gammau=0.45   ! NBQ pressure gradient scheme
      real,parameter :: gammau_2=(1./3.)*(1.-gammau)
# endif
      real
     &                Dnew(PRIVATE_2D_SCRATCH_ARRAY),
     &               rubar(PRIVATE_2D_SCRATCH_ARRAY),
     &               rvbar(PRIVATE_2D_SCRATCH_ARRAY),
     &                Drhs(PRIVATE_2D_SCRATCH_ARRAY),
     &                 UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &                 UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &                 VFx(PRIVATE_2D_SCRATCH_ARRAY),
     &                 VFe(PRIVATE_2D_SCRATCH_ARRAY),
     &                urhs(PRIVATE_2D_SCRATCH_ARRAY),  
     &                vrhs(PRIVATE_2D_SCRATCH_ARRAY),
     &                DUon(PRIVATE_2D_SCRATCH_ARRAY),
     &                DVom(PRIVATE_2D_SCRATCH_ARRAY)
      real
     &                 work(PRIVATE_2D_SCRATCH_ARRAY)
# ifdef K3FAST_C3D_UVSF
#  ifdef K3FAST_COUPLING2D
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
#  ifdef K3FAST_COUPLING3D
      real
     &       ru_int2d_nbq  (PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_int2d_nbq  (PRIVATE_2D_SCRATCH_ARRAY)
#  endif
# endif
# ifdef KNHINT_CORR
      real kzw
# endif
        real
     &   Hzw_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl  :N),
     &   Hzr_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl+1:N),
     &   Hzw_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &   Hzw_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N)
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
!$acc declare create( FX, FY )     
!$acc declare create( FC, CF )     
!$acc declare create( dZdxq_u, dZdyq_v )
# ifdef UV_COR_NT
      real
     &               ntcoru(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &               ntcorv(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N),
     &               ntcorw(PRIVATE_2D_SCRATCH_ARRAY,-N_sl:N)
# endif 
# ifdef MRL_WCI
      real vstu,ustv,dudx,dvdx,dude,dvde
# endif
# if defined INTERNAL || defined BODYTIDE
      real,parameter :: U0 = 0.02
      real,parameter :: omega = 2.*pi/(12.4*3600.)
# endif
# ifdef K3FAST_SACOUS
      real dist_d  
      integer jvar1,jvar2
# endif
# ifdef WET_DRY
      real cff1_WD,cff2_WD
# endif
# ifdef NBQ_DIAGMASS
      real*QUAD masstot,masstot2
      real*QUAD masstot0,masstot02
      common /diag_mass/ masstot0,masstot02
#  ifdef MPI
      real*QUAD allmasstot(1,NNODES)
      real*QUAD allmasstot2(1,NNODES)
      integer ierr
#  endif
# endif      
! !      
! !*********************************************************************
! ! PART A.3: Commons
! !*********************************************************************
! ! 
# ifdef OPENACC
#  include "set_obc_definitions.h"
#endif
# include "nbq.h"
# ifdef AGRIF
#  include "zoom.h"
        integer irhox, irhoy, irhot
# endif
# include "set_obc_definitions.h" 
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
# ifdef OPENACC
#  include "boundary.h"
#endif
# include "sources.h"
! !      
! !*********************************************************************
! ! PART A.4: Preprocess MASKING
! !*********************************************************************
! ! 
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
! !      
! !*********************************************************************
! ! PART A.5: Loop Bounds
! !*********************************************************************
! ! 
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
! !      
! !*********************************************************************
! ! PART A.6: OpenACC
! !*********************************************************************
! ! 
! !      
! !*********************************************************************
! ! PART A.7: KNHINT: 3 time-steps
! !*********************************************************************
! ! 
# ifndef KNHINT_3M
#  define FIRST_FAST_STEP_3M FIRST_FAST_STEP
# else
#  define FIRST_FAST_STEP_3M iif.eq.nsdtnbq 
# endif
! !
! !*********************************************************************
! ! PART A.8: Set several approximations to speed-up computation
! !*********************************************************************
! !
! !---------------------------------------------------------------------
! ! d./ds terms update frequency
! !---------------------------------------------------------------------
! !
#  ifdef NBQ_GRID_SLOW
#   define NSTEP_DS FIRST_FAST_STEP_3M 
#  else
#   define NSTEP_DS .true. 
#  endif
! !
! !---------------------------------------------------------------------
! ! Grid update frequency (Hz ...)
! !---------------------------------------------------------------------
! !
#  ifdef NBQ_GRID_SLOW
#   define NSTEP_GRID .or. iif.eq.nfast
#  else
#   define NSTEP_GRID .or. mod(iif,1).eq.0
#  endif
! !
! !*********************************************************************
! !*********************************************************************
! ! PART B:        EXTERNAL (2D) FAST MODE PROCESSING
! !
! ! Compute external 2D rhs terms rubar,rvbar of 2D fast-mode equations
! ! using Generalized Forward-Backward AB3-AM4 algorithm, and update 
! ! internal and external forcing terms for 3D fast mode equations: 
! ! ru_int, ru_ext~(rufrc+rubar).
! !
! !  1- AB3 forward step for D,ubar,vbar
! !  2- Advance zeta(m+1)
! !  3- AM4 backward step for rubar,rvbar 
! !  4- update rufrc,rvfrc
! !  5- update ru_int,rv_int
! !  6- make some backups
! !
! ! Reference for Generalized FB scheme:
! ! ----------------------------------------------
! ! Shchepetkin, A.F., and J.C. McWilliams, 2009: Computational kernel 
! ! algorithms for fine-scale, multiprocess, longtime oceanic sim. 
! ! Pp. 119–182 in Handbook of Numerical Analysis: Computational Met.
! ! for the Atmosphere and Oceans. R.M. Teman and J.J. Tribbia, eds, 
! ! Elsevier Science.
! !
! !*********************************************************************
! !*********************************************************************
! !
!         call system_clock(cput1, count_rate, count_max)
! !
! !*********************************************************************
! ! PART B.0: Prepare 2D Mode integration
! !*********************************************************************
! !
#include "k3fast_2Dmode_prep.h"
! !
! !*********************************************************************
! ! PART B.1: (AGRIF) : Load RHS values into additional
! !                     AGRIF shared array for nesting
! !*********************************************************************
! !
# ifdef AGRIF
#  include "k3fast_AGRIF0.h"
# endif /* AGRIF */
! !
! !*********************************************************************
! ! PART B.2: Compute zeta(m+1) -- first guest -- and update grid
! !*********************************************************************
! !
! !--------------------------------------------------------------------
! ! Recover Hz at first fast step (if final Hz correction needed)
! !--------------------------------------------------------------------
! !
# ifdef NBQ_HZCORRECT
!$acc kernels if(compute_on_device) default(present)
      if (iic.gt.ntstart .and. FIRST_FAST_STEP) then
        do k=-N_sl+1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              Hz(i,j,k) = Hz_correct(i,j,k)
            enddo
          enddo
        enddo
      endif
!$acc end kernels
# endif
! !
! !--------------------------------------------------------------------
! ! Get derived z grid variables at first fast step
! !--------------------------------------------------------------------
! !
# if ! defined K3FAST_COUPLING2D || defined K3FAST_UV || defined K3FAST_W
       if ( FIRST_FAST_STEP) then
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u  , Hzw_nbq_inv_v)
       endif
# endif
! !
! !--------------------------------------------------------------------
! ! zeta(m+1) is derived here from the surface vertical velocity 
! ! (surface characteristic relation) instead of depth-averaged 
! ! conservation of mass. This satisfies dynamical coupling with 
! ! the surface layer. The surface kinematic relation is:
! ! zeta(m+1) = zeta(m) + dtfast * [ wsurf_nbq(m)
! !                                 - usurf_nbq(m)*dzeta/dx(m+0.5)
! !                                 - vsurf_nbq(m)*dzeta/dy(m+0.5) ]
! !
! ! with zeta(m+0.5)=(1.5+beta)*zeta(m)-(0.5+2*beta)*zeta(m-1)
! !                                            +beta*zeta(m-2) --> AB3
! !--------------------------------------------------------------------
! !
#include "k3fast_zeta_update.h"
! !
! !*********************************************************************
! !  PART B.3: Update vertical grid
! !
! !  As soon as zeta(m+1) is known, the grid can be updated at fast
! !  step m+1. In PERF option, the grid is updated at lower frequency 
! !  and only during re-evaluation of zeta(m+1).
! !
! !  Caution: Hz_bak must be assigned only once in set_depth. The
! !  following code must thus be consistant with set_depth routine.
! !*********************************************************************
! !
# if !defined NBQ_GRID_SLOW
! !--------------------------------------------------------------------
! ! Update moving bathymetry analytically
! !--------------------------------------------------------------------
#  ifdef MVB
        call ana_mvb_tile(Istr,Iend,Jstr,Jend)
#  endif /* MVB */
        if ( (FIRST_TIME_STEP .and. FIRST_FAST_STEP) NSTEP_GRID ) then
! !--------------------------------------------------------------------
! ! Update main grid parameters:
! ! z_w, z_r, Hz, Hz_bak <= h, zeta, hinv, Hz 
! !--------------------------------------------------------------------
        call set_depth_tile(Istr,Iend,Jstr,Jend)
! !--------------------------------------------------------------------
! ! Update derived grid variables
! !--------------------------------------------------------------------
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_nbq_inv,   Hzr_nbq_inv,
     &                     Hzw_nbq_inv_u, Hzw_nbq_inv_v)
     
      endif
# endif /* !NBQ_GRID_SLOW */
! ! 
! !*********************************************************************
! ! PART B.4: AM4 backward step: compute depth-averaged RHS of  
! ! --- -------- ----            2D fast-modemomentum equations
! ! Compute here:
! ! - Surface pressure gradient
! ! - Horizontal advection terms for momentum equations (2D only)
! ! - Compute Coriolis (2D and 3D) term and advective curvilinear metric
! !    terms (2D only).
! ! - Linear and/or quadratic bottom stress.
! ! - Add 2D vortex-force terms combined with advection terms
! !*********************************************************************
! !
# include "k3fast_2Dmode_rhs.h"
! ! 
! !*********************************************************************
! ! PART B.5: Pre_step 2D: where couling with 2D mode is dealt with.
! !*********************************************************************
! !
# include "k3fast_pre_step2d.h"
! !
! !********************************************************************
! !********************************************************************
! ! PART C:         Solve fast-mode 3D compressible equations.
! !********************************************************************
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
! !********************************************************************
! ! PART C.0: Initializations and backups for 3D equations
! !********************************************************************
! !
! !--------------------------------------------------------------------
! !  Fast bottom friction
! !
! !  Set bottom stress using logarithmic or linear 
! !  and/or quadratic formulation. 
! !--------------------------------------------------------------------
! !
# if defined BSTRESS_FAST && !defined BBL
#  include "k3fast_fbf.h"
# endif
! !
! !--------------------------------------------------------------------
! !  Initializations of ru_nbq_avg2 (as qdmu_nbq) at first fast step
! !  It is used at last fast step for computation of new ru_nbq_avg2
! !--------------------------------------------------------------------
! !
# include "k3fast_init.h"
! !
! !--------------------------------------------------------------------
! !  Store boundary values of nbq variables at previous 
! !  time-step for use in radiation boundary conditions
! !--------------------------------------------------------------------
! !
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
      call nbq_bry_store_tile (Istr,Iend,Jstr,Jend)
# endif
! !
! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ! KNHINT: Fast 3D-Mode conditional integration (BEGIN)
! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! !
# ifdef KNHINT_3M
      if (mod(iif,nsdtnbq).eq.0) then  ! HEREFA
       dtfast=dtfast*float(nsdtnbq)
       dtgrid_nbq=dtgrid_nbq*float(nsdtnbq)
       ndtfast=ndtfast/nsdtnbq
# endif
! !
! !********************************************************************
! ! PART C.1: Compute Non-traditional Coriolis pseudo-force
! !
! !   ntcoru = -e W cos(a) * Hz      [m2/s2]
! !   ntcorv = +e W sin(a) * Hz
! !   ntcorw = +e (U cos(a) - V sin(a) ) * Hz
! ! 
! ! with e = 2 Omega cos(Phi)
! !      a = angle between North and meridional grid axis
! !
! !********************************************************************
!
# ifdef UV_COR_NT
#  include "k3fast_ntcoriolis.h"
# endif
! !
! !********************************************************************
! ! PART C.2: Compute "Pressure - Viscosity" component (theta)
! !********************************************************************
! ! 
# ifdef K3FAST_RHO 
!      if (FIRST_TIME_STEP.and.FIRST_FAST_STEP) then
! !
!$acc kernels default( present )    
      if (FIRST_TIME_STEP.and.FIRST_FAST_STEP_3M) then 
       do k=-N_sl+1,N
         do j=JstrV-2,Jend+1
           do i=IstrU-2,Iend+1
             thetadiv_nbq(i,j,k)=(
     &                     +rho_nbq(i,j,k) 
#  ifndef K3FAST_CSVISC2K
     &                           *soundspeed2_nbq
#  else
     &                           *soundspeed2_nbq(i,j,k)
#  endif 
     &                          )*Hzr_nbq_inv(i,j,k) 
           enddo
         enddo
       enddo   
      else
       do k=-N_sl+1,N
         do j=JstrV-2,Jend+1
           do i=IstrU-2,Iend+1
             thetadiv_nbq(i,j,k)=(
#  ifndef K3FAST_CSVISC2K 
     &                  -visc2_nbq       *thetadiv_nbq   (i,j,k)
     &                  +rho_nbq  (i,j,k)*soundspeed2_nbq
#  else
     &                  -visc2_nbq(i,j,k)*thetadiv_nbq   (i,j,k)
     &                  +rho_nbq  (i,j,k)*soundspeed2_nbq(i,j,k)
#  endif 
     &                           )*Hzr_nbq_inv(i,j,k) 
           enddo
         enddo
       enddo     
      endif
!$acc end kernels  
# endif
! !
! !********************************************************************
! ! PART C.3:  Compute fast-mode U & V momentum equations
! !
! !  dqdm/dt=dtfast*(Compressible pressure force + second viscosity
! !                  + gravity + non-NBQ forces + nudging)
! !
! !  If explicit scheme: all (x,y,z) mom components are done here
! !  If implicit scheme: z-mom component is done after
! !********************************************************************
! !
# include "k3fast_qdmuv_update.h"
! !
! !--------------------------------------------------------------------
! !      Compute horizontal part of divergence operator
! !              if implicit scheme
! !  (if so, from now on, thetadiv_nbq array is used
! !               for divergence, not theta)
! !--------------------------------------------------------------------
! !
# ifdef NBQ_IMP
#   include "k3fast_divh.h"
# endif
! !
# ifdef K3FAST_W
! !
! !********************************************************************
! ! PART C.4: Compute fast-mode W-Momentum equation: qdmw_nbq 
! !        - Explicit integration: compute qdmw_nh
! !        - Implicit integration: compute RHS of tri-diag system 
! !********************************************************************
! !
#  include "k3fast_qdmw_update.h"
! !
! !--------------------------------------------------------------------
! !   Implicit integration of W-momentum
! !--------------------------------------------------------------------
! !
! !   Initialisation FX and FYX out here
! !
! ! Compute and solve tridiag (implicit) system for W
! !
#  ifdef NBQ_IMP
#   include "k3fast_qdmw_implicit.h"
#  endif  
! !
! !********************************
! ! KNHINT numerical mode control: 
! ! remove potential surface component
! !********************************
! !
#  ifdef KNHINT_CORR 
!       kzw=2.*pi/2./sqrt(om_r(1,1)**2+on_r(1,1)**2)
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
              do k=N,Max(1,N-2*alphaNw_nbq),-1
!              cff=(alphaw_nbq-1.)
!    &               *exp(-(z_w(i,j,k)            -z_w(i,j,N))**2
!    &                    /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2)
!               cff2= cff*
!     &              qdmw_nbq(i,j,N)
!     &         *(Hz(i,j,k)+Hz(i,j,k+1))/Hz(i,j,N)
!! Linear 1:
!!!   &         *(z_w(i,j,k)+H(i,j))/(z_w(i,j,N)+H(i,j))
!! Linear 2:
!     &         *(z_w(i,j,k)-z_w(i,j,N-alphaNw_nbq))
!     &         /(z_w(i,j,N)-z_w(i,j,N-alphaNw_nbq))
!! Short waves:
!!    &         *sinh(kzw*(H(i,j)+z_w(i,j,k)))
!!    &         /sinh(kzw*H(i,j))
               qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)
     &               *(1.-alphaw_nbq
     &               *exp(-(z_w(i,j,k)            -z_w(i,j,N))**2
     &                    /(z_w(i,j,N-alphaNw_nbq)-z_w(i,j,N))**2))
!
             enddo
!            k=N
!            cff2=qdmw_nbq(i,j,N)*(alphaw_nbq-1.)
!            qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)*0.9
           enddo
         enddo
#   endif /* KNHINT_CORR  */
! !
! !--------------------------------------------------------------------
! !  W-momentum open boundary conditions
! !--------------------------------------------------------------------
! !
#  ifdef OBC_NBQ
      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
#  endif
! !
      if (LAST_FAST_STEP) then
!$acc update host( qdmw_nbq )   !! iif=last
      endif 
! !
! !--------------------------------------------------------------------
! !  Exchange qdmw (if necessary)
! !--------------------------------------------------------------------
! ! 
#  ifdef UV_COR_NT
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_w3d_sedlayb_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmw_nbq(START_2D_ARRAY,-N_sl))   
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmw_nbq(START_2D_ARRAY,0))   
#   endif
#  endif
! !
! !********************************************************************
! ! PART C.5: compute divergence operator.
! !********************************************************************
! !
! !--------------------------------------------------------------------
! ! Compute horizontal components of divergence operator
! !              if explicit scheme
! !  (if so, from now on, thetadiv_nbq array is used
! !               for divergence, not theta)
! !--------------------------------------------------------------------
! !
#  ifndef NBQ_IMP
#   include "k3fast_divh.h"
#  endif
! !
! !--------------------------------------------------------------------
! ! Vertical part of divergence
! !      div(m+1) = div(m+1) + divZ(qdmw_nbq(m+1))
! !--------------------------------------------------------------------
! !		
! !.......Computes fluxes:  
! !
#  include "k3fast_divv.h"
! !
! !--------------------------------------------------------------------
! !         Exchange divergence
! !--------------------------------------------------------------------
! !	
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   ifndef K3FAST_SEDLAYERS
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv_nbq(START_2D_ARRAY,1))
#     if defined NBQ_HZ_PROGNOSTIC  
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv2_nbq(START_2D_ARRAY,1))
#     endif
#   else
      call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv_nbq(START_2D_ARRAY,-N_sl+1))
#     if defined NBQ_HZ_PROGNOSTIC  
       call exchange_r3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv2_nbq(START_2D_ARRAY,-N_sl+1))
#     endif
#   endif
#  endif
# endif  /* K3FAST_W */
! !
! !********************************************************************
! ! PART C.6: Fast-mode conservation of mass
! !********************************************************************
! !
# include "k3fast_mass_update.h"
! !
! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ! KNHINT: Fast 3D-Mode conditional integration (END)
! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! !
# ifdef KNHINT_3M
        dtfast=dtfast/float(nsdtnbq)
        dtgrid_nbq=dtgrid_nbq/float(nsdtnbq)
        ndtfast=ndtfast*nsdtnbq
       endif 
# endif  /* KNHINT_3M */
! !
! !
! !********************************************************************
! !********************************************************************
! ! PART D:  Fast-mode post-processing...
! !********************************************************************
! !********************************************************************
! !
! !********************************************************************
! ! PART D.1: Compute time average & adjust depth-averaged variables
! !********************************************************************
! ! 
#include "k3fast_post.h"
! !
! !********************************************************************
! ! PART D.2: Adjust ZETA
! !
! !  using depth-integrated continuity equation and update grid
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
! !********************************************************************
! !
# include "k3fast_zeta_correct.h"
! !
! !********************************************************************
! ! PART D.3: AGRIF: Apply conservation requirements for nesting
! !********************************************************************
! !
# ifdef AGRIF
#  include "k3fast_AGRIF1.h"
# endif /* AGRIF */
! !
! !********************************************************************
! ! PART D.4: Adjust Hz(m+1) for internal mode
! !           by inverting internal continuity equation
! !********************************************************************
! !
# ifdef NBQ_HZCORRECT  
#  include "k3fast_hz_correct.h"
# endif  /* NBQ_HZCORRECT */   
! !
! !
! !********************************************************************
! !********************************************************************
! ! PART E:  TEST FOR CFL VIOLATION. IF SO, PRINT AND STOP
! !********************************************************************
! !********************************************************************
! !
      VMAXL=100.
      VMAX=0.
!$acc parallel loop if(compute_on_device) independent
!$acc& collapse(2)
!$ass& reduction(max:VMAX)
!$acc&      copyout( VMAX, imax, jmax )
      do j=Jstr,Jend
        do i=Istr,Iend
          cff1=ubar(i,j,knew)
          cff2=vbar(i,j,knew)
          cff=max(abs(cff1),abs(cff2))
# ifdef K3FAST_RHO
          do k=1,N
           cff1=rho_nbq(i,j,k)
           cff2=rho_nbq(i,j,k)
           cff=max(cff,abs(cff1),abs(cff2))
          enddo
# endif
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
#endif /* K3FAST */
