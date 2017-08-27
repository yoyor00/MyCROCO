! $Id: step2d.F 1615 2014-12-17 13:27:07Z rblod $
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
#define PRED_COUPLED_MODE

      subroutine step2d (tile)
      implicit none
      integer tile, trd
#include "param.h"
#include "private_scratch.h"
C$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
      trd=0
C$    trd=omp_get_thread_num()
      
      call step2D_FB_tile ( Istr,Iend,Jstr,Jend, A2d(1,1,trd)
     &                    , A2d(1, 2,trd), A2d(1, 3,trd), A2d(1, 4,trd)
     &                    , A2d(1, 5,trd), A2d(1, 6,trd), A2d(1, 7,trd)
     &                    , A2d(1, 8,trd), A2d(1, 9,trd) 
#ifdef M2_HADV_UP3
     &                    , A2d(1,10,trd), A2d(1,11,trd) 
#else
     &                    , A2d(1,10,trd), A2d(1,11,trd) 
     &                    , A2d(1,12,trd), A2d(1,13,trd)     
#endif
#if defined NBQ
#if defined NBQ_IJK
     &                    , A3d(1,1,trd), A3d(1,2,trd)
     &                    , A3d(1,3,trd), A3d(1,4,trd)
#ifdef NBQ_ZETAW
     &                    , A3d(1,5,trd), A3d(1,6,trd)
#endif
     &                    , A2d(1,14,trd)                     ! A2d(:,14-25,trd) used as working arrays for NBQ
#ifndef NBQ_ZETAW
     &                    , A3d(1,6,trd)                      ! A3d(:,6-,trd) used as working arrays for NBQ
#endif
#endif
     &                    , A2d(1,26,trd), A2d(1,27,trd)
     &                    , A2d(1,28,trd), A2d(1,29,trd)
#ifndef NBQ_ZETAW
     &                    , A2d(1,30,trd), A2d(1,31,trd)
#endif
#endif
     &                    )

      return
      end
 
      subroutine step2D_FB_tile (Istr,Iend,Jstr,Jend, zeta_new,
     &                           Dnew,rubar,rvbar, 
     &                           Drhs, UFx,UFe,
     &                           VFx,VFe
#ifdef M2_HADV_UP3
     &                          ,wrk1,wrk2
#else                           
     &                          ,urhs,vrhs 
     &                          ,DUon,DVom
#endif
#if defined NBQ
#if defined NBQ_IJK
     &                          ,Hzw_half_nbq_inv, Hzr_half_nbq_inv 
     &                          ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
#ifdef NBQ_ZETAW
     &                          ,Hzu_half_qdmu,Hzv_half_qdmv
#endif
     &                          ,work_nbq
#ifndef NBQ_ZETAW
     &                          ,work3d_nbq
#endif
#endif
     &                          ,ruext_nbq_2d_sum, rvext_nbq_2d_sum
     &                          ,ruext_nbq_2d_old, rvext_nbq_2d_old
#ifndef NBQ_ZETAW
     &                          ,Drhs_half_u,Drhs_half_v
#endif
#endif
     &                          )

!
!======================================================================
! Perform one time step for barotropic mode (free-surface and baro-
! tropic 2D momentum equations) using Generalized Forward-Backward
! AB3-AM4 algorithm. Also calculate fast-time averages to interact
! with baroclinic mode.
!
! Forward-Backward schemes for coupled equations (here free-surface zeta
! and 2D momentum ubar and vbar) are stable for longer time-steps 
! (alpha_max>1) than synchronous schemes while only needing evaluation 
! of each r.h.s. term ones. Rather than the original Euler FB schemes, 
! a generalized FB scheme with an AB3-AM4 step is used here for its added 
! robustness (for Coriolis and advection terms). It is a combination of 
! AB3-like step for zeta and AM4-like step for ubar and vbar but with 
! coefficients chosen for maximum stability (e.g. beta=5/12 in original 
! AB3 is replaced by 0.281105):
!
! AB3 Forward step:
! ----------------
! ubarrhs = (3/2+beta)*ubar(n) - (1/2+2*beta)*ubar(n-1) + beta*ubar(n-2)
! vbarrhs = (3/2+beta)*vbar(n) - (1/2+2*beta)*vbar(n-1) + beta*vbar(n-2)
! zeta(n+1) = zeta(n) + RZETA(ubarhs,vbarrhs)
!
! AM4 Backward step:
! -----------------
! zetarhs = (1/2+gamma+2*epsilon) * zeta(n+1) + 
!           (1/2-2*gamma-3*epsilon) * zeta(n) +
!           gamma * zeta(n-1) + epsilon * zeta(n-2)
! ubar(n+1) = ubar(n) + RUBAR(zetarhs,ubarrhs,vbarrhs)
! vbar(n+1) = vbar(n) + RVBAR(zetarhs,ubarrhs,vbarrhs)
!
! beta=0.281105   gamma=0.0880    epsilon=0.013 
!
! Note that the first 2D step is performed with the original Euler FB 
! scheme and the second 2D step with a AB2-AM3 FB scheme with coefficients 
! again chosen for maximum stability.
!
! Reference:
! ---------
! Shchepetkin, A.F., and J.C. McWilliams, 2009: Computational kernel 
! algorithms for fine-scale, multiprocess, longtime oceanic simulations. 
! Pp. 119–182 in Handbook of Numerical Analysis: Computational Methods 
! for the Atmosphere and Oceans. R.M. Teman and J.J. Tribbia, eds, 
! Elsevier Science.
!
!======================================================================
!
#ifdef NBQ
      use module_nh
      use module_nbq
#endif

      implicit none
#include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j, kbak, kold,
#ifdef MPI
     &         err,
#endif
#ifdef PSOURCE
     &        is,
#endif
     &        imin,imax,jmin,jmax
      real sum_c
      real    VMAX,VMAXL
      real zeta_new(PRIVATE_2D_SCRATCH_ARRAY),  cff,
     &         Dnew(PRIVATE_2D_SCRATCH_ARRAY),  cff0,
     &        rubar(PRIVATE_2D_SCRATCH_ARRAY),  cff1,
     &        rvbar(PRIVATE_2D_SCRATCH_ARRAY),  cff2,
     &         Drhs(PRIVATE_2D_SCRATCH_ARRAY),  cff3,
     &          UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &          UFe(PRIVATE_2D_SCRATCH_ARRAY),  DUnew,
     &          VFx(PRIVATE_2D_SCRATCH_ARRAY),  DVnew,
     &          VFe(PRIVATE_2D_SCRATCH_ARRAY)
#ifdef M2_HADV_UP3
      real     wrk1(PRIVATE_2D_SCRATCH_ARRAY),
     &         wrk2(PRIVATE_2D_SCRATCH_ARRAY)
      real     curvX,curvE,cffE,cffX, gamma
      parameter (gamma = -0.25)
#else
      real     urhs(PRIVATE_2D_SCRATCH_ARRAY),  
     &         vrhs(PRIVATE_2D_SCRATCH_ARRAY),
     &         DUon(PRIVATE_2D_SCRATCH_ARRAY),
     &         DVom(PRIVATE_2D_SCRATCH_ARRAY)
#endif
#if defined NBQ
#if defined NBQ_IJK
       real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N) 
       real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N) 
       real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
       real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real work_nbq(PRIVATE_2D_SCRATCH_ARRAY,11)
#ifndef NBQ_ZETAW
       real work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,N,5)
#endif

#endif
       real ruext_nbq_2d_sum(PRIVATE_2D_SCRATCH_ARRAY)
       real rvext_nbq_2d_sum(PRIVATE_2D_SCRATCH_ARRAY)
       real ruext_nbq_2d_old(PRIVATE_2D_SCRATCH_ARRAY)
       real rvext_nbq_2d_old(PRIVATE_2D_SCRATCH_ARRAY)
#ifndef NBQ_ZETAW
       real Drhs_half_u(PRIVATE_2D_SCRATCH_ARRAY)
       real Drhs_half_v(PRIVATE_2D_SCRATCH_ARRAY)
#endif
#endif
#ifdef M2FILTER_NONE
      real :: myepsilon,mybeta,myalpha,mygamma
#endif
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "coupling.h"
#include "forces.h"
#ifdef MRL_WCI
      real vstu,ustv,dudx,dvdx,dude,dvde
#endif
#include "mixing.h"
#include "climat.h"
#include "scalars.h"
#include "sources.h"
#ifdef POT_TIDES
# include "tides.h"
#endif
#ifdef AGRIF
# include "zoom.h"
        integer irhox, irhoy, irhot
#endif
#if defined INTERNAL || defined BODYTIDE
      real U0, omega
#endif
#ifdef WET_DRY
      real cff1_WD,cff2_WD
#endif
#ifdef NBQ
# include "nbq.h"
      integer k
#endif
#ifdef NBQ_ZETAW
      integer resetfromrhobar
#endif  	  

      real dthetadiv_nbqdz_u(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real dthetadiv_nbqdz_v(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real dthetadiv_nbqdz_w(PRIVATE_2D_SCRATCH_ARRAY,2) 

      real dZdxq_u(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real dZdyq_v(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real FY(PRIVATE_2D_SCRATCH_ARRAY) 

      real zwrk1(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real zwrk2(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real zwrk3(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real zwrk4(PRIVATE_2D_SCRATCH_ARRAY,2) 
      real zwrk5(PRIVATE_2D_SCRATCH_ARRAY) 
      real FC(PRIVATE_1D_SCRATCH_ARRAY,0:N) 
      real FX(PRIVATE_2D_SCRATCH_ARRAY) 
	  integer :: k1, k2, kp1

      real DC(PRIVATE_1D_SCRATCH_ARRAY,0:N+1) 
      real CF(PRIVATE_1D_SCRATCH_ARRAY,0:N+1)
      real WORK(PRIVATE_2D_SCRATCH_ARRAY)
      real dum_s
!
#ifndef NBQ_MASS
#  define Hzr_half_nbq Hz
#endif
!
#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif
!
#ifdef MPI
#include "mpi_roms.h"
      include 'mpif.h'
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif   

#include "compute_auxiliary_bounds.h"

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
!======================================================================
!======================================================================
!======================================================================
! Approximations to speed-up computations
!======================================================================
!======================================================================
!======================================================================
!
!---------------------------------------
! d./ds terms:
!---------------------------------------
!
# ifdef NBQ_NODS
#  if defined NBQ_IMP && defined NBQ_ZETAW
#   define NSTEP_DS mod(iic,1000)==1.and.iif==1
#  else
#   define NSTEP_DS mod(iic,200)==1.and.iif==1
#  endif
# endif
!
!---------------------------------------
! Grid computations:
!---------------------------------------
!
!# define NSTEP_GRID .or.mod(iif,1)==0
# define NSTEP_GRID .or.mod(iif,nfast)==0
!
!======================================================================
!======================================================================
!======================================================================
! Initializations
!======================================================================
!======================================================================
!======================================================================
!
         if (IstrU.gt.Iend) then
            do j=Jstr,Jend
	    do i=Istr,Iend+1
                FX(i,j)=0.
            enddo
            enddo
         endif
         if (JstrV.gt.Jend) then
            do j=Jstr,Jend+1
            do i=Istr,Iend 
               FY(i,j)=0.
            enddo
            enddo
         endif
!
!======================================================================
!======================================================================
!======================================================================
! AB3 Forward Step: Computes zeta(n+1) and zetarhs
!======================================================================
!======================================================================
!======================================================================
!
!---------------------------------------
!  AB3 Coef.
!---------------------------------------
!
#ifdef M2FILTER_NONE
# ifdef TANK
         myalpha   = 0.01
# else
         myalpha   = 0.5
# endif
         mybeta    = 0.281105
         myepsilon = 0.00976186 - 0.13451357*myalpha
         mygamma   = 0.08344500 - 0.51358400*myalpha
#endif
!
#ifdef NBQ_ZETAW
        if (iic==1.and.iif==1) then
         cff4= 1.
         cff5= 0.
         cff6= 0.
        elseif (iic==1.and.iif==2) then
         cff4= 1.5
         cff5= -0.5
         cff6= 0.
        else
         cff4= 1.5+mybeta
         cff5=-2.0*mybeta-0.5
         cff6= mybeta
        endif
     !  cff4= 1.
     !  cff5= 0.
     !  cff6= 0.
#endif
!
! Preliminary step: compute total depth (meters) of the water
! ----------------- column and vertically integrated mass fluxes
! which are needed to compute divergence in rhs_zeta and as input
! data to compute nonlinear advection terms for the barotropic
! momentum equations.
!
! Output: (urhs,vrhs), (DUon,DVom) --> (DU_avg2,DV_avg2)
!
#include "step2d_ab3_init.h"
!
!-----------------------------------------------------------------------
! Advance free-surface:   Compute zeta_new, which is at new time
!-------- ---- --------   step, and interpolate half-step backward
! for the subsequent computation of barotropic pressure-gradient
! terms. It should be noted that because Forward Euler step is used 
! to update zeta during the first barotropic step, the pressure term  
! must be computed via Backward step to keep it numerically stable.
! However, this would interfere with the computation of forcing terms
! "rufrc,rvfrc" because computation of pressure gradient in 3D mode
! uses exactly the initial value of "zeta", rather than value changed
! by one barotropic time step.  To resolve this conflict, the
! pressure gradient term computation during the first barotropic
! step is computed in two stages: first use just zeta(:,:,kstp) to
! insure exact consistency with 3D mode; then, after "rufrc,rvfrc"
! are finalized, add a correction term based on the difference
! zeta_new(:,:)-zeta(:,:,kstp) to "rubar,rvbar" to make them
! consistent with Backward step for pressure gradient terms.
!-----------------------------------------------------------------------
!
#ifndef NBQ_ZETAW
!
!---------------------------------------
! Zeta is computed from depth-averaged
!  equation.
!---------------------------------------
!
# include "step2d_zeta2d.h"
!
#else /* NBQ_ZETAW */
!---------------------------------------
! Zeta is computed from W_surf
! and grid can vary at dtfast.
!---------------------------------------
!
      if (iic==1.and.iif==1) then
!
# include "step2d_grid_ext.h"  
!  
      endif  
!
#endif /* ! NBQ_ZETAW*/
!
!---------------------------------------
! Computes zetarhs to use 
!   in momentum equations
!    AB3 now !!!
!---------------------------------------
!
#include "step2d_zetarhs.h"
!
!---------------------------------------
! Update Zeta (ifnot NBQ_ZETAW)
!---------------------------------------
!
#ifndef NBQ_ZETAW
!
# include "step2d_zeta2d_load.h"
!
#endif /* not NBQ_ZETAW */
!
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
# include "step2d_extint_avg12_part1.h"
!
!-----------------------------------------------------------------------
! Computes rhs for momentum equations
!-----------------------------------------------------------------------
!
# include "step2d_momentum_rhs.h"
!
!
#ifdef SOLVE3D
!
      if (FIRST_2D_STEP) then
!
!-----------------------------------------------------------------------
! First 2D step: prepares forcing and coupling
! Computes rufrc & rvfrc
! Coupling between 2D and 3D parts.
!--------- ------- -- --- -- ------
! Before the predictor step of the first barotropic time step
! arrays "rufrc" and "rvfrc" contain vertically integrals of the
! 3D right-hand-side terms for the momentum equations (including
! surface and bottom stresses, if so prescribed).
!
! During the first barotropic time step connvert them into forcing
! terms by subtracting the fast-time "rubar" and "rvbar" from them;
! These forcing terms are then extrapolated forward in time using
! optimized Adams-Bashforth weights, so that the resultant rufrc
! and rvfrc are centered effectively at time n+1/2. From now on,
! these newly computed forcing terms will remain constant during
! the fast time stepping and will added to "rubar" and "rvbar"
! during all subsequent barotropic time steps.
!-----------------------------------------------------------------------
!
# ifdef NBQ_ZETAW
#  include "step2d_zetaw_frc.h"
# else
#  include "step2d_zeta2d_frc.h"
# endif
!
      endif   !<-- FIRST_2D_STEP
!
#endif /* SOLVE3D */
!
#ifdef NBQ
!
!-----------------------------------------------------------------------
! Coupling EXT / NBQ
!-----------------------------------------------------------------------
!
# ifdef NBQ_ZETAW
#  include "step2d_extnbq_part1_zetaw.h"
# else
#  include "step2d_extnbq_part1_zeta2d.h"
# endif
!
!
!-----------------------------------------------------------------------
! Resolve 3D NBQ equation for ndtnbq time steps
! This provides NBQ-filtered rhs terms for the barotropic equation.
! Then perform 2D filter of these terms to get forcing terms for 
! the total equation.
!-----------------------------------------------------------------------
!
!
# include "step2d_nbq_prep.h"  
!
!-----------------------------------------------------------------------
! NBQ time-steps
!-----------------------------------------------------------------------
!     
# ifdef NBQ_ZETAW
#  include "step3d_nbq_zetaw.h"
# else
#  include "step3d_nbq_zeta2d.h"
#endif
!
!======================================================================
!======================================================================
!======================================================================
! AM4 Backward Step: Computes ubar,vbar at time n+1
!                    first computes rubar,rvbar
!======================================================================
!======================================================================
!======================================================================
      
!
!---------------------------------------
!  AM4 Coef.
!---------------------------------------
!
#ifdef NBQ_ZETAW
      if (iic==1.and.iif==1) then
        cff4= 1.
        cff5= 0.
        cff6= 0.
        cff7= 0.
      elseif (iic==1.and.iif==2) then
         cff4=0.5+2.*myalpha
         cff5=1.-cff4
         cff6=0.
         cff7=0.
      elseif (iic==1.and.iif==3) then
         cff4=0.5+mygamma+2.*myalpha
         cff5=1.-cff4-mygamma
         cff6=mygamma
         cff7=0.
      else
         cff4=0.5+2.*myepsilon+mygamma+2.*myalpha
         cff5=1.-cff4-mygamma-myepsilon
         cff6=mygamma
         cff7=myepsilon
      endif
   !    cff4= 1.
   !    cff5= 0.
   !    cff6= 0.
   !    cff7= 0.
#endif
!
!-----------------------------------------------------------------------
! Computes new zeta and update grid
!    (AM4)
!-----------------------------------------------------------------------
!
# ifdef NBQ_ZETAW
#  include "step2d_zetaw.h"   
# endif
!
!-----------------------------------------------------------------------
!
! Get filtered rhs terms (! M2FILTER_NONE)
! and multiply by dx*dy to get units of rho*Hz*dx*dy*ru
!
!-----------------------------------------------------------------------
!
# include "step2d_extnbq_part2.h"

!
!-----------------------------------------------------------------------
!
!  Gets rid of NBQ NBQ2EXT & NBQ2INT coupling
!
!-----------------------------------------------------------------------
!
# ifdef NBQ_NOCOUPLING 
          rubar_nbq   =0.   ! 2D
          rvbar_nbq   =0.
          ru_nbq_ext  =0.   ! 3D
          rv_nbq_ext  =0.
          rw_nbq_ext  =0.
          ru_nbq_avg2 =0.
          rv_nbq_avg2 =0.
          rw_nbq_avg2 =0.
#  ifdef NBQ_MASS
          rhobar_nbq  =1.
          rho_nbq_ext =1.
          rho_nbq_avg1=1.
#  endif
# endif
!
#endif  /* NBQ */   
!#ifdef SOLVE3D
!
!     if (FIRST_2D_STEP) then
!
!-----------------------------------------------------------------------
! First 2D step: prepares forcing and coupling
!-----------------------------------------------------------------------
!
!----------------------------------------
! Since coupling requires that pressure 
! gradient term is computed
! using zeta(:,:,kstp) instead of 
! zeta_new(:,:) needed to achieve
! numerical stability, apply compensation 
! to shift pressure gradient
! terms from "kstp" to "knew": in essense, 
! convert the fist 2D step
! from Forward Euler to Forward-Backward].
!----------------------------------------
!
!# include "step2d_zeta_update.h"
!
!      endif   !<-- FIRST_2D_STEP
!
!#endif /* SOLVE3D */
!
!
!-----------------------------------------------------------------------
! Exchanges NBQ coupling
!-----------------------------------------------------------------------
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
# ifdef NBQ_MASS
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   rhobar_nbq_avg1(START_2D_ARRAY))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                   rho_nbq_avg1(START_2D_ARRAY,1))
# endif
# ifdef NBQ
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                   ru_nbq_avg2(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                   rv_nbq_avg2(START_2D_ARRAY,1))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
     &                   rw_nbq_avg2(START_2D_ARRAY,0))
# endif
#endif
!
!-----------------------------------------------------------------------
! Perform time step for the 2D momentum equations. 

! Also compute fast-time averaged barotropic mass fluxes. 
! Doing so on the fly yields a more computationally dense code and 
! eliminates repeated multiplication by Dnew (since mass fluxes are 
! actually available as volatile variables DUnew,DVnew at this moment. 
! However, as a result of this arrangement, a special code is needed 
! to compute fast-time averages along the physical boundaries, which is 
! done below.
!-----------------------------------------------------------------------
!
! Output: (ubar,vbar), (DU_avg1,DV_avg1)
!
!#ifndef NBQ_ZETAW
# include "step2d_momentum_upd.h"
!!#else
!#
!#endif
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
#ifndef NBQ_ZETAW
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew))
#endif
#ifndef M2FILTER_NONE
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg1(START_2D_ARRAY,nnew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg1(START_2D_ARRAY,nnew))
#endif
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg2(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg2(START_2D_ARRAY))
# if defined MRL_WCI && defined WET_DRY
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ust2d(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vst2d(START_2D_ARRAY))
# endif
#endif
!
!-----------------------------------------------------------------------
!  Debugging ubar,vbar
!-----------------------------------------------------------------------
!
#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
#ifndef NBQ_ZETAW
       call check_tab2d(ubar(:,:,knew),'ubar step2d #2','u')
       call check_tab2d(vbar(:,:,knew),'vbar step2d #2','v')
#endif
C$OMP END MASTER       
#endif   
!
!-----------------------------------------------------------------------
!  Apply conservation requirements for nesting
!-----------------------------------------------------------------------
!
#ifdef AGRIF
# include "step2d_extint_avg12_part2.h"    
#endif /* AGRIF */
!
!      
!======================================================================
!======================================================================
!======================================================================
!    THE END
!======================================================================
!======================================================================
!======================================================================
!
!
!-----------------------------------------------------------------------
!  TEST FOR CFL VIOLATION. IF SO, PRINT AND STOP
!-----------------------------------------------------------------------
!
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
#ifdef MPI
            imax=i+iminmpi-1
            jmax=j+jminmpi-1
#else
            imax=i
            jmax=j
#endif
          ENDIF
        enddo
      enddo
      IF (VMAX.GT.VMAXL) THEN
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =   STEP2D:   ABNORMAL JOB END        = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '
        if (VMAX.eq.666.) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)') 
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')  
     &                                       '  IMAX JMAX  =',imax,jmax
#ifdef SOLVE3D
        write(stdout,'(A,2I6/)')
     &                                         '  IINT IEXT  =',iic,iif
#else
        write(stdout,'(A,I6/)')    '  IIC        =',iic
#endif
        may_day_flag=1
#ifdef MPI
        call mpi_abort (MPI_COMM_WORLD, err)
#else
        stop                            !-->  EXIT
#endif
      ENDIF

!       ! Fin CN
!       if (iic.eq.10) then
!#ifdef MPI       
!       call MPI_Finalize (i)
!       if (mynode.eq.0) write(6,*) 'fin step2d'
!#else
!       write(6,*) 'fin step2d'
!#endif       
!       stop 
!       endif

#if defined DEBUG_NBQ && defined MPI
!        if (iic.eq.10) then
!         call MPI_Finalize (i)
!         if (mynode.eq.0) write(6,*) 'Stop in step2d: debug NBQ'
!         stop
!        endif
#endif

  
      return
      end
