! $Id: step3d_fast_zetaw.F 1615 2014-12-17 13:27:07Z rblod $
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

#if defined NBQ && defined NBQ_ZETAW

#define PRED_COUPLED_MODE

      subroutine step3d_fast_zetaw (tile)
      implicit none
      integer tile, trd
#include "param.h"
#include "private_scratch.h"
C$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
      trd=0
C$    trd=omp_get_thread_num()
!      
!======================================================================
!      step3d_fast_zetaw
!======================================================================
!
      call step3d_fast_zetaw_tile ( Istr,Iend,Jstr,Jend, A2d(1,1,trd)
     &                    , A2d(1, 2,trd), A2d(1, 3,trd), A2d(1, 4,trd)
     &                    , A2d(1, 5,trd), A2d(1, 6,trd), A2d(1, 7,trd)
     &                    , A2d(1, 8,trd), A2d(1, 9,trd) 
#ifdef M2_HADV_UP3
     &                    , A2d(1,10,trd), A2d(1,11,trd) 
#else
     &                    , A2d(1,10,trd), A2d(1,11,trd) 
     &                    , A2d(1,12,trd), A2d(1,13,trd)     
#endif
     &                    , A3d(1,1,trd), A3d(1,2,trd)
     &                    , A3d(1,3,trd), A3d(1,4,trd)
     &                    , A3d(1,5,trd), A3d(1,6,trd)
     &                    , A2d(1,14,trd)                     ! A2d(:,14-25,trd) used as working arrays for NBQ
     &                    , A2d(1,26,trd), A2d(1,27,trd)
     &                    , A2d(1,28,trd), A2d(1,29,trd)
     &                    )

      return
      end
! 
!======================================================================
!      step3d_fast_zetaw_tile
!======================================================================
!
      subroutine step3d_fast_zetaw_tile (Istr,Iend,Jstr,Jend, zeta_new,
     &                           Dnew,rubar,rvbar, 
     &                           Drhs, UFx,UFe,
     &                           VFx,VFe
#ifdef M2_HADV_UP3
     &                          ,wrk1,wrk2
#else                           
     &                          ,urhs,vrhs 
     &                          ,DUon,DVom
#endif
     &                          ,Hzw_half_nbq_inv, Hzr_half_nbq_inv 
     &                          ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
     &                          ,Hzu_half_qdmu,Hzv_half_qdmv
     &                          ,work_nbq
     &                          ,ruext_nbq_2d_sum, rvext_nbq_2d_sum
     &                          ,ruext_nbq_2d_old, rvext_nbq_2d_old
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
       real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N) 
       real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N) 
       real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
       real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real work_nbq(PRIVATE_2D_SCRATCH_ARRAY,11)
       real ruext_nbq_2d_sum(PRIVATE_2D_SCRATCH_ARRAY)
       real rvext_nbq_2d_sum(PRIVATE_2D_SCRATCH_ARRAY)
       real ruext_nbq_2d_old(PRIVATE_2D_SCRATCH_ARRAY)
       real rvext_nbq_2d_old(PRIVATE_2D_SCRATCH_ARRAY)
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
#include "nbq.h"
      integer k
      integer resetfromrhobar
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
      integer i1_a,i2_a,j1_a,j2_a

      real DC(PRIVATE_1D_SCRATCH_ARRAY,0:N+1) 
      real CF(PRIVATE_1D_SCRATCH_ARRAY,0:N+1)
      real WORK(PRIVATE_2D_SCRATCH_ARRAY)
      real dum_s,dum_max
      real sum_r(PRIVATE_2D_SCRATCH_ARRAY)
      real,parameter :: thetaimp_nbq = 1.
      real,parameter :: gammau=0.45
      !real,parameter :: gammau=0.55
      !real,parameter :: gammau=1.0
      real,parameter :: gammau_2=(1./3.)*(1.-gammau)
      real,dimension(GLOBAL_2D_ARRAY,0:N) :: qdmw_nbq_old
!
#ifndef NBQ_MASS
#  define Hzr Hz
#endif
!
#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif
!
#ifdef MPI
# include "mpi_roms.h"
      include 'mpif.h'
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif   

#include "compute_auxiliary_bounds.h"

#ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
#else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
#endif
!
#ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
#else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
#endif

#define zwrk UFx
#define rzeta  UFe
#define rzeta2  VFe
#define rzetaSA VFx

!
!*********************************************************************
! Settings of several approximations to speed-up computations
!*********************************************************************
!
!----------------------------------------------------------------------
! d./ds terms:
!----------------------------------------------------------------------
!
#ifdef NBQ_NODS
# if defined NBQ_IMP && defined NBQ_ZETAW
#   define NSTEP_DS mod(iic,1)==1.and.iif==1
# else
#   define NSTEP_DS mod(iic,1)==1.and.iif==1
# endif
#endif
!

!----------------------------------------------------------------------
! Grid computations:
!----------------------------------------------------------------------
!
#ifdef NBQ_GRIDEXT
#  define NSTEP_GRID .or.mod(iif,1)==0
#else
#  define NSTEP_GRID .or.iif==nfast
#endif
!
!*********************************************************************
! Initializations
!*********************************************************************
!
      if (iic==ntstart.and.iif==1) then
!
!----------------------------------------------------------------------
! Initializes qdmu, qdmv, qdmw: fast-mode momentum is equal to
! internal mode values.
!
! This should be temporary.
!----------------------------------------------------------------------
!
        do k=1,N
          do j=Jstr-1,Jend+1
          do i=Istr,iend+1
#ifdef NBQ_MASS
              qdmu_nbq(i,j,k)=(1.+0.5*(rho(i,j,k)+rho(i-1,j,k))/rho0)
     &              *0.5*u(i,j,k,nrhs)*(hz(i,j,k)+hz(i-1,j,k))
#else
              qdmu_nbq(i,j,k)=0.5*u(i,j,k,nrhs)*(hz(i,j,k)+hz(i-1,j,k))
#endif
          enddo
         enddo
        enddo

        do k=1,N
          do j=Jstr,Jend+1
          do i=Istr-1,iend+1
#ifdef NBQ_MASS
              qdmv_nbq(i,j,k)=(1.+0.5*(rho(i,j,k)+rho(i,j-1,k))/rho0) 
     &              *0.5*v(i,j,k,nrhs)*(hz(i,j,k)+hz(i,j-1,k))
#else
              qdmv_nbq(i,j,k)=0.5*v(i,j,k,nrhs)*(hz(i,j,k)+hz(i,j-1,k))
#endif
          enddo
         enddo
        enddo

        do k=1,N-1
          do j=Jstr-1,Jend+1
          do i=Istr-1,iend+1
#ifdef NBQ_MASS
               qdmw_nbq(i,j,k)=(1.+0.5*(rho(i,j,k)+rho(i,j,k+1))/rho0) 
     &            *0.5*wz(i,j,k,nrhs)*(hz(i,j,k)+hz(i,j,k+1))
#else
              qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*(hz(i,j,k)+hz(i,j,k+1))
#endif
          enddo
         enddo
        enddo
        k=0 
          do j=Jstr-1,Jend+1
          do i=Istr-1,iend+1
#ifdef NBQ_MASS
            qdmw_nbq(i,j,k)=(1.+rho(i,j,k+1)/rho0) 
     &            *0.5*wz(i,j,k,nrhs)*hz(i,j,k+1)
#else
            qdmw_nbq(i,j,k)= 0.5*wz(i,j,k,nrhs)*hz(i,j,k+1)
#endif
          enddo
        enddo
        k=N 
          do j=Jstr-1,Jend+1
          do i=Istr-1,iend+1
#ifdef NBQ_MASS
            qdmw_nbq(i,j,k)=(1.+rho(i,j,k)/rho0) 
     &           *0.5*wz(i,j,k,nrhs)*hz(i,j,k)
#else
            qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)*hz(i,j,k)
#endif
          enddo
        enddo

      endif

#ifdef NBQ_TRACERS
         if (iic.ge.2.and.iif==1) then
            Hz = Hz_tra
            call grid_coef_nh(
     &       Istr,Iend,Jstr,Jend,
     &       Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &       Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &       Hzu_half_qdmu, Hzv_half_qdmv                                     
     &        )
        endif
#endif

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
!----------------------------------------------------------------------
! Constants for time-stepping:
!----------------------------------------------------------------------
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
        if (FIRST_TIME_STEP.and.iif==1) then
         cff4= 1.
         cff5= 0.
         cff6= 0.
        elseif (FIRST_TIME_STEP.and.iif==2) then
         cff4= 1.5
         cff5= -0.5
         cff6= 0.
        else
         cff4= 1.5+mybeta
         cff5=-2.0*mybeta-0.5
         cff6= mybeta
        endif
      if (FIRST_2D_STEP) then         ! Meaning of temporal indices
        kbak=kstp                     ! ------- -- -------- -------
        kold=kstp                     ! m-2     m-1      m      m+1
        cff1=1.                       ! kold    kbak    kstp    knew
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1.                  ! Logically AB2-AM3 forward-backward 
        cff2=0.                  ! scheme with coefficients chosen for
        cff3=0.                  ! maximum stability ... (see below)
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
#ifdef M2FILTER_NONE
        cff1= 1.5+mybeta
        cff2=-2.0*mybeta-0.5
        cff3= mybeta
#else
        cff1= 1.781105
        cff2=-1.06221
        cff3= 0.281105
#endif
      endif
!
!----------------------------------------------------------------------
!  AGRIF
!----------------------------------------------------------------------
!
#ifdef AGRIF
        irhox = Agrif_Irhox()
        irhoy = Agrif_Irhoy()
        irhot = Agrif_Irhot()
#endif
!
!----------------------------------------------------------------------
!  Time-step indices:
!----------------------------------------------------------------------
!
      if (FIRST_TIME_STEP.and.iif==1) then
        kbak2=4
        kold2=3
      else
        kold2=kbak2
        kbak2=kbak2+1
        if (kbak2==5) kbak2=1
      endif

#if defined M2FILTER_NONE 
      if (FIRST_2D_STEP) then
# ifdef SOLVE3D
        cff0=0.                !---> Compute pressure-gradient
        cff1=1.                !  terms using just zeta(:,:,kstp)
# else
        cff0=1.
        cff1=0.
# endif
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        cff0= 1.0833333333333    ! Logically AB2-AM3 forward-backward
        cff1=-0.1666666666666    ! scheme with coefficients chosen for
        cff2= 0.0833333333333    ! maximum stability, while maintaining
        cff3= 0.                 ! third-accuracy; alpha_max=1.73
      else
# ifdef M2FILTER_NONE
        cff0=0.5+2.*myepsilon+mygamma+2.*myalpha
        cff1=1.-cff0-mygamma-myepsilon
        cff2=mygamma
        cff3=myepsilon
# else
        cff0=0.614
        cff1=0.285
        cff2=0.088
        cff3=0.013
# endif
      endif

      if (FIRST_TIME_STEP.and.iif.eq.1) then
        cff8=1.D0
        cff9=0.D0
        cff10=0.D0
      elseif (FIRST_TIME_STEP.and.iif.eq.1+1) then
        cff8=1.5D0
        cff9=-0.5D0
        cff10=0.D0
      else
        cff8= 1.5D0+mybeta
        cff9=-2.0D0*mybeta-0.5D0
        cff10= mybeta
      endif

#endif /* M2FILTER_NONE */

#ifdef RVTK_DEBUG_ADVANCED
       call check_tab2d(zeta(:,:,kstp),'zeta step2d #0','r')
#endif  
!
!*********************************************************************
! Extrapolate (D,u,v) at m+1/2 to compute RHS 
!*********************************************************************
!
#if defined EW_PERIODIC || defined MPI || !defined M2_HADV_UP3
      imin=IstrU-2
      imax=Iend+1
# define IV_EXT_RANGE Istr-1,Iend+1
#else
      imin=max(IstrU-3,0)
      imax=min(Iend+2,Lm+1)
# define IV_EXT_RANGE max(Istr-2,0),min(Iend+2,Lm+1)
#endif
#if defined NS_PERIODIC || defined MPI || !defined M2_HADV_UP3
      jmin=JstrV-2
      jmax=Jend+1
# define JU_EXT_RANGE Jstr-1,Jend+1
#else
      jmin=max(JstrV-3,0)
      jmax=min(Jend+2,Mm+1)
# define JU_EXT_RANGE max(Jstr-2,0),min(Jend+2,Mm+1)
#endif
!
!----------------------------------------------------------------------  
! Total depth / mass  at n+1/2
!----------------------------------------------------------------------
!
      do j=jmin,jmax
        do i=imin,imax
#ifdef NBQ_MASS
          Drhs(i,j)= cff8 *(zeta(i,j,kstp)+h(i,j))
     &                    *rhobar_nbq(i,j,kstp)
     &              +cff9 *(zeta(i,j,kbak2)+h(i,j))
     &                    *rhobar_nbq(i,j,kbak2)
     &              +cff10*(zeta(i,j,kold2)+h(i,j))
     &                    *rhobar_nbq(i,j,kold2)
#else /* ! NBQ_MASS */
          Drhs(i,j)= cff8 *(zeta(i,j,kstp)+h(i,j))
     &              +cff9 *(zeta(i,j,kbak2)+h(i,j))
     &              +cff10*(zeta(i,j,kold2)+h(i,j))
#endif /* NBQ_MASS */
        enddo
      enddo

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,kstp))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,kstp))
#endif

!
!----------------------------------------------------------------------  
! Depth average velocity  at n+1/2
!----------------------------------------------------------------------
!
      do j=JU_EXT_RANGE
        do i=imin+1,imax
          urhs(i,j)=cff8*ubar(i,j,kstp) +cff9 *ubar(i,j,kbak2)
     &                                   +cff10*ubar(i,j,kold2)
#if defined MRL_WCI && defined MASKING
          urhs(i,j)=urhs(i,j)*umask(i,j)+ust2d(i,j)*(umask(i,j)-1.0)
#endif
          DUon(i,j)=0.5*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
#ifdef MRL_WCI
     &                                                   + ust2d(i,j)
#endif
     &                                                              )
        enddo
      enddo
#undef JU_EXT_RANGE

      do j=jmin+1,jmax
        do i=IV_EXT_RANGE
          vrhs(i,j)=cff8*vbar(i,j,kstp) +cff9 *vbar(i,j,kbak2)
     &                                   +cff10*vbar(i,j,kold2)
#if defined MRL_WCI && defined MASKING
          vrhs(i,j)=vrhs(i,j)*vmask(i,j)+vst2d(i,j)*(vmask(i,j)-1.0)
#endif
          DVom(i,j)=0.5*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
#ifdef MRL_WCI
     &                                                   + vst2d(i,j)
#endif
     &                                                              )
        enddo
      enddo
#undef IV_EXT_RANGE

#ifdef M2_HADV_UP3
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   urhs(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vrhs(START_2D_ARRAY))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   Duon (START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   Dvom(START_2D_ARRAY))
# endif
#endif

#ifdef OBC_VOLCONS
      call set_DUV_bc_tile (Istr,Iend,Jstr,Jend, Drhs, DUon,DVom)
#endif

#ifdef RVTK_DEBUG_ADVANCED
       call check_tab2d(zeta(:,:,kstp),'zeta step2d #1','r')
#endif 

       if (FIRST_TIME_STEP.and.iif==1) then 
!
!----------------------------------------------------------------------  
! Grid-coef to compute operators at first time_step
!----------------------------------------------------------------------
!
        call grid_coef_nh(
     &   Istr,Iend,Jstr,Jend,
     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &   Hzu_half_qdmu, Hzv_half_qdmv                                     
     &   )
      endif
!
!***********************************************************************
! Compute Zeta for rhs (pressure gradient)
!***********************************************************************
!
#ifdef NBQ_MASS
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend

       !   rhobar_nbq(i,j,knew) = 2.*rhobar_nbq(i,j,kstp)-rhobar_nbq(i,j,kbak2)
       !   zwrk(i,j)= zeta(i,j,kstp)/rhobar_nbq(i,j,kstp)-h(i,j)
           zwrk(i,j)=zeta(i,j,kstp)
!          zwrk(i,j)= 
!     &              +cff4*zeta(i,j,kstp)/rhobar_nbq(i,j,kstp)
!     &              +cff5*zeta(i,j,kbak2)/rhobar_nbq(i,j,kbak2)
!     &              +cff6*zeta(i,j,kold2)/rhobar_nbq(i,j,kold2)
!     &              -h(i,j)

# ifdef MASKING
          Dnew(i,j)=(zeta(i,j,kstp)*rmask(i,j)+h(i,j))*rhobar_nbq(i,j,kstp)  ! CAUTION: rhobar_nbq must 1 on Mask !!
          zwrk(i,j)=zwrk(i,j)*rmask(i,j)
# else
          Dnew(i,j)=(zeta(i,j,kstp)+h(i,j))*rhobar_nbq(i,j,kstp)
# endif

# if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
        enddo
      enddo
!
#else /* ! NBQ_MASS */
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
# ifndef NBQ
          zeta_new(i,j)=zeta_new(i,j) SWITCH rmask(i,j)  ! Useful ????
# endif

          zeta(i,j,knew)= zeta(i,j,knew)               ! Useful ????
     &            SWITCH rmask(i,j)                      ! Utilité ?????
#  ifdef MASKING
          Dnew(i,j)=(zeta(i,j,kstp)*rmask(i,j)+h(i,j))
#  else
          Dnew(i,j)=(zeta(i,j,kstp)+h(i,j))
#  endif
          zwrk(i,j)=zeta(i,j,kstp) SWITCH rmask(i,j)  
!         zwrk(i,j)= cff4*zeta(i,j,kstp)
!    &              +cff5*zeta(i,j,kbak2)+cff6*zeta(i,j,kold2)  
!         zwrk(i,j)= cff4*zeta(i,j,knew)+cff5*zeta(i,j,kstp)
!    &              +cff6*zeta(i,j,kbak2)+cff7*zeta(i,j,kold2)  


# if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
        enddo
      enddo

#endif /* NBQ_MASS */

!
!***********************************************************************
! Compute rhs for momentum equations
!***********************************************************************
!
!
!----------------------------------------------------------------------
! Surface pressure gradient
!----------------------------------------------------------------------
!
      cff=0.5*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(zwrk(i-1,j)-zwrk(i,j)))
#endif
#ifdef MRL_WCI
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i-1,j) )
#endif
#if defined POT_TIDES && !defined SOLVE3D
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                   *( Ptide(i,j)-Ptide(i-1,j) )
#endif
     &                                                              )    
! 
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)
#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(zwrk(i,j-1)-zwrk(i,j)))
#endif
#ifdef MRL_WCI
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i,j-1) )
#endif
#if defined POT_TIDES && !defined SOLVE3D
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                   *( Ptide(i,j)-Ptide(i,j-1) )
#endif
     &                                                              )
        enddo
      enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

#undef rzetaSA
#undef rzeta2
#undef rzeta
#undef zwrk
!
!----------------------------------------------------------------------
! Compute horizontal advection terms for momentum equations (2D only)
!-------- ---------- --------- ----- --- -------- --------- --- -----
! NOTE: mathematically necessary (minimal) index ranges for momentum-
! flux components are 
!
!      UFx(IstrU-1:Iend,Jstr:Jend)   VFx(Istr:Iend+1,JstrV:Jend)
!      UFe(IstrU:Iend,Jstr:Jend+1)   VFe(Istr,Iend,JstrV-1,Jend)
!
! however, for the purpose computational efficiency, these ranges are
! unified by suppressing U,V-suffices in order to allow fusion of the
! consecutive loops. This leads to slight increase of the redundant
! computations near western and southern boundaries in non-periodic
! directions. 
!----------------------------------------------------------------------
!
#ifdef UV_ADV
!
# if defined SOLVE3D || !defined M2_HADV_UP3

!----------------------------------------------------------------------
! Centered Second order advection scheme
!
! Numerical diffusion of momentum is implicitely added through 3D
! forcing of advection in rufrc and rvfrc (i.e., diffusion is
! at slow time scale)

!----------------------------------------------------------------------
      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))

          VFx(i+1,j)=0.25*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i+1,j)
#  endif
#  ifdef WET_DRY0
     &                                 *pmask_wet(i+1,j)
#  endif
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          VFe(i,j)=0.25*(DVom(i,j)+DVom(i,j+1))
     &                 *(vrhs(i,j)+vrhs(i,j+1))

          UFe(i,j+1)=0.25*(DVom(i,j+1)+DVom(i-1,j+1))
     &                   *(urhs(i,j+1)+urhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i,j+1)
#  endif
#  ifdef WET_DRY0
     &                                 *pmask_wet(i,j+1)
#  endif
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)

          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo !--> discard UFx,VFe,UFe,VFx, DUon,DVom

# else   
!----------------------------------------------------------------------
!  UP3 (QUICK) spatial advection scheme as in 3D part
!
!  in this case, there is no numerical requirement for explicit viscosity
!-----------------------------------------------------------------------
#  define uxx wrk1
#  define Huxx wrk2
!
#  ifdef EW_PERIODIC
#   define IU_EXT_RANGE IstrU-1,Iend+1
#  else
#   ifdef MPI
        if (WEST_INTER) then
          imin=IstrU-1
        else
          imin=max(IstrU-1,2)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
#    define IU_EXT_RANGE imin,imax
#   else
#    define IU_EXT_RANGE max(IstrU-1,2),min(Iend+1,Lm)
#   endif
#  endif

        do j=Jstr,Jend
          do i=IU_EXT_RANGE
            uxx(i,j)=urhs(i-1,j)-2.*urhs(i,j)+urhs(i+1,j)
            Huxx(i,j)=Duon(i-1,j)-2.*Duon(i,j)+Duon(i+1,j)
          enddo
        enddo
#   undef IU_EXT_RANGE
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=Jstr,Jend
            uxx(1,j)=uxx(2,j)
            Huxx(1,j)=Huxx(2,j)
          enddo
        endif
        if (EASTERN_EDGE) then
#   ifdef MPI        
          do j=Jstr,Jend
            uxx(Lmmpi+1,j)=uxx(Lmmpi,j)
            Huxx(Lmmpi+1,j)=Huxx(Lmmpi,j)
          enddo
#   else
          do j=Jstr,Jend
            uxx(Lm+1,j)=uxx(Lm,j)
            Huxx(Lm+1,j)=Huxx(Lm,j)
          enddo
#   endif          
        endif
#  endif

        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cffX=urhs(i,j)+urhs(i+1,j)
            if (cffX.gt.0.) then
              curvX=uxx(i,j)
            else
              curvX=uxx(i+1,j)
            endif
            UFx(i,j)=0.25*(cffX+gamma*curvX)*( Duon(i,j)+
     &               Duon(i+1,j)-0.125*(Huxx(i,j)+Huxx(i+1,j)))
          enddo
        enddo
#  undef Huxx
#  undef uxx
#  define vee wrk1
#  define Hvee wrk2
      
#  ifdef NS_PERIODIC
#   define JV_EXT_RANGE JstrV-1,Jend+1
#  else
#   ifdef MPI
        if (SOUTH_INTER) then
          jmin=JstrV-1
        else
          jmin=max(JstrV-1,2)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
#    define JV_EXT_RANGE jmin,jmax
#   else
#    define JV_EXT_RANGE max(JstrV-1,2),min(Jend+1,Mm)
#   endif
#  endif

        do j=JV_EXT_RANGE
          do i=Istr,Iend
            vee(i,j) =vrhs(i,j-1)-2.*vrhs(i,j)+vrhs(i,j+1)
            Hvee(i,j)=Dvom(i,j-1)-2.*Dvom(i,j)+Dvom(i,j+1)
          enddo
        enddo
#   undef JV_EXT_RANGE
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=Istr,Iend
            vee(i,1)=vee(i,2)
            Hvee(i,1)=Hvee(i,2)
          enddo
        endif
        if (NORTHERN_EDGE) then
#   ifdef MPI
          do i=Istr,Iend
            vee(i,Mmmpi+1)=vee(i,Mmmpi)
            Hvee(i,Mmmpi+1)=Hvee(i,Mmmpi)
          enddo
#   else        
          do i=Istr,Iend
            vee(i,Mm+1)=vee(i,Mm)
            Hvee(i,Mm+1)=Hvee(i,Mm)
          enddo
#   endif          
        endif
#  endif
        do j=JstrV-1,Jend
          do i=Istr,Iend
            cffE=vrhs(i,j)+vrhs(i,j+1)
            if (cffE.gt.0.) then
              curvE=vee(i,j)
            else
              curvE=vee(i,j+1)
            endif
            VFe(i,j)=0.25*(cffE+gamma*curvE)*( Dvom(i,j)+
     &               Dvom(i,j+1)-0.125*(Hvee(i,j)+Hvee(i,j+1)))
          enddo
        enddo
#  undef Hvee
#  undef vee
#  define uee wrk1
#  define Hvxx wrk2

#  ifdef NS_PERIODIC
#   define JU_EXT_RANGE Jstr-1,Jend+1
#  else
#   ifdef MPI
        if (SOUTH_INTER) then
          jmin=Jstr-1
        else
          jmin=max(Jstr-1,1)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
#    define JU_EXT_RANGE jmin,jmax
#   else
#    define JU_EXT_RANGE max(Jstr-1,1),min(Jend+1,Mm)
#   endif
#  endif

        do j=JU_EXT_RANGE
          do i=IstrU,Iend
            uee(i,j)=urhs(i,j-1)-2.*urhs(i,j)+urhs(i,j+1)
          enddo
        enddo
#  undef JU_EXT_RANGE
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=IstrU,Iend
            uee(i,0)=uee(i,1)
          enddo
        endif
        if (NORTHERN_EDGE) then
#   ifdef MPI
          do i=IstrU,Iend
            uee(i,Mmmpi+1)=uee(i,Mmmpi)
          enddo
#   else        
          do i=IstrU,Iend
            uee(i,Mm+1)=uee(i,Mm)
          enddo
#   endif          
        endif
#  endif
        do j=Jstr,Jend+1
          do i=IstrU-1,Iend
           Hvxx(i,j)=Dvom(i-1,j)-2.*Dvom(i,j)+Dvom(i+1,j)
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=IstrU,Iend
            cffX=urhs(i,j)+urhs(i,j-1)
            cffE=Dvom(i,j)+Dvom(i-1,j)
            if (cffE.gt.0.) then
              curvX=uee(i,j-1)
            else
              curvX=uee(i,j)
            endif
            UFe(i,j)=0.25*(cffX+gamma*curvX)*(cffE-0.125*(
     &                             Hvxx(i,j)+Hvxx(i-1,j) ))
          enddo
        enddo
#  undef Hvxx
#  undef uee
#  define vxx wrk1
#  define Huee wrk2

#  ifdef EW_PERIODIC 
#   define IV_EXT_RANGE Istr-1,Iend+1
#  else
#   ifdef MPI
        if (WEST_INTER) then
          imin=Istr-1
        else
          imin=max(Istr-1,1)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
#    define IV_EXT_RANGE imin,imax
#   else
#    define IV_EXT_RANGE max(Istr-1,1),min(Iend+1,Lm)
#   endif
#  endif

        do j=JstrV,Jend
          do i=IV_EXT_RANGE
            vxx(i,j)=vrhs(i-1,j)-2.*vrhs(i,j)+vrhs(i+1,j)
          enddo
        enddo
#   undef IV_EXT_RANGE
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=JstrV,Jend
            vxx(0,j)=vxx(1,j)
          enddo
        endif
        if (EASTERN_EDGE) then
#   ifdef MPI
          do j=JstrV,Jend
            vxx(Lmmpi+1,j)=vxx(Lmmpi,j)
          enddo
#   else        
          do j=JstrV,Jend
            vxx(Lm+1,j)=vxx(Lm,j)
          enddo
#   endif          
        endif
        
#  endif
        do j=JstrV-1,Jend
          do i=Istr,Iend+1
           Huee(i,j)=Duon(i,j-1)-2.*Duon(i,j)+Duon(i,j+1)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend+1
            cffE=vrhs(i,j)+vrhs(i-1,j)
            cffX=Duon(i,j)+Duon(i,j-1)
            if (cffX.gt.0.) then
              curvE=vxx(i-1,j)
            else
              curvE=vxx(i,j)
            endif
            VFx(i,j)=0.25*(cffE+gamma*curvE)*(cffX-0.125*(
     &                             Huee(i,j)+Huee(i,j-1) ))
          enddo
        enddo
#undef Huee
#undef vxx
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j)-UFx(i,j  )+UFx(i-1,j)
     &                           -UFe(i,j+1)+UFe(i  ,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j  )
     &                           -VFe(i  ,j)+VFe(i,j-1)
          enddo
        enddo
# endif /* !SOLVE3D && M2_HADV_UP3 */
#endif /* UV_ADV */
!
!-----------------------------------------------------------------------
! Compute Coriolis (2D and 3D) term and advective curvilinear metric
! terms (2D only).
!-----------------------------------------------------------------------
!
#if defined UV_COR || (defined CURVGRID && defined UV_ADV)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=Drhs(i,j)*(
# ifdef UV_COR
     &                   fomn(i,j)
# endif
# if (defined CURVGRID && defined UV_ADV)
     &          +0.5*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &                -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
# endif
     &                   )
# ifdef MRL_WCI
#  if defined CURVGRID && defined UV_ADV
          cff1 = Drhs(i,j)*(
     &    0.5*( dndx(i,j)*(vst2d(i,j)+vst2d(i,j+1))
     &          -dmde(i,j)*(ust2d(i,j)+ust2d(i+1,j)) ))
#  else
          cff1 = 0.0
#  endif
          UFx(i,j)=(cff+cff1)*(vrhs(i,j)+vrhs(i,j+1))
     &                 +cff*(vst2d(i,j)+vst2d(i,j+1))
          VFe(i,j)=(cff+cff1)*(urhs(i,j)+urhs(i+1,j))
     &                 +cff*(ust2d(i,j)+ust2d(i+1,j))
# else
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
# endif 
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar(i,j)=rubar(i,j)+0.25*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j) -0.25*(VFe(i,j)+VFe(i,j-1))
        enddo 
      enddo
#endif /* UV_COR */
!
!-----------------------------------------------------------------------
! Compute horizontal viscous stress terms (2D only).
!-----------------------------------------------------------------------
!
#ifdef SOLVE3D
# undef UV_VIS2
#endif
#ifdef UV_VIS2
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          cff=2.*Drhs(i,j)*visc2_r(i,j)
          UFx(i,j)=cff*(ubar(i+1,j,kstp)-ubar(i,j,kstp))
     &                                 *pm(i,j)*on_r(i,j)
          VFe(i,j)=cff*(vbar(i,j+1,kstp)-vbar(i,j,kstp))
     &                                 *pn(i,j)*om_r(i,j)

          cff1=0.0625*visc2_p(i+1,j+1)*( Drhs(i,j)
     &       +Drhs(i+1,j)+Drhs(i,j+1)+Drhs(i+1,j+1) )*(
     &          (pn(i+1,j+1)+pn(i,j+1)+pn(i+1,j)+pn(i,j))
     &             *(ubar(i+1,j+1,kstp)-ubar(i+1,j,kstp))
     &         +(pm(i+1,j+1)+pm(i,j+1)+pm(i+1,j)+pm(i,j))
     &             *(vbar(i+1,j+1,kstp)-vbar(i,j+1,kstp))
     &                                                  )
# ifdef MASKING
     &                     *pmask(i+1,j+1)
# endif
          UFe(i+1,j+1)=cff1*om_p(i+1,j+1)
          VFx(i+1,j+1)=cff1*on_p(i+1,j+1)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)+UFx(i,j)-UFx(i-1,j)
     &                         +UFe(i,j+1)-UFe(i,j)

          rvbar(i,j)=rvbar(i,j)+VFx(i+1,j)-VFx(i,j)
     &                         +VFe(i,j)-VFe(i,j-1)
        enddo
      enddo
#endif /* UV_VIS2 */
!
!-----------------------------------------------------------------------
! Linear and/or quadratic bottom stress.
!-----------------------------------------------------------------------
!
#ifdef BBL
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - bustr(i,j)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - bvstr(i,j)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
#else
      if (rdrg2.gt.0.) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
 
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
 
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
#endif
!
!-----------------------------------------------------------------------
! Add 2D vortex-force terms combined with advection terms
!-----------------------------------------------------------------------
!
#ifdef MRL_WCI
      do j=Jstr,Jend
        do i=IstrU,Iend
          vstu = 0.25*( vst2d(i,j) + vst2d(i,j+1)
     &                 +vst2d(i-1,j)+vst2d(i-1,j+1) )
          dudx = 0.5*( urhs(i+1,j)-urhs(i-1,j) )
          dvdx = 0.5*( vrhs(i,j) - vrhs(i-1,j)
     &                 +vrhs(i,j+1) - vrhs(i-1,j+1) )
          rubar(i,j) = rubar(i,j) + 0.5*om_u(i,j)*
     &                      ( Drhs(i-1,j)+Drhs(i,j) )
     &               *( ust2d(i,j)*dudx + vstu*dvdx )
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          ustv = 0.25*( ust2d(i,j) + ust2d(i+1,j)
     &                 +ust2d(i,j-1)+ust2d(i+1,j-1) )
          dude = 0.5*( urhs(i,j) - urhs(i,j-1)
     &                 +urhs(i+1,j) - urhs(i+1,j-1) )
ccc          dvde = 0.5*( vrhs(i,j-1)-vrhs(i,j+1) )
          dvde = 0.5*( vrhs(i,j+1)-vrhs(i,j-1) )
          rvbar(i,j) = rvbar(i,j) + 0.5*on_v(i,j)*
     &                      ( Drhs(i,j-1)+Drhs(i,j) )
     &               *( ustv*dude + vst2d(i,j)*dvde )
        enddo
      enddo
#endif

!
#ifdef SOLVE3D
!
      if (FIRST_2D_STEP) then
!
!***********************************************************************
! First 2D step: prepare forcing and coupling
! Compute rufrc & rvfrc
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
!***********************************************************************
!

# ifdef PRED_COUPLED_MODE
!
!-----------------------------------------------------------------------
!    Compute 3D coupling
!-----------------------------------------------------------------------
!
        if (FIRST_TIME_STEP) then
          cff3=0.                        ! This version is designed
          cff2=0.                        ! for coupling during 3D
          cff1=1.                        ! predictor sub-step: here
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
          cff3=0.                        ! computed as instantaneous
          cff2=-0.5                      ! value at 3D time step
          cff1=1.5                       ! "nstp" first, and then
        else                             ! extrapolated half-step
          cff3=0.281105                  ! forward using  AM3-like
          cff2=-0.5-2.*cff3              ! weights optimized for
          cff1=1.5+cff3                  ! maximum stability (with
        endif                            ! special care for startup)
        if (iic.eq.1) then
         do j=Jstr,Jend
          do i=IstrU,Iend
            cff=rufrc(i,j)
            rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                          + cff3*rufrc_bak(i,j,nstp)
            rufrc_bak(i,j,nstp)=cff
          enddo
         enddo
         do j=JstrV,Jend
          do i=Istr,Iend
            cff=rvfrc(i,j)
            rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                          + cff3*rvfrc_bak(i,j,nstp)
            rvfrc_bak(i,j,nstp)=cff
          enddo
         enddo
        else
         do j=Jstr,Jend
          do i=IstrU,Iend
            cff=rufrc(i,j)-rubar(i,j)
            rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                          + cff3*rufrc_bak(i,j,nstp)
            rufrc_bak(i,j,nstp)=cff
          enddo
         enddo
         do j=JstrV,Jend
          do i=Istr,Iend
            cff=rvfrc(i,j)-rvbar(i,j)
            rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                          + cff3*rvfrc_bak(i,j,nstp)
            rvfrc_bak(i,j,nstp)=cff
          enddo
         enddo
        endif

# else  /*  PRED_COUPLED_MODE */

        do j=Jstr,Jend                       ! This version is
          do i=Istr,Iend                     ! designed for coupling
            rufrc(i,j)=rufrc(i,j)-rubar(i,j) ! during 3D corrector
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j) ! sub-step: no forward
          enddo                              ! extrapolation needs
        enddo                                ! to be performed.

# endif /*  PRED_COUPLED_MODE */

!
      endif   !<-- FIRST_2D_STEP
!
#endif /* SOLVE3D */
!
!-----------------------------------------------------------------------
! Coupling EXT / NBQ
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
! Prepares Coupling with NBQ
!
!   ru_int_nbq       : RHS (3D) ( *mask & 2D correction)
!
!   ruext_nbq_2d     : RHS (2D)
!   ruext_nbq_2d_old : RHS (2D) at previous time-step
!   ruext_nbq_2d_sum : to compute time-averaged RHS (2D)
!   
!
!-----------------------------------------------------------------------
!
#if defined M2FILTER_NONE
!
      if (FIRST_2D_STEP) then
!
!------------------------------
! First 2D time-step only
!------------------------------
!
        do j=Jstr,Jend
          do i=IstrU,Iend
            ruext_nbq_2d_sum(i,j)=0.
	    ruext_nbq_2d_old(i,j)=0.
          enddo
        enddo

# ifdef MASKING
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)*umask(i,j)
            enddo
          enddo
        enddo
# endif
 
        do j=JstrV,Jend
          do i=Istr,Iend
            rvext_nbq_2d_sum(i,j)=0.
            rvext_nbq_2d_old(i,j)=0.
          enddo
        enddo 

# ifdef MASKING
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)*vmask(i,j)
            enddo
          enddo
        enddo
# endif

      endif
#endif /* M2FILTER_NONE */
!
!------------------------------
! All 2D time-steps 
!------------------------------
!

#define ruext_nbq_2d UFx
      do j=Jstr,Jend
        do i=IstrU,Iend
          ruext_nbq_2d(i,j)=(rufrc(i,j)+rubar(i,j))*pm_u(i,j)*pn_u(i,j)
     &
     &              /(Drhs(i,j)+Drhs(i-1,j))
#ifdef MASKING
          ruext_nbq_2d(i,j)=ruext_nbq_2d(i,j)*umask(i,j)
#endif
          ruext_nbq_2d_old(i,j)=ruext_nbq_2d(i,j)-ruext_nbq_2d_old(i,j)
#if defined M2FILTER_NONE     
          ruext_nbq_2d_sum(i,j)=ruext_nbq_2d_sum(i,j)+ruext_nbq_2d(i,j)
#endif          
        enddo
      enddo
      
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)= ru_int_nbq(i,j,k)+
     &        ruext_nbq_2d_old(i,j)
     &        *(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
     
      do j=Jstr,Jend
        do i=IstrU,Iend
          ruext_nbq_2d_old(i,j)=ruext_nbq_2d(i,j)
        enddo
      enddo
    
#undef ruext_nbq_2d

#define rvext_nbq_2d UFx
      do j=JstrV,Jend
        do i=Istr,Iend
          rvext_nbq_2d(i,j)=(rvfrc(i,j)+rvbar(i,j))*pm_v(i,j)*pn_v(i,j)
     &              /(Drhs(i,j)+Drhs(i,j-1))
#ifdef MASKING
          rvext_nbq_2d(i,j)=rvext_nbq_2d(i,j)*vmask(i,j)
#endif
          rvext_nbq_2d_old(i,j)=rvext_nbq_2d(i,j)-rvext_nbq_2d_old(i,j)
#if defined M2FILTER_NONE     
          rvext_nbq_2d_sum(i,j)=rvext_nbq_2d_sum(i,j)+rvext_nbq_2d(i,j)
#endif          
        enddo
      enddo
        
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)+
     &          rvext_nbq_2d_old(i,j)
     &         *(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          rvext_nbq_2d_old(i,j)=rvext_nbq_2d(i,j)
        enddo
      enddo
#undef rvext_nbq_2d 

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                                 ru_int_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                                 rv_int_nbq(START_2D_ARRAY,1))
#endif 

#ifdef RVTK_DEBUG
      call check_tab3d(ru_int_nbq,'ru_int_nbq (A)','u')
      call check_tab3d(rv_int_nbq,'rv_int_nbq (A)','v')
      call check_tab3d(rw_int_nbq,'rw_int_nbq (A)','wint')
#endif  

!
!
!***********************************************************************
! Prepare RHS of 3D NBQ system
!***********************************************************************
!

      if (FIRST_2D_STEP ) then	
!
!-----------------------------------------------------------------------
!  Some initializations at first 2D time-step:
!-----------------------------------------------------------------------
!
! Initializes Ru-avg2 variables
!     
#if defined M2FILTER_NONE 
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
          do k=0,N
            do j=Jstr,Jend
              do i=Istr,Iend
                rw_nbq_avg2(i,j,k)=qdmw_nbq(i,j,k)
              enddo
            enddo 
          enddo
#endif  /* M2FILTER_NONE */

      endif    ! FIRST_2D_STEP

!
!-----------------------------------------------------------------------
!  Indices & init
!-----------------------------------------------------------------------
!
      if (istr.eq.1) then
        IstrR2=Istr-1
#ifdef MPI
        IstrU2=Istr  !+1
#else
        IstrU2=Istr ! +1
#endif
      else
        IstrR2=Istr
        IstrU2=Istr
      endif

      if (iend.eq.Lm) then
        IendR2=Iend+1
      else
        IendR2=Iend
      endif

      if (jstr.eq.1) then
        JstrR2=Jstr-1
#ifdef MPI
        JstrV2=Jstr  !+1
#else
        JstrV2=Jstr !  +1
#endif
      else
        JstrR2=Jstr
        JstrV2=Jstr
      endif

      if (jend.eq.Mm) then
        JendR2=Jend+1
      else
        JendR2=Jend
      endif

#ifdef NBQ_DTDRHO2
       if (FIRST_TIME_STEP.and.iif==1) then
	  do j=JstrR2,JendR2
	    do k=1,N
	    do i=IstrR2,IendR2
	       zr_nbq(i,j,k,:)=z_r(i,j,k)
            enddo
            enddo
   	    do k=0,N
	    do i=IstrR2,IendR2
	       z_nbq (i,j,k,:)=z_w(i,j,k)
            enddo
            enddo
         enddo
       endif
#endif
!
!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
#if defined ACOUSTIC 
       if (FIRST_TIME_STEP.and.iif==1) then
          call densityijk_nbq(10)      
       endif
#endif

#undef DEBUG
!
!-------------------------------------------------------------------
!       Initialization of various test-cases
!-------------------------------------------------------------------
!       
!       if (iif==1.and.iic==1) call initial_nh_tile (3,Istr,Iend,Jstr,Jend)
!       if (iic==1.and.iif==1) thetadiv_nbq=0.
!
!-------------------------------------------------------------------
!  Get internal and external forcing terms for nbq equations:
!  ru+rubar (or lambda_ext+lambda_int)
!  dzdt*rhosurf
!-------------------------------------------------------------------
!
!
!------------------------------------------------------------------
!       Implicit part: system setup
!-------------------------------------------------------------------
!    
       do j=Jstr,Jend
         do i=Istr,Iend
	    WORK(i,j)=pm(i,j)*pn(i,j)
	 enddo
       enddo      

       rubar_nbq=0.
       rvbar_nbq=0.
       DU_nbq=0.
       DV_nbq=0.

#if defined NBQ_DTDRHO2 && defined NBQ_MASS
       if (FIRST_TIME_STEP.and.iif==1) then
        do k=1,N 
          do j=JstrR2,JendR2           
            do i=IstrR2,IendR2
	       rho_bak (i,j,k)=(   rho(i,j,k)/rho0)*Hzr(i,j,k)
	    enddo
	  enddo
	 enddo
	 z_nbq(:,:,:,1)=z_w(:,:,:)
	 z_nbq(:,:,:,2)=z_w(:,:,:)
	 z_nbq(:,:,:,3)=z_w(:,:,:)
	 z_nbq(:,:,:,4)=z_w(:,:,:)
       endif
#endif

!	
!------------------------------------------------------------------
!              Stores tendencies
!------------------------------------------------------------------
!
        do k=0,N 
          do j=Jstr,Jend             
            do i=Istr,Iend
               rw_nbq_ext (i,j,k) = qdmw_nbq(i,j,k) 
            enddo
          enddo
        enddo

#ifdef NBQ_ZETAEXP
!
!***********************************************************************
! Compute zeta and update grid.
!***********************************************************************
!
!
!-------------------------------------------------------------------
!       Computes surface mean velocities (Zeta)
!-------------------------------------------------------------------
!

! LAURENT: loop indices have to be corrected
       
        if (IstrU.le.Iend) then
         do j=Jstr,Jend
          do i=Istr,Iend+1     
               umean_nbq(i,j)=qdmu_nbq(i,j,N)  
# ifdef NBQ_MASS                      
     &            / (Hz(i,j,N)+Hz(i-1,j,N)) * 2. 
# else
     &            / (Hzr(i,j,N)+Hzr(i-1,j,N)) * 2. 
# endif
# ifdef MASKING
     &            * umask(i,j) 
# endif
          enddo 
         enddo 
        endif

        if (JstrV.le.Jend) then
         do j=Jstr,Jend+1
          do i=Istr,Iend     
               vmean_nbq(i,j)=qdmv_nbq(i,j,N)              
# ifdef NBQ_MASS     
     &            / (Hz(i,j,N)+Hz(i,j-1,N)) * 2. 
# else
     &            / (Hzr(i,j,N)+Hzr(i,j-1,N)) * 2. 
# endif
# ifdef MASKING
     &            * vmask(i,j) 
# endif
          enddo
         enddo 
        endif

        do j=Jstr,Jend
          do i=Istr,Iend
               wmean_nbq(i,j,kstp)=qdmw_nbq(i,j,N)         
# ifdef NBQ_MASS     
!    &     / (rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)+1.+rho(i,j,N)/rho0)  ! Surface BC: Rho_nbq(surf)=0 
     &     / (1.+rho(i,j,N)/rho0)  
# endif
     &             * Hzw_half_nbq_inv(i,j,N)   
# ifdef MASKING
     &             * rmask(i,j) 
# endif
          enddo
        enddo 

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
        if (IstrU.le.Iend) then
         call exchange_u2d_tile (Istr,Iend,Jstr,Jend,umean_nbq(START_2D_ARRAY))
        endif
        if (JstrV.le.Jend) then
         call exchange_v2d_tile (Istr,Iend,Jstr,Jend,vmean_nbq(START_2D_ARRAY))
        endif
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,wmean_nbq(START_2D_ARRAY,kstp))
# endif

!
!-------------------------------------------------------------------
c After this call to step2d_zetaw
c the vertical grid (z_w, z_r, Hzr, Hzw_half_nbq) will be in (n+1)
c while the derived arrays (Hzr_half_nbq_inv, Hzw_half_nbq_inv .. ... )
c are still in (n)
c since the call to grid_coef_nh is at the beginning of step2d
!-------------------------------------------------------------------
!

!
!-------------------------------------------------------------------
! Zeta is computed (first guess) and grid is updated (n+1)
!-------------------------------------------------------------------
!
# define zab3 UFx
!
!-----------------------------------------------------------------------
! Computes zeta(n+1)
!-----------------------------------------------------------------------
!
      if (FIRST_TIME_STEP.and.iif.eq.1) then
         cff8=1.D0
         cff9=0.D0
         cff10=0.D0
       elseif (FIRST_TIME_STEP.and.iif.eq.2) then
         cff8=1.5D0
         cff9=-0.5D0
        cff10=0.D0
      else
        cff8= 1.5D0+mybeta
        cff9=-2.0D0*mybeta-0.5D0
        cff10= mybeta
      endif

      do j=JstrV-2,Jend+1
        do i=IstrU-2,Iend+1
           zab3(i,j) =   cff8  * zeta(i,j,kstp)
     &                 + cff9  * zeta(i,j,kbak2)
     &                 + cff10 * zeta(i,j,kold2)
        enddo
      enddo

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
         zeta (i,j,knew)=
     &    ( zeta(i,j,kstp)
     &        +dtfast*(
# if defined NBQ_ZETAEXP
     &    +wmean_nbq(i,j,kstp)
# else
     &    +wmean_nbq(i,j,knew)
# endif
     &    -0.5*(umean_nbq(i  ,j)
     &          *(zab3(i  ,j)
     &           -zab3(i-1,j))*pm_u(i,j)
# ifdef MASKING
     &         *rmask(i,j)*rmask(i-1,j)
# endif
     &         +umean_nbq(i+1,j)
     &          *(zab3(i+1,j)
     &           -zab3(i  ,j))*pm_u(i+1,j) 
# ifdef MASKING
     &         *rmask(i,j)*rmask(i+1,j)
# endif
     &          )
# ifdef MASKING
     &         *umask(i,j)*umask(i+1,j)
# endif

     &    -0.5*(vmean_nbq(i  ,j)
     &          *(zab3(i,j  ) 
     &           -zab3(i,j-1))*pm_v(i,j)
# ifdef MASKING
     &         *rmask(i,j)*rmask(i,j-1)
# endif
     &         +vmean_nbq(i,j+1)
     &          *(zab3(i,j+1)
     &           -zab3(i,j  ))*pm_v(i,j+1) 
# ifdef MASKING
     &         *rmask(i,j)*rmask(i,j+1)
# endif
     &         )
# ifdef MASKING
     &         *vmask(i,j)*vmask(i,j+1)
# endif
     &         ) )
# ifdef MASKING
     &         *rmask(i,j)
# endif
        enddo
      enddo

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
# endif

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
           zeta_new(i,j)=zeta(i,j,knew)
        enddo
      enddo
!
!-----------------------------------------------------------------------
! Add nudging terms
!-----------------------------------------------------------------------
!
# ifdef ZNUDGING
#   ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew),zetazon)
        endif
      endif
      if (iic.eq.ntstart) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ssh(START_2D_ARRAY),sshzon)
        endif
      endif
#   endif  /* ZONAL_NUDGING */
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j) + dtfast*Znudgcof(i,j)
#   ifdef ZONAL_NUDGING
     &                                 *(sshzon(j)-zetazon(j))
#   else
     &                                 *(ssh(i,j)-zeta_new(i,j))
#   endif /* ZONAL_NUDGING */
        enddo
      enddo
# endif /* ZNUDGING */
!
!-----------------------------------------------------------------------
! Computes zetarhs to use in momentum equations
!-----------------------------------------------------------------------
!

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
 	   zeta_new(i,j)=zeta_new(i,j) SWITCH rmask(i,j)          
# if !defined NBQ_ZETAEXP
!  in case of ZETAEXP rhobar has not yet been computed
!  Dnew will be computed later
          Dnew(i,j)=(zeta_new(i,j)+h(i,j))                            ! USEFUL ????
#  if defined NBQ_MASS
     &       *rhobar_nbq(i,j,knew)
#  endif
# endif
        enddo
      enddo

!
!-----------------------------------------------------------------------
! Load new free-surface values into shared array
! Modify new free-surface to ensure that depth is > Dcrit for masked
! cells.
!-----------------------------------------------------------------------
!

# if defined WET_DRY && defined MASKING
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   (Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
        enddo
      enddo 
# endif
!
!-----------------------------------------------------------------------
! Load rhs values into additional AGRIF shared array for nesting
!-----------------------------------------------------------------------
! 
# ifdef AGRIF
      if (FIRST_2D_STEP) then
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,0)=zeta(i,j,kstp)       
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
          du_avg3(i,j,0)  = DUon(i,j)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
          dv_avg3(i,j,0)  = DVom(i,j)
          enddo
        enddo 
      endif

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   Zt_avg3(START_2D_ARRAY,0))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   du_avg3(START_2D_ARRAY,0))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   dv_avg3(START_2D_ARRAY,0))
#  endif

#  ifdef RVTK_DEBUG_ADVANCED
       if (.not.agrif_Root()) then
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(Zt_avg3(:,:,0),'Zt_avg3 (index 0) step2d','r')
       call check_tab2d(DU_avg3(:,:,0),'DU_avg3 (index 0) step2d','u')
       call check_tab2d(DV_avg3(:,:,0),'DV_avg3 (index 0) step2d','v')
C$OMP END MASTER  
       endif
#   endif  
# endif /* AGRIF */   
!
!-----------------------------------------------------------------------
! Compute wet/dry masks
!-----------------------------------------------------------------------
!
# ifdef WET_DRY
      call wetdry_tile (Istr,Iend,Jstr,Jend)
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
# endif
!
!-----------------------------------------------------------------------
! Set boundary conditions for the free-surface
!-----------------------------------------------------------------------
!
      call zetabc_tile (Istr,Iend,Jstr,Jend)

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
# ifdef SOLVE3D
       if ((iic.eq.1.and.iif==1)
     &     NSTEP_GRID
     &    ) then
!
!-----------------------------------------------------------------------
! Update Grid:
!-----------------------------------------------------------------------
!
        call set_depth_tile(Istr,Iend,Jstr,Jend
     &   ,resetfromrhobar
     &   ) 

!
!-----------------------------------------------------------------------
!     Update grid variables
!-----------------------------------------------------------------------
!
       call grid_coef_nh(
     &   Istr,Iend,Jstr,Jend,
     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &   Hzu_half_qdmu, Hzv_half_qdmv                                     
     &   )

      endif

#  undef zab3

# endif /* SOLVE3D */


#endif  /* GRIDEXT */

!
!***********************************************************************
! NH-Momemtum equations
!***********************************************************************
!
!-------------------------------------------------------------------
!       "Pressure - Viscosity" Variable (theta)
!               theta does not change
!-------------------------------------------------------------------
!      
        do k=1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              thetadiv_nbq(i,j,k)=(-visc2_nbq*(thetadiv_nbq(i,j,k)
     &                                        +thetadiv3_nbq(i,j,k))
     &                             +soundspeed2_nbq(i,j)*rho_nbq(i,j,k)) 
     &                               *Hzr_half_nbq_inv(i,j,k)  
            enddo
          enddo
        enddo
!
!-------------------------------------------------------------------
!      Horizontal Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!-------------------------------------------------------------------
!
!---------------------------
!  XI- and ETA-Directions:
!---------------------------
!
        k2 = 1
        do k=0,N
          k1=k2
	  k2=3-k1

#ifdef NBQ_NODS
          if (NSTEP_DS) then
#endif
          if (k.eq.0) then ! Bottom Boundary conditions
	
	    do j=Jstr,Jend
              do i=IstrU2,Iend
	        dthetadiv_nbqdz_u(i,j,k2)=0. 
	      enddo
	    enddo

            do j=JstrV2,Jend
              do i=Istr,Iend
	        dthetadiv_nbqdz_v(i,j,k2)=0.
	      enddo
    	    enddo 

          else

            if (k.eq.N) then ! Top Boundary conditions

              do j=JstrV2-1,Jend
                do i=IstrU2-1,Iend
#ifndef NBQ_NODS
                  dthetadiv_nbqdz(i,j)    = - thetadiv_nbq(i  ,j,k)
#else
                  dthetadiv_nbqdz(i,j,k,1)= - thetadiv_nbq(i  ,j,k)
#endif
     	        enddo
	      enddo

            else

              do j=JstrV2-1,Jend
                do i=IstrU2-1,Iend
#ifndef NBQ_NODS
                  dthetadiv_nbqdz(i,j)    =thetadiv_nbq(i  ,j,k+1) 
     &             - thetadiv_nbq(i  ,j,k)
#else
                  dthetadiv_nbqdz(i,j,k,1)=thetadiv_nbq(i  ,j,k+1)
     &              - thetadiv_nbq(i  ,j,k)
#endif
                enddo
              enddo

            endif
  
            do j=Jstr,Jend
            do i=IstrU2,Iend
#ifndef NBQ_NODS
              dthetadiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &            dthetadiv_nbqdz(i,j)
     &           +dthetadiv_nbqdz(i-1,j)) 
#else
              dthetadiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &           dthetadiv_nbqdz(i,j,k,1)
     &          +dthetadiv_nbqdz(i-1,j,k,1))    
#endif          
            enddo
            enddo
            do j=JstrV2,Jend
              do i=Istr,Iend
#ifndef NBQ_NODS
                dthetadiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(
     &         dthetadiv_nbqdz(i,j)
     &        +dthetadiv_nbqdz(i,j-1))
#else
                dthetadiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(
     &       dthetadiv_nbqdz(i,j,k,1)
     &      +dthetadiv_nbqdz(i,j-1,k,1))
#endif
              enddo
            enddo    
  
            endif    
#ifdef NBQ_NODS
          endif
#endif          

          if (k.gt.0) then
!
!---------------------------
!...........U-momentum:
!---------------------------
!
            do j=Jstr,Jend
              do i=IstrU2,Iend
                if (k.gt.1.and.k.lt.N) then 
#ifndef NBQ_NODS
	  	  dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                  *(dthetadiv_nbqdz_u(i,j,k2)+dthetadiv_nbqdz_u(i,j,k1))             ! dZdx * (d(delta p)dz)_u
#else
                  if (NSTEP_DS) then
	   	  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)
     &                    -z_r(i-1,j,k))       
     &                  *(dthetadiv_nbqdz_u(i,j,k2)+dthetadiv_nbqdz_u(i,j,k1))              ! dZdx * (d(delta p)dz)_u
                  endif 
                  dum_s=dthetadiv_nbqdz(i,j,k,1)                                        
#endif
                elseif (k.gt.1) then
#ifndef NBQ_NODS
	 	  dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                        *dthetadiv_nbqdz_u(i,j,k1)                                        ! dZdx * (d(delta p)dz)_u
     &                    +(z_w(i,j,N)-z_w(i-1,j,N))                      
     &                  *dthetadiv_nbqdz_u(i,j,k2)
#else
                  if (NSTEP_DS) then
	 	  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)
     &        -z_r(i-1,j,k))       
     &                   *dthetadiv_nbqdz_u(i,j,k1)                                       ! dZdx * (d(delta p)dz)_u
     &                  +(z_w(i,j,N)-z_w(i-1,j,N))                      
     &                  *dthetadiv_nbqdz_u(i,j,k2)                 
                  endif
	 	  dum_s=dthetadiv_nbqdz(i,j,k,1)                                          
#endif
                else
#ifndef NBQ_NODS
	  	  dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                  *2.*dthetadiv_nbqdz_u(i,j,k2)                                      ! dZdx * (d(delta p)dz)_u
#else
                  if (NSTEP_DS) then
  	  	  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)
     &                     -z_r(i-1,j,k))       
     &                    *2.*dthetadiv_nbqdz_u(i,j,k2)                                    ! dZdx * (d(delta p)dz)_u
                  endif
	  	  dum_s=dthetadiv_nbqdz(i,j,k,1)                                          
#endif
                endif
#ifdef MASKING
                if (umask(i-1,j)*umask(i+1,j)==0.) then
                dum_s=dum_s
     &   -(thetadiv_nbq(i,j,k)
     &   -thetadiv_nbq(i-1,j,k))                       ! - d(delta p)dx
                else
#endif
                dum_s=dum_s
     &   -(gammau*thetadiv_nbq(i,j,k)+gammau_2*thetadiv_nbq(i+1,j,k)
     &   -gammau*thetadiv_nbq(i-1,j,k)-gammau_2*thetadiv_nbq(i-2,j,k))                       ! - d(delta p)dx
#ifdef MASKING
                endif
#endif
                dum_s=dum_s   *Hzu_half_qdmu(i,j,k)
                qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k) + dtnbq * (
     &              dum_s + ru_int_nbq(i,j,k))  
                DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
                ru_nbq_ext (i,j,k) = dum_s / work(i,j) 
                rubar_nbq(i,j)=rubar_nbq(i,j)+ru_nbq_ext(i,j,k)
#ifdef NBQ_NUDGING
             qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*(1.-nudg_coef_nbq(i,j))
     &        +u(i,j,k,nrhs)*hzu_half_qdmu(i,j,k)
     &        * nudg_coef_nbq(i,j)
#endif

              enddo 
            enddo

!
!---------------------------
!...........V-momentum:
!---------------------------
!
            do j=JstrV2,Jend
              do i=Istr,Iend
                if (k.gt.1.and.k.lt.N) then 
#ifndef NBQ_NODS
	          dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                     *(dthetadiv_nbqdz_v(i,j,k2)
     &                  +dthetadiv_nbqdz_v(i,j,k1))    ! dZdy * (d(delta p)dz)_v
#else
                  if (NSTEP_DS) then
	          dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)
     &    -z_r(i,j-1,k)) 
     &                  *(dthetadiv_nbqdz_v(i,j,k2)+dthetadiv_nbqdz_v(i,j,k1))    ! dZdy * (d(delta p)dz)_v
                  endif
	          dum_s=dthetadiv_nbqdz(i,j,k,2)                         
#endif
                elseif (k.gt.1) then
#ifndef NBQ_NODS
	          dum_s=(z_r(i,j,k)-z_r(i,j-1,k))            
     &                  *dthetadiv_nbqdz_v(i,j,k1)                                ! dZdy * (d(delta p)dz)_v
     &                   +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                  *dthetadiv_nbqdz_v(i,j,k2)
#else
                  if (NSTEP_DS) then
     	          dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)
     &                    -z_r(i,j-1,k)) 
     &                   *dthetadiv_nbqdz_v(i,j,k1)                               ! dZdy * (d(delta p)dz)_v
     &                  +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                  *dthetadiv_nbqdz_v(i,j,k2)
                  endif
	          dum_s=dthetadiv_nbqdz(i,j,k,2) 
#endif
                else
#ifndef NBQ_NODS
	          dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                  *2.*dthetadiv_nbqdz_v(i,j,k2)                             ! dZdy * (d(delta p)dz)_v
#else
                  if (NSTEP_DS) then
	          dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)
     &                      -z_r(i,j-1,k)) 
     &                  *2.*dthetadiv_nbqdz_v(i,j,k2)                             ! dZdy * (d(delta p)dz)_v
                  endif
	          dum_s=dthetadiv_nbqdz(i,j,k,2) 
#endif
                endif
                
#ifdef MASKING
                if (vmask(i,j-1)*vmask(i,j+1)==0.) then
                dum_s=dum_s
     &   -(thetadiv_nbq(i,j,k)
     &   -thetadiv_nbq(i,j-1,k))                                                  ! - d(delta p)dy
                else
#endif
                dum_s=dum_s
     &   -(gammau*thetadiv_nbq(i,j,k)+gammau_2*thetadiv_nbq(i,j+1,k)
     &   -gammau*thetadiv_nbq(i,j-1,k)-gammau_2*thetadiv_nbq(i,j-2,k))            ! - d(delta p)dy
#ifdef MASKING
                endif
#endif
                dum_s=dum_s*Hzv_half_qdmv(i,j,k)
                qdmv_nbq(i,j,k) = qdmv_nbq(i,j,k) + dtnbq * (
     &                  dum_s + rv_int_nbq(i,j,k))    
                DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
                rv_nbq_ext (i,j,k) = dum_s / work(i,j)  
                rvbar_nbq(i,j)=rvbar_nbq(i,j)+rv_nbq_ext(i,j,k)	

#ifdef NBQ_NUDGING
             qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*(1.-nudg_coef_nbq(i,j))
     &        +v(i,j,k,nrhs)*hzv_half_qdmv(i,j,k)
     &        * nudg_coef_nbq(i,j)
#endif

              enddo
            enddo
          endif
        enddo        
        
!
!---------------------------
!  U-momentum open boundary conditions
!---------------------------
#ifdef RVTK_DEBUG
      call check_tab3d(qdmu_nbq,'qdmu_nbqint','uint')
      call check_tab3d(qdmv_nbq,'qdmv_nbqint','vint')
#endif

#ifndef NBQ_GRIDEXT    
       call u2dbc_tile      (Istr,Iend,Jstr,Jend, UFx)
       call v2dbc_tile      (Istr,Iend,Jstr,Jend, UFx)
#endif
       
#ifdef OBC_NBQ
       call unbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
       call vnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
#endif


!--------------------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!--------------------------------------------------------------------
!

#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,qdmu_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,qdmv_nbq(START_2D_ARRAY,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,DU_nbq(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,DV_nbq(START_2D_ARRAY))
#endif
#ifdef RVTK_DEBUG
      call check_tab3d(qdmu_nbq,'qdmu_nbq','u')
      call check_tab3d(qdmv_nbq,'qdmv_nbq','v')
#endif  
!-------------------------------------------------------------------
!      Explicit Vertical Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!         If implicit: (x,y)   only
!-------------------------------------------------------------------
!
#ifndef NBQ_IMP
!---------------------------
!  Z-Direction: Explicit
!---------------------------
        do j=Jstr,Jend
          do k=1,N-1
            do i=Istr,Iend                                                               
               dum_s =   thetadiv_nbq(i,j,k) - thetadiv_nbq(i,j,k+1)   
               qdmw_nbq(i,j,k)   = qdmw_nbq(i,j,k)   
     &            + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
# ifdef NBQ_GRAV
     &            -0.25*(rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)
     &                 +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
     &                *(Hzr(i,j,k)+Hzr(i,j,k+1))
     &                 *g*dtnbq
# endif
# ifdef MASKING
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
# endif            
# ifdef NBQ_NUDGING
             qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)*(1.-nudg_coef_nbq(i,j)) &
              +wz(i,j,k,nrhs)*hzw_half_nbq(i,j,k) &
              * nudg_coef_nbq(i,j)
# endif      
            enddo             
          enddo
          k=N
          do i=Istr,Iend                                                               
               dum_s =   thetadiv_nbq(i,j,N)                              
               qdmw_nbq(i,j,N) = qdmw_nbq(i,j,N)   
     &           + dtnbq * ( dum_s + rw_int_nbq(i,j,N) )
# ifdef NBQ_GRAV
     &            -rho_nbq(i,j,N)*0.5
     &                 *g*dtnbq
# endif
# ifdef MASKING
               qdmw_nbq(i,j,N) = qdmw_nbq(i,j,N) * rmask(i,j) 
# endif           
# ifdef NBQ_NUDGING
             qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)*(1.-nudg_coef_nbq(i,j)) &
              +wz(i,j,N,nrhs)*hzw_half_nbq(i,j,N) &
              * nudg_coef_nbq(i,j)
# endif    
          enddo     
        		   
! Bottom boundary:      

         enddo

!---------------------------
! Vertical momentum open boundary conditions
!---------------------------
# ifdef OBC_NBQ
        call wnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,qdmw_nbq(START_2D_ARRAY,0))
# endif
#endif  /* NBQ_IMP */

!
!-------------------------------------------------------------------
!      Mass equation (1): DX(p+1) ==> thetadiv_nbq
!-------------------------------------------------------------------
!

!--------------------------- 
! X -component 
!---------------------------

        k2 = 1
        do k=0,N
          k1=k2
	  k2=3-k1

#ifdef NBQ_NODS
          if (NSTEP_DS) then
#endif

          if (k.lt.N) then
             kp1 = k + 1
	     do j=Jstr,Jend
             do i=Istr,Iend+1
	       dZdxq_u(i,j,k2)=(z_r(i,j,kp1)-z_r(i-1,j,kp1)) 
     &                  *qdmu_nbq(i,j,kp1)    ! (dZdx * (rho u))_u
             enddo
             enddo
          endif

	  if (k.eq.0) then	! Bottom boundary conditions

#if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend
	    do i=Istr,Iend+1 
# ifndef NBQ_NODS
              dZdxq_w(i,j,k2)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                          *qdmu_nbq(i,j,1)  
     &                          / (Hzr(i,j,1)+Hzr(i-1,j,1))
# else
              dZdxq_w(i,j,k) = (z_w(i,j,0)-z_w(i-1,j,0))
     &                        *qdmu_nbq(i,j,1)  
     &                          / (Hzr(i,j,1)+Hzr(i-1,j,1))
# endif
	    enddo
	    enddo

 	    do j=Jstr,Jend
 	    do i=Istr,Iend    
# ifndef NBQ_NODS     
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k2) *pm_u(i,j) 
     &                +dZdxq_w(i+1,j,k2) *pm_u(i+1,j) ) 
     &                        * Hzr(i,j,1)      
# else
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k) *pm_u(i,j) 
     &                         +dZdxq_w(i+1,j,k) *pm_u(i+1,j) ) 
     &                        * Hzr(i,j,1)      
# endif

# if defined MASKING
              qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
# endif 
 	    enddo
 	    enddo 

#else /* NBQ_FREESLIP */

# ifndef NBQ_NODS
 	    do j=Jstr,Jend
 	    do i=Istr,Iend +1
              dZdxq_w(i,j,k2)=0.  
              qdmw_nbq(i,j,0)=0.     
 	    enddo
 	    enddo 
# endif
	  
#endif /* NBQ_FREESLIP */

          elseif (k==N) then ! Top boundary conditions
           
            do j=Jstr,Jend
	    do i=Istr,Iend+1
#ifndef NBQ_NODS
	      dZdxq_w(i,j,k2)= (z_w(i,j,N)-z_w(i-1,j,N))   
     &                 *qdmu_nbq(i,j,N)                                                            
     &                 / (Hzr(i,j,N)+Hzr(i-1,j,N)) 
#else
	      dZdxq_w(i,j,k)= (z_w(i,j,N)-z_w(i-1,j,N))   
     &                 *qdmu_nbq(i,j,N)                                                            
     &                 / (Hzr(i,j,N)+Hzr(i-1,j,N))
#endif  
            enddo
            enddo  
 
           else ! k<>0 & k<>N   ! Inner domain   

            do j=Jstr,Jend
	    do i=Istr,Iend+1
#ifndef NBQ_NODS
	       dZdxq_w(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &           dZdxq_u(i,j,k1)+dZdxq_u(i,j,k2)) 
#else
	       dZdxq_w(i,j,k )=Hzw_half_nbq_inv_u(i,j,k)*(
     &           dZdxq_u(i,j,k1)+dZdxq_u(i,j,k2)) 
#endif
            enddo 
            enddo

	   endif ! k<>0 , k<>N ,  Inner domain 

#ifdef NBQ_NODS

	  else   ! Update d./ds terms

	  if (k.eq.0) then	! Bottom boundary conditions

# if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend
	    do i=Istr,Iend+1 
              dZdxq_w(i,j,k)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                           * qdmu_nbq(i,j,1)  
     &                           / (Hzr(i,j,1)+Hzr(i-1,j,1))
	    enddo
	    enddo

 	    do j=Jstr,Jend
 	    do i=Istr,Iend    
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k) *pm_u(i,j) +dZdxq_w(i+1,j,k) 
     &                         * pm_u(i+1,j) ) 
     &                         * Hzr(i,j,1)    

#  if defined MASKING
              qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#  endif 
 	    enddo
 	    enddo 

# endif /* NBQ_FREESLIP */

          endif ! k.eq.0
          endif ! Update d./ds terms

#else  /* NBQ_NODS */

	  if (k.eq.0) then	! Bottom boundary conditions

# if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend
	    do i=Istr,Iend+1 
               dZdxq_w(i,j,k2)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                           * qdmu_nbq(i,j,1)  
     &                           / (Hzr(i,j,1)+Hzr(i-1,j,1))
	    enddo
	    enddo

 	    do j=Jstr,Jend
 	    do i=Istr,Iend    
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k2) *pm_u(i,j) +dZdxq_w(i+1,j,k2) 
     &                         * pm_u(i+1,j) ) 
     &                         * Hzr(i,j,1)    

#  if defined MASKING
              qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#  endif 
 	    enddo
 	    enddo 

# endif /* NBQ_FREESLIP */ 

          endif ! k.eq.0

#endif  /* NBQ_NODS */

          if (k.gt.0) then
           if (IstrU.le.Iend) then
            do j=Jstr,Jend
	    do i=Istr,Iend+1
#ifndef NBQ_NODS
	      FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k2)-dZdxq_w(i,j,k1))
#else
	      FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k)-dZdxq_w(i,j,k-1))
#endif

#ifdef MASKING
              FX(i,j) = FX(i,j) * umask(i,j)
#endif                
            enddo
            enddo

            do j=Jstr,Jend
            do i=Istr,Iend
	      thetadiv_nbq(i,j,k)=FX(i,j)  +FX(i+1,j)           
            enddo
            enddo
           endif 
          endif
	enddo	
	 
!---------------------------
! Y component     
!---------------------------   

        k2 = 1
	do k=0,N
	  k1=k2
	  k2=3-k1

#ifdef NBQ_NODS
          if (NSTEP_DS) then
#endif

	  if (k.lt.N) then
            kp1 = k + 1
            do j=Jstr,Jend+1
            do i=Istr,Iend
	      dZdyq_v(i,j,k2)=(z_r(i,j,kp1)-z_r(i,j-1,kp1)) 
     &                   *qdmv_nbq(i,j,kp1)    ! (dZdy * (rho v))_v
            enddo
	    enddo			 
          endif

          if (k.eq.0) then	! Bottom boundary conditions

#if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend+1
            do i=Istr,Iend
# ifndef NBQ_NODS
               dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                     *qdmv_nbq(i,j,1) 
     &                     / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
# else
               dZdyq_w(i,j,k)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                      *qdmv_nbq(i,j,1) 
     &                      / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
# endif 
	    enddo
	    enddo
 	    do j=Jstr,Jend
  	    do i=Istr,Iend   
# ifndef NBQ_NODS     
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k2)*pm_v(i,j)  
     &    +dZdyq_w(i,j+1,k2)*pm_v(i,j+1)  )    
     &                               * Hzr(i,j,1) 
# else
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k)*pm_v(i,j)  
     &       +dZdyq_w(i,j+1,k)*pm_v(i,j+1)  )    
     &                               * Hzr(i,j,1) 
# endif 
# if defined MASKING
                 qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
# endif 
  	    enddo
            enddo
#else  /* NBQ_FREESLIP */
# ifndef NBQ_NODS     
 	    do j=Jstr,Jend +1
 	    do i=Istr,Iend  
              dZdyq_w(i,j,k2)=0.  
              qdmw_nbq(i,j,0)=0.
 	    enddo
 	    enddo 
# endif
#endif /* NBQ_FREESLIP */

          elseif (k==N) then ! Top boundary conditions

            do j=Jstr,Jend+1
	    do i=Istr,Iend
#ifndef NBQ_NODS
              dZdyq_w(i,j,k2)= (z_w(i,j,N)-z_w(i,j-1,N))       
     &                        * qdmv_nbq(i,j,N)                                 
     &                        / ( Hzr(i,j,N)+Hzr(i,j-1,N) )
#else
              dZdyq_w(i,j,k)= (z_w(i,j,N)-z_w(i,j-1,N))       
     &                        * qdmv_nbq(i,j,N)                                 
     &                        / ( Hzr(i,j,N)+Hzr(i,j-1,N) )
#endif 
	    enddo
	    enddo

          else

      	    do j=Jstr,Jend+1
            do i=Istr,Iend
#ifndef NBQ_NODS
              dZdyq_w(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)
     &    *(dZdyq_v(i,j,k1)+dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
#else
              dZdyq_w(i,j,k)=Hzw_half_nbq_inv_v(i,j,k)
     &     *(dZdyq_v(i,j,k1)+dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
#endif 
            enddo 
            enddo

          endif

#ifdef NBQ_NODS

          else

           if (k.eq.0) then	! Bottom boundary conditions

# if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend+1
            do i=Istr,Iend
               dZdyq_w(i,j,k)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                   *qdmv_nbq(i,j,1) 
     &                   / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
	    enddo
	    enddo
 	    do j=Jstr,Jend
  	    do i=Istr,Iend   
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k)*pm_v(i,j)  
     &                     +dZdyq_w(i,j+1,k)*pm_v(i,j+1)  )     
     &                               * Hzr(i,j,1) 
#  if defined MASKING
                 qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#  endif 
  	    enddo
            enddo
# endif
           endif
          endif

#else /* NBQ_NODS */

          if (k.eq.0) then	! Bottom boundary conditions

# if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend+1
            do i=Istr,Iend
               dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                   *qdmv_nbq(i,j,1) 
     &                   / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
	    enddo
	    enddo
 	    do j=Jstr,Jend
  	    do i=Istr,Iend   
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k2)*pm_v(i,j)  
     &                     +dZdyq_w(i,j+1,k2)*pm_v(i,j+1)  )     
     &                               * Hzr(i,j,1) 
#  if defined MASKING
                 qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#  endif 
  	    enddo
            enddo
# endif  /* NBQ_FREESLIP */
          endif

#endif /* NBQ_NODS */

          if (k.gt.0) then
           if (JstrV.le.Jend) then
            do j=Jstr,Jend+1
            do i=Istr,Iend 
#ifndef NBQ_NODS
	      FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k2)-dZdyq_w(i,j,k1))
#else
	      FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k)-dZdyq_w(i,j,k-1))
#endif 
#ifdef MASKING
              FY(i,j) = FY(i,j) * vmask(i,j)
#endif                 
            enddo
            enddo
	    do j=Jstr,Jend
	    do i=Istr,Iend
   	      thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)+FY(i,j)+FY(i,j+1)              
            enddo
            enddo
           endif
          endif
	enddo		         
        do k=1,N

         if (IstrU.le.Iend) then
          do j=Jstr,Jend		   
          do i=Istr,Iend+1
            FX(i,j)=on_u(i,j)* qdmu_nbq(i,j,k)
#ifdef MASKING
            FX(i,j) = FX(i,j) * umask(i,j)
#endif
          enddo
          enddo	
         endif
			  
         if (JstrV.le.Jend) then
          do j=Jstr,Jend+1		   
          do i=Istr,Iend
            FY(i,j)=om_v(i,j)* qdmv_nbq(i,j,k)
#ifdef MASKING
            FY(i,j) = FY(i,j) * vmask(i,j)
#endif
          enddo
          enddo
         endif

!---------------------------
! Horizontal Divergence :
!     Dx(p+1)     
!---------------------------   
         if (IstrU.gt.Iend) then
          do j=Jstr,Jend		   
          do i=Istr,Iend			   
            thetadiv_nbq(i,j,k)=(thetadiv_nbq(i,j,k)                         
     &		            +WORK(i,j)*(FY(i,j+1)-FY(i,j))  
     &		                       ) 
#ifdef MASKING
            thetadiv_nbq(i,j,k) = thetadiv_nbq(i,j,k) * rmask(i,j)
#endif                              
          enddo
          enddo
         elseif (JstrV.gt.Jend) then
          do j=Jstr,Jend		   
          do i=Istr,Iend			   
            thetadiv_nbq(i,j,k)=(thetadiv_nbq(i,j,k)                         
     &		            +WORK(i,j)*(FX(i+1,j)-FX(i,j))  
     &		                       ) 
#ifdef MASKING
            thetadiv_nbq(i,j,k) = thetadiv_nbq(i,j,k) * rmask(i,j)
#endif                              
          enddo
          enddo
         else
          do j=Jstr,Jend		   
          do i=Istr,Iend			   
            thetadiv_nbq(i,j,k)=(thetadiv_nbq(i,j,k)                         
     &		            +WORK(i,j)*(FX(i+1,j)-FX(i,j)+FY(i,j+1)-FY(i,j))  
     &		                       ) 
#ifdef MASKING
            thetadiv_nbq(i,j,k) = thetadiv_nbq(i,j,k) * rmask(i,j)
#endif                              
          enddo
          enddo
         endif  
        enddo

!---------------------------
! Time and Bp density variations
!---------------------------  
 
#if defined NBQ_DTDRHO2 
	    z_nbq(:,:,:,knew)=z_w(:,:,:)
#endif

#if defined NBQ_DTDRHO2 

# ifndef NBQ_GRIDEXT
          if (iif==1) then
# endif
        do j=Jstr,Jend
         do i=Istr,Iend
           FC(i,0)=0.              ! Bottom boundary condition
           CF(i,0)=0.              ! Bottom boundary condition
         enddo

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=   
     &          -(z_nbq(i,j,k,knew)-z_nbq(i,j,k,kstp))/dtgrid_nbq
     &          *0.5*( (rho(i,j,k  )/rho0)
     &                  +rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k)  
     &                +(rho(i,j,k+1)/rho0)
     &                  +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
              CF(i,k)=   
     &          -(z_nbq(i,j,k,knew)-z_nbq(i,j,k,kstp))/dtgrid_nbq 
     &          *0.5*( 1.+rho(i,j,k  )/rho0
     &                  +rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k) 
     &                +1.+rho(i,j,k+1)/rho0
     &                  +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1) )

       thetadiv2_nbq(i,j,k)=(FC(i,k)-FC(i,k-1)) 
# ifdef NBQ_MASS
     &      +(Hzr(i,j,k)*(rho(i,j,k)/rho0)-rho_bak(i,j,k))/dtgrid_nbq
# endif
       thetadiv3_nbq(i,j,k)=(CF(i,k)-CF(i,k-1)) 
            enddo
          enddo

          do i=Istr,Iend
              FC(i,N)=   
     &          -(z_nbq(i,j,N,knew)-z_nbq(i,j,N,kstp))/dtgrid_nbq
     &          *0.5*( rho(i,j,N)/rho0
     &                  +rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)
     &                +rho(i,j,N)/rho0)
              CF(i,N)=   
     &          -(z_nbq(i,j,N,knew)-z_nbq(i,j,N,kstp))/dtgrid_nbq 
     &          *0.5*( 1.+rho(i,j,N)/rho0
     &                  +rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)
     &                + 1.+rho(i,j,N)/rho0)

       thetadiv2_nbq(i,j,N)=(FC(i,N)-FC(i,N-1)) 
# ifdef NBQ_MASS
     &      +(Hzr(i,j,N)*(rho(i,j,N)/rho0)-rho_bak(i,j,N))/dtgrid_nbq
# endif            
       thetadiv3_nbq(i,j,N)=(CF(i,N)-CF(i,N-1)) 
          enddo
         enddo

# ifdef NBQ_NODS
         endif
# endif

#endif  /* NBQ_DTDRHO2 */

!
!-------------------------------------------------------------------
! Implicit Vertical Momentum equation: 
!-------------------------------------------------------------------
!
#ifdef NBQ_IMP
        qdmw_nbq_old = qdmw_nbq
!  
        do j=Jstr,Jend
          do k=1,N
            do i=Istr,Iend
               FC(i,k)=  (soundspeed2_nbq(i,j)*rho_nbq(i,j,k)
     &         - thetaimp_nbq*(soundspeed2_nbq(i,j)*dtnbq+visc2_nbq) 
     &            *(thetadiv_nbq(i,j,k)+thetadiv2_nbq(i,j,k))
     &                   )
     
       FC(i,k)=FC(i,k)-thetaimp_nbq*(1.-thetaimp_nbq)*
     &   soundspeed2_nbq(i,j)*dtnbq*
     & (Hzw_half_nbq_inv(i,j,k)* qdmw_nbq(i,j,k)
     &  -Hzw_half_nbq_inv(i,j,k-1)* qdmw_nbq(i,j,k-1))
     
       FC(i,k)=FC(i,k)* Hzr_half_nbq_inv(i,j,k) 
            enddo
          enddo
     
!.........Inner layers:
          do k=1,N-1
            do i=Istr,Iend                                                               
              dum_s =   FC(i,k) - FC(i,k+1)            
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &          + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
# if defined MASKING
              qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
# endif               
            enddo             
          enddo

!.........Surface BC:
          k=N
          do i=Istr,Iend                                                               
            dum_s =   FC(i,k)                              
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &         + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
# if defined MASKING
             qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
# endif              
          enddo   

!.........Bottom BC:     

        enddo

!--------------------------- 
! Gaussian Elimination:
!---------------------------

!.......Comptuts coef.

        do j=Jstr,Jend

!..........Bottom BC:
           k=1
           do i=Istr,Iend
             cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq(i,j)
     &             *dtnbq+visc2_nbq)) 
             cff=(cff1+Hzw_half_nbq_inv(i,j,1)*(Hzr_half_nbq_inv(i,j,1)
     &                 +Hzr_half_nbq_inv(i,j,2)))
             CF(i,1)=(-Hzw_half_nbq_inv(i,j,2)*Hzr_half_nbq_inv(i,j,2))/cff
             DC(i,1)=qdmw_nbq(i,j,1)*cff1/cff   
     &                +qdmw_nbq(i,j,0)/cff*Hzw_half_nbq_inv(i,j,0)
     &                  *Hzr_half_nbq_inv(i,j,1)
           enddo

!..........Inner layers:
           do k=2,N-1
             do i=Istr,Iend
             cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq(i,j)
     &             *dtnbq+visc2_nbq)) 
               cff=(cff1+                                                                    
     &               Hzw_half_nbq_inv(i,j,k)*(Hzr_half_nbq_inv(i,j,k)
     &              +Hzr_half_nbq_inv(i,j,k+1))   
     &              +Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)
     &              *CF(i,k-1))
               CF(i,k)=(-Hzw_half_nbq_inv(i,j,k+1)
     &              *Hzr_half_nbq_inv(i,j,k+1))/cff
               DC(i,k)=(qdmw_nbq(i,j,k)*cff1+Hzw_half_nbq_inv(i,j,k-1)
     &              *Hzr_half_nbq_inv(i,j,k)*DC(i,k-1)) /cff           
             enddo            
           enddo

!..........Surface BC:
           k=N
           do i=Istr,Iend
             cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq(i,j)
     &             *dtnbq+visc2_nbq)) 
             cff=(cff1+Hzw_half_nbq_inv(i,j,N)*Hzr_half_nbq_inv(i,j,N) 
     &                     +Hzw_half_nbq_inv(i,j,N-1)
     &                     *Hzr_half_nbq_inv(i,j,N)*CF(i,N-1))  
             CF(i,N)=0. 
             DC(i,k)=(qdmw_nbq(i,j,N)*cff1+Hzw_half_nbq_inv(i,j,N-1)
     &                  *Hzr_half_nbq_inv(i,j,N)*DC(i,N-1))/cff
           enddo 

!..........Solves tri-diag system:
           do i=Istr,Iend
             qdmw_nbq(i,j,N)=DC(i,k)   
!    &            -rho_nbq(i,j,N)*0.5
!    &                 *9.81*dtnbq        

# ifdef NBQ_NUDGING
             qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)*(1.-nudg_coef_nbq(i,j))
     &        +wz(i,j,N,nrhs)*hzw_half_nbq(i,j,N)
     &        * nudg_coef_nbq(i,j)
# endif

           enddo
           do k=N-1,1,-1
             do i=Istr,Iend
               qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
!    &            -0.25*(rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)
!    &                 +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
!    &                *(Hzr(i,j,k)+Hzr(i,j,k+1))
!    &                 *9.81*dtnbq
# ifdef NBQ_NUDGING
             qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)*(1.-nudg_coef_nbq(i,j))
     &        +wz(i,j,k,nrhs)*hzw_half_nbq(i,j,k)
     &        * nudg_coef_nbq(i,j)
# endif

             enddo            
           enddo                        
        enddo    
!---------------------------
! Vertical momentum open boundary conditions
!---------------------------
# ifdef OBC_NBQ
        call wnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
                
#endif /* NBQ_IMP */


#if !defined NBQ_ZETAEXP
!
!-------------------------------------------------------------------
!       Computes surface mean velocities (Zeta
!-------------------------------------------------------------------
       
        if (IstrU.le.Iend) then
         do j=Jstr,Jend
          do i=Istr,Iend+1     
               umean_nbq(i,j)=qdmu_nbq(i,j,N)  
# ifdef NBQ_MASS                        
     &           / ( rho_nbq(i,j,N)  *Hzr_half_nbq_inv(i,j,N)
     &              +rho(i,j,N)/rho0            
     &              +rho_nbq(i-1,j,N)*Hzr_half_nbq_inv(i-1,j,N)
     &              +rho(i-1,j,N)/rho0+2.)  
     &            / (Hzr(i,j,N)+Hzr(i-1,j,N)) * 4. 
# else
     &            / (Hzr(i,j,N)+Hzr(i-1,j,N)) * 2. 
# endif
# ifdef MASKING
     &            * umask(i,j) 
# endif
          enddo 
         enddo 
        endif

        if (JstrV.le.Jend) then
         do j=Jstr,Jend+1
          do i=Istr,Iend     
               vmean_nbq(i,j)=qdmv_nbq(i,j,N)              
# ifdef NBQ_MASS     
     &            / ( rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)
     &               +rho(i,j,N)/rho0   
     &               +rho_nbq(i,j-1,N)*Hzr_half_nbq_inv(i,j-1,N)
     &               +rho(i,j-1,N)/rho0   +2.)  
     &            / (Hzr(i,j,N)+Hzr(i,j-1,N)) * 4. 
# else
     &            / (Hzr(i,j,N)+Hzr(i,j-1,N)) * 2. 
# endif
# ifdef MASKING
     &            * vmask(i,j) 
# endif
          enddo
         enddo 
        endif

        do j=Jstr,Jend
          do i=Istr,Iend
               wmean_nbq(i,j,knew)=qdmw_nbq(i,j,N)         
# ifdef NBQ_MASS     
     &     / (rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)+1.+rho(i,j,N)/rho0)  
# endif
     &             * Hzw_half_nbq_inv(i,j,N)   
# ifdef MASKING
     &             * rmask(i,j) 
# endif
          enddo
        enddo 

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
        if (IstrU.le.Iend) then
         call exchange_u2d_tile (Istr,Iend,Jstr,Jend,umean_nbq(START_2D_ARRAY))
        endif
        if (JstrV.le.Jend) then
         call exchange_v2d_tile (Istr,Iend,Jstr,Jend,vmean_nbq(START_2D_ARRAY))
        endif
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,wmean_nbq(START_2D_ARRAY,knew))
# endif

#endif /* NBQ_ZETAEXP*/
!
!-------------------------------------------------------------------
!      Mass equation (1)
!-------------------------------------------------------------------
!		
!
!.......Computes fluxes:  
!
        do j=Jstr,Jend
            
#ifdef NBQ_FREESLIP
         do i=Istr,Iend
           FC(i,0)=Hzw_half_nbq_inv(i,j,0) *
     &  (thetaimp_nbq*qdmw_nbq(i,j,0)
     &  +(1.-thetaimp_nbq)*qdmw_nbq_old(i,j,0))   ! Bottom boundary condition
         enddo
#else
         do i=Istr,Iend
           FC(i,0)=0.                                                    ! Bottom boundary condition
         enddo
#endif

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=Hzw_half_nbq_inv(i,j,k) *
     &  (thetaimp_nbq*qdmw_nbq(i,j,k)
     &  +(1.-thetaimp_nbq)*qdmw_nbq_old(i,j,k))
	      thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                           +FC(i,k)-FC(i,k-1)  
            enddo
          enddo
            do i=Istr,Iend
              FC(i,N)=Hzw_half_nbq_inv(i,j,N) *
     &  (thetaimp_nbq*qdmw_nbq(i,j,N)
     &  +(1.-thetaimp_nbq)*qdmw_nbq_old(i,j,N))
	      thetadiv_nbq(i,j,N)=thetadiv_nbq(i,j,N)
     &                           +FC(i,N)-FC(i,N-1)    
            enddo
        enddo

!.......Computes rho_nbq:
!                 
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!       call exchange_r3d_tile (Istr,Iend,Jstr,Jend,rho_nbq(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,thetadiv_nbq(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,thetadiv2_nbq(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,thetadiv3_nbq(START_2D_ARRAY,1))
#endif

         do k=1,N
         do j=JstrV-2,Jend+1
         do i=IstrU-2,Iend+1
           rho_nbq(i,j,k) = rho_nbq(i,j,k)  
     &       - dtfast*(thetadiv_nbq(i,j,k)+thetadiv2_nbq(i,j,k))
#ifdef NBQ_NUDGING
             rho_nbq(i,j,k)=rho_nbq(i,j,k)*(1.-nudg_coef_nbq(i,j))
     &        +0. * nudg_coef_nbq(i,j)
#endif

         enddo
         enddo
         enddo
!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
#if defined ACOUSTIC
       call densityijk_nbq(11)       
#endif

#ifdef NBQ_DTDRHO2
	  do j=JstrR2,JendR2
	    do k=1,N
	    do i=IstrR2,IendR2
             zr_nbq(i,j,k,knew)=z_r(i,j,k)
            enddo
            enddo
   	    do k=0,N
	    do i=IstrR2,IendR2
             z_nbq (i,j,k,knew)=z_w(i,j,k)
            enddo
            enddo
         enddo
#endif

!
!-------------------------------------------------------------------
!      rhobar and rho_avg
!-------------------------------------------------------------------
!		
#ifdef NBQ_MASS
# ifdef NBQ_ZETAREDIAG
c Compute rhobar(n+1) (used in the diagnostic of zeta via vertically integrated continuity equation)
         call densityijk_nbq(20)
# endif
#endif
!
!-------------------------------------------------------------------
!      Grid !
!-------------------------------------------------------------------
!
#if !defined NBQ_ZETAEXP
!***********************************************************************
! s-grid does not evolve with fast mode
! 
!***********************************************************************

# define zab3 UFx
!
!-----------------------------------------------------------------------
! Computes zeta(n+1)
!-----------------------------------------------------------------------
!
      if (FIRST_TIME_STEP.and.iif.eq.1) then
         cff8=1.D0
         cff9=0.D0
         cff10=0.D0
       elseif (FIRST_TIME_STEP.and.iif.eq.2) then
         cff8=1.5D0
         cff9=-0.5D0
        cff10=0.D0
      else
        cff8= 1.5D0+mybeta
        cff9=-2.0D0*mybeta-0.5D0
        cff10= mybeta
      endif

      do j=JstrV-2,Jend+1
        do i=IstrU-2,Iend+1
           zab3(i,j) =   cff8  * zeta(i,j,kstp)
     &                 + cff9  * zeta(i,j,kbak2)
     &                 + cff10 * zeta(i,j,kold2)
        enddo
      enddo

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
         zeta (i,j,knew)=
     &    ( zeta(i,j,kstp)
     &        +dtfast*(
# if defined NBQ_ZETAEXP
     &    +wmean_nbq(i,j,kstp)
# else
     &    +wmean_nbq(i,j,knew)
# endif
     &    -0.5*(umean_nbq(i  ,j)
     &          *(zab3(i  ,j)
     &           -zab3(i-1,j))*pm_u(i,j)
# ifdef MASKING
     &         *rmask(i,j)*rmask(i-1,j)
# endif
     &         +umean_nbq(i+1,j)
     &          *(zab3(i+1,j)
     &           -zab3(i  ,j))*pm_u(i+1,j) 
# ifdef MASKING
     &         *rmask(i,j)*rmask(i+1,j)
# endif
     &          )
# ifdef MASKING
     &         *umask(i,j)*umask(i+1,j)
# endif

     &    -0.5*(vmean_nbq(i  ,j)
     &          *(zab3(i,j  ) 
     &           -zab3(i,j-1))*pm_v(i,j)
# ifdef MASKING
     &         *rmask(i,j)*rmask(i,j-1)
# endif
     &         +vmean_nbq(i,j+1)
     &          *(zab3(i,j+1)
     &           -zab3(i,j  ))*pm_v(i,j+1) 
# ifdef MASKING
     &         *rmask(i,j)*rmask(i,j+1)
# endif
     &         )
# ifdef MASKING
     &         *vmask(i,j)*vmask(i,j+1)
# endif
     &         ) )
# ifdef MASKING
     &         *rmask(i,j)
# endif
        enddo
      enddo

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
# endif

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
           zeta_new(i,j)=zeta(i,j,knew)
        enddo
      enddo

!
!-----------------------------------------------------------------------
! Add nudging terms
!-----------------------------------------------------------------------
!
# ifdef ZNUDGING
#   ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew),zetazon)
        endif
      endif
      if (iic.eq.ntstart) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ssh(START_2D_ARRAY),sshzon)
        endif
      endif
#   endif  /* ZONAL_NUDGING */
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j) + dtfast*Znudgcof(i,j)
#   ifdef ZONAL_NUDGING
     &                                 *(sshzon(j)-zetazon(j))
#   else
     &                                 *(ssh(i,j)-zeta_new(i,j))
#   endif /* ZONAL_NUDGING */
        enddo
      enddo
# endif /* ZNUDGING */
!
!-----------------------------------------------------------------------
! Computes zetarhs to use in momentum equations
!-----------------------------------------------------------------------
!

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
 	   zeta_new(i,j)=zeta_new(i,j) SWITCH rmask(i,j)
# if !defined NBQ_ZETAEXP
c in case of ZETAEXP rhobar has not yet been computed
c Dnew will be computed later
          Dnew(i,j)=(zeta_new(i,j)+h(i,j))
#  if defined NBQ_MASS
     &       *rhobar_nbq(i,j,knew)
#  endif
# endif
        enddo
      enddo

!
!-----------------------------------------------------------------------
! Load new free-surface values into shared array
! Modify new free-surface to ensure that depth is > Dcrit for masked
! cells.
!-----------------------------------------------------------------------
!

# if defined WET_DRY && defined MASKING
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   (Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
        enddo
      enddo 
# endif
!
!-----------------------------------------------------------------------
! Load rhs values into additional AGRIF shared array for nesting
!-----------------------------------------------------------------------
! 
# ifdef AGRIF
      if (FIRST_2D_STEP) then
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,0)=zeta(i,j,kstp)       
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
          du_avg3(i,j,0)  = DUon(i,j)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
          dv_avg3(i,j,0)  = DVom(i,j)
          enddo
        enddo 
      endif

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   Zt_avg3(START_2D_ARRAY,0))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   du_avg3(START_2D_ARRAY,0))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   dv_avg3(START_2D_ARRAY,0))
#  endif

#  ifdef RVTK_DEBUG_ADVANCED
       if (.not.agrif_Root()) then
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(Zt_avg3(:,:,0),'Zt_avg3 (index 0) step2d','r')
       call check_tab2d(DU_avg3(:,:,0),'DU_avg3 (index 0) step2d','u')
       call check_tab2d(DV_avg3(:,:,0),'DV_avg3 (index 0) step2d','v')
C$OMP END MASTER  
       endif
#  endif  
# endif /* AGRIF */   
!
!-----------------------------------------------------------------------
! Compute wet/dry masks
!-----------------------------------------------------------------------
!
# ifdef WET_DRY
      call wetdry_tile (Istr,Iend,Jstr,Jend)
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
# endif   
!
!-----------------------------------------------------------------------
! Set boundary conditions for the free-surface
!-----------------------------------------------------------------------
!
      call zetabc_tile (Istr,Iend,Jstr,Jend)

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
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
! Update Grid:
!-----------------------------------------------------------------------
!
       if ((iic.eq.1.and.iif==1)
     &     NSTEP_GRID) then
        call set_depth_tile(Istr,Iend,Jstr,Jend
     &   ,resetfromrhobar
     &   ) 

!
!-----------------------------------------------------------------------
! s-grid for NH
! s-grid evolves with external mode (NH)
!-----------------------------------------------------------------------
!
       call grid_coef_nh(
     &   Istr,Iend,Jstr,Jend,
     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &   Hzu_half_qdmu, Hzv_half_qdmv                                     
     &   )

        endif

#  undef zab3

# endif /* SOLVE3D */


#else /* ifdef ZETAEXP */
!
!***********************************************************************
! need to update Dnew here (ZETAEXP case)
!***********************************************************************
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
# if defined NBQ_MASS
     &    *rhobar_nbq(i,j,knew)
# endif
        enddo
      enddo

#endif  /* ! ZETAEXP */

!
!***********************************************************************
! Some updates
!***********************************************************************
!
!-------------------------------------------------------------------
! Save rho_bak
!-------------------------------------------------------------------
!
#if defined NBQ_DTDRHO2
      do k=1,N 
       do j=JstrR2,JendR2             
        do i=IstrR2,IendR2
           rho_bak(i,j,k)=(    rho(i,j,k)/rho0)*Hzr(i,j,k)
        enddo
       enddo
      enddo
#endif
            
#ifdef RVTK_DEBUG
       call check_tab3d(rho_nbq,'rho_nbq','r')
#endif    

!-------------------------------------------------------------------
!......Set NBQ/EXT coupling terms
!-------------------------------------------------------------------
!
!          call ruijk_nbq(2, Istr,Iend,Jstr,Jend,WORK)


#ifdef M2FILTER_NONE
        if (LAST_2D_STEP) then
#endif 
          
        do k=0,N 
          do j=Jstr,Jend              
            do i=Istr,Iend
              rw_nbq_ext (i,j,k) = ((qdmw_nbq(i,j,k)-rw_nbq_ext(i,j,k))
     &              /dtnbq-ndtnbq*rw_int_nbq(i,j,k))/WORK(i,j)
            enddo
          enddo
        enddo
#ifdef M2FILTER_NONE        
        endif
#endif 


#ifdef RVTK_DEBUG
          call check_tab2d(rubar_nbq,'rubar_nbq step3d_nbq','uint')
          call check_tab2d(rvbar_nbq,'rvbar_nbq step3d_nbq','vint')
!         call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext step3d_nbq','r')
#endif  
!    
#ifdef NBQ_MASS
# if !defined NBQ_ZETAREDIAG

c densityijk_nbq has already been call is ZETAREDIAG is defined

         call densityijk_nbq(20)
# endif
#endif

        if (LAST_2D_STEP) then
!     
!-----------------------------------------------------------------------
! Get filtered rhs terms (M2FILTER_NONE)
! and multiply by dx*dy to get units of rho*Hz*dx*dy*ru
!
! Finalizes rho-avg1 variables
!-----------------------------------------------------------------------
!     
#ifdef NBQ_MASS
          do j=Jstr,Jend
            do i=Istr,Iend
              rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,knew)
            enddo
          enddo 
          do k=1,N
            do j=Jstr,Jend
              do i=Istr,Iend
                rho_nbq_avg1(i,j,k)=1.d0
     &        + (rho_nbq(i,j,k)/Hzr(i,j,k)+rho(i,j,k)/rho0 )   
              enddo
            enddo 
          enddo
#endif /* NBQ_MASS */
!
!-----------------------------------------------------------------------
! Finalizes Ru_avg2 variables
!-----------------------------------------------------------------------
!
          do k=1,N
            do j=Jstr,Jend
              do i=IstrU,Iend
     !    !     ru_nbq_avg1(i,j,k)=ru_nbq_ext(i,j,k)
                ru_int_nbq(i,j,k) = ru_int_nbq(i,j,k)
     &             -ruext_nbq_2d_old(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
                ru_nbq_avg2(i,j,k)=
     &             ((qdmu_nbq(i,j,k)-ru_nbq_avg2(i,j,k))/dt
     &             -ru_int_nbq(i,j,k)-(ruext_nbq_2d_sum(i,j)/nfast)*
     &             (Hz(i,j,k)+Hz(i-1,j,k)))*on_u(i,j)*om_u(i,j)
              enddo
            enddo 
          enddo          

          do k=1,N
            do j=JstrV,Jend
              do i=Istr,Iend             
          !     rv_nbq_avg1(i,j,k)=rv_nbq_ext(i,j,k) 
                rv_int_nbq(i,j,k) = rv_int_nbq(i,j,k)
     &            -rvext_nbq_2d_old(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
                rv_nbq_avg2(i,j,k)=
     &             ((qdmv_nbq(i,j,k)-rv_nbq_avg2(i,j,k))/dt
     &             -rv_int_nbq(i,j,k)-(rvext_nbq_2d_sum(i,j)/nfast)*
     &             (Hz(i,j,k)+Hz(i,j-1,k)))*on_v(i,j)*om_v(i,j)  
              enddo
            enddo 
          enddo

          do k=1,N   ! Francis: switch from k=0 to k=1 here
            do j=Jstr,Jend
              do i=Istr,Iend
         !      rw_nbq_avg1(i,j,k)=rw_nbq_ext(i,j,k)
                rw_nbq_avg2(i,j,k)=
     &             ((qdmw_nbq(i,j,k)-rw_nbq_avg2(i,j,k))/dt
     &             -rw_int_nbq(i,j,k))*on_r(i,j)*om_r(i,j)
              enddo
            enddo 
          enddo   

        endif

!
!-----------------------------------------------------------------------
!  Gets rid of NBQ NBQ2EXT & NBQ2INT coupling
!-----------------------------------------------------------------------
!
#ifdef NBQ_NOCOUPLING 
          rubar_nbq   =0.   ! 2D
          rvbar_nbq   =0.
          ru_nbq_ext  =0.   ! 3D
          rv_nbq_ext  =0.
          rw_nbq_ext  =0.
          ru_nbq_avg2 =0.
          rv_nbq_avg2 =0.
          rw_nbq_avg2 =0.
# ifdef NBQ_MASS
          rhobar_nbq  =1.
          rho_nbq_ext =1.
          rho_nbq_avg1=1.
# endif
#endif

!         rhobar_nbq  =1.
!         rho_nbq_ext =1.
!         rho_nbq_avg1=1.

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
!*********************************************************************
! 2D depth-averaged velocity & associated coupling
!*********************************************************************
!
! Output: (ubar,vbar), (DU_avg1,DV_avg1)
!

#define Dstp DUon

      cff=0.5*dtfast
#ifdef SOLVE3D
      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)
#else
      cff2=2.*dtfast
#endif

      if (iif==1) then
         DU_avg1=0.
         DV_avg1=0.
         DU_avg2=0.
         DV_avg2=0.
      endif
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew=DU_nbq(i,j) *2.
#ifdef MASKING
     &                                         *umask(i,j)
#endif
#ifdef WET_DRY
          cff1_WD=ABS(ABS(umask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DUnew)*umask_wet(i,j)
          umask_wet(i,j)=0.5*umask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.-cff1_WD)
          DUnew=DUnew*umask_wet(i,j)
# ifdef MRL_WCI
          ust2d(i,j)=ust2d(i,j)*umask_wet(i,j)
# endif
#endif
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
#ifdef SOLVE3D
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*on_u(i,j)*( DUnew
# ifdef MRL_WCI
     &                 +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
# endif
     &                                                   )
          DU_avg2(i,j)=DU_avg2(i,j)+cff2*on_u(i,j)*( DUnew 
#  ifdef MRL_WCI
     &                 +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#  endif
     &                                                   )
#endif
        enddo
      enddo 

  
      do j=JstrV,Jend
        do i=Istr,Iend
          DVnew=DV_nbq(i,j) *2.
#ifdef MASKING
     &                                         *vmask(i,j)
#endif
#ifdef WET_DRY
          cff1_WD=ABS(ABS(vmask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DVnew)*vmask_wet(i,j)
          vmask_wet(i,j)=0.5*vmask_wet(i,j)*cff1_WD
     &                        +cff2_WD*(1.-cff1_WD)
          DVnew=DVnew*vmask_wet(i,j)
# ifdef MRL_WCI
          vst2d(i,j)=vst2d(i,j)*vmask_wet(i,j)
# endif
#endif
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
#ifdef SOLVE3D
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*om_v(i,j)*(DVnew
# ifdef MRL_WCI
     &                 +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
# endif
     &                                                   )
          DV_avg2(i,j)=DV_avg2(i,j)+cff2*on_v(i,j)*( DVnew
# ifdef MRL_WCI
     &                 +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
# endif
     &                                                   )
#endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
!  Set 2D Momemtum nudging
!-----------------------------------------------------------------------
!
#if defined M2NUDGING && defined M2CLIMATOLOGY

# ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew),ubzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew),vbzon)
        endif
      endif
      if (iic.eq.ntstart) then
        if (FIRST_2D_STEP) then 
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubclm(START_2D_ARRAY),ubclmzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbclm(START_2D_ARRAY),vbclmzon)
        endif
      endif
# endif /* ZONAL_NUDGING */

      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef ZONAL_NUDGING        
          DUnew = dtfast*M2nudgcof(i,j)*(ubclmzon(j)-ubzon(j))
# else          
          DUnew = dtfast*M2nudgcof(i,j)*(ubclm(i,j)-ubar(i,j,knew))
# endif
# ifdef MASKING
     &                 * umask(i,j)
# endif 
# ifdef WET_DRY
     &                 * umask_wet(i,j)
# endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
# ifdef SOLVE3D
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# endif
        enddo
      enddo
      
      do j=JstrV,Jend
        do i=Istr,Iend
# if defined ZONAL_NUDGING  
          DVnew = dtfast*M2nudgcof(i,j)*(vbclmzon(j)-vbzon(j)) 
# else
          DVnew = dtfast*M2nudgcof(i,j)*(vbclm(i,j)-vbar(i,j,knew))    
# endif      
# ifdef MASKING
     &                 * vmask(i,j)
# endif 
# ifdef WET_DRY
     &                 * vmask_wet(i,j)
# endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
# ifdef SOLVE3D
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# endif

        enddo
      enddo
#endif /* M2NUDGING */
!
!-----------------------------------------------------------------------
!  Body force for the Internal Tide test case
!-----------------------------------------------------------------------
!
#if defined INTERNAL || defined BODYTIDE
      omega=2.*pi/(12.4*3600.)
      U0=0.02
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew = dtfast*omega*U0*cos(omega*time)
# ifdef MASKING
     &                 * umask(i,j)
# endif
# ifdef WET_DRY
     &                 * umask_wet(i,j)
# endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
# ifdef SOLVE3D
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# endif
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          DVnew = dtfast*0.5*(f(i,j)+f(i,j-1))*
     &                   U0*sin(omega*time)
# ifdef MASKING
     &                 * vmask(i,j)
# endif
# ifdef WET_DRY
     &                 * vmask_wet(i,j)
# endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
# ifdef SOLVE3D
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# endif

        enddo
      enddo
#endif /* INTERNAL */  

!
!----------------------------------------------------------------------- 
! Set boundary conditions and compute integral mass flux accross
! all open boundaries, if any.
!-----------------------------------------------------------------------
!
#if defined NBQ_GRIDEXT || !defined NBQ || !defined NBQ_ZETAW
      call u2dbc_tile (Istr,Iend,Jstr,Jend, UFx) 
      call v2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
#endif

#ifdef NBQ_ZETAREDIAG
c recompute zeta via vertically integrated continuity equation
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   DU_avg1(START_2D_ARRAY,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   DV_avg1(START_2D_ARRAY,knew))
# endif
#endif
          
#ifdef OBC_VOLCONS
      call obc_flux_tile (Istr,Iend,Jstr,Jend)
#endif

!
!-----------------------------------------------------------------------
! Compute fast-time averaged barotropic mass fluxes along physical
! boundaries.
!-----------------------------------------------------------------------
!
#ifdef SOLVE3D
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Istr-1,j)=(h(Istr-1,j)+zeta(Istr-1,j,knew))
#  if defined NBQ_MASS
     &   *rhobar_nbq(Istr-1,j,knew)
#  endif
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Iend+1,j)=(h(Iend+1,j)+zeta(Iend+1,j,knew))
#  if defined NBQ_MASS
     &   *rhobar_nbq(Iend+1,j,knew)
#  endif
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jstr-1)=(h(i,Jstr-1)+zeta(i,Jstr-1,knew))
#  if defined NBQ_MASS
     &   *rhobar_nbq(i,Jstr-1,knew)
#  endif
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jend+1)=(h(i,Jend+1)+zeta(i,Jend+1,knew))
#  if defined NBQ_MASS
     &   *rhobar_nbq(i,Jend+1,knew)
#  endif
        enddo
      endif
# endif
      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)
#  if defined NBQ_GRIDEXT || !defined NBQ || !defined NBQ_ZETAW
     &         +cff1*(Dnew(IstrU-1,j)
     &         +Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
#  else
     &         +cff1*(DU_nbq(IstrU-1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                             +ust2d(IstrU-1,j)
#  endif
     &                                             )*on_u(IstrU-1,j)

          DU_avg2(IstrU-1,j)=DU_avg2(IstrU-1,j)
#  if defined NBQ_GRIDEXT || !defined NBQ
     &         +cff2*(Dnew(IstrU-1,j)
     &         +Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
#  else
     &         +cff2*(DU_nbq(IstrU-1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                             +ust2d(IstrU-1,j)
#  endif
     &                                             )*on_u(IstrU-1,j)
        enddo
        do j=JstrV,Jend
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)
#  if defined NBQ_GRIDEXT || !defined NBQ || !defined NBQ_ZETAW
     &       +cff1*(Dnew(Istr-1,j)
     &       +Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
#  else
     &       +cff1*(DV_nbq(Istr-1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                              +vst2d(Istr-1,j)
#  endif
     &                                              )*om_v(Istr-1,j)
          DV_avg2(Istr-1,j)=DV_avg2(Istr-1,j)
#  if defined NBQ_GRIDEXT || !defined NBQ
     &       +cff2*(Dnew(Istr-1,j)
     &       +Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
#  else
     &       +cff2*(DV_nbq(Istr-1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                              +vst2d(Istr-1,j)
#  endif
     &                                              )*om_v(Istr-1,j)
        enddo
      endif
        
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)
#  if defined NBQ_GRIDEXT || !defined NBQ || !defined NBQ_ZETAW
     &            +cff1*( Dnew(Iend+1,j)
     &            +Dnew(Iend,j) )*(ubar(Iend+1,j,knew)
#  else
     &            +cff1*(DU_nbq(Iend+1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                              +ust2d(Iend+1,j)
#  endif
     &                                              )*on_u(Iend+1,j)
          DU_avg2(Iend+1,j)=DU_avg2(Iend+1,j)
#  if defined NBQ_GRIDEXT || !defined NBQ
     &            +cff2*( Dnew(Iend+1,j)
     &            +Dnew(Iend,j) )*(ubar(Iend+1,j,knew)
#  else
     &            +cff2*(DU_nbq(Iend+1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                              +ust2d(Iend+1,j)
#  endif
     &                                              )*on_u(Iend+1,j)
        enddo
        do j=JstrV,Jend
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)
#  if defined NBQ_GRIDEXT || !defined NBQ  || !defined NBQ_ZETAW
     &        +cff1*( Dnew(Iend+1,j)
     &        +Dnew(Iend+1,j-1) )*(vbar(Iend+1,j,knew)
#  else
     &        +cff1*(DV_nbq(Iend+1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                              +vst2d(Iend+1,j)
#  endif
     &                                              )*om_v(Iend+1,j)
          DV_avg2(Iend+1,j)=DV_avg2(Iend+1,j)
#  if defined NBQ_GRIDEXT || !defined NBQ
     &        +cff2*( Dnew(Iend+1,j)
     &        +Dnew(Iend+1,j-1) )*(vbar(Iend+1,j,knew)
#  else
     &        +cff2*(DV_nbq(Iend+1,j)*2.
#  endif
#  ifdef MRL_WCI
     &                                              +vst2d(Iend+1,j)
#  endif
     &                                              )*om_v(Iend+1,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=IstrU,Iend
          DU_avg1(i,Jstr-1,nnew)=DU_avg1(i,Jstr-1,nnew)
     &        +cff1*( Dnew(i,Jstr-1)
     &        +Dnew(i-1,Jstr-1) )*(ubar(i,Jstr-1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jstr-1)
#  endif
     &                                              )*on_u(i,Jstr-1)
          DU_avg2(i,Jstr-1)=DU_avg2(i,Jstr-1)
     &        +cff2*( Dnew(i,Jstr-1)
     &        +Dnew(i-1,Jstr-1) )*(ubar(i,Jstr-1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jstr-1)
#  endif
     &                                              )*on_u(i,Jstr-1)
        enddo
        do i=IstrR,IendR
          DV_avg1(i,JstrV-1,nnew)=DV_avg1(i,JstrV-1,nnew)
     &         +cff1*(Dnew(i,JstrV-1)
     &         +Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
#  ifdef MRL_WCI
     &                                              +vst2d(i,JstrV-1)
#  endif
     &                                              )*om_v(i,JstrV-1)
          DV_avg2(i,JstrV-1)=DV_avg2(i,JstrV-1)
     &         +cff2*(Dnew(i,JstrV-1)
     &         +Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
#  ifdef MRL_WCI
     &                                              +vst2d(i,JstrV-1)
#  endif
     &                                              )*om_v(i,JstrV-1)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=IstrU,Iend
          DU_avg1(i,Jend+1,nnew)=DU_avg1(i,Jend+1,nnew)
     &        +cff1*( Dnew(i,Jend+1)
     &        +Dnew(i-1,Jend+1) )*(ubar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                               +ust2d(i,Jend+1)
#  endif
     &                                               )*on_u(i,Jend+1)
          DU_avg2(i,Jend+1)=DU_avg2(i,Jend+1)
     &        +cff2*( Dnew(i,Jend+1)
     &        +Dnew(i-1,Jend+1) )*(ubar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                               +ust2d(i,Jend+1)
#  endif
     &                                               )*on_u(i,Jend+1)
        enddo
        do i=IstrR,IendR
          DV_avg1(i,Jend+1,nnew)=DV_avg1(i,Jend+1,nnew)
     &            +cff1*( Dnew(i,Jend+1)
     &            +Dnew(i,Jend) )*(vbar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                               +vst2d(i,Jend+1)
#  endif
     &                                               )*om_v(i,Jend+1)
          DV_avg2(i,Jend+1)=DV_avg2(i,Jend+1)
     &            +cff2*( Dnew(i,Jend+1)
     &            +Dnew(i,Jend) )*(vbar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                               +vst2d(i,Jend+1)
#  endif
     &                                               )*om_v(i,Jend+1)
        enddo
      endif
# endif
#endif /* SOLVE3D */

#ifdef NBQ_ZETAREDIAG
!
!*********************************************************************
!
! Adjust ZETA using vertically integrated continuity equation
!
!*********************************************************************
       
       do j=Jstr-1,Jend+1
       do i=Istr-1,Iend+1
          zeta_new(i,j)=zeta_new(i,j) SWITCH rmask(i,j)
          Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
# ifdef NBQ_MASS
     &    *rhobar_nbq(i,j,knew)
# endif
        enddo
      enddo
# ifndef NBQ_GRIDEXT
       if ((iic.eq.1.and.iif==1)
     &     .or.iif==nfast
     &     ) then
# endif
!
!-----------------------------------------------------------------------
! Update zeta
!-----------------------------------------------------------------------
!

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=((h(i,j)+zeta(i,j,kstp))
# ifdef NBQ_MASS
     &             *rhobar_nbq(i,j,kstp) 
# endif
     &   + (dtfast*pm(i,j)*pn(i,j)
     &              *0.5D0*(
     &    (Dnew(i  ,j)+Dnew(i-1,j))*(ubar(i  ,j,knew))*on_u(i  ,j)
     &   -(Dnew(i+1,j)+Dnew(i  ,j))*(ubar(i+1,j,knew))*on_u(i+1,j)
     &  + (Dnew(i,j  )+Dnew(i,j-1))*(vbar(i,j  ,knew))*om_v(i,j  )
     &   -(Dnew(i,j+1)+Dnew(i,j  ))*(vbar(i,j+1,knew))*om_v(i,j+1))))
# ifdef NBQ_MASS
     &   / rhobar_nbq(i,j,knew)-h(i,j)
# endif

        enddo
      enddo
             
      do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,knew)
          enddo
      enddo
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   Zt_avg1(START_2D_ARRAY))
# endif
!
!-----------------------------------------------------------------------
! Update some grid parameters
!-----------------------------------------------------------------------
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
!           dum_s=0.
            do k=1,N
              Hz(i,j,k)=Hz_bak2(i,j,k)
     &        - dtfast*(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k))

!             Hz(i,j,k)=
!     &        cff4*(Hz(i,j,k)
!     &        - dtfast*(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k)))
!     &        +cff5*Hz_t(i,j,k,kstp)
!     &        +cff6*Hz_t(i,j,k,kbak2)
!     &        +cff7*Hz_t(i,j,k,kold2)
!             Hz_t(i,j,k,knew)=Hz(i,j,k)
!             dum_s=dum_s+Hz(i,j,k)

              Hzr(i,j,k) =(Hz(i,j,k)-rho_nbq(i,j,k))
     &                   /(1.+rho(i,j,k)/rho0)

!            enddo
!            do k=1,N

              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
              z_r(i,j,k)=(z_w(i,j,k)+z_w(i,j,k-1))/2.D0
            enddo

        enddo
      enddo

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI      
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_w(START_2D_ARRAY,0))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_r(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hzr(START_2D_ARRAY,1))
# endif

!
!***********************************************************************
! Adjuste some grid parameters
!***********************************************************************
!
!-----------------------------------------------------------------------
! Computes grid at m
!-----------------------------------------------------------------------
!
        call grid_coef_nh(
     &   Istr,Iend,Jstr,Jend,
     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &   Hzu_half_qdmu, Hzv_half_qdmv                                     
     &   )

# ifndef NBQ_GRIDEXT
      endif
# endif

#endif /* NBQ_ZETAREDIAG */
!
!-----------------------------------------------------------------------
! Apply point sources for river runoff simulations
!-----------------------------------------------------------------------
!
#ifdef PSOURCE
      do is=1,Nsrc 
# ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
# else
        i=Isrc(is)
        j=Jsrc(is)
# endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            ubar(i,j,knew)=2.*Qbar(is)/( on_u(i,j)
     &                   *(Dnew(i-1,j)+Dnew(i,j)) )
# ifdef SOLVE3D
            DU_avg1(i,j,nnew)=Qbar(is)
# endif
          else
            vbar(i,j,knew)=2.*Qbar(is)/( om_v(i,j)
     &                   *(Dnew(i,j-1)+Dnew(i,j)) )
# ifdef SOLVE3D
            DV_avg1(i,j,nnew)=Qbar(is)
# endif
          endif
        endif
      enddo
#endif
!
!-----------------------------------------------------------------------
!  Diagnostics
!-----------------------------------------------------------------------
!
#ifndef SOLVE3D
      call diag_tile (Istr,Iend,Jstr,Jend, UFx,UFe)
#endif

!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg1(START_2D_ARRAY,nnew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg1(START_2D_ARRAY,nnew))
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
!  Apply conservation requirements for nesting
!-----------------------------------------------------------------------
!
#ifdef AGRIF
 
      if (.NOT.Agrif_Root()) THEN
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,iif)=zeta(i,j,knew)
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
            DU_avg3(i,j,iif) = 0.5*(h(i,j)+zeta(i,j,knew)+
     &          h(i-1,j)+zeta(i-1,j,knew)) *on_u(i,j)*ubar(i,j,knew)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
            DV_avg3(i,j,iif) = 0.5*(h(i,j)+zeta(i,j,knew)+
     &          h(i,j-1)+zeta(i,j-1,knew)) *om_v(i,j)*vbar(i,j,knew)
          enddo
        enddo
 
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Zt_avg3(START_2D_ARRAY,iif))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg3(START_2D_ARRAY,iif))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg3(START_2D_ARRAY,iif))
# endif

# ifdef RVTK_DEBUG_ADVANCED
       call check_tab2d(Zt_avg3(:,:,iif),'Zt_avg3 step2d','r')
       call check_tab2d(DU_avg3(:,:,iif),'DU_avg3 step2d','u')
       call check_tab2d(DV_avg3(:,:,iif),'DV_avg3 step2d','v')
# endif   
      endif

# ifdef AGRIF_CONSERV_VOL
       if (iif == nfast) then
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
# endif
 
#endif /* AGRIF */

#ifdef NBQ_TRACERS
!
!*********************************************************************
! ZETAW: correct Hz
!*********************************************************************
!
       if (iif==nfast) then

        Hz_tra = Hz

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg2(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg2(START_2D_ARRAY))
# endif

!      dum_max=0.
       do j=Jstr,Jend
       do i=Istr,Iend
           dum_s=0.
           do k=1,N
             dum_s=dum_s+
     &           hz(i,j,k)-hz_bak(i,j,k)
            enddo
            do k=1,N

c LAURENT: remove this computation of Hz_corr and its definition in nbq.h
c Kept for testing... removed soon.

             Hz_corr(i,j,k)=(dum_s
     &           + (DU_avg2(i+1,j)-DU_avg2(i,j)
     &             +DV_avg2(i,j+1)-DV_avg2(i,j)
     &               )*pm(i,j)*pn(i,j)
     &        *dt)     !/float(N)
     &                /(z_w(i,j,N)-z_w(i,j,0))
     &                *(z_w(i,j,k)-z_w(i,j,k-1))
!            dum_max=max(dum_max,abs(Hz_corr(i,j,k)))
# ifdef MASKING
     &        *rmask(i,j)
# endif

             Hz(i,j,k)=Hz(i,j,k)
     &         -(dum_s
     &           + (DU_avg2(i+1,j)-DU_avg2(i,j)
     &             +DV_avg2(i,j+1)-DV_avg2(i,j)
     &               )*pm(i,j)*pn(i,j)
     &        *dt) ! / float(N)
     &                /(z_w(i,j,N)-z_w(i,j,0))
     &                *(z_w(i,j,k)-z_w(i,j,k-1))
# ifdef MASKING
     &        *rmask(i,j)
# endif
!              Hzr(i,j,k) =(Hz(i,j,k)-rho_nbq(i,j,k))      ! Zw and Zr to be adjusted here ?
!     &                   /(1.+rho(i,j,k)/rho0)
!              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
!              z_r(i,j,k)=(z_w(i,j,k)+z_w(i,j,k-1))/2.D0
             enddo
       enddo
       enddo
!      write(6,*) iic,"Max=",dum_max

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!     call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!    &                        Hzr(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
c LAURENT: these calls for z_w and z_r can probably be removed
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_w(START_2D_ARRAY,0))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_r(START_2D_ARRAY,1))
# endif 
       endif
#endif  /* NBQ_TRACERS */

!      
!    THE END
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

#ifdef MPI
        write(stdout,'(A,I4)') '      mynode =',mynode
#endif

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

#undef zwrk 
#undef rzeta  
#undef rzeta2  
#undef rzetaSA 

      return
      end
#else
      subroutine step3d_fast_zetaw_empty
      end
#endif /* NBQ */
