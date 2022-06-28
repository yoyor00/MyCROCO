! $Id: pre_step3d.F 1568 2014-06-30 15:57:49Z gcambon $
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
#if defined SOLVE3D && defined NBQ2
# ifndef VADV_ADAPT_PRED
#  undef VADV_ADAPT_IMP
# endif

      subroutine pre_step3d_2 (tile)
!
      implicit none
      integer tile,  trd,omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "ocean3d.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call pre_step3d_tile (Istr,Iend,Jstr,Jend,
     &                A3d(1,1,trd), A3d(1,2,trd), A3d(1,3,trd),
     &                A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
# ifdef VADV_ADAPT_IMP     
     &                                            A2d(1,4,trd), 
# endif           
     &                A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd)
     &                                          , A3d(1,4,trd)         
     &                                                        )
      return
      end

      subroutine pre_step3d_tile (Istr,Iend,Jstr,Jend, ru,rv,rw,
     &                                                 FC,CF,DC, 
# ifdef VADV_ADAPT_IMP     
     &                                                       WC, !CAUTION!!!!
# endif     
     &                                               FX,FE,WORK
     &                                                 ,Hz_half
     &                                                         )  
!
!--------------------------------------------------------------------
! Preliminary step: initialize computations of the new time step
! 3D primitive variables.
!
! Since r.h.s. arrays ru,rv,rt(:,:,???[,:]), which at this moment
! contain r.h.s at time step n-2 will be overwritten by the
! subsequent routines within rhs3d driver, both [n-1 and n-2] old-
! time-step r.h.s. term in Adams-Bashforth stepping scheme are
! added at this time to the time step [n] fields and the result is
! stored as u,v,t(:,:,???[,:]).
!
! The actual time step will be completed in step3d, after the
! time step [n] r.h.s. terms and new-time step Hz will be available
! after the completion rhs3d computations and the 2D (barotropic
! mode) computations.
!--------------------------------------------------------------------
!
      implicit none
# include "param.h"
      integer Istr,Iend,Jstr,Jend, itrc, i,j,k, indx,k1
     &       ,imin,imax,jmin,jmax,nadv,iAkt
# ifdef PSOURCE
     &       ,is
# endif
      real   ru(PRIVATE_2D_SCRATCH_ARRAY,N),    cff,cff0,
     &       rv(PRIVATE_2D_SCRATCH_ARRAY,N),    cff1,
     &       rw(PRIVATE_2D_SCRATCH_ARRAY,0:N),  cff2,
     &       FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),  cff3,
     &       CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),  cff4,cff5,
     &       DC(PRIVATE_1D_SCRATCH_ARRAY,0:N),  gamma, 
     &       FX(PRIVATE_2D_SCRATCH_ARRAY),      epsil,
     &       FE(PRIVATE_2D_SCRATCH_ARRAY),        cdt,
     &       WORK(PRIVATE_2D_SCRATCH_ARRAY)
# if !defined M3FAST_COUPLING1 || !defined M3FAST_COUPLING1W
      real   cff6,cff7
# endif
# ifdef VADV_ADAPT_IMP      
      real   WC(PRIVATE_1D_SCRATCH_ARRAY,0:N) 
# endif 
      real Hz_half(PRIVATE_2D_SCRATCH_ARRAY,N)
      real cff10,cff11
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
      real   ru_int_nbq_2d(PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_int_nbq_2d(PRIVATE_2D_SCRATCH_ARRAY)
# endif
# ifdef M3FAST_DISSCOUPLING
      real alpha_uv,alpha_w
# endif
     
!# ifdef WET_DRY
!      real cff3,cff4,cff5
!#  endif
!# ifdef M3FAST
      real qdm_nstp,qdm_indx
!# endif
      parameter (gamma=1./6., epsil=1.E-16)
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "ocean2d.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
#ifdef OPENACC
# include "openacc.h"
#endif
# ifdef PSOURCE
#  include "sources.h"
# endif
# if defined M3FAST /* M3FAST_C3D_UVSF || defined M3FAST_C3D_WSF */
#  include "nbq_2.h"
# endif
#include "openacc_work.h"
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif

# ifdef M3FAST_C3D_UVFS
#  define ru_nbq_avg1 ru_nbq
#  define rv_nbq_avg1 rv_nbq
# endif
# ifdef M3FAST_C3D_WFS
#   define rw_nbq_avg1 rw_nbq
# endif

!
!--------------------------------------------------------------------
! Definition of flux operators: 1st, 2nd, 3rd, 4th, 5th or 6th order,
! used in UP5 and C6 advection schemes (and order degradation near 
! land masks). cdiff is part of laplacian diffusion in flux1 (used 
! near mask): 
!    0 --> flux1=flux2 (second order C2 advection scheme)
!    1 --> flux1 gives 1st order monotonic UP1 advection scheme
!--------------------------------------------------------------------
!
# if defined TS_HADV_UP5 ||  defined TS_HADV_C6 \
      || defined TS_HADV_WENO5 ||  defined BIO_HADV_WENO5 \
      || defined TS_VADV_WENO5

      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
      REAL    :: ua, vel, cdiff, cdif
      REAL    :: flux1, flux2, flux3, flux4, flux5, flux6
      REAL    :: flx2, flx3, flx4, flx5
      REAL    :: mask0, mask1, mask2, mask3

      flux2(q_im1, q_i, ua, cdiff) = 0.5*( q_i + q_im1 )

      flux1(q_im1, q_i, ua, cdiff) = flux2(q_im1, q_i, ua, cdiff) -
     &      0.5*cdiff*sign(1.,ua)*(q_i-q_im1)

      flux4(q_im2, q_im1, q_i, q_ip1, ua) = 
     &      ( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =
     &      flux4(q_im2, q_im1, q_i, q_ip1, ua) +
     &      sign(1.,ua)*((q_ip1 -
     &      q_im2)-3.*(q_i-q_im1))/12.0

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      ( 37.*(q_i+q_im1) - 8.*(q_ip1+q_im2)
     &      +(q_ip2+q_im3) )/60.0

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
     &      -sign(1.,ua)*( 
     &      (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )/60.0
     
# endif
#  if defined TS_HADV_WENO5 || defined BIO_HADV_WENO5 || defined TS_VADV_WENO5
        REAL    :: flux3_weno, flux5_weno
# endif

!
!--------------------------------------------------------------------
!
# include "compute_auxiliary_bounds.h"

!
!
!--------------------------------------------------------------------
! Dissipation in coupling scheme: 
!   0: no dissipation,
!   1: dissipative coupling.
!--------------------------------------------------------------------
!
# ifdef M3FAST_DISSCOUPLING
#  ifdef KH_INST
      alpha_uv=0.75
      alpha_w =0.
#  elif defined TANK
      alpha_uv=0.
      alpha_w =0.
#  elif defined MILES
      alpha_uv=0.
      alpha_w =0. 
#  else      
      alpha_uv=0.
      alpha_w =0.
#  endif
# endif
      indx=3-nstp
      nadv=  nstp   !<-- do not remove: used in compute_tracer_fluxes.h
      
      
      if (FIRST_TIME_STEP) then
        cff=0.5*dt
        cff1=1.
        cff2=0.
      else
        cff=(1.-gamma)*dt
        cff1=0.5+gamma
        cff2=0.5-gamma
      endif
      
#if defined RVTK_DEBUG_PERFRST
       call check_tab3d(Hz(:,:,:),'Hz in prestep3d','r')
       call check_tab3d(Hz_bak(:,:,:),'Hz_bak in prestep3d','r')
       call check_tab3d(Huon(:,:,:),'Huon in prestep3d','r')
#endif
!$acc kernels if(compute_on_device) default(present)
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            Hz_half(i,j,k)=cff1*Hz(i,j,k)+cff2*Hz_bak(i,j,k)
     &        -cff*pm(i,j)*pn(i,j)*( Huon(i+1,j,k)-Huon(i,j,k)
     &                              +Hvom(i,j+1,k)-Hvom(i,j,k)
     &                                  +We(i,j,k)-We(i,j,k-1)
# ifdef VADV_ADAPT_IMP     
     &                                  +Wi(i,j,k)-Wi(i,j,k-1)     
# endif     
     &                                                       )
          enddo
        enddo
      enddo
!$acc end kernels
!
!  Set extended range
!
#  ifdef EW_PERIODIC
      imin=Istr-2
      imax=Iend+2
#  else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=Jstr-2
      jmax=Jend+2
#  else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
#  endif
!
!--------------------------------------------------------------------
!  Start computation of the auxiliary tracer field
!  -----------------------------------------------
!  Once it will be completed, t(:,:,???,:) is effectively halfway
!  between time steps n and n+1. A high spatial order, centered,
!  non-conservative [but constancy preserving!] scheme is used for
!  this auxiliary step. This is done by introducing an artificial
!  continuity equation [Easter, 1993].
!
!  Since this field will be used exclussively to compute the high-
!  order fluxes during subsequent step3d_t operation, the final
!  values of t(i,j,k,??,itrc) alfer step3d_t will be computed in
!  a flux-conservative manner. The overall time step will be both 
!  conservative and constancy preserving.
!
!  This preliminary step shall be done before field t(:,:,:,???,:)
!  loses its meaningful values during the pre-step operation. 
!
!======================================================================
!
!  Compute horizontal advection
!
!======================================================================
!
      do itrc=1,NT
        do k=1,N

# ifdef BIO_HADV_WENO5
          if (itrc .le. 2) then
# endif
!
!----------------------------------------------------------
! Advection for active tracers: Compute fluxes FX and FE
!----------------------------------------------------------
!
#  define PREDICTOR
#  include "compute_horiz_tracer_fluxes.h"
#  undef PREDICTOR

!
!----------------------------------------------------------
! C6 advection for passive tracers (biology ...)
!----------------------------------------------------------
!
# ifdef BIO_HADV_WENO5
          else  ! if (itrc .gt. 2)
#  define FLUX5 flux6
#  define FLUX3 flux4
#  define FLUX2 flux2
#  undef  UP5_MASKING
!
            cdif=0.
#  include "t3dadv_order5.h"
!
#  undef FLUX5
#  undef FLUX3
#  undef FLUX2
#  undef UP5_MASKING
          endif
# endif
!
!----------------------------------------------------------
! Apply point sources for river runoff simulations
!----------------------------------------------------------
!
# ifdef PSOURCE
          do is=1,Nsrc
#  ifdef MPI
           i=Isrc_mpi(is,mynode)
           j=Jsrc_mpi(is,mynode)
#  else
           i=Isrc(is)
           j=Jsrc(is)
#  endif
            if (Istr.le.i .and. i.le.Iend+1
     &                   .and. Jstr.le.j .and. j.le.Jend+1) then
              if (Dsrc(is).eq.0) then
                if (Lsrc(is,itrc)) then
                  FX(i,j)=Huon(i,j,k)*Tsrc(is,k,itrc)
#  ifdef MASKING
                else
                  if (rmask(i,j).eq.0 .and. rmask(i-1,j).eq.1) then
                    FX(i,j)=Huon(i,j,k)*t(i-1,j,k,nstp,itrc)
                  elseif (rmask(i,j).eq.1. .and. rmask(i-1,j).eq.0) then
                    FX(i,j)=Huon(i,j,k)*t(i  ,j,k,nstp,itrc)
                  endif
#  endif
                endif
              else
                if (Lsrc(is,itrc)) then
                  FE(i,j)=Hvom(i,j,k)*Tsrc(is,k,itrc)
#  ifdef MASKING
                else
                  if (rmask(i,j).eq.0 .and. rmask(i,j-1).eq.1) then
                    FE(i,j)=Hvom(i,j,k)*t(i,j-1,k,nstp,itrc)
                  elseif (rmask(i,j).eq.1 .and. rmask(i,j-1).eq.0) then
                    FE(i,j)=Hvom(i,j,k)*t(i,j  ,k,nstp,itrc)
                  endif
#  endif
                endif
              endif
            endif
          enddo
# endif /* PSOURCE */
!
!----------------------------------------------------------
! Finalize horizontal advection: compute flux divergences
!----------------------------------------------------------
!
          if (FIRST_TIME_STEP) then
            cff=0.5*dt
            do j=Jstr,Jend
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=Hz(i,j,k)*t(i,j,k,nstp,itrc)
     &                      -cff*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                            +FE(i,j+1)-FE(i,j))
              enddo
            enddo
          else
            cff=(1.-gamma)*dt
            cff1=0.5+gamma
            cff2=0.5-gamma
            do j=Jstr,Jend
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=cff1*Hz(i,j,k)*t(i,j,k,nstp,itrc)
     &                            +cff2*Hz_bak(i,j,k)*t(i,j,k,indx,itrc)
     &                         -cff*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                               +FE(i,j+1)-FE(i,j))
              enddo
            enddo
          endif

        enddo   ! <-- k
      enddo       ! <-- itrc

!
! Continue computation of the auxiliary tracer field: auxiliary
! continuity equation (the same for all tracers): DC=1/Hz_half_new,
! where Hz_half_new is Hz at time step n+1/2 as it would be computed
! from three-dimensional divergence of volume fluxes Huon,Hvom and W.
!
!======================================================================
!
! Compute vertical advection
!
!======================================================================
!
! Finalize computation of the auxiliary tracer field: compute its
! change due to vertical advection. Computation of vertical advective
! fluxes requires interpolation of tracer values to the verical grid-
! box interfaces (W-points). This can be done by either using
! parabolic spline interpolation or, more simple local cubic
! polynomial [linear interpolation is considered obsolete].
!      
      if (FIRST_TIME_STEP) then
            cdt=0.5*dt
      else
            cdt=(1.-gamma)*dt
      endif 

# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
      do j=Jstr,Jend
        do i=IstrU,Iend
           ru_int_nbq_2d(i,j)=0.
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
           rv_int_nbq_2d(i,j)=0.
        enddo
      enddo
# endif

# if defined TS_HADV_TEST && defined SOLID_BODY_PER
      do k=1,N
        do j=JstR,JendR
          do i=IstrR,IendR
             u(i,j,k,nstp) = u(i,j,k,nnew)
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstR,IendR
             v(i,j,k,nstp) = v(i,j,k,nnew)
          enddo
        enddo
      enddo

!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!#   ifdef THREE_GHOST_POINTS_UV
!      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
!     &                             u(START_2D_ARRAY,1,nstp))
!      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
!     &                             v(START_2D_ARRAY,1,nstp))
!#   else
!      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
!     &                        u(START_2D_ARRAY,1,nstp))
!      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
!     &                        v(START_2D_ARRAY,1,nstp))
!#   endif
!#  endif
# endif

!======================================================================
! Predictor final computations
!======================================================================

        cff4=0.5+gamma
        cff5=0.5-gamma
# ifdef M3FAST
        if (FIRST_TIME_STEP) then
        cff3=0.                        ! This version is designed
        cff2=0.                        ! for coupling during 3D
        cff1=1.                        ! predictor sub-step: here
#  if !defined M3FAST_COUPLING1 || !defined M3FAST_COUPLING1W
        cff6=1.
        cff7=0.
#  endif
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
        cff3= 0.                       ! computed as instantaneous
        cff2=-0.5                      ! value at 3D time step
        cff1= 1.5                      ! "nstp" first, and then
#  if !defined M3FAST_COUPLING1 || !defined M3FAST_COUPLING1W
        cff6=2.-0.6
        cff7=-(1.-0.6)
#  endif
        else                             ! extrapolated half-step
        cff3=0.281105                 ! forward using  AB3-like
        cff2=-0.5-2.*cff3              ! weights optimized for
        cff1=1.5+cff3                 ! maximum stability (with
#  if !defined M3FAST_COUPLING1 || !defined M3FAST_COUPLING1W
        cff6=2.-0.6
        cff7=-(1.-0.6)
#  endif
        endif      
                      ! special care for startup)
# endif
 
!# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!#  ifdef THREE_GHOST_POINTS_TS
!      call exchange_w3d_3pts_tile (Istr,Iend,Jstr,Jend,
!     &            rw(START_2D_ARRAY,0))   
!#  else
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
!     &            rw(START_2D_ARRAY,0))   
!#  endif
!# endif
!!
!-----------------------------------------------------------------------
! Main j-loop 
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend

        do i=Istr,Iend
# ifndef VADV_ADAPT_IMP 
           do k=1,N
             DC(i,k)=1./Hz_half(i,j,k)
           enddo   
# endif
           DC(i,0)=cdt*pn(i,j)*pm(i,j)

        enddo
        
        do itrc=1,NT
!
!----------------------------------------------------------------------
! Compute vertical fluxes FC
!----------------------------------------------------------------------
!
# define PREDICTOR
# include "compute_vert_tracer_fluxes.h"
# undef PREDICTOR
!
# ifdef TS_HADV_TEST
      FC(istr:iend,0:N) = 0.d0
# endif
!
!----------------------------------------------------------------------    
!++ Apply vertical advection for tracers (DC corresponds to 1/Hz_half)
!----------------------------------------------------------------------
!
          do k=1,N
            do i=Istr,Iend
# ifdef VADV_ADAPT_IMP
              t(i,j,k,nnew,itrc)= t(i,j,k,nnew,itrc)  !<-- in this case division by Hz
     &               -DC(i,0)*(FC(i,k)-FC(i,k-1))     !<-- in the tridiagonal matrix problem
# else
              t(i,j,k,nnew,itrc)=DC(i,k)*( t(i,j,k,nnew,itrc)
     &               -DC(i,0)*(FC(i,k)-FC(i,k-1)))
#  ifdef CONST_TRACERS 
              t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)  
#  endif
# endif
            enddo
          enddo            !--> discard FC  
!
!----------------------------------------------------------------------
!  Compute implicit vertical advection/diffusion
!----------------------------------------------------------------------
!
# ifdef VADV_ADAPT_IMP
#  ifdef SALINITY
        iAKt=min(itrc,isalt)
#  else
        iAKt=min(itrc,itemp)
#  endif  
#  define TRIDIAG_TRA
#  undef  TRIDIAG_U
#  undef  TRIDIAG_V 
#   include "tridiag_pred.h"
#  undef TRIDIAG_TRA
# endif          
!
        enddo   !<-- itrc  !--> discard DC
!
!======================================================================
!
! Momentum equations: time stepping to time n+1/2
!
!======================================================================
!   

        do i=IstrU,Iend
            DC(i,0)=pm_u(i,j)*pn_u(i,j)
            DC(i,1)=DC(i,0)/( Zt_avg1(i  ,j)+h(i  ,j)
     &                       +Zt_avg1(i-1,j)+h(i-1,j))
        enddo
        if (FIRST_TIME_STEP) then
          do k=1,N
            do i=IstrU,Iend
# ifdef VADV_ADAPT_IMP              
              cff = 1.
# else
              cff = 2./(Hz_half(i,j,k)+Hz_half(i-1,j,k))
# endif           
              cff0=DC(i,0)*ru(i,j,k)

              u(i,j,k,nnew)=( u(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                       +cdt*cff0 )*cff
     
# ifdef M3FAST_C3D_UVSF
              ru_int_nbq(i,j,k) =cff0
#  ifndef M3FAST_COUPLING1
              ru_intt_nbq(i,j,k)=cff0
#  else
              ru_intt_nbq(i,j,k,nstp)=cff0
#  endif
# endif
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
              ru_int_nbq_2d(i,j)= ru_int_nbq_2d(i,j) ! CAUTION: M3FAST_BOTH
     &                           +ru_int_nbq(i,j,k) 
# endif
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5*( Hz(i  ,j,k)
     &                                         +Hz(i-1,j,k))
            enddo
          enddo
          
        else 
        
          do k=1,N
            do i=IstrU,Iend

# ifdef VADV_ADAPT_IMP              
              cff = 1.
# else
              cff = 2./(Hz_half(i,j,k)+Hz_half(i-1,j,k))
# endif 
              qdm_nstp = u(i,j,k,nstp)*0.5*( Hz(i  ,j,k)
     &                                      +Hz(i-1,j,k))
              qdm_indx = u(i,j,k,indx)*0.5*( Hz_bak(i  ,j,k)
     &                                      +Hz_bak(i-1,j,k))
              u(i,j,k,nnew)=( cff4*qdm_nstp
     &                       +cff5*qdm_indx
     &                       +cdt*DC(i,0)*(ru(i,j,k)
# ifdef M3FAST_C3D_UVFS
     &                                     +ru_nbq_avg1(i,j,k)
# endif     
     &                                      ))*cff 
 
# ifdef M3FAST_C3D_UVSF
#  ifndef M3FAST_COUPLING1
              cff0=((qdm_nstp-qdm_indx) /dt
     &               -DC(i,0)*ru_nbq_avg2(i,j,k))
              ru_int_nbq(i,j,k)=(cff6*cff0
     &                          +cff7*ru_intt_nbq(i,j,k))
#   ifdef MASKING
     &                          *umask(i,j)
#   endif        
              ru_intt_nbq(i,j,k)=cff0
#  else
           cff0=DC(i,0)*ru(i,j,k)! *(1.-gamma)   
#   ifndef M3FAST_ZETAW
!    &            -rufrc(i,j) *DC(i,1)  
!    &            *(Hz(i,j,k)+Hz(i-1,j,k))                    ! Check Hz_half
#   endif
           ru_int_nbq(i,j,k)=(cff1*cff0
     &                         + cff2*ru_intt_nbq(i,j,k,3-nstp)
     &                         + cff3*ru_intt_nbq(i,j,k,nstp))  
#   ifdef MASKING
     &                         *umask(i,j)
#   endif        
           ru_intt_nbq(i,j,k,nstp)=cff0
#  endif /* !M3FAST_COUPLING1 */
              
# ifdef M3FAST_DISSCOUPLING
              cff10=((qdm_nstp-qdm_indx) /dt
     &              -DC(i,0)*ru_nbq_avg2(i,j,k))
              ru_int_nbq(i,j,k) = (1.-alpha_uv) * ru_int_nbq(i,j,k) !ICI
     &                            +alpha_uv * cff10
# endif     
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
              ru_int_nbq_2d(i,j)= ru_int_nbq_2d(i,j) ! CAUTION: M3FAST_BOTH
     &                           +ru_int_nbq(i,j,k) 
# endif
# endif /* M3FAST_C3D_UVSF */
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            enddo
          enddo
          endif ! FIRST_TIME_STEP
 
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
          do i=IstrU,Iend
          do k=1,N
          ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)
     &               -ru_int_nbq_2d(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))
     &                               /( Zt_avg1(i  ,j)+h(i  ,j)
     &                                 +Zt_avg1(i-1,j)+h(i-1,j))
          enddo
          enddo
# endif     
!
c # if defined M3FAST
c #  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
c       call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
c      &                        ru_int_nbq(START_2D_ARRAY,1))
c #  endif
c # endif
c # if defined RVTK_DEBUG 
c C$OMP BARRIER
c C$OMP MASTER
c        call check_tab3d(ru_int_nbq,'ru_int_nbq (pre_step3d)','u')
c C$OMP END MASTER     
c # endif  
!
!======================================================================
!
!  Compute implicit vertical advection/viscosity
!
# ifdef VADV_ADAPT_IMP
#  undef  TRIDIAG_TRA
#  define TRIDIAG_U
#  undef  TRIDIAG_V
#   include "tridiag_pred.h"
# endif
!
!======================================================================
!
        if (j.ge.JstrV) then
          do i=Istr,Iend
            DC(i,0)=pm_v(i,j)*pn_v(i,j)
            DC(i,1)=DC(i,0)/( Zt_avg1(i,j  )+h(i,j  )
     &                       +Zt_avg1(i,j-1)+h(i,j-1))
          enddo
          if (FIRST_TIME_STEP) then
            do k=1,N
              do i=Istr,Iend
# ifdef VADV_ADAPT_IMP              
                cff = 1.
# else
                cff = 2./(Hz_half(i,j,k)+Hz_half(i,j-1,k))
# endif 
                cff0=DC(i,0)*rv(i,j,k)
   !  &                             /( Zt_avg1(i,j  )+h(i,j  )
   !  &                               +Zt_avg1(i,j-1)+h(i,j-1))

                v(i,j,k,nnew)=( v(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                         +cdt*cff0)*cff
# ifdef M3FAST_C3D_UVSF
                rv_int_nbq(i,j,k) =cff0
#  ifndef M3FAST_COUPLING1
                rv_intt_nbq(i,j,k)=cff0
#  else
                rv_intt_nbq(i,j,k,nstp)=cff0
#  endif
# endif
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
                rv_int_nbq_2d(i,j)= rv_int_nbq_2d(i,j)
     &                           +rv_int_nbq(i,j,k) 
# endif
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5*( Hz(i,j  ,k)
     &                                           +Hz(i,j-1,k))
              enddo
            enddo
          else
            do k=1,N
              do i=Istr,Iend
              
# ifdef VADV_ADAPT_IMP              
                cff = 1.
# else
                cff = 2./(Hz_half(i,j,k)+Hz_half(i,j-1,k))
# endif              
                qdm_nstp = v(i,j,k,nstp)*0.5*( Hz(i,j  ,k)
     &                                        +Hz(i,j-1,k))
                qdm_indx = v(i,j,k,indx)*0.5*( Hz_bak(i  ,j,k)
     &                                        +Hz_bak(i,j-1,k))
                v(i,j,k,nnew)=( cff4*qdm_nstp
     &                         +cff5*qdm_indx
     &                        +cdt*DC(i,0)*(rv(i,j,k)
# ifdef M3FAST_C3D_UVFS
     &                                  +rv_nbq_avg1(i,j,k)
# endif     
     &                                  ))*cff
    
# ifdef M3FAST_C3D_UVSF
#  ifndef M3FAST_COUPLING1
                cff0=((qdm_nstp-qdm_indx) /dt
     &               -DC(i,0)*rv_nbq_avg2(i,j,k))
                rv_int_nbq(i,j,k)=(cff6*cff0
     &                            +cff7*rv_intt_nbq(i,j,k))
#   ifdef MASKING
     &                           *vmask(i,j)   
#   endif     
                rv_intt_nbq(i,j,k)=cff0
#  else
               cff0=DC(i,0)*rv(i,j,k)
#   ifndef M3FAST_ZETAW
!    &            -rvfrc(i,j) *DC(i,1)  
!    &            *(Hz(i,j,k)+Hz(i,j-1,k))                    ! Check Hz_half
#   endif
               rv_int_nbq(i,j,k)= (cff1*cff0
     &                            + cff2*rv_intt_nbq(i,j,k,3-nstp)
     &                            + cff3*rv_intt_nbq(i,j,k,nstp))
!#   ifndef M3FAST_ZETAW
!     &            -rvfrc(i,j) *DC(i,1)  
!     &            *(Hz(i,j,k)+Hz(i,j-1,k))                    ! Check Hz_half
!#  endif
#   ifdef MASKING
     &                           *vmask(i,j)   
#   endif     
              rv_intt_nbq(i,j,k,nstp)=cff0
#  endif /* !M3FAST_COUPLING1 */
        
# ifdef M3FAST_DISSCOUPLING
              cff10=((qdm_nstp-qdm_indx) /dt
     &              -DC(i,0)*rv_nbq_avg2(i,j,k))
              rv_int_nbq(i,j,k) = (1.-alpha_uv) * rv_int_nbq(i,j,k)    !ICI
     &                           +alpha_uv * cff10
# endif      
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
              rv_int_nbq_2d(i,j)= rv_int_nbq_2d(i,j)
     &                           +rv_int_nbq(i,j,k) 
# endif
# endif /* M3FAST_C3D_UVSF */
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              enddo
            enddo               !--> discard DC(:,0)
          endif
     
# if defined M3FAST_COUPLING2D && defined M3FAST_C3D_UVSF
          do i=Istr,Iend
          do k=1,N
          rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)
     &               -rv_int_nbq_2d(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))
     &                               /( Zt_avg1(i,j  )+h(i,j  )
     &                                 +Zt_avg1(i,j-1)+h(i,j-1))
          enddo
          enddo
 
# endif    
        endif 
!
c # ifdef M3FAST
c #  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
c       call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
c      &                        rv_int_nbq(START_2D_ARRAY,1))
c #  endif
c # endif
c # if defined RVTK_DEBUG 
c C$OMP BARRIER
c C$OMP MASTER
c        call check_tab3d(rv_int_nbq,'rv_int_nbq (pre_step3d)','v')
c C$OMP END MASTER     
c # endif  
        
c ! #ifdef OPENACC
c ! !$acc update device(ru_int_nbq,rv_int_nbq,rw_int_nbq)
c ! !$acc& async(sync_ruv_int_nbq)
c ! #endif
c # if defined RVTK_DEBUG 
c C$OMP BARRIER
c C$OMP MASTER
c        call check_tab3d(ru_int_nbq,'ru_int_nbq (pre_step3d) 0','uint')
c        call check_tab3d(rv_int_nbq,'rv_int_nbq (pre_step3d) 0','vint')
c C$OMP END MASTER     
c # endif  
        
!  Compute implicit vertical advection/viscosity
!
# ifdef VADV_ADAPT_IMP
# undef  TRIDIAG_TRA
# undef  TRIDIAG_U
# define TRIDIAG_V
#  include "tridiag_pred.h"
# endif
!
!=====================================================================       =
!
      enddo     !<-- j
!        if (.NOT. (FIRST_TIME_STEP) ) then
!# if defined M3FAST
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
!      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
!     &                        ru_int_nbq(START_2D_ARRAY,1))
!#  endif
!# endif
!# ifdef M3FAST
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
!      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
!     &                        rv_int_nbq(START_2D_ARRAY,1))
!#  endif
!# endif
!        endif

# ifdef M3SLOW_W
!----------------------------------------------------------
!<-- at this point 
!<-- rw            contains internal 3D advection + Coriolis
!<-- rw_nbq_avg1   contains NBQ pressure gradient (+ gravity) 
!                                        + second viscosity 
!----------------------------------------------------------
      do j=Jstr,Jend
       do i=Istr,Iend
        DC(i,0)=pm(i,j)*pn(i,j)
       enddo
!**************************************
       if (FIRST_TIME_STEP) then
!**************************************
        do k=1,N-1
          do i=Istr,Iend
            cff0=DC(i,0)*rw(i,j,k)
            wz(i,j,k,nnew)=( wz(i,j,k,nstp)*(Hz(i,j,k)+Hz(i,j,k+1))
     &                              +2.*dt*cff0) 
     &                     /  (Hz_half(i,j,k)+Hz_half(i,j,k+1)) 
#  ifdef MASKING
     &                           *rmask(i,j)
#  endif    
#  ifdef M3FAST_C3D_WSF
            rw_int_nbq(i,j,k) =cff0
#   ifndef M3FAST_COUPLING1W            
            rw_intt_nbq(i,j,k)=cff0
#   else
            rw_intt_nbq(i,j,k,nstp)=cff0
#   endif
#  endif
            wz(i,j,k,indx)=wz(i,j,k,nstp)*0.5*(Hz(i,j,k  )+
     &                                         Hz(i,j,k+1))
          enddo
        enddo       
        k=N
        do i=Istr,Iend
            cff0=DC(i,0)*rw(i,j,N)
            wz(i,j,N,nnew)=(wz(i,j,N,nstp)*Hz(i,j,N)+2.*dt*cff0)
     &                     /Hz_half(i,j,N) !<-- wz has units m.s-1
#  ifdef MASKING
     &                     *rmask(i,j)
#  endif     
#  ifdef M3FAST_C3D_WSF
            rw_int_nbq(i,j,N) =cff0
#   ifndef M3FAST_COUPLING1W            
            rw_intt_nbq(i,j,N)=cff0
#   else
           rw_intt_nbq(i,j,N,nstp)=cff0
#   endif
#  endif            
            wz(i,j,N,indx)=wz(i,j,N,nstp)*0.5*Hz(i,j,N) !<-- wz(indx) has units kg.m-1.s-1

!== Special treatment for k=0 because the control volume at the top is Hz(1)/2
            rw_int_nbq (i,j,0)=0.
            
#  ifdef NBQ_FREESLIP
            wz(i,j,0,nnew)=  0.5*
     &                       ( u(i,j,1,nnew) * pm_u(i,j) 
     &                                    * ( z_w(i  ,j,0)-z_w(i-1,j,0))
     &                        +u(i+1,j,1,nnew) * pm_u(i+1,j) 
     &                                    * ( z_w(i+1,j,0)-z_w(i  ,j,0))
     &                        +v(i,j,1,nnew) * pm_v(i,j) 
     &                                    * ( z_w(i,j  ,0)-z_w(i,j-1,0)) 
     &                        +v(i,j+1,1,nnew) * pm_v(i,j+1) 
     &                                    * ( z_w(i,j+1,0)-z_w(i,j,0)))
           wz(i,j,0,indx)=wz(i,j,0,nnew)*0.5*Hz(i,j,1) !<-- wz(indx) has units kg.m-1.s-1
!          wz(i,j,0,indx)=wz(i,j,0,ndiffnew)*0.5*Hz(i,j,1) !<-- wz(indx) has units !kg.m-1.s-1
#  else
            wz(i,j,0,nnew)=0.
            wz(i,j,0,indx)=0.
#  endif
        enddo      

!**************************************
      else
!**************************************    
        do k=1,N-1
          do i=Istr,Iend
          
            qdm_nstp = wz(i,j,k,nstp)*(Hz(i,j,k  )+             
     &                                             Hz(i,j,k+1))
            qdm_indx = wz(i,j,k,indx)*(Hz_bak(i,j,k  )+     
     &                                             Hz_bak(i,j,k+1))
            wz(i,j,k,nnew)=( (cff4*qdm_nstp
     &                       +cff5*qdm_indx)  
     &                      +2.*cdt*DC(i,0)                 ! n - 1/2
     &                             *(rw(i,j,k)
#  ifdef M3FAST_C3D_WFS
     &                               +rw_nbq_avg1(i,j,k)  !
#  endif     
     &                     )) / (Hz_half(i,j,k)+Hz_half(i,j,k+1))  !<-- wz has units m.s-1
               
#  ifdef M3FAST_C3D_WSF
#   ifndef M3FAST_COUPLING1W  
           cff0 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                 - DC(i,0)*rw_nbq_avg2(i,j,k))
           rw_int_nbq(i,j,k)=(cff6*cff0
     &                       +cff7*rw_intt_nbq(i,j,k))
#    ifdef MASKING
     &                              *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,k)=cff0
#  else
            cff0=DC(i,0)*rw(i,j,k) ! wzdiff here if necessary
            rw_int_nbq(i,j,k)= (cff1*cff0
     &                        + cff2*rw_intt_nbq(i,j,k,3-nstp)
     &                        + cff3*rw_intt_nbq(i,j,k,nstp))
#   ifdef MASKING
     &                              *rmask(i,j)
#   endif     
           rw_intt_nbq(i,j,k,nstp)=cff0
#  endif /* M3FAST_COUPLING1W */

#   ifdef M3FAST_DISSCOUPLING
           cff10 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                 - DC(i,0)*rw_nbq_avg2(i,j,k))
           rw_int_nbq(i,j,k)=(1.-alpha_w)*rw_int_nbq(i,j,k)   
     &                      +alpha_w * cff10
#   endif
#  endif /* M3FAST_C3D_WSF */
            wz(i,j,k,indx)=wz(i,j,k,nstp)*0.5*(Hz(i,j,k  )+
     &                                         Hz(i,j,k+1)) 
     
          enddo
        enddo
        
        do i=Istr,Iend
!== Special treatment for k=N because the control volume at the top is Hz(N)/2  

           k=N
           qdm_nstp = wz(i,j,N,nstp)*Hz    (i,j,N)
           qdm_indx = wz(i,j,N,indx)*Hz_bak(i,j,N)
           wz(i,j,N,nnew)=((cff4*qdm_nstp
     &                     +cff5*qdm_indx)  
     &                     +2.*cdt*DC(i,0)
     &                        *(rw(i,j,N)
#  ifdef M3FAST_C3D_WFS
     &                         +rw_nbq_avg1(i,j,N)
#  endif     
     &                   ))/ Hz_half(i,j,N)          !<-- wz has units m.s-1
     
#  ifdef M3FAST_C3D_WSF
#   ifndef M3FAST_COUPLING1W  
           cff0 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                      - DC(i,0)*rw_nbq_avg2(i,j,N))
           rw_int_nbq(i,j,N)=(cff6*cff0
     &                       +cff7*rw_intt_nbq(i,j,N))
#    ifdef MASKING
     &                           *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,N)=cff0
#   else
           cff0=DC(i,0)*rw(i,j,N) ! wzdiff here if necessary
           rw_int_nbq(i,j,N)= (cff1*cff0
     &                       +cff2*rw_intt_nbq(i,j,N,3-nstp)
     &                       +cff3*rw_intt_nbq(i,j,N,nstp))
#    ifdef MASKING
     &                           *rmask(i,j)
#    endif     
           rw_intt_nbq(i,j,N,nstp)=cff0
#   endif /* M3FAST_COUPLING1W */
#    ifdef M3FAST_DISSCOUPLING
           cff10 = ((qdm_nstp-qdm_indx)*0.5/dt
     &                      - DC(i,0)*rw_nbq_avg2(i,j,N))
           rw_int_nbq(i,j,N)=  (1.-alpha_w)*rw_int_nbq(i,j,N)     ! ICI
     &                         +alpha_w*cff10
#    endif
#  endif /* M3FAST_C3D_WSF */

           wz(i,j,N,indx)=wz(i,j,N,nstp)*0.5*Hz(i,j,N)   !<-- wz(indx) has units kg.m-1.s-1
     
!== Special treatment for k=0 because the control volume at the bottom
!    is Hz(1)/2  
#  ifdef NBQ_FREESLIP
            wz(i,j,0,nnew)=  0.5*
     &                       ( u(i,j,1,nnew) * pm_u(i,j) 
     &                                    * ( z_w(i  ,j,0)-z_w(i-1,j,0))
     &                        +u(i+1,j,1,nnew) * pm_u(i+1,j) 
     &                                    * ( z_w(i+1,j,0)-z_w(i  ,j,0))
     &                        +v(i,j,1,nnew) * pm_v(i,j) 
     &                                    * ( z_w(i,j  ,0)-z_w(i,j-1,0)) 
     &                        +v(i,j+1,1,nnew) * pm_v(i,j+1) 
     &                                    * ( z_w(i,j+1,0)-z_w(i,j,0)))
          wz(i,j,0,indx)=wz(i,j,0,nnew)*0.5*Hz(i,j,1) !<-- wz(indx) has units kg.m-1.s-1
!         wz(i,j,0,indx)=wz(i,j,0,ndiffnew)*0.5*Hz(i,j,1) !<-- wz(indx) has units !kg.m-1.s-1
#  else
            wz(i,j,0,nnew)=0.
            wz(i,j,0,indx)=0.
#  endif

        enddo    
!==                  
      endif
  
      enddo  
!==  
# endif /* M3FAST_C3D_WSF */   
  
# ifdef M3FAST_C3D_UVFS
#  undef ru_nbq_avg1 
#  undef rv_nbq_avg1
# endif
# ifdef M3FAST_C3D_WFS
#  undef rw_nbq_avg1
# endif
!
!--------------------------------------------------------------------
!

# if defined TS_HADV_TEST && ( defined DIAGONAL_ADV || defined SOLID_BODY_ROT )
      do k=1,N
        do j=Jstr,Jend
          do i=IstrR,IendR
            u(i,j,k,nnew) = u(i,j,k,nstp)
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            v(i,j,k,nnew) = v(i,j,k,nstp)
          enddo
        enddo
      enddo
# endif
# if defined TS_HADV_TEST && defined SOLID_BODY_PER
      cff = cos(2.*pi*(time+dt/2.)/(float(ntimes)*dt))
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            u(i,j,k,nnew) = u(i,j,k,nstp) * cff
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            v(i,j,k,nnew) = v(i,j,k,nstp) * cff
          enddo
        enddo
      enddo
# endif
!
!--------------------------------------------------------------------
!
!======================================================================
! Set PHYSICAL lateral boundary conditions for tracers.
!======================================================================
!
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#  ifdef THREE_GHOST_POINTS_TS
        call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                               t(START_2D_ARRAY,1,nnew,itrc))
#  else
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
#  endif
# endif
      enddo

      call u3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
      call v3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
# ifdef M3SLOW_W
      call w3dbc_tile (Istr,Iend,Jstr,Jend, WORK)      
# endif
!
!======================================================================
! Coupling, include ghost points associated with PHYSICAL
! boundaries ONLY. Do not touch periodic ghost points or
! internal computational margins (MPI code). 
!======================================================================
!
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!======================================================================
! Set PHYSICAL lateral boundary conditions for momentum.
!======================================================================
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#  ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v(START_2D_ARRAY,1,nnew))
#   ifdef M3SLOW_W   
      call exchange_w3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             wz(START_2D_ARRAY,0,nnew))
#   endif
#   ifdef M3FAST_C3D_UVSF
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &            ru_int_nbq(START_2D_ARRAY,1))   
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &            rv_int_nbq(START_2D_ARRAY,1))   
#   endif
#   ifdef M3FAST_C3D_WSF
      call exchange_w3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &            rw_int_nbq(START_2D_ARRAY,0))   
#   endif
#  else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,nnew))
#   ifdef M3SLOW_W   
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        wz(START_2D_ARRAY,0,nnew))   
#   endif
#   ifdef M3FAST_C3D_UVSF
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &            ru_int_nbq(START_2D_ARRAY,1))   
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &            rv_int_nbq(START_2D_ARRAY,1))   
#   endif
#   ifdef M3FAST_C3D_WSF
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &            rw_int_nbq(START_2D_ARRAY,0))   
#   endif
#  endif

!      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
!     &            ru_int_nbq(START_2D_ARRAY,1))   
!      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
!     &            rv_int_nbq(START_2D_ARRAY,1))   
# endif

!
!======================================================================
! Prepare to start two-dimensional time stepping:
! set the initial values of the fast-time-step free-surface
! field to its fast-time-averaged values corresponding
! to the time step n (nstp).
!======================================================================
!
# ifndef M2FILTER_NONE 
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do j=JU_RANGE
        do i=IV_RANGE
          zeta(i,j,knew)=Zt_avg1(i,j)
        enddo
      enddo
       call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(START_2D_ARRAY,knew))
#  else
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta(i,j,knew)=Zt_avg1(i,j)
        enddo
      enddo     
#  endif
# endif
!
# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE
#else
      subroutine pre_step3d_2_empty
#endif
      return
      end
