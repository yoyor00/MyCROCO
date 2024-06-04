! $Id: set_depth.F 1484 2014-03-17 14:01:55Z rblod $
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
#if defined SOLVE3D || defined WET_DRY

      subroutine set_depth (tile)
      implicit none
# include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      call set_depth_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_depth_tile (Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
! Create S-coordinate system: based on model topography h(i,j),
! fast-time-averaged free-surface field and vertical coordinate
! transformation metrics compute evolving depths of of the three-
! dimensional model grid (z_r,z_w) and vertical heights of model
! grid boxes.
!----------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    cff_r,cff1_r,cff_w,cff1_w, z_r0,z_w0, zetatmp,htmp
# include "param.h"
# include "grid.h"
# ifdef AGRIF
#  include "zoom.h"
# endif
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "scalars.h"
# ifdef MRL_WCI
#  include "forces.h"
# endif
# ifdef M3FAST
#  include "nbq.h"
# endif
      real eps
      parameter (eps=1.e-14)
!
# include "compute_auxiliary_bounds.h"
!
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE JstrR,JendR
# endif
!
!----------------------------------------------------------------------
!  Initialization of free surface zeta arrays and wet/dry mask
!
!  During initialization and/or restart: copy initial free surface
!  field into array for holding fast-time averaged free surface.
!----------------------------------------------------------------------
!
      if (iic.eq.0) then
        do j=J_RANGE
          do i=I_RANGE
# ifdef WET_DRY
            wetdry(i,j)=1.
#  ifdef AGRIF
            rmask_childs(i,j)=1.
#  endif
            IF (h(i,j).eq.0.) THEN
              h(i,j)=eps
            END IF
            IF (zeta(i,j,1) .lt. Dcrit(i,j)-h(i,j)) THEN
              zeta(i,j,1)=Dcrit(i,j)-h(i,j)
              zeta(i,j,2)=zeta(i,j,1)
              zeta(i,j,3)=zeta(i,j,1)
              zeta(i,j,4)=zeta(i,j,1)
              wetdry(i,j)=0.
            ENDIF
            rmask_wet(i,j)=wetdry(i,j)
# endif
# ifdef SOLVE3D
#  ifdef NEW_S_COORD
            hinv(i,j)=1./(abs(h(i,j))+hc)
#  else
            hinv(i,j)=1./h(i,j)
#  endif
            Zt_avg1(i,j)=zeta(i,j,knew)  ! knew=1 unless modified
                                         ! in ana_initial (JET)
            hgrid(i,j)=h(i,j)  ! save in grid bathy for penalization
# endif /* SOLVE3D */
          enddo
        enddo

# ifdef NBQ_MASS
        do k=1,N
          do j=J_RANGE
            do i=I_RANGE
              rho_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, hinv)
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, hgrid)
#  ifdef WET_DRY
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, h)
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, wetdry)
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, rmask_wet)
#   ifdef AGRIF
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, rmask_childs)
#   endif
#  endif
#  if defined M3FAST || defined WET_DRY
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(START_2D_ARRAY,1))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(START_2D_ARRAY,2))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(START_2D_ARRAY,3))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(START_2D_ARRAY,4))
#  endif
# endif
      endif   !  iic==0  <-- end initialization
!
!----------------------------------------------------------------------
!  Set vertical grid
!----------------------------------------------------------------------
!
# ifdef SOLVE3D
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, Zt_avg1)
#  endif

      do j=J_RANGE              !!! WARNING: Setting must be consistent
        do i=I_RANGE            !!!          with omega.F
          z_w(i,j,0)=-hgrid(i,j)
        enddo

        do k=1,N,+1
#  ifdef NEW_S_COORD
          cff_w =hc*sc_w(k)
          cff_r =hc*sc_r(k)
          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)
#  else
          cff_w =hc*(sc_w(k)-Cs_w(k))
          cff_r =hc*(sc_r(k)-Cs_r(k))
          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)
#  endif
          do i=I_RANGE
#  ifdef M3FAST
            zetatmp=zeta(i,j,knew)
#  else
            zetatmp=Zt_avg1(i,j)
#  endif
#  if defined MASKING && !defined WET_DRY
            zetatmp=zetatmp*rmask(i,j)
#  endif
#  ifdef NEW_S_COORD
            htmp=abs(hgrid(i,j))
#  else
            htmp=hgrid(i,j)
#  endif
            z_w0=cff_w+cff1_w*htmp
            z_r0=cff_r+cff1_r*htmp
#  ifdef NEW_S_COORD
            z_w(i,j,k)=z_w0*hgrid(i,j)*hinv(i,j)+zetatmp
     &                                       *(1.+z_w0*hinv(i,j))
            z_r(i,j,k)=z_r0*hgrid(i,j)*hinv(i,j)+zetatmp
     &                                       *(1.+z_r0*hinv(i,j))
#  else
            z_w(i,j,k)=z_w0+zetatmp*(1.+z_w0*hinv(i,j))
            z_r(i,j,k)=z_r0+zetatmp*(1.+z_r0*hinv(i,j))
#  endif

#  ifdef NBQ
#   ifdef NBQ_GRID_SLOW
            if (LAST_FAST_STEP)  Hz_bak(i,j,k)=Hz(i,j,k)
#   else
            if (FIRST_FAST_STEP) Hz_bak(i,j,k)=Hz(i,j,k)
#   endif
#   ifdef NBQ_HZ_PROGNOSTIC
            Hz_bak2(i,j,k)=Hz(i,j,k)
#   endif
#  else
            Hz_bak(i,j,k)=Hz(i,j,k)
#  endif
            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)

#  ifdef NBQ_MASS
            Hzr(i,j,k)=Hz(i,j,k)
            Hz(i,j,k) =Hz(i,j,k)*(1.+rho_grd(i,j,k))
     &                              +rho_nbq(i,j,k)
#  endif
          enddo
        enddo
      enddo

#  if defined M3FAST && ! defined EXACT_RESTART
      if (iic.eq.0 .or. (FIRST_TIME_STEP .and. iif.eq.1)) then
        do k=1,N,+1
          do j=J_RANGE
            do i=I_RANGE
              Hz_bak(i,j,k)=Hz(i,j,k)
            enddo
          enddo
        enddo
      endif
#  endif
 
#  undef I_RANGE
#  undef J_RANGE

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_w(START_2D_ARRAY,0))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_r(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz_bak(START_2D_ARRAY,1))
#   ifdef NBQ_MASS
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hzr(START_2D_ARRAY,1))
#   endif
#  endif
      return
      end

#  ifdef MORPHODYN
!======================================================================
!
!                    subroutine set_depth_morphodyn
!
!======================================================================
!

      subroutine set_depth_morphodyn (tile)
      implicit none
#   include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#   include "compute_tile_bounds.h"
      call set_depth_morphodyn_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_depth_morphodyn_tile (Istr,Iend,Jstr,Jend)

      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    cff_r,cff1_r,cff_w,cff1_w, z_r0,z_w0, zetatmp,htmp
#   include "param.h"
#   include "grid.h"
#   ifdef AGRIF
#    include "zoom.h"
#   endif
#   include "ocean2d.h"
#   include "ocean3d.h"
#   include "coupling.h"
#   include "scalars.h"
#   ifdef MRL_WCI
#    include "forces.h"
#   endif
#   ifdef M3FAST
#    include "nbq.h"
#   endif
      real eps
      parameter (eps=1.e-14)
!
#   include "compute_auxiliary_bounds.h"
!
#   ifdef EW_PERIODIC
#    define I_RANGE Istr,Iend
#   else
#    define I_RANGE IstrR,IendR
#   endif
#   ifdef NS_PERIODIC
#    define J_RANGE Jstr,Jend
#   else
#    define J_RANGE JstrR,JendR
#   endif

      do j=J_RANGE              !!! WARNING: Setting must be consistent
        do i=I_RANGE            !!!          with omega.F
          h(i,j)=h(i,j)+dh(i,j)
          z_w(i,j,0)=-h(i,j)
        enddo

        do k=1,N,+1
#   ifdef NEW_S_COORD
          cff_w =hc*sc_w(k)
          cff_r =hc*sc_r(k)
          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)
#   else
          cff_w =hc*(sc_w(k)-Cs_w(k))
          cff_r =hc*(sc_r(k)-Cs_r(k))
          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)
#   endif
          do i=I_RANGE
#   ifdef M3FAST
            zetatmp=zeta(i,j,knew)
#   else
            zetatmp=Zt_avg1(i,j)
#   endif
#  ifdef NEW_S_COORD
            htmp=abs(h(i,j))
            hinv(i,j)=1./(htmp+hc)
#  else
            htmp=h(i,j)
            hinv(i,j)=1./htmp
#  endif
#   if defined MASKING && !defined WET_DRY
            zetatmp=zetatmp*rmask(i,j)
#   endif
            z_w0=cff_w+cff1_w*htmp
            z_r0=cff_r+cff1_r*htmp
#   ifdef NEW_S_COORD
            z_w(i,j,k)=z_w0*h(i,j)*hinv(i,j)+zetatmp
     &                                       *(1.+z_w0*hinv(i,j))
            z_r(i,j,k)=z_r0*h(i,j)*hinv(i,j)+zetatmp
     &                                       *(1.+z_r0*hinv(i,j))
#   else
            z_w(i,j,k)=z_w0+zetatmp*(1.+z_w0*hinv(i,j))
            z_r(i,j,k)=z_r0+zetatmp*(1.+z_r0*hinv(i,j))
#   endif
            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
#   ifdef NBQ_MASS
            Hzr(i,j,k)=Hz(i,j,k)
            Hz(i,j,k) =Hz(i,j,k)*(1.+rho_grd(i,j,k))
     &                              +rho_nbq(i,j,k)
#   endif
          enddo
        enddo
      enddo

#   undef I_RANGE
#   undef J_RANGE

#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_w(START_2D_ARRAY,0))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_r(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, h)
#    ifdef NBQ_MASS
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hzr(START_2D_ARRAY,1))
#    endif
#   endif
      return
      end
#  endif /* MORPHODYN */

#  ifdef POROSITY
!======================================================================
!
!                    subroutine set_depth_porosity
!
!======================================================================
!
      subroutine set_depth_porosity (tile)
      implicit none
#   include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#   include "compute_tile_bounds.h"
      call set_depth_porosity_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_depth_porosity_tile (Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
! Create S-coordinate system: based on model topography h(i,j),
! fast-time-averaged free-surface field and vertical coordinate
! transformation metrics compute evolving depths of of the three-
! dimensional model grid (z_r,z_w) and vertical heights of model
! grid boxes.
!----------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    cff_r,cff1_r,cff_w,cff1_w, z_r0,z_w0, zetatmp
#  include "param.h"
#  include "grid.h"
#  ifdef AGRIF
#   include "zoom.h"
#  endif
#  include "ocean2d.h"
#  include "ocean3d.h"
#  include "coupling.h"
#  include "scalars.h"
#  include "penalization.h"
      real,parameter :: smooth_crit = 2.77259
      real :: z_decal
      real,parameter :: dzpena_min=0.25
      real :: pena_ref
      real z_u,z_bot,z_top,z_v,z_h,z_right,z_r1
      real :: s1
      real z_lastindex(GLOBAL_2D_ARRAY)
      real, parameter :: alpha_poro = 0.02
      real FZ3(0:N)
      real mindz(N_poro_discrete)
      integer :: k2
      real cff4,cff5,cff6,cff7,cff8
      integer :: i1,j1,k1
      real eps
      parameter (eps=1.e-14)
!
#  include "compute_auxiliary_bounds.h"
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!----------------------------------------------------------------------
!  Compute new Hz accounting for porosity
!----------------------------------------------------------------------
!
!  First compute time-independent structure function poro_discrete
!
      if (iic.eq.0) then
        call set_depth_poro_init_tile (Istr,Iend,Jstr,Jend)
      endif
!
!  Save full Hz
!
      do j=J_RANGE
        do i=I_RANGE
          do k=1,N,+1
            Hzf(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
          enddo
        enddo
      enddo
!
!  Compute time-dependent function pororatio and new Hz
!
      do j=J_RANGE 
        do i=I_RANGE
          FZ3(0)=0.
          do k=1,N
            if (z_w(i,j,k)<0.)  then
              k2=int(N_poro_discrete*(1.+z_w(i,j,k)/hgrid(i,j)))
              FZ3(k)=z_int_poro(i,j,k2)
     &               +(z_w(i,j,k)-z_poro_discrete(i,j,k2))
     &                             *poro_discrete(i,j,k2+1)
            else
              FZ3(k)=z_int_poro(i,j,N_poro_discrete)+z_w(i,j,k)
            endif
          enddo
          do k=1,N
            Hz(i,j,k)=FZ3(k)-FZ3(k-1)
            pororatio(i,j,k)=Hz(i,j,k)/(z_w(i,j,k)-z_w(i,j,k-1))
          enddo
        enddo
      enddo

      if (iic.eq.0) then
        do k=1,N,+1
          do j=J_RANGE
            do i=I_RANGE
              Hz_bak(i,j,k)=Hz(i,j,k)
            enddo
          enddo
        enddo
      endif

!
!  Compute index arrays of penalization region
!
      do j=J_RANGE
        do i=I_RANGE
          last_index(i,j) = 0
          last_index_w(i,j)=0
        enddo
      enddo
      do j=J_RANGE 
        do i=I_RANGE
          do k=N,1,-1
            if (z_r(i,j,k)>=-h_pena_target(i,j)) then
              last_index(i,j) = k
            else
              exit
            endif
          enddo
          last_index(i,j) = last_index(i,j) - 1
        enddo
      enddo
      do j=J_RANGE 
        do i=I_RANGE
          do k=N,1,-1
            if (z_w(i,j,k)>=-h_pena_target(i,j)) then
              last_index_w(i,j) = k
            else
              exit
            endif
          enddo
        enddo
      enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hzf(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        pororatio(START_2D_ARRAY,1))
      if (iic.eq.0) then
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz_bak(START_2D_ARRAY,1))
      endif
      z_lastindex = last_index
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, z_lastindex)
      last_index = nint(z_lastindex)
      z_lastindex = last_index_w
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, z_lastindex)
      last_index_w = nint(z_lastindex)
#  endif

#  undef I_RANGE
#  undef J_RANGE

      return
      end
!======================================================================
!
!                    subroutine set_depth_poro_init
!
!======================================================================
!
      subroutine set_depth_poro_init (tile)
      implicit none
#   include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#   include "compute_tile_bounds.h"
      call set_depth_poro_init_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_depth_poro_init_tile (Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
! Create S-coordinate system: based on model topography h(i,j),
! fast-time-averaged free-surface field and vertical coordinate
! transformation metrics compute evolving depths of of the three-
! dimensional model grid (z_r,z_w) and vertical heights of model
! grid boxes.
!----------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    cff_r,cff1_r,cff_w,cff1_w, z_r0,z_w0, zetatmp
#  include "param.h"
#  include "grid.h"
#  ifdef AGRIF
#   include "zoom.h"
#  endif
#  include "ocean2d.h"
#  include "ocean3d.h"
#  include "coupling.h"
#  include "scalars.h"
#  include "penalization.h"
      real,parameter :: smooth_crit = 2.77259
      real :: z_decal
      real,parameter :: dzpena_min=0.25
      real :: pena_ref
      real z_u,z_bot,z_top,z_v,z_h,z_right,z_r1
      real :: s1
      real z_lastindex(GLOBAL_2D_ARRAY)
      real, parameter :: alpha_poro = 0.02
      real z_w_tmp(GLOBAL_2D_ARRAY,0:N)
      real FZ3(0:N)
      real mindz(N_poro_discrete)
      integer :: k2
      real cff4,cff5,cff6,cff7,cff8
      integer :: i1,j1,k1
      real eps
      parameter (eps=1.e-14)
!
#  include "compute_auxiliary_bounds.h"
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!----------------------------------------------------------------------
!  Initialization of free surface zeta arrays and wet/dry mask
!
!  During initialization and/or restart: copy initial free surface 
!  field into array for holding fast-time averaged free surface.
!----------------------------------------------------------------------
!
      z_decal = 2.19722/smooth_crit
!
! COMPUTE VERTICAL GRID WITH ZETA = 0 to get H_PORO
!
      do j=J_RANGE
        do i=I_RANGE
          z_w_tmp(i,j,0)=-hgrid(i,j)
          h(i,j)=0.
        enddo

        do k=1,N,+1
#  ifdef NEW_S_COORD
          cff_w =hc*sc_w(k)
          cff1_w=Cs_w(k)
#  else
          cff_w =hc*(sc_w(k)-Cs_w(k))
          cff1_w=Cs_w(k)
#  endif
          do i=I_RANGE
            zetatmp=0.
#  if defined MASKING && !defined WET_DRY
            zetatmp=zetatmp*rmask(i,j)
#  endif
            z_w0=cff_w+cff1_w*hgrid(i,j)        
#  ifdef NEW_S_COORD
            z_w_tmp(i,j,k)=z_w0*hgrid(i,j)*hinv(i,j)+zetatmp
     &                                     *(1.+z_w0*hinv(i,j))
#  else
            z_w_tmp(i,j,k)=z_w0+zetatmp*(1.+z_w0*hinv(i,j))
#  endif
          enddo
        enddo
      enddo

      do j=J_RANGE
        do i=I_RANGE
          z_poro_discrete(i,j,0) = -hgrid(i,j)
          z_poro_discrete(i,j,N_poro_discrete) = 0.
          do k=1,N_poro_discrete-1
            z_poro_discrete(i,j,k)=-hgrid(i,j) + k * hgrid(i,j)/N_poro_discrete
          enddo
          mindz(:) = huge(1.)
          do k=1,N_poro_discrete
            do k2=1,N
              if ((z_w_tmp(i,j,k2  ) >= z_poro_discrete(i,j,k-1)).AND.
     &            (z_w_tmp(i,j,k2-1) <= z_poro_discrete(i,j,k))) then
                mindz(k)=min(mindz(k),z_w_tmp(i,j,k2)-z_w_tmp(i,j,k2-1))
              endif
            enddo
          enddo
          do k=1,N_poro_discrete
            z_u = 0.5*(z_poro_discrete(i,j,k-1)+z_poro_discrete(i,j,k))
            z_right=(-h_pena_target(i,j)-z_u)/mindz(k)
            z_h=-hgrid(i,j)
            poro_discrete(i,j,k)=1.-0.5*(1.+tanh(smooth_crit*0.5*
     &                                        (z_right-z_decal)))
            poro_discrete(i,j,k)=max(poro_discrete(i,j,k),alpha_poro)
            poro_discrete(i,j,k)=max(poro_discrete(i,j,k),
     &                                     dzpena_min/(0.9*mindz(k)))
            poro_discrete(i,j,k)=min(poro_discrete(i,j,k),1.)
            if (abs(z_h+h_pena_target(i,j))<=0.1) then
              poro_discrete(i,j,k)=1.
            endif
          enddo

          poro_discrete(i,j,N_poro_discrete) = 1.

          z_int_poro(i,j,0) = 0.
          do k=1,N_poro_discrete
            z_int_poro(i,j,k)=z_int_poro(i,j,k-1)+
     &              (z_poro_discrete(i,j,k)-z_poro_discrete(i,j,k-1))*
     &                                            poro_discrete(i,j,k)
          enddo 

          FZ3(0)=0.         
          do k=1,N
            if (z_w_tmp(i,j,k)<0.)  then
              k2=int(N_poro_discrete*(1.+z_w_tmp(i,j,k)/hgrid(i,j)))
              FZ3(k)=z_int_poro(i,j,k2)+
     &                 (z_w_tmp(i,j,k)-z_poro_discrete(i,j,k2))
     &                                  *poro_discrete(i,j,k2+1)
            else
              FZ3(k)=z_int_poro(i,j,N_poro_discrete)+z_w_tmp(i,j,k)
            endif
          enddo  
          do k=1,N
            Hz(i,j,k)=FZ3(k)-FZ3(k-1)
            h(i,j)=h(i,j)+Hz(i,j,k)
          enddo      

        enddo ! END i   
      enddo ! END j

      


#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend, Hz(START_2D_ARRAY,1))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, h)
#  endif 

      return
      end
#  endif /* POROSITY */

#  ifdef PERMEABILITY
!======================================================================
!
!                    subroutine set_depth_permeability
!
!======================================================================
!
      subroutine set_depth_permeability (tile)
      implicit none
#   include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#   include "compute_tile_bounds.h"
      call set_depth_permeability_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_depth_permeability_tile (Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
! Create S-coordinate system: based on model topography h(i,j),
! fast-time-averaged free-surface field and vertical coordinate
! transformation metrics compute evolving depths of of the three-
! dimensional model grid (z_r,z_w) and vertical heights of model
! grid boxes.
!----------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    cff_r,cff1_r,cff_w,cff1_w, z_r0,z_w0, zetatmp
#  include "param.h"
#  include "grid.h"
#  ifdef AGRIF
#   include "zoom.h"
#  endif
#  include "ocean2d.h"
#  include "ocean3d.h"
#  include "coupling.h"
#  include "scalars.h"
#  include "penalization.h"
      real,parameter :: smooth_crit = 2.77259
      real :: z_decal
      real,parameter :: dzpena_min=0.25
      real :: pena_ref
      real z_u,z_bot,z_top,z_v,z_h,z_right,z_r1
      real :: s1
      real cff4,cff5,cff6,cff7,cff8
      integer :: i1,j1,k1
      real eps
      parameter (eps=1.e-14)
!
#  include "compute_auxiliary_bounds.h"
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!----------------------------------------------------------------------
!  Initialization of free surface zeta arrays and wet/dry mask
!
!  During initialization and/or restart: copy initial free surface 
!  field into array for holding fast-time averaged free surface.
!----------------------------------------------------------------------
!

      z_decal = 2.19722/smooth_crit
      epsilon_pena=dt/2.

      do k=N,1,-1
        do j=JstrR,JendR
          do i=IstrR,IendR
            pena_r(i,j,k)=0.
          enddo
        enddo
      enddo

      do j=JstrR,JendR
        do i=IstrR,IendR
          z_u=-h_pena_target(i,j)
          z_h=-hgrid(i,j)
          if (abs(z_h-z_u)<=-0.1) then
            pena_r(i,j,:)=0.
          else
            cff5=-HUGE(1.)
            do k=N,1,-1 ! compute velocity above floor
              cff7=(z_r(i,j,k)-z_u)/Zob(i,j)
              if (cff7>=1.) then
                cff4=0.5*(v(i  ,j,k,nrhs)+v(i  ,j+1,k,nrhs))
                cff6=0.5*(u(i  ,j,k,nrhs)+u(i+1,j ,k,nrhs))
                cff5=SQRT(cff6*cff6+cff4*cff4)
                k1 = k  ! Vertical index of last level
                cff4=min(Cdb_max,max((vonKar/LOG(cff7))**2,Cdb_min))
                pena_ref=dt*cff4*cff5/(0.5*(Hz(i-1,j,k)+Hz(i,j,k)))
                pena_ref=max(min(pena_ref,1.-1.D-4),1.D-4)
                z_bot=z_w(i,j,k-1)
                z_top=z_w(i,j,k)
                z_r1=(z_u-(z_r(i,j,k)))/(z_top-z_bot)
              else 
                exit
              endif
            enddo
            s1=2.*atanh(1.-2.*pena_ref)/(z_decal-z_r1)
            do k=N,1,-1
              z_bot=z_w(i,j,k-1)
              z_top=z_w(i,j,k)
              z_right=(z_u-(z_r(i,j,k)))/(z_top-z_bot)
              if (z_right<1.0001*z_r1) then
                pena_r(i,j,k)=0.
              else
                pena_r(i,j,k)=pena_ref+(1.-pena_ref)*
     &                       (tanh(smooth_crit*0.5*(z_right-z_r1)))
              endif
            enddo
          endif
        enddo
      enddo
       
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        pena_r(START_2D_ARRAY,1))
#  endif

      do k=1,N
        do j=JstrR,JendR
          do i=Istr,Iend+1
            pena_u(i,j,k)=max(pena_r(i-1,j,k),pena_r(i,j,k))
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend+1
          do i=IstrR,IendR
            pena_v(i,j,k)=max(pena_r(i,j-1,k),pena_r(i,j,k))
          enddo
        enddo
      enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        pena_u(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        pena_v(START_2D_ARRAY,1))
#  endif

#  undef I_RANGE
#  undef J_RANGE

      return
      end
#  endif /* PERMEABILITY */
!
!======================================================================
!
!                       subroutine set_HUV
!
!======================================================================
!
      subroutine set_HUV (tile)
      implicit none
#  include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#  include "compute_tile_bounds.h"
      call set_HUV_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_HUV_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#  include "param.h"
#  include "grid.h"
#  include "ocean3d.h"
#  include "scalars.h"
!
#  ifdef MRL_WCI
#   include "forces.h"
#  endif
#  include "compute_auxiliary_bounds.h"
!
#  ifdef EW_PERIODIC
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif

#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif

      do k=1,N
        do j=JU_RANGE
          do i=IU_RANGE
            Huon(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
     &                                       *( u(i,j,k,nrhs)
#  ifdef MRL_WCI
     &                                           + ust(i,j,k)
#  endif
     &                                                      )
#  if defined MRL_WCI && defined MASKING
            Huon(i,j,k)=Huon(i,j,k)*umask(i,j)
#  endif

#  if defined UV_VIS4 && defined UV_MIX_GEO
            z_u(i,j,k)=0.5*(z_r(i,j,k)+z_r(i-1,j,k))
            dz_u(i,j,k)=z_r(i,j,k)-z_r(i-1,j,k)
#  endif
          enddo
        enddo
        do j=JV_RANGE
          do i=IV_RANGE
            Hvom(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
     &                                       *( v(i,j,k,nrhs)
#  ifdef MRL_WCI
     &                                           + vst(i,j,k)
#  endif
     &                                                      )
#  if defined MRL_WCI && defined MASKING
            Hvom(i,j,k)=Hvom(i,j,k)*vmask(i,j)
#  endif
#  if defined UV_VIS4 && defined UV_MIX_GEO
            z_v(i,j,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
            dz_v(i,j,k)=z_r(i,j,k)-z_r(i,j-1,k)
#  endif
          enddo
        enddo
      enddo
#  undef IU_RANGE
#  undef JU_RANGE
#  undef IV_RANGE
#  undef JV_RANGE
!
! Exchange periodic boundaries, if so prescribed.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             Huon(START_2D_ARRAY,1))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             Hvom(START_2D_ARRAY,1))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        Huon(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hvom(START_2D_ARRAY,1))
#   endif
#   if defined UV_VIS4 && defined UV_MIX_GEO
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_u(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_v(START_2D_ARRAY,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        dz_u(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        dz_v(START_2D_ARRAY,1))
#   endif
#  endif
      return
      end
!
!======================================================================
!
!                       subroutine set_HUV1
!
!======================================================================
!
      subroutine set_HUV1 (tile)
      implicit none
#  include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#  include "compute_tile_bounds.h"
      call set_HUV1_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_HUV1_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#  include "param.h"
#  include "grid.h"
#  include "ocean3d.h"
#  include "scalars.h"
#  ifdef MRL_WCI
#   include "forces.h"
#  endif
!
#  include "compute_auxiliary_bounds.h"
!
#  ifdef EW_PERIODIC
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif

#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif

      do k=1,N
        do j=JU_RANGE
          do i=IU_RANGE
            Huon(i,j,k)=0.25*(3.*(Hz(i,j,k)+Hz(i-1,j,k))
     &                   -Hz_bak(i,j,k)-Hz_bak(i-1,j,k))
     &                        *on_u(i,j)*( u(i,j,k,nrhs)
#  ifdef MRL_WCI
     &                                      + ust(i,j,k)
#  endif
     &                                                 )
#  if defined MRL_WCI && defined MASKING
            Huon(i,j,k)=Huon(i,j,k)*umask(i,j)
#  endif

          enddo
        enddo
        do j=JV_RANGE
          do i=IV_RANGE
            Hvom(i,j,k)=0.25*( 3.*(Hz(i,j,k)+Hz(i,j-1,k))
     &                    -Hz_bak(i,j,k)-Hz_bak(i,j-1,k))
     &                         *om_v(i,j)*( v(i,j,k,nrhs)
#  ifdef MRL_WCI
     &                                       + vst(i,j,k)
#  endif
     &                                                  )
#  if defined MRL_WCI && defined MASKING
            Hvom(i,j,k)=Hvom(i,j,k)*vmask(i,j)
#  endif
          enddo
        enddo
      enddo
#  undef IU_RANGE
#  undef JU_RANGE
#  undef IV_RANGE
#  undef JV_RANGE
!
! Exchange periodic boundaries, if so prescribed.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             Huon(START_2D_ARRAY,1))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             Hvom(START_2D_ARRAY,1))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        Huon(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hvom(START_2D_ARRAY,1))
#   endif
#  endif
      return
      end
!
!======================================================================
!
!                       subroutine set_HUV2
!
!======================================================================
!
      subroutine set_HUV2 (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
#  include "param.h"
#  include "private_scratch.h"
#  include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call set_HUV2_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd),
     &                                         A2d(1,2,trd))
      return
      end

      subroutine set_HUV2_tile (Istr,Iend,Jstr,Jend,DC,FC)
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#  include "param.h"
      real DC(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     FC(PRIVATE_1D_SCRATCH_ARRAY,0:N)
#  include "grid.h"
#  include "ocean3d.h"
#  include "scalars.h"
#  ifdef MRL_WCI
#   include "forces.h"
#  endif
#  include "coupling.h"
!
#  include "compute_auxiliary_bounds.h"
!
#  ifdef EW_PERIODIC
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif

#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif

      do j=JU_RANGE
        do i=IU_RANGE
          DC(i,0)=0.
          FC(i,0)=0.
        enddo
        do k=1,N
          do i=IU_RANGE
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*u(i,j,k,nrhs)
          enddo
        enddo
        do i=IU_RANGE
          FC(i,0)=(FC(i,0)-DU_avg2(i,j))/DC(i,0)
#  ifdef MRL_WCI
     &                               +ust2d(i,j)
#  endif
        enddo
        do k=1,N
          do i=IU_RANGE
            u(i,j,k,nrhs)=(u(i,j,k,nrhs)-FC(i,0))
#  ifdef MASKING
     &                                *umask(i,j)
#   ifdef MRL_WCI
     &               +ust(i,j,k)*(umask(i,j)-1.0)
#   endif
#  endif
            Huon(i,j,k)=DC(i,k)*( u(i,j,k,nrhs)
#  ifdef MRL_WCI
     &                             + ust(i,j,k)
#  endif
     &                                        )
          enddo
        enddo
      enddo

      do j=JV_RANGE
        do i=IV_RANGE
          DC(i,0)=0.
          FC(i,0)=0.
        enddo
        do k=1,N
          do i=IV_RANGE
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*v(i,j,k,nrhs)
          enddo
        enddo
        do i=IV_RANGE
          FC(i,0)=(FC(i,0)-DV_avg2(i,j))/DC(i,0)
#  ifdef MRL_WCI
     &                               +vst2d(i,j)
#  endif
        enddo
        do k=1,N
          do i=IV_RANGE
            v(i,j,k,nrhs)=(v(i,j,k,nrhs)-FC(i,0))
#  ifdef MASKING
     &                                *vmask(i,j)
#   ifdef MRL_WCI
     &               +vst(i,j,k)*(vmask(i,j)-1.0)
#   endif
#  endif
            Hvom(i,j,k)=DC(i,k)*( v(i,j,k,nrhs)
#  ifdef MRL_WCI
     &                             + vst(i,j,k)
#  endif
     &                                        )
          enddo
        enddo
      enddo

#  undef IU_RANGE
#  undef JU_RANGE
#  undef IV_RANGE
#  undef JV_RANGE
!
! Exchange periodic boundaries, if so prescribed.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             Huon(START_2D_ARRAY,1))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             Hvom(START_2D_ARRAY,1))
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u(START_2D_ARRAY,1,nrhs))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v(START_2D_ARRAY,1,nrhs))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        Huon(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hvom(START_2D_ARRAY,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,nrhs))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,nrhs))
#   endif
#  endif
# endif /* SOLVE3D */
      return
      end
#else
      subroutine set_depth_empty
      return
      end
#endif /* SOLVE3D || WET_DRY */

