! $Id: cfb_stress.F $
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
#if defined SMFLUX_CFB && !defined CFB_WIND

      subroutine cfb_stress (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call cfb_stress_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine cfb_stress_tile (Istr,Iend,Jstr,Jend)

      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "scalars.h"
      integer i,j, Istr,Iend,Jstr,Jend,
     &        imin,imax,jmin,jmax
      real    cff,cff1, Wspd_min,Wstr_min, stau_ref

#  include "compute_extended_bounds.h"
!
!---------------------------------------------------------------------
!  Correct kinematic, surface wind stress (m2/s2) due to surface 
!  current feedback
!---------------------------------------------------------------------
!
!         sustr = sustr + s_tau * Uo,  (or equivalently for svstr)
!
! with  s_tau = cfb_slope  * Wspd + cfb_offset  [N.m^-3.s]  CFB_STRESS
!  or   s_tau = cfb_slope2 * Wstr + cfb_offset2 [N.m^-3.s]  CFB_STRESS2
!
!      For wind wpeed  Wspd < 3 m/s,          CFB_STRESS
!          wind stress Wstr < 0.045 N/m^2     CFB_STRESS2
!            --> take s_tau mean value of -0.0027 N.m^-3.s
!
!      Caution: sustr is in [m^2/s^2], sms_scale=1/rho0 [m^3/kg]
!      must be used for conversion
!
!---------------------------------------------------------------------
!
      Wspd_min=3.        ! [m/s]
      Wstr_min=0.045     ! [N/m2]
      stau_ref=-0.0027   ! [N.m^-3.s]

      do j=JstrR,JendR
        do i=Istr,IendR
# ifdef CFB_STRESS
          cff1=0.5*(wspd(i,j)+wspd(i-1,j))
          if (cff1 .gt. Wspd_min) then
            cff = cfb_slope*cff1 + cfb_offset
          else
            cff = stau_ref
          endif
# else 
          cff1=wstr_u(i,j)*rho0    ! CFB_STRESS2 case
          if (cff1 .gt. Wstr_min) then 
            cff = cfb_slope2*cff1 + cfb_offset2
          else
            cff = stau_ref
          endif
# endif
          sustr(i,j)=sustr(i,j) + cff*u(i,j,N,nrhs)/rho0
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
# ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
# endif
        enddo
      enddo
!
      do j=Jstr,JendR
        do i=IstrR,IendR
# ifdef CFB_STRESS
          cff1=0.5*(wspd(i,j)+wspd(i,j-1))
          if (cff1 .gt. Wspd_min) then
            cff = cfb_slope*cff1 + cfb_offset
          else
            cff = stau_ref
          endif
# else 
          cff1=wstr_v(i,j)*rho0    ! CFB_STRESS2 case
          if (cff1 .gt. Wstr_min) then
            cff = cfb_slope2*cff1 + cfb_offset2
          else
            cff = stau_ref
          endif
# endif
          svstr(i,j)=svstr(i,j) + cff*v(i,j,N,nrhs)/rho0
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# endif
        enddo
      enddo
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,sustr(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,svstr(START_2D_ARRAY))
# endif
!
      return
      end
!
!==============================================================================
!
# ifdef CFB_STRESS2
      subroutine cfb_stress_cor (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call cfb_stress_cor_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine cfb_stress_cor_tile (Istr,Iend,Jstr,Jend)

      implicit none
#  include "param.h"
#  include "scalars.h"
#  include "forces.h"
      integer i,j, Istr,Iend,Jstr,Jend,
     &        imin,imax,jmin,jmax
      real    cff

#  include "compute_extended_bounds.h"
#  ifdef MPI
#   define LOCALLM Lmmpi
#   define LOCALMM Mmmpi
#  else
#   define LOCALLM Lm
#   define LOCALMM Mm
#  endif

      if ( Iend .eq. LOCALLM ) then
#  ifdef EW_PERIODIC
        imax=Iend+1
#  else
#   ifdef MPI
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=Iend
        endif
#   else
        imax=Iend
#   endif
#  endif
      else
        imax=Iend
      endif

      if ( Jend .eq. LOCALMM ) then
#  ifdef NS_PERIODIC
        jmax=Jend+1
#  else
#   ifdef MPI
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=Jend
        endif
#   else
        jmax=Jend
#   endif
#  endif
      else
        jmax=Jend
      endif
!
!---------------------------------------------------------------------
!
!  Compute wind stress magnitude at u and v points
!
!---------------------------------------------------------------------
!
      do j=JstrR,jmax
        do i=Istr,IendR
          cff=0.25*(svstr(i-1,j  )+svstr(i,j  )
     &             +svstr(i-1,j+1)+svstr(i,j+1))
          wstr_u(i,j)=sqrt(sustr(i,j)**2+cff**2)
        enddo
      enddo
!
      do j=Jstr,JendR
        do i=IstrR,imax
          cff=0.25*(sustr(i  ,j-1)+sustr(i  ,j)
     &             +sustr(i+1,j-1)+sustr(i+1,j))
          wstr_v(i,j)=sqrt(svstr(i,j)**2+cff**2)
        enddo
      enddo
!
!  Set boundary conditions
!
      if (NORTHERN_EDGE) then
        do i=Istr,IendR
          wstr_u(i,Jend+1)=wstr_u(i,Jend)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr,JendR
          wstr_v(Iend+1,j)=wstr_v(Iend,j)
        enddo
      endif
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_u2d_tile(Istr,Iend,Jstr,Jend,wstr_u(START_2D_ARRAY))
!      call exchange_v2d_tile(Istr,Iend,Jstr,Jend,wstr_v(START_2D_ARRAY))
!#  endif

      return
      end
# endif /* CFB_STRESS2 */

#else
      subroutine cfb_stress_empty
      return
      end
#endif
