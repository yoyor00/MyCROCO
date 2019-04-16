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
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "forces.h"
#include "scalars.h"
      integer i,j, Istr,Iend,Jstr,Jend,
     &        imin,imax,jmin,jmax
      real    cff,cff1,cff2

# include "compute_extended_bounds.h"

# ifdef CFB_STRESS
      imax=IendR
      jmax=JendR
# else
      if (EASTERN_EDGE) then
        imax=Iend
      else
        imax=IendR
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend
      else
        jmax=JendR
      endif
# endif
!
!---------------------------------------------------------------------
!  Correct kinematic, surface wind stress (m2/s2) due to surface 
!  current feedback
!---------------------------------------------------------------------
!
! sustr = sustr + s_tau * Uo,  (or equivalently for svstr)
!
!   with  s_tau = cfb_slope1 * Uatm + cfb_offset1 [N.m^-3.s]  CFB_STRESS
!    or   s_tau = cfb_slope2 * Wstr + cfb_offset2 [N.m^-3.s]  CFB_STRESS2
!
!      For wind Uatm < 3 m/s, or wind stress Wstr < 0.045 N/m^2: 
!            --> take s_tau mean value of -0.0027 N.m^-3.s
!      Caution: sustr is in [m^2/s^2], sms_scale=1/rho0 [m^3/kg] must be
!      used for conversion
!
!---------------------------------------------------------------------
!
      do j=JstrR,jmax
        do i=Istr,IendR
# ifdef CFB_STRESS
          cff1=0.5*(wspd(i,j)+wspd(i-1,j))
          if (cff1 .gt. 3) then
              cff = cff1*cfb_slope + cfb_offset
          else
              cff = -0.0027
          endif
# else 
          cff2=0.25*(svstr(i-1,j  )+svstr(i,j  )  ! CFB_STRESS2 case
     &              +svstr(i-1,j+1)+svstr(i,j+1))
          cff1=sqrt(sustr(i,j)**2+cff2**2)/sms_scale
          if (cff1 .gt. 0.045) then 
            cff = cfb_slope2*cff1 + cfb_offset2
          else
            cff = -0.0027
          endif
# endif
          sustr(i,j)=sustr(i,j)+sms_scale*cff*u(i,j,N,nrhs)
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
        do i=IstrR,imax
# ifdef CFB_STRESS
          cff1=0.5*(wspd(i,j)+wspd(i,j-1))
          if (cff1 .gt. 3) then
            cff = cfb_slope*cff1 + cfb_offset
          else
            cff = -0.0027
          endif
# else 
          cff2=0.25*(sustr(i  ,j-1)+sustr(i  ,j)  ! CFB_STRESS2 case
     &              +sustr(i+1,j-1)+sustr(i+1,j))
          cff1=sqrt(svstr(i,j)**2+cff2**2)/sms_scale
          if (cff1 .gt. 0.045) then
            cff = cfb_slope2*cff1 + cfb_offset2
          else
            cff = -0.0027
          endif
# endif
          svstr(i,j)=svstr(i,j)+sms_scale*cff*v(i,j,N,nrhs)
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# endif
        enddo
      enddo
!
# ifndef CFB_STRESS
!  Set boundary conditions for CFB_STRESS2 case
      if (NORTHERN_EDGE) then
        do i=Istr,IendR
          sustr(i,Jend+1)=sustr(i,Jend)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr,JendR
          svstr(Iend+1,j)=svstr(Iend,j)
        enddo
      endif
# endif
!
      return
      end
#else
      subroutine cfb_stress_empty
      return
      end
#endif
