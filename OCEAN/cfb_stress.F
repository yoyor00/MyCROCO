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
      real cff,cff2,cff3

# ifdef EW_PERIODIC
        imin=Istr-2
        imax=Iend+2
# else
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
# endif
# ifdef NS_PERIODIC
        jmin=Jstr-2
        jmax=Jend+2
# else
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
# endif
!
!----------------------------------------------------------------------
!  Input variables:
!
!      - wspd
!      - sustr
!      - svstr
!
!----------------------------------------------------------------------
!
!
!--------------------------------------------------------------
!  Compute kinematic, surface wind stress (m2/s2).
!--------------------------------------------------------------
!
!     sustr = sustr+stau*Uo where stau=cfb_slope*Ua+cfb_offset for wind
!     larger than 3m/s, otherwise take mean value of -0.0027
!     same for svstr
!
      do j=jmin,jmax
        do i=imin+1,imax
# ifdef CFB_STRESS
          if (wspd(i+1,j) .gt. 3) then
              cff=0.5*cfb_slope*(wspd(i,j)+wspd(i+1,j))+cfb_offset
          else
              cff = -0.0027
          endif
# else 
! CFB_STRESS2 case
          cff3=0.25*(svstr(i-1,j  )+svstr(i,j  )+
     &               svstr(i-1,j+1)+svstr(i,j+1))
          if (cff3 .gt. 0.045) then 
              cff2=sqrt(sustr(i,j)**2+cff3**2)
              cff=cfb_slope2*cff2+cfb_offset2
          else
              cff = -0.0027
          endif
# endif
          sustr(i,j)=sustr(i,j)+cff*u(i,j,N,nrhs)
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
# ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
# endif
        enddo
      enddo
!
      do j=jmin+1,jmax
        do i=imin,imax
# ifdef CFB_STRESS
          if (wspd(i+1,j) .gt. 3) then
             cff=0.5*cfb_slope*(wspd(i,j)+wspd(i,j+1))+cfb_offset
          else
              cff = -0.0027
          endif
# else 
! CFB_STRESS2 case
          cff3=0.25*(sustr(i  ,j-1)+sustr(i  ,j)+
     &               sustr(i+1,j-1)+sustr(i+1,j))
          if (cff3 .gt. 0.045) then
              cff2=sqrt(svstr(i,j)**2+cff3**2)
              cff=cfb_slope2*cff2+cfb_offset2
          else
              cff = -0.0027
          endif
# endif
          svstr(i,j)=svstr(i,j)+cff*v(i,j,N,nrhs)
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# endif
        enddo
      enddo

      return
      end
#else
      subroutine cfb_stress_empty
      return
      end
#endif
