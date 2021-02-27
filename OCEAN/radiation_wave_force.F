! $Id$
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
!
#include "cppdefs.h"

# if defined RADIATION_WAVE_FORCE

      subroutine radiation_wave_force (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
#  include "param.h"
#  include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call radiation_wave_force_tile (istr,iend,jstr,jend)
      return
      end

      subroutine radiation_wave_force_tile (istr,iend,jstr,jend)
      implicit none
#  include "param.h"
      integer istr,iend,jstr,jend, i,j,k


#if defined RADIATION_WAVE_FORCE_OFFLINE
      call radiation_offline(istr,iend,jstr,jend)
#else
      call radiation2D(istr,iend,jstr,jend)
#endif

#ifdef SOLVE3D
      call vertical_radiation(istr,iend,jstr,jend)
#endif
      end subroutine radiation_wave_force_tile

#ifdef SOLVE3D
      subroutine vertical_radiation(istr,iend,jstr,jend)
      ! vertical radiation stress
      ! authors: Cristele Chevalier, Camille Mazoyer, IRD
      !          Yann Cote, Simon Huret, Cyprien Sow (2015)
      ! date: 2020
      implicit none
      integer istr,iend,jstr,jend, i,j,k
#  include "param.h"
#  include "grid.h"
#  include "ocean3d.h"
#  include "forces.h"

      do k=1,N
         do j=Jstr,Jend
            do i=Istr,Iend
               FXwave_z(i,j,k)=FXwave(i,j)*0.5*(Hz(i,j,k)+ 
     &                         Hz(i-1,j,k))/h(i,j)
               FYwave_z(i,j,k)=FYwave(i,j)*0.5*(Hz(i,j,k)+
     &                         Hz(i,j-1,k))/h(i,j)
            end do
         end do
      end do
      return
      end subroutine vertical_radiation
#endif
  

#if !defined RADIATION_WAVE_FORCE_OFFLINE
      subroutine radiation2D(istr,iend,jstr,jend)
      ! temporal radiation stress (x,y)
      ! authors: Cristele Chevalier and Camille Mazoyer,IRD
      !          Cyprien Sow, Yann Cote
      ! date: july 2020

      implicit none
#  include "param.h"
#  include "forces.h"
      integer istr,iend,jstr,jend, i,j
      
         do j=Jstr,Jend
            do i=Istr,Iend
              FXwave(i,j)=FXwave3D(i,j)*Awave(i,j)
              FYwave(i,j)=FYwave3D(i,j)*Awave(i,j)
            end do
         end do
      return
      end subroutine radiation2D
#endif
  
      
#if defined RADIATION_WAVE_FORCE_OFFLINE
      subroutine radiation_offline(istr,iend,jstr,jend)
      ! temporal radiation stress 3D (x,y,t) (offline calculus)
      ! authors: Cyprien Sow, Cristele Chevalier and Camille Mazoyer,IRD
      ! date: july 2020

      implicit none
#  include "param.h"
#  include "forces.h"
      integer istr,iend,jstr,jend, i,j
      
         do j=Jstr,Jend
            do i=Istr,Iend
              FXwave(i,j)=FXwave3D(i,j)*Awave(i,j)
              FYwave(i,j)=FYwave3D(i,j)*Awave(i,j)
            end do
         end do
      end subroutine radiation_offline
#endif

#else
      subroutine radiation_wave_force_empty
      end
#endif

