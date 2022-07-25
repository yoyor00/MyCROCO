! $Id: mld_rho.F 1458 2014-02-03 15:01:25Z gcambon $
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
!
#if defined MLD_RHO
      subroutine mld_rho (tile)
!
      implicit none
      integer tile, trd, omp_get_thread_num
#  include "param.h"
#  include "private_scratch.h"
#  include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call mld_rho_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine mld_rho_tile (Istr,Iend,Jstr,Jend)
!
!======================================================================
! Compute Mixed-layer depth using ...
!
! References:
! ----------
!
!
!====================================================================== 
!
      implicit none
#  include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j,k, kbl_tmp
#  include "grid.h"
#  include "ocean3d.h"
#  include "coupling.h"
#  include "scalars.h"
!
#  include "compute_extended_bounds.h"
!


      do j=JstrR,JendR
        do i=IstrR,IendR

          kbl_rho(i,j) = N-1
            
            do k=N-1,1,-1
              if ((rho1(i,j,k)-rho1(i,j,N)) .lt. 0.03) then
                kbl_rho(i,j) = k
              endif
            enddo   ! <-- k
            
            hbl_rho(i,j) = ( z_w(i,j,N) - z_r(i,j,kbl_rho(i,j)) )
#  ifdef MASKING
     &                                                *rmask(i,j)
#  endif     
        enddo   ! <-- I
      enddo    ! <-- j
      
      return
      end
#else
      subroutine mld_rho_empty
      end
#endif /* MLD_RHO*/
