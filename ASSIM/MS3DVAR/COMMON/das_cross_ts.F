#include "cppdefs.h"

      subroutine das_cross_ts (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_cross_ts_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_cross_ts_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real cts
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_covar.h"
#include "das_ocean.h"
#include "das_ocean9.h"
!
# include "compute_extended_bounds.h"
!
! 2D
!
!
!  Initialize 3-D primitive variables.
!
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_w(i,j,k,itemp)=t_das(i,j,k,itemp) 
     &              + cross_tsp * corr_ts(i,j,k)
     &                         * t_das(i,j,k,isalt)
              t_w(i,j,k,isalt)=t_das(i,j,k,isalt)
     &              + cross_tsp * corr_ts(i,j,k)
     &                         * t_das(i,j,k,itemp)
            enddo
          enddo
        enddo
      return
      end
