#include "cppdefs.h"

      subroutine das_adj_map (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_map_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_map_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, cfs, kdas
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_ocean_smooth.h"
# include "das_covar.h"
# include "das_lbgfs.h"
# include "das_innov.h"
!
# include "compute_auxiliary_bounds.h"

!
      do k=1,max_prf
        kdas=NDAS-k+1
        do j=JstrR,JendR
          do i=IstrR,IendR
            cft=t_sm(i,j,kdas,itemp)-t_map(i,j,k)
            cfs=t_sm(i,j,kdas,isalt)-s_map(i,j,k)
            t_sm_adj(i,j,kdas,itemp)=t_sm_adj(i,j,kdas,itemp)
     &                       +cft*map_t_mask(i,j,k)*map_ot(k)
            t_sm_adj(i,j,kdas,isalt)=t_sm_adj(i,j,kdas,isalt)
     &                       +cfs*map_s_mask(i,j,k)*map_os(k)
          enddo
        enddo
      enddo
!
!
      return
      end
