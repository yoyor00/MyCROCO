#include "cppdefs.h"

      subroutine das_adj_get_zeta_ubal_from_tot (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_get_zeta_ubal_from_tot_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_get_zeta_ubal_from_tot_tile
     &                     (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_ocean.h"
!
# include "compute_extended_bounds.h"
!
!  Initialize 3-D primitive variables.
!
        do j=JstrR,JendR
          do i=IstrR,IendR
            zeta_s_adj(i,j)= zeta_s_adj(i,j)
            zeta_h_adj(i,j)=-zeta_s_adj(i,j) + zeta_h_adj(i,j)
          enddo
        enddo
      return
      end
