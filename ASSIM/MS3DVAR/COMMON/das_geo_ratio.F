#include "cppdefs.h"

      subroutine das_geo_ratio (tile)
      implicit none
      integer tile
# include "param.h"
# include "das_param.h"
# include "das_private_scratch.h"
# include "compute_tile_bounds.h"
      call das_geo_ratio_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_geo_ratio_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real geor
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
!
!  u
!

      do k=1,NDAS

        do j=JstrR,JendR
          do i=Istr,IendR
            geor=0.5*(georatio(i-1,j)+georatio(i,j))
            u_s(i,j,k) = geor*u_s(i,j,k)
          enddo
        enddo
!
!  v
!
        do j=Jstr,JendR
          do i=IstrR,IendR
            geor=0.5*(georatio(i,j-1)+georatio(i,j))
            v_s(i,j,k) = geor*v_s(i,j,k)
          enddo
        enddo
!      
      enddo  !k
!          
      return
      end
