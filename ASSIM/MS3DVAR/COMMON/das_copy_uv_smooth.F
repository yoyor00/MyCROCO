#include "cppdefs.h"

      subroutine das_copy_uv_smooth (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_copy_uv_smooth_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_copy_uv_smooth_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_ocean.h"
#include "das_ocean_smooth.h"
!
# include "compute_extended_bounds.h"
!
!
!  Initialize 3-D primitive variables.
!
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              u_s(i,j,k)=u_sm(i,j,k)
              v_s(i,j,k)=v_sm(i,j,k)
            enddo
          enddo
        enddo

      return
      end
