#include "cppdefs.h"

      subroutine das_adj_copy_pressure_sm (tile)
      implicit none
      integer tile
#include "param.h"
#include "ms3dvar_analysis_param.h"
#include "compute_tile_bounds.h"
      call das_adj_copy_pressure_sm_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_copy_pressure_sm_tile
     &                     (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#include "param.h"
#include "ms3dvar_analysis_param.h"
#include "scalars.h"
#include "ms3dvar_ocean.h"
#include "ms3dvar_ocean_smooth.h"
!
# include "compute_extended_bounds.h"
!
!  Initialize 3-D primitive variables.
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            p_sm_adj(i,j,k)=p_sm_adj(i,j,k) + p_s_adj(i,j,k)
            p_s_adj(i,j,k)=0.0
          enddo
        enddo
      enddo
      return
      end
