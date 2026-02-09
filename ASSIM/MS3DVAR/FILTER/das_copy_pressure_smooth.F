#include "cppdefs.h"

      subroutine das_copy_pressure_sm (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_copy_pressure_sm_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_copy_pressure_sm_tile
     &                  (Istr,Iend,Jstr,Jend)
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
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              pp_s(i,j,k)=p_sm(i,j,k)
            enddo
          enddo
        enddo
      return
      end
