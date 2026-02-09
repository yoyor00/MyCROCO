#include "cppdefs.h"

      subroutine das_copy_zetah_sm (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_copy_zetah_sm_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_copy_zetah_sm_tile
     &                  (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_ocean.h"
#include "das_ocean_smooth.h"
!
# include "compute_extended_bounds.h"
!
!  Initialize 3-D primitive variables.
!
        do j=JstrR,JendR
          do i=IstrR,IendR
            zeta_h(i,j)=zeta_sm(i,j)
            zeta_s(i,j)=zeta_s(i,j)+zeta_sm(i,j)
          enddo
        enddo
      return
      end
