#include "cppdefs.h"
#if defined DAS_SWOTSSH 

      subroutine das_innov_swot (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_swot_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_swot_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
      real cft, crt
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
# include "das_ocean_smooth.h"
!
!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
      crt=0.60
      do j=JstrR,JendR
        do i=IstrR,IendR
          cft=ssh_swot(i,j)-(zeta_sm(i,j)-ssh_ref(i,j))
          if (abs(cft) .gt. crt ) swot_mask(i,j)=0.0
          swot_mask(i,j)=swot_mask(i,j)*ref_mask(i,j)
          ssh_swot(i,j)=cft * swot_mask(i,j)
          swot_oin(i,j)=swot_ossh
        enddo
      enddo
!
!
      return
      end
#else
      subroutine das_innov_swot_empty
      return
      end
#endif /* DAS_SWOTSSH */
