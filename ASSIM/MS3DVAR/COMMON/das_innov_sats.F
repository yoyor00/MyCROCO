#include "cppdefs.h"
#if defined DAS_SATSSSS

      subroutine das_innov_sats (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_sats_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_sats_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real cft, crt
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
!
!
# include "compute_auxiliary_bounds.h"
!
# ifdef DAS_SATSSSS
      crt=0.5
      do j=JstrR,JendR
        do i=IstrR,IendR
          cft=sss_sats(i,j)-t_das(i,j,NDAS,isalt)
          if (abs(cft) .gt. crt ) sats_mask(i,j)=0.0
          sss_sats(i,j)=cft * sats_mask(i,j)
          sats_oin(i,j)=sats_osss
        enddo
      enddo
# endif
      return
      end
#else
      subroutine das_innov_sats_empty
      return
      end
#endif /* DAS_SATSSSS  */
