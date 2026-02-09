#include "cppdefs.h"

      subroutine das_zeta_avg_tide_bias (tile)
!
      implicit none
      integer i,j
      real cft, crt,cbias, gall
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
      crt=0.6
      cbias=0.0
      gall=0.0
      do j=0,Mm+1
        do i=0,Lm+1
          cft=ssh_swot(i,j)-(zeta_sm(i,j)-ssh_ref(i,j))
          if (abs(cft) .gt. crt ) swot_mask(i,j)=0.0
          swot_mask(i,j)=swot_mask(i,j)*ref_mask(i,j)
          cbias=cbias+cft * swot_mask(i,j)
          gall=gall + swot_mask(i,j)
        enddo
      enddo
      if (gall .tt. 50.0) then
        cbias=cbias/gall
      else
        
!
!
      return
      end
