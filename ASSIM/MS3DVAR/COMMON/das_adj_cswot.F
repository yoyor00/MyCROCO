#include "cppdefs.h"
#if defined DAS_SWOTSSH 

      subroutine das_adj_cswot (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_cswot_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_cswot_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, cfm
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

# ifdef EW_PERIODIC
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE  Istr,IendR
# endif

# ifdef NS_PERIODIC
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE  Jstr,JendR
# endif
!
      do j=JR_RANGE
         do i=IR_RANGE
           cft=zeta_sm(i,j)-ssh_swot(i,j)
           zeta_sm_adj(i,j)=zeta_sm_adj(i,j)
     &            +cft*swot_mask(i,j)*swot_oin(i,j)
         enddo
       enddo
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
!
      return
      end
#else
      subroutine das_adj_cswot_empty
      return
      end
#endif /* DAS_SWOTSSH  */

