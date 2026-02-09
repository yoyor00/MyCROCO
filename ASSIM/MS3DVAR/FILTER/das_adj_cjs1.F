#include "cppdefs.h"
#if defined DAS_JASONSSH || defined DAS_TPSSH \
  || defined DAS_ERS2SSH || defined DAS_GFOSSH

      subroutine das_adj_cjs1 (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_cjs1_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_cjs1_tile (Istr,Iend,Jstr,Jend) 
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
# ifdef DAS_JASONSSH
           cft=zeta_sm(i,j)-ssh_js1(i,j)
           zeta_sm_adj(i,j)=zeta_sm_adj(i,j)
     &            +cft*js1_mask(i,j)*js1_oin(i,j)
# endif
# ifdef DAS_TPSSH
           cfm=zeta_sm(i,j)-ssh_tp(i,j)
           zeta_sm_adj(i,j)=zeta_sm_adj(i,j)
     &             +cfm*tp_mask(i,j)*tp_oin(i,j)
# endif
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
      subroutine das_adj_cjs1_empty
      return
      end
#endif /* DAS_JASONSSH DAS_TPSSH */

