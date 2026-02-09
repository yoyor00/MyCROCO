#include "cppdefs.h"
#if defined DAS_MCSST || defined DAS_TMISST

      subroutine das_adj_csst_lin (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_csst_lin_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_csst_lin_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, cfm
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
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
# ifdef DAS_TMISST
           cft=t_s(i,j,NDAS,itemp)-sst_tmi(i,j)
           t_s_adj(i,j,NDAS,itemp)=t_s_adj(i,j,NDAS,itemp)
     &           +cft*tmi_mask(i,j)*tmi_oin(i,j)*costf_adj
# endif
# ifdef DAS_MCSST
           cfm=t_s(i,j,NDAS,itemp)-sst_mc(i,j)
           t_s_adj(i,j,NDAS,itemp)=t_s_adj(i,j,NDAS,itemp)
     &           +cfm* mc_mask(i,j)* mc_oin(i,j)*costf_adj
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
      subroutine das_adj_csst_lin_empty
      return
      end
#endif /* DAS_TMISST DAS_MCSST */

