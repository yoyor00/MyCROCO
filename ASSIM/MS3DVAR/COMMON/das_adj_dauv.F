#include "cppdefs.h"

      subroutine das_adj_dauv (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_dauv_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_dauv_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, cfs, kdas
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
!      dor_ot=1.0/(0.4*0.4)
!      dor_os=1.0/(0.25*0.25)
!
      do k=1,max_prf_dor
        kdas=NDAS-k+1
        do j=JR_RANGE
          do i=IR_RANGE
            cft=t_s(i,j,kdas,itemp)-t_dor(i,j,k)
            cfs=t_s(i,j,kdas,isalt)-s_dor(i,j,k)
            t_s_adj(i,j,kdas,itemp)=t_s_adj(i,j,kdas,itemp)
     &                       +cft*dor_mask(i,j,k)*dor_ot
            t_s_adj(i,j,kdas,isalt)=t_s_adj(i,j,kdas,isalt)
     &                       +cfs*dor_mask(i,j,k)*dor_os
          enddo
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
