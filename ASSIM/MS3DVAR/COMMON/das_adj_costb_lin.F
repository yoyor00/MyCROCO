#include "cppdefs.h"

      subroutine das_adj_costb_lin (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_costb_lin_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_costb_lin_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
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
      do k=1,NDAS
        do j=JR_RANGE
          do i=IR_RANGE
            t_adj(i,j,k,itemp)=t_das(i,j,k,itemp)
     &              * rmask_das(i,j,k)*costf_adj
            t_adj(i,j,k,isalt)=t_das(i,j,k,isalt)
     &              * rmask_das(i,j,k)*costf_adj
          enddo
        enddo
      enddo
!
# if !defined DAS_HYDRO_STRONG & !defined DAS_HYDRO_WEAK
      do j=JR_RANGE
        do i=IR_RANGE
          zeta_adj(i,j)=zeta_das(i,j)
     &              * rmask_das(i,j,k)*costf_adj
        enddo
      enddo
# endif
#if !defined DAS_GEOS_STRONG & !defined DAS_GEOS_WEAK
!
! u
!
      do k=1,NDAS
        do j=JR_RANGE
          do i=IU_RANGE
            u_adj(i,j,k)=u_das(i,j,k)
     &                * umask_das(i,j,k)*costf_adj
          enddo
        enddo
      enddo
!
! v
!
      do k=1,NDAS
        do j=JV_RANGE
          do i=IR_RANGE
            v_adj(i,j,k)=v_das(i,j,k)
     &                * vmask_das(i,j,k)*costf_adj
          enddo
        enddo
      enddo
#endif
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
!
      return
      end

