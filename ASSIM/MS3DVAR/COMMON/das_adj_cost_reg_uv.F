#include "cppdefs.h"

      subroutine das_adj_cost_reg_uv (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_cost_reg_uv_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_cost_reg_uv_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
       real nozero
!
# include "compute_auxiliary_bounds.h"

!
! u
!
      do k=1,NDAS
        do j=Jstr,Jend
          do i=Istr,IendR
            u_s_adj(i,j,k)=u_s_adj(i,j,k)
     &             + u_s(i,j,k) * reguv
          enddo
        enddo
      enddo
!
! v
!
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,Iend
            v_s_adj(i,j,k)=v_s_adj(i,j,k)
     &             + v_s(i,j,k) * reguv
          enddo
        enddo
      enddo
!
      return
      end

