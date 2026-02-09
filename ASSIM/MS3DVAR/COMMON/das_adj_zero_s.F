#include "cppdefs.h"

      subroutine das_adj_zero_s (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_zero_s_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_zero_s_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_ocean.h"
#include "das_ocean_smooth.h"

      real init    !!!!  0xFFFA5A5A ==> NaN

      parameter (init=0.)

!
# include "compute_extended_bounds.h"
!
!
!  2-D variables.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_s_adj(i,j)=init
        enddo
      enddo
!
!  Initialize 3-D primitive variables.
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            u_s_adj(i,j,k)=init
            v_s_adj(i,j,k)=init
          enddo
        enddo
      enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_s_adj(i,j,k,itemp)=init
              t_sm_adj(i,j,k,itemp)=init
              t_w_adj(i,j,k,itemp)=init
            enddo
          enddo
        enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_s_adj(i,j,k,isalt)=init
              t_sm_adj(i,j,k,isalt)=init
              t_w_adj(i,j,k,isalt)=init
            enddo
          enddo
        enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              rho_s_adj(i,j,k)=init
              p_s_adj(i,j,k)=init
            enddo
          enddo
        enddo
      return
      end
