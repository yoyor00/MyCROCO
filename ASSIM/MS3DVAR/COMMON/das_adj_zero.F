#include "cppdefs.h"

      subroutine das_adj_zero (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_zero_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_zero_tile (Istr,Iend,Jstr,Jend)
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
          zeta_adj(i,j)=init
          zeta_s_adj(i,j)=init
          zeta_h_adj(i,j)=init
          zeta_sm_adj(i,j)=init
#if defined DAS_PWTHZETATOT
          zeta_h_sm_adj(i,j)=init
#endif
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
            psi_adj(i,j,k)=init
            psi_s_adj(i,j,k)=init
            chi_adj(i,j,k)=init
            chi_s_adj(i,j,k)=init
          enddo
        enddo
      enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_adj(i,j,k,itemp)=init
              t_s_adj(i,j,k,itemp)=init
              t_sm_adj(i,j,k,itemp)=init
              t_w_adj(i,j,k,itemp)=init
            enddo
          enddo
        enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_adj(i,j,k,isalt)=init
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
              p_sm_adj(i,j,k)=init
            enddo
          enddo
        enddo
      return
      end
