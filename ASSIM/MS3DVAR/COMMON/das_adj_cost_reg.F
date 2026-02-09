#include "cppdefs.h"

      subroutine das_adj_cost_reg (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_cost_reg_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_cost_reg_tile (Istr,Iend,Jstr,Jend) 
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

      nozero=0.0001
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_s_adj(i,j,k,itemp)=t_s_adj(i,j,k,itemp)
     &        +t_s(i,j,k,itemp) 
     &        * regp/(bt_das(i,j,k,itemp)*bt_das(i,j,k,itemp) + nozero)
            t_s_adj(i,j,k,isalt)=t_s_adj(i,j,k,isalt)
     &        +t_s(i,j,k,isalt) 
     &        * regp/(bt_das(i,j,k,isalt)*bt_das(i,j,k,isalt) + nozero)
          enddo
        enddo
      enddo
!
# if !defined DAS_HYDRO_STRONG
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_s_adj(i,j)=zeta_s_adj(i,j)
     &          +zeta_s(i,j) 
     &          * regp/(bz_das(i,j)*bz_das(i,j) + nozero)
        enddo
      enddo
# endif

#if !defined DAS_GEOS_STRONG
# if defined DAS_PSICHI_COSTREG
!
! psi
!
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,IendR
            psi_s_adj(i,j,k)=psi_s_adj(i,j,k)
     &             + psi_s(i,j,k) 
     &             * regp/(bpsi_das(i,j,k)*bpsi_das(i,j,k)+ nozero)
          enddo
        enddo
      enddo
!
! chi
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            chi_s_adj(i,j,k)=chi_s_adj(i,j,k)
     &             + chi_s(i,j,k) 
     &             * regp/(bchi_das(i,j,k)*bchi_das(i,j,k)+ nozero)
          enddo
        enddo
      enddo
# endif
#endif
!
      return
      end

