#include "cppdefs.h"

      subroutine das_adj_costb (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_costb_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_costb_tile (Istr,Iend,Jstr,Jend) 
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

      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_adj(i,j,k,itemp)=t_w(i,j,k,itemp) 
            t_adj(i,j,k,isalt)=t_w(i,j,k,isalt)
          enddo
        enddo
      enddo
!
# if !defined DAS_HYDRO_STRONG
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_adj(i,j)=zeta_das(i,j)
     &                         
        enddo
      enddo
# endif
#if !defined DAS_GEOS_STRONG
!
! psi
!
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,IendR
            psi_adj(i,j,k)=psi_das(i,j,k)
          enddo
        enddo
      enddo
!
! chi
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            chi_adj(i,j,k)=chi_das(i,j,k)
          enddo
        enddo
      enddo
# endif
!
      return
      end

