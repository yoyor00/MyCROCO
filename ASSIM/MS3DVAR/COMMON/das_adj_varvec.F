#include "cppdefs.h"

      subroutine das_adj_varvec (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_varvec_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_varvec_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k,ie,je,ke
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
!
#ifdef MPI
      include 'mpif.h'
      integer size, step, status(MPI_STATUS_SIZE), ierr
      real buff
#endif

# include "compute_auxiliary_bounds.h"
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_s_adj(i,j,k,itemp)=t_s_adj(i,j,k,itemp)
     &                 * bt_das(i,j,k,itemp)
            t_s_adj(i,j,k,isalt)=t_s_adj(i,j,k,isalt)
     &                 * bt_das(i,j,k,isalt)
#if !defined DAS_GEOS_STRONG
            chi_s_adj(i,j,k)=chi_s_adj(i,j,k) * bchi_das(i,j,k)
#endif
          enddo
        enddo
      enddo
#if !defined DAS_GEOS_STRONG
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,IendR
            psi_s_adj(i,j,k)=psi_s_adj(i,j,k) * bpsi_das(i,j,k)
          enddo
        enddo
      enddo
#endif
#if !defined DAS_HYDRO_STRONG
        do j=JstrR,JendR
          do i=IstrR,IendR
            zeta_s_adj(i,j)=zeta_s_adj(i,j) * bz_das(i,j)
          enddo
        enddo
#endif
!
      return
      end

