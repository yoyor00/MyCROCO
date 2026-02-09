#include "cppdefs.h"

      subroutine das_varvec (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_varvec_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_varvec_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k, ie,je,ke
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
!
#ifdef MPI
      write (*,'(6x,A)')
     &     'The routine das_matvec is not ready for MPI'
      stop
#endif

# include "compute_auxiliary_bounds.h"
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_s(i,j,k,itemp)=t_s(i,j,k,itemp) 
     &                                     * bt_das(i,j,k,itemp)
            t_s(i,j,k,isalt)=t_s(i,j,k,isalt) 
     &                                     * bt_das(i,j,k,isalt)
#if !defined DAS_GEOS_STRONG
            chi_s(i,j,k)=chi_s(i,j,k) * bchi_das(i,j,k)
#endif
          enddo
        enddo
      enddo
!
!psi
!
#if !defined DAS_GEOS_STRONG
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,IendR
            psi_s(i,j,k)=psi_s(i,j,k) * bpsi_das(i,j,k)
          enddo
        enddo
      enddo
#endif
!
! zeta
!
#if !defined DAS_HYDRO_STRONG
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_s(i,j)=zeta_s(i,j) * bz_das(i,j)
        enddo
      enddo
#endif
!
      return
      end

