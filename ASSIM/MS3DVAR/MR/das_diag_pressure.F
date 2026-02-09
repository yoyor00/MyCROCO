#include "cppdefs.h"

      subroutine das_diag_pressure (tile)
      implicit none
      integer tile,trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_diag_pressure_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_diag_pressure_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
! compute geostropic current from sea surface height increments
! zeta_s and density increments rho_s. called by ROMS-DAS
!
! It is not applicable near the equator(-1.5 - 1.5). 
!
! In consistence with ROMS that computes XI-component (ETA-component)
! of pressure at u-grids (v-grids), v-component (u-component) of
! geostrophic flow is thus computed at u-grids (v-grids) and then interpolated
! to v-grids (u-grids). 
!
! In the geometry vertical coordinates, there is a possibilty of
! a single column or row. The increment in those areas is zero in
! the following pragram.
!
! Boundary condition: sliding, and velocity normal to coastline zero
! 
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real cff, cff1, crit_lat, crit, rho02
      real a0
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean_smooth.h"
# include "das_ocean9.h"
!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
! ... p=P/g
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            pp_s(i,j,k)=0.
          enddo
        enddo
      enddo

      do j=JstrR,JendR
        do i=IstrR,IendR
#if defined DAS_PWTHZETATOT
! ... Add contribution of the smoothed unbalanced SSH
# if defined DAS_BKGZETA
          pp_s(i,j,NDAS) = rho_s(i,j,NDAS)*zetab_das9(i,j) 
# else
          pp_s(i,j,NDAS) = rho_s(i,j,NDAS)*zeta_das9(i,j) 
# endif
     &           + (rho0+rho_das9(i,j,NDAS))*( 
     &             (1.-alph0)  *  zeta_h(i,j)  + 
     &              beta0      * zeta_sm(i,j)  +
     &             (gama0-beta0)* zeta_h_sm(i,j))
#else
# if defined DAS_BKGZETA
          pp_s(i,j,NDAS) = rho_s(i,j,NDAS)*zetab_das9(i,j) 
# else
          pp_s(i,j,NDAS) = rho_s(i,j,NDAS)*zeta_das9(i,j) 
# endif
     &           + (rho0+rho_das9(i,j,NDAS))*zeta_h(i,j)
#endif
        enddo
      enddo

      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=NDAS-1, 1, -1
            pp_s(i,j,k)=pp_s(i,j,k+1)
     &            + 0.5*(rho_s(i,j,k) + rho_s(i,j,k+1))
     &             * (z_das(k+1) - z_das(k))
          enddo
        enddo
      enddo
!          
      return
      end
