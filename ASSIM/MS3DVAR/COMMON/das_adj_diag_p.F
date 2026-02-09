#include "cppdefs.h"

      subroutine das_adj_diag_p (tile)
      implicit none
      integer tile,trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_adj_diag_p_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_diag_p_tile (Istr,Iend,Jstr,Jend)
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

      do j=JstrR,JendR                                                                                                                   
        do i=IstrR,IendR                                                                                                                 
          do k=1,NDAS-1
            rho_s_adj(i,j,k)=rho_s_adj(i,j,k)
     &           +0.5*p_s_adj(i,j,k)* (z_das(k+1) - z_das(k))
            rho_s_adj(i,j,k+1)=rho_s_adj(i,j,k+1)
     &           +0.5*p_s_adj(i,j,k)* (z_das(k+1) - z_das(k))
            p_s_adj(i,j,k+1)=p_s_adj(i,j,k+1)+p_s_adj(i,j,k)
            p_s_adj(i,j,k)=0.
          enddo                                                                                                                          
        enddo                                                                                                                            
      enddo
!
      do j=JstrR,JendR                                                                                                                   
        do i=IstrR,IendR                                                                                                                 
          rho_s_adj(i,j,NDAS)=rho_s_adj(i,j,NDAS)
# if defined DAS_BKGZETA
     &             +p_s_adj(i,j,NDAS)*zetab_das9(i,j)
# else
     &             +p_s_adj(i,j,NDAS)*zeta_das9(i,j)
# endif
#if defined DAS_PWTHZETATOT
! ... Add adjoint contribution of the smoothed unbalanced SSH
          zeta_h_adj(i,j)=zeta_h_adj(i,j)
     &             +(rho0+rho_das9(i,j,NDAS))*(1.-alph0)
     &                                       *p_s_adj(i,j,NDAS)
          zeta_sm_adj(i,j)=zeta_sm_adj(i,j)
     &             +(rho0+rho_das9(i,j,NDAS))*beta0
     &                                       *p_s_adj(i,j,NDAS)
          zeta_h_sm_adj(i,j)=zeta_h_sm_adj(i,j)
     &             +(rho0+rho_das9(i,j,NDAS))*(gama0-beta0)
     &                                       *p_s_adj(i,j,NDAS)
#else
          zeta_h_adj(i,j)=zeta_h_adj(i,j)
     &             +(rho0+rho_das9(i,j,NDAS))*p_s_adj(i,j,NDAS)
#endif
          p_s_adj(i,j,NDAS)=0.
        enddo                                                                                                                            
      enddo

      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            p_s_adj(i,j,k)=0.
          enddo
        enddo
      enddo
!          
      return
      end
