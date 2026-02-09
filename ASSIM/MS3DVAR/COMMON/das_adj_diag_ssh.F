#include "cppdefs.h"

      subroutine das_adj_diag_ssh (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_adj_diag_ssh_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_diag_ssh_tile (Istr,Iend,Jstr,Jend) 
!
! compute ssh increments by integration of the
! linearized hydrostatic equation
!
      implicit none
# include "param.h"
# include "das_param.h"
# include "das_ocean.h"
# include "das_ocean9.h"
# include "scalars.h"
!
      integer Istr,Iend,Jstr,Jend, i,j,k 
!
# include "compute_extended_bounds.h"
!
! rho is the deviation from the Boussinesque approximation mean density
! total density is rho+rho0 (rho0=1000kg/m^3)
!
        do j=JstrR,JendR
          do i=IstrR,IendR
!            zeta_h_adj(i,j)=zeta_h_adj(i,j) + zeta_s_adj(i,j)

#if defined DAS_BKGZETA
            rho_s_adj(i,j,NDAS)=rho_s_adj(i,j,NDAS)
     &          - zeta_h_adj(i,j) * zetab_das9(i,j)
     &                    /(rho0+rho_das9(i,j,NDAS))
#else
            rho_s_adj(i,j,NDAS)=rho_s_adj(i,j,NDAS)
     &          - zeta_h_adj(i,j) * zeta_das9(i,j)
     &                    /(rho0+rho_das9(i,j,NDAS))
#endif
            zeta_h_adj(i,j)=zeta_h_adj(i,j)
     &                    /(rho0+rho_das9(i,j,NDAS))
#if defined DAS_ZETAHADJCORR
            do k=NDAS-1,1,-1
#else
            do k=1,NDAS-1
#endif
              rho_s_adj(i,j,k)=rho_s_adj(i,j,k)
     &            - 0.5*(z_das(k+1) - z_das(k))*zeta_h_adj(i,j)
              rho_s_adj(i,j,k+1)=rho_s_adj(i,j,k+1)
     &            - 0.5*(z_das(k+1) - z_das(k))*zeta_h_adj(i,j)
            enddo
            zeta_h_adj(i,j)=0.0
          enddo
        enddo
      return
      end

