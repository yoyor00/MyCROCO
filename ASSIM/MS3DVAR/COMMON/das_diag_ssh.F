#include "cppdefs.h"

      subroutine das_diag_ssh (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_diag_ssh_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_diag_ssh_tile (Istr,Iend,Jstr,Jend) 
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
      integer Istr,Iend,Jstr,Jend, i,j,k 

# include "compute_extended_bounds.h"
!
! rho is the deviation from the Boussinesque approximation mean density
! total density is rho+rho0 (rho0=1000kg/m^3)
!
        do j=JstrR,JendR
          do i=IstrR,IendR       
            zeta_h(i,j)=0.0
            do k=1,NDAS-1                              
              zeta_h(i,j)=zeta_h(i,j)
     &            -  0.5*(rho_s(i,j,k) + rho_s(i,j,k+1))
     &             * (z_das(k+1) - z_das(k))
            enddo
#if defined DAS_BKGZETA
            zeta_h(i,j)=(zeta_h(i,j) - rho_s(i,j,NDAS)*zetab_das9(i,j))
     &                    /(rho0+rho_das9(i,j,NDAS))
#else
            zeta_h(i,j)=(zeta_h(i,j) - rho_s(i,j,NDAS)*zeta_das9(i,j))
     &                    /(rho0+rho_das9(i,j,NDAS))
#endif
!            zeta_s(i,j)=zeta_s(i,j) + zeta_h(i,j)
          enddo
        enddo
      return
      end

