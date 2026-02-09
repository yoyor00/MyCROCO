#include "cppdefs.h"
#if defined DAS_BKGZETA
      subroutine das_diag_steric_sshb (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_diag_steric_sshb_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_diag_steric_sshb_tile (Istr,Iend,Jstr,Jend) 
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
            zetab_das9(i,j)=0.0
            do k=1,NDAS-1                              
! Oops quite large steric ssh ! 
!              zetab_das9(i,j)=zetab_das9(i,j)
!     &         - (rho0+0.5*(rho_das9(i,j,k) + rho_das9(i,j,k+1)))
!     &         * (z_das(k+1) - z_das(k))
              zetab_das9(i,j)=zetab_das9(i,j)
     &         - (0.5*(rho_das9(i,j,k) + rho_das9(i,j,k+1)))
     &         * (z_das(k+1) - z_das(k))
            enddo
            zetab_das9(i,j)=zetab_das9(i,j)/(rho0+rho_das9(i,j,NDAS))
          enddo
        enddo
      return
      end
#endif /* DAS_BKGZETA */
