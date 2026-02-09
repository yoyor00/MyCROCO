#include "cppdefs.h"

      subroutine das_psichi_uv (tile)
      implicit none
      integer tile,trd, mp_my_threadnum
# include "param.h"
# include "das_param.h"
# include "compute_tile_bounds.h"
      call das_psichi_uv_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_psichi_uv_tile (Istr,Iend,Jstr,Jend)
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
      real cff
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
# include "das_ocean_diag.h"

!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
! NOTE: u/v at the physical boundaries are computed,
!       but the values MUST NOT be used
!
! compute ageostrophic u/v from agetrophic psi/chi,
!
      cff=0.5
!
!  psi -> u
!
      do k=1,NDAS
!
        do j=Jstr,Jend
          do i=Istr,IendR
#if defined DAS_UVAGEO_DIAG
# if defined DAS_ADIMPSICHI
            u_ap(i,j,k) = 
     &              -cff*(psi_s(i,j+1,k) - psi_s(i,j,k))
     &               *adp*(pn(i-1,j)+pn(i,j))
            u_s(i,j,k) = u_s(i,j,k) + u_ap(i,j,k)
# else
            u_ap(i,j,k) = 
     &              -cff*(psi_s(i,j+1,k) - psi_s(i,j,k))
     &               *(pn(i-1,j)+pn(i,j))
            u_s(i,j,k) = u_s(i,j,k) + u_ap(i,j,k)
# endif
#else
# if defined DAS_ADIMPSICHI
            u_s(i,j,k) = u_s(i,j,k)
     &              -cff*(psi_s(i,j+1,k) - psi_s(i,j,k))
     &               *adp*(pn(i-1,j)+pn(i,j))
# else
            u_s(i,j,k) = u_s(i,j,k)
     &              -cff*(psi_s(i,j+1,k) - psi_s(i,j,k))
     &               *(pn(i-1,j)+pn(i,j))
# endif
#endif
          enddo
        enddo
!
! psi -> v
!
        do j=Jstr,JendR
          do i=Istr,Iend
#if defined DAS_UVAGEO_DIAG
# if defined DAS_ADIMPSICHI
            v_ap(i,j,k) =
     &             +cff * (psi_s(i+1,j,k) -psi_s(i,j,k))
     &              *adp*(pm(i,j-1)+pm(i,j))

            v_s(i,j,k) = v_s(i,j,k) + v_ap(i,j,k)
# else
            v_ap(i,j,k) =
     &             +cff * (psi_s(i+1,j,k) -psi_s(i,j,k))
     &              *(pm(i,j-1)+pm(i,j))

            v_s(i,j,k) = v_s(i,j,k) + v_ap(i,j,k)
# endif
#else
# if defined DAS_ADIMPSICHI
            v_s(i,j,k) = v_s(i,j,k)
     &             +cff * (psi_s(i+1,j,k) -psi_s(i,j,k))
     &              *adp*(pm(i,j-1)+pm(i,j))
# else
            v_s(i,j,k) = v_s(i,j,k)
     &             +cff * (psi_s(i+1,j,k) -psi_s(i,j,k))
     &              *(pm(i,j-1)+pm(i,j))

# endif
#endif
          enddo
        enddo
!
!........
!
! chi -> u
!
        do j=Jstr,Jend
          do i=Istr,IendR
#if defined DAS_UVAGEO_DIAG
# if defined DAS_ADIMPSICHI
            u_ac(i,j,k)=
     &               cff * (chi_s(i,j,k) -chi_s(i-1,j,k))
     &              *adp*(pm(i-1,j)+pm(i,j))
            u_s(i,j,k)= u_s(i,j,k) + u_ac(i,j,k)
# else
            u_ac(i,j,k)=
     &               cff * (chi_s(i,j,k) -chi_s(i-1,j,k))
     &              *(pm(i-1,j)+pm(i,j))
            u_s(i,j,k)= u_s(i,j,k) + u_ac(i,j,k)
# endif
#else
# if defined DAS_ADIMPSICHI
            u_s(i,j,k)= u_s(i,j,k)+
     &               cff * (chi_s(i,j,k) -chi_s(i-1,j,k))
     &              *adp*(pm(i-1,j)+pm(i,j))
# else
            u_s(i,j,k)= u_s(i,j,k)+
     &               cff * (chi_s(i,j,k) -chi_s(i-1,j,k))
     &              *(pm(i-1,j)+pm(i,j))
# endif
#endif
          enddo
        enddo
!
! chi -> v
!
        do j=Jstr,JendR
          do i=Istr,Iend
#if defined DAS_UVAGEO_DIAG
# if defined DAS_ADIMPSICHI
            v_ac(i,j,k)=
     &               cff*(chi_s(i,j,k) - chi_s(i,j-1,k))
     &              *adp*(pn(i,j-1)+pn(i,j))
            v_s(i,j,k)=v_s(i,j,k) + v_ac(i,j,k)
# else
            v_ac(i,j,k)=
     &               cff*(chi_s(i,j,k) - chi_s(i,j-1,k))
     &              *(pn(i,j-1)+pn(i,j))
            v_s(i,j,k)=v_s(i,j,k) + v_ac(i,j,k)
# endif
#else
# if defined DAS_ADIMPSICHI
            v_s(i,j,k)=v_s(i,j,k)+
     &               cff*(chi_s(i,j,k) - chi_s(i,j-1,k))
     &              *adp*(pn(i,j-1)+pn(i,j))
# else
            v_s(i,j,k)=v_s(i,j,k)+
     &               cff*(chi_s(i,j,k) - chi_s(i,j-1,k))
     &              *(pn(i,j-1)+pn(i,j))
# endif
#endif
          enddo
        enddo
!      
      enddo  !k
!          
      return
      end
