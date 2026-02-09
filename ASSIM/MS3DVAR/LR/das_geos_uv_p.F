#include "cppdefs.h"

      subroutine das_geos_uv_p (tile)
      implicit none
      integer tile,trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "das_private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_geos_uv_p_tile (Istr,Iend,Jstr,Jend,
     &                     A2d(1,1,trd),A2d(1,2,trd),
     &                     A2d(1,3,trd))
      return
      end

      subroutine das_geos_uv_p_tile (Istr,Iend,Jstr,Jend,ug,vg)
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
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
# include "das_ocean_smooth.h"
      real ug(PRIVATE_2D_SCRATCH_ARRAY),
     &     vg(PRIVATE_2D_SCRATCH_ARRAY)
      real cff, ff0, crit_lat, crit, rho02, fcr
!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
!
      cff=2.0*g
      rho02=2.0*rho0
      
      crit_lat=2.5
      crit=2.0*1.e-5*sin(crit_lat*3.1415926/180.0)    ! 5 -> 1.7431e-06
                                                      ! 3 -> 1.0467e-06
      ff0=1.5e-7    ! 0.45
!      print*,'crit_lat=',crit_lat
!
! NOTE: p=P/g
!  u
!
      do k=1,NDAS
        do j=Jstr,Jend+1
          do i=Istr-1,IendR
            if (abs(f(i,j-1)+f(i,j)) .lt. crit) then
              fcr=f(i,j-1)+f(i,j)
              if (abs(fcr) .lt. ff0) then     ! 0.15 degrees
                ug(i,j)=0.0
              else
                ug(i,j)=-cff*(p_s(i,j,k) -p_s(i,j-1,k))
     &            *(pn(i,j-1)+pn(i,j))/fcr
     &             /(rho02+rho_das9(i,j,k)+rho_das9(i,j-1,k))
                ug(i,j)=ug(i,j)*abs(fcr)/crit
              endif
            else
              ug(i,j)=-cff*(p_s(i,j,k) -p_s(i,j-1,k))
     &            *(pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
     &             /(rho02+rho_das9(i,j,k)+rho_das9(i,j-1,k))
            endif
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,IendR
            u_s(i,j,k)=0.25*(
     &        ug(i-1,j)+ug(i-1,j+1)+ug(i,j)+ug(i,j+1) )
          enddo
        enddo
        if (SOUTHERN_EDGE) then
          do i=Istr,IendR
            u_s(i,0,k) = 0.5*(ug(i-1,1)+ug(i,1))
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=Istr,IendR
            u_s(i,Mm+1,k) = 0.5*(ug(i-1,Mm+1)+ug(i,Mm+1))
          enddo
        endif
!
!  v
!
        do j=Jstr-1,JendR
          do i=Istr,Iend+1
            if (abs(f(i-1,j)+f(i,j)) .lt. crit) then
              fcr=f(i-1,j)+f(i,j)
              if (abs(fcr) .lt. ff0) then     ! 0.15 degrees
                vg(i,j)=0.0
              else
                vg(i,j)= cff * (p_s(i,j,k) -p_s(i-1,j,k))
     &            *(pm(i-1,j)+pm(i,j))/fcr
     &            /(rho02+rho_das9(i-1,j,k)+rho_das9(i,j,k))
                vg(i,j)=vg(i,j)*abs(fcr)/crit
               endif
            else
              vg(i,j)= cff * (p_s(i,j,k) -p_s(i-1,j,k))
     &            *(pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
     &            /(rho02+rho_das9(i-1,j,k)+rho_das9(i,j,k))
            endif
          enddo
        enddo
        do j=Jstr,JendR
          do i=Istr,Iend
            v_s(i,j,k)=0.25* (
     &        vg(i,j-1)+vg(i,j)+vg(i+1,j-1)+vg(i+1,j)  )
          enddo
        enddo
        if (WESTERN_EDGE) then
          do j=Jstr,JendR
            v_s(0,j,k) = 0.5*(vg(1,j-1)+vg(1,j))
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=Jstr,JendR
            v_s(Lm+1,j,k) = 0.5*(vg(Lm+1,j-1)+vg(Lm+1,j))
          enddo
        endif
      enddo   !k
     
      return
      end
