#include "cppdefs.h"

      subroutine das_adj_geos_uv_p (tile)
      implicit none
      integer tile,trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "das_private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_adj_geos_uv_p_tile (Istr,Iend,Jstr,Jend,
     &                     A2d(1,1,trd),A2d(1,2,trd),
     &                     A3d(1,1,trd))
      return
      end

      subroutine das_adj_geos_uv_p_tile (Istr,Iend,Jstr,Jend,
     &               ug_adj,vg_adj,ps_adj)
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
      real cff, ff0, crit_lat, crit, rho02, fcr
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
# include "das_ocean_smooth.h"
      real ug_adj(PRIVATE_2D_SCRATCH_ARRAY),
     &     vg_adj(PRIVATE_2D_SCRATCH_ARRAY),
     &     ps_adj(PRIVATE_2D_SCRATCH_ARRAY,NDAS)
!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
!
      cff=2.0*g
      rho02=2.0*rho0
      
      crit_lat=2.5
      crit=2.0*1.e-5*sin(crit_lat*3.1415926/180.0)

      ff0=1.2e-7    ! 0.16
!
! adjoint
!
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          ug_adj(i,j)=0.0
          vg_adj(i,j)=0.0
        enddo
      enddo
      do k=1,NDAS
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            ps_adj(i,j,k)=0.0
          enddo
        enddo
      enddo
!
! NOTE: p=P/g
!  u
!
      do k=1,NDAS

        if (NORTHERN_EDGE) then                                                                                                          
          do i=Istr,IendR                                                                                                                
            ug_adj(i-1,Mm+1)=ug_adj(i-1,Mm+1)
     &           +0.5*u_s_adj(i,Mm+1,k) 
            ug_adj(i,Mm+1)=ug_adj(i,Mm+1)
     &           +0.5*u_s_adj(i,Mm+1,k) 
            u_s_adj(i,Mm+1,k) = 0.
          enddo                                                                                                                          
        endif
        if (SOUTHERN_EDGE) then
          do i=Istr,IendR                                                                                                                
            ug_adj(i-1,1)=ug_adj(i-1,1)
     &              +0.5*u_s_adj(i,0,k)
            ug_adj(i,1)=ug_adj(i,1)
     &              +0.5*u_s_adj(i,0,k)
            u_s_adj(i,0,k) = 0.
          enddo                                                                                                                          
        endif
        do j=Jstr,Jend                                                                                                                   
          do i=Istr,IendR                                                                                                                
            ug_adj(i-1,j)=ug_adj(i-1,j)
     &              +0.25*u_s_adj(i,j,k)
            ug_adj(i-1,j+1)=ug_adj(i-1,j+1)
     &              +0.25*u_s_adj(i,j,k)
            ug_adj(i,j)=ug_adj(i,j)
     &              +0.25*u_s_adj(i,j,k)
            ug_adj(i,j+1)=ug_adj(i,j+1)
     &              +0.25*u_s_adj(i,j,k)
            u_s_adj(i,j,k)=0.
          enddo                                                                                                                          
        enddo 
        do j=Jstr,Jend+1
          do i=Istr-1,IendR
            fcr=f(i,j-1)+f(i,j)
            if (abs(f(i,j-1)+f(i,j)) .lt. crit) then
              if (abs(fcr) .lt. ff0) then     ! 0.15 degrees
                ug_adj(i,j)=0.0
              else
                ug_adj(i,j)=ug_adj(i,j)*abs(fcr)/crit
                ps_adj(i,j,k)=ps_adj(i,j,k)  
     &            -cff*ug_adj(i,j)*(pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
     &            /(rho02+rho_das9(i,j,k)+rho_das9(i,j-1,k))
                ps_adj(i,j-1,k)=ps_adj(i,j-1,k)
     &             +cff*ug_adj(i,j)*(pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
     &            /(rho02+rho_das9(i,j,k)+rho_das9(i,j-1,k))
                ug_adj(i,j)=0.
              endif
            else
              ps_adj(i,j,k)=ps_adj(i,j,k)
     &              -cff*ug_adj(i,j)*(pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
     &            /(rho02+rho_das9(i,j,k)+rho_das9(i,j-1,k))
              ps_adj(i,j-1,k)=ps_adj(i,j-1,k)
     &              +cff*ug_adj(i,j)*(pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
     &            /(rho02+rho_das9(i,j,k)+rho_das9(i,j-1,k))
              ug_adj(i,j)=0.
            endif
          enddo
        enddo
!
!  v
!
        if (EASTERN_EDGE) then                                                                                                           
          do j=Jstr,JendR                                                                                                                
            vg_adj(Lm+1,j-1)=vg_adj(Lm+1,j-1)
     &             +0.5*v_s_adj(Lm+1,j,k)
            vg_adj(Lm+1,j)=vg_adj(Lm+1,j)
     &             +0.5*v_s_adj(Lm+1,j,k)
            v_s_adj(Lm+1,j,k) = 0.
          enddo                                                                                                                          
        endif
        if (WESTERN_EDGE) then                                                                                                           
          do j=Jstr,JendR                                                                                                                
            vg_adj(1,j-1)=vg_adj(1,j-1)
     &             +0.5*v_s_adj(0,j,k)
            vg_adj(1,j)=vg_adj(1,j)
     &             +0.5*v_s_adj(0,j,k)
            v_s_adj(0,j,k) = 0.
          enddo                                                                                                                          
        endif 
        do j=Jstr,JendR                                                                                                                  
          do i=Istr,Iend                                                                                                                 
            vg_adj(i,j-1)=vg_adj(i,j-1)
     &             +0.25* v_s_adj(i,j,k)
            vg_adj(i,j)=vg_adj(i,j)
     &             +0.25* v_s_adj(i,j,k)
            vg_adj(i+1,j-1)=vg_adj(i+1,j-1)
     &             +0.25* v_s_adj(i,j,k)
            vg_adj(i+1,j)=vg_adj(i+1,j)
     &             +0.25* v_s_adj(i,j,k)
            v_s_adj(i,j,k)=0.
          enddo                                                                                                                          
        enddo 
        do j=Jstr-1,JendR
          do i=Istr,Iend+1
            if (abs(f(i-1,j)+f(i,j)) .lt. crit ) then
              fcr=f(i-1,j)+f(i,j)
              if (abs(fcr) .lt. ff0) then     ! 0.15 degrees
                vg_adj(i,j)=0.0
              else
                vg_adj(i,j)=vg_adj(i,j)*abs(fcr)/crit
                ps_adj(i,j,k)=ps_adj(i,j,k)
     &            +cff * vg_adj(i,j)*(pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
     &            /(rho02+rho_das9(i-1,j,k)+rho_das9(i,j,k))
                ps_adj(i-1,j,k)=ps_adj(i-1,j,k)
     &            -cff * vg_adj(i,j)*(pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
     &            /(rho02+rho_das9(i-1,j,k)+rho_das9(i,j,k))
                vg_adj(i,j)= 0.
              endif
            else
              ps_adj(i,j,k)=ps_adj(i,j,k)
     &          +cff * vg_adj(i,j)*(pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
     &            /(rho02+rho_das9(i-1,j,k)+rho_das9(i,j,k))
              ps_adj(i-1,j,k)=ps_adj(i-1,j,k)
     &          -cff * vg_adj(i,j)*(pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
     &            /(rho02+rho_das9(i-1,j,k)+rho_das9(i,j,k))
              vg_adj(i,j)= 0.
            endif
          enddo
        enddo
      enddo   !k
C$OMP CRITICAL (adj_geosp_cr)
      do k=1,NDAS
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            p_s_adj(i,j,k)=p_s_adj(i,j,k)+ps_adj(i,j,k)
            ps_adj(i,j,k)=0.0
          enddo
        enddo
      enddo
C$OMP END CRITICAL (adj_geosp_cr)
!     
      return
      end
