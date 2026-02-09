#include "cppdefs.h"

      subroutine das_adj_psichi_uv (tile)
      implicit none
      integer tile,trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "das_private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_adj_psichi_uv_tile (Istr,Iend,Jstr,Jend,
     &                     A3d(1,1,trd),A3d(1,2,trd))
      return
      end

      subroutine das_adj_psichi_uv_tile (Istr,Iend,Jstr,Jend,
     &                                     psis_adj,chis_adj)
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
     &           ,NSUB
      real cff
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
      real 
     &     psis_adj(PRIVATE_2D_SCRATCH_ARRAY,NDAS),
     &     chis_adj(PRIVATE_2D_SCRATCH_ARRAY,NDAS)
!
# include "compute_auxiliary_bounds.h"
!
!------------------------------------------------------------------
!
      cff=0.5
!
      do k=1,NDAS
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            psis_adj(i,j,k)=0.0
            chis_adj(i,j,k)=0.0
          enddo
        enddo
      enddo
!
      do k=1,NDAS

!        do j=Jstr,Jend                                                                                       
!          do i=Istr,IendR                                                                                    
!            u_s_adj(i,j,k) = u_s_adj(i,j,k)                                                                          
!     &               *pmask_das(i,j,k)*pmask_das(i,j+1,k) !exclude one grid at physical boundary             
!          enddo                                                                                              
!        enddo                                                                                                
!!                                                                                                            
!! psi -> v                                                                                                   
!!                                                                                                            
!        do j=Jstr,JendR                                                                                      
!          do i=Istr,Iend                                                                                     
!            v_s_adj(i,j,k) = v_s_adj(i,j,k)                                                                          
!     &              *pmask_das(i,j,k)*pmask_das(i+1,j,k)  !exclude one grid at physical boundary             
!          enddo                                                                                              
!        enddo 
!
!
! chi -> v
!
        do j=Jstr,JendR
          do i=Istr,Iend
# if defined DAS_ADIMPSICHI
            chis_adj(i,j,k)=chis_adj(i,j,k)+cff*v_s_adj(i,j,k)
     &              *adp*(pn(i,j-1)+pn(i,j))
            chis_adj(i,j-1,k)=chis_adj(i,j-1,k)-cff*v_s_adj(i,j,k)
     &              *adp*(pn(i,j-1)+pn(i,j))
# else
            chis_adj(i,j,k)=chis_adj(i,j,k)+cff*v_s_adj(i,j,k)
     &              *(pn(i,j-1)+pn(i,j))
            chis_adj(i,j-1,k)=chis_adj(i,j-1,k)-cff*v_s_adj(i,j,k)
     &              *(pn(i,j-1)+pn(i,j))
# endif
          enddo
        enddo
!
! chi -> u
!
        do j=Jstr,Jend
          do i=Istr,IendR
# if defined DAS_ADIMPSICHI
            chis_adj(i,j,k)=chis_adj(i,j,k)+cff *u_s_adj(i,j,k)
     &              *adp*(pm(i-1,j)+pm(i,j))
            chis_adj(i-1,j,k)=chis_adj(i-1,j,k)-cff *u_s_adj(i,j,k)
     &              *adp*(pm(i-1,j)+pm(i,j))
# else
            chis_adj(i,j,k)=chis_adj(i,j,k)+cff *u_s_adj(i,j,k)
     &              *(pm(i-1,j)+pm(i,j))
            chis_adj(i-1,j,k)=chis_adj(i-1,j,k)-cff *u_s_adj(i,j,k)
     &              *(pm(i-1,j)+pm(i,j))
# endif
          enddo
        enddo
!
! psi -> v
!
        do j=Jstr,JendR
          do i=Istr,Iend
# if defined DAS_ADIMPSICHI
            psis_adj(i+1,j,k)=psis_adj(i+1,j,k)+ cff * v_s_adj(i,j,k)
     &                 *adp*(pm(i,j-1)+pm(i,j))
            psis_adj(  i,j,k)=psis_adj(  i,j,k)- cff * v_s_adj(i,j,k)
     &                 *adp*(pm(i,j-1)+pm(i,j))
# else
            psis_adj(i+1,j,k)=psis_adj(i+1,j,k)+ cff * v_s_adj(i,j,k)
     &                 *(pm(i,j-1)+pm(i,j))
            psis_adj(  i,j,k)=psis_adj(  i,j,k)- cff * v_s_adj(i,j,k)
     &                 *(pm(i,j-1)+pm(i,j))
# endif
          enddo
        enddo
!
!  psi -> u
!
!
        do j=Jstr,Jend
          do i=Istr,IendR
# if defined DAS_ADIMPSICHI
            psis_adj(i,j+1,k)=psis_adj(i,j+1,k)-cff*u_s_adj(i,j,k)
     &                  *adp*(pn(i-1,j)+pn(i,j))
            psis_adj(i,j,k)=psis_adj(i,j,k)    +cff*u_s_adj(i,j,k)
     &                  *adp*(pn(i-1,j)+pn(i,j))
# else
            psis_adj(i,j+1,k)=psis_adj(i,j+1,k)-cff*u_s_adj(i,j,k)
     &                  *(pn(i-1,j)+pn(i,j))
            psis_adj(i,j,k)=psis_adj(i,j,k)    +cff*u_s_adj(i,j,k)
     &                  *(pn(i-1,j)+pn(i,j))
# endif
          enddo
        enddo
!      
      enddo  !k
!
! serailized
!
C$OMP CRITICAL (adj_psichi_cr)
!
      do k=1,NDAS
        do j=Jstr,Jend+1
          do i=Istr,Iend+1
            psi_s_adj(i,j,k)=psi_s_adj(i,j,k)+psis_adj(i,j,k)
            psis_adj(i,j,k)=0.0
          enddo
        enddo
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            chi_s_adj(i,j,k)=chi_s_adj(i,j,k)+chis_adj(i,j,k)
            chis_adj(i,j,k)=0.0
          enddo
        enddo
      enddo
!      do k=1,NDAS
!        do j=Jstr,Jend+1
!          do i=Istr,Iend+1
!            psi_s_adj(i,j,k)=psi_s_adj(i,j,k)*pmask_das(i,j,k)
!          enddo
!        enddo
!        do j=Jstr-1,Jend+1
!          do i=Istr-1,Iend+1
!            chi_s_adj(i,j,k)=chi_s_adj(i,j,k)*rmask_das(i,j,k)
!          enddo
!        enddo
!      enddo
!
C$OMP END CRITICAL (adj_psichi_cr)      
!          
      return
      end
