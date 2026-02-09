#include "cppdefs.h"

      subroutine das_adj_chfradar (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_chfradar_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_chfradar_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, fl_ot
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
# include "das_innov.h"
!
# include "compute_auxiliary_bounds.h"

!
! u
!
      do j=Jstr,Jend
         do i=Istr,IendR
           cft=u_s(i,j,NDAS)-u_hf(i,j)
           u_s_adj(i,j,NDAS)=u_s_adj(i,j,NDAS)
     &                     +cft*hf_umask(i,j)*hf_ouv
         enddo
       enddo
!
! v
!
      do j=Jstr,JendR
         do i=Istr,Iend
           cft=v_s(i,j,NDAS)-v_hf(i,j)
           v_s_adj(i,j,NDAS)=v_s_adj(i,j,NDAS)
     &                     +cft*hf_vmask(i,j)*hf_ouv
         enddo
       enddo
!
!6KM
!
! u
!
      do j=Jstr,Jend
         do i=Istr,IendR
           cft=u_s(i,j,NDAS)-u_hf6(i,j)
           u_s_adj(i,j,NDAS)=u_s_adj(i,j,NDAS)
     &                     +cft*hf_umask6(i,j)*hf_uierr6(i,j)
         enddo
       enddo
!
! v
!
      do j=Jstr,JendR
         do i=Istr,Iend
           cft=v_s(i,j,NDAS)-v_hf6(i,j)
           v_s_adj(i,j,NDAS)=v_s_adj(i,j,NDAS)
     &                     +cft*hf_vmask6(i,j)*hf_vierr6(i,j)
         enddo
       enddo
!
      return
      end
