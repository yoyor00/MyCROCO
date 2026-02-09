#include "cppdefs.h"

      subroutine das_add_lr_s (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_add_lr_s_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_add_lr_s_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#include "param.h"
#include "das_param.h"
#include "das_param_lr.h"
#include "scalars.h"
#include "das_ocean.h"
#include "das_ocean_lr.h"

      real init    !!!!  0xFFFA5A5A ==> NaN

      parameter (init=0.)

!
# include "compute_extended_bounds.h"
!
!
!      do j=JstrR,JendR
!        do i=IstrR,IendR
!          zeta_s(i,j) = zeta_s(i,j) + zeta_lr_s(i,j)
!        enddo
!      enddo

      do k=1,NDAS
        do j=JstrR,JendR
          do i=Istr,IendR
            u_s(i,j,k)=u_s(i,j,k) + u_lr_s(i,j,k)
          enddo
        enddo
        do j=Jstr,JendR
          do i=IstrR,IendR
            v_s(i,j,k)=v_s(i,j,k) + v_lr_s(i,j,k)
          enddo
        enddo

        do j=JstrR,JendR
          do i=IstrR,IendR
            t_s(i,j,k,itemp)=t_s(i,j,k,itemp) + t_lr_s(i,j,k,itemp)
            t_s(i,j,k,isalt)=t_s(i,j,k,isalt) + t_lr_s(i,j,k,isalt)
          enddo
        enddo
      enddo    !k
      return
      end
