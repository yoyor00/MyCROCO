#include "cppdefs.h"

      subroutine das_add_lr (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_add_lr_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_add_lr_tile (Istr,Iend,Jstr,Jend)
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
!  2-D variables.
!
!      do j=JstrR,JendR
!        do i=IstrR,IendR
!          zeta_das(i,j)=zeta_das(i,j) + zeta_lr_s(i,j)
!        enddo
!      enddo
!
!  Initialize 3-D primitive variables.
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=Istr,IendR
            u_das(i,j,k)=u_das(i,j,k) + u_lr_s(i,j,k)
          enddo
        enddo
        do j=Jstr,JendR
          do i=IstrR,IendR
            v_das(i,j,k)=v_das(i,j,k) + v_lr_s(i,j,k)
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_das(i,j,k,itemp)=t_das(i,j,k,itemp) + t_lr_s(i,j,k,itemp)
            t_das(i,j,k,isalt)=t_das(i,j,k,isalt) + t_lr_s(i,j,k,isalt)
          enddo
        enddo
      enddo    !k
      return
      end
