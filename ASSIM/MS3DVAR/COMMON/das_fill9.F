#include "cppdefs.h"

      subroutine das_fill9 (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_fill9_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_fill9_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_ocean.h"
#include "das_ocean9.h"
!
# include "compute_extended_bounds.h"
!
! The vertical interpolation may generate bad apples in T/S
! A protective check is needed
!
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
             if ( t_das(i,j,k,itemp) .lt. 1.99999 )
     &            t_das(i,j,k,itemp)=1.99999
             if ( t_das(i,j,k,isalt) .lt. 5.99999 )
     &            t_das(i,j,k,isalt)=5.99999
            enddo
          enddo
        enddo
!
! 2D
!
        do j=JstrR,JendR
          do i=IstrR,IendR
            zeta_das(i,j)=zeta_das(i,j)*rmask_das(i,j,NDAS)
          enddo
        enddo
!
!  Initialize 3-D primitive variables.
!
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_das(i,j,k,itemp)=t_das(i,j,k,itemp)
     &                  +(1.0-rmask_das(i,j,k))*tc1d_das(k)
              t_das(i,j,k,isalt)=t_das(i,j,k,isalt)
     &                  +(1.0-rmask_das(i,j,k))*sc1d_das(k)
            enddo
          enddo
        enddo
      return
      end
