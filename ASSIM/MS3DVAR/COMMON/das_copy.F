#include "cppdefs.h"

      subroutine das_copy (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_copy_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_copy_tile (Istr,Iend,Jstr,Jend)
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
! 2D
!
        do j=JstrR,JendR
          do i=IstrR,IendR
            zeta_das(i,j)=zeta_s(i,j)
          enddo
        enddo
!
!  Initialize 3-D primitive variables.
!
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_das(i,j,k,itemp)=t_s(i,j,k,itemp)
              t_das(i,j,k,isalt)=t_s(i,j,k,isalt)
            enddo
          enddo
        enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=Istr,IendR
              u_das(i,j,k)=u_s(i,j,k)
            enddo
          enddo
        enddo
        do k=1,NDAS
          do j=Jstr,JendR
            do i=IstrR,IendR
              v_das(i,j,k)=v_s(i,j,k)
            enddo
          enddo
        enddo

      return
      end
