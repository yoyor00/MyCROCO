#include "cppdefs.h"

      subroutine das_ca2zero (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_ca2zero_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_ca2zero_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "das_ocean.h"

      real init    !!!!  0xFFFA5A5A ==> NaN

      parameter (init=0.)

!
# include "compute_extended_bounds.h"
!
!
!  2-D variables.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_das(i,j)=init
        enddo
      enddo
!
!  Initialize 3-D primitive variables.
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            u_das(i,j,k)=init
            v_das(i,j,k)=init
          enddo
        enddo
      enddo
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_das(i,j,k,itemp)=init
              t_das(i,j,k,isalt)=init
            enddo
          enddo
        enddo
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            psi_das(i,j,k)=init
            chi_das(i,j,k)=init
          enddo
        enddo
      enddo
      return
      end
