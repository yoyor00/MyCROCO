#include "cppdefs.h"

      subroutine das_setup1 (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_setup1_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_setup1_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc, kmin
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "ocean3d.h"
#include "grid.h"
#include "das_ocean.h"
!
      real init    !!!!  0xFFFA5A5A ==> NaN
      parameter (init=0.)
      real r1
!
# include "compute_extended_bounds.h"
!
!  Initialize 2-D primitive variables.
!
! for z -> s
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          kmin=ndas
          do k=1,ndas
            if (z_das(kmin) .ge. z_r(i,j,1)) then
              kmin=kmin-1
            else
              nzr_das(i,j)=kmin+1
              goto 10
            endif
          enddo
          nzr_das(i,j)=1
  10      continue
        enddo
      enddo
!
! deepest level above nzr_das
! ns_das >= 1
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          kmin=N
          do k=1,N
            if (z_r(i,j,kmin) .ge. z_das(nzr_das(i,j)) ) then
              kmin=kmin-1
            else
              nsr_das(i,j)=kmin+1
              goto 20
            endif
          enddo
          nsr_das(i,j)=1
  20      continue
        enddo
      enddo
!
!  Initialize 3-D primitive variables.
!
!  rmask_das, umask_das, and vmask_das have been
!  initialize to be zero
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=nzr_das(i,j),ndas
            rmask_das(i,j,k)=rmask(i,j)
          enddo
        enddo
      enddo

!
      return
      end

