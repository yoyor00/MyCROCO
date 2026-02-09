#include "cppdefs.h"

      subroutine das_setup2 (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_setup2_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_setup2_tile (Istr,Iend,Jstr,Jend)
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
#ifdef EW_PERIODIC
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE  Istr,IendR
# ifdef MPI                            
      if (WEST_INTER) IstrR=Istr          ! computational boundary
      if (EAST_INTER) IendR=Iend          ! are filled during
# endif                               
#endif                               
#ifdef NS_PERIODIC
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE  Jstr,JendR
# ifdef MPI
      if (SOUTH_INTER) JstrR=Jstr         ! same as above.
      if (NORTH_INTER) JendR=Jend         !
# endif
#endif
!
! psi
!
!NOTE, pmask_das is DIFFERENT FROM pmask even in the horizontal
!
      do j=JV_RANGE
        do i=IU_RANGE
          do k=1,ndas
            pmask_das(i,j,k) = rmask_das(i,j,k) * rmask_das(i-1,j,k)
     &         * rmask_das(i,j-1,k) * rmask_das(i-1,j-1,k)
!            if (r1 .lt. 1.5) then
!              pmask_das(i,j,k)=0.0
!            else
!              pmask_das(i,j,k)=1.0
!            endif
          enddo
        enddo
      enddo
!
! u
!
! 1. mask
!
      do j=JR_RANGE
        do i=IU_RANGE
          do k=1,ndas
            umask_das(i,j,k)=rmask_das(i,j,k)*rmask_das(i-1,j,k)
          enddo
        enddo
      enddo
!
! 2. nzu_das 
!
      do j=JR_RANGE
        do i=IU_RANGE
          do k=1,ndas
            if (umask_das(i,j,k) .gt. 0.5) then
              nzu_das(i,j)=k
              goto 30
            else
              nzu_das(i,j)=ndas
            endif
          enddo
  30      continue
        enddo
      enddo
!
! 3. nsu_das >=1 
!
      do j=JR_RANGE
        do i=IU_RANGE
          kmin=N
          do k=1,N
            r1=0.5*(z_r(i-1,j,kmin)+z_r(i,j,kmin))
            if (r1 .ge. z_das(nzu_das(i,j)) ) then
              kmin=kmin-1
            else
              nsu_das(i,j)=kmin+1
              goto 40
            endif
          enddo
          nsu_das(i,j)=1
  40      continue
        enddo
      enddo
!...
!
! v
!
! 1. mask
!
      do j=JV_RANGE
        do i=IR_RANGE
          do k=1,ndas
            vmask_das(i,j,k)=rmask_das(i,j,k)*rmask_das(i,j-1,k)
          enddo
        enddo
      enddo
!
! 2. nzv_das
!
      do j=JV_RANGE
        do i=IR_RANGE
          do k=1,ndas
            if (vmask_das(i,j,k) .gt. 0.5) then
              nzv_das(i,j)=k
              goto 50
            else
              nzv_das(i,j)=ndas
            endif
          enddo
  50      continue
        enddo
      enddo
!
! nsv_das >= 1
!
      do j=JV_RANGE
        do i=IR_RANGE
          kmin=N
          do k=1,N
            r1=0.5*(z_r(i,j-1,kmin)+z_r(i,j,kmin))
            if (r1 .ge. z_das(nzv_das(i,j)) ) then
              kmin=kmin-1
            else
              nsv_das(i,j)=kmin+1
              goto 60
            endif
          enddo
          nsv_das(i,j)=1
  60      continue
        enddo
      enddo
!
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE
!
      return
      end

