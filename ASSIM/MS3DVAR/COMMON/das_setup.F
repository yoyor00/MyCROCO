#include "cppdefs.h"

      subroutine das_setup (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_setup_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_setup_tile (Istr,Iend,Jstr,Jend)
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
      if (WEST_INTER) IstrR=Istr   
      if (EAST_INTER) IendR=Iend  
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
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE
!
      return
      end

