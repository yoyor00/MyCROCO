#include "cppdefs.h"

      subroutine das_init_arrays (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_init_arrays_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_init_arrays_tile (Istr,Iend,Jstr,Jend)
!
! Initialize (first touch) globally accessable arrays. Most of them
! are assigned to zeros, vertical mixing coefficients are assinged
! to their background values. These will remain unchenged if no
! vertical mixing scheme is applied. Because of the "first touch"
! default data distribution policy, this operation actually performs
! distribution of the shared arrays accross the cluster, unless
! another distribution policy is specified to override the default.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "grid.h"
#include "das_ocean.h"

      real init    !!!!  0xFFFA5A5A ==> NaN

      parameter (init=0.)

!
# include "compute_extended_bounds.h"
!
!
!  Initialize 2-D primitive variables.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_das(i,j)=init
          zeta_s(i,j)=init
          nsr_das(i,j)=init
          nsu_das(i,j)=init
          nsv_das(i,j)=init
          nzr_das(i,j)=init
          nzu_das(i,j)=init
          nzv_das(i,j)=init
        enddo
      enddo
#ifdef SOLVE3D
!
!  Initialize 3-D primitive variables.
!

      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            u_das(i,j,k)=init
            v_das(i,j,k)=init
            u_s(i,j,k)=init
            v_s(i,j,k)=init
            rho_s(i,j,k)=init
# ifdef MASKING
            rmask_das(i,j,k)=init
            umask_das(i,j,k)=init
            vmask_das(i,j,k)=init
#  endif
          enddo
        enddo
      enddo
      do itrc=1,NT
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_das(i,j,k,itrc)=init
              t_s(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
#endif /* SOLVE3D */
!
! initialize those not innitialized in ROMS
!
#ifdef MASKING
      do j=JstrR,JendR
        do i=IstrR,IendR
          rmask(i,j)=init
          umask(i,j)=init
          vmask(i,j)=init
        enddo
      enddo
#endif
      return
      end
