
      subroutine cost_fun(ad_x, cost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      double precision cost
      double precision xcost

      double precision Zobt_bck

      integer k, ierr, i

      Zobt_bck = Zobt
      call set_state(ad_x)

      cost = 0

C$AD BINOMIAL-CKP ad_nt+1 40 1
      do k=1,ad_nt
         call step()
         call cost_fun_step(cost, k)
      end do

c     cost = cost +
c     &     (x(1)-Zob_bck)*(x(1)-Zob_bck)

      write (*,*) 'Zob=',Zobt

      write (*,*) 'cost=',cost

      return
      end subroutine cost_fun

C$AD NOCHECKPOINT
      subroutine cost_fun_step(icost, k)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"
      double precision icost
      integer k
      integer ntrds,trd,my_first,my_last,tile,range

#ifdef OPENMP
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
#else
      ntrds=1
      trd=0
#endif
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      do tile=my_first,my_last
         call cost_fun_step_2d(tile, icost, k)
      end do

      end subroutine cost_fun_step

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d(tile, icost, k)

      implicit none
# include "param.h"
      integer k
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, k)

      end subroutine cost_fun_step_2d

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, k)
      implicit none
#include "param.h"
#include "grid.h"
#include "scalars.h"
#include "ocean2d.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend,i,j
      integer k, kk
      double precision delta_z, icost

      do j=Jstr,Jend
         do i=Istr,Iend
c            if ((zeta(i,j,knew)-ad_obs(i,j,k+1)) > 1D-10) then
c               write(*,*) i,j,knew,sim_iicroot,k,sim_iicroot+k,
c     &              zeta(i,j,knew),ad_obs(i,j,k+1)
c            end if
            delta_z = (zeta(i,j,knew)-ad_obs(i,j,k+1))
            icost = icost + delta_z*delta_z
         end do
      end do
      
      end subroutine cost_fun_step_2d_tile

C$AD NOCHECKPOINT
      subroutine set_state(ad_x)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"

      integer ntrds,trd,my_first,my_last,tile,range

#ifdef OPENMP
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
#else
      ntrds=1
      trd=0
#endif
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      do tile=my_first,my_last
         call set_state_2d(tile, ad_x)
      end do

      end subroutine set_state


C$AD NOCHECKPOINT
      subroutine set_state_2d(tile, ad_x)

      implicit none
# include "param.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, ad_x)

      end subroutine set_state_2d

C$AD NOCHECKPOINT
      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend,ad_x)
      implicit none

#include "param.h"
#include "grid.h"
#include "ocean2d.h"
#include "scalars.h"
#include "mpi_roms.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend, i,j
      double precision Zobt_bck, h0(GLOBAL_2D_ARRAY)

      integer bi, bj
#include "compute_extended_bounds.h"

      bi = IendR-IstrR+1
      bj = JendR-JstrR+1

c      write(*,*) 'set_state'
c      write(*,*) 'krhs=',krhs

c      write(*,*) bi, bj, ad_array_size

c      do j=JstrR,JendR
c         do i=IstrR,IendR
c            write (*,*) '###',i,j,h(i,j),
c     &           1+(i-IstrR+1)+bi*(j-JstrR),
c     &           x(1+(i-IstrR+1)+bi*(j-JstrR))
c             h(i,j) = x(1+(i-IstrR+1)+bi*(j-JstrR))
c         end do
c      end do

      Zobt = real(ad_x(1))
      zob(:,:) = Zobt
      
      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif
      write(*, *) 'Zob=',Zobt
      write(*, *) '*******************'

 900  format(1PE16.9)

      return
      end subroutine set_state_2d_tile
