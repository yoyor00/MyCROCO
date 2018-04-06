      subroutine cost_fun(x, cost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "adparam.h"

# ifdef MPI
# include "mpi_roms.h"
# ifdef AMPI
# include "ampi/ampif.h"
# else
      include 'mpif.h'
# endif
# endif

      double precision cost

      double precision acost(NNODES)
      double precision rdrg2_bck

      integer k, ierr, i

      rdrg2_bck = rdrg2
      call set_state(x)

      cost = 0
      do k=1,120
C$AD NOCHECKPOINT
         call step()
         call cost_fun_step(cost, k)
      end do

      write (*,*) NNODES

# ifdef MPI
# ifdef AMPI
c     issue with Reduce (bws_ampi_reduce is present but signature
c     different, association by name, from what tapenade writes) call
c     call AMPI_Reduce(icost, cost, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, & 0,
c     MPI_COMM_WORLD, ierr);
      if (mynode.eq.0) then
         do i = 1, NNODES
            acost(i) = 0.
         end do
      end if
      call AMPI_Gather(cost, 1, AMPI_ADOUBLE_PRECISION, acost, 1,
     &     AMPI_ADOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

c     cost = 0
      if (NNODES.ge.2) then
         if (mynode.eq.0) then
            do i = 2, NNODES
               cost = cost + acost(i)
            end do
         end if
      end if

# endif
# else

# endif

c      MPI_master_only cost = cost +
c     &     (x(1)-rdrg2_bck)*(x(1)-rdrg2_bck)

      write (*,*) 'rdrg2=',rdrg2

      write (*,*) 'cost=',cost

      return
      end subroutine cost_fun

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

      subroutine cost_fun_step_2d(tile, icost, k)

      implicit none
# include "param.h"
      integer k
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, k)

      end subroutine cost_fun_step_2d

      subroutine cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, k)
      implicit none
#include "param.h"
#include "grid.h"
#include "scalars.h"
#include "ocean2d.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend,i,j
      integer k
      double precision delta_z, icost
#include "compute_extended_bounds.h"

      do j=Jstr,Jend
         do i=Istr,Iend
c            write(*,*) i,j,knew,sim_iicroot+k,
c     &           zeta(i,j,knew),obs(i,j,sim_iicroot+k)

            delta_z = (zeta(i,j,knew)-obs(i,j,sim_iicroot+k))
            icost = icost + delta_z*delta_z
         end do
      end do

      end subroutine cost_fun_step_2d_tile


      subroutine set_state(x)
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
         call set_state_2d(tile, x)
      end do

      end subroutine set_state


      subroutine set_state_2d(tile, x)

      implicit none
# include "param.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, x)

      end subroutine set_state_2d


      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend,x)
      implicit none

#include "param.h"
#include "grid.h"
#include "ocean2d.h"
#include "scalars.h"
#include "mpi_roms.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend, i,j
      double precision rdrg2_bck, h0(GLOBAL_2D_ARRAY)

      integer bi, bj
#include "compute_extended_bounds.h"

      bi = IendR-IstrR+1
      bj = JendR-JstrR+1

      write(*,*) 'set_state'
      write(*,*) 'krhs=',krhs

      write(*,*) bi, bj, ad_array_size

c      do j=JstrR,JendR
c         do i=IstrR,IendR
c            write (*,*) '###',i,j,h(i,j),
c     &           1+(i-IstrR+1)+bi*(j-JstrR),
c     &           x(1+(i-IstrR+1)+bi*(j-JstrR))
c             h(i,j) = x(1+(i-IstrR+1)+bi*(j-JstrR))
c         end do
c      end do

      rdrg2 = real(x(1))

      write(*, *) 'rdrg2=',rdrg2
      write(*, *) 'IstrR=',IstrR
      write(*, *) 'JstrR=',JstrR
      write(*, *) 'h(IstrR,JstrR)=',h(IstrR,JstrR)

 900  format(1PE16.9)

      return
      end subroutine set_state_2d_tile
