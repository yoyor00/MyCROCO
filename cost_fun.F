      subroutine cost_fun(x, cost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
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
      real cost, icost

      real rdrg2_bck, h0(GLOBAL_2D_ARRAY), acost(NNODES)
      integer k, ierr, i

      call set_state(x, h0, rdrg2_bck)

      cost = 0
      do k=1,120
         call step()

         call cost_fun_step(cost)
         
      end do

# ifdef MPI
# ifdef AMPI
c     issue with Reduce (bws_ampi_reduce is present but signature
c     different, association by name, from what tapenade writes) call
c     call AMPI_Reduce(icost, cost, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, & 0,
c     MPI_COMM_WORLD, ierr);
      call AMPI_Gather(cost, 1, AMPI_ADOUBLE_PRECISION, acost, NNODES,
     &     AMPI_ADOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

c      cost = 0
      if (mynode.eq.0) then
         do i = 1, NNODES
            cost = cost + acost(i)
         end do
      end if
# endif
# else
      cost = icost
# endif
      
      MPI_master_only cost = cost +
     &     (x(1) - rdrg2_bck)*(x(1) - rdrg2_bck)

      return
      end subroutine cost_fun

      subroutine cost_fun_step(icost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"
      real icost
      integer ntrds,trd,my_first,my_last,tile,range,k

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
      integer tile, k
      real icost
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
      integer Istr,Iend,Jstr,Jend,i,j,k
      real delta_z, icost
#include "compute_extended_bounds.h"
      
      do j=Jstr,Jend
         do i=Istr,Iend
!     print*, '>>>',i,j,sim_iicroot+k, zeta(i,j,knew),
!     &              obsz(i,j,sim_iicroot+k)
            delta_z = (zeta(i,j,knew)-obs(i,j,sim_iicroot+k))
            icost = icost + delta_z*delta_z

         end do
      end do

      end subroutine cost_fun_step_2d_tile
      

      subroutine set_state(x, h0, rdrg2_bck)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"      

      real rdrg2_bck, h0(GLOBAL_2D_ARRAY)
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
         call set_state_2d(tile, x, h0, rdrg2_bck)
      end do

      end subroutine set_state


      subroutine set_state_2d(tile, x, h0, rdrg2_bck)

      implicit none
# include "param.h"
# include "adparam.h"
      integer tile
      real rdrg2_bck, h0(GLOBAL_2D_ARRAY)
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, x, h0, rdrg2_bck)

      end subroutine set_state_2d


      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend,x,h0,rdrg2_bck)
      implicit none

#include "param.h"
#include "grid.h"
#include "ocean2d.h" 
#include "scalars.h"m
#include "mpi_roms.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend, i,j
      real rdrg2_bck, h0(GLOBAL_2D_ARRAY)
#include "compute_extended_bounds.h"
      
      write(*,*) 'set_state'
      write(*,*) 'krhs=',krhs

      do j=JstrR,JendR
         do i=IstrR,IendR
            h0(i,j) = h(i,j)
            h(i,j) = real(x(1+i+(IstrR-IendR+1)*(j-1)))
         end do
      end do

      rdrg2_bck = rdrg2
      rdrg2 = real(x(1))

c      write(777, 900) rdrg2
c      write(778, 900) h

      write(*, *) 'rdrg2=',rdrg2
      write(*, *) 'IstrR=',IstrR
      write(*, *) 'JstrR=',JstrR
      write(*, *) 'h(IstrR,JstrR)=',h(IstrR,JstrR)

 900  format(1PE16.9)

      return
      end subroutine set_state_2d_tile
