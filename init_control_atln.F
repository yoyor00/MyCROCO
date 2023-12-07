      subroutine init_control(ad_x,ad_g)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer ierr
      integer k
      
      do k=1,ad_array_size/nnodes
         ad_g(k)=0.
         ad_x(k)=0.
      end do
      
      call init_local_arrays(ad_x)
      
c      MPI_master_only write(*,*) 'B0##>',ad_g_f

      end subroutine init_control

      subroutine init_local_arrays(ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k
      integer ierr
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
         call init_local_arrays_2d(tile,ad_x)
      end do

      end subroutine init_local_arrays

      subroutine init_local_arrays_2d(tile,ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer k
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call init_local_arrays_2d_tile(Istr, Iend, Jstr, Jend,ad_x)

      end subroutine init_local_arrays_2d

      subroutine init_local_arrays_2d_tile(Istr,Iend,Jstr,Jend,ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,i,j,k
      
      k=0
      ad_x(1) = zobt
c      do j=Jstr,Jend
c         do i=Istr,Iend
c            if (rmask(i,j).ne.0) then
c               k=k+1
c               ad_x(k)=zobt
c            end if
c         end do
c      end do
c      ad_array_node_size=k
      ad_array_real_node_size=1
      
      write(*,*) 'node,array_node_size, add_array_size/nnodes:',mynode,
     &     ad_array_node_size,
     &     ad_array_size/nnodes

c      write(*,*) 'minad_x,maxad_x',minval(ad_x),maxval(ad_x)
c      write(*,*) 'minh,maxh',minval(h),maxval(h)
c      write(*,*) 'minobs,maxobs',minval(ad_obs),maxval(ad_obs)
      
      end subroutine init_local_arrays_2d_tile
      
