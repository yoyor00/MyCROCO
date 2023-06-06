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

      integer Istr,Iend,Jstr,Jend,i,j,k,ik
      double precision ridge_width, ridge_height,depth,cff1,cff,rnd

# include "compute_extended_bounds.h"

#ifdef INTERNAL      
c set a lower depth
      depth=2100.
      ridge_width = 30.0D3
      ridge_height = 1600

      cff1=1./(ridge_width*ridge_width)
      do j=JstrR,JendR
        ik=0
        k=0
        do i=IstrR,IendR
          if ((xr(i,j).gt.ridge_width).or.
     &       (xr(i,j).lt.(-ridge_width))) then
            cff=0.
          else
            cff=1-(xr(i,j)*xr(i,j)*cff1);
          endif
          h(i,j)=depth-ridge_height*cff*cff;
          ad_tab0(i,j) = 0.
#   ifdef INTERNALSHELF
          if (xr(i,j).ge.0.) then
            h(i,j)=depth-ridge_height
          endif
#   endif
       enddo
      enddo

c      ad_x(1) =  2.9525087633157719E-007
c      ad_x(2) = -1.3858541365846966E-320
c      ad_x(3) = 5.1097759790190960E-008
      ad_x(1) = 0
c      ad_x(2) = 1D-10 * ridge_width
c      ad_x(3) = 1D-10 * ridge_height
      ad_array_real_node_size=1


      write(*,*) 'node,array_node_size, add_array_size/nnodes:',mynode,
     &     ad_array_real_node_size,
     &     ad_array_size/nnodes

      write(*,*) 'minad_x,maxad_x',minval(ad_x),maxval(ad_x)
      write(*,*) 'minh,maxh',minval(h),maxval(h)
      write(*,*) 'minobs,maxobs',minval(ad_obs),maxval(ad_obs)

#elif defined BASIN
#endif      

      end subroutine init_local_arrays_2d_tile
