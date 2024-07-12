      subroutine set_state(ad_x)
      implicit none
# include "adinc.h"
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

      subroutine set_state_2d(tile, ad_x)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, ad_x)

      end subroutine set_state_2d

      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k,cpoint
      double precision lzob
      character*6 counter_str

      double precision zobmin,zobmax
      
#include "compute_extended_bounds.h"
      zobmin=1.e99
      zobmax=0.
 
      k=0
      do j=Jstr,Jend
         do i=Istr,Iend
c            if (rmask(i,j).ne.0) then
            k=k+1
            zob(i,j)= ad_x(1)
c            write(*,*) 'set zob(i,j)=',  zob(i,j)
c            end if
            if (zobmin.gt.zob(i,j)) then
               zobmin=zob(i,j)
            end if
            if (zobmax.lt.zob(i,j)) then
               zobmax=zob(i,j)
            end if
         end do
      end do
               
      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif

      write(*, *) 'k, ad_array_node_size', k, ad_array_node_size
c      write(*, *) 'min ad_x', minval(ad_x(1:k))
c      write(*, *) 'max ad_x', maxval(ad_x(1:k))
      write(*, *) 'ad_counter=',ad_counter
      write(*, *) 'min zob=',zobmin
      write(*, *) 'max zob=',zobmax
      write(*, *) '*******************'


      write(*,*) 'Istr,Iend',Istr,Iend
      write(*,*) 'Jstr,Jend',Jstr,Jend

      write(counter_str, 910) mynode,ad_counter
      open(unit=99,file='z0b.'//counter_str)
      write(99, *) Iend-Istr+1,Jend-Jstr+1,0
      do j=Jstr,Jend
         do i=Istr,Iend
#ifdef SPHERICAL
            write(99, *) lonr(i,j),latr(i,j),zob(i,j)
#else
            write(99, *) xr(i,j),yr(i,j),zob(i,j)
#endif
         end do
      end do
      close(99)


      return
 900  format(1PE16.9)
 910  format(I2.2,'-',I3.3)
      end subroutine set_state_2d_tile
