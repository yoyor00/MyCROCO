      subroutine set_state(ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer ntrds,trd,my_first,my_last,tile,subs,range,ierr
      integer iif_bak

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

c$$$      call init_scalars(ierr)
c$$$
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call start_timers()
c$$$        call init_arrays (tile)
c$$$      enddo
c$$$
c$$$c      call check_tab2d(h(:,:),'h initialisation #1','r')
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call setup_grid1 (tile)
c$$$      enddo
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call setup_grid2 (tile)
c$$$      enddo
c$$$
c$$$      call set_scoord
c$$$      call set_weights
c$$$
c$$$      iif_bak = iif
c$$$      iif = -1
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call set_depth (tile)
c$$$      enddo
c$$$
c$$$      iif = iif_bak
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call grid_stiffness (tile)
c$$$      enddo
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$c#ifdef INTERNAL
c$$$c#endif
      do tile=0,NSUB_X*NSUB_E-1
         call set_state_2d(tile, ad_x)
      enddo
c      do tile=0,NSUB_X*NSUB_E-1
c        call ana_initial (tile)
c      enddo

      
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call set_HUV (tile)
c$$$# ifdef RESET_RHO0
c$$$c        call reset_rho0 (tile)
c$$$# endif
c$$$      enddo
c$$$CR      write(*,*)  ' -4' MYID
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$         call omega (tile)
c$$$         call rho_eos (tile)
c$$$      enddo
c$$$CR      write(*,*)  ' -3' MYID


      end subroutine set_state

      subroutine set_state_2d(tile, ad_x)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile, iif_bak
# include "compute_tile_bounds.h"

      call set_state_2d_1_tile(Istr, Iend, Jstr, Jend, ad_x)

c      iif_bak=iif
c     iif=-1
c      call set_scoord
c      call set_weights

c      call set_depth_bathy_modified(tile)
c      iif=iif_bak

c      call grid_stiffness(tile)

c      call set_state_2d_2_tile(Istr, Iend, Jstr, Jend, ad_x)

      end subroutine set_state_2d

      subroutine set_state_2d_1_tile(Istr,Iend,Jstr,Jend, ad_x)
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
            if (ad_obs(i,j,1).ne.spval) then
               k=k+1
               zob(i,j)= ad_x(1)
            end if
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
            write(99, *) i,j,zob(i,j)
         end do
      end do
      close(99)


      return
 900  format(1PE16.9)
 910  format(I2.2,'-',I3.3)

      end subroutine set_state_2d_1_tile
