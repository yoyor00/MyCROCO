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

      integer Istr,Iend,Jstr,Jend, i,j,k,ik,cpoint,itrc
      double precision cff1,cff2,cff,depth,ridge_width,ridge_height,d1,
     &     d2,d3
      character*6 counter_str
c #include "compute_extended_bounds.h"
#include "compute_auxiliary_bounds.h"

      d1=1D8/(200*200)
      d2=1D8/(4D3*4D3)
      d3=1D8/(200*200)

      depth = 2100         -  ad_x(1)
      ridge_width  = 30.0D3
      ridge_height = 1600
c      ridge_width = 28.0D3 +  4D3*(1- exp(-d2*(ad_x(2)*ad_x(2))))
c      ridge_height = 1600  +  200*(1- exp(-d3*(ad_x(3)*ad_x(3))))

      
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
          pena_h(i,j)=depth-ridge_height*cff*cff;
#   ifdef INTERNALSHELF
          if (xr(i,j).ge.0.) then
            pena_h(i,j)=depth-ridge_height
          endif
#   endif
       enddo
      enddo

      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif

      write(*, *) 'k, ad_array_real_node_size', k,
     &     ad_array_real_node_size
      write(*, *) 'ad_x', ad_x(:)
      write(*, *) 'ad_counter=',ad_counter
      write(*, *) 'depth =', depth
      write(*, *) 'width =', ridge_width
      write(*, *) 'height=', ridge_height
      write(*, *) '*******************'


      write(*,*) 'Istr,Iend',Istr,Iend
      write(*,*) 'Jstr,Jend',Jstr,Jend

      write(counter_str, 910) mynode,ad_counter
      open(unit=99,file='pena_h.'//counter_str)
      write(99, *) Iend-Istr+1,Jend-Jstr+1,0
      do j=Jstr,Jend
         do i=Istr,Iend
            write(99, *) i,j,h(i,j),pena_h(i,j)
         end do
      end do
      close(99)


      return
 900  format(1PE16.9)
 910  format(I2.2,'-',I3.3)

      end subroutine set_state_2d_1_tile
