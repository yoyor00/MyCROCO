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

      call init_scalars(ierr)

      do tile=0,NSUB_X*NSUB_E-1
        call start_timers()
        call init_arrays (tile)
      enddo
      do tile=my_first,my_last
         call set_state_2d(tile, ad_x)
      end do

c      call check_tab2d(h(:,:),'h initialisation #1','r')

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid1 (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid2 (tile)
      enddo

      call set_scoord
      call set_weights

      iif_bak = iif
      iif = -1
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo

      iif = iif_bak

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call grid_stiffness (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_initial (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
c      do tile=0,NSUB_X*NSUB_E-1
c        call set_HUV (tile)
# ifdef RESET_RHO0
c        call reset_rho0 (tile)
# endif
c      enddo
CR      write(*,*)  ' -4' MYID
C$OMP PARALLEL DO PRIVATE(tile)
c      do tile=0,NSUB_X*NSUB_E-1
c         call omega (tile)
c         call rho_eos (tile)
c      enddo
CR      write(*,*)  ' -3' MYID


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

      subroutine set_state_2d_2_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cff
      integer Istr,Iend,Jstr,Jend,i,j,k,itrc

      real t_bckp(GLOBAL_2D_ARRAY,N,3,NT)
#include "compute_auxiliary_bounds.h"
!
#if defined EW_PERIODIC && !defined MPI
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
#endif
#if defined NS_PERIODIC && !defined MPI
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#endif

c     cf ana_initial
      cff=((2.e-3)**2)/g
      do k=1,N
         do j=JR_RANGE
            do i=IR_RANGE
               t_bckp(i,j,k,1,itemp)=t(i,j,k,1,itemp)
               t_bckp(i,j,k,2,itemp)=t(i,j,k,2,itemp)
               t(i,j,k,1,itemp)=rho0*(1-cff*z_r(i,j,k)) - 1000.
               t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
               if (t_bckp(i,j,k,1,itemp) .ne. t(i,j,k,1,itemp)) then
                  write(*,*) '!!!!!!',t_bckp(i,j,k,1,itemp),
     &                 '!=',t(i,j,k,1,itemp)
               end if
               if (t_bckp(i,j,k,2,itemp) .ne. t(i,j,k,2,itemp)) then
                  write(*,*) '!!!!!!',t_bckp(i,j,k,2,itemp),
     &                 '!=',t(i,j,k,2,itemp)
               end if
            enddo
         enddo
      enddo

      do itrc=1,NT
         if (.not.got_tini(itrc)) then
            call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &           t(START_2D_ARRAY,1,1,itrc))
            call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &           t(START_2D_ARRAY,1,2,itrc))
         end if
      end do

      end subroutine set_state_2d_2_tile

      subroutine set_state_2d_1_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k,ik,cpoint
      double precision cff1,cff,depth,ridge_width,ridge_height,d1,d2,d3
      character*6 counter_str
#include "compute_extended_bounds.h"

#ifdef INTERNAL
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
#elif defined BASIN
#endif
      end subroutine set_state_2d_1_tile
