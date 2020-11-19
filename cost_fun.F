      subroutine ad_step()
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k

C$AD BINOMIAL-CKP ad_ns+1 80 1
      do k=1,ad_ns
         call step()
      end do

      end subroutine ad_step

      subroutine cost_fun(ad_x, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost,icost
      integer ta,k,xk

      integer i,j
      double precision delta_z

      call set_state(ad_x)

      icost=0
      cost=0

      do ta=1,ad_nt
         call ad_step()
         call cost_fun_step(ad_x,icost,ta*ad_ns,3)
         write(*,*) '==icost3==>',icost
c         call cost_fun_step(ad_x,icost,ta*ad_ns,2)
      end do

      cost = icost
      ad_irms = icost

      ad_cost = cost
      ad_ta = ad_nt*ad_ns
      write (*,*) 'cost=',cost

      return
      end subroutine cost_fun

      subroutine cost_fun_step(ad_x,icost,ta,mode)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision icost
      integer ta,mode
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
         call cost_fun_step_2d(tile,ad_x,icost,ta,mode)
      end do

      end subroutine cost_fun_step

      subroutine cost_fun_step_2d(tile,ad_x,icost,ta,mode)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer ta,mode
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call cost_fun_step_2d_tile(Istr,Iend,Jstr,Jend,ad_x,icost,ta,
     &     mode)

      end subroutine cost_fun_step_2d

      subroutine cost_fun_step_2d_tile(Istr,Iend,Jstr,Jend,
     &     ad_x,icost,ta,mode)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer ta,mode
      integer Istr,Iend,Jstr,Jend,i,j,k
      double precision delta_z, delta_h, icost

c      do j=Jstr,Jend
c         do i=Istr,Iend
c            if (ad_obs(i,j,1).ne.ad_spval) then
c               delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+2))
c               delta_h = (h(i,j) - h_bck(i,j))
c               icost = icost + delta_z*delta_z + 1.0D-2*
c     &              delta_h*delta_h/
c     &              (h(i,j)*h(i,j))
c               ad_cost_counter = ad_cost_counter + 1
c            end if
c         end do
c      end do

c      if (mode.eq.1) then
c         do k=1,ncpoints
c            i=ad_i(k)
c            j=ad_j(k)
c            delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+2))
c         delta_h = (h(i,j) - h_bck(i,j))
c         end do
      if (mode.eq.2) then
         do j=Jstr,Jend
            do i=Istr,Iend
               delta_h = (h(i,j) - h_bck(i,j))
               icost = icost + delta_h*delta_h
            end do
         end do
      else if (mode.eq.3) then
         k=0
         do j=Jstr,Jend
            do i=Istr,Iend
               if (rmask(i,j).ne.0) then
                  k=k+1
                  delta_z = zeta(i,j,knew) - ad_obs(i,j,ta+2)
                  ad_dz(k) = delta_z
                  icost = icost + delta_z*delta_z
c                  write(*,*) 'AAAAA',ad_cost_counter,time,
c     &                 ad_obs_time(ta),zeta(i,j,knew),
c     &                 ad_obs(i,j,ta+2),
c     &                 ad_obs(i,j,ta+1),
c     &                 ad_obs(i,j,ta)
c$$$                  write(*,*) 'AAAAA',mynode,i,j,
c$$$     &                 knew,ad_sim_iicroot,k,
c$$$     &                 ad_sim_iicroot+k,
c$$$     &                 zeta(i,j,knew),ad_obs(i,j,ta+2),
c$$$     &                 h(i,j),h_bck(i,j),
c$$$     &                 delta_z,
c$$$     &                 delta_h
c$$$                  write(*,*) 'BBBBB',mynode,zeta(i,j,knew),
c$$$     &                 ad_obs(i,j,ta+3),
c$$$     &                 ad_obs(i,j,ta+2),
c$$$     &                 ad_obs(i,j,ta+1),
c$$$     &                 ad_obs(i,j,ta),
c$$$     &                 ad_obs(i,j,ta-1),
c$$$     &                 ad_obs(i,j,ta-2)
                  ad_cost_counter = ad_cost_counter + 1
               end if
            end do
         end do
      end if

      write(*,*) '===cost_fun_step==='
      write(*,*) 'min delta_z',minval(ad_dz)
      write(*,*) 'max delta_z',maxval(ad_dz)
      write(*,*) '==================='

      end subroutine cost_fun_step_2d_tile

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

      call set_depth(tile)

      call grid_stiffness(tile)

      end subroutine set_state_2d

      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k,cpoint
      double precision cff1,cff,depth,ridge_width,ridge_height
      character*6 counter_str

#include "compute_extended_bounds.h"

      k=ad_array_node_size
      depth=ad_x(3)
      ridge_width=ad_x(2)
      ridge_height=ad_x(1)

      write(*,*) '########',ad_x(1),'...',ad_x(2),'...',ad_x(3)


      cff1=1./(ridge_width*ridge_width)
      do j=JstrR,JendR
        do i=IstrR,IendR
          if ((xr(i,j).gt.ridge_width).or.
     &       (xr(i,j).lt.(-ridge_width))) then
            cff=0.
          else
            cff=1-(xr(i,j)*xr(i,j)*cff1);
          endif
          h(i,j)=depth-ridge_height*cff*cff;
#   ifdef INTERNALSHELF
          if (xr(i,j).ge.0.) then
            h(i,j)=depth-ridge_height
          endif
#   endif
        enddo
      enddo


      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif

      write(*, *) 'k, ad_array_node_size', k, ad_array_node_size
      write(*, *) 'min ad_x', minval(ad_x(1:k))
      write(*, *) 'max ad_x', maxval(ad_x(1:k))
      write(*, *) 'ad_counter=',ad_counter
      write(*, *) '*******************'


      write(*,*) 'Istr,Iend',Istr,Iend
      write(*,*) 'Jstr,Jend',Jstr,Jend

      write(counter_str, 910) mynode,ad_counter
      open(unit=99,file='z0b.'//counter_str)
      write(99, *) Iend-Istr+1,Jend-Jstr+1,0
      do j=Jstr,Jend
         do i=Istr,Iend
#ifdef SPHERICAL
            write(99, *) lonr(i,j),latr(i,j),h(i,j)
#else
            write(99, *) xr(i,j),yr(i,j),h(i,j)
#endif
         end do
      end do
      close(99)


      return
 900  format(1PE16.9)
 910  format(I2.2,'-',I3.3)
      end subroutine set_state_2d_tile


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

      subroutine save_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"

      ad_ubar_bck(:,:,:) = ubar(:,:,:)
      ad_vbar_bck(:,:,:) = vbar(:,:,:)
      ad_zeta_bck(:,:,:) = zeta(:,:,:)
      h_bck(:,:) = h(:,:)

      kstp_bck = kstp
      krhs_bck = krhs
      knew_bck = knew
      iic_bck = iic

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) 'minh',minval(h_bck)
      write(*,*) 'maxh',maxval(h_bck)
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
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

      ubar(:,:,:) = ad_ubar_bck(:,:,:)
      vbar(:,:,:) = ad_vbar_bck(:,:,:)
      zeta(:,:,:) = ad_zeta_bck(:,:,:)
      h(:,:) = h_bck(:,:)

      kstp = kstp_bck
      krhs = krhs_bck
      knew = knew_bck
      iic = iic_bck

      do tile=my_first,my_last
         call restore_croco_state_2d(tile, ad_x)
      end do

      write(*,*) 'restore+++'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '++++++++++'

      end subroutine restore_croco_state

      subroutine restore_croco_state_2d(tile)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_depth(tile)

      call grid_stiffness(tile)

      end subroutine restore_croco_state_2d


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

# include "compute_extended_bounds.h"
      k=0

c     ridge_height
      ad_x(1)=800
c     ridge_width=30.e3
      ad_x(2)=30.e3
c     depth
      ad_x(3)=2000

      ad_array_node_size=3

      write(*,*) 'node,array_node_size, add_array_size/nnodes:',mynode,
     &     ad_array_node_size,
     &     ad_array_size/nnodes

      write(*,*) 'minad_x,maxad_x',minval(ad_x),maxval(ad_x)
      write(*,*) 'minh,maxh',minval(h),maxval(h)
      write(*,*) 'minobs,maxobs',minval(ad_obs),maxval(ad_obs)

      end subroutine init_local_arrays_2d_tile
