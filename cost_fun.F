#ifndef MASKING
#define rmask(X,Y) 1
#endif
      
      subroutine ad_step()
      implicit none
# include "adinc.h"
# include "adparam.h"
      
      integer k

C$AD BINOMIAL-CKP ad_ns+1 95 1
      do k=1,ad_ns
         call step()
      end do

      end subroutine ad_step
      
      subroutine cost_fun(ad_x, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost,icost
      double precision sigma2
      integer ta,k,xk

      integer i,j
      double precision delta_z

      call set_state(ad_x)

c      sigma2 = ad_cost/ad_cost_counter
      icost=0
      cost=0
      do ta=1,ad_nt
         call ad_step()
         call cost_fun_step(ad_x,icost,ta,3)
         write(*,*) '==icost3==>',icost
         ad_irms = cost
      
         cost=cost+icost
      end do

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
      double precision delta_z, delta_zob, icost

c      do j=Jstr,Jend
c         do i=Istr,Iend
c            if (ad_obs(i,j,1).ne.ad_spval) then
c               delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+2))
c               delta_zob = (zob(i,j) - zob_bck(i,j))
c               icost = icost + delta_z*delta_z + 1.0D-2*
c     &              delta_zob*delta_zob/
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
c         delta_zob = (zob(i,j) - zob_bck(i,j))
c         end do
      if (mode.eq.2) then
         do j=Jstr,Jend
            do i=Istr,Iend
               delta_zob = (zob(i,j) - zob_bck(i,j))
               icost = icost + delta_zob*delta_zob
            end do
         end do
      else if (mode.eq.3) then
         k=0
         do j=Jstr,Jend
            do i=Istr,Iend
               if (rmask(i,j).ne.0) then
                  k=k+1
                  delta_z = zeta(i,j,knew) - ad_obs(i,j,ta+2)
                  ad_dz(k)=delta_z
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
c$$$     &                 zob(i,j),zob_bck(i,j),
c$$$     &                 delta_z,
c$$$     &                 delta_zob
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
            if (rmask(i,j).ne.0) then
               k=k+1
               zob(i,j)= ad_x(k)
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

      write(*, *) 'k, ad_array_node_size', k, ad_array_node_size
      write(*, *) 'min ad_x', minval(ad_x(1:k))
      write(*, *) 'max ad_x', maxval(ad_x(1:k))
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
      zob_bck(:,:) = zob(:,:)
      zobt_bck = zobt

      kstp_bck = kstp
      krhs_bck = krhs
      knew_bck = knew
      iic_bck = iic

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) 'minzob',minval(zob_bck)
      write(*,*) 'maxzob',maxval(zob_bck)
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"

      ubar(:,:,:) = ad_ubar_bck(:,:,:)
      vbar(:,:,:) = ad_vbar_bck(:,:,:)
      zeta(:,:,:) = ad_zeta_bck(:,:,:)
      zob(:,:) = zob_bck(:,:)
      Zobt = Zobt_bck

      kstp = kstp_bck
      krhs = krhs_bck
      knew = knew_bck
      iic = iic_bck

      write(*,*) 'restore+++'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '++++++++++'

      end subroutine restore_croco_state



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
      do j=Jstr,Jend
         do i=Istr,Iend
            if (rmask(i,j).ne.0) then
               k=k+1
               ad_x(k)=zob(i,j)
            end if
         end do
      end do
      ad_array_node_size=k
      
      write(*,*) 'node,array_node_size, add_array_size/nnodes:',mynode,
     &     ad_array_node_size,
     &     ad_array_size/nnodes

      write(*,*) 'minad_x,maxad_x',minval(ad_x),maxval(ad_x)
      write(*,*) 'minh,maxh',minval(h),maxval(h)
      write(*,*) 'minobs,maxobs',minval(ad_obs),maxval(ad_obs)
      
      end subroutine init_local_arrays_2d_tile


