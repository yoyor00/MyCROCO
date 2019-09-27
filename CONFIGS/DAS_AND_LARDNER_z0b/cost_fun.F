
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
      icost = 0

C$AD BINOMIAL-CKP ad_nt*ad_ns+1 80 1
      do ta=1,ad_nt*ad_ns
         call step()
         if (mod(ta,ad_ns).eq.0) then
            call cost_fun_step(ad_x,icost,ta,3)
         end if
      end do

      cost=icost
      icost=0
      call cost_fun_step(ad_x,icost,ad_nt,2)

      write(*,*) '==cost==>',cost,sqrt(icost)
      cost=cost
      
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
c            if (ad_obs(i,j,1).ne.spval) then
c               delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+2))
c               delta_zob = (zob(i,j) - zob_bck(i,j))
c               icost = icost + delta_z*delta_z + 1.0D-2*
c     &              delta_zob*delta_zob/
c     &              (h(i,j)*h(i,j))
c               ad_cost_counter = ad_cost_counter + 1
c            end if
c         end do
c      end do

      if (mode.eq.1) then
         do k=1,ncpoints
            i=ad_i(k)
            j=ad_j(k)
            delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+2))
c         delta_zob = (zob(i,j) - zob_bck(i,j))
         end do
      else if (mode.eq.2) then
         do j=Jstr,Jend
            do i=Istr,Iend
               delta_zob = (zob(i,j) - zob_bck(i,j))
               icost = icost + delta_zob*delta_zob
            end do
         end do
      else if (mode.eq.3) then
         do j=Jstr,Jend
            do i=Istr,Iend
               if (ad_obs(i,j,1).ne.spval) then
                  delta_z = zeta(i,j,knew) - ad_obs(i,j,ta+2)
                  icost = icost + delta_z*delta_z
c$$$                  write(*,*) 'AAAAA',mynode,i,j,knew,sim_iicroot,k,
c$$$     &                 sim_iicroot+k,
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

      end subroutine cost_fun_step_2d_tile

      subroutine set_state(ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"
      double precision rms
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
      double precision control_to_param
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
               
            
c            k=ad_colloc(i,j)
c            if (k.lt.1.or.k.gt.ad_array_size) then
c               write (*,*) 'AD BIG PROBLEM k=',k
c               stop
c     endif
c            if (ad_obs(i,j,1).ne.spval) then
c               zob(i,j)=0.0
c               cpoint=ad_colloc(i,j)
c               zob(i,j)=control_to_param(
c     &              ad_x_f(cpoint))
c               cpoint=0
c               do k=1,npcpoints
c                  if(ad_cpoint_node(k).eq.mynode) then
c                     cpoint=cpoint+1
c                     lzob=control_to_param(
c     &                    ad_x_f(mynode*(ad_array_size/nnodes)+cpoint))
c                     zob(i,j) = zob(i,j) + lzob*W(i,j,cpoint)/SkW(i,j)
c                  endif
c     write(*,*) '%%%%',mynode,k,ad_x_f(k),
c     &                 lzob,zob(i,j)
c               end do
c            end if
c         end do
c      end do

c      do k=1,ncpoints
c         zob(ad_i(k), ad_j(k))=
c     &        control_to_param(ad_x_f(mynode*(ad_array_size/nnodes)+k))
c      end do

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
      end subroutine set_state_2d_tile


      subroutine init_control(ad_x,ad_g)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision param_to_control
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

      ubar_bck(:,:,:) = ubar(:,:,:)
      vbar_bck(:,:,:) = vbar(:,:,:)
      zeta_bck(:,:,:) = zeta(:,:,:)
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

      ubar(:,:,:) = ubar_bck(:,:,:)
      vbar(:,:,:) = vbar_bck(:,:,:)
      zeta(:,:,:) = zeta_bck(:,:,:)
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
      double precision param_to_control
      
      k=0
      do j=Jstr,Jend
         do i=Istr,Iend
            if (ad_obs(i,j,1).ne.spval) then
               k=k+1
               ad_x(k)=zob(i,j)
            end if
         end do
      end do
      ad_array_node_size=k
      
      write(*,*) 'node,array_node_size, add_array_size/nnodes:',mynode,
     &     ad_array_node_size,
     &     ad_array_size/nnodes

      
      end subroutine init_local_arrays_2d_tile


      double precision function control_to_param(y)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision y

      double precision a,b,c,d
      parameter(a=0.0002,b=0.1,c=-1000000,d=1000000)
c     (c,d) -> (a,b)
      control_to_param = (-a*d + a*y + b*c - b*y)/(c - d)
c      control_to_param = 1/(1+exp(y))

      end function control_to_param

      double precision function param_to_control(x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision x

      double precision a,b,c,d
      parameter(a=0.0002,b=0.1,c=-1000000,d=1000000)
      
c     (a,b) -> (c,d)
      param_to_control = c+((d-c)*(x-a)/(b-a))
      
c      param_to_control = log(1./(x)  - 1.)
c      param_to_control = 1/(a - x) + 1/(b - x)

      end function param_to_control
