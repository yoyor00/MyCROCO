
      subroutine cost_fun(ad_x_f, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost
      double precision xcost
      double precision sigma2
      integer ta,k,xk

      integer i,j
      double precision delta_z

      call set_state(ad_x_f)

c      sigma2 = ad_cost/ad_cost_counter
      cost = 0

C$AD BINOMIAL-CKP ad_nt+1 40 1
      do ta=1,ad_nt
         do xk=1,ad_ns
            call step()
         end do
         call cost_fun_step(ad_x_f,cost,ta)
      end do

c      if (sigma2.gt.0) then
c         cost = cost / sigma2
c      end if
      
      ad_cost = cost
      write (*,*) 'cost=',cost

      return
      end subroutine cost_fun

C$AD NOCHECKPOINT
      subroutine cost_fun_step(ad_x_f,icost,ta)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision icost
      integer ta
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
         call cost_fun_step_2d(tile,ad_x_f,icost,ta)
      end do

      end subroutine cost_fun_step

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d(tile,ad_x_f,icost,ta)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer ta
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call cost_fun_step_2d_tile(Istr,Iend,Jstr,Jend,ad_x_f,icost,ta)

      end subroutine cost_fun_step_2d

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d_tile(Istr,Iend,Jstr,Jend,
     &     ad_x_f,icost,ta)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,i,j
      integer ta, k
      double precision delta_z, delta_zob, icost

      do k=1,ncpoints
         i=ad_i(k)
         j=ad_j(k)
         delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+2))
         delta_zob = (zob(i,j) - zob_bck(i,j))
c         write(*,*) 'AAAAA',mynode,i,j,knew,sim_iicroot,k,sim_iicroot+k,
c     &        zeta(i,j,knew),ad_obs(i,j,ta+1),zob(i,j),zob_bck(i,j)
c         write(*,*) 'BBBBB',mynode,delta_z,zeta(i,j,knew),
c     &        ad_obs(i,j,ta+1),
c     &        ad_obs(i,j,ta),ad_obs(i,j,ta+2),ad_obs(i,j,ta-1),
c     &        ad_obs(i,j,ta-2)
         icost = icost + delta_z*delta_z + 0.0001*delta_zob*delta_zob/
     &        (h(i,j)*h(i,j))
         ad_cost_counter = ad_cost_counter + 1
      end do

      end subroutine cost_fun_step_2d_tile

C$AD NOCHECKPOINT
      subroutine set_state(ad_x_f)
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
         call set_state_2d(tile, ad_x_f)
      end do

      end subroutine set_state


C$AD NOCHECKPOINT
      subroutine set_state_2d(tile, ad_x_f)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, ad_x_f)

      end subroutine set_state_2d

C$AD NOCHECKPOINT
      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend, ad_x_f)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k
      double precision lzob
      double precision control_to_param
      character*6 counter_str

      double precision zobmin,zobmax
      
#include "compute_extended_bounds.h"
      zobmin=1000.
      zobmax=0.
      
      do j=Jstr,Jend
         do i=Istr,Iend
c            k=ad_colloc(i,j)
c            if (k.lt.1.or.k.gt.ad_array_size) then
c               write (*,*) 'AD BIG PROBLEM k=',k
c               stop
c     endif
            if (ad_obs(i,j,1).ne.spval) then
               zob(i,j)=0.0
               do k=1,ad_array_size
                  lzob=control_to_param(ad_x_f(k))
                  zob(i,j) = zob(i,j) + lzob*W(i,j,k)/SkW(i,j)
c                  write(*,*) '%%%%',mynode,k,ad_x_f(k),
c     &                 lzob,zob(i,j)
               end do
               if (zobmin.gt.zob(i,j)) then
                  zobmin=zob(i,j)
               end if
               if (zobmax.lt.zob(i,j)) then
                  zobmax=zob(i,j)
               end if
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
            write(99, *) lonr(i,j),latr(i,j),zob(i,j)
         end do
      end do
      close(99)


      return
 900  format(1PE16.9)
 910  format(I2.2,'-',I3.3)
      end subroutine set_state_2d_tile


      subroutine init_control(ad_x_f,ad_g_f)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision param_to_control
      integer ierr
      integer k
      
      call init_collocation_points()
      
      do k=1,ad_array_size
         ad_x_f(k)=param_to_control(5.0D-3)
         ad_g_f(k)=0.
      end do



c      MPI_master_only write(*,*) 'A0##>',ad_x_f
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



      subroutine init_collocation_points()
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

      if (mynode.eq.0) then
         open(unit=99,file="cpoints.dat")
         do k=1,ad_array_size
            read(99,*) ad_lonr_f(k), ad_latr_f(k)
         end do
         close(unit=99)
      end if
      call AMPI_Barrier(MPI_COMM_WORLD,ierr)
      call AMPI_Bcast(ad_lonr_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &     0,MPI_COMM_WORLD,ierr)
      call AMPI_Bcast(ad_latr_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &     0,MPI_COMM_WORLD,ierr)

      do tile=my_first,my_last
         call init_collocation_points_2d(tile)
      end do

      end subroutine init_collocation_points

      subroutine init_collocation_points_2d(tile)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer k
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call init_collocation_points_2d_tile(Istr, Iend, Jstr, Jend)

      end subroutine init_collocation_points_2d

      subroutine init_collocation_points_2d_tile(Istr,Iend,Jstr,Jend)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,i,j,k,ierr
      logical out
      double precision mindist,dist,rk,distance,xr,yr
      character*2 node_str

      do k=1,ad_array_size
         write(*,*) '###',mynode,ad_latr_f(k),ad_lonr_f(k)
      end do

      ncpoints=0
      do k=1,ad_array_size
         call get_ij(ad_lonr_f(k),ad_latr_f(k),i,j,out)
         if (.not.out.and.(ncpoints.lt.(ad_array_size))) then
            ncpoints=ncpoints+1
            ad_i(ncpoints)=i
            ad_j(ncpoints)=j
         end if
      end do

      write(node_str, 930) mynode
      open(unit=99, file='cpoints-'//node_str//'.dat')
      do k=1,ncpoints
         write(99,*) lonr(ad_i(k),ad_j(k)), latr(ad_i(k),ad_j(k))
      end do
      close(99)
 930  format(I2.2)
!     sum of ncpoints should be ad_array_size

      do j=Jstr,Jend
         do i=Istr,Iend
            SkW(i,j)=0.
            rk=max(0.01*sqrt(g*h(i,j))*TM2,10000.)
            do k=1,ad_array_size
               dist=distance(lonr(i,j),latr(i,j),
     &              ad_lonr_f(k),ad_latr_f(k))
               W(i,j,k)=exp(-dist/rk)
               SkW(i,j)=SkW(i,j)+W(i,j,k)
            end do
         end do
      end do

      do j=Jstr,Jend
         do i=Istr,Iend
            mindist=1.D230
            do k=1,ad_array_size
               dist=distance(lonr(i,j),latr(i,j),
     &              ad_lonr_f(k),ad_latr_f(k))
               if (dist.lt.mindist) then
                  ad_colloc(i,j)=k
               endif
            enddo
         enddo
      enddo

      end subroutine init_collocation_points_2d_tile


      double precision function control_to_param(y)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision y

c     R -> (a,b)
      control_to_param = 1/(1+exp(y))

      end function control_to_param

      double precision function param_to_control(x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision x

c     (a,b) -> R
      param_to_control = log(1./(x)  - 1.)
c      param_to_control = 1/(a - x) + 1/(b - x)

      end function param_to_control
