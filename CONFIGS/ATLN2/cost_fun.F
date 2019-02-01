
      subroutine cost_fun(ad_x_f, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost
      double precision xcost
      integer ta

      call set_state(ad_x_f)

      cost = 0

C$AD BINOMIAL-CKP ad_nt 40 1
      do ta=1,ad_nt-1
         call step()
         call cost_fun_step(ad_x_f,cost,ta)
      end do

      write (*,*) 'Zob=',Zobt

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
         delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+1))
         delta_zob = (zob(i,j) - zob_bck(i,j))
c         write(*,*) 'AAAAA',mynode,i,j,knew,sim_iicroot,k,sim_iicroot+k,
c     &        zeta(i,j,knew),ad_obs(i,j,ta+1),zob(i,j),zob_bck(i,j)
c         write(*,*) 'BBBBB',mynode,delta_z,delta_zob
         icost = icost + delta_z*delta_z + 1.D-5*delta_zob*delta_zob
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

#include "compute_extended_bounds.h"
      do j=Jstr,Jend
         do i=Istr,Iend
c            k=ad_colloc(i,j)
c            if (k.lt.1.or.k.gt.ad_array_size) then
c               write (*,*) 'AD BIG PROBLEM k=',k
c               stop
c            endif
            zob(i,j)=0
            do k=1,ad_array_size
               lzob=control_to_param(k,ad_x_f)
c               write(*,*) '###',mynode,i,j,k,ad_x_f(k),W(i,j,k)
               zob(i,j) = zob(i,j) + lzob*W(i,j,k)/SkW(i,j)
            end do
c            zob(i,j)=1.D-7 + 1/(1+exp(ad_x_f(k)))
         end do
      end do

      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif
      write(*, *) 'min zob=',minval(zob)
      write(*, *) 'max zob=',maxval(zob)
      write(*, *) '*******************'

 900  format(1PE16.9)

      return
      end subroutine set_state_2d_tile


      subroutine init_control(ad_x_f,ad_g_f)
      implicit none
# include "adinc.h"
# include "adparam.h"
      
      double precision param_to_control
      integer ierr
      integer k

      do k=1,ad_array_size
         ad_x_f(k)=param_to_control(k)
         ad_g_f(k)=0.
      end do

      call init_collocation_points()

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

      do k=1,ad_array_size
         write(*,*) '###',mynode,ad_latr_f(k),ad_lonr_f(k)
      end do

      ncpoints=0
      do k=1,ad_array_size
         call get_ij(ad_lonr_f(k),ad_latr_f(k),i,j,out)
         if (.not.out.and.(ncpoints.lt.(ad_array_size/nnodes))) then
            ncpoints=ncpoints+1
            ad_i(ncpoints)=i
            ad_j(ncpoints)=j
         end if
      end do

      write(*,*) mynode,'ncpoints=',ncpoints
      do k=1,ncpoints
         write(*,*) 'cpoints',mynode,ad_i(k),ad_j(k)
      end do
      
!     sum of ncpoints should be ad_array_size      
      
      do j=Jstr,Jend
         do i=Istr,Iend
            SkW(i,j)=0.
            rk=max(0.01*sqrt(g*h(i,j))*TM2,10000.)
            do k=1,ad_array_size
               dist=distance(latr(i,j),lonr(i,j),
     &              ad_latr_f(k),ad_lonr_f(k))
               W(i,j,k)=exp(-dist/rk)
               SkW(i,j)=SkW(i,j)+W(i,j,k)
            end do
         end do
      end do

      do j=Jstr,Jend
         do i=Istr,Iend
            mindist=1.D230
            do k=1,ad_array_size
               dist=distance(latr(i,j),lonr(i,j),
     &              ad_latr_f(k),ad_lonr_f(k))
               if (dist.lt.mindist) then
                  ad_colloc(i,j)=k
               endif
            enddo
         enddo
      enddo

      end subroutine init_collocation_points_2d_tile


      double precision function control_to_param(k,ad_x_f)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k
      double precision y, a, b
      parameter (a = 1D-6)
      parameter (b = 1D-2)

      y = ad_x_f(k)

c     R -> (a,b)
      control_to_param = 1.D-7 + 1/(1+exp(y))
c      if (y.le.0) then 
c         control_to_param = (a*y + b*y
c     &        + sqrt(a*a*y*y - 2*a*b*y*y + b*b*y*y + 4) - 2)/(2*y)
c      else
c         control_to_param = (a*y + b*y
c     &        - sqrt(a*a*y*y - 2*a*b*y*y + b*b*y*y + 4) - 2)/(2*y)
c      end if

c      if (control_to_param.ge.b) then
c         write(*,*) '>>>>>', control_to_param,'>=',b
c         stop
c      end if
c      if (control_to_param.le.a) then
c         write(*,*) '>>>>>', control_to_param, '=<',a
c         stop
         
c      end if
      
      end function control_to_param

      double precision function param_to_control(k)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k
      double precision x, a, b
      parameter (a = 1D-6)
      parameter (b = 1D-2)
      
      x = Zobt

c     (a,b) -> R      
      param_to_control = log(1./x  - 1.)
c      param_to_control = 1/(a - x) + 1/(b - x)
      
      end function param_to_control
