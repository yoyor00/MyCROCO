
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
      call cost_fun_step(ad_x_f,cost,ta)

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

      integer Istr,Iend,Jstr,Jend,i,j,ix,jx
      integer ta, k
! unused out of domain flag (1 proc only)      
      logical out
      external get_ij
      double precision delta, xlat,xlon,xxl,xxw,icost
      double precision distance,xx,xxmin,mmin

      do j=Jstr,Jend
         do i=Istr,Iend
            xxmin = 1.D38
            if (h(i,j).gt.10.and.ad_obs(i,j,1).ne.spval) then
               do k=1,ad_array_size/2
                  xlat=ad_x_f(2*k-1)
                  xlon=ad_x_f(2*k)
                  call get_ij(xlon,xlat,ix,jx,out)
                  xxw=sqrt(g*(min(h(i,j),h(ix,jx))))*TM2*3600
                  xxl=distance(lonr(i,j),latr(i,j),xlon,xlat)
                  xx=xxl/(xxw*xxw)
                  xx=xx*xx
                  xxmin=min(xxmin,xx)
c     write(*,*) '#####>',i,j,xlat,xlon,latr(i,j),lonr(i,j),k,
c     &              icost
               end do
               icost = icost + xxmin
            end if
         end do
      end do

      end subroutine cost_fun_step_2d_tile

      function mmin(x,y) result (r)
      implicit none
      double precision x,y,r
      r = 0.5*(x+y-abs(x-y))
      end function mmin

      
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
      subroutine set_state_2d(tile,ad_x_f)

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
      write(*, *) '**** set_state ****'
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


      call init_collocation_points()
      
      do k=1,ad_array_size/2
         ad_x_f(2*k-1)=ad_latr_f(k)
         ad_x_f(2*k)  =ad_lonr_f(k)
         ad_g_f(2*k-1)=0.
         ad_g_f(2*k)  =0.
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
      kstp_bck = kstp
      krhs_bck = krhs
      knew_bck = knew
      iic_bck = iic

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"

      ubar(:,:,:) = ubar_bck(:,:,:)
      vbar(:,:,:) = vbar_bck(:,:,:)
      zeta(:,:,:) = zeta_bck(:,:,:)
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
      double precision mindist,dist,rk,xr,yr,distance

      k=1
      do while (k.le.(ad_array_size/2))
         call random_number(xr)
         call random_number(yr)
         i=Istr+floor((Iend+1-Istr)*xr)
         j=Jstr+floor((Jend+1-Jstr)*yr)
         if (h(i,j).gt.10.and.ad_obs(i,j,1).ne.spval) then
            write(*,*) 'new collocation point:',mynode,i,j,latr(i,j),
     &           lonr(i,j),h(i,j),
     &           ad_obs(i,j,1)
            ad_latr(k)=latr(i,j)
            ad_lonr(k)=lonr(i,j)
            ad_h(k)=h(i,j)
            write(778,*) lonr(i,j),latr(i,j)
            k=k+1
         endif
            
      enddo
            
#ifdef AMPI
      call AMPI_Barrier(MPI_COMM_WORLD,ierr)
      call AMPI_Gather(ad_latr,ad_array_size/(2*nnodes),
     &     AMPI_ADOUBLE_PRECISION,ad_latr_f,ad_array_size/(2*nnodes),
     &     AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call AMPI_Gather(ad_lonr,ad_array_size/(2*nnodes),
     &     AMPI_ADOUBLE_PRECISION,ad_lonr_f,ad_array_size/(2*nnodes),
     &     AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call AMPI_Gather(ad_h,ad_array_size/(2*nnodes),
     &     AMPI_ADOUBLE_PRECISION,ad_h_f,ad_array_size/(2*nnodes),
     &     AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call AMPI_Bcast(ad_latr_f,ad_array_size/2,AMPI_ADOUBLE_PRECISION,
     &     0,MPI_COMM_WORLD,ierr)
      call AMPI_Bcast(ad_lonr_f,ad_array_size/2,AMPI_ADOUBLE_PRECISION,
     &     0,MPI_COMM_WORLD,ierr)
      call AMPI_Bcast(ad_h_f,ad_array_size/2,AMPI_ADOUBLE_PRECISION,
     &     0,MPI_COMM_WORLD,ierr)
#else
      ad_latr_f = ad_latr
      ad_lonr_f = ad_lonr
#endif
      do k=1,(ad_array_size/2)
         write(*,*) '###',ad_latr_f(k),ad_lonr_f(k),ad_h_f(k)
      end do

      do j=Jstr,Jend
         do i=Istr,Iend
            SkW(i,j)=0.
            rk=max(0.01*sqrt(g*h(i,j))*TM2,10000.)
            do k=1,ad_array_size/2
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
            do k=1,ad_array_size/2
               dist=distance(lonr(i,j),latr(i,j),
     &              ad_lonr_f(k),ad_latr_f(k))
               if (dist.lt.mindist) then
                  ad_colloc(i,j)=k
               endif
            enddo
         enddo
      enddo

      end subroutine init_collocation_points_2d_tile

