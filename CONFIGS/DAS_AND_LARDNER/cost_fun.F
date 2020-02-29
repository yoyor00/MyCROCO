
      subroutine cost_fun(ad_x, cost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      double precision cost
      double precision xcost

      integer k, ierr, i

      rdrg2_bck = rdrg2
      call set_state(ad_x)

      cost = 0

C$AD BINOMIAL-CKP 120 10 1
      do k=1,119
         call step()
         call cost_fun_step(cost, k)
      end do

c     cost = cost +
c     &     (ad_x(1)-control_var_bck)*(ad_x(1)-control_var_bck)

      write (*,*) 'rdrg2=',rdrg2
      write (*,*) 'cost=',cost

      return
      end subroutine cost_fun

C$AD NOCHECKPOINT
      subroutine cost_fun_step(icost, k)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"
      double precision icost
      integer k
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
         call cost_fun_step_2d(tile, icost, k)
      end do

      end subroutine cost_fun_step

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d(tile, icost, k)

      implicit none
# include "param.h"
      integer k
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, k)

      end subroutine cost_fun_step_2d

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, k)
      implicit none
#include "param.h"
#include "grid.h"
#include "scalars.h"
#include "ocean2d.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend,i,j
      integer k
      double precision delta_z, icost

      do j=Jstr,Jend
         do i=Istr,Iend
c            write(*,*) i,j,knew,ad_sim_iicroot+k,
c     &           zeta(i,j,knew),ad_obs(i,j,ad_sim_iicroot+k)

            delta_z = (zeta(i,j,knew)-ad_obs(i,j,ad_sim_iicroot+k))
            icost = icost + delta_z*delta_z
         end do
c         write(*,*) 'icost=',j,icost
      end do

      end subroutine cost_fun_step_2d_tile

C$AD NOCHECKPOINT
      subroutine set_state(ad_x)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
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


C$AD NOCHECKPOINT
      subroutine set_state_2d(tile, ad_x)

      implicit none
# include "param.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, ad_x)

      end subroutine set_state_2d

C$AD NOCHECKPOINT
      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend,ad_x)
      implicit none

#include "param.h"
#include "grid.h"
#include "ocean2d.h"
#include "scalars.h"
#include "mpi_roms.h"
#include "adparam.h"
      integer Istr,Iend,Jstr,Jend, i,j

      integer bi, bj
#include "compute_extended_bounds.h"

      bi = IendR-IstrR+1
      bj = JendR-JstrR+1

c      write(*,*) 'set_state'
c      write(*,*) 'krhs=',krhs

c      write(*,*) bi, bj, ad_array_size

c      do j=JstrR,JendR
c         do i=IstrR,IendR
c            write (*,*) '###',i,j,h(i,j),
c     &           1+(i-IstrR+1)+bi*(j-JstrR),
c     &           ad_x(1+(i-IstrR+1)+bi*(j-JstrR))
c             h(i,j) = ad_x(1+(i-IstrR+1)+bi*(j-JstrR))
c         end do
c      end do

      rdrg2 = real(ad_x(1))

      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif
      write(*, *) 'rdrg2=',rdrg2
      write(*, *) '*******************'

 900  format(1PE16.9)

      return
      end subroutine set_state_2d_tile

      subroutine save_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"

      ad_ubar_bck(:,:,:) = ubar(:,:,:)
      ad_vbar_bck(:,:,:) = vbar(:,:,:)
      ad_zeta_bck(:,:,:) = zeta(:,:,:)
      rdrg2_bck = rdrg2

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

      ubar(:,:,:) = ad_ubar_bck(:,:,:)
      vbar(:,:,:) = ad_vbar_bck(:,:,:)
      zeta(:,:,:) = ad_zeta_bck(:,:,:)
      rdrg2 = rdrg2_bck

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
      
      subroutine init_control(ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer i
      
      do i=1,ad_array_size
         ad_x(i)=rdrg2
      end do
      end subroutine init_control

      double precision function control_to_param(y)

      control_to_param = y

      end function control_to_param

      double precision function param_to_control(x)

      param_to_control = x

      end function param_to_control


