
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
c               delta_zob = (zob(i,j) - zob_bck(i,j))
c               icost = icost + delta_zob*delta_zob
            end do
         end do
      else if (mode.eq.3) then
         do j=Jstr,Jend
            do i=Istr,Iend
               if (ad_obs(i,j,1).ne.spval) then
                  delta_z = zeta(i,j,knew) - ad_obs(i,j,ta+2)
                  icost = icost + delta_z*delta_z
c$$$                  write(*,*) 'AAAAA',mynode,i,j,
c$$$     &                 zeta(i,j,knew),ad_obs(i,j,ta),
c$$$     &                 ad_obs(i,j,ta+1), ad_obs(i,j,ta+2)
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



      subroutine cost_fun_full_state(ad_x, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost

      call cost_fun(ad_x, cost)

      end subroutine cost_fun_full_state
