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
      cost=0
      do ta=1,ad_nt
         icost=0
         call ad_step()
         call cost_fun_step(ad_x,icost,ad_ast-1+ta*ad_ns,3)
         write(*,*) '==icost3==>',icost

         cost=cost+icost
      end do

      ad_irms = cost
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
      integer Istr,Iend,Jstr,Jend,i,j,k,l
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
c     if (mode.eq.2.or.mode.eq.3) then
         do l=1,ad_nobs
            i=ad_obs_i(l)
            j=ad_obs_j(l)
            delta_z = zeta(i,j,knew) - ad_obs(l,ta)
            icost = icost + delta_z*delta_z
            write(*,*) 'AAAAA',mynode,ta,i,j,
     &                 delta_z,
     &                 icost,
     &                 zeta(i,j,knew),
     &                 ad_obs(l,ta),
     &                 ad_obs(l,ta+1),
     &                 ad_obs(l,ta+2)
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
c$$$  &                 ad_obs(i,j,ta-2)
         end do
         ad_cost_counter = ad_cost_counter + 1
c      end if

c      write(*,*) '===cost_fun_step==='
c      write(*,*) 'min delta_z',minval(ad_dz)
c      write(*,*) 'max delta_z',maxval(ad_dz)
c      write(*,*) '==================='
      
      end subroutine cost_fun_step_2d_tile


      subroutine cost_fun_full_state(ad_x, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost

      call cost_fun(ad_x, cost)

      end subroutine cost_fun_full_state
