      subroutine ad_step(nstep)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k,iifroot,nstep

      write(*,*) 'start ad_step'
C$AD BINOMIAL-CKP nstep+1 72 1
      do k=1,nstep
         do iifroot = 0,nfast+2
            call step()
         end do
c     write (*,*) 'k,iifroot,kstp,knew:',k,iifroot,kstp,knew
         ad_step_counter = ad_step_counter + 1
      end do
      write(*,*) 'end ad_step'

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
      ad_irms=0

c      next_kstp=kstp
c      time_start=time
c      iif = -1
c      nbstep3d = 0
c      iic = ntstart

      if (ad_main_st.gt.0) then
         call ad_step(ad_main_st)
      end if

      do ta=1,ad_nt
         call ad_step(ad_ns)
         call cost_fun_step(ad_x,icost,ad_ast-1+ta,3)
         write(*,*) '==icost3==>',icost
c         call cost_fun_step(ad_x,icost,ta*ad_ns,2)
      end do

      cost = icost / ad_nt

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
      double precision delta_z, delta_h, dpena, icost

      write(*,*) '===cost_fun_step==='
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
                  delta_z = zeta(i,j,knew) - ad_obs(i,j,ta)
c                  dpena =
c     &                 (pena_h(i+1,j)-ad_tab0(i+1,j) +
c     &                  pena_h(i-1,j)-ad_tab0(i-1,j)
c     &                  -2*(pena_h(i,j)-ad_tab0(i,j)))/
c     &                 (om_u(i,j)*om_u(i,j))
                  icost = icost + delta_z*delta_z
c     &                 + 1.D4*dpena*dpena

                  ad_irms = ad_irms + delta_z*delta_z
c                  write(*,*) '=====>',mynode,i,j,ta, ad_step_counter,
c     &                 zeta(i,j,knew),
c     &                 ad_obs(i,j,ta),ad_obs(i,j,ta+1),ad_obs(i,j,ta+2),
c     &                 ad_obs(i,j,ta+3)


c$$$                  write(*,*) 'AAAAA',ad_cost_counter,time,
c$$$     &                 ad_obs_time(ta),zeta(i,j,knew),
c$$$     &                 ad_obs(i,j,ta+2),
c$$$     &                 ad_obs(i,j,ta+1),
c$$$     &                 ad_obs(i,j,ta)
c$$$                  write(*,*) 'AAAAA',mynode,i,j,
c$$$     &                 knew,ad_sim_iicroot,k,
c$$$     &                 ad_sim_iicroot+k,
c$$$     &                 zeta(i,j,knew),ad_obs(i,j,ta+2),
c$$$     &                 h(i,j),h_bck(i,j),
c$$$     &                 delta_z,
c$$$  &                 delta_h
c                  if (zeta(i,j,knew).gt.0) then
c                     if ((ta.gt.3).and.(ta.lt.ad_ns*ad_nt)) then
c                  write(*,*) 'BBBBB',ad_sim_iicroot,ta,mynode,i,j,
c     &                 zeta(i,j,knew),
c     &                 ad_obs(i,j,ta+2),
c     &                 ad_obs(i,j,ta+1),
c     &                 ad_obs(i,j,ta),
c     &                 ad_obs(i,j,ta-1),
c     &                 ad_obs(i,j,ta-2)
c                     end if
c                  end if
                  ad_cost_counter = ad_cost_counter + 1
               end if
            end do
         end do
         ad_irms = ad_irms / k
      end if

      write(*,*) '==================='

      end subroutine cost_fun_step_2d_tile

      subroutine cost_fun_full_state(ad_x, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost

      call cost_fun(ad_x, cost)

      end subroutine cost_fun_full_state
