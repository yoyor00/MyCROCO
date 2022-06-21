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
c         write (*,*) 'k,iifroot,kstp,knew:',k,iifroot,kstp,knew
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

      next_kstp=kstp
      time_start=time
      iif = -1
      nbstep3d = 0
      iic = ntstart

      if (ad_main_st.gt.0) then
         call ad_step(ad_main_st)
      end if

      do ta=1,ad_nt
         call ad_step(ad_ns)
         call cost_fun_step(ad_x,icost,ad_ast-1+ta,3)
         write(*,*) '==icost3==>',icost
c         call cost_fun_step(ad_x,icost,ta*ad_ns,2)
      end do

      cost = icost

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
c     end do
c      do j=Jstr,Jend
c         do i=Istr,Iend
c            write(*,*) ta,i,j,zeta(i,j,knew),ad_obs(i,j,ta+1),
c     &           ad_obs(i,j,ta+2)
c         end do
c      end do

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
                  dpena = (pena_h(i+1,j) + pena_h(i-1,j)-2*pena_h(i,j))/
     &                 (om_u(i,j)*om_u(i,j))
                  icost = icost + delta_z*delta_z
c     &                 + 1.D4*dpena*dpena

                  ad_irms = ad_irms + delta_z*delta_z
c                  write(*,*) '=====>',zeta(i,j,knew)
c     &                 ,ad_obs(i,j,ta-1),ad_obs(i,j,ta),ad_obs(i,j,ta+1)
c     &                 ,ad_obs(i,j,ta+2),ad_obs(i,j,ta+3)


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

      subroutine set_state(ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer ntrds,trd,my_first,my_last,tile,subs,range,ierr
      integer iif_bak

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

      call init_scalars(ierr)

      do tile=0,NSUB_X*NSUB_E-1
        call start_timers()
        call init_arrays (tile)
      enddo
      do tile=my_first,my_last
         call set_state_2d(tile, ad_x)
      end do

c      call check_tab2d(h(:,:),'h initialisation #1','r')

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid1 (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid2 (tile)
      enddo

      call set_scoord
      call set_weights

c      iif_bak = iif
c      iif = -1
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo

c      iif = iif_bak

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call grid_stiffness (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_initial (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
c      do tile=0,NSUB_X*NSUB_E-1
c        call set_HUV (tile)
# ifdef RESET_RHO0
c        call reset_rho0 (tile)
# endif
c      enddo
CR      write(*,*)  ' -4' MYID
C$OMP PARALLEL DO PRIVATE(tile)
c      do tile=0,NSUB_X*NSUB_E-1
c         call omega (tile)
c         call rho_eos (tile)
c      enddo
CR      write(*,*)  ' -3' MYID


      end subroutine set_state

      subroutine set_state_2d(tile, ad_x)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile, iif_bak
# include "compute_tile_bounds.h"

      call set_state_2d_1_tile(Istr, Iend, Jstr, Jend, ad_x)

c      iif_bak=iif
c     iif=-1
c      call set_scoord
c      call set_weights

c      call set_depth_bathy_modified(tile)
c      iif=iif_bak

c      call grid_stiffness(tile)

c      call set_state_2d_2_tile(Istr, Iend, Jstr, Jend, ad_x)

      end subroutine set_state_2d

      subroutine set_state_2d_2_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cff
      integer Istr,Iend,Jstr,Jend,i,j,k,itrc

      real t_bckp(GLOBAL_2D_ARRAY,N,3,NT)
#include "compute_auxiliary_bounds.h"
!
#if defined EW_PERIODIC && !defined MPI
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
#endif
#if defined NS_PERIODIC && !defined MPI
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#endif

c     cf ana_initial
      cff=((2.e-3)**2)/g
      do k=1,N
         do j=JR_RANGE
            do i=IR_RANGE
               t_bckp(i,j,k,1,itemp)=t(i,j,k,1,itemp)
               t_bckp(i,j,k,2,itemp)=t(i,j,k,2,itemp)
               t(i,j,k,1,itemp)=rho0*(1-cff*z_r(i,j,k)) - 1000.
               t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
               if (t_bckp(i,j,k,1,itemp) .ne. t(i,j,k,1,itemp)) then
                  write(*,*) '!!!!!!',t_bckp(i,j,k,1,itemp),
     &                 '!=',t(i,j,k,1,itemp)
               end if
               if (t_bckp(i,j,k,2,itemp) .ne. t(i,j,k,2,itemp)) then
                  write(*,*) '!!!!!!',t_bckp(i,j,k,2,itemp),
     &                 '!=',t(i,j,k,2,itemp)
               end if
            enddo
         enddo
      enddo

      do itrc=1,NT
         if (.not.got_tini(itrc)) then
            call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &           t(START_2D_ARRAY,1,1,itrc))
            call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &           t(START_2D_ARRAY,1,2,itrc))
         end if
      end do

      end subroutine set_state_2d_2_tile

      subroutine set_state_2d_1_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k,cpoint
      double precision cff1,cff,depth,ridge_width,ridge_height
      character*6 counter_str

#include "compute_extended_bounds.h"

      do j=JstrR,JendR
         k = 1
         do i=IstrR,IendR
            pena_h(i,j)= ad_x(k)
            k=k+1
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
      open(unit=99,file='pena_h.'//counter_str)
      write(99, *) Iend-Istr+1,Jend-Jstr+1,0
      do j=Jstr,Jend
         do i=Istr,Iend
            write(99, *) i,j,h(i,j),pena_h(i,j)
         end do
      end do
      close(99)


      return
 900  format(1PE16.9)
 910  format(I2.2,'-',I3.3)
      end subroutine set_state_2d_1_tile

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

      CALL PUSHINTEGER4(ad_cost_counter)
      CALL PUSHREAL8ARRAY(srflx, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(btflx, (lm+padd_x+4)*(mm+padd_e+4)*nt)
      CALL PUSHREAL8ARRAY(stflx, (lm+padd_x+4)*(mm+padd_e+4)*nt)
      CALL PUSHREAL8ARRAY(bvstr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(bustr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(svstr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(sustr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(ad_x_rv, n3d*npp)
      CALL PUSHREAL8ARRAY(ad_x_ru, n3d*npp)
      CALL PUSHREAL8ARRAY(a3d, n3d*7*npp)
      CALL PUSHREAL8ARRAY(a2d, n2d*nsa*npp)
      CALL PUSHREAL8ARRAY(dv_avg2, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(du_avg2, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(dv_avg1, (lm+padd_x+4)*(mm+padd_e+4)*5)
      CALL PUSHREAL8ARRAY(du_avg1, (lm+padd_x+4)*(mm+padd_e+4)*5)
      CALL PUSHREAL8ARRAY(zt_avg1, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(rvfrc_bak, (lm+padd_x+4)*(mm+padd_e+4)*2)
      CALL PUSHREAL8ARRAY(rufrc_bak, (lm+padd_x+4)*(mm+padd_e+4)*2)
      CALL PUSHREAL8ARRAY(rvfrc, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(rufrc, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(rhos, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(rhoa, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(akt, (lm+padd_x+4)*(mm+padd_e+4)*(n+1)*2)
      CALL PUSHREAL8ARRAY(akv, (lm+padd_x+4)*(mm+padd_e+4)*(n+1))
      CALL PUSHREAL8ARRAY(rho, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL PUSHREAL8ARRAY(we, (lm+padd_x+4)*(mm+padd_e+4)*(n+1))
      CALL PUSHREAL8ARRAY(z_r, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL PUSHREAL8ARRAY(hz, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL PUSHREAL8ARRAY(hvom, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL PUSHREAL8ARRAY(huon, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL PUSHREAL8ARRAY(z_w, (lm+padd_x+4)*(mm+padd_e+4)*(n+1))
      CALL PUSHREAL8ARRAY(hz_bak, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL PUSHREAL8ARRAY(t, (lm+padd_x+4)*(mm+padd_e+4)*n*3*nt)
      CALL PUSHREAL8ARRAY(v, (lm+padd_x+4)*(mm+padd_e+4)*n*3)
      CALL PUSHREAL8ARRAY(u, (lm+padd_x+4)*(mm+padd_e+4)*n*3)
      CALL PUSHREAL8ARRAY(vbar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL PUSHREAL8ARRAY(ubar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL PUSHREAL8ARRAY(zeta, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL PUSHINTEGER4(may_day_flag)
      CALL PUSHBOOLEAN(synchro_flag)
      CALL PUSHREAL8(rdrg2)
      CALL PUSHINTEGER4(nnew)
      CALL PUSHINTEGER4(nrhs)
      CALL PUSHINTEGER4(nstp)
      CALL PUSHINTEGER4(knew)
      CALL PUSHREAL8(time)
      CALL PUSHREAL8ARRAY(hinv, (lm+padd_x+4)*(mm+padd_e+4))

      write(*,*) 'save------'
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer ntrds,trd,my_first,my_last,tile,range

      CALL POPREAL8ARRAY(hinv, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8(time)
      CALL POPINTEGER4(knew)
      CALL POPINTEGER4(nstp)
      CALL POPINTEGER4(nrhs)
      CALL POPINTEGER4(nnew)
      CALL POPREAL8(rdrg2)
      CALL POPBOOLEAN(synchro_flag)
      CALL POPINTEGER4(may_day_flag)
      CALL POPREAL8ARRAY(zeta, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL POPREAL8ARRAY(ubar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL POPREAL8ARRAY(vbar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL POPREAL8ARRAY(u, (lm+padd_x+4)*(mm+padd_e+4)*n*3)
      CALL POPREAL8ARRAY(v, (lm+padd_x+4)*(mm+padd_e+4)*n*3)
      CALL POPREAL8ARRAY(t, (lm+padd_x+4)*(mm+padd_e+4)*n*3*nt)
      CALL POPREAL8ARRAY(hz_bak, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL POPREAL8ARRAY(z_w, (lm+padd_x+4)*(mm+padd_e+4)*(n+1))
      CALL POPREAL8ARRAY(huon, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL POPREAL8ARRAY(hvom, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL POPREAL8ARRAY(hz, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL POPREAL8ARRAY(z_r, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL POPREAL8ARRAY(we, (lm+padd_x+4)*(mm+padd_e+4)*(n+1))
      CALL POPREAL8ARRAY(rho, (lm+padd_x+4)*(mm+padd_e+4)*n)
      CALL POPREAL8ARRAY(akv, (lm+padd_x+4)*(mm+padd_e+4)*(n+1))
      CALL POPREAL8ARRAY(akt, (lm+padd_x+4)*(mm+padd_e+4)*(n+1)*2)
      CALL POPREAL8ARRAY(rhoa, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(rhos, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(rufrc, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(rvfrc, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(rufrc_bak, (lm+padd_x+4)*(mm+padd_e+4)*2)
      CALL POPREAL8ARRAY(rvfrc_bak, (lm+padd_x+4)*(mm+padd_e+4)*2)
      CALL POPREAL8ARRAY(zt_avg1, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(du_avg1, (lm+padd_x+4)*(mm+padd_e+4)*5)
      CALL POPREAL8ARRAY(dv_avg1, (lm+padd_x+4)*(mm+padd_e+4)*5)
      CALL POPREAL8ARRAY(du_avg2, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(dv_avg2, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(a2d, n2d*nsa*npp)
      CALL POPREAL8ARRAY(a3d, n3d*7*npp)
      CALL POPREAL8ARRAY(ad_x_ru, n3d*npp)
      CALL POPREAL8ARRAY(ad_x_rv, n3d*npp)
      CALL POPREAL8ARRAY(sustr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(svstr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(bustr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(bvstr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPREAL8ARRAY(stflx, (lm+padd_x+4)*(mm+padd_e+4)*nt)
      CALL POPREAL8ARRAY(btflx, (lm+padd_x+4)*(mm+padd_e+4)*nt)
      CALL POPREAL8ARRAY(srflx, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPINTEGER4(ad_cost_counter)

      write(*,*) 'restore+++'
      write(*,*) '++++++++++'

      end subroutine restore_croco_state

      subroutine restore_croco_state_2d(tile)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

c      call set_depth(tile)

c      call grid_stiffness(tile)

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
      double precision ridge_width, ridge_height,depth,cff1,cff,rnd

# include "compute_extended_bounds.h"

c set a lower depth
      depth=2100
      ridge_width = 28.0D3
      ridge_height = 1400

      cff1=1./(ridge_width*ridge_width)
      do j=JstrR,JendR
        k=1
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

      ridge_width = 28.0D3
      ridge_height = 1400

      cff1=1./(ridge_width*ridge_width)
      j=JstrR
      k=1
      do i=IstrR,IendR
         if ((xr(i,j).gt.ridge_width).or.
     &        (xr(i,j).lt.(-ridge_width))) then
            cff=0.
         else
            cff=1-(xr(i,j)*xr(i,j)*cff1);
         endif
         ad_x(k) = depth-ridge_height*cff*cff
c     &        - 100;
         k=k+1
      enddo

      ad_array_node_size=k

      write(*,*) 'node,array_node_size, add_array_size/nnodes:',mynode,
     &     ad_array_node_size,
     &     ad_array_size/nnodes

      write(*,*) 'minad_x,maxad_x',minval(ad_x),maxval(ad_x)
      write(*,*) 'minh,maxh',minval(h),maxval(h)
      write(*,*) 'minobs,maxobs',minval(ad_obs),maxval(ad_obs)

      end subroutine init_local_arrays_2d_tile


      subroutine cost_fun_full_state(ad_x, cost)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision cost

      call cost_fun(ad_x, cost)

      end subroutine cost_fun_full_state
