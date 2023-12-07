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

c$$$      call init_scalars(ierr)
c$$$
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call start_timers()
c$$$        call init_arrays (tile)
c$$$      enddo
c$$$
c$$$c      call check_tab2d(h(:,:),'h initialisation #1','r')
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call setup_grid1 (tile)
c$$$      enddo
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call setup_grid2 (tile)
c$$$      enddo
c$$$
c$$$      call set_scoord
c$$$      call set_weights
c$$$
c$$$      iif_bak = iif
c$$$      iif = -1
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call set_depth (tile)
c$$$      enddo
c$$$
c$$$      iif = iif_bak
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call grid_stiffness (tile)
c$$$      enddo
c$$$
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$c#ifdef INTERNAL
c$$$c#endif
      do tile=0,NSUB_X*NSUB_E-1
         call set_state_2d(tile, ad_x)
      enddo
c      do tile=0,NSUB_X*NSUB_E-1
c        call ana_initial (tile)
c      enddo

      
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$        call set_HUV (tile)
c$$$# ifdef RESET_RHO0
c$$$c        call reset_rho0 (tile)
c$$$# endif
c$$$      enddo
c$$$CR      write(*,*)  ' -4' MYID
c$$$C$OMP PARALLEL DO PRIVATE(tile)
c$$$      do tile=0,NSUB_X*NSUB_E-1
c$$$         call omega (tile)
c$$$         call rho_eos (tile)
c$$$      enddo
c$$$CR      write(*,*)  ' -3' MYID


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

      subroutine set_state_2d_1_tile(Istr,Iend,Jstr,Jend, ad_x)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k,ik,cpoint,itrc
      double precision cff1,cff2,cff
      character*6 counter_str

c      rdrg = ad_x(1)
      if (nrrec.eq.0) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            zeta(i,j,1)=0.D0
            zeta(i,j,2)=0.D0
            ubar(i,j,1)=0.D0
            vbar(i,j,1)=0.D0
            ubar(i,j,2)=0.D0
            vbar(i,j,2)=0.D0
          enddo
        enddo
      else
      endif
      if (nrrec.eq.0) then
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              u(i,j,k,1)=0.D0
              v(i,j,k,1)=0.D0
              u(i,j,k,2)=0.D0
              v(i,j,k,2)=0.D0
            enddo
          enddo
        enddo
      else
      endif
      
c      if (nrrec.eq.0) then
      cff1=ad_x(1)
      write(*,*) 'cff1:', cff1
      cff2=cff1*(rho0*800.D0/g)*(5.0D-5/((42.689D0/44.690D0)**2))
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              t(i,j,k,1,itemp)=cff2*exp(z_r(i,j,k)/800.D0)
     &                   *(0.6D0-0.4D0*tanh(z_r(i,j,k)/800.D0))
              t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
             enddo
          enddo
        enddo
c     endif
#if defined(MPI)        
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(:, : ,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(:, :,1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(:, :,1))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(:, :,2))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(:, :,2))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &     vbar(:, :,2))
c      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
c     &     u(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(:, :,:,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(-1:,-1:,1:,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(-1:,-1:,1:,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(-1:,-1:,1:,2))
      do itrc=1,NT
         if (.not.got_tini(itrc)) then
          call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                            t(-1:,-1:,1:,1,itrc))
          call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                            t(-1:,-1:,1:,2,itrc))
         endif
      enddo
#endif     
      return

#endif
      end subroutine set_state_2d_1_tile
