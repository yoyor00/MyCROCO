c     this file is not differentiated

      subroutine state_control(iicroot)

      implicit none
#include "adinc.h"
#include "adparam.h"

      integer iicroot
      integer ierr

      integer indic,imp,io,imode(3),omode,niter,nsim,iz(5),ndz,
     &     reverse

c     ndz must be >= 6*ad_array_size+1 in DIS mode
      parameter(ndz=6*ad_array_size+1)

      external simul_rc,euclid,ctonbe,ctcabe
      external simul
      character*3 normtype

      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost, cost_f

      double precision dx, df1,epsrel

      double precision dz(ndz)

      double precision control_to_param

      integer i

      ad_sim_iicroot = iicroot
      ad_cost_counter = 1
      ad_cost = 1
      ad_counter = 0


      write(*,*) 'ad_array_size=', ad_array_size
      write(*,*) 'nx:',2+lm+2+padd_x
      write(*,*) 'ne:',2+mm+2+padd_e

      call init_control(ad_x,ad_g)



c      MPI_master_only write(*,*) 'A1##>',ad_x_f
c      MPI_master_only write(*,*) 'B1##>',ad_g_f

      indic = 4
      call simul (indic,ad_array_size,ad_x,cost,ad_g,izs,rzs,dzs)

      ad_x_f(:) = 0.
      ad_g_f(:) = 0.

#ifdef MPI
#ifdef AMPI
      call AMPI_Gather(ad_x,ad_array_size/nnodes,
     &     AMPI_ADOUBLE_PRECISION,ad_x_f,ad_array_size/nnodes,
     &     AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call AMPI_Gather(ad_g,ad_array_size/nnodes,
     &     AMPI_ADOUBLE_PRECISION,ad_g_f,ad_array_size/nnodes,
     &     AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#else
      ad_x_f=ad_x
      ad_g_f=ad_g
#endif

#ifdef MPI
#ifdef AMPI
      call AMPI_Reduce(cost, cost_f, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &     0, MPI_COMM_WORLD, ierr)

      if (mynode.eq.0) then
         cost = cost_f
      end if
#endif
#endif
      ad_counter = ad_counter + 1

      dx=1.e-10
      df1=cost
      epsrel=1.e-4
      niter=200
      nsim=200
      normtype = 'dfn'
      io=6
      imp=5
c     imode(1)=0 ==> DIS mode
      imode(1)=0
c     imode(2)=0 ==> cold start
      imode(2)=0
c     imode(3)=0 ==> simulator never called with indic=1
      imode(3)=0

      write(*,*) 'start state_control'
      open(unit=777, file="coords.dat", action="write",
     &     form='FORMATTED')
      open(unit=778, file="icoords.dat", action="write",
     &     form='FORMATTED')
      open(unit=779, file="ctrl.dat", action="write",
     &     form='FORMATTED')
      open(unit=780, file="grad.dat", action="write",
     &     form='FORMATTED')


      reverse = 1
c     reverse mode
      do while (reverse.ge.1)
c         write(*,*) '+++>', mynode, ad_x(1), cost, gg(1)

c         MPI_master_only call optim_driver
c     &        (x,cost,gg,reverse,indic)

c         do i=1,ad_array_size
c            write (*,*) 'before: mynode, x, g',mynode,i,ad_x_f(i),
c     &           ad_g_f(i)
c         end do

c         MPI_master_only write(*,*) '1##>',ad_x_f
c         MPI_master_only write(*,*) '2##>',ad_sg_f


C        get ad_x_f, the new ad state on mpi master
         MPI_master_only call m1qn3 (simul_rc,euclid,ctonbe,ctcabe,
     &        ad_array_size,ad_x_f,
     &        cost,ad_g_f,dx,df1,epsrel,
     &        normtype,imp,io,imode,omode,niter,nsim,iz,dz,ndz,
     &        reverse,indic,izs,rzs,dzs)

c         do i=1,ad_array_size
c            write (*,*) 'after: mynode, x, g',mynode,i,
c     &           control_to_param(ad_x_f(i)),
c     &           ad_g_f(i)
c         end do
#ifdef MPI
#ifdef AMPI
         call AMPI_Scatter(ad_x_f,ad_array_size/nnodes,
     &        AMPI_ADOUBLE_PRECISION,ad_x,ad_array_size/nnodes,
     &        AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
c         call AMPI_Scatter(ad_g_f,ad_array_size/nnodes,
c     &        AMPI_ADOUBLE_PRECISION,ad_g,ad_array_size/nnodes,
c     &        AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#else
         ad_x=ad_x_f
         ad_g=ad_g_f
#endif
#ifdef MPI
#ifdef AMPI
         call AMPI_Bcast(reverse,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(indic,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)
#endif
#endif


C        get cost & ad_g, the gradient on each proc (& ad_rms)
         call simul(indic,ad_array_size,ad_x,cost,ad_g,izs,rzs,dzs)

         ad_g_f(:) = 0
#ifdef MPI
#ifdef AMPI
         call AMPI_Gather(ad_g,ad_array_size/nnodes,
     &        AMPI_ADOUBLE_PRECISION,ad_g_f,ad_array_size/nnodes,
     &        AMPI_ADOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#else
         ad_g_f=ad_g
#endif

#ifdef MPI
#ifdef AMPI
         call AMPI_Reduce(cost, cost_f, 1, AMPI_ADOUBLE_PRECISION,
     &        MPI_SUM,
     &        0, MPI_COMM_WORLD, ierr)

         if (mynode.eq.0) then
            cost = cost_f
         end if
#endif
#endif


         if (mynode.eq.0) then
            write (779,*) ad_counter,ad_rms,cost,norm2(ad_g_f),
     &           ad_dir_time,ad_adj_time
            call flush(779)
            write (780,*) ad_counter,ad_g_f
            call flush(780)
         end if

         df1 = cost

         ad_counter = ad_counter+1
      end do

C     ...
      call set_state(ad_x_f)

      close (unit=777)
      close (unit=778)
      close (unit=779)
      close (unit=780)
      write(*,*) 'end state_control'

      return
      end subroutine state_control



      subroutine handle_err(errcode)
      implicit none
# include "netcdf.inc"
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end


      function rms()
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision rms,irms,crms
      integer i,j,ne,cne,ierr

      ne = 0
      irms = 0
      rms = 0
      cne = 0
      crms = 0

      call rms_fun_step(irms,ne)

#ifdef MPI
#ifdef AMPI
      call AMPI_Reduce(irms, crms, 1, AMPI_ADOUBLE_PRECISION,
     &     MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call AMPI_Reduce(ne, cne, 1, MPI_INTEGER,
     &     MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif
#else
      crms = irms
      cne = ne
#endif

      if (mynode.eq.0) then
         if (cne.gt.0) then
            rms = sqrt(crms)/cne
         else
            rms = 0
         end if
      end if

      end function rms

      subroutine rms_fun_step(irms,ne)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision irms
      integer ne
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
         call rms_fun_step_2d(tile,irms,ne)
      end do

      end subroutine rms_fun_step

      subroutine rms_fun_step_2d(tile,irms,ne)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile
      double precision irms
      integer ne
# include "compute_tile_bounds.h"

      call rms_fun_step_2d_tile(Istr,Iend,Jstr,Jend,irms,ne)

      end subroutine rms_fun_step_2d

      subroutine rms_fun_step_2d_tile(Istr,Iend,Jstr,Jend,
     &     irms,ne)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,i,j
      integer k,ne
      double precision delta_z, delta_zob, irms

      do j=Jstr,Jend
         do i=Istr,Iend
            if (ad_obs(i,j,1).ne.spval) then
               delta_z = (zeta(i,j,knew)-ad_obs(i,j,ad_nt+2))
               irms = irms + delta_z*delta_z
               ne = ne+1
            end if
         end do
      end do

      end subroutine rms_fun_step_2d_tile
