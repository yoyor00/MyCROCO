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

      double precision cost, xf, rms, crms, rms0

      double precision dx, df1,epsrel

      double precision dz(ndz)

      integer i
      
      sim_iicroot = iicroot
      ad_cost_counter = 1
      ad_cost = 1
      ad_counter = 0

      
      write(*,*) 'ad_array_size=', ad_array_size
      write(*,*) 'nx:',2+lm+2+padd_x
      write(*,*) 'ne:',2+mm+2+padd_e
      
      call init_control(ad_x_f,ad_g_f)

      MPI_master_only write(*,*) 'A1##>',ad_x_f
      MPI_master_only write(*,*) 'B1##>',ad_g_f
      
      indic = 4
      call simul (indic,ad_array_size,ad_x_f,cost,ad_g_f,izs,rzs,dzs)
#ifdef AMPI
      call AMPI_Barrier(MPI_COMM_WORLD,ierr)
      call AMPI_Reduce(cost, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)
      call AMPI_Reduce(ad_g_f,ad_sg_f,ad_array_size,
     &     AMPI_ADOUBLE_PRECISION,
     &     MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)

      ad_counter = ad_counter + 1

      if (mynode.eq.0) then
         cost = xf
      end if
#endif
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


      reverse = 1
c     reverse mode
      do while (reverse.ge.1)
c         write(*,*) '+++>', mynode, ad_x(1), cost, gg(1)

c         MPI_master_only call optim_driver
c     &        (x,cost,gg,reverse,indic)

         do i=1,ad_array_size
            write (*,*) 'before: mynode, x, g',mynode,i,ad_x_f(i),
     &           ad_g_f(i)
         end do

         MPI_master_only write(*,*) '1##>',ad_x_f
         MPI_master_only write(*,*) '2##>',ad_sg_f
         

C        get ad_x_f, the new ad state on mpi master
         MPI_master_only call m1qn3 (simul_rc,euclid,ctonbe,ctcabe,
     &        ad_array_size,ad_x_f,
     &        cost,ad_sg_f,dx,df1,epsrel,
     &        normtype,imp,io,imode,omode,niter,nsim,iz,dz,ndz,
     &        reverse,indic,izs,rzs,dzs)


#ifdef AMPI
C        broadcast the new ad state
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(ad_x_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(ad_sg_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(reverse,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(indic,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)
#endif
         ad_g_f = ad_sg_f
         
C        get ad_g_f the gradient on each proc (& ad_rms)
         call simul(indic,ad_array_size,ad_x_f,cost,ad_g_f,izs,rzs,dzs)

#ifdef AMPI
C        sum all proc gradients & rms on master
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Reduce(cost, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Reduce(ad_rms, rms0, 1, AMPI_ADOUBLE_PRECISION,
     &        MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call AMPI_Reduce(ad_g_f,ad_sg_f,ad_array_size,
     &     AMPI_ADOUBLE_PRECISION,
     &     MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)
#endif
         if (mynode.eq.0) then
            cost = xf
            write (779,*) ad_counter,rms0,cost,norm2(ad_sg_f)
            call flush(779)
         end if
         
         df1 = cost

         do i=1,ad_array_size
            write (*,*) 'after: mynode,x,g',mynode,i,ad_x_f(i),ad_g_f(i)
         end do
         
         ad_counter = ad_counter+1
      end do

C     ...
      call set_state(ad_x_f)

      close (unit=777)
      close (unit=778)
      close (unit=779)
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

      double precision rms,irms
      integer i,j
      double precision delta_z

      irms = 0

      call rms_fun_step(irms)

      rms=sqrt(irms)
      
      end function rms

      subroutine rms_fun_step(irms)
      implicit none
# include "adinc.h"
# include "adparam.h"

      double precision irms
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
         call rms_fun_step_2d(tile,irms)
      end do

      end subroutine rms_fun_step

      subroutine rms_fun_step_2d(tile,irms)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile
      double precision irms
# include "compute_tile_bounds.h"

      call rms_fun_step_2d_tile(Istr,Iend,Jstr,Jend,irms)

      end subroutine rms_fun_step_2d

      subroutine rms_fun_step_2d_tile(Istr,Iend,Jstr,Jend,
     &     irms)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,i,j
      integer k
      double precision delta_z, delta_zob, irms

      do j=Jstr,Jend
         do i=Istr,Iend
            if (ad_obs(i,j,1).ne.spval) then
               delta_z = (zeta(i,j,knew)-ad_obs(i,j,ad_nt+2))
               irms = irms + delta_z*delta_z
            end if
         end do
      end do

      end subroutine rms_fun_step_2d_tile
