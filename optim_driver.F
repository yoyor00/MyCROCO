c     this file is not differentiated
      
# include "cppdefs.h"

      subroutine state_control(iicroot)

      implicit none
# include "cppdefs.h"
# include "param.h"
# include "adparam.h"
# include "scalars.h"
# ifdef MPI
# include "mpi_roms.h"
# ifdef AMPI
# include "ampi/ampif.h"
# else
      include 'mpif.h'
# endif
# endif
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

      double precision f, xf

      double precision dx, df1,epsrel

      double precision dz(ndz)

      integer i

      sim_iicroot = iicroot

      call init_control()
      
      indic = 4
      call simul (indic,ad_array_size,ad_x_f,f,ad_g_f,izs,rzs,dzs)
#ifdef AMPI
      call AMPI_Barrier(MPI_COMM_WORLD,ierr)
      call AMPI_Reduce(f, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)
      call AMPI_Reduce(ad_g_f,ad_sg_f,ad_array_size,
     &     AMPI_ADOUBLE_PRECISION,
     &     MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)
      
      if (mynode.eq.0) then
         f = xf
      end if
#endif
      dx=1.e-10
      df1=f
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
      open(unit=777, file="zob.dat", action="write",
     &     form='FORMATTED')
      open(unit=778, file="h.dat", action="write",
     &     form='FORMATTED')


      reverse = 1
c     reverse mode
      do while (reverse.ge.1)
c         write(*,*) '+++>', mynode, ad_x(1), f, gg(1)

c         MPI_master_only call optim_driver
c     &        (x,f,gg,reverse,indic)

         do i=1,ad_array_size
            write (*,*) 'before: mynode, x, g',mynode,i,ad_x_f(i),ad_g_f(i)
         end do
         
         MPI_master_only call m1qn3 (simul_rc,euclid,ctonbe,ctcabe,
     &        ad_array_size,ad_x_f,
     &        f,ad_g_f,dx,df1,epsrel,
     &        normtype,imp,io,imode,omode,niter,nsim,iz,dz,ndz,
     &        reverse,indic,izs,rzs,dzs)

#ifdef AMPI
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(ad_x_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(ad_g_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &        0,MPI_COMM_WORLD,ierr)
c        call AMPI_Bcast(ad_x_f,ad_array_size,AMPI_ADOUBLE_PRECISION,
c     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(reverse,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(indic,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)
#endif
         call simul(indic,ad_array_size,ad_x_f,f,ad_g_f,izs,rzs,dzs)

#ifdef AMPI
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Reduce(f, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, 0,
     &        MPI_COMM_WORLD, ierr)
         call AMPI_Reduce(ad_g_f,ad_sg_f,ad_array_size,
     &     AMPI_ADOUBLE_PRECISION,
     &     MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)
#endif
         if (mynode.eq.0) then
            f = xf
         end if

         df1 = f

         do i=1,ad_array_size
            write (*,*) 'after: mynode,x,g',mynode,i,ad_x_f(i),ad_g_f(i)
         end do
         
      end do

C ...
      call set_state(ad_x_f)

      close (unit=777)
      close (unit=778)
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
      
