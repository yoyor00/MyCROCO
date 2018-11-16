c     this file is not differentiated
      
# include "cppdefs.h"

      subroutine optim_driver(ad_x,f,gg,reverse,indic)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "adparam.h"
# ifdef MPI
# include "mpi_roms.h"
# ifdef AMPI
# include "ampi/ampif.h"
# else
      include 'mpif.h'
# endif
# endif

      integer indic,imp,io,imode(3),omode,niter,nsim,iz(5),ndz,
     &     reverse

      parameter(ndz=10*ad_array_size*ad_array_size+1)

      external simul_rc,euclid,ctonbe,ctcabe
      external simul
      character*3 normtype

      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision f
      double precision gg(ad_array_size)

      double precision dx, df1,epsrel,
     &     z(5*ad_array_size+1)

      double precision dz(ndz)

      integer i

      dx=1.e-10
      df1=f
      epsrel=1.e-2
      niter=200
      nsim=200
      normtype = 'dfn'
      io=6
      imp=5
      imode(1)=0
      imode(2)=0
      imode(3)=0

      do i=1,ad_array_size
         write (*,*) 'before: x, gg', i, ad_x(i), gg(i)
      end do

      call m1qn3 (simul,euclid,ctonbe,ctcabe,ad_array_size,ad_x,
     &     f,gg,dx,df1,epsrel,
     &     normtype,imp,io,imode,omode,niter,nsim,iz,dz,ndz,
     &     reverse,indic,izs,rzs,dzs)


      do i=1,ad_array_size
         write (*,*) 'after: x, gg', i, ad_x(i), gg(i)
      end do


      end subroutine optim_driver



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

      parameter(ndz=10*ad_array_size*ad_array_size+1)

      external simul_rc,euclid,ctonbe,ctcabe
      external simul
      character*3 normtype

      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision f, xf
      double precision gg(ad_array_size), xgg(ad_array_size)

      double precision dx, df1,epsrel,
     &     z(5*ad_array_size+1)

      double precision dz(ndz)

      integer i

      sim_iicroot = iicroot

      do i=1,ad_array_size
         gg(i) = 0.
      end do

      call init_control(ad_x)
      
      indic = 4
      call simul (indic,ad_array_size,ad_x,f,gg,izs,rzs,dzs)
#ifdef AMPI
      call AMPI_Barrier(MPI_COMM_WORLD,ierr)
      call AMPI_Reduce(f, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)
      call AMPI_Reduce(gg, xgg, ad_array_size,AMPI_ADOUBLE_PRECISION,
     &     MPI_SUM, 0,
     &     MPI_COMM_WORLD, ierr)

      if (mynode.eq.0) then
         write(*,*) 'NODE 0:', xf, xgg(1)
         f = xf
         gg(:) = xgg(:)
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
      imode(1)=0
      imode(2)=0
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

         MPI_master_only call m1qn3 (simul_rc,euclid,ctonbe,ctcabe,
     &        ad_array_size,ad_x,
     &        f,gg,dx,df1,epsrel,
     &        normtype,imp,io,imode,omode,niter,nsim,iz,dz,ndz,
     &        reverse,indic,izs,rzs,dzs)

c         write(*,*) '===>', mynode, ad_x(1), f, gg(1)
#ifdef AMPI
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(ad_x,ad_array_size,AMPI_ADOUBLE_PRECISION,
     &        0,MPI_COMM_WORLD,ierr)
         call AMPI_Bcast(reverse,1,MPI_INTEGER,
     &        0,MPI_COMM_WORLD,ierr)   
#endif
         call simul(indic,ad_array_size,ad_x,f,gg,izs,rzs,dzs)

#ifdef AMPI
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Reduce(f, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, 0,
     &        MPI_COMM_WORLD, ierr)
         call AMPI_Reduce(gg, xgg, ad_array_size,AMPI_ADOUBLE_PRECISION,
     &        MPI_SUM, 0,
     &        MPI_COMM_WORLD, ierr)

         if (mynode.eq.0) then
            write(*,*) 'NODE 0:', xf, xgg(1)
            f = xf
            gg(:) = xgg(:)
         end if
#endif

         df1 = f

      end do
C ...
      call set_state(ad_x)

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
      
