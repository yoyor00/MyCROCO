c     this file is not differentiated
      
      subroutine check_driver(iicroot)

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
      integer k
      double precision f,xf,f0
      double precision fp,xfp,costb,cost
      double precision ad_eps,ad_check

      
      sim_iicroot = iicroot

      call init_control(ad_x_f,ad_g_f)

      write(*,*) 'start gradient ad_check'
      open(unit=777, file="zob.dat", action="write",
     &     form='FORMATTED')
      open(unit=778, file="h.dat", action="write",
     &     form='FORMATTED')

      ad_eps = 1D-5
      do while (ad_eps.gt.1D-10)

         call save_croco_state()
         call cost_fun(ad_x_f, f)
         call restore_croco_state()
         costb=1
         cost=f
         call cost_fun_b(ad_x_f,ad_g_f,cost,costb)
         call restore_croco_state()
#ifdef AMPI
         call AMPI_Barrier(MPI_COMM_WORLD,ierr)
         call AMPI_Reduce(f, xf, 1, AMPI_ADOUBLE_PRECISION, MPI_SUM, 0,
     &        MPI_COMM_WORLD, ierr)
         call AMPI_Reduce(ad_g_f,ad_sg_f,ad_array_size,
     &        AMPI_ADOUBLE_PRECISION,
     &        MPI_SUM, 0,
     &        MPI_COMM_WORLD, ierr)
         
         ad_g_f = ad_sg_f
         if (mynode.eq.0) then
            f = xf
         end if
#endif
c     perturbations
         do k=1,ad_array_size
            ad_x_f(k) = ad_x_f(k)+ad_eps
            call save_croco_state()
            call cost_fun(ad_x_f,fp)
            call restore_croco_state()
            ad_x_f(k) = ad_x_f(k)-ad_eps
#ifdef AMPI
            call AMPI_Barrier(MPI_COMM_WORLD,ierr)
            call AMPI_Reduce(fp, xfp, 1, AMPI_ADOUBLE_PRECISION,
     &           MPI_SUM, 0,
     &           MPI_COMM_WORLD, ierr)
            if (mynode.eq.0) then
               fp = xfp
            end if
#endif
            MPI_master_only ad_check=
     &           abs(ad_sg_f(k)-((fp-f)/ad_eps))/abs(ad_sg_f(k))
            
            MPI_master_only write(*,*) 'k,ad_eps,ad_check='
     &           ,k,ad_eps,ad_check
         
         end do
         ad_eps =ad_eps/2.
      end do
      close (unit=777)
      close (unit=778)
      write(*,*) 'end state_control'

      return
      end subroutine check_driver



