c     this file is not differentiated
      
      subroutine check_driver(iicroot)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer iicroot
      integer ierr
      integer k
      double precision icost,costb,cost,cost_p
      double precision ad_eps,ad_check,ad_check_f,s_ad_g

c     perturbated ad_x
      double precision ad_x_p(ad_array_size/nnodes)
      
      sim_iicroot = iicroot

      open(unit=789, file="ad_check.dat", action="write",
     &     form='FORMATTED')
      
      call init_control(ad_x,ad_g)

      write(*,*) 'start gradient ad_check'

c     1. compute gradient with adjoint
      call save_croco_state()
      icost=0
      call cost_fun(ad_x, icost)
      call restore_croco_state()
      costb=1
      cost=icost
      call cost_fun_b(ad_x,ad_g,cost,costb)
      call restore_croco_state()

c     ad_g . (1 1 1 ... )
      s_ad_g = sum(ad_g)
      
c     perturbations      
      ad_eps = 1D-2
      do while (ad_eps.gt.1D-10)   
         
         do k=1,ad_array_size/nnodes
            ad_x_p(k) = ad_x(k)+ad_eps
         end do
         
         call save_croco_state()
         icost=0
         call cost_fun(ad_x_p,icost)
         call restore_croco_state()

         cost_p = icost
         
         ad_check=(s_ad_g-(cost_p-cost)/ad_eps)/s_ad_g

         write(*,*) '==check==',mynode,ad_check

#if defined(MPI) && defined(AMPI)
         call AMPI_Reduce(ad_check, ad_check_f, 1,
     &        AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &        0, MPI_COMM_WORLD, ierr)
#endif

         if (mynode.eq.0) then
            ad_check_f = ad_check_f/nnodes
            write(789,*)
     &           ad_eps,ad_check_f
            call flush(789)
         end if
         ad_eps =ad_eps/2.
      end do
      close (unit=789)
      write(*,*) 'end state_control'

      return
      end subroutine check_driver



