c     this file is not differentiated
      
      subroutine check_driver(iicroot)

      implicit none
# include "adinc.h"
# include "adparam.h"
      integer iicroot
      integer ierr
      integer k
      double precision icost,costb,cost,cost_p
      double precision ad_eps,ad_check_tgt,n_ad_g
      double precision ad_check_div
      double precision n_ad_g_f,ad_d_f,cost_f,cost_p_f
      external cost_fun_d

      
c     perturbated ad_x
      double precision ad_d(1)
      double precision ad_x_p(ad_array_size/nnodes)

      double precision control_to_param,rand
      
      ad_sim_iicroot = iicroot

      open(unit=789, file="ad_check.dat", action="write",
     &     form='FORMATTED')
      
      call init_control(ad_x,ad_g)

      write(*,*) 'start gradient ad_check'

c     1. compute gradient with adjoint
      icost=0
      call save_croco_state()
      call cost_fun(ad_x, icost)
      call restore_croco_state()
      costb=1
      cost=icost
      call save_croco_state()
      call cost_fun_b(ad_x,ad_g,cost,costb)
      call restore_croco_state()

         
c     2. compute tangent          
      ad_d(1)=0

      call save_croco_state()
      call cost_fun_d(ad_x, ad_g, cost, ad_d)
      call restore_croco_state()
      
c     n_ad_g = < ad_g . ad_g >
      n_ad_g = 0.
      do k=1,ad_array_node_size
         n_ad_g = n_ad_g + ad_g(k)*ad_g(k)
      end do

      if (ad_d(1).gt.epsilon(ad_d(1))) then
         ad_check_tgt = (n_ad_g / ad_d(1))
      else
         ad_check_tgt = 0
      end if

#if defined(MPI) && defined(AMPI)
      call AMPI_Reduce(n_ad_g, n_ad_g_f, 1,
     &     AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &     0, MPI_COMM_WORLD, ierr)
      call AMPI_Reduce(ad_d, ad_d_f, 1,
     &     AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &     0, MPI_COMM_WORLD, ierr)
      call AMPI_Reduce(cost, cost_f, 1,
     &     AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &     0, MPI_COMM_WORLD, ierr)
#else
      n_ad_g_f = n_ad_g
      ad_d_f = ad_d(1)
      cost_f = cost
#endif

      ad_eps = 0.
      if (mynode.eq.0) then
         if (abs(ad_d_f).gt.1D-20) then
            ad_check_tgt = n_ad_g_f/ad_d_f
         else
            ad_check_tgt = 0
         end if
         write(789,*)
     &        ad_eps,ad_check_tgt
         call flush(789)
      end if
      
c     3. finite difference estimation
      ad_eps = 1.D-4
      do while (ad_eps.gt.1D-10)   
         
         do k=1,ad_array_node_size
            ad_x_p(k) = ad_x(k)+ad_eps*ad_g(k)
         end do

         icost=0
         call save_croco_state()
         call cost_fun(ad_x_p,icost)
         call restore_croco_state()
         cost_p = icost

         
         write(*,*) 'n_ad_g:',mynode,n_ad_g
         write(*,*) 'cost:',mynode,cost
         write(*,*) 'cost_p:',mynode,cost_p

         write(*,*) '==check==',mynode,ad_check_div

#if defined(MPI) && defined(AMPI)
         call AMPI_Reduce(cost_p, cost_p_f, 1,
     &        AMPI_ADOUBLE_PRECISION, MPI_SUM,
     &        0, MPI_COMM_WORLD, ierr)
#else
         cost_p_f = cost_p
#endif

         if (mynode.eq.0) then
            if (abs(cost_p_f - cost_f).gt.1D-48) then
               ad_check_div=(n_ad_g_f / ((cost_p_f-cost_f)/ad_eps))
            else
               ad_check_div=0
            end if
            
            write(789,*)
     &           ad_eps,ad_check_div
            call flush(789)
         end if
         ad_eps =ad_eps/2.
      end do
      close (unit=789)
      write(*,*) 'end state_control'

      return
      end subroutine check_driver



