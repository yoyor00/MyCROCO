#include "cppdefs.h"
#ifdef NBQ
!
!======================================================================
!                      NBQ-Mode for NH-modeling
!                            Main Routine
!======================================================================
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief SNBQ driver : Non-hydrostatic algorithm with the 
!                       Non-boussinesq solver.
!
!> @details NBQ time step. See SNBQ web pages :
!  http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Html_pages
!    Algorithme_NBQ.htm      --> SNBQ algorithm:               
!    Algebrique_SNBQ.htm     --> SNBQ algebric representation
!    Couplage_Numerique.htm  --> Numerical coupling
!    Couplage_Modes_SNBQ.htm --> Coupling:
!    Couplage_Split_NBQ.htm  --> Coupling Splitting
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!======================================================================
!
      subroutine step3d_fb_nbq (Istr,Iend,Jstr,Jend, WORK)

      use module_nh 
      use module_nbq

      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"
# ifdef MPI
      include 'mpif.h'
# endif
      integer Istr,Iend,Jstr,Jend,i,j,k
      real    WORK(PRIVATE_2D_SCRATCH_ARRAY)
      real :: dum_s
# ifndef MPI
      integer mynode
      mynode=0
# endif
# undef DEBUG

!
!-------------------------------------------------------------------
!       Initialization of various test-cases
!-------------------------------------------------------------------
!       
        if (iif==1.and.iic==1) call initial_nh (3)
!
!-------------------------------------------------------------------
!  Get internal and external forcing terms for nbq equations:
!  ru+rubar (or lambda_ext+lambda_int)
!  dzdt*rhosurf
!-------------------------------------------------------------------
!
      call ru_nbq(1)
      call density_nbq(1)
!
!------------------------------------------------------------------
!       Implicit part: system setup
!-------------------------------------------------------------------
!
# ifdef NBQ_IMP
      if (iif.eq.1.and.ifl_imp_nbq.eq.1) call implicit_nbq (1)
# endif
!
!*******************************************************************
!*******************************************************************
!              NBQ mode iteration (main loop)
!*******************************************************************
!*******************************************************************

      do iteration_nbq=1,iteration_nbq_max

!.......1st iteration NBQ: treats NBQ
        if (iteration_nbq>1) then
!        Momentum equation: switch time indices (move forward)
         call ru_nbq(7)
        endif
!
!-------------------------------------------------------------------
!      Compute pressure gradient and gravity terms (AMUX)
!                rhp_nbq_a  ==> rhs1_nbq
!-------------------------------------------------------------------
!
        call ru_nbq(6)
!
!-------------------------------------------------------------------
!      Horizontal Momentum equation: leapfrog time stepping
!         If explicit: (x,y,z) is dealt with here
!-------------------------------------------------------------------
!
!  XI-Direction:
!
        do l_nbq = 1,nequ_nh(7)
          dum_s =             soundspeed2_nbq  * rhs1_nbq (l_nbq)                      &
                            - visc2_nbq * rhsd2_nbq(l_nbq) 
          qdm_nbq_a(l_nbq,vnnew_nbq) = qdm_nbq_a(l_nbq,vnrhs_nbq)  + dtnbq*(           &
                                       dum_s                                           &
                                  + dqdmdt_nbq_a(l_nbq)                     )  
          rhssum_nbq_a(l_nbq) = rhssum_nbq_a(l_nbq)  +  dum_s        
        enddo 
!
!  U-momentum open boundary conditions
!
# ifdef OBC_NBQ
          call unbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
          
!
!  Message passing: Send U (51) 
!
        call parallele_nbq(51)
!       call parallele_nbq(151)
!
!  ETA-Direction:
!
         do l_nbq = nequ_nh(7)+1,neqv_nh(7)  
           dum_s =             soundspeed2_nbq  * rhs1_nbq (l_nbq)                      &
                             - visc2_nbq * rhsd2_nbq(l_nbq) 
           qdm_nbq_a(l_nbq,vnnew_nbq) = qdm_nbq_a(l_nbq,vnrhs_nbq)  + dtnbq*(            &
                                dum_s                                                    &
                              + dqdmdt_nbq_a(l_nbq)                        )  
           rhssum_nbq_a(l_nbq) = rhssum_nbq_a(l_nbq)  +  dum_s    
         enddo 
!
!  V-momentum open boundary conditions
!
# ifdef OBC_NBQ
            call vnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!  Message passing: Send V (52) 
!
         call parallele_nbq(52)
!        call parallele_nbq(152)
!
!-------------------------------------------------------------------
!      Vertical Momentum equation: leapfrog time stepping
!         If explicit: (x,y,z) is dealt with here
!         If implicit: (x,y)   only
!-------------------------------------------------------------------
!
# ifndef NBQ_IMP
!
!  Z-Direction: Explicit
!
           do l_nbq = neqv_nh(7)+1,neqw_nh(7)
             dum_s =             soundspeed2_nbq  * rhs1_nbq (l_nbq)                   &
                               - visc2_nbq * rhsd2_nbq(l_nbq)
             qdm_nbq_a(l_nbq,vnnew_nbq) = qdm_nbq_a(l_nbq,vnrhs_nbq)  + dtnbq*(        &
                                  dum_s                                                &
                                + dqdmdt_nbq_a(l_nbq)                          )  
             rhssum_nbq_a(l_nbq) = rhssum_nbq_a(l_nbq)  +  dum_s
           enddo 
# else
           call parallele_nbq(151)  ! u only 
           call parallele_nbq(152)  ! v only 
!          Receive:
           call implicit_nbq (2)
# endif
!
!      Vertical momentum open boundary conditions
!
# ifdef OBC_NBQ
!       call wnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!-------------------------------------------------------------------
!      Message passing 
!-------------------------------------------------------------------
!
!  Send
        call parallele_nbq(53)    ! w only 
!       call parallele_nbq(153) 
!
!  Receive
# ifndef NBQ_IMP
        call parallele_nbq(151)     
        call parallele_nbq(152)    
# endif 
        call parallele_nbq(153)     

!
!-------------------------------------------------------------------
!      Compute divergence term (AMUX):
!          qdm_nbq_a ==> div_nbq_a
!-------------------------------------------------------------------
!
        call density_nbq(60)
!
!-------------------------------------------------------------------
!      Message passing
!-------------------------------------------------------------------
!
!  Send
        call parallele_nbq(7) 

!  Receive
        call parallele_nbq(17) 

!
!-------------------------------------------------------------------
!      Compute Second viscosity (product mat*vect):
!            div_nbq_a ==> rhsd2_nbq
!-------------------------------------------------------------------
!
        call viscous_nbq (2)
!
!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
#  ifdef ACOUSTIC
        call density_nbq(11)
#  endif
!
!-------------------------------------------------------------------
!      Mass equation: leapfrog time stepping:
!-------------------------------------------------------------------
!
        do l_nbq=1,neqcont_nh
          rhp_nbq_a(l_nbq,rnnew_nbq) = rhp_nbq_a(l_nbq,rnrhs_nbq)         &
                             - div_nbq_a(l_nbq,dnrhs_nbq) * dtnbq 
  !        i     = l2iq_nh (l_nbq)
  !        j     = l2jq_nh (l_nbq)
  !        k     = l2kq_nh (l_nbq)
  !        if (i==13.and.j==6) rhp_nbq_a(l_nbq,rnnew_nbq) = 0.
  !        if (i==14.and.j==6) rhp_nbq_a(l_nbq,rnnew_nbq) = 0.
          
        enddo

!
!-------------------------------------------------------------------
!      Density open boundary conditions
!-------------------------------------------------------------------
!
# ifdef OBC_NBQ
!        call rnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!-------------------------------------------------------------------
!       Mass equation: switch time indices (move forward)
!-------------------------------------------------------------------
!
        call density_nbq(7)
!
!*******************************************************************
!*******************************************************************
      enddo    ! NBQ loop

!*******************************************************************
!*******************************************************************
!
!
!-------------------------------------------------------------------
!......Move forward: momentum
!-------------------------------------------------------------------
!
       call ru_nbq(7)
!
!-------------------------------------------------------------------
!......Set NBQ/EXT coupling terms
!-------------------------------------------------------------------
!
      call ru_nbq(2)

#ifdef RVTK_DEBUG
!       call check_tab2d(rubar_nbq(:,:),'rubar_nbq step3d_nbq','u')
!       call check_tab2d(rvbar_nbq(:,:),'rvbar_nbq step3d_nbq','v')
!       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext step3d_nbq','r')
#endif      
      call density_nbq(20)

   
      end subroutine step3d_fb_nbq

#else
      subroutine step3d_fb_nbq_empty
      end subroutine step3d_fb_nbq_empty
#endif
