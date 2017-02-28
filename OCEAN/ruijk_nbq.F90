#include "cppdefs.h"
#if defined NBQ && defined NBQ_IJK

      subroutine ruijk_nbq(icall,Istr,Iend,Jstr,Jend)

!**********************************************************************
!
!                      Various Computations related to
!                            NBQ momentum
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!
! DESCRIPTION: 
!
!> @brief NBQ momentum routine.
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!**********************************************************************

      use module_nh 
      use module_nbq
      implicit none

      
# include "param_F90.h"
# include "scalars_F90.h"
# include "grid.h"
# include "ocean3d.h"
# include "nbq.h"

      integer :: i,j,k
      integer :: icall,ncp
      real    :: cff
      integer :: Istr,Iend,Jstr,Jend
      real WORK(PRIVATE_2D_SCRATCH_ARRAY)      


      if (icall.eq.1) then      ! TBD
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!      call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
!                                       ruint_nbq(START_2D_ARRAY,1))
!      call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
!                        ruext_nbq(START_2D_ARRAY,1))
!      call exchange_u3d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &
!                                       rvint_nbq(START_2D_ARRAY,1))
!      call exchange_u3d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &
!                        rvext_nbq(START_2D_ARRAY,1))
!      call exchange_u3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,  &
!                                       rwint_nbq(START_2D_ARRAY,0))
#endif

      elseif (icall.eq.2) then
!
!*******************************************************************
!  Prepare feedback of NBQ rhs terms to external equations
!*******************************************************************
!
        cff=1./(rho0*real(ndtnbq))
     
!        
! X-direction:
!
        rubar_nbq(:,:)=0.
        do k=1,N         
          do j=JstrU_nh,JendU_nh   
            do i=IstrU_nh,IendU_nh
              ru_nbq_ext (i,j,k) = cff*rhssumu_nbq(i,j,k)*on_u(i,j)*om_u(i,j)
              rhssumu_nbq(i,j,k) = 0.
              rubar_nbq(i,j)     = rubar_nbq(i,j)+ru_nbq_ext(i,j,k)             
            enddo
          enddo
        enddo        

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI            
!      call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &      ! TBD
!                                       ru_nbq_ext(START_2D_ARRAY,1))
!      call exchange_u2d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
!                        rubar_nbq(START_2D_ARRAY))
#endif

#ifdef RVTK_DEBUG
!       call check_tab3d(ru_nbq_ext(:,:,1:N),'ru_nbq_ext (ru_nbq)','u')
!       call check_tab2d(rubar_nbq(:,:),'rubar_nbq (ru_nbq)','u')
#endif    
!
! Y-direction:
!
        rvbar_nbq(:,:)=0.
        do k=1,N         
          do j=JstrV_nh,JendV_nh   
            do i=IstrV_nh,IendV_nh
              rv_nbq_ext (i,j,k) = cff*rhssumv_nbq(i,j,k)*om_v(i,j)*on_v(i,j)
              rhssumv_nbq(i,j,k) = 0.
              rvbar_nbq(i,j)     = rvbar_nbq(i,j)+rv_nbq_ext(i,j,k)
            enddo
          enddo
        enddo
        

#ifdef M2FILTER_NONE
        if (LAST_2D_STEP) then
#endif 
        
!    
! Z-direction:
!

        do j=Jstr_nh,Jend_nh             
          do i=Istr_nh,Iend_nh
            WORK(i,j) = cff*on_r(i,j)*om_r(i,j)
          enddo
        enddo
          
        do k=0,N 
          do j=Jstr_nh,Jend_nh             
            do i=Istr_nh,Iend_nh
              rw_nbq_ext (i,j,k) = ((qdmw_nbq(i,j,k)-rw_nbq_ext(i,j,k))/dtnbq-ndtnbq*rho0*rwint_nbq(i,j,k))*WORK(i,j)
            enddo
          enddo
        enddo
#ifdef M2FILTER_NONE        
        endif
#endif        

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!      call exchange_r3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,  &        ! TBD
!                                       rw_nbq_ext(START_2D_ARRAY,0)) 
#endif
#ifdef RVTK_DEBUG
!       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext (ru_nbq)','v')
#endif    

      elseif (icall.eq.6) then
!
!*******************************************************************
!  Increment momentum :
!*******************************************************************
!
       
       endif  ! icall

       return
       end
#else
      subroutine ruijk_nbq_empty 
      return
      end
#endif

