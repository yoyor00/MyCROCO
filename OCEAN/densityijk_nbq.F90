#include "cppdefs.h"
#if defined NBQ && defined NBQ_IJK

      subroutine densityijk_nbq(icall)

!======================================================================
!
!                      Various Computations related to
!                            NBQ density
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief Density anomaly driver 
!
!> @details 
! - icall=0 vertical averaging density to get the external 
!           mode density (rhpbar_nbq_t).
! - icall=1 internal mode density (rhpio_nbq_t).
! - icall=5 computes the NH pressure gradient for the internal mode.
! - icall=6 computes the divergence of the momentum (div_nbq_a) 
!      + used in the momentum equation (second viscosity term with 
!        gradient operator)
!      + used in the continuity equation
! - icall=7 time incrementation of NBQ mode density anomaly.
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!======================================================================
!
!      use module_principal , only : kount0, iteration3d, rhp_t, rho, mask_t, &
!			    dz_t, iteration2d_max_now, hz_w, iteration2d,    &
!			    imax, jmax, kmax
!     use module_parallele                ! #MPI#
      use module_nh                       ! #NH#
      use module_nbq                      ! #NBQ#
      implicit none
      integer :: icall, i,j ,k, k1,indm_d
      integer Istr,Iend,Jstr,Jend

# include "param_F90.h"
# include "scalars_F90.h"
# include "work.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "nbq.h"

      integer :: ncp
      real    :: dist_d,rho_sum

      double precision :: t1_d,t2_d

!#include "compute_tile_bounds.h"

      if (icall.eq.20) then
# ifdef NBQ_MASS
!
!**********************************************************************
!  Transfer density field to i,j,k array 
!  and time filter, ready for external mode
!**********************************************************************
#ifdef NBQ_ZETAW
        rhobar_nbq(istrq_nh-1:iendq_nh+1,jstrq_nh-1:jendq_nh+1,knew2)=0.
#else
        rhobar_nbq(istrq_nh-1:iendq_nh+1,jstrq_nh-1:jendq_nh+1,knew)=0.
#endif 

         do k=1,N
           do j=jstrq_nh-1,jendq_nh+1
             do i=istrq_nh-1,iendq_nh+1
#ifdef NBQ_ZETAW
               rhobar_nbq(i,j,knew2)= rhobar_nbq(i,j,knew2) + rho_nbq(i,j,k) &  !XXX1
                             +rho(i,j,k)/rho0*hzr(i,j,k)
#else
               rhobar_nbq(i,j,knew)= rhobar_nbq(i,j,knew) + rho_nbq(i,j,k)  !XXX1
#endif
               rho_nbq_ext(i,j,k)  = 1.+rho_nbq(i,j,k)/Hzr(i,j,k) &
                                       +rho(i,j,k)/rho0
             enddo  
           enddo  
         enddo

!         rho_sum=0.

         do j=jstrq_nh-1,jendq_nh+1
           do i=istrq_nh-1,iendq_nh+1
#  ifdef NBQ_ZETAW
             rhobar_nbq(i,j,knew2) = 1.+(rhobar_nbq(i,j,knew2)) & 
                / (z_w(i,j,N)-z_w(i,j,0))
#  else
             rhobar_nbq(i,j,knew) = 1.+(rhobar_nbq_int(i,j)+rhobar_nbq(i,j,knew)) & 
                / (zw_half_nbq(i,j,N)-zw_half_nbq(i,j,0))
#  endif
           enddo
         enddo
 
!          rhobar_nbq=1.
!          rho_nbq_ext=1.       

!
!      

!c LAURENT: the two next exchanges should not be needed
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
         call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,rhobar_nbq(START_2D_ARRAY,knew2))   ! TBD
         call exchange_r3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh &
                                ,rho_nbq_ext(START_2D_ARRAY,1))
# endif
 
  
# ifdef RVTK_DEBUG
!      call check_tab3d(rho_nbq_ext(:,:,1:N),'rho_nbq_ext (density_nbq)','r')
!      call check_tab2d(rhobar_nbq(:,:,knew),'rhobar_nbq (density_nbq)','r')
# endif    

#endif

      elseif (icall.eq.6) then
      elseif (icall.eq.60) then
      elseif (icall.eq.61) then
      elseif (icall.eq.7) then

# if defined ACOUSTIC  && defined NBQ_IJK
      elseif (icall.eq.10) then
!
!*******************************************************************
!......Acoustic waves: Initialization
!*******************************************************************
          period_exp = 0.025/2.
          for_a_exp  = 2.5
          amp_exp    = 1.e-2
          hmax_exp   = 128.
          dg_exp     = 2.
          time_nbq = 0. 

      elseif (icall.eq.11) then
!
!*******************************************************************
!......Acoustic waves: Forcing
!*******************************************************************
!
        time_nbq = time_nbq + dtfast

        do k=1,N
      !  do j=Jstr-1,Jend+1
      !  do i=Istr-1,Iend+1
         do j=jstrq_nh-1,jendq_nh+1
           do i=istrq_nh-1,iendq_nh+1

             dist_d=sqrt((xr(i,j)-xl/2.)**2+(0.*(yr(i,j)-el/2.))**2    & 
                               +(abs(z_w(i,j,k))-hmax_exp/2.)**2)

             rho_nbq(i,j,k)=rho_nbq(i,j,k)+amp_exp                       &
                                          *sin(2*acos(-1.)*time_nbq/period_exp)   &  
                                          *exp(-dist_d**2/for_a_exp**2)  &
#  ifdef NBQ_MASS
                        *Hzr(i,j,k)
#  else
                        *Hz(i,j,k)
#  endif

            enddo
          enddo
        enddo

# endif /* ACOUSTIC */

      endif  ! icall

        return
        end
#else
        subroutine densityijk_nbq_empty
        return
        end
#endif
