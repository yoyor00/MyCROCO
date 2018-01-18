 
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
!           dum_s=0.
            do k=1,N
              Hz(i,j,k)=Hz_bak2(i,j,k)
     &        - dtfast*(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k))

!             Hz(i,j,k)=
!     &        cff4*(Hz(i,j,k)
!     &        - dtfast*(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k)))
!     &        +cff5*Hz_t(i,j,k,kstp2)
!     &        +cff6*Hz_t(i,j,k,kbak2)
!     &        +cff7*Hz_t(i,j,k,kold2)
!             Hz_t(i,j,k,knew2)=Hz(i,j,k)
!             dum_s=dum_s+Hz(i,j,k)

              Hzr(i,j,k) =(Hz(i,j,k)-rho_nbq(i,j,k))
     &                   /(1.+rho(i,j,k)/rho0)

!            enddo
!            do k=1,N

              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
              z_r(i,j,k)=(z_w(i,j,k)+z_w(i,j,k-1))/2.D0
            enddo

        enddo
      enddo

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI      
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew2))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_w(START_2D_ARRAY,0))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        z_r(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hzr(START_2D_ARRAY,1))
#endif

!       if ((iic.eq.1.and.iif==1)
!     &     NSTEP_GRID
!     &     .or.iif.eq.nfast) then
!        flag_grid=1
!        call set_depth_tile(Istr,Iend,Jstr,Jend
!     &   ,resetfromrhobar
!     &   ) 
!        endif

#include "step2d_grid_ext.h"

!       call grid_coef_nh(
!    &   Istr,Iend,Jstr,Jend,
!    &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
!    &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
!    &   Hzu_half_qdmu, Hzv_half_qdmv
!    &   )

