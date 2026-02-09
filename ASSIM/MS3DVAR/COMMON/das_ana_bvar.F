#include "cppdefs.h"
#ifdef DAS_ANA_BVAR

      subroutine das_ana_bvar (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_ana_bvar_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_ana_bvar_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
!
      real btk,bsk,buk,bvk,bz,vr0,expz(NDAS)
!
# include "compute_auxiliary_bounds.h"
!
# ifdef EW_PERIODIC
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE  Istr,IendR
# endif

# ifdef NS_PERIODIC
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE  Jstr,JendR
# endif
!
!  forecast error variance
!------------------------------------------------------------------
!
!      btk=2.5
!      bsk=3.5
      btk=5.0
      bsk=0.15
      buk=0.15
      bvk=0.15
      bz=0.05
      vr0=280.0    ! e-fold depth
      do k=1,NDAS
        expz(k)=exp(z_das(k)/vr0)
      enddo
!
! zeta
!
      do j=JR_RANGE
        do i=IR_RANGE
          bz_das(i,j)=(bz+0.002*bz)
        enddo
      enddo
!
! temp and salt
!
      do k=1,NDAS
        do j=JR_RANGE
          do i=IR_RANGE
            bt_das(i,j,k,itemp)=(btk*expz(k)+0.002*btk)
            bt_das(i,j,k,isalt)=(bsk*expz(k)+0.002*bsk)
          enddo
        enddo
      enddo
      do k=1,NDAS
        do j=JR_RANGE
          do i=IU_RANGE
            bu_das(i,j,k)=(buk*expz(k)+0.002*buk)
          enddo
        enddo
      enddo
      do k=1,NDAS
        do j=JV_RANGE
          do i=IR_RANGE
            bv_das(i,j,k)=(bvk*expz(k)+0.002*bvk)
          enddo
        enddo
      enddo

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        bz_das(START_2D_ARRAY))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                         bu_das(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                         bv_das(START_2D_ARRAY,1))
      do itrc=1,NT
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                   bt_das(START_2D_ARRAY,1,itrc))
      enddo
#endif
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end

#else
      subroutine das_ana_bvar_empty
      return
      end
#endif /* DAS_ANA_BVAR */
