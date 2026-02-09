#include "cppdefs.h"

      subroutine das_matvec (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_matvec_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_matvec_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k, ie,je,ke
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
!
#ifdef MPI
      write (*,'(6x,A)')
     &     'The routine das_matvec is not ready for MPI'
      stop
#endif

# include "compute_auxiliary_bounds.h"
!
! compute the entire model domain
! but the lower triangular matrices !?
!
! see das_param.h Local_len correlation set to zero
! below this Local_len=24
! BLXD das_matvec/das_adj_matvec introduces bdy spurious values 
!      in the increment through the asymetric interplay of L and L^T 
! why calculating the increment at bdy data unless 
! if ms3dvar would support MPI subdomain decomposition
! => temporary fix to remove bdy L/L^T artefacts : 
! removing the calculation at indices 0/0 Lm+1/Mm+1

# ifdef EW_PERIODIC
#  define IR_ENTIRE max(1,I-Local_len),I
#  define IU_ENTIRE max(1,I-Local_len),I
# else
#  ifdef MPI
#   define IR_ENTIRE max(0,I-Local_len),I
#   define IU_ENTIRE max(1,I-Local_len),I
#  else 
#   define IR_ENTIRE max(1,I-Local_len),I
#   define IU_ENTIRE max(1,I-Local_len),I
#  endif
# endif

# ifdef NS_PERIODIC
#  define JR_ENTIRE max(1,J-Local_len),J
#  define JV_ENTIRE max(1,J-Local_len),J
# else
#  ifdef MPI
#   define JR_ENTIRE max(0,J-Local_len),J
#   define JV_ENTIRE max(1,J-Local_len),J
#  else 
#   define JR_ENTIRE max(1,J-Local_len),J
#   define JV_ENTIRE max(1,J-Local_len),J
#  endif
# endif
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
!
#if !defined DAS_BVAR_CORR
!
            t_s(i,j,k,itemp)=t_w(i,j,k,itemp)
            t_s(i,j,k,isalt)=t_w(i,j,k,isalt)
# if !defined DAS_GEOS_STRONG
            chi_s(i,j,k)=chi_das(i,j,k)
# endif
#else
            t_s(i,j,k,itemp)=0.0
            t_s(i,j,k,isalt)=0.0
# if !defined DAS_GEOS_STRONG
            chi_s(i,j,k)=0.0
# endif
!
            do ke=1,k
              do je=JR_ENTIRE
                do ie=IR_ENTIRE
                  t_s(i,j,k,itemp)=t_s(i,j,k,itemp)
     &                  +ctX_das(i,ie,itemp)*ctE_das(j,je,itemp)
     &                   *ctZ_das(k,ke,itemp)*t_w(ie,je,ke,itemp)
                  t_s(i,j,k,isalt)=t_s(i,j,k,isalt)
     &                  +ctX_das(i,ie,isalt)*ctE_das(j,je,isalt)
     &                   *ctZ_das(k,ke,isalt)*t_w(ie,je,ke,isalt)
# if !defined DAS_GEOS_STRONG
                  chi_s(i,j,k)=chi_s(i,j,k)
     &                  +cchiX_das(i,ie)*cchiE_das(j,je)
     &                   *cchiZ_das(k,ke)*chi_das(ie,je,ke)
# endif
                enddo
              enddo
            enddo
#endif
!             
          enddo
        enddo
      enddo
!
#if !defined DAS_GEOS_STRONG
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,IendR
!
# if !defined DAS_BVAR_CORR
            psi_s(i,j,k)=psi_das(i,j,k)
# else
            psi_s(i,j,k)=0.0
            do ke=1,k
              do je=JV_ENTIRE
                do ie=IU_ENTIRE
                  psi_s(i,j,k)=psi_s(i,j,k)
     &                  +cpsiX_das(i,ie)*cpsiE_das(j,je)
     &                   *cpsiZ_das(k,ke)*psi_das(ie,je,ke)
                enddo
              enddo
            enddo
# endif
          enddo
        enddo
      enddo
#endif

!
! zeta
!
#if !defined DAS_HYDRO_STRONG
      do j=JstrR,JendR
        do i=IstrR,IendR
!
# if !defined DAS_BVAR_CORR
!
          zeta_s(i,j)=zeta_das(i,j)
# else
          zeta_s(i,j)=0.0
!
          do je=JR_ENTIRE
            do ie=IR_ENTIRE
              zeta_s(i,j)=zeta_s(i,j)
     &                +czX_das(i,ie)*czE_das(j,je)
     &                          *zeta_das(ie,je)
            enddo
          enddo
# endif
        enddo
      enddo
#endif
!
# undef IR_ENTIRE
# undef IU_ENTIRE
# undef JR_ENTIRE
# undef JV_ENTIRE
!
      return
      end

