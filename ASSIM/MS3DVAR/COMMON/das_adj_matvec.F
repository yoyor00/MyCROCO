#include "cppdefs.h"

      subroutine das_adj_matvec (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_matvec_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_matvec_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k,ie,je,ke
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
!
#ifdef MPI
      include 'mpif.h'
      integer size, step, status(MPI_STATUS_SIZE), ierr
      real buff
#endif

# include "compute_auxiliary_bounds.h"
!
! compute the entire model domain
! but the lower triangular matrices 
!
! BLXD das_matvec/das_adj_matvec introduces bdy spurious values 
!      in the increment through the asymetric interplay of L and L^T 
! why calculating the increment at bdy data unless 
! if ms3dvar would support MPI subdomain decomposition
! => temporary fix to remove bdy L/L^T artefacts : 
! removing the calculation at indices 0/0 Lm+1/Mm+1

# ifdef EW_PERIODIC
#  define IR_ENTIRE I,min(I+Local_len,Lm)
#  define IU_ENTIRE I,min(I+Local_len,Lm)
# else
#  ifdef MPI
#   define IR_ENTIRE I,min(I+Local_len,Lm+1)
#   define IU_ENTIRE I,min(I+Local_len,Lm+1)
#  else 
#   define IR_ENTIRE I,min(I+Local_len,Lm)
#   define IU_ENTIRE I,min(I+Local_len,Lm)
#  endif
# endif

# ifdef NS_PERIODIC
#  define JR_ENTIRE J,min(J+Local_len,Mm)
#  define JV_ENTIRE J,min(J+Local_len,Mm)
# else
#  ifdef MPI
#   define JR_ENTIRE J,min(J+Local_len,Mm+1)
#   define JV_ENTIRE J,min(J+Local_len,Mm+1)
#  else 
#   define JR_ENTIRE J,min(J+Local_len,Mm)
#   define JV_ENTIRE J,min(J+Local_len,Mm)
#  endif

# endif
!
      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
!
#if !defined DAS_BVAR_CORR
!
            t_w_adj(i,j,k,itemp)=t_w_adj(i,j,k,itemp)
     &           + t_s_adj(i,j,k,itemp)
            t_w_adj(i,j,k,isalt)=t_w_adj(i,j,k,isalt)
     &           + t_s_adj(i,j,k,isalt)
            chi_adj(i,j,k)=chi_adj(i,j,k)
     &            + chi_s_adj(i,j,k)
! Missing line continuation undef DAS_BVAR_CORR never tested ?
#else
!
            do ke=k,NDAS
              do je=JR_ENTIRE
                do ie=IR_ENTIRE
                  t_w_adj(i,j,k,itemp)=t_w_adj(i,j,k,itemp)
     &                  +ctX_das(i,ie,itemp)*ctE_das(j,je,itemp)
     &                   *ctZ_das(k,ke,itemp)*t_s_adj(ie,je,ke,itemp)
                  t_w_adj(i,j,k,isalt)=t_w_adj(i,j,k,isalt)
     &                  +ctX_das(i,ie,isalt)*ctE_das(j,je,isalt)
     &                   *ctZ_das(k,ke,isalt)*t_s_adj(ie,je,ke,isalt)
                  chi_adj(i,j,k)=chi_adj(i,j,k)
     &                  +cchiX_das(i,ie)*cchiE_das(j,je)
     &                   *cchiZ_das(k,ke)*chi_s_adj(ie,je,ke)
                enddo
              enddo
            enddo
#endif
!
          enddo
        enddo
      enddo
!
!psi
!
      do k=1,NDAS
        do j=Jstr,JendR
          do i=Istr,IendR
!
#if !defined DAS_BVAR_CORR
            psi_adj(i,j,k)=psi_adj(i,j,k)
     &                 + psi_s_adj(i,j,k)
! Missing line continuation undef DAS_BVAR_CORR never tested ?
#else
!
            do ke=k,NDAS
              do je=JV_ENTIRE                         
                do ie=IU_ENTIRE                      
                  psi_adj(i,j,k)=psi_adj(i,j,k)     
     &                  +cpsiX_das(i,ie)*cpsiE_das(j,je) 
     &                   *cpsiZ_das(k,ke)*psi_s_adj(ie,je,ke) 
                enddo                           
              enddo                                     
            enddo                                      
#endif
!                               
          enddo                   
        enddo                    
      enddo
!
! zeta
!
#if !defined DAS_HYDRO_STRONG
      do j=JstrR,JendR
        do i=IstrR,IendR
!
# if !defined DAS_BVAR_CORR
!
          zeta_adj(i,j)=zeta_adj(i,j) + zeta_s_adj(i,j)
# else
!
          do je=JR_ENTIRE
            do ie=IR_ENTIRE
              zeta_adj(i,j)=zeta_adj(i,j)
     &                +czX_das(i,ie)*czE_das(j,je)
     &                 *zeta_s_adj(ie,je)
            enddo
          enddo
# endif
!
        enddo
      enddo
#endif
# undef IR_ENTIRE
# undef IU_ENTIRE
# undef JR_ENTIRE
# undef JV_ENTIRE
!
      return
      end

