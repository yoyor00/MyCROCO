#include "cppdefs.h"

      subroutine das_increment (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_increment_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_increment_tile (Istr,Iend,Jstr,Jend) 
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
! compute the entire model domain
! but the lower triangular matrices 
! 
!
# ifdef EW_PERIODIC
#  define IR_ENTIRE 1,I
#  define IU_ENTIRE 1,I
# else
#  define IR_ENTIRE 0,I
#  define IU_ENTIRE 1,I
# endif

# ifdef NS_PERIODIC
#  define JR_ENTIRE 1,J
#  define JV_ENTIRE 1,J
# else
#  define JR_ENTIRE 0,J
#  define JV_ENTIRE 1,J
# endif
!
      do k=1,NDAS
        do j=JR_RANGE
          do i=IR_RANGE
!
#if !defined DAS_BVAR_CORR
            t_s(i,j,k,itemp)=t_das(i,j,k,itemp) 
            t_s(i,j,k,isalt)=t_das(i,j,k,isalt) 
#else
            t_s(i,j,k,itemp)=0.0
            t_s(i,j,k,isalt)=0.0
!
            do ke=1,k
              do je=JR_ENTIRE
                do ie=IR_ENTIRE
                  t_s(i,j,k,itemp)=t_s(i,j,k,itemp) 
     &                    +ctX_das(i,ie,itemp)*ctE_das(j,je,itemp)
     &                     *ctZ_das(k,ke,itemp)*t_das(ie,je,ke,itemp)
                  t_s(i,j,k,isalt)=t_s(i,j,k,isalt) 
     &                    +ctX_das(i,ie,isalt)*ctE_das(j,je,isalt)
     &                     *ctZ_das(k,ke,isalt)*t_das(ie,je,ke,isalt)
                enddo
              enddo
            enddo
#endif
!
            t_s(i,j,k,itemp)=t_s(i,j,k,itemp) 
     &                       * rmask_das(i,j,k) * bt_das(i,j,k,itemp)
            t_s(i,j,k,isalt)=t_s(i,j,k,isalt) 
     &                       * rmask_das(i,j,k) * bt_das(i,j,k,isalt)
          enddo
        enddo
      enddo
!
      return
      end

