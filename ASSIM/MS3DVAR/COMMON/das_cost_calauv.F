#include "cppdefs.h"

      subroutine das_cost_calauv (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_cost_calauv_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_cost_calauv_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, NSUB 
#include "param.h"
# include "das_param.h"
      real cost_tile_my
      real cost_my_t,cost_my_s
      real cft, cfs, kdas
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
# include "das_innov.h"
!
#ifdef MPI
      include 'mpif.h'
      integer size, step, status(MPI_STATUS_SIZE), ierr
      real buff
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

!      cal_ot=1.0/(0.4*0.4)
!      cal_os=1.0/(0.25*0.25)
!
! tile summation
!
      cost_my_t=0.
      cost_my_s=0.
      do k=1,max_prf_cal
        kdas=NDAS-k+1
        do j=JR_RANGE
          do i=IR_RANGE
            cft=t_s(i,j,kdas,itemp)-t_cal(i,j,k)  
            cfs=t_s(i,j,kdas,isalt)-s_cal(i,j,k)  
            cost_my_t=cost_my_t + cft*cft
     &                 *cal_mask(i,j,k)*cal_ot
            cost_my_s=cost_my_s + cfs*cfs
     &                 *cal_mask(i,j,k)*cal_os
          enddo
        enddo
      enddo
!
      cost_tile_my=cost_my_t+cost_my_s
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
!
      if (SINGLE_TILE_MODE) then
        NSUB=1
      else
        NSUB=NSUB_X*NSUB_E
      endif
!
! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out. 
!
! NOTE: costf has the part of the background
!
C$OMP CRITICAL (cost_auv_cr)
      costf=costf + 0.5 * cost_tile_my            ! summation among
                                                  ! the threads
      tile_count=tile_count+1         ! This counter identifies
      if (tile_count.eq.NSUB) then    ! the last thread, whoever
        tile_count=0                  ! it is, not always master.
#ifdef MPI
        if (NNODES.gt.1) then         ! Perform global summation 
          size=NNODES                 ! among MPI processes
   1      step=(size+1)/2 
          if (mynode.ge.step .and. mynode.lt.size) then
            buff=costf         ! This is MPI_Reduce
            call MPI_Send (buff,  1, MPI_DOUBLE_PRECISION,
     &             mynode-step, 17, MPI_COMM_WORLD,      ierr)
            elseif (mynode .lt. size-step) then
              call MPI_Recv (buff,  1, MPI_DOUBLE_PRECISION,
     &             mynode+step, 17, MPI_COMM_WORLD, status, ierr)
              costf=costf+buff
            endif
            size=step
            if (size.gt.1) goto 1
        endif
        if (mynode.eq.0) then
#endif
!
! Raise may_day_flag to stop computations in the case of blowing up.
! [Criterion for blowing up here is the numerical overflow, so that
! avgkp is 'INF' or 'NAN' (any mix of lover and uppercase letters),
! therefore it is sufficient to check for the presence of letter 'N'.
!
# ifdef MPI
        endif    ! <-- mynode.eq.0
# endif
      endif
C$OMP END CRITICAL (cost_auv_cr)
      return
      end

