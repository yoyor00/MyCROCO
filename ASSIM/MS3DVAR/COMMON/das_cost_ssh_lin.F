#include "cppdefs.h"
#if defined DAS_JASONSSH || defined DAS_TPSSH

      subroutine das_cost_ssh_lin (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_cost_ssh_lin_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_cost_ssh_lin_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j, NSUB 
#include "param.h"
# include "das_param.h"
      real cost_tile_my, cost_my_js1, cost_my_tp
      real cft, cfm
      real cost_tile_my_adj, cost_my_js1_adj, cost_my_tp_adj
      real cft_adj, cfm_adj
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

!
! tile summation
!
# if defined DAS_JASONSSH
      cost_my_js1_adj=0.
      cost_my_js1=0.
# endif
# if defined DAS_TPSSH
      cost_my_tp_adj=0.
      cost_my_tp=0.
# endif
      do j=JR_RANGE
        do i=IR_RANGE
# if defined DAS_JASONSSH
          cft_adj=zeta_s_adj(i,j)
          cft=zeta_s(i,j)-ssh_js1(i,j)  
          cost_my_js1_adj=cost_my_js1_adj + 2.*cft*cft_adj
     &                     *js1_mask(i,j)*js1_oin(i,j)
          cost_my_js1=cost_my_js1 + cft*cft
     &                     *js1_mask(i,j)*js1_oin(i,j)
# endif
# if defined DAS_TPSSH
          cfm_adj=zeta_s_adj(i,j)
          cfm=zeta_s(i,j)-ssh_tp(i,j)
          cost_my_tp_adj=cost_my_tp_adj + 2.*cfm*cfm_adj
     &                    *tp_mask(i,j)*tp_oin(i,j)
          cost_my_tp=cost_my_tp + cfm*cfm
     &                    *tp_mask(i,j)*tp_oin(i,j)
# endif
        enddo
      enddo
!
      cost_tile_my_adj=
# if defined DAS_JASONSSH
     &               + cost_my_js1_adj
# endif
# if defined DAS_TPSSH
     &               + cost_my_tp_adj
# endif
      cost_tile_my=
# if defined DAS_JASONSSH
     &               + cost_my_js1
# endif
# if defined DAS_TPSSH
     &               + cost_my_tp
# endif
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
C$OMP CRITICAL (lin_csst_cr)
      costf_adj=costf_adj + 0.5 * cost_tile_my_adj            ! summation among
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
C$OMP END CRITICAL (lin_csst_cr)
      return
      end

#else
      subroutine das_cost_ssh_lin_empty
      return
      end
#endif /* DAS_JASONSSH DAS_TPSSH */
