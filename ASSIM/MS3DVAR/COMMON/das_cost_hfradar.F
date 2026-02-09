#include "cppdefs.h"

      subroutine das_cost_hfradar (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_cost_hfradar_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_cost_hfradar_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j, NSUB 
#include "param.h"
# include "das_param.h"
      real cost_tile_my, cost_my_hfu,cost_my_hfv
      real cft
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
!
! tile summation
!
! u
!
      cost_my_hfu=0.0
      do j=Jstr,Jend
        do i=Istr,IendR
          cft=u_s(i,j,NDAS)-u_hf(i,j)  
          cost_my_hfu=cost_my_hfu + cft*cft
     &                 *hf_umask(i,j)*hf_ouv
        enddo
      enddo
!
! v
!
      cost_my_hfv=0.0
      do j=Jstr,JendR
        do i=Istr,Iend
          cft=v_s(i,j,NDAS)-v_hf(i,j)
          cost_my_hfv=cost_my_hfv + cft*cft
     &                 *hf_vmask(i,j)*hf_ouv
        enddo
      enddo
!
! 6KM
!
!
! u
!
      do j=Jstr,Jend
        do i=Istr,IendR
          cft=u_s(i,j,NDAS)-u_hf6(i,j)
          cost_my_hfu=cost_my_hfu + cft*cft
     &                 *hf_umask6(i,j)*hf_uierr6(i,j)
        enddo
      enddo
!
! v
!
      do j=Jstr,JendR
        do i=Istr,Iend
          cft=v_s(i,j,NDAS)-v_hf6(i,j)
          cost_my_hfv=cost_my_hfv + cft*cft
     &                 *hf_vmask6(i,j)*hf_vierr6(i,j)
        enddo
      enddo
!
      cost_tile_my=cost_my_hfu + cost_my_hfv
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
C$OMP CRITICAL (cost_fl_cr)
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
C$OMP END CRITICAL (cost_fl_cr)
      return
      end

