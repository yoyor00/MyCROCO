#include "cppdefs.h"

      subroutine das_costb (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
#include "param.h"
#include "das_param.h"
#include "das_private_scratch.h"
#include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_costb_tile (Istr,Iend,Jstr,Jend, 
     &          A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd))
      return
      end

      subroutine das_costb_tile (Istr,Iend,Jstr,Jend, 
     &                                cft2d,cfs2d,cfz2d)
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k, NSUB 
#include "param.h"
# include "das_param.h"
      real cft2d(PRIVATE_2D_SCRATCH_ARRAY),
     &     cfs2d(PRIVATE_2D_SCRATCH_ARRAY),
     &     cfz2d(PRIVATE_2D_SCRATCH_ARRAY)
      real cost_tile_my, costt_my, costs_my, costz_my,
     &                   costu_my,costv_my
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
! Compute and report volume averaged kinetic, potential and total
!-----------------------------------------------------------------
! energy densities for either two- (shallow water) or three-
! dimensional versions of the model.
!
! At first, compute kinetic and potential energies, as well as total
! volume within the tile [subdomain of indices (Istr:Iend,Jstr:Jend)]
! by individual threads. In the case of three dimensions also perform
! verical summation at this stage.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          cft2d(i,j)=0.
          cfs2d(i,j)=0.
          cfz2d(i,j)=0.
        enddo
        do k=NDAS,1,-1
          do i=IstrR,IendR
            cft2d(i,j)=cft2d(i,j)+t_das(i,j,k,itemp)
     &                   * t_w(i,j,k,itemp) 
            cfs2d(i,j)=cfs2d(i,j)+t_das(i,j,k,isalt)
     &                   * t_w(i,j,k,isalt) 
          enddo
        enddo
#if !defined DAS_HYDRO_STRONG
        do i=IstrR,IendR
          cfz2d(i,j)=zeta_das(i,j)
     &                   * zeta_das(i,j) 
        enddo
#endif
      enddo
!
! tile summation
!
      costt_my=0.
      costs_my=0.
#if !defined DAS_HYDRO_STRONG
      costz_my=0.
#endif
      do j=JstrR,JendR
        do i=IstrR,IendR
          costt_my=costt_my+cft2d(i,j)
          costs_my=costs_my+cfs2d(i,j)
#if !defined DAS_HYDRO_STRONG
          costz_my=costz_my+cfz2d(i,j)
#endif
        enddo
      enddo
!
#if !defined DAS_GEOS_STRONG
!
!
! psi and chi
!
      do j=Jstr,JendR
        do i=Istr,IendR
          cft2d(i,j)=0.
        enddo
        do k=NDAS,1,-1
          do i=Istr,IendR
            cft2d(i,j)=cft2d(i,j)+psi_das(i,j,k)                                                                               
     &                   * psi_das(i,j,k) 
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          cfs2d(i,j)=0.
        enddo
        do k=NDAS,1,-1
          do i=IstrR,IendR
            cfs2d(i,j)=cfs2d(i,j)+chi_das(i,j,k)                                                                               
     &                   * chi_das(i,j,k) 
          enddo
        enddo
      enddo
      costu_my=0.
      costv_my=0.
      do j=Jstr,JendR
        do i=Istr,IendR
          costu_my=costu_my+cft2d(i,j)
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          costv_my=costv_my+cfs2d(i,j)
        enddo
      enddo
#endif
       cost_tile_my=    costt_my
     &               + costs_my
#if !defined DAS_HYDRO_STRONG
     &               + costz_my
#endif
#if !defined DAS_GEOS_STRONG
     &               + costu_my+costv_my
#endif
!
      if (SINGLE_TILE_MODE) then
        NSUB=1
      else
        NSUB=NSUB_X*NSUB_E
      endif
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
C$OMP CRITICAL (costb_cr)
!      if (tile_count.eq.0) then
!        costf=0.0
!      endif                                       ! Perform global
      costf=costf + 0.5*cost_tile_my                  ! summation among
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
#ifdef MPI
        endif    ! <-- mynode.eq.0
#endif
      endif
C$OMP END CRITICAL (costb_cr)
      return
      end

