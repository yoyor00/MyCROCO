#include "cppdefs.h"

      subroutine das_cost_reg (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
#include "param.h"
#include "das_param.h"
#include "das_private_scratch.h"
#include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_cost_reg_tile (Istr,Iend,Jstr,Jend, 
     &          A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd))
      return
      end

      subroutine das_cost_reg_tile (Istr,Iend,Jstr,Jend, 
     &                                cft2d,cfs2d,cfz2d)
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k, NSUB 
#include "param.h"
# include "das_param.h"
      real cft2d(PRIVATE_2D_SCRATCH_ARRAY),
     &     cfs2d(PRIVATE_2D_SCRATCH_ARRAY),
     &     cfz2d(PRIVATE_2D_SCRATCH_ARRAY)
      real cost_tile_my, costt_my, costs_my, costz_my,
     &                   costu_my,costv_my, nozero
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

#if defined DAS_COUNT_PSICHI
      integer :: count_non_null_psi
      integer :: count_non_null_chi
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
#if defined DAS_COUNT_PSICHI
      count_non_null_psi=0
      count_non_null_chi=0
#endif

      nozero=0.0001
      do j=JstrR,JendR
        do i=IstrR,IendR
          cft2d(i,j)=0.
          cfs2d(i,j)=0.
          cfz2d(i,j)=0.
        enddo
        do k=NDAS,1,-1
          do i=IstrR,IendR
            cft2d(i,j)=cft2d(i,j)+t_s(i,j,k,itemp)
     &           * t_s(i,j,k,itemp) 
     &           * regp/(bt_das(i,j,k,itemp)*bt_das(i,j,k,itemp) + nozero)
            cfs2d(i,j)=cfs2d(i,j)+t_s(i,j,k,isalt)
     &           * t_s(i,j,k,isalt) 
     &           * regp/(bt_das(i,j,k,isalt)*bt_das(i,j,k,isalt) + nozero)
          enddo
        enddo
#if !defined DAS_HYDRO_STRONG
        do i=IstrR,IendR
          cfz2d(i,j)=cfz2d(i,j)+zeta_s(i,j)
     &                   * zeta_s(i,j) 
     &                   * regp/(bz_das(i,j)*bz_das(i,j) + nozero)
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
# if defined DAS_PSICHI_COSTREG
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
            cft2d(i,j)=cft2d(i,j)+psi_s(i,j,k)                                                                               
     &          * psi_s(i,j,k) 
     &          * regp/(bpsi_das(i,j,k)*bpsi_das(i,j,k)+ nozero)
#if defined DAS_COUNT_PSICHI
          if (abs(psi_s(i,j,k)) .gt. 0. ) then
            count_non_null_psi=count_non_null_psi+1
          endif
#endif
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          cfs2d(i,j)=0.
        enddo
        do k=NDAS,1,-1
          do i=IstrR,IendR
            cfs2d(i,j)=cfs2d(i,j)+chi_s(i,j,k)                                                                               
     &          * chi_s(i,j,k) 
     &          * regp/(bchi_das(i,j,k)*bchi_das(i,j,k)+ nozero)
#if defined DAS_COUNT_PSICHI
          if (abs(chi_s(i,j,k)) .gt. 0. ) then
            count_non_null_chi=count_non_null_chi+1
          endif
#endif
          enddo
        enddo
      enddo
#if defined DAS_COUNT_PSICHI
      if ( count_non_null_psi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'Non null sig*C*psi_s = ',count_non_null_psi
      endif 
      if ( count_non_null_chi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'Non null sig*C*chi_s = ',count_non_null_chi
      endif 
#endif
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
# endif
#endif
       cost_tile_my=    costt_my
     &               + costs_my
#if !defined DAS_HYDRO_STRONG
     &               + costz_my
#endif
#if !defined DAS_GEOS_STRONG
# if defined DAS_PSICHI_COSTREG
     &               + costu_my+costv_my
# endif
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
C$OMP CRITICAL (costreg_cr)
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
C$OMP END CRITICAL (costreg_cr)
      return
      end

