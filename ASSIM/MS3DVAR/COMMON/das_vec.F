#include "cppdefs.h"
      subroutine das_ca2vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! control variable arrays to vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      implicit none
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "das_innov.h"
#include "das_lbgfs.h"
      integer i,j,k,l0

!
      l0=0
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            xmin(l0)=t_das(i,j,k,itemp)
          enddo 
        enddo 
      enddo 
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            xmin(l0)=t_das(i,j,k,isalt)
          enddo 
        enddo 
      enddo 
#if !defined DAS_GEOS_STRONG
      do k=1,NDAS
        do j=1,Mm+1
          do i=1,Lm+1
            l0=l0+1
            xmin(l0)=psi_das(i,j,k)
          enddo
        enddo
      enddo
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            xmin(l0)=chi_das(i,j,k)
          enddo
        enddo
      enddo
#endif
#if !defined DAS_HYDRO_STRONG
      do j=0,Mm+1
        do i=0,Lm+1
          l0=l0+1
          xmin(l0)=zeta_das(i,j)
        enddo
      enddo
#endif
      if ( l0 .ne . NDIM ) then
        write(*, '(6x,A)') 
     &     '!!! das_a2vec: ERROR in dimenision'
        stop
      endif 
!      
      return
      end
!
      subroutine das_ga2vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! control variable gradient arrays to vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      implicit none
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "das_innov.h"
#include "das_lbgfs.h"
      integer i,j,k,l0
#if defined DAS_COUNT_PSICHI
      integer :: count_non_null_psi
      integer :: count_non_null_chi
      count_non_null_psi=0
      count_non_null_chi=0
#endif
!
      l0=0
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            xgrd(l0)=t_adj(i,j,k,itemp)
          enddo 
        enddo 
      enddo 
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            xgrd(l0)=t_adj(i,j,k,isalt)
          enddo 
        enddo 
      enddo 
#if !defined DAS_GEOS_STRONG
      do k=1,NDAS
        do j=1,Mm+1
          do i=1,Lm+1
            l0=l0+1
            xgrd(l0)=psi_adj(i,j,k)
#if defined DAS_COUNT_PSICHI
          if (abs(psi_adj(i,j,k)) .gt. 0. ) then
            count_non_null_psi=count_non_null_psi+1
          endif
#endif
          enddo
        enddo
      enddo
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            xgrd(l0)=chi_adj(i,j,k)
#if defined DAS_COUNT_PSICHI
          if (abs(chi_adj(i,j,k)) .gt. 0. ) then
            count_non_null_chi=count_non_null_chi+1
          endif
#endif
          enddo
        enddo
      enddo
#endif
#if !defined DAS_HYDRO_STRONG
      do j=0,Mm+1
        do i=0,Lm+1
          l0=l0+1
          xgrd(l0)=zeta_adj(i,j)
        enddo
      enddo
#endif
#if defined DAS_COUNT_PSICHI
      if ( count_non_null_psi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null psi_adj = ',count_non_null_psi
      endif 
      if ( count_non_null_chi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null chi_adj = ',count_non_null_chi
      endif 
#endif
      if ( l0 .ne . NDIM ) then
        write(*, '(6x,A)') 
     &     '!!! das_ga2vec: ERROR in dimenision'
        stop
      endif 
!      
      return
      end

      subroutine das_vec2ca
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vec to control variable arrays 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      implicit none
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "das_innov.h"
#include "das_lbgfs.h"
      integer i,j,k,l0
#if defined DAS_COUNT_PSICHI
      integer :: count_non_null_psi
      integer :: count_non_null_chi
      count_non_null_psi=0
      count_non_null_chi=0
#endif
!
      l0=0
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            t_das(i,j,k,itemp)=xmin(l0)
          enddo 
        enddo 
      enddo 
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            t_das(i,j,k,isalt)=xmin(l0)
          enddo 
        enddo 
      enddo 
#if !defined DAS_GEOS_STRONG
      do k=1,NDAS
        do j=1,Mm+1
          do i=1,Lm+1
            l0=l0+1
            psi_das(i,j,k)=xmin(l0)
#if defined DAS_COUNT_PSICHI
          if (abs(psi_das(i,j,k)) .gt. 0. ) then
            count_non_null_psi=count_non_null_psi+1
          endif
#endif
          enddo
        enddo
      enddo
      do k=1,NDAS
        do j=0,Mm+1
          do i=0,Lm+1
            l0=l0+1
            chi_das(i,j,k)=xmin(l0)
#if defined DAS_COUNT_PSICHI
          if (abs(chi_das(i,j,k)) .gt. 0. ) then
            count_non_null_chi=count_non_null_chi+1
          endif
#endif
          enddo
        enddo
      enddo
#endif
#if !defined DAS_HYDRO_STRONG
      do j=0,Mm+1
        do i=0,Lm+1
          l0=l0+1
          zeta_das(i,j)=xmin(l0)
        enddo
      enddo
#endif
#if defined DAS_COUNT_PSICHI
      if ( count_non_null_psi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null psi = ',count_non_null_psi
      endif 
      if ( count_non_null_chi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null chi = ',count_non_null_chi
      endif 
#endif
      if ( l0 .ne . NDIM ) then
        write(*, '(6x,A)') 
     &     '!!! das_vec2a: ERROR in dimenision'
        stop
      endif 
!      
      return
      end

