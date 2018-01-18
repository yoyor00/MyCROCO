!***********************************************************************
! s-grid evolves with external mode (NH)
!    step 2: grid at n+1
!***********************************************************************

#  define zab3 UFx
!
!-----------------------------------------------------------------------
! Computes zeta(n+1)
!-----------------------------------------------------------------------
!
      if (FIRST_TIME_STEP.and.iif.eq.1) then
         cff8=1.D0
         cff9=0.D0
         cff10=0.D0
       elseif (FIRST_TIME_STEP.and.iif.eq.2) then
         cff8=1.5D0
         cff9=-0.5D0
        cff10=0.D0
      else
        cff8= 1.5D0+mybeta
        cff9=-2.0D0*mybeta-0.5D0
        cff10= mybeta
      endif

      do j=JstrV-2,Jend+1
        do i=IstrU-2,Iend+1
           zab3(i,j) =   cff8  * zeta(i,j,kstp2)
     &                 + cff9  * zeta(i,j,kbak2)
     &                 + cff10 * zeta(i,j,kold2)
        enddo
      enddo

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
         zeta (i,j,knew2)=
     &    ( zeta(i,j,kstp2)
     &        +dtfast*(
#if defined NBQ_ZETAEXP
     &    +wmean_nbq(i,j,kstp2)
#else
     &    +wmean_nbq(i,j,knew2)
#endif
     &    -0.5*(umean_nbq(i  ,j)
     &          *(zab3(i  ,j)
     &           -zab3(i-1,j))*pm_u(i,j)
#ifdef MASKING
     &         *rmask(i,j)*rmask(i-1,j)
#endif
     &         +umean_nbq(i+1,j)
     &          *(zab3(i+1,j)
     &           -zab3(i  ,j))*pm_u(i+1,j) 
#ifdef MASKING
     &         *rmask(i,j)*rmask(i+1,j)
#endif
     &          )
#ifdef MASKING
     &         *umask(i,j)*umask(i+1,j)
#endif

     &    -0.5*(vmean_nbq(i  ,j)
     &          *(zab3(i,j  ) 
     &           -zab3(i,j-1))*pm_v(i,j)
#ifdef MASKING
     &         *rmask(i,j)*rmask(i,j-1)
#endif
     &         +vmean_nbq(i,j+1)
     &          *(zab3(i,j+1)
     &           -zab3(i,j  ))*pm_v(i,j+1) 
#ifdef MASKING
     &         *rmask(i,j)*rmask(i,j+1)
#endif
     &         )
#ifdef MASKING
     &         *vmask(i,j)*vmask(i,j+1)
#endif
     &         ) )
#ifdef MASKING
     &         *rmask(i,j)
#endif
        enddo
      enddo

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew2))
# endif

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
           zeta_new(i,j)=zeta(i,j,knew2)
        enddo
      enddo

!
!-----------------------------------------------------------------------
! Add nudging terms
!-----------------------------------------------------------------------
!
#  ifdef ZNUDGING
#   ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew2),zetazon)
        endif
      endif
      if (iic.eq.ntstart) then
        if (FIRST_2D_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ssh(START_2D_ARRAY),sshzon)
        endif
      endif
#   endif  /* ZONAL_NUDGING */
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j) + dtfast*Znudgcof(i,j)
#   ifdef ZONAL_NUDGING
     &                                 *(sshzon(j)-zetazon(j))
#   else
     &                                 *(ssh(i,j)-zeta_new(i,j))
#   endif /* ZONAL_NUDGING */
        enddo
      enddo
#  endif /* ZNUDGING */
!
!-----------------------------------------------------------------------
! Computes zetarhs to use in momentum equations
!-----------------------------------------------------------------------
!

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
 	   zeta_new(i,j)=zeta_new(i,j) SWITCH rmask(i,j)
#if !defined NBQ_ZETAEXP
c in case of ZETAEXP rhobar has not yet been computed
c Dnew will be computed later
          Dnew(i,j)=(zeta_new(i,j)+h(i,j))
#if defined NBQ_MASS
     &       *rhobar_nbq(i,j,knew2)
#endif
#endif
        enddo
      enddo

!
!-----------------------------------------------------------------------
! Load new free-surface values into shared array
! Modify new free-surface to ensure that depth is > Dcrit for masked
! cells.
!-----------------------------------------------------------------------
!

#  if defined WET_DRY && defined MASKING
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew2)=zeta(i,j,knew2)+ 
     &                   (Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
        enddo
      enddo 
#  endif
!
!-----------------------------------------------------------------------
! Load rhs values into additional AGRIF shared array for nesting
!-----------------------------------------------------------------------
! 
#  ifdef AGRIF
      if (FIRST_2D_STEP) then
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,0)=zeta(i,j,kstp)       
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
          du_avg3(i,j,0)  = DUon(i,j)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
          dv_avg3(i,j,0)  = DVom(i,j)
          enddo
        enddo 
      endif

#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   Zt_avg3(START_2D_ARRAY,0))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   du_avg3(START_2D_ARRAY,0))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   dv_avg3(START_2D_ARRAY,0))
#   endif

#   ifdef RVTK_DEBUG_ADVANCED
       if (.not.agrif_Root()) then
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(Zt_avg3(:,:,0),'Zt_avg3 (index 0) step2d','r')
       call check_tab2d(DU_avg3(:,:,0),'DU_avg3 (index 0) step2d','u')
       call check_tab2d(DV_avg3(:,:,0),'DV_avg3 (index 0) step2d','v')
C$OMP END MASTER  
       endif
#   endif  
#  endif /* AGRIF */   
!
!-----------------------------------------------------------------------
! Compute wet/dry masks
!-----------------------------------------------------------------------
!
#  ifdef WET_DRY
      call wetdry_tile (Istr,Iend,Jstr,Jend)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew2))
#  endif      
!
!-----------------------------------------------------------------------
! Debug zeta
!-----------------------------------------------------------------------
!
!#  ifdef RVTK_DEBUG_ADVANCED
!C$OMP BARRIER
!C$OMP MASTER
!       call check_tab2d(zeta(:,:,knew),'zeta step2d #2','r')
!C$OMP END MASTER       
!#  endif       
!
!-----------------------------------------------------------------------
! Set boundary conditions for the free-surface
!-----------------------------------------------------------------------
!
      call zetabc_tile (Istr,Iend,Jstr,Jend)

!
!----------------------------------------------------------------------
! Compute time averaged fields over all short timesteps.
!
! Reset/initialise arrays for averaged fields during the first
! barotropic time step; Accumulate averages after that. Include
! physical boundary points, but not periodic ghost points or
! computation  MPI computational margins.
!----------------------------------------------------------------------
!
#  ifdef SOLVE3D
c LAURENT: unnecessary lines with M2FILTER_NONE (only at last time step): REMOVED
!      cff1=weight(1,iif)
!      cff2=weight(2,iif)
!      if (FIRST_2D_STEP) then
!        do j=JstrR,JendR
!          do i=IstrR,IendR
!            Zt_avg1(i,j)=cff1*zeta(i,j,knew2)
!          enddo
!        enddo 
!      else
!        do j=JstrR,JendR
!          do i=IstrR,IendR
!            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew2)
!          enddo
!        enddo
!      endif

!
!-----------------------------------------------------------------------
! Update Grid:
!-----------------------------------------------------------------------
!
       if ((iic.eq.1.and.iif==1)
     &     NSTEP_GRID) then
c setting flag_grid to 1
c enforces the recomputation of vertical grid derived arrays (Hz(r,w)_xxx_inv) at next time step
        flag_grid=1
        call set_depth_tile(Istr,Iend,Jstr,Jend
     &   ,resetfromrhobar
     &   ) 

#include "step2d_grid_ext.h"

        endif

#  undef zab3

#  endif /* SOLVE3D */


