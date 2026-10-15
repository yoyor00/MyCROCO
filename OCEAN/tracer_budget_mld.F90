!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================

#include "cppdefs.h" 

      module tracer_budget_mld
!============================================================================
!!                   ***  MODULE  tracer_budget_mld  ***
!!
!! This module manages tracer budget averaged over MLD
!============================================================================

      implicit none
 
      private 
 
      public :: idx_mld_density_and_temp

      contains

!==========================================================================
       subroutine idx_mld_density_and_temp(Istr,Iend,Jstr,Jend,itrc)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixed layer depth
!              for both density (kbld) and temperature (kblt) criteria, using the reference
!              depth (mld_depth_ref)
!            - compute the depth of MLD for both density and temperature criteria 
!              (both stored in Tcrit_mld)
!           
! Called by : step3d_t.F, t3dmix_ISO.F, t3dmix_S.F, t3dmix_spg.F
!----------------------------------------------------------------------
        use croco_namelist, only: mld_crit_D, mld_crit_T, mld_depth_ref
        use param
        use scalars
#include "diagnostics.h" 
#include "mixing.h"  
#include "grid.h"  
!#include "param.h" 
#include "ocean3d.h"
     
        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend,itrc
     
        ! Local declarations
        integer :: kblt(PRIVATE_2D_SCRATCH_ARRAY),indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k 
     
        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth mld_depth_ref
        ! -----------------------------------------------------------------
          if (z_w(i,j,0) .gt. -mld_depth_ref) then
        ! if bathymetry shallower than reference depth, clamp to bottom
                      indx_depth = 1
          else
                      indx_depth = N
                      nearest_diff = abs(z_r(i,j,N)+mld_depth_ref)
           do k=N-1,1,-1
            diff_k = abs(z_r(i,j,k)+mld_depth_ref)
            if (diff_k .le. nearest_diff) then
               nearest_diff = diff_k
               indx_depth = k
            else
               exit
            endif
           enddo
          endif
        ! ---------------------------------------------------
        ! Compute kbld (density) and kblt (temperature)
        ! both starting from indx_depth downward
        ! ---------------------------------------------------
          kbld(i,j)=indx_depth
          kblt(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff.eq.0.and.rho(i,j,k).le. &
     &        rho(i,j,indx_depth)+mld_crit_D) then
            kbld(i,j)=k
           endif
           if (rho(i,j,k).gt. &
     &        rho(i,j,indx_depth)+mld_crit_D)  cff=1
        
           if (cff1.eq.0.and.t(i,j,k,nnew,itemp).ge. &
     &        t(i,j,indx_depth,nnew,itemp)-mld_crit_T) then
             kblt(i,j)=k
           endif
           if (t(i,j,k,nnew,itemp).lt. &
     &        t(i,j,indx_depth,nnew,itemp)-mld_crit_T) cff1=1
          enddo
        ! ----------------------------------------------------
        ! Store MLD: itrc=1 -> density, itrc=2 -> temperature
        ! ----------------------------------------------------
          if (itrc.eq.1) then
#  ifdef MASKING
           Tcrit_mld(i,j,itrc)=z_w(i,j,N)-z_r(i,j,kbld(i,j)) &
     &                     *rmask(i,j)
#  else
           Tcrit_mld(i,j,itrc)=z_w(i,j,N)-z_r(i,j,kbld(i,j)) 
     
#  endif
        
          elseif (itrc.eq.2) then
#  ifdef MASKING
           Tcrit_mld(i,j,itrc)=z_w(i,j,N)-z_r(i,j,kblt(i,j)) &
     &                      *rmask(i,j)
#  else
           Tcrit_mld(i,j,itrc)=z_w(i,j,N)-z_r(i,j,kblt(i,j)) 
#  endif
          endif
         enddo
        enddo
       end subroutine idx_mld_density_and_temp
     
      end module tracer_budget_mld
