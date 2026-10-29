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

      module module_mld
!============================================================================
!!                   ***  MODULE  module_mld  ***
!!
!! This module computes MLD and index depending on the wanted criterion
!! (turbulence, density, temperature)
!============================================================================

      implicit none
 
      private 
 
      public :: idx_mld_density,idx_mld_temperature,idx_mld_turbulence

      contains
!==========================================================================
       subroutine idx_reference_depth(i,j,mld_depth_ref,indx_depth)
!----------------------------------------------------------------------
! Purposes : - find vertical index of reference depth (mld_depth_ref)
!
!----------------------------------------------------------------------
        use param
        use scalars
!#include "grid.h"
#include "ocean3d.h"

        ! Arguments 
        integer, intent(in) :: i,j
        real, intent(in) :: mld_depth_ref
        ! Outputs
        integer, intent(out) :: indx_depth
        ! Local declarations
        integer :: k
        real(kind=8) :: nearest_diff, diff_k
       
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

       end subroutine idx_reference_depth

!==========================================================================
       subroutine mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kbl_local,mld_local)
!----------------------------------------------------------------------
! Purposes : compute the depth of MLD from vertical index
!
!----------------------------------------------------------------------
        use param
#include "grid.h"
#include "ocean3d.h"

        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        integer, intent(inout) :: kbl_local(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld_local(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: i,j,k

        do j=Jstr,Jend
         do i=Istr,Iend
#  ifdef MASKING
          mld_local(i,j)=z_w(i,j,N)-z_w(i,j,kbl_local(i,j)-1) &
     &                      *rmask(i,j)
#  else
          mld_local(i,j)=z_w(i,j,N)-z_w(i,j,kbl_local(i,j)-1)
#  endif
         enddo
        enddo
       end subroutine mld_from_vertical_idx


!==========================================================================
       subroutine idx_mld_density(Istr,Iend,Jstr,Jend, &
                       mld_crit_D,mld_depth_ref, &
                       kbl_dens,mld)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixed layer depth
!              for density (kbl_dens) criteria, using the reference
!              depth (mld_depth_ref)
!            - compute the depth of MLD for density criteria 
!              (stored in mld)
!
! Called by : step3d_t.F, t3dmix_ISO.F, t3dmix_S.F, t3dmix_spg.F
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"
     
        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in) :: mld_crit_D,mld_depth_ref
        integer, intent(inout) :: kbl_dens(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k 
     
        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth mld_depth_ref
        ! -----------------------------------------------------------------
          call idx_reference_depth(i,j,mld_depth_ref,indx_depth)
        ! ---------------------------------------------------
        ! Compute kbl_dens (density)
        ! starting from indx_depth downward
        ! ---------------------------------------------------
          kbl_dens(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff.eq.0.and.rho(i,j,k).le. &
     &        rho(i,j,indx_depth)+mld_crit_D) then
            kbl_dens(i,j)=k
           endif
           if (rho(i,j,k).gt. &
     &        rho(i,j,indx_depth)+mld_crit_D)  cff=1
          enddo
         enddo
        enddo

        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
        call mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kbl_dens,mld)
       end subroutine idx_mld_density


!==========================================================================
       subroutine idx_mld_temperature(Istr,Iend,Jstr,Jend, &
                       mld_crit_T,mld_depth_ref, &
                       kbl_temp,mld)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixed layer depth
!              for temperature (kbl_temp) criteria, using the reference
!              depth (mld_depth_ref)
!            - compute the depth of MLD for temperature criteria 
!              (stored in mld)
!
! Called by : step3d_t.F, t3dmix_ISO.F, t3dmix_S.F, t3dmix_spg.F
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"

        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in) :: mld_crit_T,mld_depth_ref
        integer, intent(inout) :: kbl_temp(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k
        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth mld_depth_ref
        ! -----------------------------------------------------------------
          call idx_reference_depth(i,j,mld_depth_ref,indx_depth)
        ! ---------------------------------------------------
        ! Compute kbl_temp (temperature)
        ! starting from indx_depth downward
        ! ---------------------------------------------------
          kbl_temp(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff1.eq.0.and.t(i,j,k,nnew,itemp).ge. &
     &        t(i,j,indx_depth,nnew,itemp)-mld_crit_T) then
             kbl_temp(i,j)=k
           endif
           if (t(i,j,k,nnew,itemp).lt. &
     &        t(i,j,indx_depth,nnew,itemp)-mld_crit_T) cff1=1
          enddo
         enddo
        enddo

        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
        call mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kbl_temp,mld)
      
       end subroutine idx_mld_temperature


!==========================================================================
       subroutine idx_mld_turbulence(Istr,Iend,Jstr,Jend, &
                       mld_crit_AKt,mld_depth_ref, &
                       kbl_turb, mld)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixing layer depth
!              for both turbulence (kbl_turb) criteria, using the reference
!              depth (mld_depth_ref)
!            - compute the depth of MLD for turbulence criterion
!              (stored in hbl)
!
! Called by : gls_mixing.F
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"
#include "mixing.h"

        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in) :: mld_crit_AKt,mld_depth_ref
        integer,intent(inout) :: kbl_turb(GLOBAL_2D_ARRAY)
        real,intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k


        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth mld_depth_ref
        ! -----------------------------------------------------------------
          call idx_reference_depth(i,j,mld_depth_ref,indx_depth)
        ! ---------------------------------------------------
        ! Compute kbl_turb
        ! both starting from indx_depth downward
        ! ---------------------------------------------------
          kbl_turb(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff.eq.0.and.Akt(i,j,k,itemp).gt. &
     &        mld_crit_AKt) then
            kbl_turb(i,j)=k
           endif
           if (Akt(i,j,k,itemp).le. &
     &        mld_crit_AKt)  cff=1

          enddo
         enddo
        enddo

        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
        call mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kbl_turb,mld)

       end subroutine idx_mld_turbulence

      end module module_mld
