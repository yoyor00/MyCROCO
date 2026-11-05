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

#  ifdef MASKING
#   define SWITCH *
#  else
#   define SWITCH !
#  endif

      implicit none
 
      private 
 
      public :: idx_mld_density,idx_mld_temperature,idx_mld_turbulence

      contains
!==========================================================================
       subroutine idx_reference_depth(i,j,depth_ref,indx_depth)
!----------------------------------------------------------------------
! Purposes : - find vertical index of reference depth (depth_ref)
!
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"

        ! Arguments 
        integer, intent(in) :: i,j
        real, intent(in) :: depth_ref
        ! Outputs
        integer, intent(out) :: indx_depth
        ! Local declarations
        integer :: k
        real(kind=8) :: nearest_diff, diff_k
       
          if (z_w(i,j,0) .gt. -depth_ref) then
        ! if bathymetry shallower than reference depth, clamp to bottom
                      indx_depth = 1
          else
                      indx_depth = N
                      nearest_diff = abs(z_r(i,j,N)+depth_ref)
           do k=N-1,1,-1
            diff_k = abs(z_r(i,j,k)+depth_ref)
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
                       kmld,mld)
!----------------------------------------------------------------------
! Purposes : compute the depth of MLD from vertical index
!
!----------------------------------------------------------------------
        use param
#include "grid.h"
#include "ocean3d.h"

        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        integer, intent(inout) :: Kmld(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: i,j,k

        do j=Jstr,Jend
         do i=Istr,Iend
          mld(i,j)=(z_w(i,j,N)-z_w(i,j,Kmld(i,j)-1)) SWITCH rmask(i,j)
         enddo
        enddo
       end subroutine mld_from_vertical_idx


!==========================================================================
       subroutine idx_mld_density(Istr,Iend,Jstr,Jend, &
                       threshold_dens,depth_ref, &
                       kmld_dens,mld)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixed layer depth
!              for density (kmld_dens) criteria, using the reference
!              depth (depth_ref)
!            - compute the depth of MLD for density criteria 
!
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"
     
        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in) :: threshold_dens,depth_ref
        integer, intent(inout) :: kmld_dens(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k 
     
        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth depth_ref
        ! -----------------------------------------------------------------
          call idx_reference_depth(i,j,depth_ref,indx_depth)
        ! ---------------------------------------------------
        ! Compute kmld_dens (density)
        ! starting from indx_depth downward
        ! ---------------------------------------------------
          kmld_dens(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff.eq.0.and.rho(i,j,k).le. &
     &        rho(i,j,indx_depth)+threshold_dens) then
            kmld_dens(i,j)=k
           endif
           if (rho(i,j,k).gt. &
     &        rho(i,j,indx_depth)+threshold_dens)  cff=1
          enddo
         enddo
        enddo

        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
        call mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kmld_dens,mld)
       end subroutine idx_mld_density


!==========================================================================
       subroutine idx_mld_temperature(Istr,Iend,Jstr,Jend, &
                       threshold_temp,depth_ref, &
                       kmld_temp,mld)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixed layer depth
!              for temperature (kmld_temp) criteria, using the reference
!              depth (depth_ref)
!            - compute the depth of MLD for temperature criteria 
!
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"

        ! Arguments 
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in) :: threshold_temp,depth_ref
        integer, intent(inout) :: kmld_temp(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k
        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth depth_ref
        ! -----------------------------------------------------------------
          call idx_reference_depth(i,j,depth_ref,indx_depth)
        ! ---------------------------------------------------
        ! Compute kmld_temp (temperature)
        ! starting from indx_depth downward
        ! ---------------------------------------------------
          kmld_temp(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff1.eq.0.and.t(i,j,k,nnew,itemp).ge. &
     &        t(i,j,indx_depth,nnew,itemp)-threshold_temp) then
             kmld_temp(i,j)=k
           endif
           if (t(i,j,k,nnew,itemp).lt. &
     &        t(i,j,indx_depth,nnew,itemp)-threshold_temp) cff1=1
          enddo
         enddo
        enddo

        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
        call mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kmld_temp,mld)
      
       end subroutine idx_mld_temperature


!==========================================================================
       subroutine idx_mld_turbulence(Istr,Iend,Jstr,Jend, &
                       threshold_turb,depth_ref, &
                       kmld_turb, mld)
!----------------------------------------------------------------------
! Purposes : - compute the vertical indices corresponding to the mixing layer depth
!              for both turbulence (kmld_turb) criteria, using the reference
!              depth (depth_ref)
!            - compute the depth of MLD for turbulence criterion
!
!----------------------------------------------------------------------
        use param
        use scalars
#include "ocean3d.h"
#include "mixing.h"

        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in) :: threshold_turb,depth_ref
        integer,intent(inout) :: kmld_turb(GLOBAL_2D_ARRAY)
        real,intent(inout) :: mld(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: indx_depth, cff1,cff
        integer :: i,j,k
        real(kind=8) :: nearest_diff, diff_k


        do j=Jstr,Jend
         do i=Istr,Iend
        ! -----------------------------------------------------------------
        ! Find the index corresponding to the reference depth depth_ref
        ! -----------------------------------------------------------------
          call idx_reference_depth(i,j,depth_ref,indx_depth)
        ! ---------------------------------------------------
        ! Compute kmld_turb
        ! both starting from indx_depth downward
        ! ---------------------------------------------------
          kmld_turb(i,j)=indx_depth
          cff=0
          cff1=0
          do k=indx_depth,1,-1
           if (cff.eq.0.and.Akt(i,j,k,itemp).gt. &
     &        threshold_turb) then
            kmld_turb(i,j)=k
           endif
           if (Akt(i,j,k,itemp).le. &
     &        threshold_turb)  cff=1

          enddo
         enddo
        enddo

        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
        call mld_from_vertical_idx(Istr,Iend,Jstr,Jend, &
                       kmld_turb,mld)

       end subroutine idx_mld_turbulence

      end module module_mld

