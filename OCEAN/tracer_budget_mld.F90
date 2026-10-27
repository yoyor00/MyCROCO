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
 
      public :: idx_mld_density,idx_mld_temperature, idx_mld_turbulence, &
              budget_entrainment_term, budget_average_over_mld

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
#include "grid.h"  
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
        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
#  ifdef MASKING
          mld(i,j)=z_w(i,j,N)-z_w(i,j,kbl_dens(i,j)-1) &
     &                     *rmask(i,j)
#  else
          mld(i,j)=z_w(i,j,N)-z_w(i,j,kbl_dens(i,j)-1) 
     
#  endif
         enddo
        enddo
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
#include "grid.h"
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
        ! ----------------------------------------------------
        ! Store MLD
        ! ----------------------------------------------------
#  ifdef MASKING
          mld(i,j)=z_w(i,j,N)-z_w(i,j,kbl_temp(i,j)-1) &
     &                      *rmask(i,j)
#  else
          mld(i,j)=z_w(i,j,N)-z_w(i,j,kbl_temp(i,j)-1)
#  endif
         enddo
        enddo

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
#include "grid.h"
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
        ! ----------------------------------------------------
        ! Store MLD for turbulence
        ! ----------------------------------------------------
#  ifdef MASKING
          mld(i,j)=z_w(i,j,N)-z_w(i,j,kbl_turb(i,j)-1) &
     &                     *rmask(i,j)
#  else
          mld(i,j)=z_w(i,j,N)-z_w(i,j,kbl_turb(i,j)-1)

#  endif
         enddo
        enddo
       end subroutine idx_mld_turbulence


!==========================================================================
       subroutine budget_average_over_mld(Istr,Iend,Jstr,Jend, &
                       itrc, T_mld_nnew, T_mld_nstp,mode, & 
                       TXadv_mld_local,TYadv_mld_local,TVadv_mld_local, &
                       TVmix_mld_local,THmix_mld_local,TForc_mld_local, &
                       Tentr_mld_local,Trate_mld_local,Taver_mld_local, &
                       TXadv_local,TYadv_local,TVadv_local, & 
                       TVmix_local,THmix_local,TForc_local, &
                       kbl_local)
!----------------------------------------------------------------------
! Purposes : - average budget terms of TS diagnostic over MLD
!
! Called by : step3d_t.F, t3dmix_ISO.F, t3dmix_S.F, t3dmix_spg.F
!----------------------------------------------------------------------
        use param
        use scalars
        use croco_namelist, only: dt
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "diagnostics.h"

        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend,itrc
        integer, intent(in) :: mode 
        integer, intent(in) :: kbl_local(GLOBAL_2D_ARRAY)
        real :: T_mld_nnew(PRIVATE_2D_SCRATCH_ARRAY), &
                T_mld_nstp(PRIVATE_2D_SCRATCH_ARRAY)
        real,intent(inout) :: TXadv_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TYadv_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TVadv_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TVmix_mld_local(GLOBAL_2D_ARRAY,NT), &
                              THmix_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TForc_mld_local(GLOBAL_2D_ARRAY,NT), &
                              Tentr_mld_local(GLOBAL_2D_ARRAY,NT), &
                              Trate_mld_local(GLOBAL_2D_ARRAY,NT), &
                              Taver_mld_local(GLOBAL_2D_ARRAY,NT)
        real,intent(in)    :: TXadv_local(GLOBAL_2D_ARRAY,N,NT), &
                              TYadv_local(GLOBAL_2D_ARRAY,N,NT), &
                              TVadv_local(GLOBAL_2D_ARRAY,N,NT), &
                              TVmix_local(GLOBAL_2D_ARRAY,N,NT), &
                              THmix_local(GLOBAL_2D_ARRAY,N,NT), &
                              TForc_local(GLOBAL_2D_ARRAY,N,NT)

        ! Local declarations
        integer :: i,j,k,kmld
        integer :: mode_step3d=1
        integer :: mode_t3dmix_S=2
        integer :: mode_t3dmix_ISO=3
        integer :: mode_t3dmix_spg=4
        real :: cff

        do j=Jstr,Jend
          do i=Istr,Iend
            if (mode == mode_step3d) then
              TXadv_mld_local(i,j,itrc)=0.
              TYadv_mld_local(i,j,itrc)=0.
              TVadv_mld_local(i,j,itrc)=0.
              TVmix_mld_local(i,j,itrc)=0.
              THmix_mld_local(i,j,itrc)=0.
              TForc_mld_local(i,j,itrc)=0.
              Trate_mld_local(i,j,itrc)=0.
              Tentr_mld_local(i,j,itrc)=0.
            else if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_spg .or. &
                      mode == mode_t3dmix_ISO) then
              THmix_mld_local(i,j,itrc)=0.
              Taver_mld_local(i,j,itrc)=0.
              if (mode == mode_t3dmix_ISO) then
#  ifdef TS_MIX_IMP
                TVmix_mld_local(i,j,itrc)=0.
#  endif
              endif
            endif
            T_mld_nnew(i,j)=0.
            T_mld_nstp(i,j)=0.
          enddo
        enddo

        do j=Jstr,Jend
          do i=Istr,Iend
            kmld=kbl_local(i,j)
            do k=N,kmld,-1
              cff=Hz(i,j,k)/(z_w(i,j,N)-z_w(i,j,kmld-1))
              if (mode == mode_step3d) then
                TXadv_mld_local(i,j,itrc)=TXadv_mld_local(i,j,itrc)+ &
                                   TXadv_local(i,j,k,itrc)*cff
                TYadv_mld_local(i,j,itrc)=TYadv_mld_local(i,j,itrc)+ &
                                   TYadv_local(i,j,k,itrc)*cff
                TVadv_mld_local(i,j,itrc)=TVadv_mld_local(i,j,itrc)+ &
                                   TVadv_local(i,j,k,itrc)*cff
                TVmix_mld_local(i,j,itrc)=TVmix_mld_local(i,j,itrc)+ &
                                   TVmix_local(i,j,k,itrc)*cff
                TForc_mld_local(i,j,itrc)=TForc_mld_local(i,j,itrc)+ &
                                   TForc_local(i,j,k,itrc)*cff
              else if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_spg .or. &
                      mode == mode_t3dmix_ISO) then
                THmix_mld_local(i,j,itrc)=THmix_mld_local(i,j,itrc)+ &
                                   THmix_local(i,j,k,itrc)*cff
                if (mode == mode_t3dmix_ISO) then
#  ifdef TS_MIX_IMP
                  TVmix_mld_local(i,j,itrc)=TVmix_mld_local(i,j,itrc)+&
                                   TVmix_local(i,j,k,itrc)*cff
#  endif
                endif             
              endif             
              T_mld_nnew(i,j)=T_mld_nnew(i,j)+ &
                                   t(i,j,k,nnew,itrc)*cff
#  ifdef MASKING
              if (mode == mode_step3d) then
                TXadv_mld_local(i,j,itrc)= TXadv_mld_local(i,j,itrc)*rmask(i,j)
                TYadv_mld_local(i,j,itrc)= TYadv_mld_local(i,j,itrc)*rmask(i,j)
                TVadv_mld_local(i,j,itrc)= TVadv_mld_local(i,j,itrc)*rmask(i,j)
                TVmix_mld_local(i,j,itrc)= TVmix_mld_local(i,j,itrc)*rmask(i,j)
                TForc_mld_local(i,j,itrc)= TForc_mld_local(i,j,itrc)*rmask(i,j)
              else if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_spg .or. &
                      mode == mode_t3dmix_ISO) then
                THmix_mld_local(i,j,itrc)= THmix_mld_local(i,j,itrc)*rmask(i,j)
                if (mode == mode_t3dmix_ISO) then
#   ifdef TS_MIX_IMP
                TVmix_mld_local(i,j,itrc)= TVmix_mld_local(i,j,itrc)*rmask(i,j)
#   endif
#  endif
                endif             
              endif             
            enddo 
            if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_spg .or. &
                    mode == mode_t3dmix_ISO) then
              Taver_mld_local(i,j,itrc)=T_mld_nnew(i,j) ! mean over mld
#  ifdef MASKING
              Taver_mld_local(i,j,itrc) = Taver_mld_local(i,j,itrc)*rmask(i,j)
#  endif
            endif
          enddo
        enddo


       end subroutine budget_average_over_mld

!==========================================================================
       subroutine budget_entrainment_term(Istr,Iend,Jstr,Jend,itrc, &
                                  T_mld_nnew, T_mld_nstp,Trate_mld_local, &
                                     TXadv_mld_local,TYadv_mld_local,TVadv_mld_local, &
                                     TVmix_mld_local,THmix_mld_local,TForc_mld_local, &
                                     Tentr_mld_local,kbl_nstp_local,kbl_local)
!----------------------------------------------------------------------
! Purposes : - compute the entrainment term for TS diagnostique averaged on MLD 
!
! Called by : step3d_t.F, t3dmix_ISO.F, t3dmix_S.F, t3dmix_spg.F
!----------------------------------------------------------------------
        use param
        use scalars
        use croco_namelist, only: dt
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "diagnostics.h"

        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend,itrc
        real :: T_mld_nnew(PRIVATE_2D_SCRATCH_ARRAY), &
                T_mld_nstp(PRIVATE_2D_SCRATCH_ARRAY)
        integer,intent(in) :: kbl_local(GLOBAL_2D_ARRAY)
        real,intent(inout) :: TXadv_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TYadv_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TVadv_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TVmix_mld_local(GLOBAL_2D_ARRAY,NT), &
                              THmix_mld_local(GLOBAL_2D_ARRAY,NT), &
                              TForc_mld_local(GLOBAL_2D_ARRAY,NT), &
                              Tentr_mld_local(GLOBAL_2D_ARRAY,NT), &
                              Trate_mld_local(GLOBAL_2D_ARRAY,NT)
        integer,intent(inout) :: kbl_nstp_local(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: i,j,k,kmld
        real :: cff

        do j=Jstr,Jend
          do i=Istr,Iend
            if (kbl_nstp_local(i,j).eq.0) kbl_nstp_local(i,j)=kbl_local(i,j)
            kmld=kbl_nstp_local(i,j)
            do k=N,kmld,-1
              cff=Hz_bak(i,j,k)/(z_w(i,j,N)-z_w(i,j,kmld-1))
              T_mld_nstp(i,j)=T_mld_nstp(i,j)+ &
                             t(i,j,k,nstp,itrc)*cff
            enddo
            if (itrc .eq. NT) kbl_nstp_local(i,j)=kbl_local(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            Trate_mld_local(i,j,itrc)=(T_mld_nnew(i,j)-T_mld_nstp(i,j))/dt 
#     if defined MASKING
            Trate_mld_local(i,j,itrc)=Trate_mld_local(i,j,itrc)*rmask(i,j)
#     endif
            Tentr_mld_local(i,j,itrc)=Trate_mld_local(i,j,itrc)- &
                                      TXadv_mld_local(i,j,itrc)-  &
                                      TYadv_mld_local(i,j,itrc)-  &
                                      TVadv_mld_local(i,j,itrc)-  &
                                      TVmix_mld_local(i,j,itrc)-  &
#     if (defined TS_DIF2 || defined TS_DIF4 || defined SPONGE_DIF2)
                                THmix_mld_local(i,j,itrc)-  &
#     endif
                                TForc_mld_local(i,j,itrc) 
#     if defined MASKING
             Tentr_mld_local(i,j,itrc)=Tentr_mld_local(i,j,itrc)*rmask(i,j)
#     endif
          enddo
        enddo


       end subroutine budget_entrainment_term


      end module tracer_budget_mld
