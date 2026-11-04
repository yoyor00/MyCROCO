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
      module module_tracer_budget_mld
!============================================================================
!!                   ***  MODULE  module_tracer_budget_mld  ***
!!
!! This module manages tracer budget averaged over MLD
!============================================================================
#include "cppdefs.h" 

# if defined DIAGNOSTICS_TS_MLD

#  ifdef MASKING
#   define SWITCH *
#  else
#   define SWITCH !
#  endif

      implicit none
 
      private 
 
      public :: get_diagTS_over_mld

      contains
!==========================================================================
!==========================================================================
       subroutine get_diagTS_over_mld(Istr,Iend,Jstr,Jend, &
                  itrc,T_mld_nnew,T_mld_nstp,mode, &
                  kbl_turb_local,kbl_dens_local,kbl_temp_local, &
                  mld_turb_local,mld_dens_local,mld_temp_local)
!--------------------------------------------------------------------------
! Purposes : main routine to get TS diagnostics integrated over MLD
!            3 steps : (1) get mld and vertical index of mld
!                      (2) compute the integral of TS diagnostics over MLD
!                      (3) compute rate and entrainement/detrainement terms
!
! Called by : step3d_t.F, t3dmix_ISO.F, t3dmix_S.F, t3dmix_spg.F
!--------------------------------------------------------------------------
        use param
        use croco_namelist, only: mld_crit_type

        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend,itrc
        integer, intent(in) :: mode
        integer, intent(inout) :: kbl_turb_local(GLOBAL_2D_ARRAY), &
                                  kbl_dens_local(GLOBAL_2D_ARRAY), &
                                  kbl_temp_local(GLOBAL_2D_ARRAY)
        real, intent(inout) :: mld_turb_local(GLOBAL_2D_ARRAY), &
                               mld_dens_local(GLOBAL_2D_ARRAY), &
                               mld_temp_local(GLOBAL_2D_ARRAY)
        real :: T_mld_nnew(PRIVATE_2D_SCRATCH_ARRAY), &
                T_mld_nstp(PRIVATE_2D_SCRATCH_ARRAY)

        ! Local declarations
        integer :: mode_step3d=1
        integer :: mode_t3dmix_S=2
        integer :: mode_t3dmix_ISO=3
        integer :: mode_t3dmix_spg=4
        integer :: kbl_local(GLOBAL_2D_ARRAY)


        call get_mld_idx(Istr,Iend,Jstr,Jend, &
                kbl_turb_local,kbl_dens_local,kbl_temp_local, &
                mld_turb_local,mld_dens_local,mld_temp_local)

        if (mld_crit_type == 'turbulence') then
                kbl_local = kbl_turb_local
        else if (mld_crit_type == 'density') then
                kbl_local = kbl_dens_local
        else if (mld_crit_type == 'temperature') then
                kbl_local = kbl_temp_local
        endif

        call integrate_over_mld(Istr,Iend,Jstr,Jend, &
                       itrc, T_mld_nnew, T_mld_nstp,mode,kbl_local)

        if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_ISO &
                   .or. mode == mode_t3dmix_spg) then 
          call get_rate_and_entrainement_terms(Istr,Iend,Jstr,Jend,itrc, &
                                  T_mld_nnew,T_mld_nstp,kbl_local)
#  if (!defined TS_DIF2 && !defined TS_DIF4 && !defined SPONGE_DIF2)
        else if (mode == mode_step3d) then
          call get_rate_and_entrainement_terms(Istr,Iend,Jstr,Jend,itrc, &
                                  T_mld_nnew,T_mld_nstp,kbl_local)
#  endif
        endif
       end subroutine get_diagTS_over_mld

!==========================================================================
!==========================================================================
       subroutine get_mld_idx(Istr,Iend,Jstr,Jend, &
                       kbl_turb_local,kbl_dens_local,kbl_temp_local, &
                       mld_turb_local,mld_dens_local,mld_temp_local)
!--------------------------------------------------------------------------
! Purposes : - get the mld and mld vertical index using the three 
!              criteria: turbulence, density and temperature
!--------------------------------------------------------------------------
        use module_mld, only: idx_mld_density, &
                           idx_mld_temperature, &
                           idx_mld_turbulence
        use croco_namelist, only: mld_crit_D, mld_crit_T, &
                                  mld_crit_AKt, mld_depth_ref
        use param
        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        integer, intent(inout) :: kbl_turb_local(GLOBAL_2D_ARRAY), &
                                  kbl_dens_local(GLOBAL_2D_ARRAY), &
                                  kbl_temp_local(GLOBAL_2D_ARRAY) 
        real, intent(inout) :: mld_turb_local(GLOBAL_2D_ARRAY), &
                               mld_dens_local(GLOBAL_2D_ARRAY), &
                               mld_temp_local(GLOBAL_2D_ARRAY)

        call idx_mld_density(Istr,Iend,Jstr,Jend, &
            mld_crit_D,mld_depth_ref, &
            kbl_dens_local,mld_dens_local)
        call idx_mld_temperature(Istr,Iend,Jstr,Jend, &
            mld_crit_T,mld_depth_ref, &
            kbl_temp_local,mld_temp_local)
#    if defined LMD_SKPP || defined GLS_MIXING || defined ANA_VMIX
        call idx_mld_turbulence(Istr,Iend,Jstr,Jend, & 
            mld_crit_AKt,mld_depth_ref, &
            kbl_turb_local,mld_turb_local) 
#    endif
       end subroutine get_mld_idx

!==========================================================================
!==========================================================================
       subroutine integrate_over_mld(Istr,Iend,Jstr,Jend, &
                       itrc, T_mld_nnew, T_mld_nstp,mode, &
                       kbl_local)
!--------------------------------------------------------------------------
! Purposes : - integrate budget terms of TS diagnostics over MLD
!--------------------------------------------------------------------------
        use param
        use scalars
#include "diagnostics.h"
#include "grid.h"
#include "ocean3d.h"

        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend,itrc
        integer, intent(in) :: mode
        integer, intent(in) :: kbl_local(GLOBAL_2D_ARRAY)
        real :: T_mld_nnew(PRIVATE_2D_SCRATCH_ARRAY), &
                T_mld_nstp(PRIVATE_2D_SCRATCH_ARRAY)

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
              TXadv_mld(i,j,itrc)=0.
              TYadv_mld(i,j,itrc)=0.
              TVadv_mld(i,j,itrc)=0.
              TVmix_mld(i,j,itrc)=0.
              THmix_mld(i,j,itrc)=0.
              TForc_mld(i,j,itrc)=0.
              Trate_mld(i,j,itrc)=0.
              Tentr_mld(i,j,itrc)=0.
              Taver_mld(i,j,itrc)=0.
            else if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_ISO &
                   .or. mode == mode_t3dmix_spg) then 
              THmix_mld(i,j,itrc)=0.
              Taver_mld(i,j,itrc)=0.
#  ifdef TS_MIX_IMP
              if (mode == mode_t3dmix_ISO) then 
                TVmix_mld(i,j,itrc)=0.
              endif
#  endif
#  if defined TS_SPONGE_IMP || defined TS_MIX_IMP
              if (mode == mode_t3dmix_spg) then 
                TVmix_mld(i,j,itrc)=0.
              endif
#  endif
            endif
            T_mld_nnew(i,j)=0.
            T_mld_nstp(i,j)=0.
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            kmld=kbl_local(i,j)
            do k=N,kmld,-1
              if (mode == mode_step3d) then
                TXadv_mld(i,j,itrc)=(TXadv_mld(i,j,itrc)+TXadv(i,j,k,itrc)) SWITCH rmask(i,j)
                TYadv_mld(i,j,itrc)=(TYadv_mld(i,j,itrc)+TYadv(i,j,k,itrc)) SWITCH rmask(i,j)
                TVadv_mld(i,j,itrc)=(TVadv_mld(i,j,itrc)+TVadv(i,j,k,itrc)) SWITCH rmask(i,j)
                TVmix_mld(i,j,itrc)=(TVmix_mld(i,j,itrc)+TVmix(i,j,k,itrc)) SWITCH rmask(i,j)
                THmix_mld(i,j,itrc)=(THmix_mld(i,j,itrc)+THmix(i,j,k,itrc)) SWITCH rmask(i,j)
                TForc_mld(i,j,itrc)=(TForc_mld(i,j,itrc)+TForc(i,j,k,itrc)) SWITCH rmask(i,j)
              else if (mode == mode_t3dmix_S .or. mode == mode_t3dmix_spg .or. &
                      mode == mode_t3dmix_ISO) then
                THmix_mld(i,j,itrc)=(THmix_mld(i,j,itrc)+THmix(i,j,k,itrc)) SWITCH rmask(i,j)
#  ifdef TS_MIX_IMP
                if (mode == mode_t3dmix_ISO) then 
                  TVmix_mld(i,j,itrc)=(TVmix_mld(i,j,itrc)+TVmix(i,j,k,itrc)) SWITCH rmask(i,j)
                endif
#  endif
#  if defined TS_SPONGE_IMP || defined TS_MIX_IMP
!
! TVmix was updated by the implicit sponge (or isopycnal) mixing
!
                if (mode == mode_t3dmix_spg) then 
                  TVmix_mld(i,j,itrc)=(TVmix_mld(i,j,itrc)+TVmix(i,j,k,itrc)) SWITCH rmask(i,j)
                endif
#  endif
              endif
              cff=Hz(i,j,k)/(pm(i,j)*pn(i,j))
              T_mld_nnew(i,j)=T_mld_nnew(i,j)+t(i,j,k,nnew,itrc)*cff
            enddo
!
! Mean tracer over the MLD 
! (updated in t3dmix routines if called)
!
            Taver_mld(i,j,itrc)=T_mld_nnew(i,j) &
                              /(z_w(i,j,N)-z_w(i,j,kmld-1)) &
                              *(pm(i,j)*pn(i,j)) SWITCH rmask(i,j)
          enddo
        enddo

       end subroutine integrate_over_mld

!==========================================================================
!==========================================================================
       subroutine get_rate_and_entrainement_terms(Istr,Iend,Jstr,Jend,itrc, &
                                      T_mld_nnew,T_mld_nstp,kbl_local)
!--------------------------------------------------------------------------
! Purposes : - compute the rate and entrainement/detrainement terms for 
!              TS diagnostics integrated over MLD 
!--------------------------------------------------------------------------
        use croco_namelist, only: dt
        use scalars
        use param
#include "diagnostics.h"
#include "grid.h"
#include "ocean3d.h"
        ! Arguments
        integer, intent(in) :: Istr,Iend,Jstr,Jend,itrc
        real :: T_mld_nnew(PRIVATE_2D_SCRATCH_ARRAY), &
                T_mld_nstp(PRIVATE_2D_SCRATCH_ARRAY)
        integer,intent(in) :: kbl_local(GLOBAL_2D_ARRAY)
        ! Local declarations
        integer :: i,j,k,kmld
        real :: cff

        do j=Jstr,Jend
          do i=Istr,Iend
            if (kbl_nstp(i,j).eq.0) kbl_nstp(i,j)=kbl_local(i,j)
            kmld=kbl_nstp(i,j)
            do k=N,kmld,-1
              cff=Hz_bak(i,j,k)/(pm(i,j)*pn(i,j))
              T_mld_nstp(i,j)=T_mld_nstp(i,j)+ &
                             t(i,j,k,nstp,itrc)*cff
            enddo
            if (itrc .eq. NT) kbl_nstp(i,j)=kbl_local(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            Trate_mld(i,j,itrc)=(T_mld_nnew(i,j)-T_mld_nstp(i,j))/dt SWITCH rmask(i,j)
            Tentr_mld(i,j,itrc)=(Trate_mld(i,j,itrc)- &
                                 TXadv_mld(i,j,itrc)- &
                                 TYadv_mld(i,j,itrc)- &
                                 TVadv_mld(i,j,itrc)- &
                                 TVmix_mld(i,j,itrc)- &
                                 THmix_mld(i,j,itrc)- &
                                 TForc_mld(i,j,itrc)) SWITCH rmask(i,j)
          enddo
        enddo

       end subroutine get_rate_and_entrainement_terms
# endif /* DIAGNOSTICS_TS_MLD */
      end module module_tracer_budget_mld

