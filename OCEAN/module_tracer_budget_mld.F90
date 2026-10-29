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

      module module_tracer_budget_mld
!============================================================================
!!                   ***  MODULE  module_tracer_budget_mld  ***
!!
!! This module manages tracer budget averaged over MLD
!============================================================================

      implicit none
 
      private 
 
      public :: budget_entrainment_term, budget_average_over_mld

      contains
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
#  ifdef TS_MIX_IMP
              if (mode == mode_t3dmix_ISO) then
                TVmix_mld_local(i,j,itrc)=0.
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
#  ifdef TS_MIX_IMP
                if (mode == mode_t3dmix_ISO) then
                  TVmix_mld_local(i,j,itrc)=TVmix_mld_local(i,j,itrc)+&
                                   TVmix_local(i,j,k,itrc)*cff
                endif             
#  endif
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
#   ifdef TS_MIX_IMP
                if (mode == mode_t3dmix_ISO) then
                TVmix_mld_local(i,j,itrc)= TVmix_mld_local(i,j,itrc)*rmask(i,j)
                endif             
#   endif
              endif             
#  endif
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


      end module module_tracer_budget_mld
