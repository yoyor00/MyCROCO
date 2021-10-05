! $Id: cpl_prism_get.F 1574 2014-07-01 15:13:38Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
      SUBROUTINE cpl_prism_get (rcvtime)
!     !---------------------------------------------------------------------
!     !              ***  ROUTINE cpl_prism_get  ***
!     !
!     ! ** Purpose : - 
!     !           
!     !
!     !----------------------------------------------------------------------
#include "cppdefs.h"
#if defined OA_COUPLING || defined OW_COUPLING
      USE mod_prism             ! OASIS3 prism module
# ifdef AGRIF
      use Agrif_Util
# endif
      IMPLICIT NONE

#include "cppdefs.h"

#include "param.h"
#include "forces.h"
#include "mpi_cpl.h"
#include "ocean3d.h"
#include "scalars.h"
#include "grid.h"

      INTEGER    :: rcvtime, i, j, io, jo
      INTEGER    :: iok_taux, iok_tauy, iok_taum
      INTEGER    :: iok_srflx, iok_itemp , iok_isalt
      INTEGER    :: iok_wdry                                   ! OW_COUPLING 
      INTEGER    :: iok_wt0m1, iok_whs, iok_wdir               ! OW_COUPLING 
      INTEGER    :: iok_wbhd, iok_wubr, iok_wfoc, iok_wlm      ! OW_COUPLING
      INTEGER    :: iok_wubrx, iok_wubry 
      INTEGER    :: iok_wtwox, iok_wtwoy, iok_wtaux, iok_wtauy ! OW_COUPLING 
      INTEGER    :: imin_halo, imax_halo, jmin_halo, jmax_halo
      INTEGER    :: level
      CHARACTER*1 clevel

      real, parameter :: eps=1.e-20
      real,dimension(GLOBAL_2D_ARRAY) :: tab_t_tmp
      real,dimension(GLOBAL_2D_ARRAY) :: wave_dir 
#if (!defined OA_GRID_UV)
      real,dimension(GLOBAL_2D_ARRAY) :: tab_u_tmp
      real,dimension(GLOBAL_2D_ARRAY) :: tab_v_tmp
#endif

#ifdef AGRIF
      level=Agrif_Fixed()
#else
      level=0
#endif
      write(clevel,'(i1.1)') level
!      
!!$      if (level==0) then
!!$         MPI_master_only 
!!$     & write(*,*) 'CPL-CROCO : Get for parent at time', rcvtime
!!$      else		
!!$         MPI_master_only 
!!$     & write(*,*) 'CPL-CROCO : Get for child:', level,
!!$     &  'at time', rcvtime
!!$      endif
!
      imin_halo = 1  
      imax_halo = Lmmpi  
      jmin_halo = 1  
      jmax_halo = Mmmpi  

      if (WEST_INTER) then
         imin_halo = imin_halo - 1
      endif 
      if (EAST_INTER) then
         imax_halo = imax_halo + 1
      endif 
      if (SOUTH_INTER) then     
         jmin_halo = jmin_halo - 1
      endif
      if (NORTH_INTER) then
         jmax_halo = jmax_halo + 1
      endif 
!
      iok_srflx = 0
      iok_isalt = 0
      iok_itemp = 0
      iok_taux = 0
      iok_tauy = 0
      iok_taum = 0
!
!     add for OW_COUPLING
      iok_wdry = 0
      iok_wt0m1 = 0
      iok_whs = 0
      iok_wdir = 0
      iok_wbhd = 0
      iok_wubr = 0
      iok_wubrx = 0
      iok_wubry = 0
      iok_wfoc = 0
      iok_wlm = 0
      iok_wtwox = 0
      iok_wtwoy = 0
      iok_wtaux = 0
      iok_wtauy = 0
!     add end for OW_COUPLING
!
!      print *,'CROCO, number of received fields = ',krcv
      DO i = 1, nmaxfld
#ifdef OA_COUPLING
!-------------------------------------------------
         if (srcv_clname(i)=='RRMSRFL'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_srflx) 
            if (iok_srflx /= 0) then
!              Transmit to CROCO field
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               srflx = tab_t_tmp / (rho0*Cp) 
            endif
         end if  
!--------------------------------------------------
#ifdef SALINITY
         if (srcv_clname(i)=='RRMEVPR'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_isalt) 
            if (iok_isalt /= 0) then
!              Transmit to CROCO field : stflx(:,:,isalt)
!              Changing Unit from [kg/m^2/s] (<--> [mm/s]) to [PSU m/s]  --> coeff= 1/1000.
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               stflx(:,:,isalt) = tab_t_tmp * t(:,:,N,nstp,isalt) / 1000.
            endif
         end if
#endif
!-------------------------------------------------
         if (srcv_clname(i)=='RRMSTFL'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_itemp) 
            if (iok_itemp /= 0)  then
!              Transmit to CROCO field : stflx(:,:,itemp)
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               stflx(:,:,itemp) = tab_t_tmp / (rho0*Cp)
            endif     
         end if
!-------------------------------------------------
         if (srcv_clname(i)=='RRMTAUX'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_taux) 
            if (iok_taux /= 0) then  
#if (!defined OA_GRID_UV)  
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TAUX on grid T --> U
               do jo=jmin_halo-1,jmax_halo+1
               do io=imin_halo,imax_halo+1
                  tab_u_tmp(io,jo)=0.5*(tab_t_tmp(io-1,jo)+tab_t_tmp(io,jo))
               enddo
               enddo
!              Transmit to CROCO field : sustr
               sustr  = tab_u_tmp / rho0
               if (iok_taum == 0) tab_u_tmp = tab_t_tmp ! temporary storage...
#else
!              Transmit to CROCO field : sustr
               call exchange_u2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               sustr  = tab_t_tmp / rho0
#endif   
            endif
         end if
!-------------------------------------------------
         if (srcv_clname(i)=='RRMTAUY'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_tauy) 
            if (iok_tauy /= 0) then 
#if (!defined OA_GRID_UV) 
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TAUY on grid T --> U
               do jo=jmin_halo,jmax_halo+1
               do io=imin_halo-1,imax_halo+1
                  tab_v_tmp(io,jo)=0.5*(tab_t_tmp(io,jo-1)+tab_t_tmp(io,jo))
               enddo
               enddo
!              Transmit to CROCO field : svstr
               svstr = tab_v_tmp(:,:) / rho0
               if (iok_taum == 0) tab_v_tmp = tab_t_tmp ! temporary storage...
#else
!              Transmit to CROCO field : svstr
               call exchange_v2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               svstr = tab_t_tmp(:,:) / rho0                
#endif 
            endif
         end if
!-------------------------------------------------
         if (srcv_clname(i)=='RRMTAUM'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_taum) 
            if (iok_taum /= 0) then  
!              Transmit to CROCO field : smstr
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               smstr = tab_t_tmp / rho0
!     avoid slight negative values (after interpolation for example)...
               where ( smstr < 0. )  
                  smstr = 0.
               end where
            endif
         end if
#endif
#ifdef OW_COUPLING
!        print *,'entering get ow_coupling'
!        print *,'srcv_clname(i)=',srcv_clname(i)
!-------------------------------------------------
!------ OW_COUPLING - WAVE VARIABLES -------------
!-------------------------------------------------
!        mask to manage wet-drying
!        not tested/performed yet
         if (srcv_clname(i)=='RRM_DRY'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wdry)
!            if (iok_wdry /= 0)  then
!!              Transmit to CROCO field : stflx(:,:,itemp)
!               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!               stflx(:,:,itemp) = tab_t_tmp / (rho0*Cp)
!            endif
         end if
!-------------------------------------------------
!        wave_t0m1 : mean period (s)
         if (srcv_clname(i)=='RRMT0M1'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wt0m1)
            if (iok_wt0m1 /= 0)  then
!              Transmit to CROCO field : wfrq(:,:)
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               wfrq = rmask * (2*pi) / max(tab_t_tmp, eps)
!     gc
               where ( wfrq > 1000. )  
                 wfrq = 0.
               end where
!     End gc
            endif
         end if
!-------------------------------------------------
!        wave_hs : HS significant wave height (m) 
         if (srcv_clname(i)=='RRM__HS'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_whs)
            if (iok_whs /= 0)  then
!              Transmit to CROCO field : whrm(:,:)
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               whrm = rmask * tab_t_tmp * 0.70710678
            endif
         end if
!-------------------------------------------------
!     wave_dir :  mean direction + then we apply croco grid rotation (-angler)
!                                  and projection on xi/eta trigo croco_grid
         if (srcv_clname(i)=='RRM_DIR'//clevel) then
           call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &          tab_t_tmp, rcvtime, iok_wdir)
           if (iok_wdir /= 0)  then
!     Transmit to CROCO field : wdr then compute wdrx(:,:) + wdre(:,:)
             call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
             wdrx = rmask * cos(tab_t_tmp - angler)
             wdre = rmask * sin(tab_t_tmp - angler)
             wave_dir = tab_t_tmp 
           endif
         endif
!-------------------------------------------------
!        wave_bhd : wave-induced Bernoulli head pressure (bhd in N.m-1)
!        not tested/performed yet
         if (srcv_clname(i)=='RRM_BHD'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wbhd)
!            if (iok_wbhd /= 0)  then
!!              Transmit to CROCO field : stflx(:,:,itemp)
!               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!               stflx(:,:,itemp) = tab_t_tmp / (rho0*Cp)
!            endif
         endif
!-------------------------------------------------
# ifdef OW_COUPLING_FULL
!        ORBITAL BOTTOM VELOCITY      
!        wave_ubr : rms amplitude of orbital bottom velocity of the waves (m/s)
!!         if (srcv_clname(i)=='RRM_UBR'//clevel) then
!!            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
!!     &  tab_t_tmp, rcvtime, iok_wubr)
!!            if (iok_wubr /= 0)  then
!!!              Transmit to CROCO field : ubr
!!               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!!               ubr = tab_t_tmp
!!            endif
!!         endif
!        wave_ubrx : x-component rms amplitude of wave orbital bottom velocity(m/s)
         if (srcv_clname(i)=='RRM_UBRX'//clevel) then
             call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wubrx)
            if (iok_wubrx /= 0)  then
!              Transmit to CROCO field : ubr
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               ubr = tab_t_tmp
            endif
         endif
!        wave_ubry : y-component rms amplitude of wave orbital bottom velocity(m/s)
         if (srcv_clname(i)=='RRM_UBRX'//clevel) then
             call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wvbry)
            if (iok_wvbry /= 0)  then
!              Transmit to CROCO field : vbr
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               ubr = sqrt(ubr**2+tab_t_tmp**2)
            endif
         endif
# endif         
!
!-------------------------------------------------
!        wave_phioc : Wave-to-ocean TKE flux (phioc in W.m-2)
!        not tested/performed yet
         if (srcv_clname(i)=='RRM_FOC'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wfoc)
!            if (iok_wfoc /= 0)  then
!!              Transmit to CROCO field : stflx(:,:,itemp)
!               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!               stflx(:,:,itemp) = tab_t_tmp / (rho0*Cp)
!            endif
         endif
!-------------------------------------------------
!        wave_wlm : mean length wave (m)
!        not tested/performed yet
         if (srcv_clname(i)=='RRM__LM'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wlm)
!            if (iok_wlm /= 0)  then
!!              Transmit to CROCO field : stflx(:,:,itemp)
!               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!               stflx(:,:,itemp) = tab_t_tmp / (rho0*Cp)
!            endif
         endif
!-------------------------------------------------
!        wave_twox : x-component of the wave-to-ocean momentum flux (twox in m2.s-2)
         if (srcv_clname(i)=='RRMTWOX'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wtwox)
            if (iok_wtwox /= 0) then
#if (!defined OA_GRID_UV)  
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TWOX on grid T --> U
               do jo=jmin_halo-1,jmax_halo+1
               do io=imin_halo,imax_halo+1
                  twox(io,jo)=0.5*(tab_t_tmp(io-1,jo)+tab_t_tmp(io,jo))
               enddo
               enddo
#else
               call exchange_u2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               twox  = tab_t_tmp
#endif   
            endif
         endif
!-------------------------------------------------
!        wave_twoy : y-component of the wave-to-ocean momentum flux (twoy in m2.s-2)
         if (srcv_clname(i)=='RRMTWOY'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wtwoy)
            if (iok_wtwoy /= 0) then
#if (!defined OA_GRID_UV) 
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TWOY on grid T --> U
               do jo=jmin_halo,jmax_halo+1
               do io=imin_halo-1,imax_halo+1
                  twoy(io,jo)=0.5*(tab_t_tmp(io,jo-1)+tab_t_tmp(io,jo))
               enddo
               enddo
#else
               call exchange_v2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               twoy = tab_t_tmp
#endif 
            endif
         endif
!-------------------------------------------------
!        wave_tawx : Net wave-supported stress, x-component (tawx in m2.s-2)
!        (ie atmosphere-to-wave stress)
         if (srcv_clname(i)=='RRMTAWX'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wtaux)
            if (iok_wtaux /= 0) then
#if (!defined OA_GRID_UV)  
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TWOX on grid T --> U
               do jo=jmin_halo-1,jmax_halo+1
               do io=imin_halo,imax_halo+1
                  tawx(io,jo)=0.5*(tab_t_tmp(io-1,jo)+tab_t_tmp(io,jo))
               enddo
               enddo
#else
               call exchange_u2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               tawx = tab_t_tmp
#endif   
            endif
         endif
!-------------------------------------------------
!        wave_tawy : Net wave-supported stress, y-component (tawy in m2.s-2)
!        (ie atmosphere-to-wave stress)
         if (srcv_clname(i)=='RRMTAWY'//clevel) then
            call cpl_prism_getvar(srcv_nid(:,i), cplmsk, nmaxatm,
     &  tab_t_tmp, rcvtime, iok_wtauy)
            if (iok_wtauy /= 0) then
#if (!defined OA_GRID_UV) 
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TWOY on grid T --> U
               do jo=jmin_halo,jmax_halo+1
               do io=imin_halo-1,imax_halo+1
                  tawy(io,jo)=0.5*(tab_t_tmp(io,jo-1)+tab_t_tmp(io,jo))
               enddo
               enddo
#else
               call exchange_v2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               tawy = tab_t_tmp
#endif 
            endif
         endif
#endif
      END DO

#ifdef OW_COUPLING
!-------------------------------------------------
!    Update CROCO wind stress accounting for wave bulk
!-------------------------------------------------
# if defined WAVE_SMFLUX
         sustr = twox
         svstr = twoy
# else
         sustr = sustr - tawx + twox
         svstr = svstr - tawy + twoy
# endif
#endif 
!-------------------------------------------------
!
#ifdef OA_COUPLING
      if  (iok_itemp /= 0) then
         stflx(:,:,itemp) = stflx(:,:,itemp) + srflx
      endif
      IF  (iok_taum == 0 .and. iok_taux /= 0 .and. iok_tauy /= 0) THEN
#if (!defined OA_GRID_UV) 
         smstr = sqrt(tab_u_tmp**2 + tab_v_tmp**2)  / rho0
#else
         case not coded...
#endif 
      ENDIF
#endif
      
!	
#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
#endif 
      END SUBROUTINE cpl_prism_get
