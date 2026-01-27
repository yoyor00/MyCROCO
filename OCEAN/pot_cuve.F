! $Id: clm_tides.F 1476 2014-02-17 08:55:42Z rblod $
!
!======================================================================
! CROCO is +a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined POT_CUVE
      subroutine pot_cuve (tile)
!
!================================================== Robert Hetland ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine adds tidal elevation (m) and tidal currents (m/s) to !
!  sea surface height and 2D momentum climatologies, respectively.   !
!                                                                    !
!=====================================================================
!
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
!
      trd=omp_get_thread_num()
      call pot_cuve_tile (Istr,Iend,Jstr,Jend)
      return
      end
!
!*********************************************************************
      subroutine pot_cuve_tile (Istr,Iend,Jstr,Jend)
!*********************************************************************
!
      implicit none
# include "param.h"
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "tides.h"
# include "boundary.h"
!
      integer
     &     Iend, Istr, Jend, Jstr, i, j, itide, a,
     &     N_plunger
      real
     &     angle, phase, omega_loc, ramp, cff, xrc,
     &     yrc, amp_p, Ttide, rad_p, t_stop, rad_p2,
     &     Ttide2, xrc2, yrc2, amp_p2
      logical :: P_stop
!
# include "compute_auxiliary_bounds.h"
# ifdef MPI
#  define IRANGE Istr-1,Iend+1
#  define JRANGE Jstr-1,Jend+1
# else
#  define IRANGE IstrR,IendR
#  define JRANGE JstrR,JendR
# endif
!
      angle=0.
      phase=0.
      ramp=1.
      
# ifdef TIDERAMP
      ramp=TANH(dt*float(iic-1)/35.)
# endif
      cff=2.*pi*time
!
c     cdate = tool_sectodat(time)
c     CALL tool_decompdate(cdate,ijour,imois,ian,iheure,iminu,isec)
c     tempis =iheure*3600.0+iminu*60.0+isec
!
!  If appropriate, load tidal forcing into boundary arrays.
!
# if defined POT_CUVE
!
!-----------------------------------------------------------------------
!  Compute tidal potential (m) 
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
!
      ! ==================
      ! --- Parameters ---
      ! ==================
      N_plunger = 1       ! Number of plunger(s)
      if (N_plunger == 1) then
         Ttide = 34.49     ! Tide period (in seconds) (a priori, same for EXP09)
         rad_p = 1.        ! Radius of the circle (in meters)
         xrc   = -6.2      !-5.70   ! x-position of the center (in meters)
         yrc   =  3.6      ! 5.52   ! y-position of the center (in meters)
         amp_p =  0.03*ramp     ! Plunger amplitude in excel file (in meters)
                           ! [Factor 2 if two plungers (e.g, spring tide)]
      else if (N_plunger == 2) then
         ! --- Plunger B1:
         Ttide = 35.77
         rad_p = 0.25
         xrc   = -5.720
         yrc   =  0.385
         amp_p =  0.03*ramp
         ! --- Plunger B2:
         Ttide2 = 35.77
         rad_p2 = 0.25
         xrc2   = -6.940
         yrc2   =  4.475 
         amp_p2 =  0.03*ramp
      end if
      ! =================
      ! --- Plunger stop:
      ! =================
      P_stop = .False.             ! P_stop: True if you want to stop the plunger
                                  !         after a given time set by T_stop
      if (P_stop == .True.) then
        t_stop = 394.             ! in seconds
      end if
      ! -----------------
      omega_loc=cff/Ttide
      ! --------------
      do j=JstrR,JendR
        a=0.0
        do i=IstrR,IendR
!     Force de gradient de pression apliquée sous la forme d'un cercle
          if (N_plunger == 1) then
             if (P_stop == .False.) then
                Ptide(i,j)=((amp_p*(1.+pi/4.))/(pi*rad_p**2.))
     &                      *exp(-((xr(i,j)-xrc)**2.+(yr(i,j)-
     &                      yrc)**2.)/rad_p**2.)*SIN(omega_loc)
             else
                if (time < t_stop) then
                   Ptide(i,j)=((amp_p*(1.+pi/4.))/(pi*rad_p**2.))
     &                        *exp(-((xr(i,j)-xrc)**2.+(yr(i,j)-
     &                        yrc)**2.)/rad_p**2.)*SIN(omega_loc)
               else if (time >= t_stop) then
                  Ptide(i,j) =0.
               end if
             end if
          else if (N_plunger == 2) then
             if (P_stop == .False.) then
                Ptide(i,j)=((amp_p*(1.+pi/4.))/(pi*rad_p**2.))
     &                      *exp(-((xr(i,j)-xrc)**2.+(yr(i,j)-
     &                      yrc)**2.)/rad_p**2.)*SIN(omega_loc)+
     &                      ((amp_p2*(1.+pi/4.))/(pi*rad_p2**2.))
     &                      *exp(-((xr(i,j)-xrc2)**2.+(yr(i,j)-
     &                      yrc2)**2.)/rad_p2**2.)*SIN(omega_loc)
             else
                if (time < t_stop) then
                   Ptide(i,j)=((amp_p*(1.+pi/4.))/(pi*rad_p**2.))
     &                        *exp(-((xr(i,j)-xrc)**2.+(yr(i,j)-
     &                        yrc)**2.)/rad_p**2.)*SIN(omega_loc)+
     &                        ((amp_p2*(1.+pi/4.))/(pi*rad_p2**2.))
     &                        *exp(-((xr(i,j)-xrc2)**2.+(yr(i,j)-
     &                        yrc2)**2.)/rad_p2**2.)*SIN(omega_loc)
               else if (time >= t_stop) then 
                  Ptide(i,j) =0.
               end if
             end if
          end if
!     Force de gradient de pression appliqué sous la forme de deux demi-cercle séparé par un rectangle orientée selon l'axe x
!         if (xr(i,j) .le. 2.6) then
!          Ptide(i,j)=0.01
!     &     *exp(-((xr(i,j)-2.6)**2+(yr(i,j)-9.62)**2)/(0.5)**2)
!     &     *COS(omega_loc)
!         else if (xr(i,j) .le. (2.6+1)) then
!          Ptide(i,j)=0.01
!     &     *exp(-((xr(i,j)-(2.6+a*0.1))**2+(yr(i,j)-(9.62))**2)/(0.5)**2)
!     &     *COS(omega_loc)
!          a=a+1.0
!         else
!          Ptide(i,j)=0.01
!     &     *exp(-((xr(i,j)-(2.6+1))**2+(yr(i,j)-(9.62))**2)/(0.5)**2)
!     &     *COS(omega_loc)
!         end if
!
#  ifdef MASKING
              Ptide(i,j)=Ptide(i,j)*rmask(i,j)
#  endif
        enddo
      enddo
c     do itide=1,Ntides
c       if (Tperiod(itide).gt.0.) then
c         omega_loc=cff/Tperiod(itide)
c         do j=JstrR,JendR
c           do i=IstrR,IendR
c             Ptide(i,j)=Ptide(i,j)+ 
c    &                   ramp*POT_Tamp(i,j,itide)*
c    &                   COS(omega_loc-POT_Tphase(i,j,itide))
c           enddo
c         enddo
c       endif
c     enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Ptide(START_2D_ARRAY))
#  endif

# endif /* POT_CUVE */

#else
      subroutine pot_cuve_empty
#endif /* SSH_TIDES || UV_TIDES */
      return
      end
