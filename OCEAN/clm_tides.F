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
!
#include "cppdefs.h"
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      subroutine clm_tides (tile)
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
      call clm_tides_tile (Istr,Iend,Jstr,Jend,
     &                     A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd),
     &                     A2d(1,4,trd),A2d(1,5,trd),A2d(1,6,trd),
     &                     A2d(1,7,trd),A2d(1,8,trd))
      return
      end
!
!*********************************************************************
      subroutine clm_tides_tile (Istr,Iend,Jstr,Jend,Cangle,Sangle,
     &                           Cphase,Sphase,Etide,Utide,Vtide,
     &                           Ptide_loc)
!*********************************************************************
!
# ifdef TIDES_MAS
      USE mod_tides_mas
# endif
      USE croco_namelist, only:dt
      implicit none
# ifndef TIDES_MAS
#  include "param.h"
# endif
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "tides.h"
# include "boundary.h"
!
# ifdef TIDES_MAS
      INTEGER               :: ian,imois,ijour
      real                  :: tempis
      logical               :: new_day
# endif
      integer
     &     Iend, Istr, Jend, Jstr, i, j, itide
      real
     &     angle, phase, omega_loc, ramp, cff
      real
     &        Cangle(PRIVATE_2D_SCRATCH_ARRAY),
     &        Cphase(PRIVATE_2D_SCRATCH_ARRAY),
     &        Sangle(PRIVATE_2D_SCRATCH_ARRAY),
     &        Sphase(PRIVATE_2D_SCRATCH_ARRAY),
     &        Etide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Utide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Vtide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Ptide_loc(PRIVATE_2D_SCRATCH_ARRAY)
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
      ramp=TANH(dt*sec2day*float(iic-1))
# endif
      cff=2.*pi*time
!
# ifdef TIDES_MAS
!-------------
! MAS INIT
!-------------
      IF (FIRST_RST_TIME_STEP) THEN
!$acc kernels if(compute_on_device) default(present)
#  if defined SSH_TIDES
        call mas_init_ssh(Istr,Iend,Jstr,Jend,IstrR,IendR,JstrR,JendR,
     &                     Ntides,k_of,nn_of,
     &                     SSH_Tamp,SSH_Tphase,mrssh,mgssh,mhssh)
#  endif
#  if defined UV_TIDES
        call mas_init_uv(Istr,Iend,Jstr,Jend,IstrR,IendR,JstrR,JendR,
     &                    Ntides,k_of,nn_of,
     &                    UV_Tmajor,UV_Tminor,UV_Tangle,UV_Tphase,
     &                    angler,mru,mgu,mrv,mgv,mhu,mhv)
#  endif
#  if defined POT_TIDES
        call mas_init_pot(Istr,Iend,Jstr,Jend,IstrR,IendR,JstrR,JendR,
     &                     Ntides,k_of,nn_of,
     &                     POT_Tamp,POT_Tphase,mrpot,mgpot,mhpot)
#  endif
!$acc end kernels
      ENDIF ! FIRST TIME STEP
# endif /* MAS INIT */

# if defined TIDES_MAS
      ijour  = mas_ijour
      imois  = mas_imois
      ian    = mas_ian
      tempis = mas_tempis
      new_day = mas_is_new_day
!$acc kernels if(compute_on_device) default(present)
#  if defined SSH_TIDES
      call mas_daily_update(IstrR,IendR,JstrR,JendR,new_day,
     &                       ijour,imois,ian,mrssh,mgssh,mhssh)
      call mas_interpolate(IstrR,IendR,JstrR,JendR,
     &                      Istr,Iend,Jstr,Jend,tempis,mhssh,Etide)
#  endif
#  if defined UV_TIDES
      call mas_daily_update(Istr,IendR,JstrR,JendR,new_day,
     &                       ijour,imois,ian,mru,mgu,mhu)
      call mas_daily_update(IstrR,IendR,Jstr,JendR,new_day,
     &                       ijour,imois,ian,mrv,mgv,mhv)
      call mas_interpolate(IstrR,IendR,JstrR,JendR,
     &                      Istr,Iend,Jstr,Jend,tempis,mhu,Utide)
      call mas_interpolate(IstrR,IendR,JstrR,JendR,
     &                      Istr,Iend,Jstr,Jend,tempis,mhv,Vtide)
#  endif
#  if defined POT_TIDES
      call mas_daily_update(IstrR,IendR,JstrR,JendR,new_day,
     &                       ijour,imois,ian,mrpot,mgpot,mhpot)
      call mas_interpolate(IstrR,IendR,JstrR,JendR,
     &                      Istr,Iend,Jstr,Jend,tempis,mhpot,Ptide_loc)
#  endif
#  if defined SSH_TIDES && defined MASKING
      do j=JRANGE
        do i=IRANGE
          Etide(i,j)=Etide(i,j)*ramp*rmask(i,j)
        enddo
      enddo
#  endif
#  if defined UV_TIDES && defined MASKING
      do j=JstrR,JendR
        do i=Istr,IendR
          Utide(i,j)=Utide(i,j)*ramp*umask(i,j)
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          Vtide(i,j)=Vtide(i,j)*ramp*vmask(i,j)
        enddo
      enddo
#  endif
#  if defined POT_TIDES && defined MASKING
      do j=JstrR,JendR
        do i=IstrR,IendR
          Ptide_loc(i,j)=Ptide_loc(i,j)*ramp*rmask(i,j)
        enddo
      enddo
#  endif
!$acc end kernels
# endif /* TIDES_MAS */

# if defined SSH_TIDES
#  ifndef TIDES_MAS
!$acc kernels if(compute_on_device) default(present)
!---------------------------------------------------------------------
!  Add tidal elevation (m) to sea surface height climatology.
!---------------------------------------------------------------------
!
      do j=JRANGE
        do i=IRANGE
          Etide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=JRANGE
            do i=IRANGE
              Etide(i,j)=Etide(i,j)+ramp*
     &                   SSH_Tamp(i,j,itide)*
     &                   COS(omega_loc-SSH_Tphase(i,j,itide))
#   ifdef MASKING
              Etide(i,j)=Etide(i,j)*rmask(i,j)
#   endif
            enddo
          enddo
        endif
      enddo
!$acc end kernels
#  endif 
!$acc kernels if(compute_on_device) default(present)
!----------------------------------------------------------------
!  Add sub-tidal forcing and adjust climatology to include tides.
!----------------------------------------------------------------
#  if defined ZCLIMATOLOGY
      do j=JRANGE
        do i=IRANGE
          ssh(i,j)=ssh(i,j)+Etide(i,j)
#   ifdef WET_DRY
          if (ssh(i,j) .lt. Dcrit(i,j)-h(i,j)) then
            ssh(i,j)=Dcrit(i,j)-h(i,j)
          endif
#   endif
        enddo
      enddo
#  endif /* ZCLIMATOLOGY */
!----------------------------------------------------------------
!  If appropriate, load tidal forcing into boundary arrays.
!----------------------------------------------------------------
#  ifdef Z_FRC_BRY
#   if defined OBC_WEST
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          zetabry_west(j)=zetabry_west(j)+Etide(Istr-1,j)
#    ifdef WET_DRY
          if (zetabry_west(j) .lt. Dcrit(Istr-1,j)-h(Istr-1,j)) then
            zetabry_west(j)=Dcrit(Istr-1,j)-h(Istr-1,j)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_EAST
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          zetabry_east(j)=zetabry_east(j)+Etide(Iend+1,j)
#    ifdef WET_DRY
          if (zetabry_east(j) .lt. Dcrit(Iend+1,j)-h(Iend+1,j)) then
            zetabry_east(j)=Dcrit(Iend+1,j)-h(Iend+1,j)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
          zetabry_south(i)=zetabry_south(i)+Etide(i,Jstr-1)
#    ifdef WET_DRY
          if (zetabry_south(i) .lt. Dcrit(i,Jstr-1)-h(i,Jstr-1)) then
            zetabry_south(i)=Dcrit(i,Jstr-1)-h(i,Jstr-1)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=IstrR,IendR
          zetabry_north(i)=zetabry_north(i)+Etide(i,Jend+1)
#    ifdef WET_DRY
          if (zetabry_north(i) .lt. Dcrit(i,Jend+1)-h(i,Jend+1)) then
            zetabry_north(i)=Dcrit(i,Jend+1)-h(i,Jend+1)
          endif
#    endif
        enddo
      endif
#   endif
#  endif /* Z_FRC_BRY */
!$acc end kernels
# endif /* SSH_TIDES */
!
# if defined UV_TIDES
!$acc kernels if(compute_on_device) default(present)
!
!---------------------------------------------------------------------
!  Add tidal currents (m/s) to 2D momentum climatologies.
!---------------------------------------------------------------------
!
#  ifndef TIDES_MAS
      do j=JstrR,JendR
        do i=Istr,IendR
          Utide(i,j)=0.
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          Vtide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=Jstr-1,JendR
            do i=Istr-1,IendR
              angle=UV_Tangle(i,j,itide)
#  ifdef CURVGRID
     &                      -angler(i,j)
#  endif
              phase=omega_loc-UV_Tphase(i,j,itide)
              Cangle(i,j)=COS(angle)
              Cphase(i,j)=COS(phase)
              Sangle(i,j)=SIN(angle)
              Sphase(i,j)=SIN(phase)
            enddo
          enddo
          do j=JstrR,JendR
            do i=Istr,IendR
              Utide(i,j)=Utide(i,j)+ramp*
     &                        0.125*((UV_Tmajor(i-1,j,itide)+
     &                                UV_Tmajor(i  ,j,itide))*
     &                               (Cangle(i-1,j)+Cangle(i,j))*
     &                               (Cphase(i-1,j)+Cphase(i,j))-
     &                               (UV_Tminor(i-1,j,itide)+
     &                                UV_Tminor(i  ,j,itide))*
     &                               (Sangle(i-1,j)+Sangle(i,j))*
     &                               (Sphase(i-1,j)+Sphase(i,j)))
#  ifdef MASKING
             Utide(i,j)=Utide(i,j)*umask(i,j)
#  endif
            enddo
          enddo
          do j=Jstr,JendR
            do i=IstrR,IendR
              Vtide(i,j)=Vtide(i,j)+ramp*
     &                        0.125*((UV_Tmajor(i,j-1,itide)+
     &                                UV_Tmajor(i,j  ,itide))*
     &                               (Sangle(i,j-1)+Sangle(i,j))*
     &                               (Cphase(i,j-1)+Cphase(i,j))+
     &                               (UV_Tminor(i,j-1,itide)+
     &                                UV_Tminor(i,j  ,itide))*
     &                               (Cangle(i,j-1)+Cangle(i,j))*
     &                               (Sphase(i,j-1)+Sphase(i,j)))
#  ifdef MASKING
              Vtide(i,j)=Vtide(i,j)*vmask(i,j)
#  endif
            enddo
          enddo
        endif
      enddo
#  endif 

#  if defined M2CLIMATOLOGY
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
      do j=JstrR,JendR
        do i=Istr,IendR
          ubclm(i,j)=ubclm(i,j)+Utide(i,j)
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          vbclm(i,j)=vbclm(i,j)+Vtide(i,j)
        enddo
      enddo
#  endif
!
!  If appropriate, load tidal forcing into boundary arrays.
!
#  ifdef M2_FRC_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) THEN
        do j=JstrR,JendR
          ubarbry_west(j)=ubarbry_west(j)+Utide(Istr,j)
        enddo
        do j=Jstr,JendR
          vbarbry_west(j)=vbarbry_west(j)+Vtide(Istr-1,j)
        enddo
      endif
#   endif
#   ifdef OBC_EAST
      if (EASTERN_EDGE) THEN
        do j=JstrR,JendR
          ubarbry_east(j)=ubarbry_east(j)+Utide(Iend+1,j)
        enddo
        do j=Jstr,JendR
          vbarbry_east(j)=vbarbry_east(j)+Vtide(Iend+1,j)
        enddo
      endif
#   endif
#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) THEN
        do i=Istr,IendR
          ubarbry_south(i)=ubarbry_south(i)+Utide(i,Jstr-1)
        enddo
        do i=IstrR,IendR
          vbarbry_south(i)=vbarbry_south(i)+Vtide(i,Jstr)
        enddo
      endif
#   endif
#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) THEN
        do i=Istr,IendR
          ubarbry_north(i)=ubarbry_north(i)+Utide(i,Jend+1)
        enddo
        do i=IstrR,IendR
          vbarbry_north(i)=vbarbry_north(i)+Vtide(i,Jend+1)
        enddo
      endif
#   endif
#  endif /* M2_FRC_BRY */
!$acc end kernels
# endif /* UV_TIDES */

# if defined POT_TIDES
!$acc kernels if(compute_on_device) default(present)
#  ifndef TIDES_MAS
      do j=JstrR,JendR
        do i=IstrR,IendR
          Ptide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=JstrR,JendR
            do i=IstrR,IendR
              Ptide(i,j)=Ptide(i,j)+
     &                   ramp*POT_Tamp(i,j,itide)*
     &                   COS(omega_loc-POT_Tphase(i,j,itide))
#   ifdef MASKING
              Ptide(i,j)=Ptide(i,j)*rmask(i,j)
#   endif
            enddo
          enddo
        endif
      enddo
#  else
      do j=JstrR,JendR
        do i=IstrR,IendR
          Ptide(i,j)=Ptide_loc(i,j)
        enddo
      enddo
#  endif

!$acc end kernels
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Ptide(START_2D_ARRAY))
#  endif

# endif /* POT_TIDES */

#else
      subroutine clm_tides_empty
#endif /* SSH_TIDES || UV_TIDES */
      return
      end
