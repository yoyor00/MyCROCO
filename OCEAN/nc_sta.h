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
! This is include file "nc_sta.h".
! ==== == ======= ==== ============
!
! stafield     Number of station fields for output
! wrtsta       Logical vector with flags for output
! indxsta[...] Index of logical flag to output several fields
!       Grd  - grid level
!       Temp - temp
!       Salt - Salt
!       Rho  - Density
!       Vel  - u and v components
! ncidsta      id of station output file
! nrecsta      step to output station data
! sta[...]     several reference names of netcdf output
! staname      station output filename
! staposname   station input data filename

      integer stafield
      parameter(stafield=5)
      integer indxstaGrd, indxstaTemp, indxstaSalt
      integer indxstaRho, indxstaVel
      parameter (indxstaGrd=1, indxstaTemp=2, indxstaSalt=3)
      parameter (indxstaRho=4, indxstaVel=5)


      integer ncidsta,    nrecsta,    staGlevel
      integer staTstep,   staTime,    staXgrd,   staYgrd
      integer staZgrd,    staZeta,    staU,      staV
#ifdef SPHERICAL
      integer staLon,     staLat
#else
      integer staX,       staY
#endif
#ifdef SOLVE3D
      integer staDepth,   staDen
# ifdef TEMPERATURE
      integer staTemp
# endif
# ifdef SALINITY
      integer staSal
# endif
# ifdef MUSTANG
      integer staMUS(NT-2)
# endif
#endif
      logical wrtsta(stafield)

      common/incscrum_sta/ncidsta,    nrecsta,    staGlevel
      common/incscrum_sta/staTstep,   staTime,    staXgrd,   staYgrd
      common/incscrum_sta/staZgrd,    staZeta,    staU,      staV
#ifdef SPHERICAL
      common/incscrum_sta/staLon,     staLat
#else
      common/incscrum_sta/staX,       staY
#endif
#ifdef SOLVE3D
      common/incscrum_sta/staDepth,   staDen
# ifdef TEMPERATURE
      common/incscrum_sta/staTemp
# endif
# ifdef SALINITY
      common/incscrum_sta/staSal
# endif
# ifdef MUSTANG
      common/incscrum_sta/staMUS
# endif
#endif
      common/incscrum_sta/wrtsta

