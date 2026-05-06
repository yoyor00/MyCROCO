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
      integer indxstaGrd, indxstaTemp, indxstaSalt,
     & indxstaRho, indxstaVel
      parameter (     indxstaGrd=1, indxstaTemp=2,
     & indxstaSalt=3, indxstaRho=4,  indxstaVel=5)


      integer ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
#ifdef SPHERICAL
     &      , staLon,     staLat
#else
     &      , staX,       staY
#endif
#ifdef SOLVE3D
     &      , staDepth,   staDen
# ifdef TEMPERATURE
     &      , staTemp
# endif
# ifdef SALINITY
     &      , staSal
# endif
# ifdef MUSTANG
     &      , staMUS(NT-2)
# endif
#endif
      logical wrtsta(stafield)

      common/incscrum_sta/
     &        ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
#ifdef SPHERICAL
     &      , staLon,     staLat
#else
     &      , staX,       staY
#endif
#ifdef SOLVE3D
     &      , staDepth,   staDen
# ifdef TEMPERATURE
     &      ,   staTemp
# endif
# ifdef SALINITY
     &      , staSal
# endif
# ifdef MUSTANG
     &      , staMUS
# endif

#endif
     &      , wrtsta


      character*80  staname,   staposname
      common /cncscrum_sta/ staname,   staposname
