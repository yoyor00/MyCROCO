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
! Include file "init_sta.h".
! ==============================
!
! statitle         Stations application title.
! stat0,stax0, ... input start time and positions from sta.in
! stacoor          type of coordinates in input(lat,lon or x,y)a
! STgrd            input grid level station in nested applications

      real stat0(Msta), stax0(Msta), stay0(Msta), staz0(Msta)
      common /ncrealsta/ stat0, stax0, stay0, staz0

      integer  stacoor(Msta), STgrd(Msta)
      common /ncintsta/ stacoor,STgrd

! to be tested if STgrd is really necessary ALVARO

      character*80 statitle
      common /nccharsta/ statitle
