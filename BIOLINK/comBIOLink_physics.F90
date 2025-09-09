!---------------------------------------------------------------------------
!
                     MODULE comBIOLink_physics
!
!---------------------------------------------------------------------------
#include "cppdefs.h"

#if defined SUBSTANCE && defined BIOLink

  !!======================================================================================!!
  !!                   ***  MODULE  comBIOLink  ***                                       !!
  !! Ocean dynamics Bio :  declare and initialize all common variables                    !!
  !!                       related to physical part of BIOLink module                     !!
  !!                       ( water column, velocities and density )                       !!
  !!                                                                                      !!
  !!   History : ! February 2022 : Creation from comBIOLink (G. Koenig)                   !!                                                      !!
  !!                                                                                      !!
  !!======================================================================================!!

#include "coupleur_define_BIOLink.h"

      !! * Modules used
      USE comsubstance

     !*************************************************************************!
     !*************************************************************************!
     !******************** Variables from the hydro model *********************!
     !*************************************************************************!
     !*************************************************************************!

     !=====================================================================
     !  Height of the water column variables
     !=====================================================================
#if !defined ECO3M
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE  :: THICKLAYERWC 
      ! Thickness of cells (dz), centered around C,T,Sal
#endif
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE  :: THICKLAYERWW 
      ! Thickness of cells (dz), centered around W, interface
#  if ! defined MUSTANG
      REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: TOTAL_WATER_HEIGHT 
      ! Total water height of the column
#  endif /* MUSTANG */

     !=====================================================================
     !  Bottom current variables
     !=====================================================================

#  if defined BLOOM && defined key_benthos
      REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: BOTTCURRENTBIO 
      ! Current in the bottom layer
#  endif /* BLOOM && key_benthos */

     !=====================================================================
     !  Temperature and salinity variables
     !=====================================================================
#  if !defined ECO3M
      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:)               :: SAL_BIOLink,TEMP_BIOLink
      ! Salinity and temperature in the water column
#  endif /* ECO3M */


     CONTAINS

#endif /* SUBSTANCE && BIOLink */

END MODULE
