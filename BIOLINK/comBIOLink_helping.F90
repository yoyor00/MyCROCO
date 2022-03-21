!---------------------------------------------------------------------------
!
                     MODULE comBIOLink_helping
!
!---------------------------------------------------------------------------
#include "cppdefs.h"

#if defined SUBSTANCE && defined BIOLink

  !!======================================================================================!!
  !!                   ***  MODULE  comBIOLink  ***                                       !!
  !! Ocean dynamics Bio :  Declare and initialize all common variables                    !!
  !!                        for BIOLink helping functions                                 !!
  !!                                                                                      !!
  !!   History : ! February 2022 : Creation from comBIOLink ( G. Koenig)                  !!          !!                                                                                      !!
  !!                                                                                      !!
  !!======================================================================================!!

#include "coupleur_define_BIOLink.h"

      !! * Modules used
      USE comsubstance

      IMPLICIT NONE

     !*************************************************************************!
     !*************************************************************************!
     !*************** Variables for the model helping routines ****************!
     !********************* Or sediment model *********************************!
     !*************************************************************************!
     !*************************************************************************!

     !=====================================================================
     !  Sinking velocity of tracers
     !=====================================================================

      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:,:,:)             :: WS_BIOLink 
      ! Settling velocities, it depends on particles density and vertical currents

     !=====================================================================
     !  Radiation variables
     !=====================================================================

     ! Photosynthetic available radiations

#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR
      ! Photosynthetic Available Radiation to be read in CROCO
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer
      ! Photosynthetic Available Radiation in the top layer
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE :: PAR_avg_layer_phyto 
      ! Photosynthetic Available Radiation in a lyaer concerned with phytoplankton,
      ! I am not sure
#    if defined PEPTIC
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: PAR_top_layer_day 
      ! Photosynthetic Available Radiation in the top layer per day 
!      REAL(KIND=rsh),DIMENSION(:,:)  ,ALLOCATABLE :: light_ave_daily,light_integ
      ! Light averaged per day and total light ( per day ? I am not sure)
#    endif /* PEPTIC */
#  endif /* BIOLink_PAR_eval */
     
     ! Extinction and attenuation

#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: EXTINCTION_RAD
      ! Coefficient of light extinction of the water column, I am not sure
#  endif /* BIOLink_PAR_eval */

#  ifdef BLOOM
      REAL(KIND=rsh),DIMENSION(:,:,:,:,:)  ,ALLOCATABLE :: extinction_tab
      ! Table of light extinction, but I do not entirely understand the dimensions
      REAL(KIND=rsh),DIMENSION(:,:,:)  ,ALLOCATABLE :: extinction_ave4d,extinction_aveh
      ! Extinction averaged on four days and per hour
      REAL(KIND=rsh)                                :: t_cum_extinctionh
      ! Extinction cumulated on time
#  endif /* BLOOM */

     ! Variables that affect the extinction

      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: SPMTOT_MGL
      ! Concentration of suspended matter

#  if defined  BIOLink_PAR_eval
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: BIOLink_chloro
      ! Chlorophyll in the water column, I am not sure
#  endif /* BIOLink_PAR_eval */

#  ifdef BLOOM
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_diat,effetlumiere_day_dino,effetlumiere_day_nano
      ! Integration over one day of the effect of light on the phytoplankton growth ( from midnight
      !  to midnight) of diatoms, dinoflagellates and nanoplankton
#    ifdef key_phaeocystis
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_phaeocystis
      ! Effect of light on the growth of phaeocystis
#    endif /* key_phaeocystis */
#    ifdef key_karenia
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_karenia
      ! Effect of light on the growth of karenia
#    endif /* key_karenia */
#    ifdef key_psnz
      REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:,:)  ,PUBLIC :: effetlumiere_day_psnz
      ! Effect of light on the growth of psnz
#    endif /* key_psnz */
#  endif /* BLOOM */

     !=====================================================================
     !  Estimation of the suspended matter via satellite
     !=====================================================================

#  ifdef key_messat
      CHARACTER(LEN=lchain)                       :: filemessat_clim,filemessat_obs
      ! Files of climatology and observation of suspended matter 
      LOGICAL                                     :: l_messat_clim,l_messat_obs 
      ! Boolean for the activation of climatology of suspended mattter climatology and observations
      INTEGER, PUBLIC                             :: idateinf_messat,idatesup_messat
      ! Date of beginning of climatology and satellite observations
      INTEGER, PUBLIC                             :: date_s_annee_messat,idimt_messat,idimt_messat_obs
      ! Temporal variables for the satellite observations, I am note entirely sure
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat
      ! Concentration of suspended matter observed by satellite
      REAL(KIND=rlg),DIMENSION(1:51)              :: t_clim_messat
      ! Time of the satellite observations, it is one per week I guess
      REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: messat_obs
      ! Observation of suspended matter via satellite ?
      REAL(KIND=rlg),DIMENSION(:),ALLOCATABLE      :: t_obs_messat
      ! Time of observation of the satellite measurements ?
      CHARACTER(LEN=lchain)                       :: file_mes_obs_path
      ! Path to the file of observations measurements ?
#  endif /* key_messat */

     !=====================================================================
     !  Wind variables
     !=====================================================================

#  if ! defined BULK_FLUX
      REAL(kind=rsh),ALLOCATABLE,DIMENSION(:,:)                :: WIND_SPEED 
      ! Wind velocity
#  endif /* BULK_FLUX */


CONTAINS

#endif /* SUBSTANCE && BIOLink */

END MODULE
