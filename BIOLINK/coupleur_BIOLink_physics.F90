#include "cppdefs.h"
!---------------------------------------------------------------------------
!
                     MODULE coupleur_BIOLink_physics
!
!---------------------------------------------------------------------------
   

#if defined BIOLink

   !&E======================================================================
   !&E                   ***  MODULE  BIOLink_physics  ***
   !&E
   !&E ** Purpose : Functions for the necessary physics that is not handled
   !&E              by the hydro model
   !&E
   !&E ** Description :
   !&E     subroutine BIOLink_water_column     ! evaluation of total water height and 
   !&E                                          vertical meshes thickness - 
   !&E                                          called by BIOLink_update
   !&E
   !&E   History :
   !&E    !  2022 (G. Koenig) Creation from coupleur_BIOLink.F90
   !&E
   !&E=========================================================================================

!*****************************************************************************************!
!                                                                                         !
!                                      Interface                                          !
!                                                                                         !
!*****************************************************************************************!

#include "coupleur_define_BIOLink.h" ! Equivalence of names between the hydro and the different
                                     ! Tracer and biological model. Also contains a function
                                     ! For the allocation of the main tables

#ifdef key_MARS
#include "coupleur_dimhydro_BIOLink.h"
   USE sflxatm,      ONLY : rad
#else
   USE module_BIOLink ! script that groups all the files of BIOLink together and allows the 
                      ! access to their functions/subroutines/variables
   USE comsubstance   ! Access to the functions/variables of SUBSTANCE ( from the MUSTANG
                      ! sediment model)
#endif /* key_MARS */

   USE comBIOLink     ! Common variables of the BIOLink coupler

   USE comBIOLink_helping ! Common variables for the helping functions
                          ! of BIOLink

   USE comBIOLink_physics ! Common variables for the physical functions
                          ! of BIOLink

#if defined MUSTANG
   USE comMUSTANG ,  ONLY : htot ! Height of the water column
#endif /*MUSTANG*/

#if defined ECO3M   
   USE COUPLEUR_PHY_BIO ! Internal functions of coupling of ECO3M
#endif /* ECO3M */

 
   IMPLICIT NONE
  
!*****************************************************************************************!
!                                                                                         !
!                                  Accessibility                                          !
!                                                                                         !
!*****************************************************************************************!

     !====================================================================
     ! Declarations of functions
     !====================================================================

     PUBLIC BIOLink_water_column ! Function to determine the thickness of 
                                 ! the water columns and of its cells

  CONTAINS

  SUBROUTINE BIOLink_water_column(ifirst,ilast,jfirst,jlast)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_water_column ***
  !&E
  !&E ** Purpose : Computation of total water thickness and thickness of 
  !&E              each vertical layer
  !&E
  !&E ** Description : We read the depth and water elevation from the 
  !&E                  hydro model and use it to compute the total water
  !&E                  thickness. Then we use the difference of height
  !&E                  of each cell to compute the thickness of each cell
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : limitation of the subgrid given by MPI ifirst,
  !&E                     ilast, jfirst and jlast
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : Creation by B. Thouvenin ( date unknown)
  !&E              Commenting by G. Koenig ( february 2022)
  !&E
  !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined key_MARS

   USE toolgeom,     ONLY : f_dzu,f_dzw
   USE comvars2d,    ONLY : ig,id,jb,jh,hm

#endif /* key_MARS */

     !====================================================================
     ! External arguments
     !====================================================================

   INTEGER, INTENT(IN)      :: ifirst,ilast,jfirst,jlast
 
     !====================================================================
     ! Local declarations of variables
     !====================================================================

   INTEGER                  :: i,j,k,kmaxmod

    
     !====================================================================
     ! Execution of the function
     !====================================================================

      !******************* Water height computation ************************!

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod)

   DO j=jfirst,jlast ! The loop is on the entire meridional direction

#if defined key_MARS

      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1) ! I do not know

        IF(j.GE.jb(i)+1 .AND. j .LE. jh(i)-1) THEN
#else

      DO i=ifirst,ilast ! And on the entire zonal direction

#endif /* key_MARS */

#if ! defined MUSTANG

#  if defined key_MARS

          TOTAL_WATER_HEIGHT(i,j)=BATHY_H0(i,j)+                 & ! The to
                                  WATER_ELEVATION(i,j) ! tal water thickness 
                                                       ! is computed from 
                                                       ! the depth
          
#  else

           TOTAL_WATER_HEIGHT(i,j)=z_w(i,j,N)+h(i,j) ! The total water 
                                                     ! thickness is computed
                                                     ! from the top of the 
                                                     ! last vertical layer


#  endif /* key_MARS */

#endif /* MUSTANG */

#if defined key_MARS

            IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN

              kmaxmod=1

            ELSE
 
              kmaxmod=NB_LAYER_WAT

            ENDIF

            IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN

                THICKLAYERWC(1,i,j)=TOTAL_WATER_HEIGHT(i,j)
                THICKLAYERWW(1,i,j)=TOTAL_WATER_HEIGHT(i,j)*0.5_rsh
                THICKLAYERWC(2:NB_LAYER_WAT,i,j)=0.0_rsh
                THICKLAYERWW(2:NB_LAYER_WAT,i,j)=THICKLAYERWW(1,i,j)

            ELSE

              DO k=1,kmaxmod          

                THICKLAYERWC(k,i,j)=f_dzu(BATHY_H0(i,j),WATER_ELEVATION(i,j),k,i,j)
                THICKLAYERWW(k,i,j)=f_dzw(BATHY_H0(i,j),WATER_ELEVATION(i,j),k,i,j)

              ENDDO

            ENDIF

        ENDIF

#else

          DO k=1,NB_LAYER_WAT-1          

                THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1) ! The thickness
                THICKLAYERWW(k,i,j)=z_r(i,j,k+1)-z_r(i,j,k) ! of each cell
                                                            ! is computed by
                                                            ! the height 
                                                            ! difference of 
                                                            ! the cells,
                                                            ! both from top
                                                            ! and center

          ENDDO

          k=NB_LAYER_WAT ! Last layer
          THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1)
          THICKLAYERWW(k,i,j)=0._rsh

#endif /* key_MARS */

      ENDDO
   ENDDO
!$OMP END DO

 
  END SUBROUTINE  BIOLink_water_column

#endif /* BIOLink */

END MODULE coupleur_BIOLink_physics

