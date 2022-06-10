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

#include "coupleur_define_BIOLink.h" 
                                     ! Equivalence of names between the hydro and the different
                                     ! Tracer and biological model. Also contains a function
                                     ! For the allocation of the main tables

   USE module_BIOLink ! script that groups all the files of BIOLink together and allows the 
                      ! access to their functions/subroutines/variables
   USE comsubstance   ! Access to the functions/variables of SUBSTANCE ( from the MUSTANG
                      ! sediment model)

   USE comBIOLink     ! Common variables of the BIOLink coupler

   USE comBIOLink_helping ! Common variables for the helping functions
                          ! of BIOLink

   USE comBIOLink_physics ! Common variables for the physical functions
                          ! of BIOLink

#if defined MUSTANG

   USE comMUSTANG ,  ONLY : htot ! Height of the water column

#endif /*MUSTANG*/

#if defined ECO3M   
   USE mod_eco3m, ONLY : THICKLAYERWC

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

      DO i=ifirst,ilast ! And on the entire zonal direction

#if ! defined MUSTANG

         
           TOTAL_WATER_HEIGHT(i,j)=z_w(i,j,N)+h(i,j) ! The total water 
                                                     ! thickness is computed
                                                     ! from the top of the 
                                                     ! last vertical layer


#endif /* MUSTANG */
#if !defined ECO3M  /* No need to calculate WW for ECO3M */ 
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
#else
          DO k=1,NB_LAYER_WAT-1          
                THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1) 
          ENDDO
          k=NB_LAYER_WAT ! Last layer
          THICKLAYERWC(k,i,j)=z_w(i,j,k)-z_w(i,j,k-1)
#endif /* ECO3M */

      END DO

   END DO
!$OMP END DO

 
  END SUBROUTINE  BIOLink_water_column

#endif /* BIOLink */

END MODULE coupleur_BIOLink_physics

