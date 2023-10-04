MODULE OBSTRUCTIONS

#include "cppdefs.h"

#ifdef OBSTRUCTION
   !!==========================================================================
   !!                   ***  MODULE  OBSTRUCTIONS  ***
   !!
   !! ** Purpose : obstructions interactions with hydrodynamics,
   !!              link from hydro model to OBSTRUCTIONS1DV module
   !!
   !!              Strong dependancy with hydro model for mesh size, 
   !!              position of u,v, name of variable, order of indexes i,j,k...
   !!
   !! ** Description :
   !!     
   !!     subroutine OBSTRUCTIONS_update   ! Main subroutine controlling updates 
   !!                                       of obstructions characteristics 
   !!==========================================================================

   !! * Modules used
   USE comOBSTRUCTIONS
   USE OBSTRUCTIONS1DV

   IMPLICIT NONE

   !! * Accessibility
   PUBLIC OBSTRUCTIONS_update

   PRIVATE

   CONTAINS

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_update(limin, limax, ljmin, ljmax)
   !!---------------------------------------------------------------------
   !!                 ***  ROUTINE OBSTRUCTIONS_update  ***
   !!
   !! ** Purpose : Update of obstruction parameters at each time step
   !!---------------------------------------------------------------------
   !! * Modules used
   USE initOBSTRUCTIONS, ONLY : OBSTRUCTIONS_readfile_char
   USE module_OBSTRUCTIONS  ! needed for zob, h, u,v, cm0, Zt_avg1
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

   !! * Local declaration
   TYPE(output_obst) :: obst_output_ij 
   REAL(KIND=rsh)            :: hwat
   REAL(KIND=rsh),DIMENSION(obst_kmax)   :: obst_uz
   REAL(KIND=rsh),DIMENSION(obst_kmax)   :: obst_vz
   REAL(KIND=rsh),DIMENSION(obst_kmax)   :: obst_zc
   REAL(KIND=rsh),DIMENSION(obst_kmax)   :: obst_dz
   INTEGER :: i, j, k

   !! * Executable part

   ! ******************************
   ! * READING CHARACTERISTICS FILE
   ! ******************************
   CALL OBSTRUCTIONS_readfile_char(limin, limax, ljmin, ljmax)


   DO j = ljmin, ljmax
      DO i = limin, limax
         ! **************************************
         ! * COMPUTES 3D VELOCITIES AT RHO POINT
         ! **************************************
         obst_uz(:) = (u(i,j,:,nstp) + u(i+1,j,:,nstp)) / 2.
         obst_vz(:) = (v(i,j,:,nstp) + v(i,j+1,:,nstp)) / 2.
         ! ***********************************************
         ! * COMPUTES CELL THICKNESS AND HEIGHT AT CENTER
         ! ***********************************************
         obst_dz(:) = Hz(i,j,:)
         obst_zc(1) = Hz(i,j,1) / 2.
         DO k = 2, obst_kmax
            obst_zc(k) = obst_zc(k-1) + Hz(i,j,k-1) / 2. + Hz(i,j,k) / 2.
         ENDDO

         hwat = h(i,j) + Zt_avg1(i,j)

         ! ************************
         ! * UPDATING OBSTRUCTION 
         ! ************************
         CALL o1dv_main(hwat, cm0, zob(i,j), obst_uz, obst_vz, obst_dz, obst_zc, &
            obst_height_inst(:,i,j), obst_dens_inst(:,i,j), &
            obst_thick_inst(:,i,j), obst_width_inst(:,i,j), &
            obst_position(:,i,j), obst_height(:,i,j), &
            obst_output_ij)
         
         ! ************************
         ! * SAVE FOR OUTPUT
         ! ************************
         obst_height(:,i,j) = obst_output_ij%height(:)
         obst_dens3d(:,:,i,j) = obst_output_ij%dens(:,:)
         obst_width3d(:,:,i,j) = obst_output_ij%width(:,:)
         obst_thick3d(:,:,i,j) = obst_output_ij%thick(:,:)
         obst_fuz(i,j,:) = obst_output_ij%fuz(:)
         obst_fvz(i,j,:) = obst_output_ij%fvz(:)
         obst_a3d(:,i,j,:) = obst_output_ij%a3d(:,:)
         obst_tau(i,j,:) = obst_output_ij%tau(:)
         obst_t(i,j,:) = obst_output_ij%t(:)
      ENDDO
   ENDDO

   END SUBROUTINE OBSTRUCTIONS_update 

   !!==========================================================================================================

#endif
END MODULE OBSTRUCTIONS