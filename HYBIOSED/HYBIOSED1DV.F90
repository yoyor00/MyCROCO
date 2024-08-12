MODULE HYBIOSED1DV
   !
   ! This module manage the interactions and feedbacks between hydrodynamics,
   ! biology and sediment dynamics (in one cell "i,j").
   !
   ! No use, no commons, no include this module is independant.
   !
   ! Authors : F.Ganthy
   ! adaptation from MARS version to an independant module S.Le Gac, 2024
   !
   ! PUBLIC SUBROUTINES, VARIABLES and TYPES names of this module are precede
   ! by "hbs1dv_"
   ! All variables are private except ***TODO***
   !

   IMPLICIT NONE
   PRIVATE

   ! VARIABLE
   PUBLIC :: hbs1dv_param
   ! SUBROUTINE
   PUBLIC :: hbs1dv_alloc
   PUBLIC :: hbs1dv_init
   PUBLIC :: hbs1dv_update_root_level
   PUBLIC :: hbs1dv_ws_trapp
   PUBLIC :: hbs1dv_ws_block
   PUBLIC :: hbs1dv_update_ws_coeff
   PUBLIC :: hbs1dv_zroot_troot_ero
   PUBLIC :: hbs1dv_comp_erosion

   ! TYPE
   PUBLIC :: hbs1dv_param_type

   ! Declaration

   ! parameters
   INTEGER, PARAMETER :: lchain = 200 ! max length for string
   INTEGER, PARAMETER :: rsh = 8 ! working precision
   REAL(KIND=rsh), PARAMETER  :: pi = atan(1.D0)*4.D0

   TYPE hbs1dv_param_type
      ! characteristic of an HYBIOSED veriable
      LOGICAL :: l_root_accomodation
      LOGICAL :: l_root_erosion
      LOGICAL :: l_bedsedstab
      INTEGER ::   root_accomod_type
      REAL(KIND=rsh) :: c_tauceroot_x0
      REAL(KIND=rsh) :: c_tauceroot_x1
      REAL(KIND=rsh) :: c_E0root_x0
      REAL(KIND=rsh) :: c_E0root_x1
      REAL(KIND=rsh) :: c_root_accomod_vel_max
      REAL(KIND=rsh) :: c_root_accomod_day
      REAL(KIND=rsh) :: c_opt_root_depth
      REAL(KIND=rsh) :: c_opt_root_thick
   END TYPE

   ! variables of the module

   INTEGER :: hbs1dv_nbvar    ! Total number of hybiosed variables
   INTEGER :: hbs1dv_kmax    ! TODO
   LOGICAL :: hbs1dv_l_suspsed_trapp
   LOGICAL :: hbs1dv_l_suspsed_block
   REAL(KIND=rsh) :: hbs1dv_h0fond
   REAL(KIND=rsh) :: hbs1dv_c_suspsed_trapp_exp
   REAL(KIND=rsh) :: hbs1dv_c_suspsed_trapp_max
   REAL(KIND=rsh) :: hbs1dv_c_suspsed_block_exp
   REAL(KIND=rsh) :: hbs1dv_c_suspsed_block_max
   TYPE(hbs1dv_param_type), DIMENSION(:), ALLOCATABLE :: hbs1dv_param
   ! Hybiosed parameters

CONTAINS
   !============================================================================
   SUBROUTINE hbs1dv_alloc(nbvar, kmax)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_alloc  ***
   !!
   !! ** Purpose : Initialize the dimensions
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      INTEGER, INTENT(IN) :: nbvar ! number of hybiosed variables
      INTEGER, INTENT(IN) :: kmax  ! number of layers in water
      !! * Executable part
      hbs1dv_nbvar = nbvar
      hbs1dv_kmax = kmax
      IF (hbs1dv_nbvar > 0) THEN
         ALLOCATE (hbs1dv_param(nbvar))
      END IF

      RETURN
   END SUBROUTINE hbs1dv_alloc

   !============================================================================

   SUBROUTINE hbs1dv_init(h0fond, hbs_l_suspsed_trapp, &
                          hbs_c_suspsed_trapp_max, hbs_c_suspsed_trapp_exp, &
                          hbs_l_suspsed_block, hbs_c_suspsed_block_max, &
                          hbs_c_suspsed_block_exp, &
                          hbs_c_tauceroot_x0, hbs_c_tauceroot_x1, &
                          hbs_c_E0root_x0, hbs_c_E0root_x1, &
                          hbs_c_root_accomod_vel_max, hbs_c_root_accomod_day, &
                          hbs_c_opt_root_depth, hbs_c_opt_root_thick, &
                          hbs_c_root_accomod_type, l_hbs_root_accomodation, &
                          l_hbs_root_erosion, l_hbs_bedsedstab)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_init  ***
   !!
   !! ** Purpose : Allocate the global arrays of the module.
   !!---------------------------------------------------------------------------
      IMPLICIT NONE

      REAL, INTENT(IN) :: h0fond
      LOGICAL, INTENT(IN) :: hbs_l_suspsed_trapp
      LOGICAL, INTENT(IN) :: hbs_l_suspsed_block
      REAL(KIND=rsh), INTENT(IN) :: hbs_c_suspsed_trapp_exp
      REAL(KIND=rsh), INTENT(IN) :: hbs_c_suspsed_trapp_max
      REAL(KIND=rsh), INTENT(IN) :: hbs_c_suspsed_block_exp
      REAL(KIND=rsh), INTENT(IN) :: hbs_c_suspsed_block_max
      LOGICAL, DIMENSION(hbs1dv_nbvar), INTENT(IN) :: l_hbs_root_accomodation
      LOGICAL, DIMENSION(hbs1dv_nbvar), INTENT(IN) :: l_hbs_root_erosion
      LOGICAL, DIMENSION(hbs1dv_nbvar), INTENT(IN) :: l_hbs_bedsedstab
      INTEGER, DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_root_accomod_type
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_tauceroot_x0
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_tauceroot_x1
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_E0root_x0
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_E0root_x1
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_root_accomod_vel_max
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_root_accomod_day
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_opt_root_depth
      REAL(KIND=rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: hbs_c_opt_root_thick

      INTEGER :: iv

      hbs1dv_h0fond = h0fond
      hbs1dv_l_suspsed_trapp = hbs_l_suspsed_trapp
      hbs1dv_l_suspsed_block = hbs_l_suspsed_block
      hbs1dv_c_suspsed_trapp_exp = hbs_c_suspsed_trapp_exp
      hbs1dv_c_suspsed_trapp_max = hbs_c_suspsed_trapp_max
      hbs1dv_c_suspsed_block_exp = hbs_c_suspsed_block_exp
      hbs1dv_c_suspsed_block_max = hbs_c_suspsed_block_max

      DO iv = 1, hbs1dv_nbvar
         hbs1dv_param(iv)%l_root_accomodation = l_hbs_root_accomodation(iv)
         hbs1dv_param(iv)%l_root_erosion = l_hbs_root_erosion(iv)
         hbs1dv_param(iv)%l_bedsedstab = l_hbs_bedsedstab(iv)
         hbs1dv_param(iv)%root_accomod_type = hbs_c_root_accomod_type(iv)
         hbs1dv_param(iv)%c_tauceroot_x0 = hbs_c_tauceroot_x0(iv)
         hbs1dv_param(iv)%c_tauceroot_x1 = hbs_c_tauceroot_x1(iv)
         hbs1dv_param(iv)%c_E0root_x0 = hbs_c_E0root_x0(iv)
         hbs1dv_param(iv)%c_E0root_x1 = hbs_c_E0root_x1(iv)
         hbs1dv_param(iv)%c_root_accomod_vel_max = hbs_c_root_accomod_vel_max(iv)
         hbs1dv_param(iv)%c_root_accomod_day = hbs_c_root_accomod_day(iv)
         hbs1dv_param(iv)%c_opt_root_depth = hbs_c_opt_root_depth(iv)
         hbs1dv_param(iv)%c_opt_root_thick = hbs_c_opt_root_thick(iv)
      END DO

      RETURN
   END SUBROUTINE hbs1dv_init
   !============================================================================

   SUBROUTINE hbs1dv_update_ws_coeff(uz, vz, s3d, a3d, ws_trapp, ws_block)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_update_ws_coeff  ***
   !!
   !! ** Purpose : TODO
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN)   :: uz
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN)   :: vz
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN)   :: s3d
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN)   :: a3d
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(INOUT) :: ws_trapp
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(INOUT) :: ws_block
      !! * Executable part
      ! *******************************
      ! * COMPUTES CORRECTION FACTOR FOR
      ! * SETTLING VELOCITIES
      ! *******************************
      ws_trapp = hbs1dv_ws_trapp(uz, vz, s3d)
      ws_block = hbs1dv_ws_block(a3d)
      RETURN
   END SUBROUTINE hbs1dv_update_ws_coeff
   !============================================================================

   FUNCTION hbs1dv_ws_trapp(uz, vz, s3d) RESULT(ws_trapp)
   !!---------------------------------------------------------------------------
   !!                 ***  SUBROUTINE hbs1dv_comp_ws_trapp  ***
   !!
   !! ** Description : Computes correction factor for settling velocities
   !!                  due to trapping by obstructions
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN) ::  uz
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN) ::  vz
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN) ::  s3d
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax) ::  ws_trapp
      !! * Local declaration
      REAL(KIND=rsh) :: uzvz
      INTEGER        :: k
      !! * Executable part
      ws_trapp(:) = 1.0_rsh
      IF (hbs1dv_l_suspsed_trapp) THEN
         DO k = 1, hbs1dv_kmax
            uzvz = SQRT(uz(k)**2.0_rsh + vz(k)**2.0_rsh)
            IF (hbs1dv_c_suspsed_trapp_exp == 0.0_rsh) THEN
               ws_trapp(k) = 1.0_rsh + hbs1dv_c_suspsed_trapp_max*s3d(k)*uzvz
            ELSE
               ws_trapp(k) = 1.0_rsh + hbs1dv_c_suspsed_trapp_max* &
                             (EXP(hbs1dv_c_suspsed_trapp_exp* &
                                  s3d(k)*uzvz) - 1.0_rsh)/ &
                             (EXP(hbs1dv_c_suspsed_trapp_exp) - 1.0_rsh)
            END IF
            ws_trapp(k) = MIN(ws_trapp(k), hbs1dv_c_suspsed_trapp_max)
         END DO
      END IF
   END FUNCTION hbs1dv_ws_trapp
   !!===========================================================================

   FUNCTION hbs1dv_ws_block(a3d) RESULT(ws_block)
   !!---------------------------------------------------------------------------
   !!                 ***  SUBROUTINE hbs1dv_comp_ws_block ***
   !!
   !! ** Description : Computes correction factor for settling velocities
   !!                  due to blockage by obstructions
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax), INTENT(IN) ::  a3d
      REAL(KIND=rsh), DIMENSION(hbs1dv_kmax) :: ws_block
      !! * Local declaration
      INTEGER        :: k
      !! * Executable part
      ws_block(:) = 1.0_rsh
      IF (hbs1dv_l_suspsed_block) THEN
         DO k = 1, hbs1dv_kmax
            IF (a3d(k) > 0.0_rsh) THEN
               IF (a3d(k) >= hbs1dv_c_suspsed_block_max) THEN
                  ws_block(k) = 0.0_rsh
               ELSE
                  IF (hbs1dv_c_suspsed_block_exp .EQ. 0.0_rsh) THEN
                     ws_block(k) = 1.0_rsh + &
                                   a3d(k)*(-1.0_rsh/hbs1dv_c_suspsed_block_max)
                  ELSE
                     ws_block(k) = 1.0_rsh - &
                                   (EXP(hbs1dv_c_suspsed_block_exp* &
                                        (a3d(k)/hbs1dv_c_suspsed_block_max)) &
                                    - 1.0_rsh) &
                                   /(EXP(hbs1dv_c_suspsed_block_exp) - 1.0_rsh)
                  END IF
               END IF ! END TEST A3D>Amax
            END IF ! END TEST A3D
         END DO ! END LOOP kmax
      END IF
   END FUNCTION hbs1dv_ws_block
   !!===========================================================================

   SUBROUTINE hbs1dv_update_root_level(zup_root0, thick_root0, &
                                       zup_root, thick_root, dz_root, &
                                       dthick_root, dhsed, &
                                       position_wat, position_bed, dt, jjulien)
   !!---------------------------------------------------------------------------
   !!                 ***  SUBROUTINE hbs1dv_update_root_level ***
   !!
   !! ** Description : Update root level depth and thickness depending on
   !!                  erosion/deposition and root accomodation
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(rsh), INTENT(IN) :: dhsed ! dhsed : effective change of sediment
      ! height, dhsed>0 for accretion , dhsed<0 for erosion
      REAL(rsh), INTENT(IN) :: dt
      INTEGER, INTENT(IN) :: jjulien
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: position_wat
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: position_bed
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: zup_root0
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: thick_root0
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: zup_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: thick_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: dz_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: dthick_root

      !! * Local declaration
      INTEGER          :: iv
      REAL(KIND=rsh)   :: dz_acc, dz_bot, dz_up, vel_tmp, dherod
      CHARACTER(LEN=lchain) :: CaseS
      TYPE(hbs1dv_param_type) :: var

      dz_root(:) = 0.0_rsh
      dthick_root(:) = 0.0_rsh
      dz_acc = 0.0_rsh
      dz_bot = 0.0_rsh
      dz_up = 0.0_rsh
      IF (dhsed .LT. 0.0_rsh) THEN
         dherod = -dhsed
      ELSE
         dherod = 0.0_rsh
      END IF
      ! ***************************************************
      ! * COMPUTES CHANGE OF ROOT LEVEL DEPTH AND THICKNESS
      ! ***************************************************
      DO iv = 1, hbs1dv_nbvar
         var = hbs1dv_param(iv)

         IF ((var%l_bedsedstab) .AND. &
             ((position_wat(iv) .GT. 0.0_rsh) .OR. &
              (position_bed(iv) .GT. 0.0_rsh))) THEN

            dz_acc = hbs1dv_accomodation_potential(var, &
                                                   position_wat(iv), &
                                                   position_bed(iv), &
                                                   dt, jjulien)

            ! *************************************** !
            IF (dhsed .LT. 0.0_rsh) THEN
               ! **************************
               ! * CASE-1 : EROSION OCCURED
               ! **************************
               IF (zup_root0(iv) .EQ. 0.0_rsh) THEN
                  ! ***************************************************** !
                  ! * CASE-1.1 : ROOT LEVEL ALREADY AT SEDIMENT SURFACE * !
                  ! ***************************************************** !
                  ! Update depth of root level related to bed erosion
                  CaseS = '1.1'
                  zup_root(iv) = zup_root0(iv) ! No change, already = 0.0
                  IF (var%l_root_erosion) THEN
                     ! *************************** !
                     ! * ROOT EROSION IS ALLOWED * !
                     ! *************************** !
                     ! Update thickness of root level related to erosion
                     thick_root(iv) = MAX(thick_root0(iv) - dherod, 0.0_rsh)
                     ! Update thickness of root level related to
                     ! downward accomodation
                     thick_root(iv) = MIN(thick_root(iv) + dz_acc, &
                                          var%c_opt_root_thick)
                  END IF
               ELSE
                  ! ***************************************************** !
                  ! * CASE-1.2 : ROOT LEVEL NOT YET AT SEDIMENT SURFACE * !
                  ! ***************************************************** !
                  ! Update first related to bed erosion
                  IF (zup_root0(iv) .GE. dherod) THEN
                     ! ****************************************************** !
                     ! * Thickness of eroded layer is lower than depth of   * !
                     ! * root level                                         * !
                     ! ****************************************************** !
                     CaseS = '1.2a'
                     zup_root(iv) = zup_root0(iv) - dherod
                     thick_root(iv) = thick_root0(iv)! No change of thickness
                  ELSE
                     ! ****************************************************** !
                     ! * Thickness of eroded layer is higher than depth of  * !
                     ! * root level                                         * !
                     ! * All surficial layer is eroded + part of root level * !
                     ! ****************************************************** !
                     CaseS = '1.2b'
                     dz_bot = dherod - zup_root0(iv) ! Thickness of roots which
                     ! will be eroded
                     zup_root(iv) = 0.0_rsh ! Roots reache sediment surface
                     IF (var%l_root_erosion) THEN
                        ! *************************** !
                        ! * ROOT EROSION IS ALLOWED * !
                        ! *************************** !
                        ! Update thickness of root level related to erosion
                        thick_root(iv) = MAX(thick_root0(iv) + dz_bot, 0.0_rsh)
                     END IF
                  END IF ! * END TEST ON hbs_zup_root0.GE.dhsed
                  ! Now update due to accomodation
                  IF (zup_root(iv) .GT. var%c_opt_root_depth) THEN
                     ! ******************************************************* !
                     ! * CASE-1.2.1 :                                        * !
                     ! * ROOT LEVEL IS DEEPER THAN OPTIMAL DEPTH             * !
                     ! * IT MAY ACCOMODATE UPWARD (AND MAYBE DOWNWARD)       * !
                     ! ******************************************************* !
                     dz_up = -(var%c_opt_root_depth - zup_root(iv))
                     ! Authorized upward accomodation
                     IF (dz_up .LE. dz_acc) THEN
                        ! ******************************************************* !
                        ! * CASE-1.2.1.1 :                                      * !
                        ! * ACCOMODATION POTENTIAL IS HIGHER THAN AUTHORIZED    * !
                        ! * UPWARD ACCOMODATION                                 * !
                        ! * THERE WILL HAVE BOTH UPWARD & DOWNWARD ACCOMODATION * !
                        ! ******************************************************* !
                        CaseS = '1.2.1.1'
                        zup_root(iv) = zup_root(iv) - dz_up
                        ! Root level moves upward (dz_up)
                        thick_root(iv) = MIN(thick_root(iv) + dz_acc, &
                                             var%c_opt_root_thick)
                        ! Update root level thickness
                     ELSE
                        ! **************************************************** !
                        ! * CASE-1.2.1.2 :                                   * !
                        ! * ACCOMODATION POTENTIAL IS LOWER THAN AUTHORIZED  * !
                        ! * UPWARD ACCOMODATION                              * !
                        ! * THERE WILL HAVE ONLY UPWARD ACCOMODATION         * !
                        ! **************************************************** !
                        CaseS = '1.2.1.2'
                        zup_root(iv) = zup_root(iv) - dz_acc
                        ! Root level moves upward (dz_up)
                        thick_root(iv) = MIN(thick_root(iv) + dz_acc, &
                                             var%c_opt_root_thick)
                        ! Update root level thickness
                     END IF ! * END TEST ON dz_up.LE.dz_acc
                  ELSE
                     ! ************************************************** !
                     ! * CASE-1.2.2 : ROOT LEVEL IS ABOVE OPTIMAL DEPTH * !
                     ! *              IT MAY ONLY ACCOMODATE DOWNWARD   * !
                     ! ************************************************** !
                     CaseS = '1.2.2'
                     thick_root(iv) = MIN(thick_root(iv) + dz_acc, &
                                          var%c_opt_root_thick)
                  END IF ! * END TEST ON hbs_zup_root.GT.hbs_c_opt_root_depth
               END IF ! * END TEST ON hbs_zup_root0.EQ.0.0
            ELSE
               ! *****************************
               ! * CASE-2 : DEPOSITION OCCURED
               ! *****************************
               ! Update first related to bed accretion
               zup_root(iv) = zup_root0(iv) + dhsed ! Roots level moves downward
               ! Now update due to accomodation
               IF (zup_root(iv) .LE. var%c_opt_root_depth) THEN
                  ! ********************************************************** !
                  ! * CASE-2.1 :                                             * !
                  ! * ROOT LEVEL IS ROOT LEVEL IS ABOVE OPTIMAL DEPTH        * !
                  ! * IT MAY ONLY ACCOMODATE DOWNWARD                        * !
                  ! ********************************************************** !
                  CaseS = '2.1'
                  thick_root(iv) = MIN(thick_root0(iv) + dz_acc, &
                                       var%c_opt_root_thick)
               ELSE
                  ! ********************************************************** !
                  ! * CASE-2.2 : ROOT LEVEL IS DEEPER THAN OPTIMAL DEPTH     * !
                  ! * IT MAY ACCOMODATE UPWARD (AND MAYBE DOWNWARD)          * !
                  ! ********************************************************** !
                  dz_up = -(var%c_opt_root_depth - zup_root(iv))
                  ! Authorized upward accomodation
                  IF (dz_up .LE. dz_acc) THEN
                     ! ******************************************************* !
                     ! * CASE-2.2.1 :                                        * !
                     ! * ACCOMODATION POTENTIAL IS HIGHER THAN AUTHORIZED    * !
                     ! * UPWARD ACCOMODATION                                 * !
                     ! * THERE WILL HAVE BOTH UPWARD & DOWNWARD ACCOMODATION * !
                     ! ******************************************************* !
                     CaseS = '2.2.1'
                     zup_root(iv) = zup_root(iv) - dz_up
                     ! Root level moves upward (dz_up)
                     thick_root(iv) = MIN(thick_root0(iv) + dz_acc, &
                                          var%c_opt_root_thick)
                     ! Update root level thickness
                  ELSE
                     ! ******************************************************* !
                     ! * CASE-2.2.2 :                                        * !
                     ! * ACCOMODATION POTENTIAL IS LOWER THAN AUTHORIZED     * !
                     ! * UPWARD ACCOMODATION                                 * !
                     ! * THERE WILL HAVE ONLY UPWARD ACCOMODATION            * !
                     ! ******************************************************* !
                     CaseS = '2.2.2'
                     zup_root(iv) = zup_root(iv) - dz_acc
                     ! Root level moves upward (dz_up)
                     thick_root(iv) = MIN(thick_root0(iv) + dz_acc, &
                                          var%c_opt_root_thick)
                     ! Update root level thickness
                  END IF ! * END TEST ON dz_up.LE.dz_acc
               END IF ! * END TEST ON hbs_zup_root.LE.hbs_c_opt_root_depth
            END IF ! * END TEST ON EROSION/DEPOSITION

            ! *********************** !
            ! * UPDATE TOTAL CHANGE * !
            ! *********************** !
            zup_root(iv) = MIN(MAX(zup_root(iv), 0.0_rsh), 1000.0_rsh)
            ! Limiter (0-1000 m) to debug
            thick_root(iv) = MIN(MAX(thick_root(iv), 0.0_rsh), 1000.0_rsh)
            ! Limiter (0-1000 m) to debug
            dz_root(iv) = zup_root0(iv) - zup_root(iv)
            ! Positive if root level moves upward (0 is interface with water)
            dthick_root(iv) = thick_root(iv) - thick_root0(iv)
            ! Positive if thickness increases
         END IF ! * END TEST ON l_hbs_bedsedstab AND position
      END DO ! * END LOOP ON iv

   END SUBROUTINE hbs1dv_update_root_level
   !============================================================================

   FUNCTION hbs1dv_accomodation_potential(var, position_wat, position_bed, &
                                          dt, jjulien) RESULT(dz_acc)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_accomodation_potential ***
   !!
   !! ** Purpose : COMPUTE ROOT ACCOMODATION POTENTIAL
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      TYPE(hbs1dv_param_type), INTENT(IN) :: var
      REAL(KIND=rsh), INTENT(IN) ::  position_wat
      REAL(KIND=rsh), INTENT(IN) ::  position_bed
      REAL(KIND=rsh), INTENT(IN) ::  dt
      INTEGER, INTENT(IN) :: jjulien
      REAL(KIND=rsh) :: dz_acc

      REAL(KIND=rsh) :: vel_tmp

      dz_acc = 0.0_rsh
      IF (position_wat > 0.0_rsh) THEN !! There are alived roots
         IF (var%l_root_accomodation) THEN
            IF (var%root_accomod_type .EQ. 0) THEN
               ! Constant accomodation velocity is used
               dz_acc = dt* &
                        hbs1dv_conv_mmmonth_2_msec(var%c_root_accomod_vel_max)
            ELSEIF (var%root_accomod_type .EQ. 1) THEN
               ! Parameterized accomodation velocity is used
               vel_tmp = hbs1dv_conv_mmmonth_2_msec(var%c_root_accomod_vel_max)
               vel_tmp = (vel_tmp - (0.05*vel_tmp))*0.5_rsh
               dz_acc = dt*(vel_tmp + vel_tmp*COS(2.0_rsh*pi* &
                                                  (REAL(jjulien, rsh) - &
                                                   var%c_root_accomod_day)/ &
                                                  365.0_rsh))
            END IF
         END IF
      END IF

   END FUNCTION hbs1dv_accomodation_potential
   !============================================================================

   FUNCTION hbs1dv_conv_mmmonth_2_msec(V) RESULT(Vconv)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_conv_mmmonth_2_msec ***
   !!
   !! ** Purpose : Convert velocity given in mm.month-1 into m.sec-1
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: V
      REAL(KIND=rsh) :: Vconv
      Vconv = (V/1000.0_rsh)/(30.0_rsh*86400.0_rsh)

   END FUNCTION hbs1dv_conv_mmmonth_2_msec
   !============================================================================

   FUNCTION hbs1dv_conv_mmday_2_msec(V) RESULT(Vconv)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_conv_mmday_2_msec ***
   !!
   !! ** Purpose : Convert velocity given in mm.day-1 into m.sec-1
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: V
      REAL(KIND=rsh) :: Vconv
      Vconv = (V/1000.0_rsh)/(86400.0_rsh)

   END FUNCTION hbs1dv_conv_mmday_2_msec
   !============================================================================

   FUNCTION hbs1dv_comp_eros_flx(tauskin, toce, excespowr, xeros) RESULT(ero)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_conv_mmmonth_2_msec ***
   !!
   !! ** Purpose : Computes raw erosion flux, in the same way than in
   !!              MUSTANG module
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: tauskin
      REAL(KIND=rsh), INTENT(IN) :: toce
      REAL(KIND=rsh), INTENT(IN) :: excespowr
      REAL(KIND=rsh), INTENT(IN) :: xeros
      REAL(KIND=rsh) :: ero
      IF (tauskin > toce) THEN
         ero = xeros*(tauskin/toce - 1.0_rsh)**excespowr
      ELSE
         ero = 0.0_rsh
      END IF

   END FUNCTION hbs1dv_comp_eros_flx
   !============================================================================

   FUNCTION hbs1dv_E0_coef(position_bed, zup_root, thick_root, root_biomass) &
      RESULT(E0_coef)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_E0_coef ***
   !!
   !! ** Purpose : Computes correction factor for erosion rate
   !!              due to presence of root system
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: position_bed
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: zup_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: thick_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: root_biomass
      REAL(KIND=rsh) :: E0_coef
      !! * Local declaration
      INTEGER        :: iv
      REAL(KIND=rsh) :: ntot, coef_tmp
      TYPE(hbs1dv_param_type) :: var

      !! * Executable part
      E0_coef = 0.0_rsh
      ntot = 0.0_rsh
      DO iv = 1, hbs1dv_nbvar
         var = hbs1dv_param(iv)
         coef_tmp = 0.0_rsh
         IF (var%l_bedsedstab) THEN
            IF (position_bed(iv) > 0.0_rsh) THEN
               IF ((zup_root(iv) == 0.0_rsh) .AND. &
                   (thick_root(iv) > 0.0_rsh)) THEN
                  coef_tmp = (position_bed(iv)* &
                              ((EXP((root_biomass(iv)/var%c_E0root_x0)* &
                                    var%c_E0root_x1) - 1.0_rsh) &
                               /(EXP(var%c_E0root_x1) - 1.0_rsh))) + &
                             position_bed(iv)
               END IF
            END IF
         END IF
         E0_coef = E0_coef + (1.0_rsh - coef_tmp)
         ntot = ntot + 1.0_rsh
      END DO
      E0_coef = MAX(E0_coef/ntot, 1.0E-6)
   END FUNCTION hbs1dv_E0_coef
   !============================================================================

   FUNCTION hbs1dv_tauce_coef(position_bed, zup_root, thick_root, &
                              root_biomass) RESULT(tauce_coef)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_tauce_coef ***
   !!
   !! ** Purpose : Computes correction factor for critical bottom shear
   !!              stress for erosion due to presence of root system
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: position_bed
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: zup_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: thick_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN) :: root_biomass
      REAL(KIND=rsh) :: tauce_coef
      !! * Local declaration
      INTEGER        :: iv
      REAL(KIND=rsh) :: ntot, coef_tmp
      TYPE(hbs1dv_param_type) :: var

      !! * Executable part
      tauce_coef = 0.0_rsh
      ntot = 0.0_rsh
      DO iv = 1, hbs1dv_nbvar
         var = hbs1dv_param(iv)
         coef_tmp = 0.0_rsh
         IF (var%l_bedsedstab) THEN
            IF (position_bed(iv) > 0.0_rsh) THEN
               IF ((zup_root(iv) == 0.0_rsh) .AND. &
                   (thick_root(iv) > 0.0_rsh)) THEN
                  coef_tmp = position_bed(iv)* &
                             (var%c_tauceroot_x0* &
                              root_biomass(iv)**var%c_tauceroot_x1)
               END IF
            END IF
         END IF
         tauce_coef = tauce_coef + (1.0_rsh + coef_tmp)
         ntot = ntot + 1.0_rsh
      END DO
      tauce_coef = tauce_coef/ntot
   END FUNCTION hbs1dv_tauce_coef
   !============================================================================

   LOGICAL FUNCTION hbs1dv_isthere_roots(position_bed, thick_root)
   !!---------------------------------------------------------------------------
   !!                 *** FUNCTION hbs1dv_isthere_roots ***
   !!
   !! ** Purpose : Check if there is possibility to have roots in current cell
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: position_bed
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: thick_root
      !! * Local declaration
      INTEGER   :: iv
      !! * Executable part
      hbs1dv_isthere_roots = .FALSE.
      DO iv = 1, hbs1dv_nbvar
         IF ((hbs1dv_param(iv)%l_bedsedstab) .AND. &
             (position_bed(iv) .GT. 0.0_rsh) .AND. &
             (thick_root(iv) .GT. 0.0_rsh)) THEN
            hbs1dv_isthere_roots = .TRUE.
         END IF
      END DO
   END FUNCTION hbs1dv_isthere_roots
   !============================================================================

   SUBROUTINE hbs1dv_comp_erosion(position_bed, thick_root, &
                                  zup_root, root_biomass, &
                                  tauskin, toce, csanmud, excespowr, xeros, &
                                  dzs, MF, dt, erosi)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_comp_erosion ***
   !!
   !! ** Purpose : Computes total sediment mass that can be eroded during
   !!              a time-step (dt)
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: position_bed
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: thick_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: zup_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: root_biomass
      REAL(KIND=rsh), INTENT(IN) :: tauskin, csanmud, excespowr, MF, dt, dzs
      REAL(KIND=rsh), INTENT(IN) :: toce, xeros
      REAL(KIND=rsh), INTENT(INOUT) :: erosi
      !! * Local declaration
      LOGICAL        :: Lchange
      INTEGER        :: iv
      REAL(KIND=rsh) :: dt_rest, toce_tmp, xeros_tmp, erosi_tmp, &
                        ero_tmp, dz_erod, dz_erod_max, Zchange, dt_used
      REAL(KIND=rsh), DIMENSION(1:hbs1dv_nbvar) :: zroot_tmp, troot_tmp
      !! * Executable part
      IF (csanmud > 0.0_rsh) THEN
         IF (hbs1dv_isthere_roots(position_bed, thick_root)) THEN
          !! *********************************
          !! * THERE ARE ROOTS IN CURRENT CELL
          !! *********************************
            !----------------------------
            ! Preliminary initializations
            !----------------------------
            erosi = 0.0_rsh
            dt_rest = dt
            zroot_tmp(:) = zup_root(:)
            troot_tmp(:) = thick_root(:)
            ! *****************
            ! * START PROCEDURE
            ! *****************
            DO WHILE (dt_rest > 0.0_rsh)
               !--------------------------------------------
               ! Computation of corrected erosion parameters
               !--------------------------------------------
               toce_tmp = toce*hbs1dv_tauce_coef(position_bed, zroot_tmp, &
                                                 troot_tmp, root_biomass)
               xeros_tmp = xeros*hbs1dv_E0_coef(position_bed, zroot_tmp, &
                                                troot_tmp, root_biomass)
               !----------------------
               ! Erosion mass [kg.m-2]
               !----------------------
               ero_tmp = hbs1dv_comp_eros_flx(tauskin, toce_tmp, excespowr, &
                                              xeros_tmp)
               erosi_tmp = ero_tmp*MF*dt_rest
               !----------------------------------------------------------------
               ! Maximum possible erosion with in these sediment/root conditions
               !----------------------------------------------------------------
               dz_erod_max = erosi_tmp/csanmud
               !--------------------------------------------
               ! Search if there may have a change of regime
               !--------------------------------------------
               Lchange = .FALSE.
               Zchange = 999.0_rsh
               CALL hbs1dv_search_next_change(zroot_tmp, troot_tmp, &
                                              Lchange, Zchange)
               IF ((.NOT.Lchange) .OR. (Zchange .GE. dz_erod_max)) THEN
                  !----------------------------------------------------------
                  ! No more change in bed characteristics during the
                  ! time-step or the current sediment layer
                  ! OR 
                  ! erosion does not reach the change 
                  !----------------------------------------------------------
                  erosi = erosi + erosi_tmp
                  dt_rest = 0.0_rsh
               ELSE
                  !----------------------------------------------------------
                  ! There is a change in bed characteristics during the
                  ! time-step or the current sediment layer
                  !----------------------------------------------------------
                  erosi_tmp = Zchange*csanmud  ! Mass to be eroded
                  erosi = erosi + erosi_tmp
                  dt_used = ((Zchange*csanmud)/(ero_tmp*MF)) ! Time need to erod
                  ! from the top to the depth where change occured
                  dt_rest = dt_rest - dt_used
                  ! Update Zroot and Troot
                  CALL hbs1dv_update_zroot_troot(Zchange, troot_tmp, zroot_tmp)
               END IF
            END DO ! do while dt_rest > 0.
         ELSE ! hbs1dv_isthere_roots
          !! ************************************
          !! * THERE ARE NO ROOTS IN CURRENT CELL
          !! ************************************
            erosi = hbs1dv_comp_eros_flx(tauskin, toce, excespowr, xeros)*MF*dt
         END IF ! hbs1dv_isthere_roots
      ELSE ! csanmud > 0
         erosi = 0.0_rsh
      END IF
   END SUBROUTINE hbs1dv_comp_erosion
   !============================================================================

   SUBROUTINE hbs1dv_search_next_change(Zroot, Troot, Lchange, Zchange)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_search_next_change ***
   !!
   !! ** Purpose : Search next change in bottom conditions
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(KIND=rsh), DIMENSION(1:hbs1dv_nbvar), INTENT(IN) :: Zroot, Troot
      LOGICAL, INTENT(OUT)        :: Lchange
      REAL(KIND=rsh), INTENT(OUT) :: Zchange
      !! * Local declaration
      INTEGER        :: i, j, iv
      REAL(KIND=rsh) :: Zr, Tr
      LOGICAL, DIMENSION(1:hbs1dv_nbvar)        :: Lch_iv
      REAL(KIND=rsh), DIMENSION(1:hbs1dv_nbvar) :: Zch_iv
      !! * Executable part
      Lch_iv(:) = .FALSE.
      Zch_iv(:) = 999.0_rsh
      DO iv = 1, hbs1dv_nbvar
         IF ((Zroot(iv) .GT. 0.0_rsh) .AND. (Troot(iv) .GT. 0.0_rsh)) THEN
            !--------------------------------------
            ! Roots exist and are below the surface
            !--------------------------------------
            Lch_iv(iv) = .TRUE.
            Zch_iv(iv) = Zroot(iv)
         ELSEIF ((Zroot(iv) .EQ. 0.0_rsh) .AND. (Troot(iv) .GT. 0.0_rsh)) THEN
            !---------------------------------------------------
            ! Root exist and are at the sediment-water interface
            !---------------------------------------------------
            Lch_iv(iv) = .TRUE.
            Zch_iv(iv) = Troot(iv)
         END IF
      END DO
      Zchange = MINVAL(Zch_iv(:))
      IF (Zchange .EQ. 999.0_rsh) THEN
         Lchange = .FALSE.
      ELSE
         Lchange = .TRUE.
      END IF
   END SUBROUTINE hbs1dv_search_next_change
   !============================================================================

   SUBROUTINE hbs1dv_zroot_troot_ero(position_bed, dz_root, &
                                     thick_root, zup_root)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_zroot_troot_ero ***
   !!
   !! ** Purpose : Update root level depth and thickness when erosion occured
   !!              during erosion procedure
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(IN)  :: position_bed
      REAL(KIND=rsh), INTENT(IN) :: dz_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: thick_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: zup_root
      !! * Local declaration
      INTEGER        :: iv
      LOGICAL        :: Lroots
      REAL(KIND=rsh) :: dz_rest
      !! * Executable part
      Lroots = hbs1dv_isthere_roots(position_bed, thick_root)
      IF ((Lroots) .AND. (dz_root .GT. 0.0_rsh)) THEN
         CALL hbs1dv_update_zroot_troot(dz_root, thick_root, zup_root)
      END IF

   END SUBROUTINE hbs1dv_zroot_troot_ero
   !============================================================================

   SUBROUTINE hbs1dv_update_zroot_troot(dz_root, thick_root, zup_root)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs1dv_zroot_troot_ero ***
   !!
   !! ** Purpose : Update root level depth and thickness when erosion occured
   !!              during erosion procedure
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(KIND=rsh), INTENT(IN) :: dz_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: thick_root
      REAL(rsh), DIMENSION(hbs1dv_nbvar), INTENT(INOUT)  :: zup_root
      !! * Local declaration
      INTEGER        :: iv
      REAL(KIND=rsh) :: dz_rest
      !! * Executable part
      DO iv = 1, hbs1dv_nbvar
         IF (zup_root(iv) .EQ. 0.0_rsh) THEN
            !-------------------------------------------------------------
            ! Roots are already at sediment surface, no change on zup_root
            !-------------------------------------------------------------
            IF (hbs1dv_param(iv)%l_root_erosion) THEN
               thick_root(iv) = MAX(thick_root(iv) - dz_root, 0.0_rsh)
            END IF
         ELSE
            IF (zup_root(iv) .GE. dz_root) THEN
               !------------------------------------------------------
               ! Roots are not at sediment surface, no change on troot
               !------------------------------------------------------
               zup_root(iv) = zup_root(iv) - dz_root
            ELSE
               dz_rest = dz_root - zup_root(iv)
               zup_root(iv) = 0.0_rsh
               IF (hbs1dv_param(iv)%l_root_erosion) THEN
                  thick_root(iv) = MAX(thick_root(iv) - dz_rest, 0.0_rsh)
               END IF
            END IF
         END IF
      END DO

   END SUBROUTINE hbs1dv_update_zroot_troot
   !============================================================================

END MODULE HYBIOSED1DV
