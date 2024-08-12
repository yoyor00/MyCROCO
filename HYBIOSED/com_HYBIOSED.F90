MODULE com_HYBIOSED

#include "cppdefs.h"

#if defined HYBIOSED
   !!===========================================================================
   !!                   ***  MODULE com_HYBIOSED  ***
   !!
   !! ** Purpose : Declaration of all variables used by HYBIOSED module
   !!
   !!===========================================================================

   IMPLICIT NONE

   ! default
   PUBLIC

   ! ---------------------------------------------------------------------------
   ! Definition of rsh, rlg, riosh, riolg, lchain
   ! ---------------------------------------------------------------------------
   INTEGER, PARAMETER :: riosh = 8, riolg = 8, rlg = 8, rsh = 8 ! working precision
   INTEGER, PARAMETER :: lchain = 200

   INTEGER :: hbs_kmax ! number of layers in water
   INTEGER :: ierrorlog, iwarnlog, iscreenlog ! logging IO

   ! Variables nomenclature :
   ! hbs_*    : variables header
   ! hbs_i_*  : variables used for initialization
   ! hbs_c_*  : constant parameters
   ! hbs_fn_* : filename variables
   ! hbs_l_* or l_hbs* : logical variables

   ! * Variables for inputs purpose
   !-------------------------------
   CHARACTER(LEN=lchain) :: hbs_fn_position
   INTEGER :: hbs_nbvar

   ! * Interaction suspended sediment/canopy
   !----------------------------------------
   LOGICAL :: l_hbs_suspsed_trapp              ! Activation of effects of obstructions to flow on suspended particles trapping (FORCE Substance/MUSTANG module)
   LOGICAL :: l_hbs_suspsed_block              ! Activation of effects of obstructions to flow on suspended particles blockage (FORCE Substance/MUSTANG module)

   REAL(KIND=rsh) :: hbs_c_suspsed_trapp_max   ! Maximum correction factor on settling velocity due to particles trapping
   REAL(KIND=rsh) :: hbs_c_suspsed_trapp_exp   ! Exponential coefficient for correction factor on settling velocity due to particles trapping [-Inf[0]+Inf]
   REAL(KIND=rsh) :: hbs_c_suspsed_block_max   ! Maximum correction factor on settling velocity due to particles blockage [0:1]
   REAL(KIND=rsh) :: hbs_c_suspsed_block_exp   ! Exponential coefficient for correction factor on settling velocity due to particles blockage [-Inf[0]+Inf]

   ! * Variables on (k,i,j)
   !-----------------------
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_ws_trapp    ! The coefficient for increased settling velocity due to frontal collision (k,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_ws_block    ! The coefficient for decreased settling velocity due to vertical blackage (k,i,j), [-]

   ! * Variables on (iv,i,j)
   !-----------------------
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_zup_root0 ! The depth of root level (from seabed) at the beginning of time-step (iv,i,j,) [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_thick_root0 ! The thickness of root level at the beginning of time-step (iv,i,j) [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_zup_root ! The depth of root level (from seabed) for change in erosion parameters (iv,i,j,) [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_thick_root ! The thickness of root level for change in erosion parameters (iv,i,j) [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_dz_root ! The change of root level depth for change in erosion parameters (iv,i,j,) [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_dthick_root ! The change of root level thickness (iv,i,j) [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_position_wat ! The table where the position of the differents wat variables is defined through occupation rate (iv,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_position_bed ! The table where the position of the differents bed variables is defined through occupation rate (iv,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: hbs_root_biomass ! The biomass of roots (iv,i,j) [gDW.m-2]

   ! * Variables on (i,j)
   !-----------------------
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: hbs_hsed_prev

   ! * Variables on (iv)
   CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: hbs_fn_var             ! Name of parameters file
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_tauceroot_x0              ! First coefficient for tauce coefficient due to root presence
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_tauceroot_x1              ! Second coefficient for tauce coefficient due to root presence
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_E0root_x0                 ! First coefficient for E0 coefficient due to root presence
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_E0root_x1                 ! Second coefficient for E0 coefficient due to root presence
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_root_accomod_vel_max      ! Maximum vertical accomodation of root system (mm.month-1)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_root_accomod_day          ! Day of the year for maximum accomodation velocity
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_opt_root_depth            ! Optimal root level depth (m)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_c_opt_root_thick            ! Optimal (maximum) root level thickness (m)
   INTEGER, DIMENSION(:), ALLOCATABLE :: hbs_c_root_accomod_type                ! Type of root accomodation (0:constant, 1:parameterized depending on day of the year, 2:from biological model)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: l_hbs_root_accomodation                ! Take account (or not) for vertical root accomodation
   LOGICAL, DIMENSION(:), ALLOCATABLE :: l_hbs_root_erosion                     ! Take account (or not) for root erosion
   CHARACTER(len=lchain), DIMENSION(:), ALLOCATABLE :: hbs_varname              ! Name of variable (iv)
   CHARACTER(len=lchain), DIMENSION(:), ALLOCATABLE :: hbs_vardesc              ! Short description
   LOGICAL, DIMENSION(:), ALLOCATABLE :: l_hbs_bedsedstab                       ! If this variable as an impact on bed sediment stabilization (USE Mustang module)

   LOGICAL, DIMENSION(:), ALLOCATABLE :: l_hbs_initfromfile                     ! If initialization from a restart file
   CHARACTER(len=lchain), DIMENSION(:), ALLOCATABLE :: hbs_fn_initfile          ! Name of file for root biomass spatial initialization
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_i_rbiom                     ! Initial root biomass (if l_hbs_rbiom_hom=.TRUE. AND l_hbs_rbiom_stat=.TRUE.)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_i_zroot                     ! Initial root level depth (if l_hbs_zroot_hom=.TRUE.)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: hbs_i_thickroot                 ! Initial root level thickness (if l_hbs_thickroot_hom=.TRUE.)

   LOGICAL :: l_hbsout_pos
   LOGICAL :: l_hbsout_susp_trapp
   LOGICAL :: l_hbsout_susp_block
   LOGICAL :: l_hbsout_rbiom
   LOGICAL :: l_hbsout_zroot
   LOGICAL :: l_hbsout_thickroot

   INTEGER, DIMENSION(:), ALLOCATABLE :: hisHbs                    ! Output identifier
   INTEGER, DIMENSION(:), ALLOCATABLE :: avgHbs                    ! Output identifier
   LOGICAL, DIMENSION(:), ALLOCATABLE :: outHbs                    ! To choose which variable is outputed
   LOGICAL, DIMENSION(:), ALLOCATABLE :: out2DHbs                  ! To indicate if the output variable is 2D or 3D
   CHARACTER(LEN=75), DIMENSION(:, :), ALLOCATABLE :: vname_hbs



   CHARACTER(LEN=lchain) :: hbs_nout_pos_bed  ! Name of in output file
   CHARACTER(LEN=lchain) :: hbs_nout_pos_wat ! Name of in output file
   CHARACTER(LEN=lchain) :: hbs_nout_susp_trapp  ! Name variable correction factor on settling velocity due to trapping
   CHARACTER(LEN=lchain) :: hbs_nout_susp_block  ! Name variable correction factor on settling velocitydue to blockage
   CHARACTER(LEN=lchain) :: hbs_nout_rbiom ! Name variable root biomass
   CHARACTER(LEN=lchain) :: hbs_nout_zroot ! Name variable depth of root level
   CHARACTER(LEN=lchain) :: hbs_nout_thickroot ! Name variable thickness of root level
   CHARACTER(LEN=lchain) :: hbs_nout_dzroot ! Name variable change of depth of root level
   CHARACTER(LEN=lchain) :: hbs_nout_dthickroot ! Name variable change of thickness of root level
   CHARACTER(LEN=lchain) :: hbs_nout_tauce_coef  ! Name variable correction factor for critical shears stress for erosion
   CHARACTER(LEN=lchain) :: hbs_nout_E0_coef ! Name variable correction factor for erosion rate

CONTAINS

   !!===========================================================================

#endif /* HYBIOSED */
END MODULE com_HYBIOSED
