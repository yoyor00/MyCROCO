 MODULE comOBSTRUCTIONS

#ifdef key_OBSTRUCTIONS
   !&E======================================================================
   !&E                   ***  MODULE comOBSTRUCTIONS  ***
   !&E
   !&E ** Purpose : Declaration of all variables used by obstructions module
   !&E
   !&E ** Description :
   !&E
   !&E     subroutine OBSTRUCTIONS_alloc_nbvar        ! Allocates variables depending on
   !&E                                                ! number of obstructions variables
   !&E
   !&E     subroutine OBSTRUCTIONS_alloc_xyz          ! Allocates spatial tables for obstructions
   !&E
   !&E     subroutine OBSTRUCTIONS_alloc_other        ! Allocates other tables for obstructions
   !&E
   !&E     subroutine OBSTRUCTIONS_dealloc            ! Deallocates variables for obstructions
   !&E
   !&E ** History :
   !&E     ! 2018-04-13 (F. Ganthy) Original code : moved from obstructions.F90
   !&E     ! 2022-01-11 (A. Le Pevedic) Introduction of averaged parameters to send to WW3 
   !&E
   !&E======================================================================

#include "toolcpp.h"

   !! * Module used
   USE parameters, ONLY : lchain,rsh,riosh

   !! * Declaration
   IMPLICIT NONE

   !! * Shared or public module variables

   ! Variables nomenclature :
   ! obst_*    : variables for obstructions
   ! obst_i_*  : variables used for initialization
   ! obst_c_*  : constant parameters
   ! obst_fn_* : filename variables
   ! l_obst_*  : logical variables

   !--------------------------------------------------------------------------
   ! * VARIABLES WITHOUT DEPENDANCE ON THE NUMBER OF OBSTRUCTIONS VARIABLES  :
   !--------------------------------------------------------------------------

   ! * Variables for inputs purpose
   CHARACTER(LEN=lchain),PUBLIC :: obst_fn_position                         ! Name of the input file for the obstruction position within the domain

   ! * Variables for outputs purpose
   CHARACTER(LEN=lchain),PUBLIC :: obst_fn_out                              ! Name of the output file for obstructions

   LOGICAL,PUBLIC :: l_obstout_pos                                          ! Write obstruction position (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_height_f                                     ! Write 2D obstruction height (forcing) (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_height_e                                     ! Write 2D obstruction height (effective) (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_height_mean                                  ! Write mean obstruction height (over all obstructions) (i,j)
   LOGICAL,PUBLIC :: l_obstout_dens_f                                       ! Write 2D obstruction density (forcing) (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_dens_e                                       ! Write 3D obstruction density (3D effective) (iv,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_dens_2d                                      ! Write 2D obstruction density (2D effective) (i,j)
   LOGICAL,PUBLIC :: l_obstout_width_f                                      ! Write 2D obstruction width (forcing) (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_width_e                                      ! Write 3D obstruction width (3D effective) (iv,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_width_2d                                     ! Write 2D obstruction width (2D effective) (i,j)
   LOGICAL,PUBLIC :: l_obstout_thick_f                                      ! Write 2D obstruction thick (forcing )(iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_thick_e                                      ! Write 3D obstruction thick (3D effective) (iv,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_oai                                          ! Write 2D obstruction area index, (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_theta                                        ! Write 3D obstruction bending angle, (iv,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_cover                                        ! Write 2D obstruction coverage, (iv,i,j)
   LOGICAL,PUBLIC :: l_obstout_frac_z                                       ! Write 3D obstruction fraction of sigma layer, (iv,k,i,j)

   LOGICAL,PUBLIC :: l_obstout_fuv                                          ! Write 2D obstruction friction force (i,j)
   LOGICAL,PUBLIC :: l_obstout_fuzvz                                        ! Write 3D obstruction friction force (k,i,j)
   LOGICAL,PUBLIC :: l_obstout_a2d                                          ! Write 2D obstruction horizontal area (iv+3,i,j)
   LOGICAL,PUBLIC :: l_obstout_a3d                                          ! Write 3D obstruction horizontal area (iv+3,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_s2d                                          ! Write 2D obstruction vertical area (iv+3,i,j)
   LOGICAL,PUBLIC :: l_obstout_s3d                                          ! Write 3D obstruction vertical area (iv+3,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_drag                                         ! Write 3D obstruction drag coefficient (iv,k,i,j)
   LOGICAL,PUBLIC :: l_obstout_tau                                          ! Write 3D obstruction turbulence dissipation scale (k,i,j)
   LOGICAL,PUBLIC :: l_obstout_z0bed                                        ! Write 2D bottom roughness length for bed (i,j)
   LOGICAL,PUBLIC :: l_obstout_z0obst                                       ! Write 2D bottom roughness length for obstructions (iv+3,i,j)
   LOGICAL,PUBLIC :: l_obstout_z0bstress                                    ! Write 2D bottom roughness length used for bottom shear stress computation (i,j)
   LOGICAL,PUBLIC :: l_obstout_bstress                                      ! Write 2D total bottom shear stress (i,j)
   LOGICAL,PUBLIC :: l_obstout_bstressc                                     ! Write 2D current bottom shear stress (i,j)
   LOGICAL,PUBLIC :: l_obstout_bstressw                                     ! Write 2D wave bottom shear stress (i,j)

   CHARACTER(LEN=lchain),PUBLIC :: name_out_pos                             ! Name obstruction position (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_height_f                        ! Name 2D obstruction height (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_height_e                        ! Name 2D obstruction height (effective) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_dens_f                          ! Name 2D obstruction density (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_dens_e                          ! Name 3D obstruction density (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_width_f                         ! Name 2D obstruction width (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_width_e                         ! Name 3D obstruction width (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_thick_f                         ! Name 2D obstruction thick (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_thick_e                         ! Name 3D obstruction thick (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_oai                             ! Name 2D obstruction area index, (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_theta                           ! Name 3D obstruction bending angle (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_cover                           ! Name 2D obstruction coverage, (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_frac_z                          ! Name 3D obstruction fraction of sigma layer, (iv,i,j)


   CHARACTER(LEN=lchain),PUBLIC :: name_out_fuv                             ! Name 2D obstruction friction force (i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_fuzvz                           ! Name 3D obstruction friction force (k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_a2d                             ! Name 2D obstruction horizontal area (iv+3,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_a3d                             ! Name 3D obstruction horizontal area (iv+3,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_s2d                             ! Name 2D obstruction vertical area (iv+3,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_s3d                             ! Name 3D obstruction vertical area (iv+3,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_drag                            ! Name 3D obstruction drag coefficient (k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_tau                             ! Name 3D obstruction turbulence dissipation scale (k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_z0bed                           ! Name 2D bottom roughness length for bed (i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_z0obst                          ! Name 2D bottom roughness length for obstruction (iv+3,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_z0bstress                       ! Name 2D bottom roughness length used for bottom shear stress computation (i,j)
   CHARACTER(LEN=lchain),PUBLIC :: name_out_bstress                         ! Name 2D total bottom shear stress within output file
   CHARACTER(LEN=lchain),PUBLIC :: name_out_bstressc                        ! Name 2D current bottom shear stress within output file
   CHARACTER(LEN=lchain),PUBLIC :: name_out_bstressw                        ! Name 2D wave bottom shear stress within output file

   REAL(KIND=riosh),PUBLIC :: riog_valid_min_pos,riog_valid_max_pos         ! Valid minimum and maximum for output obstruction position
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_height,riog_valid_max_height   ! Valid minimum and maximum for output obstruction height
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_dens,riog_valid_max_dens       ! Valid minimum and maximum for output obstruction density
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_width,riog_valid_max_width     ! Valid minimum and maximum for output obstruction width
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_thick,riog_valid_max_thick     ! Valid minimum and maximum for output obstruction tickness
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_oai,riog_valid_max_oai         ! Valid minimum and maximum for output obstruction area index
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_theta,riog_valid_max_theta     ! Valid minimum and maximum for output obstruction bending angle
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_cover,riog_valid_max_cover     ! Valid minimum and maximum for output obstruction coverage
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_fracz,riog_valid_max_fracz     ! Valid minimum and maximum for output obstruction fraction in z
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_fuv,riog_valid_max_fuv         ! Valid minimum and maximum for output obstructin friction force
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_a,riog_valid_max_a             ! Valid minimum and maximum for output obstruction horizontal area
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_s,riog_valid_max_s             ! Valid minimum and maximum for output obstruction vertical area
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_drag,riog_valid_max_drag       ! Valid minimum and maximum for output obstruction drag
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_tau,riog_valid_max_tau         ! Valid minimum and maximum for output obstruction tau
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_z0,riog_valid_max_z0           ! Valid minimum and maximum for output obstruction roughness length
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_bstress,riog_valid_max_bstress ! Valid minimum and maximum for output bottom shear stress
   REAL(KIND=riosh),PUBLIC :: riog_valid_min_zroot,riog_valid_max_zroot     ! Valid minimum and maximum for output obstruction root depth

   ! Other variables/parameters
   INTEGER,PUBLIC  :: obst_iv
   INTEGER,PUBLIC  :: obst_nbvar                       ! The total number of obstruction variables
   INTEGER,PUBLIC  :: obst_nv_up                       ! Number of variable for upward obstructions
   INTEGER,PUBLIC  :: obst_nv_do                       ! Number of variable for downward obstruction
   INTEGER,PUBLIC  :: obst_nv_3d                       ! Number of variable for full 3d obstructions
   INTEGER,PUBLIC  :: obst_nv_turb                     ! Number of variable for obstructions using full turbulence procedure
   INTEGER,PUBLIC  :: obst_nv_noturb                   ! Number of variable for obstructions using simplified (roughness length) procedure
   INTEGER,PUBLIC  :: obst_nv_rigid_up                 ! Number of variable for rigid and upward (from sea-bed) obstructions
   INTEGER,PUBLIC  :: obst_nv_rigid_do                 ! Number of variable for rigid and downward (from sea-surface) obstructions
   INTEGER,PUBLIC  :: obst_nv_flexi_up                 ! Number of variable for flexible and upward (from sea-bed) obstructions
   INTEGER,PUBLIC  :: obst_nv_flexi_do                 ! Number of variable for flexible and downward (from sea-surface) obstructions
   INTEGER,PUBLIC  :: obst_kmax                        ! Number of sigma layer effectively used (=kmax for 3D modele, =obst_kmax2d for 2dmodele)

   LOGICAL,PUBLIC :: l_obst_z0bstress_tot              ! IF ONLY ONE obstruction VARIABLE USED Z0SED

   REAL(KIND=rsh),PUBLIC  :: obst_c_paramhuv           ! The coefficient of obstruction height for computation of velocity
   REAL(KIND=rsh),PUBLIC  :: obst_c_imp3d              ! Implicitation coefficient for obstructions formulations in 3D
   REAL(KIND=rsh),PUBLIC  :: obst_c_exp3d              ! Explicitation coefficient for obstructions formulations in 3D
   REAL(KIND=rsh),PUBLIC  :: obst_c_imp2d              ! Implicitation coefficient for obstructions formulations in 2D
   REAL(KIND=rsh),PUBLIC  :: obst_c_exp2d              ! Explicitation coefficient for obstructions formulations in 2D

!FG settling velocity term for zostera (should be defined as iv dependant)
#if defined key_casobstflume_ganthy2015_sedim
   INTEGER,PUBLIC         :: obst_c_flumetrapp_opt   ! Option for flume trapping factor (0 : linear, 1 : power, 2 : exp, 3 : linear LOG10)
   REAL(KIND=rsh),PUBLIC  :: obst_c_flumetrapp_x0    ! First term for flume trapping factor
   REAL(KIND=rsh),PUBLIC  :: obst_c_flumetrapp_x1    ! Second term for flume trapping factor
#endif
   REAL(KIND=rsh),PUBLIC  :: obst_c_settletrapp_x0   ! first term for parameterized obstructions trapping factor
   REAL(KIND=rsh),PUBLIC  :: obst_c_settletrapp_x1   ! second term for parameterized obstructions trapping factor 
   REAL(KIND=rsh),PUBLIC  :: obst_c_settleblock_x0   ! first term for parameterized obstructions blocking factor
   REAL(KIND=rsh),PUBLIC  :: obst_c_settleblock_x1   ! second term for parameterized obstructions blocking factor
!END FG

   REAL(KIND=rsh),PUBLIC  :: obst_i_z0bstress          ! Roughness length for bottom shear stress without obstruction (= z0seduni if key_sedim)
   REAL(KIND=rsh),PUBLIC  :: fricwav,fws2              ! Wave related friction factor for bottom shear stress 

   INTEGER,PARAMETER,PUBLIC         :: obst_kmax2d = 20       ! Number of sigma layer used for modele2D
   REAL(KIND=rsh),PARAMETER,PUBLIC  :: obst_p_hmin = 0.05_rsh ! Minimum coefficient value for obstruction height under bending  

   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ONLY ON THE NUMBER OF OBSTRUCTIONS VARIABLES (iv) :
   !--------------------------------------------------------------------------

   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_varnum                   ! The number allocated to each variables (iv)
   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_fracxy_type              ! The type of correction for horizontal coverage ((grid cell not completely fill with obstructions) (iv)
   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_nbhnorm                  ! Number of vertical steps from the distribution file (iv)
   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_abdel_nmax             ! Number of segments for Abdlerhman method (bending) (iv) 

   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_filechar               ! For reading time-series of obstructions characteristics from a file (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_filedistri             ! For reading the vertical distribution of obstructions from a file (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_init_spatial           ! For reading spatial file for density, height, width and thickness (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_flexible               ! For obstructions flexibility (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_cylindre               ! For obstruction shape (cylinder, ellispe / parallelepipeds), (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_downward               ! For downward obstructions (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_3dobst                 ! For fully 3D obstructions (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_noturb                 ! For use of simplified formulation (roughness length) instead of full turbulent formulation (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_abdelrough_cste        ! For use constant Abdelrhman 2003 coefficient (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_fracxy                 ! For horizontal coverage correction (grid cell not completely fill with obstructions) (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_abdelposture           ! For computation of obstructions posture following abdelrhman 2007 (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_param_height           ! For use a parameterization (from velocity) for obstruction height (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_drag_cste              ! For constant/variable drag coefficient (or corrected from bending angle) (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: l_obst_z0bstress              ! For using a roughness length for obstructions (iv)

   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_z0bstress_option         ! For using various parameterizations of bottom roughness


   CHARACTER(len=lchain),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_varname    ! Name of obstructions variables (iv)
   CHARACTER(len=lchain),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_fn_vardat  ! Name of the input file for obstruction variables
   CHARACTER(len=lchain),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_fn_char    ! Name of the file for times-series of obstructions characteristics (iv)
   CHARACTER(len=lchain),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_fn_distrib ! Name of the file for the vertical distribution of obstructions (iv)

   ! Initialization variables
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_i_height          ! Initial height of obstructions (iv), [m]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_i_width           ! Initial width of obstructions element (iv), [m]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_i_thick           ! Initial thickness of obstructions element (iv), [m]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_i_dens            ! Initial density of obstructions (iv), [m-2]

   ! Constant parameters
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_rho             ! Initial density of obstructions (iv), [kg.m-3]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_drag            ! Drag coefficient for obstructions (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_lift            ! Lift coefficient for obstructions (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_fracxy_k0       ! First parameter for horizontal coverage correction (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_fracxy_k1       ! Second parameter for horizontal coverage correction (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_fracxy_l        ! Third parameter for horizontal coverage correction (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_lz              ! obstructions spacing coefficient (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_shelter         ! sheltering coefficient (iv) [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_crough_x0       ! First coefficient (or cst value) for Abdelrhman (2003) c value [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_crough_x1       ! Second coefficient for Abdelrhman (2003) c value [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_height_x0       ! First parameter for obstruction or Zostera noltii height (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_height_x1       ! Second parameter for obstruction or Zostera noltii height (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_z0bstress       ! roughness length for obstructions (iv), [m]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_z0bstress_x0    ! First parameter for roughness length parameterization (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_z0bstress_x1    ! Second parameter for roughness length parameterization (iv), [-]
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_z0bstress_x2    ! Third parameter for roughness length parameterization (correction for 2Dsmall-depth) (iv), [-]


   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ONLY ON GRID (i,j) or (k,i,j)                     :
   !--------------------------------------------------------------------------
   ! Variables on (k)
   !-----------------
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC   :: obst_sig             ! Sigma levels (same as sig (comvars3d) for 3D modele, or allocated on obst_kmax2d for 2D modele)
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC   :: obst_dsig            ! Sigme layer thickness (same as sig (comvars3d) for 3D modele, or allocated on obst_kmax2d for 2D modele)
   ! Variables on (i,j)
   !----------------------
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_roswat_bot      ! Water density at the bottom (i,j)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_fu_i            ! Implicit Sink term in the momentum equation for 2D formulation (x direction) (i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_fv_i            ! Implicit Sink term in the momentum equation for 2D formulation (y direction) (i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_fu_e            ! Explicit Sink term in the momentum equation for 2D formulation (x direction) (i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_fv_e            ! Explicit Sink term in the momentum equation for 2D formulation (y direction) (i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_z0bed           ! z0 for bed without obstructions from source code (used where no obstructions are present) (i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_bstress         ! Total bottom shear stress (i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_bstressc        ! Current bottom shear stress (i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_bstressw        ! Wave bottom shear stress (i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_z0bstress       ! z0Sed from sedimento (i,j) but modified by obstructions, [m]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_raphbx
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_raphby
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_frofonx
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_frofony
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_dens_mean       ! Mean obstruction density (i,j), [n.m-2]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_width_mean      ! Mean obstruction width (i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_height_mean     ! Mean obstruction height (i,j) [m]


   ! Variables on (k,i,j)
   !---------------------
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_zc            ! Height z at center of sigma layer (k,i,j), (m)
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_dz            ! Dz (k,i,j), (m)
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_uz            ! 3D (or pseudo-3D) velocities used within the whole obstructions procedure (x direction) (k,i,j), (m.s-1)
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_vz            ! 3D (or pseudo-3D) velocities used within the whole obstructions procedure (y direction) (k,i,j), (m.s-1)
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_fuz_i         ! Implicit Sink term in the momentum equation for 3D formulation (x direction) (k,i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_fvz_i         ! Implicit Sink term in the momentum equation for 3D formulation (y direction) (k,i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_fuz_e         ! Explicit Sink term in the momentum equation for 3D formulation (x direction) (k,i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_fvz_e         ! Explicit Sink term in the momentum equation for 3D formulation (y direction) (k,i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_t             ! Source term T(z) in y equation (k,i,j), [N.m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_tau           ! Source term tau_veg in y equation (k,i,j), [N.m-2]

   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ONLY ON BOTH (iv,i,j) or (iv,k,i,j) or other      :
   !--------------------------------------------------------------------------
   ! Variables on (iv,i,j)
   !----------------------
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_dens_inst       ! Instantaneous obstruction real density (iv,i,j), [m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_width_inst      ! Instantaneous obstruction real width (iv,i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_thick_inst      ! Instantaneous obstruction real thickness (iv,i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_height_inst     ! Instantaneous obstruction real unbend height (iv,i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_area_index_inst ! Instantaneous obstruction real area index (iv,i,j), [-]

   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_position      ! The table where the position of the differents variables is defined through occupation rate (iv,i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_height        ! Obstruction height within the domain (iv,i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_oai           ! Obstruction area index (iv,i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_fracxy        ! Obstruction correction factor for horizontal coverage (iv,i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_a2d           ! Horizontal area occupied by obstructions per unit area (iv+3,i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_s2d           ! Vertical area occupied by obstructions per unit area (iv+3,i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC :: obst_z0obst        ! z0 for each obstructions + for no turb (iv+3,i,j), [m]


   ! Variables on (iv,k,i,j)
   !------------------------
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_dens3d      ! Obstructions density for use with 3D formulations (iv,k,i,j), [m-2]
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_width3d     ! Obstructions width (in flow direction) for use with 3D formulations (iv,k,i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_thick3d     ! Obstructions thickness (perpendicular to flow direction) for use with 3D formulations (iv,k,i,j), [m]
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_theta3d     ! Obstruction bending angle (from vertical) (iv,k,i,j), [rad]
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_fracz3d     ! Fraction of sigma layer occupied by obstruction (iv,k,i,j) [0:1]
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_drag3d      ! Obstructions drag coefficient used by 3D formulations (iv,k,i,j), [-]

   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_a3d         ! Horizontal area occupied by obstructions at height k per unit area (iv+3,k,i,j), [-]
   REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE,PUBLIC :: obst_s3d         ! Vertical area occupied by obstructions at height k per unit area (iv+3,k,i,j), [-]

   ! Variables on (iv,kk)
   !---------------------
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_dens_norm       ! Normalized density from the distribution file (iv,kk), allocated within obst_readfile_char subroutine, [-]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_height_norm     ! Normalized heigth from the distribution file (iv,kk), allocated within obst_readfile_char subroutine, [-]

   ! Variables on (iv,t)
   !--------------------
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_dens_t          ! Instantaneous obstruction density (iv,t), allocated within obst_readfile_char subroutine, [-]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_width_t         ! Instantaneous obstruction width (iv,t), allocated within obst_readfile_char subroutine, [m]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_thick_t         ! Instantaneous obstruction thickness (iv,t), allocated within obst_readfile_char subroutine, [m]
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_height_t        ! Instantaneous obstruction heights (iv,t), allocated within obst_readfile_char subroutine, [m]
   ! Variables on (iv,other)
   !------------------------
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_fx         ! Horizontal force acting on leaves seagment (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_fz         ! Vertical force acting on leaves seagment (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_zcent      ! Z position of segment centres (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_t0cent     ! First bending angles of segment centres (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_t1cent     ! Second bending angles of segment centres (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_dtheta     ! Difference in bending angles of segment centres between two iterations (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_uvcent     ! Local velocity at segment centres (Abdelrhman method for bending) (iv,obst_c_abdel_nmax)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_zn         ! Height of end of each segment (Abdelrhman method for bending) (iv,0:obst_abdel_c_nmax+1)
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: obst_abdel_tn         ! Bending angle of end of each segment (Abdelrhman method for bending) (iv,0:obst_c_abdel_nmax+1)

   !! * Private variables

   CONTAINS

   !!==========================================================================================================

SUBROUTINE OBSTRUCTIONS_alloc_nbvar

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_alloc_nbvar  ***
   !&E
   !&E ** Purpose : Allocation of tables depending on number of obstructions variables
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : casxxx or obst_readvar
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08    (F. Ganthy) Add variables pre-initialization
   !&E       ! 2014-10    (F. Ganthy) More modifications + computation of obstructions
   !&E                                posture following Abdelrhman 2007
   !&E       ! 2016-08    (F. Ganthy) Optimization on Abdelrhman 2007 method
   !&E       ! 2017-02-16 (F. Ganthy) Some modifications:
   !&E                                - Allowing multiple obstructions type in a single grid cell
   !&E                                  --> allowed multispecific computation.
   !&E                                  This imply that some tables must be allocated depending on (iv,k,i,j) or (iv,i,j)
   !&E                                  These allocations are done within obst_alloc_xyz (because tables allocated within
   !&E                                  obst_alloc_nbvar are only those read within namelist).
   !&E                                - Changes on instantaneous obstruction state variables for future coupling with Zostera growth module
   !&E                                - Differenciation of cylindric / parallelepipedic structures
   !&E                                - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell)
   !&E                                - Cleaning (removing useless parameters and tests)
   !&E       ! 2017-04 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh

   IMPLICIT NONE

   !! * Local declaration

   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBST_ALLOC_NBVAR'
   !-------------------------
   ! Variables on (iv)
   !--------------------
   ALLOCATE(obst_varnum               (1:obst_nbvar))
   obst_varnum(:)                     = 0
   ALLOCATE(obst_fracxy_type          (1:obst_nbvar))
   obst_fracxy_type(:)                = 0
   ALLOCATE(obst_nbhnorm              (1:obst_nbvar))
   obst_nbhnorm(:)                    = 0
   ALLOCATE(obst_c_abdel_nmax         (1:obst_nbvar))
   obst_c_abdel_nmax(:)               = 1

   ALLOCATE(l_obst_filechar           (1:obst_nbvar))
   l_obst_filechar(:)                 = .FALSE.
   ALLOCATE(l_obst_filedistri         (1:obst_nbvar))
   l_obst_filedistri(:)               = .FALSE.
   ALLOCATE(l_obst_init_spatial       (1:obst_nbvar))
   l_obst_init_spatial(:)             = .FALSE.
   ALLOCATE(l_obst_flexible           (1:obst_nbvar))
   l_obst_flexible(:)                 = .FALSE.
   ALLOCATE(l_obst_cylindre           (1:obst_nbvar))
   l_obst_cylindre(:)                 = .FALSE.
   ALLOCATE(l_obst_downward           (1:obst_nbvar))
   l_obst_downward(:)                 = .FALSE.
   ALLOCATE(l_obst_3dobst             (1:obst_nbvar))
   l_obst_3dobst(:)                   = .FALSE.
   ALLOCATE(l_obst_noturb             (1:obst_nbvar))
   l_obst_noturb(:)                   = .FALSE.
   ALLOCATE(l_obst_abdelrough_cste    (1:obst_nbvar))
   l_obst_abdelrough_cste             = .FALSE.
   ALLOCATE(l_obst_fracxy             (1:obst_nbvar))
   l_obst_fracxy(:)                    = .FALSE.
   ALLOCATE(l_obst_abdelposture       (1:obst_nbvar))
   l_obst_abdelposture(:)             = .FALSE.
   ALLOCATE(l_obst_param_height       (1:obst_nbvar))
   l_obst_param_height(:)             = .FALSE.
   ALLOCATE(l_obst_drag_cste          (1:obst_nbvar))
   l_obst_drag_cste(:)                = .TRUE.
   ALLOCATE(l_obst_z0bstress          (1:obst_nbvar))
   l_obst_z0bstress(:)                = .FALSE.

   ALLOCATE(obst_z0bstress_option     (1:obst_nbvar))
   obst_z0bstress_option(:)           = 0

   ALLOCATE(obst_varname              (1:obst_nbvar))
   obst_varname(:)                    = '.'
   ALLOCATE(obst_fn_vardat            (1:obst_nbvar))
   obst_fn_vardat(:)                  = '.'
   ALLOCATE(obst_fn_char              (1:obst_nbvar))
   obst_fn_char(:)                    = '.'
   ALLOCATE(obst_fn_distrib           (1:obst_nbvar))
   obst_fn_distrib(:)                 = '.'

   ALLOCATE(obst_i_height             (1:obst_nbvar))
   obst_i_height(:)                   = 0.0_rsh
   ALLOCATE(obst_i_width              (1:obst_nbvar))
   obst_i_width(:)                    = 0.0_rsh
   ALLOCATE(obst_i_thick              (1:obst_nbvar))
   obst_i_thick(:)                    = 0.0_rsh
   ALLOCATE(obst_i_dens               (1:obst_nbvar))
   obst_i_dens(:)                     = 0.0_rsh

   ALLOCATE(obst_c_rho                (1:obst_nbvar))
   obst_c_rho(:)                      = 0.0_rsh
   ALLOCATE(obst_c_drag               (1:obst_nbvar))
   obst_c_drag(:)                     = 0.0_rsh
   ALLOCATE(obst_c_lift               (1:obst_nbvar))
   obst_c_lift(:)                     = 0.0_rsh
   ALLOCATE(obst_c_z0bstress          (1:obst_nbvar))
   obst_c_z0bstress(:)                = 0.0_rsh
   ALLOCATE(obst_c_fracxy_k0          (1:obst_nbvar))
   obst_c_fracxy_k0(:)                = 0.0_rsh
   ALLOCATE(obst_c_fracxy_k1          (1:obst_nbvar))
   obst_c_fracxy_k1(:)                = 0.0_rsh
   ALLOCATE(obst_c_fracxy_l           (1:obst_nbvar))
   obst_c_fracxy_l(:)                 = 0.0_rsh
   ALLOCATE(obst_c_crough_x0          (1:obst_nbvar))
   obst_c_crough_x0(:)                = 0.0_rsh
   ALLOCATE(obst_c_crough_x1          (1:obst_nbvar))
   obst_c_crough_x1(:)                = 0.0_rsh
   ALLOCATE(obst_c_lz                 (1:obst_nbvar))
   obst_c_lz(:)                       = 0.0_rsh
   ALLOCATE(obst_c_shelter            (1:obst_nbvar))
   obst_c_shelter(:)                  = 1.0_rsh
   ALLOCATE(obst_c_height_x0          (1:obst_nbvar))
   obst_c_height_x0                   = 0.0_rsh
   ALLOCATE(obst_c_height_x1          (1:obst_nbvar))
   obst_c_height_x1(:)                = 0.0_rsh
   ALLOCATE(obst_c_z0bstress_x0       (1:obst_nbvar))
   obst_c_z0bstress_x0(:)             = 0.0_rsh
   ALLOCATE(obst_c_z0bstress_x1       (1:obst_nbvar))
   obst_c_z0bstress_x1(:)             = 0.0_rsh
   ALLOCATE(obst_c_z0bstress_x2       (1:obst_nbvar))
   obst_c_z0bstress_x2(:)             = 0.0_rsh
   !-------------------------
   PRINT_DBG*, 'END OBSTRUCTIONS_ALLOC_NBVAR'
END SUBROUTINE OBSTRUCTIONS_alloc_nbvar

   !!==========================================================================================================

SUBROUTINE OBSTRUCTIONS_alloc_xyz

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_alloc_xyz  ***
   !&E
   !&E ** Purpose : Allocation of spatial obstruction tables
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : casxxx or obst_init
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08    (F. Ganthy) Add variables pre-initialization
   !&E       ! 2014-10    (F. Ganthy) More modifications + computation of obstructions
   !&E                                posture following Abdelrhman 2007
   !&E       ! 2016-03    (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2017-02-16 (F. Ganthy) Some modifications:
   !&E                                - Allowing multiple obstructions type in a single grid cell
   !&E                                  --> allowed multispecific computation.
   !&E                                  This imply that some tables must be allocated depending on (iv,k,i,j) or (iv,i,j)
   !&E                                  These allocations are done within obst_alloc_xyz (because tables allocated within
   !&E                                  obst_alloc_nbvar are only those read within namelist).
   !&E                                - Changes on instantaneous obstruction state variables for future coupling with Zostera growth module
   !&E                                - Differenciation of cylindric / parallelepipedic structures
   !&E                                - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell)
   !&E                                - Cleaning (removing useless parameters and tests)
   !&E       ! 2017-04    (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04    (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-11    (F. Ganthy) Change order of initializations
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh,limin,liminm1,limax,ljmin,ljminm1,ljmax,kmax
   USE comvars2d,  ONLY  : l_modele2d

   IMPLICIT NONE

   !! * Local declaration

   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_ALLOC_XYZ'
   !------------------------------------
   ! Definition of effective kmax to use
   !------------------------------------
   IF((l_modele2d).OR.(kmax.EQ.1))THEN
       obst_kmax = obst_kmax2d
   ELSE
       obst_kmax = kmax
   ENDIF
   !----------------------
   ! Variables on (k)
   !----------------------
   ALLOCATE(obst_sig             (0:obst_kmax+1))
   obst_sig(:)                   = 0.0_rsh
   ALLOCATE(obst_dsig            (obst_kmax))
   obst_dsig(:)                  = 0.0_rsh
   !----------------------
   ! Variables on (iv,i,j)
   !----------------------
   ALLOCATE(obst_dens_inst       (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_dens_inst(:,:,:)         = 0.0_rsh
   ALLOCATE(obst_width_inst      (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_width_inst(:,:,:)        = 0.0_rsh
   ALLOCATE(obst_thick_inst      (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_thick_inst(:,:,:)        = 0.0_rsh
   ALLOCATE(obst_height_inst     (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_height_inst(:,:,:)       = 0.0_rsh
   ALLOCATE(obst_area_index_inst (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_area_index_inst(:,:,:)   = 0.0_rsh

   ALLOCATE(obst_position        (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_position(:,:,:)          = 0.0_rsh
   ALLOCATE(obst_height          (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_height(:,:,:)            = 0.0_rsh
   ALLOCATE(obst_oai             (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_oai(:,:,:)               = 0.0_rsh
   ALLOCATE(obst_fracxy          (1:obst_nbvar,limin:limax,ljmin:ljmax))
   obst_fracxy(:,:,:)            = 0.0_rsh

   ALLOCATE(obst_a2d             (1:obst_nbvar+3,limin:limax,ljmin:ljmax))
   obst_a2d(:,:,:)               = 0.0_rsh
   ALLOCATE(obst_s2d             (1:obst_nbvar+3,limin:limax,ljmin:ljmax))
   obst_s2d(:,:,:)               = 0.0_rsh
   ALLOCATE(obst_z0obst          (1:obst_nbvar+3,limin:limax,ljmin:ljmax))
   obst_z0obst(:,:,:)            = 0.0_rsh

   !-------------------
   ! Variables on (i,j)
   !-------------------
   ALLOCATE(obst_roswat_bot      (limin:limax,ljmin:ljmax))
   obst_roswat_bot               = 0.0_rsh
   ALLOCATE(obst_fu_i            (limin:limax,ljmin:ljmax))
   obst_fu_i(:,:)                = 0.0_rsh
   ALLOCATE(obst_fv_i            (limin:limax,ljmin:ljmax))
   obst_fv_i(:,:)                = 0.0_rsh
   ALLOCATE(obst_fu_e            (limin:limax,ljmin:ljmax))
   obst_fu_e(:,:)                = 0.0_rsh
   ALLOCATE(obst_fv_e            (limin:limax,ljmin:ljmax))
   obst_fv_e(:,:)                = 0.0_rsh
   ALLOCATE(obst_z0bed           (limin:limax,ljmin:ljmax))
   obst_z0bed(:,:)               = 0.0_rsh
   ALLOCATE(obst_bstress         (limin:limax,ljmin:ljmax))
   obst_bstress(:,:)             = 0.0_rsh
   ALLOCATE(obst_bstressc        (limin:limax,ljmin:ljmax))
   obst_bstressc(:,:)            = 0.0_rsh
   ALLOCATE(obst_bstressw        (limin:limax,ljmin:ljmax))
   obst_bstressw(:,:)            = 0.0_rsh
   ALLOCATE(obst_z0bstress       (limin:limax,ljmin:ljmax))
   obst_z0bstress(:,:)           = 0.0_rsh
   ALLOCATE(obst_raphbx          (liminm1:limax,ljmin:ljmax))
   obst_raphbx(:,:)              = 0.0_rsh
   ALLOCATE(obst_raphby          (limin:limax,ljminm1:ljmax))
   obst_raphby(:,:)              = 0.0_rsh
    ALLOCATE(obst_frofonx        (liminm1:limax,ljmin:ljmax))
   obst_frofonx(:,:)             = 0.0_rsh
   ALLOCATE(obst_frofony         (limin:limax,ljminm1:ljmax))
   obst_frofony(:,:)             = 0.0_rsh
   ALLOCATE(obst_dens_mean       (limin:limax,ljminm1:ljmax))
   obst_dens_mean(:,:)           = 0.0_rsh
   ALLOCATE(obst_width_mean      (limin:limax,ljminm1:ljmax))
   obst_width_mean(:,:)          = 0.0_rsh
   ALLOCATE(obst_height_mean     (limin:limax,ljminm1:ljmax))
   obst_height_mean(:,:)         = 0.0_rsh
   !------------------------
   ! Variables on (iv,k,i,j)
   !------------------------
   ALLOCATE(obst_dens3d          (1:obst_nbvar,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_dens3d(:,:,:,:)          = 0.0_rsh
   ALLOCATE(obst_width3d         (1:obst_nbvar,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_width3d(:,:,:,:)         = 0.0_rsh
   ALLOCATE(obst_thick3d         (1:obst_nbvar,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_thick3d(:,:,:,:)         = 0.0_rsh
   ALLOCATE(obst_theta3d         (1:obst_nbvar,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_theta3d(:,:,:,:)         = 0.0_rsh
   ALLOCATE(obst_fracz3d         (1:obst_nbvar,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_fracz3d(:,:,:,:)         = 0.0_rsh
   ALLOCATE(obst_drag3d          (1:obst_nbvar,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_drag3d(:,:,:,:)          = 0.0_rsh

   ALLOCATE(obst_a3d             (1:obst_nbvar+3,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_a3d(:,:,:,:)             = 0.0_rsh
   ALLOCATE(obst_s3d             (1:obst_nbvar+3,1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_s3d(:,:,:,:)             = 0.0_rsh
   !---------------------
   ! Variables on (k,i,j)
   !---------------------
   ALLOCATE(obst_zc              (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_zc(:,:,:)                = 0.0_rsh
   ALLOCATE(obst_dz              (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_dz(:,:,:)                = 0.0_rsh
   ALLOCATE(obst_uz              (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_uz(:,:,:)                = 0.0_rsh
   ALLOCATE(obst_vz              (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_vz(:,:,:)                = 0.0_rsh
   ALLOCATE(obst_fuz_i           (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_fuz_i(:,:,:)             = 0.0_rsh
   ALLOCATE(obst_fvz_i           (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_fvz_i(:,:,:)             = 0.0_rsh
   ALLOCATE(obst_fuz_e           (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_fuz_e(:,:,:)             = 0.0_rsh
   ALLOCATE(obst_fvz_e           (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_fvz_e(:,:,:)             = 0.0_rsh
   ALLOCATE(obst_t               (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_t(:,:,:)                 = 0.0_rsh
   ALLOCATE(obst_tau             (1:obst_kmax,limin:limax,ljmin:ljmax))
   obst_tau(:,:,:)               = 0.0_rsh
   !------------------------------
   PRINT_DBG*, 'END OBSTRUCTIONS_ALLOC_XYZ'
END SUBROUTINE OBSTRUCTIONS_alloc_xyz

   !!==========================================================================================================

SUBROUTINE OBSTRUCTIONS_alloc_other

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_alloc_other  ***
   !&E
   !&E ** Purpose : Allocation of other obstruction tables
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : casxxx or obst_init
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2016-08-16 (F. Ganthy) Original code
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh

   IMPLICIT NONE

   !! * Local declaration

   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_ALLOC_OTHER'
   !--------------------------
   ! OTHER : Abdelrhman method
   ALLOCATE(obst_abdel_fx     (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_fx (:,:)        = 0.0_rsh
   ALLOCATE(obst_abdel_fz     (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_fz(:,:)         = 0.0_rsh
   ALLOCATE(obst_abdel_zcent  (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_zcent(:,:)      = 0.0_rsh
   ALLOCATE(obst_abdel_t0cent (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_t0cent(:,:)     = 0.0_rsh
   ALLOCATE(obst_abdel_t1cent (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_t1cent(:,:)     = 0.0_rsh
   ALLOCATE(obst_abdel_dtheta (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_dtheta(:,:)     = 0.0_rsh
   ALLOCATE(obst_abdel_uvcent (1:obst_nbvar,1:MAXVAL(obst_c_abdel_nmax)))
   obst_abdel_uvcent(:,:)     = 0.0_rsh
   ALLOCATE(obst_abdel_zn     (1:obst_nbvar,0:MAXVAL(obst_c_abdel_nmax)+1))
   obst_abdel_zn(:,:)         = 0.0_rsh
   ALLOCATE(obst_abdel_tn     (1:obst_nbvar,0:MAXVAL(obst_c_abdel_nmax)+1))
   obst_abdel_tn(:,:)         = 0.0_rsh
   !--------------------------
   PRINT_DBG*, 'END OBSTRUCTIONS_ALLOC_OTHER'
END SUBROUTINE OBSTRUCTIONS_alloc_other

   !!==========================================================================================================

SUBROUTINE OBSTRUCTIONS_dealloc

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_dealloc  ***
   !&E
   !&E ** Purpose : Deallocation of variables
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : main.F90
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2016-03-11 (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2016-08    (F. Ganthy) Optimization on Abdelrhman 2007 method
   !&E       ! 2017-02-16 (F. Ganthy) Some modifications:
   !&E                                - Allowing multiple obstructions type in a single grid cell
   !&E                                  --> allowed multispecific computation.
   !&E                                  This imply that some tables must be allocated depending on (iv,k,i,j) or (iv,i,j)
   !&E                                  These allocations are done within obst_alloc_xyz (because tables allocated within
   !&E                                  obst_alloc_nbvar are only those read within namelist).
   !&E                                - Changes on instantaneous obstruction state variables for future coupling with Zostera growth module
   !&E                                - Differenciation of cylindric / parallelepipedic structures
   !&E                                - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell)
   !&E                                - Cleaning (removing useless parameters and tests)
   !&E       ! 2017-04    (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04    (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh,limin,limax,jmin,jmax,kmax,liminm1,ljminm1

   IMPLICIT NONE

   !! * Local declaration

   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_DEALLOC'
   !-------------------------
   ! Variables on (iv)
   !--------------------
   DEALLOCATE(obst_varnum)
   DEALLOCATE(obst_nbhnorm)
   DEALLOCATE(obst_fracxy_type)
   DEALLOCATE(obst_c_abdel_nmax)

   DEALLOCATE(l_obst_filechar)
   DEALLOCATE(l_obst_filedistri)
   DEALLOCATE(l_obst_init_spatial)
   DEALLOCATE(l_obst_flexible)
   DEALLOCATE(l_obst_cylindre)
   DEALLOCATE(l_obst_downward)
   DEALLOCATE(l_obst_3dobst)
   DEALLOCATE(l_obst_noturb)
   DEALLOCATE(l_obst_abdelrough_cste)
   DEALLOCATE(l_obst_fracxy)
   DEALLOCATE(l_obst_abdelposture)
   DEALLOCATE(l_obst_param_height)
   DEALLOCATE(l_obst_drag_cste)
   DEALLOCATE(l_obst_z0bstress)

   DEALLOCATE(obst_z0bstress_option)

   DEALLOCATE(obst_varname)
   DEALLOCATE(obst_fn_vardat)
   DEALLOCATE(obst_fn_char)
   DEALLOCATE(obst_fn_distrib)

   DEALLOCATE(obst_i_height)
   DEALLOCATE(obst_i_width)
   DEALLOCATE(obst_i_thick)
   DEALLOCATE(obst_i_dens)

   DEALLOCATE(obst_c_rho)
   DEALLOCATE(obst_c_drag)
   DEALLOCATE(obst_c_lift)
   DEALLOCATE(obst_c_z0bstress)
   DEALLOCATE(obst_c_fracxy_k0)
   DEALLOCATE(obst_c_fracxy_k1)
   DEALLOCATE(obst_c_fracxy_l)
   DEALLOCATE(obst_c_crough_x0)
   DEALLOCATE(obst_c_crough_x1)
   DEALLOCATE(obst_c_lz)
   DEALLOCATE(obst_c_shelter)
   DEALLOCATE(obst_c_height_x0)
   DEALLOCATE(obst_c_height_x1)
   DEALLOCATE(obst_c_z0bstress_x0)
   DEALLOCATE(obst_c_z0bstress_x1)
   DEALLOCATE(obst_c_z0bstress_x2)

   !----------------------
   ! Variables on (k)
   !----------------------
   DEALLOCATE(obst_sig)
   DEALLOCATE(obst_dsig)
   !----------------------
   ! Variables on (iv,i,j)
   !----------------------
   DEALLOCATE(obst_position)

   DEALLOCATE(obst_dens_inst)
   DEALLOCATE(obst_width_inst)
   DEALLOCATE(obst_thick_inst)
   DEALLOCATE(obst_height_inst)
   DEALLOCATE(obst_area_index_inst)

   DEALLOCATE(obst_height)
   DEALLOCATE(obst_oai)
   DEALLOCATE(obst_fracxy)

   DEALLOCATE(obst_a2d)
   DEALLOCATE(obst_s2d)
   DEALLOCATE(obst_z0obst)

   !-------------------
   ! Variables on (i,j)
   !-------------------
   DEALLOCATE(obst_fu_i)
   DEALLOCATE(obst_fv_i)
   DEALLOCATE(obst_fu_e)
   DEALLOCATE(obst_fv_e)
   DEALLOCATE(obst_z0bed)
   DEALLOCATE(obst_bstress)
   DEALLOCATE(obst_bstressc)
   DEALLOCATE(obst_bstressw)
   DEALLOCATE(obst_z0bstress)

   DEALLOCATE(obst_raphbx)
   DEALLOCATE(obst_raphby)
   DEALLOCATE(obst_frofonx)
   DEALLOCATE(obst_frofony)

   DEALLOCATE(obst_dens_mean)
   DEALLOCATE(obst_width_mean)
   DEALLOCATE(obst_height_mean)
   !------------------------
   ! Variables on (iv,k,i,j)
   !------------------------
   DEALLOCATE(obst_dens3d)
   DEALLOCATE(obst_width3d)
   DEALLOCATE(obst_thick3d)
   DEALLOCATE(obst_theta3d)
   DEALLOCATE(obst_fracz3d)
   DEALLOCATE(obst_drag3d)

   DEALLOCATE(obst_a3d)
   DEALLOCATE(obst_s3d)
   !---------------------
   ! Variables on (k,i,j)
   !---------------------
   DEALLOCATE(obst_zc)
   DEALLOCATE(obst_dz)
   DEALLOCATE(obst_uz)
   DEALLOCATE(obst_vz)
   DEALLOCATE(obst_fuz_i)
   DEALLOCATE(obst_fvz_i)
   DEALLOCATE(obst_fuz_e)
   DEALLOCATE(obst_fvz_e)
   DEALLOCATE(obst_t)
   DEALLOCATE(obst_tau)
   !---------------------
   ! Variables on (other)
   !---------------------
   DEALLOCATE(obst_abdel_fx)
   DEALLOCATE(obst_abdel_fz)
   DEALLOCATE(obst_abdel_zcent)
   DEALLOCATE(obst_abdel_t0cent)
   DEALLOCATE(obst_abdel_t1cent)
   DEALLOCATE(obst_abdel_dtheta)
   DEALLOCATE(obst_abdel_uvcent)
   DEALLOCATE(obst_abdel_zn)
   DEALLOCATE(obst_abdel_tn)

   DEALLOCATE(obst_dens_norm)
   DEALLOCATE(obst_height_norm)
   DEALLOCATE(obst_dens_t)
   DEALLOCATE(obst_width_t)
   DEALLOCATE(obst_thick_t)
   DEALLOCATE(obst_height_t)
   !-------------------------
   PRINT_DBG*, 'END OBSTRUCTIONS_DEALLOC'
END SUBROUTINE OBSTRUCTIONS_dealloc

   !!==========================================================================================================



#endif
 END MODULE comOBSTRUCTIONS
