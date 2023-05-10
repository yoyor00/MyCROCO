#include "cppdefs.h"
 
MODULE comobstructions

#ifdef OBSTRUCTION
   !&E======================================================================
   !&E                   ***  MODULE comOBSTRUCTIONS  ***
   !&E
   !&E ** Purpose : Declaration of all variables used by obstructions module
   !&E
   !&E ** History :
   !&E     ! 2018-04-13 (F. Ganthy) Original code : moved from obstructions.F90
   !&E     ! 2022-01-11 (A. Le Pevedic) Introduction of averaged parameters to send to WW3 
   !&E     ! 2022-12 (S. Le Gac) Adaptation for CROCO, move alloc/dealloc to initOBSTRUCTIONS.F90
   !&E     ! 2023-04 (S. Le Gac) Adaptation for CROCO, rename ALL variables with obst_ at the begining
   !&E
   !&E======================================================================


   !! * Declaration
   IMPLICIT NONE

   ! default
   public


   ! -------------------------------------------------------------------------
   ! Definition of rsh, rlg, riosh, riolg, lchain
   ! -------------------------------------------------------------------------
   INTEGER,PARAMETER :: riosh = 8, riolg = 8, rlg = 8, rsh = 8
   INTEGER,PARAMETER :: lchain = 200
   INTEGER,PARAMETER :: pi = 3.14159265358979323846
   REAL,PARAMETER :: obst_c2turb = 1.92  ! beta2 in gls kepsilon TODO : add compatibility test to be sure GLS_KEPSILON is used and take beta2 value at initialisation

   INTEGER,PUBLIC :: imin, imax, jmin, jmax, kmax 
   INTEGER,PUBLIC :: ierrorlog, iwarnlog, iscreenlog 
   REAL,PUBLIC :: obst_h0fond
   REAL,PUBLIC :: obst_cmu


   !! * Shared or public module variables

   ! Variables nomenclature :
   ! obst_*    : variables for obstructions
   ! obst_i_*  : variables used for initialization
   ! obst_c_*  : constant parameters
   ! obst_fn_* : filename variables
   ! obst_l_*  : logical variables

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

   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_pos                             ! Name obstruction position (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_height_f                        ! Name 2D obstruction height (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_height_e                        ! Name 2D obstruction height (effective) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_dens_f                          ! Name 2D obstruction density (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_dens_e                          ! Name 3D obstruction density (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_width_f                         ! Name 2D obstruction width (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_width_e                         ! Name 3D obstruction width (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_thick_f                         ! Name 2D obstruction thick (forcing) (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_thick_e                         ! Name 3D obstruction thick (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_oai                             ! Name 2D obstruction area index, (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_theta                           ! Name 3D obstruction bending angle (iv,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_cover                           ! Name 2D obstruction coverage, (iv,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_frac_z                          ! Name 3D obstruction fraction of sigma layer, (iv,i,j)


   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_fuv                             ! Name 2D obstruction friction force (i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_fuzvz                           ! Name 3D obstruction friction force (k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_a2d                             ! Name 2D obstruction horizontal area (iv+3,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_a3d                             ! Name 3D obstruction horizontal area (iv+3,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_s2d                             ! Name 2D obstruction vertical area (iv+3,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_s3d                             ! Name 3D obstruction vertical area (iv+3,k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_drag                            ! Name 3D obstruction drag coefficient (k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_tau                             ! Name 3D obstruction turbulence dissipation scale (k,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_z0bed                           ! Name 2D bottom roughness length for bed (i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_z0obst                          ! Name 2D bottom roughness length for obstruction (iv+3,i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_z0bstress                       ! Name 2D bottom roughness length used for bottom shear stress computation (i,j)
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_bstress                         ! Name 2D total bottom shear stress within output file
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_bstressc                        ! Name 2D current bottom shear stress within output file
   CHARACTER(LEN=lchain),PUBLIC :: obst_nout_bstressw                        ! Name 2D wave bottom shear stress within output file

   REAL(KIND=riosh),PUBLIC :: obst_riog_min_pos,obst_riog_max_pos         ! Valid minimum and maximum for output obstruction position
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_height,obst_riog_max_height   ! Valid minimum and maximum for output obstruction height
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_dens,obst_riog_max_dens       ! Valid minimum and maximum for output obstruction density
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_width,obst_riog_max_width     ! Valid minimum and maximum for output obstruction width
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_thick,obst_riog_max_thick     ! Valid minimum and maximum for output obstruction tickness
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_oai,obst_riog_max_oai         ! Valid minimum and maximum for output obstruction area index
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_theta,obst_riog_max_theta     ! Valid minimum and maximum for output obstruction bending angle
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_cover,obst_riog_max_cover     ! Valid minimum and maximum for output obstruction coverage
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_fracz,obst_riog_max_fracz     ! Valid minimum and maximum for output obstruction fraction in z
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_fuv,obst_riog_max_fuv         ! Valid minimum and maximum for output obstructin friction force
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_a,obst_riog_max_a             ! Valid minimum and maximum for output obstruction horizontal area
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_s,obst_riog_max_s             ! Valid minimum and maximum for output obstruction vertical area
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_drag,obst_riog_max_drag       ! Valid minimum and maximum for output obstruction drag
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_tau,obst_riog_max_tau         ! Valid minimum and maximum for output obstruction tau
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_z0,obst_riog_max_z0           ! Valid minimum and maximum for output obstruction roughness length
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_bstress,obst_riog_max_bstress ! Valid minimum and maximum for output bottom shear stress
   REAL(KIND=riosh),PUBLIC :: obst_riog_min_zroot,obst_riog_max_zroot     ! Valid minimum and maximum for output obstruction root depth

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

   LOGICAL,PUBLIC :: obst_l_z0bstress_tot              ! IF ONLY ONE obstruction VARIABLE USED Z0SED

   REAL(KIND=rsh),PUBLIC  :: obst_c_paramhuv           ! The coefficient of obstruction height for computation of velocity
   REAL(KIND=rsh),PUBLIC  :: obst_c_imp3d              ! Implicitation coefficient for obstructions formulations in 3D
   REAL(KIND=rsh),PUBLIC  :: obst_c_exp3d              ! Explicitation coefficient for obstructions formulations in 3D
   REAL(KIND=rsh),PUBLIC  :: obst_c_imp2d              ! Implicitation coefficient for obstructions formulations in 2D
   REAL(KIND=rsh),PUBLIC  :: obst_c_exp2d              ! Explicitation coefficient for obstructions formulations in 2D

!FG settling velocity term for zostera (should be defined as iv dependant)
   REAL(KIND=rsh),PUBLIC  :: obst_c_settletrapp_x0   ! first term for parameterized obstructions trapping factor
   REAL(KIND=rsh),PUBLIC  :: obst_c_settletrapp_x1   ! second term for parameterized obstructions trapping factor 
   REAL(KIND=rsh),PUBLIC  :: obst_c_settleblock_x0   ! first term for parameterized obstructions blocking factor
   REAL(KIND=rsh),PUBLIC  :: obst_c_settleblock_x1   ! second term for parameterized obstructions blocking factor
!END FG

   REAL(KIND=rsh),PUBLIC  :: obst_i_z0bstress          ! Roughness length for bottom shear stress without obstruction (= z0seduni if key_sedim)

   INTEGER,PARAMETER,PUBLIC         :: obst_kmax2d = 20       ! Number of sigma layer used for modele2D
   REAL(KIND=rsh),PARAMETER,PUBLIC  :: obst_p_hmin = 0.05_rsh ! Minimum coefficient value for obstruction height under bending  

   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ONLY ON THE NUMBER OF OBSTRUCTIONS VARIABLES (iv) :
   !--------------------------------------------------------------------------

   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_varnum                   ! The number allocated to each variables (iv)
   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_fracxy_type              ! The type of correction for horizontal coverage ((grid cell not completely fill with obstructions) (iv)
   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_nbhnorm                  ! Number of vertical steps from the distribution file (iv)
   INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_c_abdel_nmax             ! Number of segments for Abdlerhman method (bending) (iv) 

   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_filechar               ! For reading time-series of obstructions characteristics from a file (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_filedistri             ! For reading the vertical distribution of obstructions from a file (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_init_spatial           ! For reading spatial file for density, height, width and thickness (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_flexible               ! For obstructions flexibility (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_cylindre               ! For obstruction shape (cylinder, ellispe / parallelepipeds), (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_downward               ! For downward obstructions (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_3dobst                 ! For fully 3D obstructions (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_noturb                 ! For use of simplified formulation (roughness length) instead of full turbulent formulation (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_abdelrough_cste        ! For use constant Abdelrhman 2003 coefficient (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_fracxy                 ! For horizontal coverage correction (grid cell not completely fill with obstructions) (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_abdelposture           ! For computation of obstructions posture following abdelrhman 2007 (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_param_height           ! For use a parameterization (from velocity) for obstruction height (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_drag_cste              ! For constant/variable drag coefficient (or corrected from bending angle) (iv)
   LOGICAL,DIMENSION(:),ALLOCATABLE,PUBLIC :: obst_l_z0bstress              ! For using a roughness length for obstructions (iv)

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
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE,PUBLIC   :: obst_dsig            ! Sigma layer thickness (same as sig (comvars3d) for 3D modele, or allocated on obst_kmax2d for 2D modele)
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

#endif

END MODULE comOBSTRUCTIONS
