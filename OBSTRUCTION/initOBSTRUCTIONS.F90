MODULE initobstructions

#include "cppdefs.h"

#ifdef OBSTRUCTION
   !&E==========================================================================
   !&E                   ***  MODULE  initOBSTRUCTIONS  ***
   !&E
   !&E
   !&E ** Purpose : concerns all subroutines related to obstructions interactions
   !&E    with hydrodynamics
   !&E
   !&E ** Description :
   !&E     
   !&E     subroutine OBSTRUCTIONS_init               ! Defines obstructions outputs parameters
   !&E                                                ! and initialize some variables
   !&E
   !&E     subroutine OBSTRUCTIONS_readvar            ! Read the file (variables_obstructions.dat)
   !&E                                                ! for initialization
   !&E
   !&E     subroutine OBSTRUCTIONS_write_summary      ! Perform some checks and write a summary 
   !&E                                                ! of obstructions variables
   !&E
   !&E     subroutine OBSTRUCTIONS_compatibility      ! Check compatibility of parameterization used
   !&E
   !&E     subroutine OBSTRUCTIONS_readfile_char      ! Reads input files (time-series) for 
   !&E                                                ! obstructions characteristics
   !&E
   !&E     subroutine OBSTRUCTIONS_readfile_pos       ! Reads input file for obstruction position
   !&E                                                ! within the domain
   !&E
   !&E     subroutine OBSTRUCTIONS_readfile_distri    ! Read input file for obstructions vertical
   !&E                                                ! distribution (normalized in height and density)
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
   !&E     ! 2021-10 (F. Ganthy) Original code extracted from OBSTRUCTION.F90
   !&E     ! 2022-12 (S. Le Gac) Adaptation for CROCO, extraction of alloc/dealloc from comOBSTRUCTIONS.F90
   !&E
   !&E==========================================================================
   !! * Modules used
   USE comOBSTRUCTIONS

   IMPLICIT NONE

   !! * Accessibility

   ! function & routines of this module, called outside :
   ! PUBLIC functions
   PUBLIC OBSTRUCTIONS_init_dimension
   PUBLIC OBSTRUCTIONS_init
   PUBLIC OBSTRUCTIONS_readfile_char ! used in OBSTRUCTIONS

   PRIVATE

  CONTAINS

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_init_dimension(imin_in, imax_in, jmin_in, jmax_in, N, stdout)

    INTEGER,INTENT(IN) :: imin_in, imax_in, jmin_in, jmax_in, N, stdout

    imin = imin_in
    imax = imax_in
    jmin = jmin_in
    jmax = jmax_in
    kmax = N
    iscreenlog = stdout
    ierrorlog = stdout
    iwarnlog = stdout

    !!#TODO check mynode in MPI

   END SUBROUTINE OBSTRUCTIONS_init_dimension

   SUBROUTINE OBSTRUCTIONS_init(z0b, h0fond_in)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_init  ***
   !&E
   !&E ** Purpose : Initialization of obstruction parameters
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08    (F. Ganthy) Routine simplification, removing useless
   !&E                                initializations
   !&E       ! 2014-10    (F. Ganthy) Add bending angle
   !&E       ! 2015-09    (F. Ganthy) Minor correction (obst_cmu)
   !&E       ! 2016-03    (F. Ganthy) Add logical to choose turbulence coefficient
   !&E                                values (default or user-defined)
   !&E       ! 2017-04    (F. Ganthy) Change initialization of turbulence coefficients cmu and c2turb
   !&E       ! 2017-11    (F. Ganthy) Remove useless variables related to turbulence coefficients
   !&E       ! 2021-10    (F. Ganthy) Few modification/organization on paraOBSTRUCTIONS
   !&E       ! 2022-01    (A. Le Pevedic) Added wave friction factor (used for shear stress computation)
   !&E---------------------------------------------------------------------


   IMPLICIT NONE

   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: z0b 
   REAL(KIND=rsh),INTENT(IN) :: h0fond_in


   !! * Local declaration
   CHARACTER(len=lchain)                   :: filepc
   INTEGER                                 :: iv,k
   CHARACTER(len=lchain) :: obst_fn_var1,obst_fn_var2,obst_fn_var3,      &
                            obst_fn_var4,obst_fn_var5,obst_fn_var6
   REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE :: sw
   !! * Namelist
   NAMELIST/obst_main/obst_nbvar,obst_fn_position,             & 
                      obst_fn_var1,obst_fn_var2,obst_fn_var3,  &
                      obst_fn_var4,obst_fn_var5,obst_fn_var6
   NAMELIST/obst_numerics/obst_i_z0bstress,obst_c_paramhuv,obst_c_imp3d,obst_c_imp2d,obst_fricwav
   NAMELIST/obst_output/obst_l_out_pos,obst_l_out_height_f,obst_l_out_height_e,obst_l_out_dens_f,     &
                        obst_l_out_dens_e,obst_l_out_width_f,obst_l_out_width_e,                      &
                        obst_l_out_thick_f,obst_l_out_thick_e,obst_l_out_oai,                         &
                        obst_l_out_theta,obst_l_out_cover,obst_l_out_frac_z,                          &
                        obst_l_out_fuv,obst_l_out_fuzvz,obst_l_out_a2d,obst_l_out_a3d,obst_l_out_s2d, &
                        obst_l_out_s3d,obst_l_out_drag,obst_l_out_tau,                                &
                        obst_l_out_z0bed,obst_l_out_z0obst,obst_l_out_z0bstress,obst_l_out_bstress,   &
                        obst_l_out_bstressc,obst_l_out_bstressw
   !!----------------------------------------------------------------------
   !! * Executable part
    obst_h0fond = h0fond_in

   ! ************************
   ! * READING NAMELIST
   ! ************************
   filepc = './paraOBSTRUCTIONS.txt'
   OPEN(50,file=filepc,status='old',form='formatted',access='sequential')
   READ(50,obst_main)
   READ(50,obst_numerics)
   READ(50,obst_output)
   CLOSE(50)
   MPI_master_only  WRITE(iscreenlog,*) ' '
   MPI_master_only  WRITE(iscreenlog,*) ' '
   MPI_master_only  WRITE(iscreenlog,*) '*****************************************************'
   MPI_master_only  WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBST_INIT *****'
   MPI_master_only  WRITE(iscreenlog,*) '*****************************************************'
   MPI_master_only  WRITE(iscreenlog,*) ' Reading file ',TRIM(filepc)
   MPI_master_only  WRITE(iscreenlog,*) '***********************'


   ! ************************************
   ! * ALLOCATION OF VARIABLES PARAMETERS
   ! ************************************
   CALL OBSTRUCTIONS_alloc_nbvar
   !=========================
   ! Names of parameter files
   !=========================
   DO iv=1,obst_nbvar
     IF(iv.EQ.1)THEN
       obst_fn_vardat(iv) = obst_fn_var1
     ELSEIF(iv.EQ.2)THEN
       obst_fn_vardat(iv) = obst_fn_var2
     ELSEIF(iv.EQ.3)THEN
       obst_fn_vardat(iv) = obst_fn_var3
     ELSEIF(iv.EQ.4)THEN
       obst_fn_vardat(iv) = obst_fn_var4
     ELSEIF(iv.EQ.5)THEN
       obst_fn_vardat(iv) = obst_fn_var5
     ELSEIF(iv.EQ.6)THEN
       obst_fn_vardat(iv) = obst_fn_var6
     ELSE
        MPI_master_only  WRITE(ierrorlog,*) ' '
        MPI_master_only  WRITE(ierrorlog,*) '*****************************************************************'
        MPI_master_only  WRITE(ierrorlog,*) '******* module OBSTRUCTIONS, subroutine OBSTRUCTIONS_INIT *******'
        MPI_master_only  WRITE(ierrorlog,*) '*****************************************************************'
        MPI_master_only  WRITE(ierrorlog,*) '!!! ERROR : Maximum number of variables >6                    !!!'
        MPI_master_only  WRITE(ierrorlog,*) '!!! Few changes must be performed in the code                 !!!'
        MPI_master_only  WRITE(ierrorlog,*) '!!! to increase number of allowed variables :                 !!!'
        MPI_master_only  WRITE(ierrorlog,*) '!!! --> Namelist contains (subroutine OBSTRUCTIONS_INIT)      !!!'
        MPI_master_only  WRITE(ierrorlog,*) '!!! --> Allocation of names of parameters files               !!!'
        MPI_master_only  WRITE(ierrorlog,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     ENDIF
   ENDDO
   ! ***********************
   ! * OUTPUTS VARIABLES
   ! ***********************
#ifdef MUSTANG
   obst_l_out_bstress      = .FALSE.
   obst_l_out_bstressc     = .FALSE.
   obst_l_out_bstressw     = .FALSE.
#endif
   obst_nout_pos           = 'pos'       ! Name obstruction position
   obst_nout_height_f      = 'height_f'  ! Name 2D obstruction height (forcing) (iv,i,j)
   obst_nout_height_e      = 'height_e'  ! Name 2D obstruction height (effective) (iv,i,j)
   obst_nout_dens_f        = 'dens_f'    ! Name 2D obstruction density (forcing) (iv,i,j)
   obst_nout_dens_e        = 'dens_e'    ! Name 3D obstruction density (3D effective) (iv,k,i,j)
   obst_nout_width_f       = 'width_f'   ! Name 2D obstruction width (forcing) (iv,i,j)
   obst_nout_width_e       = 'width_e'   ! Name 3D obstruction width (3D effective) (iv,k,i,j)
   obst_nout_thick_f       = 'thick_f'   ! Name 2D obstruction thick (forcing) (iv,i,j)
   obst_nout_thick_e       = 'thick_e'   ! Name 3D obstruction thick (3D effective) (iv,k,i,j)
   obst_nout_oai           = 'oai'       ! Name 2D obstruction area index (iv,i,j)
   obst_nout_theta         = 'theta'     ! Name 3D obstruction bending angle (iv,k,i,j)
   obst_nout_cover         = 'frac_xy'   ! Name 2D obstruction coverage (iv,i,j)
   obst_nout_frac_z        = 'frac_z'    ! Name 2D obstruction fraction of sigma layer (iv,i,j)

   obst_nout_fuv           = 'fuv'       ! Name 2D obstruction friction force (i,j)
   obst_nout_fuzvz         = 'fuzvz'     ! Name 3D obstruction friction force (k,i,j)
   obst_nout_a3d           = 'a2d'       ! Name 2D obstruction horizontal area (iv+3,i,j)
   obst_nout_a3d           = 'a3d'       ! Name 3D obstruction horizontal area (iv+3,k,i,j)
   obst_nout_s3d           = 's2d'       ! Name 2D obstruction vertical area (iv+3,i,j)
   obst_nout_s3d           = 's3d'       ! Name 3D obstruction vertical area (iv+3,k,i,j)
   obst_nout_drag          = 'cd3d'      ! Name 3D obstruction drag coefficient (iv,k,i,j)
   obst_nout_tau           = 'tau3d'     ! Name 3D obstruction turbulence dissipation scale (k,i,j)
   obst_nout_z0bed         = 'z0bed'     ! Name 2D bottom roughness length of bed (i,j)
   obst_nout_z0obst        = 'z0obst'    ! Name 2D bottom roughness length of obstructions (i,j)
   obst_nout_z0bstress     = 'z0bstress' ! Name 2D bottom roughness length used for bottom shear stress computation (i,j)
   obst_nout_bstress       = 'taub'      ! Name 2D total bottom shear stress
   obst_nout_bstressc      = 'taubc'     ! Name 2D current bottom shear stress
   obst_nout_bstressw      = 'taubw'     ! Name 2D wave bottom shear stress

   obst_riog_min_pos     =  0.0_riosh
   obst_riog_max_pos     =  1.0_riosh
   obst_riog_min_height  =  0.0_riosh
   obst_riog_max_height  =  10000.0_riosh
   obst_riog_min_dens    =  0.0_riosh
   obst_riog_max_dens    =  1000000.0_riosh
   obst_riog_min_width   =  0.0_riosh
   obst_riog_max_width   =  1000.0_riosh
   obst_riog_min_thick   =  0.0_riosh
   obst_riog_max_thick   =  1000.0_riosh
   obst_riog_min_oai     =  0.0_riosh
   obst_riog_max_oai     =  1000.0_riosh
   obst_riog_min_theta   =  0.0_riosh
   obst_riog_max_theta   =  180._riosh
   obst_riog_min_cover   =  0.0_riosh
   obst_riog_max_cover   =  1.0_riosh
   obst_riog_min_fracz   =  0.0_riosh
   obst_riog_max_fracz   =  1.0_riosh
   obst_riog_min_fuv     = -10000.0_riosh
   obst_riog_max_fuv     =  10000.0_riosh
   obst_riog_min_a       =  0.0_riosh
   obst_riog_max_a       =  1.0_riosh
   obst_riog_min_s       =  0.0_riosh
   obst_riog_max_s       =  10000.0_riosh
   obst_riog_min_drag    =  0.0_riosh
   obst_riog_max_drag    =  10.0_riosh
   obst_riog_min_tau     = -10000.0_riosh
   obst_riog_max_tau     =  10000.0_riosh
   obst_riog_min_z0      =  0.0_riosh
   obst_riog_max_z0      =  1.0_riosh
   obst_riog_min_bstress =  0.0_riosh
   obst_riog_max_bstress =  100.0_riosh
   obst_riog_min_zroot   =  0.0_riosh
   obst_riog_max_zroot   =  10.0_riosh
   ! ***********************
   ! * READING VARIABLES
   ! ***********************
   DO iv=1,obst_nbvar
     CALL OBSTRUCTIONS_readvar(iv)
   ENDDO
   ! **********************
   ! * TABLES ALLOCATION
   ! ***********************
   CALL OBSTRUCTIONS_alloc_xyz
   CALL OBSTRUCTIONS_alloc_other
   ! **********************
   ! * CHECKS AND SUMMARY
   ! ***********************
   CALL OBSTRUCTIONS_write_summary
   CALL OBSTRUCTIONS_compatibility
   ! ***********************
   ! * READING POSITION FILE
   ! ***********************
   CALL OBSTRUCTIONS_readfile_pos
   ! **********************
   ! * READING DISTRIBUTION FILE
   ! ***********************
   CALL OBSTRUCTIONS_readfile_distri
   ! ***************************************
   ! * INITIALIZATION OF OBT_SIG / OBST_DSIG
   ! ***************************************
   ! Dynamic allocation of sw
   ALLOCATE(sw(0:obst_kmax+1))
   ! Allocation of virtual sigma for 2D
   ! if 3D, obst_sig and obst_dsig compute at each time step !TODO
   
   ! Dynamic deallocation
   DEALLOCATE(sw)
   ! ***********************************
   ! * INITIALIZATION OF ROUGHNESS SEDIM
   ! ***********************************
   obst_l_z0bstress_tot = .FALSE.
   DO iv=1,obst_nbvar
     IF(obst_l_z0bstress(iv)) THEN
       obst_l_z0bstress_tot = .TRUE. ! Only one variable used z0sed
     ENDIF
   ENDDO
   IF(.NOT.obst_l_z0bstress_tot) THEN
     obst_z0bstress(:,:) = obst_i_z0bstress
   ENDIF
   obst_fws2=obst_fricwav*0.5_rsh
   ! **********************
   ! * OTHER INITIALIZATIONS
   ! ***********************
   obst_z0bed(:,:)       = z0b(:,:) ! Saving the surface z0 tables from hydraulical initialization (without obstructions)
   obst_c_exp3d = 1.0_rsh-obst_c_imp3d
   obst_c_exp2d = 1.0_rsh-obst_c_imp2d
   ! ***********************

   END SUBROUTINE OBSTRUCTIONS_init

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readvar(iv)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readvar  ***
   !&E
   !&E ** Purpose : read "obstruction_variables.dat" file describing obstructions
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
   !&E       ! 2014-08    (F. Ganthy) Routine simplification, removing useless
   !&E                                initializations
   !&E       ! 2014-10    (F. Ganthy) More modifications + computation of obstructions
   !&E                                posture following Abdelrhman 2007
   !&E       ! 2016-08    (F. Ganthy) Optimization on Abdelrhman 2007 method :
   !&E                                - added number of segment within namelist for each obstruction variable
   !&E                                - added sheltering coefficient witinh namelist
   !&E                                - update number according with these changes
   !&E       ! 2017-02    (F. Ganthy) Some modifications:
   !&E                                - Differenciation of cylindric / parallelepipedic structures
   !&E                                - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell) 
   !&E                                - Cleaning (removing useless parameters)
   !&E       ! 2017-10    (F. Ganthy) Add count for different types of variables
   !&E       ! 2021-10    (F. Ganthy) Convert into namelist reading
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Arguments
   INTEGER,INTENT(IN) :: iv
   !! * Local declaration
   CHARACTER(len=lchain) :: filepc
   ! For obst_var_main
   INTEGER               :: r_obst_varnum
   CHARACTER(len=lchain) :: r_obst_varname
   ! For obst_var_option
   LOGICAL               :: r_obst_l_cylindre,r_obst_l_flexible,r_obst_l_downward, &
                            r_obst_l_3dobst,r_obst_l_noturb
   ! For obst_var_init
   LOGICAL               :: r_obst_l_filechar,r_obst_l_init_spatial,r_obst_l_filedistri
   CHARACTER(len=lchain) :: r_obst_fn_char,r_obst_fn_distrib
   REAL(KIND=rsh)        :: r_obst_i_height,r_obst_i_width,r_obst_i_thick,r_obst_i_dens
   ! For obst_var_flexibility
   LOGICAL               :: r_obst_l_abdelposture,r_obst_l_param_height
   INTEGER               :: r_obst_c_abdel_nmax
   REAL(KIND=rsh)        :: r_obst_c_rho,r_obst_c_lift,r_obst_c_shelter,r_obst_c_height_x0,r_obst_c_height_x1
   ! For obst_var_roughdrag
   LOGICAL               :: r_obst_l_drag_cste,r_obst_l_abdelrough_cste
   REAL(KIND=rsh)        :: r_obst_c_crough_x0,r_obst_c_crough_x1,r_obst_c_drag,r_obst_c_lz
   ! For obst_var_fracxy
   LOGICAL               :: r_obst_l_fracxy
   INTEGER               :: r_obst_fracxy_type
   REAL(KIND=rsh)        :: r_obst_c_fracxy_k0,r_obst_c_fracxy_k1,r_obst_c_fracxy_l
   ! For obst_var_bstress
   LOGICAL               :: r_obst_l_z0bstress
   INTEGER               :: r_obst_z0bstress_option
   REAL(KIND=rsh)        :: r_obst_c_z0bstress,r_obst_c_z0bstress_x0,r_obst_c_z0bstress_x1,r_obst_c_z0bstress_x2
   !! * Namelists
   NAMELIST/obst_var_main/r_obst_varnum,r_obst_varname
   NAMELIST/obst_var_option/r_obst_l_cylindre,r_obst_l_flexible,r_obst_l_downward,          &
                            r_obst_l_3dobst,r_obst_l_noturb
   NAMELIST/obst_var_init/r_obst_l_filechar,r_obst_l_init_spatial,r_obst_l_filedistri,      &
                          r_obst_fn_char,r_obst_fn_distrib,                                 &
                          r_obst_i_height,r_obst_i_width,r_obst_i_thick,r_obst_i_dens
   NAMELIST/obst_var_flexibility/r_obst_l_abdelposture,r_obst_l_param_height,               &
                                 r_obst_c_abdel_nmax,r_obst_c_rho,r_obst_c_lift,            &
                                 r_obst_c_shelter,r_obst_c_height_x0,r_obst_c_height_x1
   NAMELIST/obst_var_roughdrag/r_obst_l_drag_cste,r_obst_l_abdelrough_cste,                 &
                               r_obst_c_crough_x0,r_obst_c_crough_x1,r_obst_c_drag,r_obst_c_lz
   NAMELIST/obst_var_fracxy/r_obst_l_fracxy,r_obst_fracxy_type,r_obst_c_fracxy_k0,          &
                            r_obst_c_fracxy_k1,r_obst_c_fracxy_l
   NAMELIST/obst_var_bstress/r_obst_l_z0bstress,r_obst_z0bstress_option,r_obst_c_z0bstress, &
                             r_obst_c_z0bstress_x0,r_obst_c_z0bstress_x1,r_obst_c_z0bstress_x2
   !!----------------------------------------------------------------------
   !! * Executable part

   ! save into simu.log
   !-------------------
    MPI_master_only  WRITE(iscreenlog,*) ' '
    MPI_master_only  WRITE(iscreenlog,*) ' '
    MPI_master_only  WRITE(iscreenlog,*) '****************************************************************'
    MPI_master_only  WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READVAR *****'
    MPI_master_only  WRITE(iscreenlog,*) '****************************************************************'
    MPI_master_only  WRITE(iscreenlog,*) ' Reading file ',TRIM(obst_fn_vardat(iv))
    MPI_master_only  WRITE(iscreenlog,*) ' defining obstructions parameters'
    MPI_master_only  WRITE(iscreenlog,*) '********************************************************'
   ! ******************************
   ! * Start reading variables file
   ! ******************************
   filepc = obst_fn_vardat(iv)
   OPEN(55,file=filepc,status='old',form='formatted',access='sequential')
   READ(55,obst_var_main)
   READ(55,obst_var_option)
   READ(55,obst_var_init)
   READ(55,obst_var_flexibility)
   READ(55,obst_var_roughdrag)
   READ(55,obst_var_fracxy)
   READ(55,obst_var_bstress)
   CLOSE(55)
   ! ************************************
   ! * Alocate to corresponding parameter
   ! ************************************
   ! * For namelist obst_var_main
   obst_varnum(iv)            = r_obst_varnum
   obst_varname(iv)           = r_obst_varname
   ! * For namelist obst_var_option
   obst_l_cylindre(iv)        = r_obst_l_cylindre
   obst_l_flexible(iv)        = r_obst_l_flexible
   obst_l_downward(iv)        = r_obst_l_downward
   obst_l_3dobst(iv)          = r_obst_l_3dobst
   obst_l_noturb(iv)          = r_obst_l_noturb
   ! * For namelist obst_var_init
   obst_l_filechar(iv)        = r_obst_l_filechar
   obst_l_init_spatial(iv)    = r_obst_l_init_spatial
   obst_l_filedistri(iv)      = r_obst_l_filedistri
   obst_fn_char(iv)           = r_obst_fn_char
   obst_fn_distrib(iv)        = r_obst_fn_distrib
   obst_i_height(iv)          = r_obst_i_height
   obst_i_width(iv)           = r_obst_i_width
   obst_i_thick(iv)           = r_obst_i_thick
   obst_i_dens(iv)            = r_obst_i_dens
   ! * For namelist obst_var_flexibility
   obst_l_abdelposture(iv)    = r_obst_l_abdelposture
   obst_l_param_height(iv)    = r_obst_l_param_height
   obst_c_abdel_nmax(iv)      = r_obst_c_abdel_nmax
   obst_c_rho(iv)             = r_obst_c_rho
   obst_c_lift(iv)            = r_obst_c_lift
   obst_c_shelter(iv)         = r_obst_c_shelter
   obst_c_height_x0(iv)       = r_obst_c_height_x0
   obst_c_height_x1(iv)       = r_obst_c_height_x1
   ! * For namelist obst_var_roughdrag
   obst_l_drag_cste(iv)       = r_obst_l_drag_cste
   obst_l_abdelrough_cste(iv) = r_obst_l_abdelrough_cste
   obst_c_crough_x0(iv)       = r_obst_c_crough_x0
   obst_c_crough_x1(iv)       = r_obst_c_crough_x1
   obst_c_drag(iv)            = r_obst_c_drag
   obst_c_lz(iv)              = r_obst_c_lz
   ! * For namelist obst_var_fracxy
   obst_l_fracxy(iv)          = r_obst_l_fracxy
   obst_fracxy_type(iv)       = r_obst_fracxy_type
   obst_c_fracxy_k0(iv)       = r_obst_c_fracxy_k0
   obst_c_fracxy_k1(iv)       = r_obst_c_fracxy_k1
   obst_c_fracxy_l(iv)        = r_obst_c_fracxy_l
   ! * For namelist obst_var_bstress
   obst_l_z0bstress(iv)       = r_obst_l_z0bstress
   obst_z0bstress_option(iv)  = r_obst_z0bstress_option
   obst_c_z0bstress(iv)       = r_obst_c_z0bstress
   obst_c_z0bstress_x0(iv)    = r_obst_c_z0bstress_x0
   obst_c_z0bstress_x1(iv)    = r_obst_c_z0bstress_x1
   obst_c_z0bstress_x2(iv)    = r_obst_c_z0bstress_x2
   !!********************************** 
   END SUBROUTINE OBSTRUCTIONS_readvar

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_write_summary
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_write_summary  ***
   !&E
   !&E ** Purpose : Write a summary of obstructions parameters
   !&E
   !&E ** Called by : OBSTRUCTIONS_init
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code (initally within obst_read_var)
   !&E       ! 2014-02    (F. Ganthy) Mooved into independant subroutine
   !&E                                to permit it use for test-cases
   !&E       ! 2014-02    (F. Ganthy) Add test on turbulence closure scheme
   !&E       ! 2014-10    (F. Ganthy) Some modifications + computation of obstructions
   !&E                                posture following Abdelrhman 2007
   !&E       ! 2015-10    (F. Ganthy) Add test on OMP parallelization
   !&E       ! 2016-08    (F. Ganthy) Optimization on Abdelrhman 2007 method :
   !&E                                updated summary according with
   !&E       ! 2017-02    (F. Ganthy) Some modifications:
   !&E                                - Differenciation of cylindric / parallelepipedic structures
   !&E                                - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell) 
   !&E                                - Cleaning (removing useless parameters)
   !&E       ! 2017-04    (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04    (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10    (F. Ganthy) Add vairables types and test of validity
   !&E       ! 2021-10    (F. Ganthy) Convert into namelist reading
   !&E       ! 2021-10    (F. Ganthy) Extract all check into OBSTRUCTIONS_compatibility
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Local declaration
   INTEGER          :: iv
   !!----------------------------------------------------------------------
   !! * Executable part
   !***********************************************
MPI_master_only  WRITE(iscreenlog,*) ' '
MPI_master_only  WRITE(iscreenlog,*) ' '
MPI_master_only  WRITE(iscreenlog,*) '***********************************************************************'
MPI_master_only  WRITE(iscreenlog,*) '********************** module OBSTRUCTIONS ****************************'
MPI_master_only  WRITE(iscreenlog,*) '************** subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(iscreenlog,*) '***********************************************************************'
MPI_master_only  WRITE(iscreenlog,*) 'LISTING OF OBSTRUCTION VARIABLES :'
MPI_master_only  WRITE(iscreenlog,*) '------------------------------------------------'
MPI_master_only  WRITE(iscreenlog,*) 'TOTAL NUMBER OF OBSTRUCTION VARIABLES : ',obst_nbvar
MPI_master_only  WRITE(iscreenlog,*) 'Number of RIGID_UP                    : ',obst_nv_rigid_up
MPI_master_only  WRITE(iscreenlog,*) 'Number of RIGID_DO                    : ',obst_nv_rigid_do
MPI_master_only  WRITE(iscreenlog,*) 'Number of FLEXI_UP                    : ',obst_nv_flexi_up
MPI_master_only  WRITE(iscreenlog,*) 'Number of FLEXI_DO                    : ',obst_nv_flexi_do
MPI_master_only  WRITE(iscreenlog,*) 'Number of 3DVARS                      : ',obst_nv_3d
MPI_master_only  WRITE(iscreenlog,*) '------------------------------------------------'
MPI_master_only  WRITE(iscreenlog,*) 'File for obstruction position is : ',TRIM(obst_fn_position)
     DO iv = 1,obst_nbvar
MPI_master_only  WRITE(iscreenlog,*) '***********************************************************************************'
MPI_master_only  WRITE(iscreenlog,*) '!===========================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_main  !'
MPI_master_only  WRITE(iscreenlog,*) '!===========================!'
MPI_master_only  WRITE(iscreenlog,*) '  - Number (identifier) of the variable  : ',obst_varnum(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Name (identifier) of the variable    : ',obst_varname(iv)
MPI_master_only  WRITE(iscreenlog,*) '!=============================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_option  !'
MPI_master_only  WRITE(iscreenlog,*) '!=============================!'
MPI_master_only  WRITE(iscreenlog,*) '  - If the current variable is cylinder-like                : ',obst_l_cylindre(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - If the current variable is flexible                     : ',obst_l_flexible(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - If the current variable is downvard                     : ',obst_l_downward(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - If the current variable is full 3D                      : ',obst_l_3dobst(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - If the current variable is considered a macro-roughness : ',obst_l_noturb(iv)
MPI_master_only  WRITE(iscreenlog,*) '!===========================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_init  !'
MPI_master_only  WRITE(iscreenlog,*) '!===========================!'
MPI_master_only  WRITE(iscreenlog,*) '  - Use a time-series file for obstructions characteristics                : ',obst_l_filechar(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Use a spatially variable file obstructions charcateristics             : ',obst_l_init_spatial(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Use a file describing the vertical distribution of obstruction density : ',obst_l_filedistri(iv)
       IF(obst_l_filechar(iv).OR.obst_l_init_spatial(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Name of temporal or spatial file for obstructions charcateristics      : ',obst_fn_char(iv)
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Name of temporal or spatial file for obstructions charcateristics      : NOT USED'
       ENDIF
       IF(obst_l_filedistri(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Name of file for the vertical distribution of obstruction density      : ',obst_fn_distrib(iv)
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Name of file for the vertical distribution of obstruction density      : NOT USED'
       ENDIF
       IF(obst_l_filechar(iv).OR.obst_l_init_spatial(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Initial height (unbent, eg. leaf-length for segrasses) of obstructions : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Initial width (or diameter for cylindric obstructions) of obstructions : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Initial thick of obstructions (along the flow)                         : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Initial density of obstructions (or maximum density)                   : NOT USED'
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Initial height (unbent, eg. leaf-length for segrasses) of obstructions : ',obst_i_height(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Initial width (or diameter for cylindric obstructions) of obstructions : ',obst_i_width(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Initial thick of obstructions (along the flow)                         : ',obst_i_thick(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Initial density of obstructions (or maximum density)                   : ',obst_i_dens(iv)
       ENDIF
MPI_master_only  WRITE(iscreenlog,*) '!==================================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_flexibility  !'
MPI_master_only  WRITE(iscreenlog,*) '!==================================!'
       IF(obst_l_flexible(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Use Abdhelhrmans (2007) procedure to compute bending                   : ',obst_l_abdelposture(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Use empirical (exponential desrease) formulation to compute bending    : ',obst_l_param_height(iv)
         IF(obst_l_abdelposture(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : ',obst_c_abdel_nmax(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : ',obst_c_rho(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : ',obst_c_lift(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : ',obst_c_shelter(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : NOT USED'
         ELSEIF(obst_l_param_height(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : ',obst_c_height_x0(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : ',obst_c_height_x1(iv)
         ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : NOT USED'
         ENDIF
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Use Abdhelhrmans (2007) procedure to compute bending                   : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Use empirical (exponential desrease) formulation to compute bending    : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : NOT USED'
       ENDIF
MPI_master_only  WRITE(iscreenlog,*) '!================================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_roughdrag  !'
MPI_master_only  WRITE(iscreenlog,*) '!================================!'
MPI_master_only  WRITE(iscreenlog,*) '  - Use a constant drag coefficient for obstructions in hydrodynamics      : ',obst_l_drag_cste(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Use a constant drag during reconstruction of velocity profiles         : ',obst_l_abdelrough_cste(iv)
       IF(obst_l_abdelrough_cste(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - First coefficient for CD during velocity profiles reconstruction       : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second coefficient for CD during velocity profiles reconstruction      : NOT USED'
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - First coefficient for CD during velocity profiles reconstruction       : ',obst_c_crough_x0(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Second coefficient for CD during velocity profiles reconstruction      : ',obst_c_crough_x1(iv)
       ENDIF
MPI_master_only  WRITE(iscreenlog,*) '  - Drag coefficient (max value if not constant) for obstructions elements : ',obst_c_drag(iv)
       IF(obst_l_noturb(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Coef. turbulent dissipation time-scale between obstructions elements   : NOT USED'
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Coef. turbulent dissipation time-scale between obstructions elements   : ',obst_c_lz(iv)
       ENDIF
MPI_master_only  WRITE(iscreenlog,*) '!=============================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_fracxy  !'
MPI_master_only  WRITE(iscreenlog,*) '!=============================!'
MPI_master_only  WRITE(iscreenlog,*) '  - Take account for non-linear patchiness correction                      : ',obst_l_fracxy(iv)
       IF(.NOT.obst_l_fracxy(iv))THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Kind of non-linear correction method                                   : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Coefficient for the exponential correction                             : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for correction of the exponential coefficient          : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for correction of the exponential coefficient         : NOT USED'
       ELSE
         IF(obst_fracxy_type(iv).EQ.1)THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Kind of non-linear correction method                                   : Exponential'
MPI_master_only  WRITE(iscreenlog,*) '  - Coefficient for the exponential correction                             : ',obst_c_fracxy_k0(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for correction of the exponential coefficient          : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for correction of the exponential coefficient         : NOT USED'
         ELSEIF(obst_fracxy_type(iv).EQ.2)THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Kind of non-linear correction method                                   : Exponential K0-variable'
MPI_master_only  WRITE(iscreenlog,*) '  - Coefficient for the exponential correction                             : ',obst_c_fracxy_k0(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for correction of the exponential coefficient          : ',obst_c_fracxy_k1(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for correction of the exponential coefficient         : ',obst_c_fracxy_l(iv)
         ENDIF
       ENDIF
MPI_master_only  WRITE(iscreenlog,*) '!==============================!'
MPI_master_only  WRITE(iscreenlog,*) '! NAMELIST : obst_var_bstress  !'
MPI_master_only  WRITE(iscreenlog,*) '!==============================!'
MPI_master_only  WRITE(iscreenlog,*) '  - To activate the impact of obstruction on Z0 (BSS computation)          : ',obst_l_z0bstress(iv)
       IF(obst_l_z0bstress(iv))THEN
         IF(obst_z0bstress_option(iv).EQ.0)THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Option to compute the obstruction induced roughness length             : Constant'
MPI_master_only  WRITE(iscreenlog,*) '  - Constant (corrected value of roughness length)                         : ',obst_c_z0bstress(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for rouhgness length computation (in 3D)               : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for rouhgness length computation (in 3D)              : NOT USED'       
         ELSEIF(obst_z0bstress_option(iv).EQ.1)THEN
MPI_master_only  WRITE(iscreenlog,*) '  - Option to compute the obstruction induced roughness length             : Parameterized'
MPI_master_only  WRITE(iscreenlog,*) '  - Constant (corrected value of roughness length)                         : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for rouhgness length computation (in 3D)               : ',obst_c_z0bstress_x0(iv)
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for rouhgness length computation (in 3D)              : ',obst_c_z0bstress_x1(iv)
         ENDIF
MPI_master_only  WRITE(iscreenlog,*) '  - Coefficient to correct 3D roughness length into 2D roughness length    : ',obst_c_z0bstress_x2(iv)
       ELSE
MPI_master_only  WRITE(iscreenlog,*) '  - Option to compute the obstruction induced roughness length             : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Constant (corrected value of roughness length)                         : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - First parameter for rouhgness length computation (in 3D)               : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Second parameter for rouhgness length computation (in 3D)              : NOT USED'
MPI_master_only  WRITE(iscreenlog,*) '  - Coefficient to correct 3D roughness length into 2D roughness length    : NOT USED'
       ENDIF
     ENDDO
     MPI_master_only  WRITE(iscreenlog,*) '***************************************************************'
   !-------------------------------------------
   END SUBROUTINE OBSTRUCTIONS_write_summary

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_compatibility
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_compatibility ***
   !&E
   !&E ** Purpose : Check compatibility of parameterization used
   !&E
   !&E ** Called by : OBSTRUCTIONS_init
   !&E
   !&E ** History :
   !&E       ! 2021-10    (F. Ganthy) Original code (extracted from OBSTRUCTIONS_write_summary)
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Local declaration
   INTEGER          :: iv,nv_tot
   LOGICAL          :: l_turb

   !!----------------------------------------------------------------------
   !! * Executable part
   ! ********************************************
   ! COUNT NUMBER OF VARIABLES OF DIFFERENT TYPES
   ! ********************************************
   obst_nv_up       = 0
   obst_nv_do       = 0
   obst_nv_3d       = 0
   obst_nv_rigid_up = 0
   obst_nv_rigid_do = 0
   obst_nv_flexi_up = 0
   obst_nv_flexi_do = 0
   obst_nv_turb     = 0
   obst_nv_noturb   = 0
   DO iv=1,obst_nbvar
     IF (obst_l_noturb(iv)) THEN
         obst_nv_noturb     = obst_nv_noturb+1
     ELSE
         obst_nv_turb       = obst_nv_turb+1
     ENDIF
     IF (obst_l_3dobst(iv)) THEN
       obst_nv_3d           = obst_nv_3d+1
     ELSE
       IF(obst_l_downward(iv))THEN
         obst_nv_do         = obst_nv_do+1
         IF(obst_l_flexible(iv))THEN
           obst_nv_flexi_do = obst_nv_flexi_do+1
         ELSE
           obst_nv_rigid_do = obst_nv_rigid_do+1
         ENDIF
       ELSE
         obst_nv_up         = obst_nv_up+1
         IF(obst_l_flexible(iv))THEN
           obst_nv_flexi_up = obst_nv_flexi_up+1
         ELSE
           obst_nv_rigid_up = obst_nv_rigid_up+1
         ENDIF
       ENDIF
     ENDIF
   ENDDO
   ! *************************************
   ! TESTS AND WRITE ON SIMULOG WARNLOG...
   ! *************************************
   !--------------------------------
   ! TESTING THE NUMBER OF VARIABLES
   !--------------------------------
   nv_tot = obst_nv_rigid_up + obst_nv_rigid_do + obst_nv_flexi_up + obst_nv_flexi_do + obst_nv_3d
   IF (nv_tot /= obst_nbvar) THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : The total number of obstruction variables read from'
MPI_master_only  WRITE(ierrorlog,*) '         file',TRIM(obst_fn_vardat(iv))
MPI_master_only  WRITE(ierrorlog,*) '         is DIFFERENT from the number of variables defined'
MPI_master_only  WRITE(ierrorlog,*) '         nbvar',obst_nbvar
MPI_master_only  WRITE(ierrorlog,*) '         nv_tot',nv_tot
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
   ENDIF
   !------------------------
   ! TESTING VARIABLES ORDER
   !------------------------
   DO iv = 1,obst_nbvar
     IF(obst_varnum(iv) /= iv) THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Wrong obstruction variable order'
MPI_master_only  WRITE(ierrorlog,*) '         Variable number ', obst_varnum(iv)
MPI_master_only  WRITE(ierrorlog,*) '         corresponds to the ', iv, ' variable read'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
   ENDDO
   
   !---------------------------
   ! TEST FOR 2D and GLS_KEPSILON
   !---------------------------
#ifndef SOLVE3D
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : OBSTRUCTION only available with SOLVE3D'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
#endif
   !---------------------------------------------------
   ! TEST FOR VARIABLES USING FULL TURBULENCE PROCEDURE
   !---------------------------------------------------
#ifndef GLS_KEPSILON
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '************************* module OBSTRUCTIONS ************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : OBSTRUCTION only available with GLS_KEPSILON'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
#endif
   !--------------------------------------------------------------
   ! TEST CONSISTENCY FOR SCHEME AND DESCRIPTION FOR EACH VARIALBE
   !--------------------------------------------------------------
   DO iv=1,obst_nbvar
     IF((obst_l_downward(iv)).AND.(obst_l_noturb(iv)))THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
MPI_master_only  WRITE(ierrorlog,*) '         This variable is defined as a downward one, while'
MPI_master_only  WRITE(ierrorlog,*) '         it uses simplified (obst_l_noturb) procedure'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
     IF((obst_l_3dobst(iv)).AND.(obst_l_noturb(iv)))THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '************************* module OBSTRUCTIONS ************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
MPI_master_only  WRITE(ierrorlog,*) '         This variable is defined as a full 3d one, while'
MPI_master_only  WRITE(ierrorlog,*) '         it uses simplified (obst_l_noturb) procedure'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
     IF((obst_l_3dobst(iv)).AND.(obst_l_flexible(iv)))THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
MPI_master_only  WRITE(ierrorlog,*) '         This variable is defined as a full 3d one, while'
MPI_master_only  WRITE(ierrorlog,*) '         it is also defined as a flexible one'
MPI_master_only  WRITE(ierrorlog,*) '         This kind of variable is not available yet...'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
     IF((obst_l_3dobst(iv)).AND.(obst_l_downward(iv)))THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
MPI_master_only  WRITE(ierrorlog,*) '         This variable is defined as a full 3d one, while'
MPI_master_only  WRITE(ierrorlog,*) '         it is also defined as a downward one'
MPI_master_only  WRITE(ierrorlog,*) '         This kind of variable is not available yet...'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
     IF((obst_l_downward(iv)).AND.(obst_l_flexible(iv)).AND.(obst_l_abdelposture(iv)))THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
MPI_master_only  WRITE(ierrorlog,*) '         This variable is defined as a flexible and downward one'
MPI_master_only  WRITE(ierrorlog,*) '         and uses Abdelhrman procedure'
MPI_master_only  WRITE(ierrorlog,*) '         This procedure is not available yet for downward variables...'
MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
     IF(obst_l_init_spatial(iv).AND.obst_l_filechar(iv))THEN
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) ' '
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
MPI_master_only  WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
MPI_master_only  WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
MPI_master_only  WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
MPI_master_only  WRITE(ierrorlog,*) '         Spatial initialization is defined with also spatial forcing'
MPI_master_only  WRITE(ierrorlog,*) '         for obstruction characteristics : This is not ready yet...'
MPI_master_only  WRITE(ierrorlog,*) '         Choose between : spatial initialization (but no time-varying characteristics)'
MPI_master_only  WRITE(ierrorlog,*) '         and homogeneous initialization (but with time-varying characteristics)'
MPI_master_only  WRITE(ierrorlog,*) '**********************************************************************'
        STOP
     ENDIF
   ENDDO
   !-------------------------------------------
   END SUBROUTINE OBSTRUCTIONS_compatibility

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readfile_char(limin, limax, ljmin, ljmax)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readfile_char  ***
   !&E
   !&E ** Purpose : Readinf .dat file for time-varying obstructions characteristics
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-10-16 (F. Ganthy) Replacing diameter by width and thickness  
   !&E       ! 2015-01-15 (F. Ganthy) MPI parallelization
   !&E       ! 2016-09-20 (F. Ganthy) Change computation of obstruction area index to be more 
   !&E       !                        consistent with biological aspects
   !&E       ! 2017-02-16 (F. Ganthy) Changes on instantaneous obstruction state variables for 
   !&E                                future coupling with Zostera growth module
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax


   !! * Local declaration
   LOGICAL               :: ex
   INTEGER               :: i,j,eof,kk,iv,numfile
   CHARACTER(LEN=lchain) :: rec
   REAL(KIND=rlg)        :: tool_datosec,tint1,tint2,dt1,dt2,t1,tdb,tfi
   REAL(KIND=rsh)        :: height1,height2,width1,width2,thick1,thick2,dens1,dens2

   !!--------------------------------------------------------------------------
   !! * Executable part
   !-------------------------
   ! **** Initialization ****
   !-------------------------
        !TODO : add reading of temporal file (careful, no tchrono in croco)
   !-------------------------------------
   ! **** End of initialization part ****
   !-------------------------------------
   !
   !-------------------------
   ! *** Current progress ***
   !-------------------------
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv = 1,obst_nbvar
         IF(.NOT.obst_l_init_spatial(iv))THEN
           IF(obst_position(iv,i,j).GT.0.0_rsh) THEN
             obst_height_inst(iv,i,j)      = obst_i_height(iv) !TODO : add reading of temporal file (careful, no tchrono in croco)
             obst_width_inst(iv,i,j)       = obst_i_width(iv)  !TODO : add reading of temporal file (careful, no tchrono in croco)
             obst_thick_inst(iv,i,j)       = obst_i_thick(iv)  !TODO : add reading of temporal file (careful, no tchrono in croco)
             obst_dens_inst(iv,i,j)        = obst_i_dens(iv)   !TODO : add reading of temporal file (careful, no tchrono in croco)
             IF(obst_l_cylindre(iv))THEN ! Cylindric/Ellipse obstruction
               obst_area_index_inst(iv,i,j)  = obst_dens_inst(iv,i,j) * obst_height_inst(iv,i,j) * &
                                                (2.0_rsh*pi*SQRT(0.5_rsh*(obst_width_inst(iv,i,j)**2.0_rsh + & 
                                                obst_thick_inst(iv,i,j)**2.0_rsh)))
             ELSE ! Parallelepipedic obstruction
               obst_area_index_inst(iv,i,j)  = 2.0_rsh * obst_dens_inst(iv,i,j) * obst_width_inst(iv,i,j) * &
                                                obst_height_inst(iv,i,j)
             ENDIF
           ELSE
             obst_height_inst(iv,i,j)      = 0.0_rsh
             obst_width_inst(iv,i,j)       = 0.0_rsh
             obst_thick_inst(iv,i,j)       = 0.0_rsh
             obst_dens_inst(iv,i,j)        = 0.0_rsh
             obst_area_index_inst(iv,i,j)  = 0.0_rsh
           ENDIF
         ELSE
           IF(obst_position(iv,i,j).GT.0.0_rsh) THEN
             IF(obst_l_cylindre(iv))THEN ! Cylindric/Ellipse obstruction
               obst_area_index_inst(iv,i,j)  = obst_dens_inst(iv,i,j) * obst_height_inst(iv,i,j) * &
                                                (2.0_rsh*pi*SQRT(0.5_rsh*(obst_width_inst(iv,i,j)**2.0_rsh + &
                                                obst_thick_inst(iv,i,j)**2.0_rsh)))
             ELSE ! Parallelepipedic obstruction
               obst_area_index_inst(iv,i,j)  = 2.0_rsh * obst_dens_inst(iv,i,j) * obst_width_inst(iv,i,j) * &
                                                obst_height_inst(iv,i,j)
             ENDIF
           ELSE
             obst_area_index_inst(iv,i,j)  = 0.0_rsh
           ENDIF
         ENDIF
       ENDDO
      ENDDO
   ENDDO
   ! *********************************
   END SUBROUTINE OBSTRUCTIONS_readfile_char

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readfile_pos
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readfile_pos  ***
   !&E
   !&E ** Purpose : read the obstruction position file
   !&E              in the historical DEL/AO format
   !&E
   !&E ** Called by : obst_init
   !&E
   !&E ** Modified variables : obst_position
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   ! to do for croco with appropriate netcdf reading

  
   END SUBROUTINE OBSTRUCTIONS_readfile_pos

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readfile_distri
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readfile_distri  ***
   !&E
   !&E ** Purpose : Initialization of obstruction parameters
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Local declaration
   INTEGER :: kk,iv,nb_max_hnorm
   LOGICAL :: ex
   !!----------------------------------------------------------------------
   !! * Executable part

   MPI_master_only  WRITE(iscreenlog,*) ' '
   MPI_master_only  WRITE(iscreenlog,*) ' '
   MPI_master_only  WRITE(iscreenlog,*) '************************************************************************'
   MPI_master_only  WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_DISTRI *****'
   MPI_master_only  WRITE(iscreenlog,*) '************************************************************************'
   DO iv = 1,obst_nbvar
     IF(obst_l_filedistri(iv)) THEN
        MPI_master_only  WRITE(iscreenlog,*) ' '
        MPI_master_only  WRITE(iscreenlog,*) 'Obstruction variable : ',obst_varname(iv)
        MPI_master_only  WRITE(iscreenlog,*) 'file defining vertical distribution of obstuctions :'
        MPI_master_only  WRITE(iscreenlog,*) TRIM(obst_fn_distrib(iv))
       INQUIRE(file=obst_fn_distrib(iv),exist=ex)
       IF(ex) THEN
         OPEN(53,file=obst_fn_distrib(iv),form='formatted')
         READ(53,*) ! Filename
         READ(53,*) ! Title
         READ(53,*) obst_nbhnorm(iv)
         MPI_master_only  WRITE(iscreenlog,*) ' '
         MPI_master_only  WRITE(iscreenlog,*) 'Number of vertical discretization : ',obst_nbhnorm(iv)
         CLOSE(53)
       ELSE
    MPI_master_only  WRITE(ierrorlog,*) ' '
    MPI_master_only  WRITE(ierrorlog,*) ' '
    MPI_master_only  WRITE(ierrorlog,*) '************************************************************************'
    MPI_master_only  WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_DISTRI *****'
    MPI_master_only  WRITE(ierrorlog,*) '************************************************************************'
    MPI_master_only  WRITE(ierrorlog,*) ' ERROR : File ',TRIM(obst_fn_distrib(iv)),'does not exist'
    MPI_master_only  WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
    MPI_master_only  WRITE(ierrorlog,*) '************************************************************************'
        STOP
       ENDIF
     ELSE
       obst_nbhnorm(iv) = 11
    MPI_master_only  WRITE(iscreenlog,*) ' '
    MPI_master_only  WRITE(iscreenlog,*) 'Obstruction variable : ',obst_varname(iv)
    MPI_master_only  WRITE(iscreenlog,*) 'An homogene distribution is applied'
    MPI_master_only  WRITE(iscreenlog,*) ' '
    MPI_master_only  WRITE(iscreenlog,*) 'Number of vertical discretization : ',obst_nbhnorm(iv)
     ENDIF
   ENDDO

   nb_max_hnorm = 0
   IF (obst_nbvar.EQ.1)THEN
     nb_max_hnorm = obst_nbhnorm(1)
   ELSE
     DO iv = 1,obst_nbvar-1
       nb_max_hnorm = MAX(obst_nbhnorm(iv),obst_nbhnorm(iv+1))
     ENDDO
   ENDIF
   ALLOCATE(obst_dens_norm   (1:obst_nbvar,1:nb_max_hnorm))
   ALLOCATE(obst_height_norm (1:obst_nbvar,1:nb_max_hnorm))

   DO iv = 1,obst_nbvar
     IF(obst_l_filedistri(iv)) THEN
       OPEN(53,file=obst_fn_distrib(iv),form='formatted')
       READ(53,*) ! Filename
       READ(53,*) ! Title
       READ(53,*) !NBhnorm
       READ(53,*) ! Column titles
       kk=1
       DO WHILE (kk.LE.obst_nbhnorm(iv))
         READ(53,*) obst_height_norm(iv,kk),obst_dens_norm(iv,kk)
    MPI_master_only  WRITE(iscreenlog,*) 'kk',kk
    MPI_master_only  WRITE(iscreenlog,*) 'height : ',obst_height_norm(iv,kk), ' dens : ',obst_dens_norm(iv,kk)
         kk = kk + 1
       ENDDO
       CLOSE(53)
     ELSE
       ! Homogene distribution of obstructions distribution
       ! Homogene distribution of obstructions distribution
       obst_dens_norm(iv,:) = 100.0_rsh
       obst_height_norm(iv,:) = 0.001_rsh
       DO kk=1,obst_nbhnorm(iv)
         obst_height_norm(iv,kk) = (REAL(kk-1))*10.0_rsh
         obst_dens_norm(iv,kk) = 100.0_rsh
       ENDDO
       obst_height_norm(iv,obst_nbhnorm(iv))=100.01_rsh
     ENDIF
   ENDDO

END SUBROUTINE OBSTRUCTIONS_readfile_distri

   !!==========================================================================================================


SUBROUTINE OBSTRUCTIONS_alloc_nbvar

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE OBSTRUCTIONS_alloc_nbvar  ***
    !&E
    !&E ** Purpose : Allocation of tables depending on number of obstructions variables
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

    IMPLICIT NONE
 
    !!----------------------------------------------------------------------
    !! * Executable part
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
 
    ALLOCATE(obst_l_filechar           (1:obst_nbvar))
    obst_l_filechar(:)                 = .FALSE.
    ALLOCATE(obst_l_filedistri         (1:obst_nbvar))
    obst_l_filedistri(:)               = .FALSE.
    ALLOCATE(obst_l_init_spatial       (1:obst_nbvar))
    obst_l_init_spatial(:)             = .FALSE.
    ALLOCATE(obst_l_flexible           (1:obst_nbvar))
    obst_l_flexible(:)                 = .FALSE.
    ALLOCATE(obst_l_cylindre           (1:obst_nbvar))
    obst_l_cylindre(:)                 = .FALSE.
    ALLOCATE(obst_l_downward           (1:obst_nbvar))
    obst_l_downward(:)                 = .FALSE.
    ALLOCATE(obst_l_3dobst             (1:obst_nbvar))
    obst_l_3dobst(:)                   = .FALSE.
    ALLOCATE(obst_l_noturb             (1:obst_nbvar))
    obst_l_noturb(:)                   = .FALSE.
    ALLOCATE(obst_l_abdelrough_cste    (1:obst_nbvar))
    obst_l_abdelrough_cste             = .FALSE.
    ALLOCATE(obst_l_fracxy             (1:obst_nbvar))
    obst_l_fracxy(:)                    = .FALSE.
    ALLOCATE(obst_l_abdelposture       (1:obst_nbvar))
    obst_l_abdelposture(:)             = .FALSE.
    ALLOCATE(obst_l_param_height       (1:obst_nbvar))
    obst_l_param_height(:)             = .FALSE.
    ALLOCATE(obst_l_drag_cste          (1:obst_nbvar))
    obst_l_drag_cste(:)                = .TRUE.
    ALLOCATE(obst_l_z0bstress          (1:obst_nbvar))
    obst_l_z0bstress(:)                = .FALSE.
 
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
 END SUBROUTINE OBSTRUCTIONS_alloc_nbvar
 
    !!==========================================================================================================
 
 SUBROUTINE OBSTRUCTIONS_alloc_xyz
 
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE OBSTRUCTIONS_alloc_xyz  ***
    !&E
    !&E ** Purpose : Allocation of spatial obstruction tables
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
 
    IMPLICIT NONE
 
    !! * Local declaration
 
    !!----------------------------------------------------------------------
    !! * Executable part
    !------------------------------------
    ! Definition of effective kmax to use
    !------------------------------------
    obst_kmax = kmax
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
    ALLOCATE(obst_dens_inst       (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_dens_inst(:,:,:)         = 0.0_rsh
    ALLOCATE(obst_width_inst      (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_width_inst(:,:,:)        = 0.0_rsh
    ALLOCATE(obst_thick_inst      (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_thick_inst(:,:,:)        = 0.0_rsh
    ALLOCATE(obst_height_inst     (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_height_inst(:,:,:)       = 0.0_rsh
    ALLOCATE(obst_area_index_inst (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_area_index_inst(:,:,:)   = 0.0_rsh
 
    ALLOCATE(obst_position        (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_position(:,:,:)          = 0.0_rsh
    ALLOCATE(obst_height          (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_height(:,:,:)            = 0.0_rsh
    ALLOCATE(obst_oai             (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_oai(:,:,:)               = 0.0_rsh
    ALLOCATE(obst_fracxy          (1:obst_nbvar,imin:imax, jmin:jmax))
    obst_fracxy(:,:,:)            = 0.0_rsh
 
    ALLOCATE(obst_a2d             (1:obst_nbvar+3,imin:imax, jmin:jmax))
    obst_a2d(:,:,:)               = 0.0_rsh
    ALLOCATE(obst_s2d             (1:obst_nbvar+3,imin:imax, jmin:jmax))
    obst_s2d(:,:,:)               = 0.0_rsh
    ALLOCATE(obst_z0obst          (1:obst_nbvar+3,imin:imax, jmin:jmax))
    obst_z0obst(:,:,:)            = 0.0_rsh
 
    !-------------------
    ! Variables on (i,j)
    !-------------------
    ALLOCATE(obst_roswat_bot      (imin:imax, jmin:jmax))
    obst_roswat_bot               = 0.0_rsh
    ALLOCATE(obst_fu_i            (imin:imax, jmin:jmax))
    obst_fu_i(:,:)                = 0.0_rsh
    ALLOCATE(obst_fv_i            (imin:imax, jmin:jmax))
    obst_fv_i(:,:)                = 0.0_rsh
    ALLOCATE(obst_fu_e            (imin:imax, jmin:jmax))
    obst_fu_e(:,:)                = 0.0_rsh
    ALLOCATE(obst_fv_e            (imin:imax, jmin:jmax))
    obst_fv_e(:,:)                = 0.0_rsh
    ALLOCATE(obst_z0bed           (imin:imax, jmin:jmax))
    obst_z0bed(:,:)               = 0.0_rsh
    ALLOCATE(obst_bstress         (imin:imax, jmin:jmax))
    obst_bstress(:,:)             = 0.0_rsh
    ALLOCATE(obst_bstressc        (imin:imax, jmin:jmax))
    obst_bstressc(:,:)            = 0.0_rsh
    ALLOCATE(obst_bstressw        (imin:imax, jmin:jmax))
    obst_bstressw(:,:)            = 0.0_rsh
    ALLOCATE(obst_z0bstress       (imin:imax, jmin:jmax))
    obst_z0bstress(:,:)           = 0.0_rsh
    ALLOCATE(obst_raphbx          (imin:imax, jmin:jmax))
    obst_raphbx(:,:)              = 0.0_rsh
    ALLOCATE(obst_raphby          (imin:imax, jmin:jmax))
    obst_raphby(:,:)              = 0.0_rsh
     ALLOCATE(obst_frofonx        (imin:imax, jmin:jmax))
    obst_frofonx(:,:)             = 0.0_rsh
    ALLOCATE(obst_frofony         (imin:imax, jmin:jmax))
    obst_frofony(:,:)             = 0.0_rsh
    ALLOCATE(obst_dens_mean       (imin:imax, jmin:jmax))
    obst_dens_mean(:,:)           = 0.0_rsh
    ALLOCATE(obst_width_mean      (imin:imax, jmin:jmax))
    obst_width_mean(:,:)          = 0.0_rsh
    ALLOCATE(obst_height_mean     (imin:imax, jmin:jmax))
    obst_height_mean(:,:)         = 0.0_rsh
    !------------------------
    ! Variables on (iv,k,i,j)
    !------------------------
    ALLOCATE(obst_dens3d          (1:obst_nbvar,1:obst_kmax,imin:imax, jmin:jmax))
    obst_dens3d(:,:,:,:)          = 0.0_rsh
    ALLOCATE(obst_width3d         (1:obst_nbvar,1:obst_kmax,imin:imax, jmin:jmax))
    obst_width3d(:,:,:,:)         = 0.0_rsh
    ALLOCATE(obst_thick3d         (1:obst_nbvar,1:obst_kmax,imin:imax, jmin:jmax))
    obst_thick3d(:,:,:,:)         = 0.0_rsh
    ALLOCATE(obst_theta3d         (1:obst_nbvar,1:obst_kmax,imin:imax, jmin:jmax))
    obst_theta3d(:,:,:,:)         = 0.0_rsh
    ALLOCATE(obst_fracz3d         (1:obst_nbvar,1:obst_kmax,imin:imax, jmin:jmax))
    obst_fracz3d(:,:,:,:)         = 0.0_rsh
    ALLOCATE(obst_drag3d          (1:obst_nbvar,1:obst_kmax,imin:imax, jmin:jmax))
    obst_drag3d(:,:,:,:)          = 0.0_rsh
 
    ALLOCATE(obst_a3d             (1:obst_nbvar+3,1:obst_kmax,imin:imax, jmin:jmax))
    obst_a3d(:,:,:,:)             = 0.0_rsh
    ALLOCATE(obst_s3d             (1:obst_nbvar+3,1:obst_kmax,imin:imax, jmin:jmax))
    obst_s3d(:,:,:,:)             = 0.0_rsh
    !---------------------
    ! Variables on (k,i,j)
    !---------------------
    ALLOCATE(obst_zc              (1:obst_kmax,imin:imax, jmin:jmax))
    obst_zc(:,:,:)                = 0.0_rsh
    ALLOCATE(obst_dz              (1:obst_kmax,imin:imax, jmin:jmax))
    obst_dz(:,:,:)                = 0.0_rsh
    ALLOCATE(obst_uz              (1:obst_kmax,imin:imax, jmin:jmax))
    obst_uz(:,:,:)                = 0.0_rsh
    ALLOCATE(obst_vz              (1:obst_kmax,imin:imax, jmin:jmax))
    obst_vz(:,:,:)                = 0.0_rsh
    ALLOCATE(obst_fuz_i           (1:obst_kmax,imin:imax, jmin:jmax))
    obst_fuz_i(:,:,:)             = 0.0_rsh
    ALLOCATE(obst_fvz_i           (1:obst_kmax,imin:imax, jmin:jmax))
    obst_fvz_i(:,:,:)             = 0.0_rsh
    ALLOCATE(obst_fuz_e           (1:obst_kmax,imin:imax, jmin:jmax))
    obst_fuz_e(:,:,:)             = 0.0_rsh
    ALLOCATE(obst_fvz_e           (1:obst_kmax,imin:imax, jmin:jmax))
    obst_fvz_e(:,:,:)             = 0.0_rsh
    ALLOCATE(obst_t               (1:obst_kmax,imin:imax, jmin:jmax))
    obst_t(:,:,:)                 = 0.0_rsh
    ALLOCATE(obst_tau             (1:obst_kmax,imin:imax, jmin:jmax))
    obst_tau(:,:,:)               = 0.0_rsh
    !------------------------------
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
 
    IMPLICIT NONE
    !!----------------------------------------------------------------------
    !! * Executable part
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
 END SUBROUTINE OBSTRUCTIONS_alloc_other
 
    !!==========================================================================================================
 
 SUBROUTINE OBSTRUCTIONS_dealloc
 
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE OBSTRUCTIONS_dealloc  ***
    !&E
    !&E ** Purpose : Deallocation of variables
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
    IMPLICIT NONE
 
    !! * Local declaration
 
    !!----------------------------------------------------------------------
    !! * Executable part
    !-------------------------
    ! Variables on (iv)
    !--------------------
    DEALLOCATE(obst_varnum)
    DEALLOCATE(obst_nbhnorm)
    DEALLOCATE(obst_fracxy_type)
    DEALLOCATE(obst_c_abdel_nmax)
 
    DEALLOCATE(obst_l_filechar)
    DEALLOCATE(obst_l_filedistri)
    DEALLOCATE(obst_l_init_spatial)
    DEALLOCATE(obst_l_flexible)
    DEALLOCATE(obst_l_cylindre)
    DEALLOCATE(obst_l_downward)
    DEALLOCATE(obst_l_3dobst)
    DEALLOCATE(obst_l_noturb)
    DEALLOCATE(obst_l_abdelrough_cste)
    DEALLOCATE(obst_l_fracxy)
    DEALLOCATE(obst_l_abdelposture)
    DEALLOCATE(obst_l_param_height)
    DEALLOCATE(obst_l_drag_cste)
    DEALLOCATE(obst_l_z0bstress)
 
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
 END SUBROUTINE OBSTRUCTIONS_dealloc
 
!!==========================================================================================================
 


#endif
END MODULE initOBSTRUCTIONS
