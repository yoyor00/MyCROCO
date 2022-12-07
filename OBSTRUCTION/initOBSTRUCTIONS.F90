MODULE initOBSTRUCTIONS

#ifdef key_OBSTRUCTIONS
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
   !&E ** History :
   !&E     ! 2021-10 (F. Ganthy) Original code extracted from OBSTRUCTION.F90
   !&E
   !&E==========================================================================

#include "toolcpp.h"
   !! * Modules used
   USE comOBSTRUCTIONS
   USE parameters,    ONLY : lchain,rsh,rlg,imin,imax,jmin,jmax,kmax,       &
                             limin,limax,ljmin,ljmax,riosh,liminm1,ljminm1, &
                             liminp1,limaxm1,ljminp1,ljmaxm1

   IMPLICIT NONE

   !! * Accessibility

   ! function & routines of this module, called outside :
   ! PUBLIC functions
   PUBLIC OBSTRUCTIONS_init
   PUBLIC OBSTRUCTIONS_readvar
   PUBLIC OBSTRUCTIONS_write_summary
   PUBLIC OBSTRUCTIONS_compatibility
   PUBLIC OBSTRUCTIONS_readfile_char
   PUBLIC OBSTRUCTIONS_readfile_pos
   PUBLIC OBSTRUCTIONS_readfile_distri

   PRIVATE

   !! * Shared or public module variables

   !! * Private variables

  CONTAINS


   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_init
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_init  ***
   !&E
   !&E ** Purpose : Initialization of obstruction parameters
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
   !&E       ! 2014-10    (F. Ganthy) Add bending angle
   !&E       ! 2015-09    (F. Ganthy) Minor correction (obst_cmu)
   !&E       ! 2016-03    (F. Ganthy) Add logical to choose turbulence coefficient
   !&E                                values (default or user-defined)
   !&E       ! 2017-04    (F. Ganthy) Change initialization of turbulence coefficients cmu and c2turb
   !&E       ! 2017-11    (F. Ganthy) Remove useless variables related to turbulence coefficients
   !&E       ! 2021-10    (F. Ganthy) Few modification/organization on paraOBSTRUCTIONS
   !&E       ! 2022-01    (A. Le Pevedic) Added wave friction factor (used for shear stress computation)
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters,  ONLY  : rsh,riosh
   USE comvars2d,   ONLY  : iscreenlog,ierrorlog,z0b,l_modele2d
   USE comvars3d,   ONLY  : sig,dsigu
   USE_MPI toolmpi, ONLY : MASTER

   IMPLICIT NONE

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
   NAMELIST/obst_numerics/obst_i_z0bstress,obst_c_paramhuv,obst_c_imp3d,obst_c_imp2d,fricwav
   NAMELIST/obst_output/l_obstout_pos,l_obstout_height_f,l_obstout_height_e,l_obstout_dens_f,     &
                        l_obstout_dens_e,l_obstout_width_f,l_obstout_width_e,                     &
                        l_obstout_thick_f,l_obstout_thick_e,l_obstout_oai,                        &
                        l_obstout_theta,l_obstout_cover,l_obstout_frac_z,                         &
                        l_obstout_fuv,l_obstout_fuzvz,l_obstout_a2d,l_obstout_a3d,l_obstout_s2d,  &
                        l_obstout_s3d,l_obstout_drag,l_obstout_tau,                               &
                        l_obstout_z0bed,l_obstout_z0obst,l_obstout_z0bstress,l_obstout_bstress,   &
                        l_obstout_bstressc,l_obstout_bstressw
   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBST_INIT'
   ! ************************
   ! * READING NAMELIST
   ! ************************
   filepc = './paraOBSTRUCTIONS.txt'
   OPEN(50,file=filepc,status='old',form='formatted',access='sequential')
   READ(50,obst_main)
   READ(50,obst_numerics)
   READ(50,obst_output)
   CLOSE(50)
   IF_MPI(MASTER) THEN
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) '*****************************************************'
     WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBST_INIT *****'
     WRITE(iscreenlog,*) '*****************************************************'
     WRITE(iscreenlog,*) ' Reading file ',TRIM(filepc)
     WRITE(iscreenlog,*) '***********************'
   ENDIF_MPI
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
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '*****************************************************************'
         WRITE(ierrorlog,*) '******* module OBSTRUCTIONS, subroutine OBSTRUCTIONS_INIT *******'
         WRITE(ierrorlog,*) '*****************************************************************'
         WRITE(ierrorlog,*) '!!! ERROR : Maximum number of variables >6                    !!!'
         WRITE(ierrorlog,*) '!!! Few changes must be performed in the code                 !!!'
         WRITE(ierrorlog,*) '!!! to increase number of allowed variables :                 !!!'
         WRITE(ierrorlog,*) '!!! --> Namelist contains (subroutine OBSTRUCTIONS_INIT)      !!!'
         WRITE(ierrorlog,*) '!!! --> Allocation of names of parameters files               !!!'
         WRITE(ierrorlog,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       ENDIF_MPI
     ENDIF
   ENDDO
   ! ***********************
   ! * OUTPUTS VARIABLES
   ! ***********************
#ifdef key_sedim_MUSTANG
   l_obstout_bstress      = .FALSE.
   l_obstout_bstressc     = .FALSE.
   l_obstout_bstressw     = .FALSE.
#endif
   name_out_pos           = 'pos'       ! Name obstruction position
   name_out_height_f      = 'height_f'  ! Name 2D obstruction height (forcing) (iv,i,j)
   name_out_height_e      = 'height_e'  ! Name 2D obstruction height (effective) (iv,i,j)
   name_out_dens_f        = 'dens_f'    ! Name 2D obstruction density (forcing) (iv,i,j)
   name_out_dens_e        = 'dens_e'    ! Name 3D obstruction density (3D effective) (iv,k,i,j)
   name_out_width_f       = 'width_f'   ! Name 2D obstruction width (forcing) (iv,i,j)
   name_out_width_e       = 'width_e'   ! Name 3D obstruction width (3D effective) (iv,k,i,j)
   name_out_thick_f       = 'thick_f'   ! Name 2D obstruction thick (forcing) (iv,i,j)
   name_out_thick_e       = 'thick_e'   ! Name 3D obstruction thick (3D effective) (iv,k,i,j)
   name_out_oai           = 'oai'       ! Name 2D obstruction area index (iv,i,j)
   name_out_theta         = 'theta'     ! Name 3D obstruction bending angle (iv,k,i,j)
   name_out_cover         = 'frac_xy'   ! Name 2D obstruction coverage (iv,i,j)
   name_out_frac_z        = 'frac_z'    ! Name 2D obstruction fraction of sigma layer (iv,i,j)

   name_out_fuv           = 'fuv'       ! Name 2D obstruction friction force (i,j)
   name_out_fuzvz         = 'fuzvz'     ! Name 3D obstruction friction force (k,i,j)
   name_out_a3d           = 'a2d'       ! Name 2D obstruction horizontal area (iv+3,i,j)
   name_out_a3d           = 'a3d'       ! Name 3D obstruction horizontal area (iv+3,k,i,j)
   name_out_s3d           = 's2d'       ! Name 2D obstruction vertical area (iv+3,i,j)
   name_out_s3d           = 's3d'       ! Name 3D obstruction vertical area (iv+3,k,i,j)
   name_out_drag          = 'cd3d'      ! Name 3D obstruction drag coefficient (iv,k,i,j)
   name_out_tau           = 'tau3d'     ! Name 3D obstruction turbulence dissipation scale (k,i,j)
   name_out_z0bed         = 'z0bed'     ! Name 2D bottom roughness length of bed (i,j)
   name_out_z0obst        = 'z0obst'    ! Name 2D bottom roughness length of obstructions (i,j)
   name_out_z0bstress     = 'z0bstress' ! Name 2D bottom roughness length used for bottom shear stress computation (i,j)
   name_out_bstress       = 'taub'      ! Name 2D total bottom shear stress
   name_out_bstressc      = 'taubc'     ! Name 2D current bottom shear stress
   name_out_bstressw      = 'taubw'     ! Name 2D wave bottom shear stress

   riog_valid_min_pos     =  0.0_riosh
   riog_valid_max_pos     =  1.0_riosh
   riog_valid_min_height  =  0.0_riosh
   riog_valid_max_height  =  10000.0_riosh
   riog_valid_min_dens    =  0.0_riosh
   riog_valid_max_dens    =  1000000.0_riosh
   riog_valid_min_width   =  0.0_riosh
   riog_valid_max_width   =  1000.0_riosh
   riog_valid_min_thick   =  0.0_riosh
   riog_valid_max_thick   =  1000.0_riosh
   riog_valid_min_oai     =  0.0_riosh
   riog_valid_max_oai     =  1000.0_riosh
   riog_valid_min_theta   =  0.0_riosh
   riog_valid_max_theta   =  180._riosh
   riog_valid_min_cover   =  0.0_riosh
   riog_valid_max_cover   =  1.0_riosh
   riog_valid_min_fracz   =  0.0_riosh
   riog_valid_max_fracz   =  1.0_riosh
   riog_valid_min_fuv     = -10000.0_riosh
   riog_valid_max_fuv     =  10000.0_riosh
   riog_valid_min_a       =  0.0_riosh
   riog_valid_max_a       =  1.0_riosh
   riog_valid_min_s       =  0.0_riosh
   riog_valid_max_s       =  10000.0_riosh
   riog_valid_min_drag    =  0.0_riosh
   riog_valid_max_drag    =  10.0_riosh
   riog_valid_min_tau     = -10000.0_riosh
   riog_valid_max_tau     =  10000.0_riosh
   riog_valid_min_z0      =  0.0_riosh
   riog_valid_max_z0      =  1.0_riosh
   riog_valid_min_bstress =  0.0_riosh
   riog_valid_max_bstress =  100.0_riosh
   riog_valid_min_zroot   =  0.0_riosh
   riog_valid_max_zroot   =  10.0_riosh
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
   ! Allocation of virtual or reel sigma
   IF((l_modele2d).OR.(kmax.EQ.1))THEN
     DO k=1,obst_kmax
       obst_sig(k) = (REAL(k-obst_kmax,rsh)-0.5_rsh)/REAL(obst_kmax,rsh)
     ENDDO
     obst_sig(obst_kmax+1)=-obst_sig(obst_kmax)
     obst_sig(0)=ABS(obst_sig(1))-2.0_rsh
     DO k=1,obst_kmax
       sw(k)=0.5_rsh*(obst_sig(k)+obst_sig(k+1))
     ENDDO
     sw(obst_kmax+1)=-sw(obst_kmax)
     sw(0)=-1.0_rsh
     DO k=1,obst_kmax
       obst_dsig(k) = sw(k)-sw(k-1)
     ENDDO
   ELSE
     DO k=1,obst_kmax
       obst_sig(k)  = sig(k)
       obst_dsig(k) = dsigu(k)
     ENDDO
   ENDIF
   ! Dynamic deallocation
   DEALLOCATE(sw)
   ! ***********************************
   ! * INITIALIZATION OF ROUGHNESS SEDIM
   ! ***********************************
   l_obst_z0bstress_tot = .FALSE.
   DO iv=1,obst_nbvar
     IF(l_obst_z0bstress(iv)) THEN
       l_obst_z0bstress_tot = .TRUE. ! Only one variable used z0sed
     ENDIF
   ENDDO
   IF(.NOT.l_obst_z0bstress_tot) THEN
     obst_z0bstress(:,:) = obst_i_z0bstress
   ENDIF
   fws2=fricwav*0.5_rsh
   ! **********************
   ! * OTHER INITIALIZATIONS
   ! ***********************
   obst_z0bed(:,:)       = z0b(:,:) ! Saving the surface z0 tables from hydraulical initialization (without obstructions)
   obst_c_exp3d = 1.0_rsh-obst_c_imp3d
   obst_c_exp2d = 1.0_rsh-obst_c_imp2d
   ! ***********************
   PRINT_DBG*, 'END OBSTRUCTIONS_INIT'

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
   !! * Modules used
   USE comvars2d, ONLY  : iscreenlog,iwarnlog,ierrorlog
   USE_MPI toolmpi,ONLY : MASTER

   IMPLICIT NONE

   !! * Arguments
   INTEGER,INTENT(IN) :: iv
   !! * Local declaration
   CHARACTER(len=lchain) :: filepc
   ! For obst_var_main
   INTEGER               :: r_obst_varnum
   CHARACTER(len=lchain) :: r_obst_varname
   ! For obst_var_option
   LOGICAL               :: r_l_obst_cylindre,r_l_obst_flexible,r_l_obst_downward, &
                            r_l_obst_3dobst,r_l_obst_noturb
   ! For obst_var_init
   LOGICAL               :: r_l_obst_filechar,r_l_obst_init_spatial,r_l_obst_filedistri
   CHARACTER(len=lchain) :: r_obst_fn_char,r_obst_fn_distrib
   REAL(KIND=rsh)        :: r_obst_i_height,r_obst_i_width,r_obst_i_thick,r_obst_i_dens
   ! For obst_var_flexibility
   LOGICAL               :: r_l_obst_abdelposture,r_l_obst_param_height
   INTEGER               :: r_obst_c_abdel_nmax
   REAL(KIND=rsh)        :: r_obst_c_rho,r_obst_c_lift,r_obst_c_shelter,r_obst_c_height_x0,r_obst_c_height_x1
   ! For obst_var_roughdrag
   LOGICAL               :: r_l_obst_drag_cste,r_l_obst_abdelrough_cste
   REAL(KIND=rsh)        :: r_obst_c_crough_x0,r_obst_c_crough_x1,r_obst_c_drag,r_obst_c_lz
   ! For obst_var_fracxy
   LOGICAL               :: r_l_obst_fracxy
   INTEGER               :: r_obst_fracxy_type
   REAL(KIND=rsh)        :: r_obst_c_fracxy_k0,r_obst_c_fracxy_k1,r_obst_c_fracxy_l
   ! For obst_var_bstress
   LOGICAL               :: r_l_obst_z0bstress
   INTEGER               :: r_obst_z0bstress_option
   REAL(KIND=rsh)        :: r_obst_c_z0bstress,r_obst_c_z0bstress_x0,r_obst_c_z0bstress_x1,r_obst_c_z0bstress_x2
   !! * Namelists
   NAMELIST/obst_var_main/r_obst_varnum,r_obst_varname
   NAMELIST/obst_var_option/r_l_obst_cylindre,r_l_obst_flexible,r_l_obst_downward,          &
                            r_l_obst_3dobst,r_l_obst_noturb
   NAMELIST/obst_var_init/r_l_obst_filechar,r_l_obst_init_spatial,r_l_obst_filedistri,      &
                          r_obst_fn_char,r_obst_fn_distrib,                                 &
                          r_obst_i_height,r_obst_i_width,r_obst_i_thick,r_obst_i_dens
   NAMELIST/obst_var_flexibility/r_l_obst_abdelposture,r_l_obst_param_height,               &
                                 r_obst_c_abdel_nmax,r_obst_c_rho,r_obst_c_lift,            &
                                 r_obst_c_shelter,r_obst_c_height_x0,r_obst_c_height_x1
   NAMELIST/obst_var_roughdrag/r_l_obst_drag_cste,r_l_obst_abdelrough_cste,                 &
                               r_obst_c_crough_x0,r_obst_c_crough_x1,r_obst_c_drag,r_obst_c_lz
   NAMELIST/obst_var_fracxy/r_l_obst_fracxy,r_obst_fracxy_type,r_obst_c_fracxy_k0,          &
                            r_obst_c_fracxy_k1,r_obst_c_fracxy_l
   NAMELIST/obst_var_bstress/r_l_obst_z0bstress,r_obst_z0bstress_option,r_obst_c_z0bstress, &
                             r_obst_c_z0bstress_x0,r_obst_c_z0bstress_x1,r_obst_c_z0bstress_x2
   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_READVAR'

   ! save into simu.log
   !-------------------
   IF_MPI(MASTER) THEN
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) '****************************************************************'
     WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READVAR *****'
     WRITE(iscreenlog,*) '****************************************************************'
     WRITE(iscreenlog,*) ' Reading file ',TRIM(obst_fn_vardat(iv))
     WRITE(iscreenlog,*) ' defining obstructions parameters'
     WRITE(iscreenlog,*) '********************************************************'
   ENDIF_MPI
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
   l_obst_cylindre(iv)        = r_l_obst_cylindre
   l_obst_flexible(iv)        = r_l_obst_flexible
   l_obst_downward(iv)        = r_l_obst_downward
   l_obst_3dobst(iv)          = r_l_obst_3dobst
   l_obst_noturb(iv)          = r_l_obst_noturb
   ! * For namelist obst_var_init
   l_obst_filechar(iv)        = r_l_obst_filechar
   l_obst_init_spatial(iv)    = r_l_obst_init_spatial
   l_obst_filedistri(iv)      = r_l_obst_filedistri
   obst_fn_char(iv)           = r_obst_fn_char
   obst_fn_distrib(iv)        = r_obst_fn_distrib
   obst_i_height(iv)          = r_obst_i_height
   obst_i_width(iv)           = r_obst_i_width
   obst_i_thick(iv)           = r_obst_i_thick
   obst_i_dens(iv)            = r_obst_i_dens
   ! * For namelist obst_var_flexibility
   l_obst_abdelposture(iv)    = r_l_obst_abdelposture
   l_obst_param_height(iv)    = r_l_obst_param_height
   obst_c_abdel_nmax(iv)      = r_obst_c_abdel_nmax
   obst_c_rho(iv)             = r_obst_c_rho
   obst_c_lift(iv)            = r_obst_c_lift
   obst_c_shelter(iv)         = r_obst_c_shelter
   obst_c_height_x0(iv)       = r_obst_c_height_x0
   obst_c_height_x1(iv)       = r_obst_c_height_x1
   ! * For namelist obst_var_roughdrag
   l_obst_drag_cste(iv)       = r_l_obst_drag_cste
   l_obst_abdelrough_cste(iv) = r_l_obst_abdelrough_cste
   obst_c_crough_x0(iv)       = r_obst_c_crough_x0
   obst_c_crough_x1(iv)       = r_obst_c_crough_x1
   obst_c_drag(iv)            = r_obst_c_drag
   obst_c_lz(iv)              = r_obst_c_lz
   ! * For namelist obst_var_fracxy
   l_obst_fracxy(iv)          = r_l_obst_fracxy
   obst_fracxy_type(iv)       = r_obst_fracxy_type
   obst_c_fracxy_k0(iv)       = r_obst_c_fracxy_k0
   obst_c_fracxy_k1(iv)       = r_obst_c_fracxy_k1
   obst_c_fracxy_l(iv)        = r_obst_c_fracxy_l
   ! * For namelist obst_var_bstress
   l_obst_z0bstress(iv)       = r_l_obst_z0bstress
   obst_z0bstress_option(iv)  = r_obst_z0bstress_option
   obst_c_z0bstress(iv)       = r_obst_c_z0bstress
   obst_c_z0bstress_x0(iv)    = r_obst_c_z0bstress_x0
   obst_c_z0bstress_x1(iv)    = r_obst_c_z0bstress_x1
   obst_c_z0bstress_x2(iv)    = r_obst_c_z0bstress_x2
   !!********************************** 
   PRINT_DBG*, 'END OBSTRUCTIONS_READVAR'
   END SUBROUTINE OBSTRUCTIONS_readvar

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_write_summary
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_write_summary  ***
   !&E
   !&E ** Purpose : Write a summary of obstructions parameters
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : OBSTRUCTIONS_init
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
   !! * Modules used
   USE parameters,  ONLY  : rsh
   USE comvars2d,   ONLY  : iscreenlog
   USE_MPI toolmpi, ONLY  : MASTER

   IMPLICIT NONE

   !! * Local declaration
   INTEGER          :: iv
   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_WRITE_SUMMARY'
   !***********************************************
   IF_MPI(MASTER) THEN
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) '***********************************************************************'
     WRITE(iscreenlog,*) '********************** module OBSTRUCTIONS ****************************'
     WRITE(iscreenlog,*) '************** subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
     WRITE(iscreenlog,*) '***********************************************************************'
     WRITE(iscreenlog,*) 'LISTING OF OBSTRUCTION VARIABLES :'
     WRITE(iscreenlog,*) '------------------------------------------------'
     WRITE(iscreenlog,*) 'TOTAL NUMBER OF OBSTRUCTION VARIABLES : ',obst_nbvar
     WRITE(iscreenlog,*) 'Number of RIGID_UP                    : ',obst_nv_rigid_up
     WRITE(iscreenlog,*) 'Number of RIGID_DO                    : ',obst_nv_rigid_do
     WRITE(iscreenlog,*) 'Number of FLEXI_UP                    : ',obst_nv_flexi_up
     WRITE(iscreenlog,*) 'Number of FLEXI_DO                    : ',obst_nv_flexi_do
     WRITE(iscreenlog,*) 'Number of 3DVARS                      : ',obst_nv_3d
     WRITE(iscreenlog,*) '------------------------------------------------'
     WRITE(iscreenlog,*) 'File for obstruction position is : ',TRIM(obst_fn_position)
     DO iv = 1,obst_nbvar
       WRITE(iscreenlog,*) '***********************************************************************************'
       WRITE(iscreenlog,*) '!===========================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_main  !'
       WRITE(iscreenlog,*) '!===========================!'
       WRITE(iscreenlog,*) '  - Number (identifier) of the variable  : ',obst_varnum(iv)
       WRITE(iscreenlog,*) '  - Name (identifier) of the variable    : ',obst_varname(iv)
       WRITE(iscreenlog,*) '!=============================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_option  !'
       WRITE(iscreenlog,*) '!=============================!'
       WRITE(iscreenlog,*) '  - If the current variable is cylinder-like                : ',l_obst_cylindre(iv)
       WRITE(iscreenlog,*) '  - If the current variable is flexible                     : ',l_obst_flexible(iv)
       WRITE(iscreenlog,*) '  - If the current variable is downvard                     : ',l_obst_downward(iv)
       WRITE(iscreenlog,*) '  - If the current variable is full 3D                      : ',l_obst_3dobst(iv)
       WRITE(iscreenlog,*) '  - If the current variable is considered a macro-roughness : ',l_obst_noturb(iv)
       WRITE(iscreenlog,*) '!===========================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_init  !'
       WRITE(iscreenlog,*) '!===========================!'
       WRITE(iscreenlog,*) '  - Use a time-series file for obstructions characteristics                : ',l_obst_filechar(iv)
       WRITE(iscreenlog,*) '  - Use a spatially variable file obstructions charcateristics             : ',l_obst_init_spatial(iv)
       WRITE(iscreenlog,*) '  - Use a file describing the vertical distribution of obstruction density : ',l_obst_filedistri(iv)
       IF(l_obst_filechar(iv).OR.l_obst_init_spatial(iv))THEN
         WRITE(iscreenlog,*) '  - Name of temporal or spatial file for obstructions charcateristics      : ',obst_fn_char(iv)
       ELSE
         WRITE(iscreenlog,*) '  - Name of temporal or spatial file for obstructions charcateristics      : NOT USED'
       ENDIF
       IF(l_obst_filedistri(iv))THEN
         WRITE(iscreenlog,*) '  - Name of file for the vertical distribution of obstruction density      : ',obst_fn_distrib(iv)
       ELSE
         WRITE(iscreenlog,*) '  - Name of file for the vertical distribution of obstruction density      : NOT USED'
       ENDIF
       IF(l_obst_filechar(iv).OR.l_obst_init_spatial(iv))THEN
         WRITE(iscreenlog,*) '  - Initial height (unbent, eg. leaf-length for segrasses) of obstructions : NOT USED'
         WRITE(iscreenlog,*) '  - Initial width (or diameter for cylindric obstructions) of obstructions : NOT USED'
         WRITE(iscreenlog,*) '  - Initial thick of obstructions (along the flow)                         : NOT USED'
         WRITE(iscreenlog,*) '  - Initial density of obstructions (or maximum density)                   : NOT USED'
       ELSE
         WRITE(iscreenlog,*) '  - Initial height (unbent, eg. leaf-length for segrasses) of obstructions : ',obst_i_height(iv)
         WRITE(iscreenlog,*) '  - Initial width (or diameter for cylindric obstructions) of obstructions : ',obst_i_width(iv)
         WRITE(iscreenlog,*) '  - Initial thick of obstructions (along the flow)                         : ',obst_i_thick(iv)
         WRITE(iscreenlog,*) '  - Initial density of obstructions (or maximum density)                   : ',obst_i_dens(iv)
       ENDIF
       WRITE(iscreenlog,*) '!==================================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_flexibility  !'
       WRITE(iscreenlog,*) '!==================================!'
       IF(l_obst_flexible(iv))THEN
         WRITE(iscreenlog,*) '  - Use Abdhelhrmans (2007) procedure to compute bending                   : ',l_obst_abdelposture(iv)
         WRITE(iscreenlog,*) '  - Use empirical (exponential desrease) formulation to compute bending    : ',l_obst_param_height(iv)
         IF(l_obst_abdelposture(iv))THEN
           WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : ',obst_c_abdel_nmax(iv)
           WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : ',obst_c_rho(iv)
           WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : ',obst_c_lift(iv)
           WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : ',obst_c_shelter(iv)
           WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : NOT USED'
           WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : NOT USED'
         ELSEIF(l_obst_param_height(iv))THEN
           WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
           WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
           WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
           WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
           WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : ',obst_c_height_x0(iv)
           WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : ',obst_c_height_x1(iv)
         ELSE
           WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
           WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
           WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
           WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
           WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : NOT USED'
           WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : NOT USED'
         ENDIF
       ELSE
         WRITE(iscreenlog,*) '  - Use Abdhelhrmans (2007) procedure to compute bending                   : NOT USED'
         WRITE(iscreenlog,*) '  - Use empirical (exponential desrease) formulation to compute bending    : NOT USED'
         WRITE(iscreenlog,*) '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
         WRITE(iscreenlog,*) '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
         WRITE(iscreenlog,*) '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
         WRITE(iscreenlog,*) '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
         WRITE(iscreenlog,*) '  - First parameter for empirical formulation                              : NOT USED'
         WRITE(iscreenlog,*) '  - Second parameter for empirical formulation                             : NOT USED'
       ENDIF
       WRITE(iscreenlog,*) '!================================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_roughdrag  !'
       WRITE(iscreenlog,*) '!================================!'
       WRITE(iscreenlog,*) '  - Use a constant drag coefficient for obstructions in hydrodynamics      : ',l_obst_drag_cste(iv)
       WRITE(iscreenlog,*) '  - Use a constant drag during reconstruction of velocity profiles         : ',l_obst_abdelrough_cste(iv)
       IF(l_obst_abdelrough_cste(iv))THEN
         WRITE(iscreenlog,*) '  - First coefficient for CD during velocity profiles reconstruction       : NOT USED'
         WRITE(iscreenlog,*) '  - Second coefficient for CD during velocity profiles reconstruction      : NOT USED'
       ELSE
         WRITE(iscreenlog,*) '  - First coefficient for CD during velocity profiles reconstruction       : ',obst_c_crough_x0(iv)
         WRITE(iscreenlog,*) '  - Second coefficient for CD during velocity profiles reconstruction      : ',obst_c_crough_x1(iv)
       ENDIF
       WRITE(iscreenlog,*) '  - Drag coefficient (max value if not constant) for obstructions elements : ',obst_c_drag(iv)
       IF(l_obst_noturb(iv))THEN
         WRITE(iscreenlog,*) '  - Coef. turbulent dissipation time-scale between obstructions elements   : NOT USED'
       ELSE
         WRITE(iscreenlog,*) '  - Coef. turbulent dissipation time-scale between obstructions elements   : ',obst_c_lz(iv)
       ENDIF
       WRITE(iscreenlog,*) '!=============================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_fracxy  !'
       WRITE(iscreenlog,*) '!=============================!'
       WRITE(iscreenlog,*) '  - Take account for non-linear patchiness correction                      : ',l_obst_fracxy(iv)
       IF(.NOT.l_obst_fracxy(iv))THEN
         WRITE(iscreenlog,*) '  - Kind of non-linear correction method                                   : NOT USED'
         WRITE(iscreenlog,*) '  - Coefficient for the exponential correction                             : NOT USED'
         WRITE(iscreenlog,*) '  - First parameter for correction of the exponential coefficient          : NOT USED'
         WRITE(iscreenlog,*) '  - Second parameter for correction of the exponential coefficient         : NOT USED'
       ELSE
         IF(obst_fracxy_type(iv).EQ.1)THEN
           WRITE(iscreenlog,*) '  - Kind of non-linear correction method                                   : Exponential'
           WRITE(iscreenlog,*) '  - Coefficient for the exponential correction                             : ',obst_c_fracxy_k0(iv)
           WRITE(iscreenlog,*) '  - First parameter for correction of the exponential coefficient          : NOT USED'
           WRITE(iscreenlog,*) '  - Second parameter for correction of the exponential coefficient         : NOT USED'
         ELSEIF(obst_fracxy_type(iv).EQ.2)THEN
           WRITE(iscreenlog,*) '  - Kind of non-linear correction method                                   : Exponential K0-variable'
           WRITE(iscreenlog,*) '  - Coefficient for the exponential correction                             : ',obst_c_fracxy_k0(iv)
           WRITE(iscreenlog,*) '  - First parameter for correction of the exponential coefficient          : ',obst_c_fracxy_k1(iv)
           WRITE(iscreenlog,*) '  - Second parameter for correction of the exponential coefficient         : ',obst_c_fracxy_l(iv)
         ENDIF
       ENDIF
       WRITE(iscreenlog,*) '!==============================!'
       WRITE(iscreenlog,*) '! NAMELIST : obst_var_bstress  !'
       WRITE(iscreenlog,*) '!==============================!'
       WRITE(iscreenlog,*) '  - To activate the impact of obstruction on Z0 (BSS computation)          : ',l_obst_z0bstress(iv)
       IF(l_obst_z0bstress(iv))THEN
         IF(obst_z0bstress_option(iv).EQ.0)THEN
           WRITE(iscreenlog,*) '  - Option to compute the obstruction induced roughness length             : Constant'
           WRITE(iscreenlog,*) '  - Constant (corrected value of roughness length)                         : ',obst_c_z0bstress(iv)
           WRITE(iscreenlog,*) '  - First parameter for rouhgness length computation (in 3D)               : NOT USED'
           WRITE(iscreenlog,*) '  - Second parameter for rouhgness length computation (in 3D)              : NOT USED'       
         ELSEIF(obst_z0bstress_option(iv).EQ.1)THEN
           WRITE(iscreenlog,*) '  - Option to compute the obstruction induced roughness length             : Parameterized'
           WRITE(iscreenlog,*) '  - Constant (corrected value of roughness length)                         : NOT USED'
           WRITE(iscreenlog,*) '  - First parameter for rouhgness length computation (in 3D)               : ',obst_c_z0bstress_x0(iv)
           WRITE(iscreenlog,*) '  - Second parameter for rouhgness length computation (in 3D)              : ',obst_c_z0bstress_x1(iv)
         ENDIF
         WRITE(iscreenlog,*) '  - Coefficient to correct 3D roughness length into 2D roughness length    : ',obst_c_z0bstress_x2(iv)
       ELSE
         WRITE(iscreenlog,*) '  - Option to compute the obstruction induced roughness length             : NOT USED'
         WRITE(iscreenlog,*) '  - Constant (corrected value of roughness length)                         : NOT USED'
         WRITE(iscreenlog,*) '  - First parameter for rouhgness length computation (in 3D)               : NOT USED'
         WRITE(iscreenlog,*) '  - Second parameter for rouhgness length computation (in 3D)              : NOT USED'
         WRITE(iscreenlog,*) '  - Coefficient to correct 3D roughness length into 2D roughness length    : NOT USED'
       ENDIF
     ENDDO
     WRITE(iscreenlog,*) '***************************************************************'
   ENDIF_MPI
   !-------------------------------------------
   PRINT_DBG*, 'END OBSTRUCTIONS_WRITE_SUMMARY'
   END SUBROUTINE OBSTRUCTIONS_write_summary

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_compatibility
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_compatibility ***
   !&E
   !&E ** Purpose : Check compatibility of parameterization used
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : OBSTRUCTIONS_init
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
   !&E       ! 2021-10    (F. Ganthy) Original code (extracted from OBSTRUCTIONS_write_summary)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters,  ONLY  : rsh
   USE comvars2d,   ONLY  : iscreenlog,iwarnlog,ierrorlog,l_modele2d
   USE comvars3d,   ONLY  : turb_nbeq,turb_2eq_option
   USE_MPI toolmpi, ONLY  : MASTER

   IMPLICIT NONE

   !! * Local declaration
   INTEGER          :: iv,nv_tot
   LOGICAL          :: IERR_MPI,l_turb

   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_COMPATIBILITY'
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
     IF (l_obst_noturb(iv)) THEN
         obst_nv_noturb     = obst_nv_noturb+1
     ELSE
         obst_nv_turb       = obst_nv_turb+1
     ENDIF
     IF (l_obst_3dobst(iv)) THEN
       obst_nv_3d           = obst_nv_3d+1
     ELSE
       IF(l_obst_downward(iv))THEN
         obst_nv_do         = obst_nv_do+1
         IF(l_obst_flexible(iv))THEN
           obst_nv_flexi_do = obst_nv_flexi_do+1
         ELSE
           obst_nv_rigid_do = obst_nv_rigid_do+1
         ENDIF
       ELSE
         obst_nv_up         = obst_nv_up+1
         IF(l_obst_flexible(iv))THEN
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
     IF_MPI(MASTER) THEN
       WRITE(ierrorlog,*) ' '
       WRITE(ierrorlog,*) ' '
       WRITE(ierrorlog,*) '**********************************************************************'
       WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
       WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
       WRITE(ierrorlog,*) ' ERROR : The total number of obstruction variables read from'
       WRITE(ierrorlog,*) '         file',TRIM(obst_fn_vardat(iv))
       WRITE(ierrorlog,*) '         is DIFFERENT from the number of variables defined'
       WRITE(ierrorlog,*) '         nbvar',obst_nbvar
       WRITE(ierrorlog,*) '         nv_tot',nv_tot
       WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
       WRITE(ierrorlog,*) '**********************************************************************'
       CALL_MPI MPI_FINALIZE(IERR_MPI)
       STOP
     ENDIF_MPI
   ENDIF
   !------------------------
   ! TESTING VARIABLES ORDER
   !------------------------
   DO iv = 1,obst_nbvar
     IF(obst_varnum(iv) /= iv) THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Wrong obstruction variable order'
         WRITE(ierrorlog,*) '         Variable number ', obst_varnum(iv)
         WRITE(ierrorlog,*) '         corresponds to the ', iv, ' variable read'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
   ENDDO
   !---------------------------------------------------
   ! TEST FOR VARIABLES USING FULL TURBULENCE PROCEDURE
   !---------------------------------------------------
   IF(obst_nv_turb /=0) THEN
     IF(l_modele2d)THEN
       IF_MPI(MASTER) THEN
         WRITE(iwarnlog,*) ' '
         WRITE(iwarnlog,*) ' '
         WRITE(iwarnlog,*) '**********************************************************************'
         WRITE(iwarnlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(iwarnlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(iwarnlog,*) ' WARNING : At least one obstruction variable use turbulence'
         WRITE(iwarnlog,*) '           procedure while you are running a 2D model'
         WRITE(iwarnlog,*) '**********************************************************************'
       ENDIF_MPI
     ENDIF
     IF(turb_nbeq /= 2) THEN
       IF_MPI(MASTER) THEN
         WRITE(iwarnlog,*) ' '
         WRITE(iwarnlog,*) ' '
         WRITE(iwarnlog,*) '**********************************************************************'
         WRITE(iwarnlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(iwarnlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(iwarnlog,*) ' WARNING : At least one obstruction variable use turbulence'
         WRITE(iwarnlog,*) '           procedure, while the wrong number of equation for'
         WRITE(iwarnlog,*) '           turbulence should be : turb_nbeq=2'
         WRITE(iwarnlog,*) '**********************************************************************'
       ENDIF_MPI
     ENDIF
     IF(turb_2eq_option /= 2) THEN
       IF_MPI(MASTER) THEN
         WRITE(iwarnlog,*) ' '
         WRITE(iwarnlog,*) ' '
         WRITE(iwarnlog,*) '**********************************************************************'
         WRITE(iwarnlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(iwarnlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(iwarnlog,*) ' WARNING : At least one obstruction variable use turbulence'
         WRITE(iwarnlog,*) '           procedure, while the wrong 2 equations turbulence'
         WRITE(iwarnlog,*) '           option should be : turb_2eq_option=2'
         WRITE(iwarnlog,*) '**********************************************************************'
       ENDIF_MPI
     ENDIF
   ENDIF
   !---------------------------
   ! TEST FOR FULL 3D VARIABLES
   !---------------------------
   IF(obst_nv_3d/=0) THEN
     IF(l_modele2d)THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : At least one obstruction variable is a full 3D'
         WRITE(ierrorlog,*) '         variable while you are running a 2D model'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF(turb_nbeq /= 2) THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : At least one obstruction variable is a full 3D'
         WRITE(ierrorlog,*) '         variable while the wrong number of equation for'
         WRITE(ierrorlog,*) '         turbulence should be : turb_nbeq=2'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF(turb_2eq_option /= 2) THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : At least one obstruction variable is a full 3D'
         WRITE(ierrorlog,*) '         variable while the wrong 2 equations turbulence'
         WRITE(ierrorlog,*) '         option should be : turb_2eq_option=2'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
   ENDIF
   !---------------------------------------------------
   ! TEST FOR VARIABLES USING FULL TURBULENCE PROCEDURE
   !---------------------------------------------------
   IF(obst_nv_do /=0) THEN
     IF(l_modele2d)THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '************************* module OBSTRUCTIONS ************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : At least one obstruction variable is a downward'
         WRITE(ierrorlog,*) '         variable while you are running a 2D model'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
   ENDIF
   !--------------------------------------------------------------
   ! TEST CONSISTENCY FOR SCHEME AND DESCRIPTION FOR EACH VARIALBE
   !--------------------------------------------------------------
   DO iv=1,obst_nbvar
     IF((l_obst_downward(iv)).AND.(l_obst_noturb(iv)))THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
         WRITE(ierrorlog,*) '         This variable is defined as a downward one, while'
         WRITE(ierrorlog,*) '         it uses simplified (l_obst_noturb) procedure'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF((l_obst_3dobst(iv)).AND.(l_obst_noturb(iv)))THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '************************* module OBSTRUCTIONS ************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
         WRITE(ierrorlog,*) '         This variable is defined as a full 3d one, while'
         WRITE(ierrorlog,*) '         it uses simplified (l_obst_noturb) procedure'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF((l_obst_3dobst(iv)).AND.(l_obst_flexible(iv)))THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
         WRITE(ierrorlog,*) '         This variable is defined as a full 3d one, while'
         WRITE(ierrorlog,*) '         it is also defined as a flexible one'
         WRITE(ierrorlog,*) '         This kind of variable is not available yet...'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF((l_obst_3dobst(iv)).AND.(l_obst_downward(iv)))THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
         WRITE(ierrorlog,*) '         This variable is defined as a full 3d one, while'
         WRITE(ierrorlog,*) '         it is also defined as a downward one'
         WRITE(ierrorlog,*) '         This kind of variable is not available yet...'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF((l_obst_downward(iv)).AND.(l_obst_flexible(iv)).AND.(l_obst_abdelposture(iv)))THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
         WRITE(ierrorlog,*) '         This variable is defined as a flexible and downward one'
         WRITE(ierrorlog,*) '         and uses Abdelhrman procedure'
         WRITE(ierrorlog,*) '         This procedure is not available yet for downward variables...'
         WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
     IF(l_obst_init_spatial(iv).AND.l_obst_filechar(iv))THEN
       IF_MPI(MASTER) THEN
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '**********************************************************************'
         WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
         WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
         WRITE(ierrorlog,*) ' ERROR : Inconsistency for variable :',obst_varname(iv)
         WRITE(ierrorlog,*) '         Spatial initialization is defined with also spatial forcing'
         WRITE(ierrorlog,*) '         for obstruction characteristics : This is not ready yet...'
         WRITE(ierrorlog,*) '         Choose between : spatial initialization (but no time-varying characteristics)'
         WRITE(ierrorlog,*) '         and homogeneous initialization (but with time-varying characteristics)'
         WRITE(ierrorlog,*) '**********************************************************************'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     ENDIF
   ENDDO
   !----------------------------
   ! TEST ON OMP PARALLELIZATION
   !----------------------------
#ifdef key_MPIOMP
   IF_MPI(MASTER) THEN
     WRITE(ierrorlog,*) ' '
     WRITE(ierrorlog,*) ' '
     WRITE(ierrorlog,*) '**********************************************************************'
     WRITE(ierrorlog,*) '********************** module OBSTRUCTIONS ***************************'
     WRITE(ierrorlog,*) '************* subroutine OBSTRUCTIONS_WRITE_SUMMARY ******************'
     WRITE(ierrorlog,*) 'ERROR : OMP parallelization not available yet for OBSTRUCTIONS'
     WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!!'
     WRITE(ierrorlog,*) '**********************************************************************'
     CALL_MPI MPI_FINALIZE(IERR_MPI)
     STOP
   ENDIF_MPI
#endif
   !-------------------------------------------
   PRINT_DBG*, 'END OBSTRUCTIONS_COMPATIBILITY'
   END SUBROUTINE OBSTRUCTIONS_compatibility

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readfile_char
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readfile_char  ***
   !&E
   !&E ** Purpose : Readinf .dat file for time-varying obstructions characteristics
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : obst_init & obst_update
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
   !&E       ! 2014-10-16 (F. Ganthy) Replacing diameter by width and thickness  
   !&E       ! 2015-01-15 (F. Ganthy) MPI parallelization
   !&E       ! 2016-09-20 (F. Ganthy) Change computation of obstruction area index to be more 
   !&E       !                        consistent with biological aspects
   !&E       ! 2017-02-16 (F. Ganthy) Changes on instantaneous obstruction state variables for 
   !&E                                future coupling with Zostera growth module
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh,limin,limax,jmin,jmax,kmax,pi
   USE comvars2d
   USE_MPI toolmpi, ONLY : MASTER

   IMPLICIT NONE

   !! * Local declaration
   LOGICAL               :: ex,IERR_MPI
   INTEGER               :: i,j,eof,kk,iv,numfile
   CHARACTER(LEN=lchain) :: rec
   REAL(KIND=rlg)        :: tool_datosec,tint1,tint2,dt1,dt2,t1,tdb,tfi
   REAL(KIND=rsh)        :: height1,height2,width1,width2,thick1,thick2,dens1,dens2

   !!--------------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_READFILE_CHAR'
   !-------------------------
   ! **** Initialization ****
   !-------------------------
   IF(nbouc.EQ.0)THEN
     !--------------------
     ! dynamic allocation
     !--------------------
     ALLOCATE(obst_dens_t   (obst_nbvar,nvalmax))
     ALLOCATE(obst_width_t  (obst_nbvar,nvalmax))
     ALLOCATE(obst_thick_t  (obst_nbvar,nvalmax))
     ALLOCATE(obst_height_t (obst_nbvar,nvalmax))
     !--------------------------
     !*** Start reading file ***
     !--------------------------
     DO iv = 1,obst_nbvar
       IF (l_obst_filechar(iv))THEN
         IF_MPI(MASTER) THEN
           WRITE(iscreenlog,*) ' '
           WRITE(iscreenlog,*) ' '
           WRITE(iscreenlog,*) '**********************************************************************'
           WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_CHAR *****'
           WRITE(iscreenlog,*) '**********************************************************************'
           WRITE(iscreenlog,*) ' Reading characteristics file ',obst_fn_char(iv)
           WRITE(iscreenlog,*) ' Corresponding to obstructions variable ',obst_varname(iv)
           WRITE(iscreenlog,*) '**********************************************************************'
         ENDIF_MPI
         !---------------------------
         ! TEST FOR FILE AVAILABILITY
         !---------------------------
         eof=0
         INQUIRE(file=obst_fn_char(iv),exist=ex)
         IF(ex) THEN
           numfile=1770+iv
           OPEN(numfile,file = obst_fn_char(iv),form = 'formatted')
           i=1
           !-----------------------------------------------
           ! Check the first available date within the file
           !-----------------------------------------------
           READ(numfile,*)
           READ(numfile,*)
           READ(numfile,*)
           READ(numfile,*)
           READ(numfile,*)
           READ(numfile,*)
           READ(numfile,*)
           READ(numfile,'(a)',END=999) rec
           kk = INDEX(rec,'/')-2
           READ(rec(kk+19:kk+27),*) height1
           READ(rec(kk+28:kk+36),*) width1
           READ(rec(kk+37:kk+45),*) thick1
           READ(rec(kk+46:),*) dens1
           tint1 = tool_datosec(rec(kk:kk+18))
           t1 = tchrono(i)
           DO WHILE (t1 < tint1.AND.i<=nvalmax)
             obst_height_t(iv,i) = height1
             obst_width_t(iv,i) = width1
             obst_thick_t(iv,i) = thick1
             obst_dens_t(iv,i) = dens1
             i=i+1
             IF(i <= nvalmax) t1=tchrono(i)
           ENDDO
           !---------------------------
           ! Time interpolation of data
           !---------------------------
           IF (i < nvalmax) THEN
             eof = 0
             DO WHILE (eof==0)
               READ(numfile,'(a)',END=999) rec
               kk = INDEX(rec,'/')-2
               READ(rec(kk+19:kk+27),*) height2
               READ(rec(kk+28:kk+36),*) width2
               READ(rec(kk+37:kk+45),*) thick2
               READ(rec(kk+46:),*) dens2
               tint2 = tool_datosec(rec(kk:kk+18))
               DO WHILE ((tchrono(i)-tint1)*(tchrono(i)-tint2) <= 0.0_rlg)
                 dt1 = (tchrono(i)-tint1)/(tint2-tint1)
                 dt2 = (tint2-tchrono(i))/(tint2-tint1)
                 obst_height_t(iv,i) = height1*dt2+height2*dt1
                 obst_width_t(iv,i)  = width1*dt2+width2*dt1
                 obst_thick_t(iv,i)  = thick1*dt2+thick2*dt1
                 obst_dens_t(iv,i)   = dens1*dt2+dens2*dt1
                 i=i+1
                 IF (i > nvalmax) GOTO 999
               ENDDO
               tint1   = tint2
               height1 = height2
               width1  = width2
               thick1  = thick2
               dens1   = dens2
             ENDDO
           ENDIF

   999 CONTINUE

           IF (i < nvalmax) THEN
             obst_height_t(iv,i) = height1
             obst_width_t(iv,i)  = width1
             obst_thick_t(iv,i)  = thick1
             obst_dens_t(iv,i)   = dens1
           ENDIF
          CLOSE(numfile)
        ELSE
          IF_MPI (MASTER) THEN
            WRITE(ierrorlog,*) ' '
            WRITE(ierrorlog,*) ' '
            WRITE(ierrorlog,*) '**************************************************************'
            WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBST_READFILE_CHAR *****'
            WRITE(ierrorlog,*) ' ERROR : File ', TRIM(obst_fn_char(iv)),' does not exist !!! '
            WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
            WRITE(ierrorlog,*) '**************************************************************'
          CALL_MPI MPI_FINALIZE(IERR_MPI)
          STOP
          ENDIF_MPI
        ENDIF ! end test of file availability
      ELSE ! NO FILE FOR TIME-VARYING CHARACTERISTICS APPLYING INITIAL PARAMETERS
        obst_height_t(iv,:) = obst_i_height(iv)
        obst_width_t(iv,:)  = obst_i_width(iv)
        obst_thick_t(iv,:)  = obst_i_thick(iv)
        obst_dens_t(iv,:)   = obst_i_dens(iv)
      ENDIF ! END TEST ON CHARACTERISTICS FROM A FILE
     ENDDO ! END LOOP ON OBSTRUCTION VARIABLES
   ENDIF ! end test of nbouc
   !-------------------------------------
   ! **** End of initialization part ****
   !-------------------------------------
   !
   !-------------------------
   ! *** Current progress ***
   !-------------------------
   dt1 = (t-tchrono(indchrono))/(tchrono(indchrono+1)-tchrono(indchrono))
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv = 1,obst_nbvar
         IF(.NOT.l_obst_init_spatial(iv))THEN
           IF(obst_position(iv,i,j).GT.0.0_rsh) THEN
             obst_height_inst(iv,i,j)      = dt1*obst_height_t(iv,indchrono+1) + (1.d0-dt1)*obst_height_t(iv,indchrono)
             obst_width_inst(iv,i,j)       = dt1*obst_width_t(iv,indchrono+1) + (1.d0-dt1)*obst_width_t(iv,indchrono)
             obst_thick_inst(iv,i,j)       = dt1*obst_thick_t(iv,indchrono+1) + (1.d0-dt1)*obst_thick_t(iv,indchrono)
             obst_dens_inst(iv,i,j)        = dt1*obst_dens_t(iv,indchrono+1) + (1.d0-dt1)*obst_dens_t(iv,indchrono)
             IF(l_obst_cylindre(iv))THEN ! Cylindric/Ellipse obstruction
               obst_area_index_inst(iv,i,j)  = obst_dens_inst(iv,i,j) * obst_height_inst(iv,i,j) * &
                                               (2.0_rsh*pi*SQRT(0.5_rsh*(obst_width_inst(iv,i,j)**2.0_rsh + obst_thick_inst(iv,i,j)**2.0_rsh)))
             ELSE ! Parallelepipedic obstruction
               obst_area_index_inst(iv,i,j)  = 2.0_rsh * obst_dens_inst(iv,i,j) * obst_width_inst(iv,i,j) * obst_height_inst(iv,i,j)
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
             IF(l_obst_cylindre(iv))THEN ! Cylindric/Ellipse obstruction
               obst_area_index_inst(iv,i,j)  = obst_dens_inst(iv,i,j) * obst_height_inst(iv,i,j) * &
                                               (2.0_rsh*pi*SQRT(0.5_rsh*(obst_width_inst(iv,i,j)**2.0_rsh + obst_thick_inst(iv,i,j)**2.0_rsh)))
             ELSE ! Parallelepipedic obstruction
               obst_area_index_inst(iv,i,j)  = 2.0_rsh * obst_dens_inst(iv,i,j) * obst_width_inst(iv,i,j) * obst_height_inst(iv,i,j)
             ENDIF
           ELSE
             obst_area_index_inst(iv,i,j)  = 0.0_rsh
           ENDIF
         ENDIF
       ENDDO
      ENDDO
   ENDDO
   ! *********************************
   PRINT_DBG*, 'END OBSTRUCTIONS_READFILE_CHAR'
   END SUBROUTINE OBSTRUCTIONS_readfile_char

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readfile_pos
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readfile_pos  ***
   !&E
   !&E ** Purpose : read the obstruction position file
   !&E              in the historical DEL/AO format
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : obst_init
   !&E
   !&E ** External calls : 
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables : obst_position
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E
   !&E---------------------------------------------------------------------

#include "toolcpp.h"

   !! * Modules used
   USE parameters
   USE ionc4
   USE netcdf
   USE typeSizes
   USE comionc4
   USE comvars2d
   USE comvars3d, ONLY : sig
#ifdef key_tide_saverotatedobc
   USE tidesaveobcrotated, ONLY : lon2d_child,lat2d_child,    &
                                  lon2du_child,lon2dv_child,  &
                                  lat2du_child,lat2dv_child,  &
                                  hx_child,hy_child,h0_child, &
                                  ig_child,id_child,jh_child,jb_child
#endif
#if defined key_siggen || defined key_gencoord
   USE comsiggen, ONLY : hc,hc_sig_g,hcx_sig_g,hcy_sig_g
#endif
   USE_MPI toolmpi,   ONLY : MASTER


   IMPLICIT NONE

   !! * Local declarations
   LOGICAL                                                  :: IERR_MPI
   INTEGER                                                  :: i,j,iv
   REAL(KIND=rlg)                                           :: latmin_b,latmax_b,lonmin_b, &
                                                               lonmax_b,latmin_o,latmax_o, &
                                                               lonmin_o,lonmax_o
   REAL(KIND=rlg),DIMENSION(imin:imax,jmin:jmax)            :: lat_b,lon_b,lat_o,lon_o
   REAL(KIND=riolg),DIMENSION(imin:imax,jmin:jmax)          :: tmp2drlg
   REAL(KIND=riosh),DIMENSION(imin:imax,jmin:jmax)          :: zlit
   REAL(kind=rsh),DIMENSION(obst_nbvar,imin:imax,jmin:jmax) :: posobst_tmp

   INTEGER :: nc_err,nc_id,ndims,ngatts,xdimid,var_id,var_type,var_ndims,var_natts,nvars
   INTEGER,DIMENSION(10):: var_ndim
   CHARACTER(len=30) :: var_name
   LOGICAL :: test


   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_READFILE_POS'

   ! READING LAT LON FROM FILE_BATHY
   CALL ionc4_openr(file_bathy)
   CALL ionc4_read_xy(file_bathy,name_in_longitude,tmp2drlg,imin,imax,jmin,jmax)
   lon_b(:,:)=tmp2drlg(:,:)
   CALL ionc4_read_xy(file_bathy,name_in_latitude,tmp2drlg,imin,imax,jmin,jmax)
   lat_b(:,:)=tmp2drlg(:,:)
   CALL ionc4_close(file_bathy)

   latmax_b=maxval(maxval(lat_b,dim=2))
   latmin_b=minval(minval(lat_b,dim=2))
   lonmin_b=minval(minval(lon_b,dim=1))
   lonmax_b=maxval(maxval(lon_b,dim=1))

   ! READING LAT LON FROM FILE_OBSTRUCTIONS
   CALL ionc4_openr(obst_fn_position)
   CALL ionc4_read_xy(obst_fn_position,name_in_longitude,tmp2drlg,imin,imax,jmin,jmax)
   lon_o(:,:)=tmp2drlg(:,:)
   CALL ionc4_read_xy(obst_fn_position,name_in_latitude,tmp2drlg,imin,imax,jmin,jmax)
   lat_o(:,:)=tmp2drlg(:,:)
   CALL ionc4_close(obst_fn_position)

   latmax_o=maxval(maxval(lat_o,dim=2))
   latmin_o=minval(minval(lat_o,dim=2))
   lonmin_o=minval(minval(lon_o,dim=1))
   lonmax_o=maxval(maxval(lon_o,dim=1))

   ! TEST FOR VALID DOMAIN
   IF((ABS(latmax_b-latmax_o).GE.1.0E-3).OR.(ABS(latmin_b-latmin_o).GE.1.0E-3).OR.(ABS(lonmax_b-lonmax_o).GE.1.0E-3).OR.(ABS(lonmin_b-lonmin_o).GE.1.0E-3)) THEN
     WRITE(ierrorlog,*) ' '
     WRITE(ierrorlog,*) ' '
     WRITE(ierrorlog,*) '*********************************************************************'
     WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_POS *****'
     WRITE(ierrorlog,*) ' ERROR : When reading the file : ',TRIM(obst_fn_position)
     WRITE(ierrorlog,*) '         Latitude or longitude from the obstructions file'
     WRITE(ierrorlog,*) '         didn t match Latitude or longitude defined of the'
     WRITE(ierrorlog,*) '         bathymetric file ',TRIM(file_bathy)
     WRITE(ierrorlog,*) '         ---------------------'
     WRITE(ierrorlog,*) '         FROM BATHYMETRY FILE:'
     WRITE(ierrorlog,*) '         lonmin = ',lonmin_b
     WRITE(ierrorlog,*) '         lonmax = ',lonmax_b
     WRITE(ierrorlog,*) '         latmax = ',latmax_b
     WRITE(ierrorlog,*) '         latmin = ',latmin_b
     WRITE(ierrorlog,*) ' '
     WRITE(ierrorlog,*) '         FROM OBSTRUCTION FILE :'
     WRITE(ierrorlog,*) '         lonmin = ',lonmin_o
     WRITE(ierrorlog,*) '         lonmax = ',lonmax_o
     WRITE(ierrorlog,*) '         latmax = ',latmax_o
     WRITE(ierrorlog,*) '         latmin = ',latmin_o
     WRITE(ierrorlog,*) '         ---------------------'
     WRITE(ierrorlog,*) ' --> SIMUMLATION IS STOPPED'
     CALL_MPI MPI_FINALIZE(IERR_MPI)
     STOP
   ENDIF

   ! RETRIEVE AND CHECK VARIABLES NAME
   CALL ionc4_openr(obst_fn_position)
   CALL ionc4_corres(obst_fn_position,nc_id)
   nc_err = nf90_inquire(nc_id, ndims, nvars, ngatts, xdimid)
   test = .FALSE.
   DO iv=1,obst_nbvar
     ! READING OBSTRUCTIONS POSITIONS
     DO var_id=1,nvars
       nc_err = nf90_inquire_variable(nc_id, var_id, var_name,var_type, var_ndims,var_ndim, var_natts)
       IF(TRIM('Pos_')//TRIM(obst_varname(iv)).EQ.var_name) test=.TRUE.
     ENDDO
     IF (test) THEN
       CALL ionc4_read_xy(obst_fn_position,TRIM('Pos_')//TRIM(obst_varname(iv)),zlit,imin,imax,jmin,jmax)
       posobst_tmp(iv,:,:) = zlit(:,:)
       DO j=ljmin,ljmax
         DO i=limin,limax
           IF(posobst_tmp(iv,i,j).GT.0.0_rsh) obst_position(iv,i,j) = posobst_tmp(iv,i,j)
         ENDDO
       ENDDO
     ELSE
       WRITE(ierrorlog,*) ' '
       WRITE(ierrorlog,*) ' '
       WRITE(ierrorlog,*) '*********************************************************************'
       WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_POS *****'
       WRITE(ierrorlog,*) ' ERROR : When reading the file : ',TRIM(obst_fn_position)
       WRITE(ierrorlog,*) '         Names of obstructions variables within position file'
       WRITE(ierrorlog,*) '         didn t match names obstructions defined within file'
       WRITE(ierrorlog,*) '         ',TRIM(obst_fn_vardat(iv))
       WRITE(ierrorlog,*) ' --> SIMUMLATION IS STOPPED'
       CALL_MPI MPI_FINALIZE(IERR_MPI)
       STOP
     ENDIF
     IF(l_obst_init_spatial(iv))THEN
       ! READING OBSTRUCTIONS DENSITY
       DO var_id=1,nvars
         nc_err = nf90_inquire_variable(nc_id, var_id, var_name,var_type, var_ndims,var_ndim, var_natts)
         IF(TRIM('Dens_')//TRIM(obst_varname(iv)).EQ.var_name) test=.TRUE.
       ENDDO
       IF (test) THEN
         CALL ionc4_read_xy(obst_fn_position,TRIM('Dens_')//TRIM(obst_varname(iv)),zlit,imin,imax,jmin,jmax)
         DO j=ljmin,ljmax
           DO i=limin,limax
             IF(posobst_tmp(iv,i,j).GT.0.0_rsh) obst_dens_inst(iv,i,j) = zlit(i,j)
           ENDDO
         ENDDO
       ELSE
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '*********************************************************************'
         WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_POS *****'
         WRITE(ierrorlog,*) ' ERROR : When reading the file : ',TRIM(obst_fn_position)
         WRITE(ierrorlog,*) '         When reading spatial obstruction density'
         WRITE(ierrorlog,*) '         Names of obstructions variables within position file'
         WRITE(ierrorlog,*) '         didn t match names obstructions defined within file'
         WRITE(ierrorlog,*) '         ',TRIM(obst_fn_vardat(iv))
         WRITE(ierrorlog,*) ' --> SIMUMLATION IS STOPPED'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF
       ! READING OBSTRUCTIONS HEIGHT
       DO var_id=1,nvars
         nc_err = nf90_inquire_variable(nc_id, var_id, var_name,var_type, var_ndims,var_ndim, var_natts)
         IF(TRIM('Height_')//TRIM(obst_varname(iv)).EQ.var_name) test=.TRUE.
       ENDDO
       IF (test) THEN
         CALL ionc4_read_xy(obst_fn_position,TRIM('Height_')//TRIM(obst_varname(iv)),zlit,imin,imax,jmin,jmax)
         DO j=ljmin,ljmax
           DO i=limin,limax
             IF(posobst_tmp(iv,i,j).GT.0.0_rsh) obst_height_inst(iv,i,j) = zlit(i,j)
           ENDDO
         ENDDO
       ELSE
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '*********************************************************************'
         WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_POS *****'
         WRITE(ierrorlog,*) ' ERROR : When reading the file : ',TRIM(obst_fn_position)
         WRITE(ierrorlog,*) '         When reading spatial obstruction height'
         WRITE(ierrorlog,*) '         Names of obstructions variables within position file'
         WRITE(ierrorlog,*) '         didn t match names obstructions defined within file'
         WRITE(ierrorlog,*) '         ',TRIM(obst_fn_vardat(iv))
         WRITE(ierrorlog,*) ' --> SIMUMLATION IS STOPPED'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF
       ! READING OBSTRUCTIONS WIDTH
       DO var_id=1,nvars
         nc_err = nf90_inquire_variable(nc_id, var_id, var_name,var_type, var_ndims,var_ndim, var_natts)
         IF(TRIM('Width_')//TRIM(obst_varname(iv)).EQ.var_name) test=.TRUE.
       ENDDO
       IF (test) THEN
         CALL ionc4_read_xy(obst_fn_position,TRIM('Width_')//TRIM(obst_varname(iv)),zlit,imin,imax,jmin,jmax)
         DO j=ljmin,ljmax
           DO i=limin,limax
             IF(posobst_tmp(iv,i,j).GT.0.0_rsh) obst_width_inst(iv,i,j) = zlit(i,j)
           ENDDO
         ENDDO
       ELSE
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '*********************************************************************'
         WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_POS *****'
         WRITE(ierrorlog,*) ' ERROR : When reading the file : ',TRIM(obst_fn_position)
         WRITE(ierrorlog,*) '         When reading spatial obstruction width'
         WRITE(ierrorlog,*) '         Names of obstructions variables within position file'
         WRITE(ierrorlog,*) '         didn t match names obstructions defined within file'
         WRITE(ierrorlog,*) '         ',TRIM(obst_fn_vardat(iv))
         WRITE(ierrorlog,*) ' --> SIMUMLATION IS STOPPED'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF
       ! READING OBSTRUCTIONS THICKNESS
       DO var_id=1,nvars
         nc_err = nf90_inquire_variable(nc_id, var_id, var_name,var_type, var_ndims,var_ndim, var_natts)
         IF(TRIM('Thick_')//TRIM(obst_varname(iv)).EQ.var_name) test=.TRUE.
       ENDDO
       IF (test) THEN
         CALL ionc4_read_xy(obst_fn_position,TRIM('Thick_')//TRIM(obst_varname(iv)),zlit,imin,imax,jmin,jmax)
         DO j=ljmin,ljmax
           DO i=limin,limax
             IF(posobst_tmp(iv,i,j).GT.0.0_rsh) obst_thick_inst(iv,i,j) = zlit(i,j)
           ENDDO
         ENDDO
       ELSE
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) ' '
         WRITE(ierrorlog,*) '*********************************************************************'
         WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_POS *****'
         WRITE(ierrorlog,*) ' ERROR : When reading the file : ',TRIM(obst_fn_position)
         WRITE(ierrorlog,*) '         When reading spatial obstruction thickness'
         WRITE(ierrorlog,*) '         Names of obstructions variables within position file'
         WRITE(ierrorlog,*) '         didn t match names obstructions defined within file'
         WRITE(ierrorlog,*) '         ',TRIM(obst_fn_vardat(iv))
         WRITE(ierrorlog,*) ' --> SIMUMLATION IS STOPPED'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF
     ENDIF
   ENDDO

   CALL ionc4_close(obst_fn_position)

   PRINT_DBG*, 'OBSTRUCTIONS_readfile_pos'
   END SUBROUTINE OBSTRUCTIONS_readfile_pos

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_readfile_distri
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_readfile_distri  ***
   !&E
   !&E ** Purpose : Initialization of obstruction parameters
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
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh
   USE comvars2d, ONLY  : iscreenlog,iwarnlog,ierrorlog
   USE_MPI toolmpi, ONLY : MASTER

   IMPLICIT NONE

   !! * Local declaration
   INTEGER :: kk,iv,nb_max_hnorm
   LOGICAL :: ex,IERR_MPI
   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBSTRUCTIONS_READFILE_DISTRI'

   IF_MPI(MASTER) THEN
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) ' '
     WRITE(iscreenlog,*) '************************************************************************'
     WRITE(iscreenlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_DISTRI *****'
     WRITE(iscreenlog,*) '************************************************************************'
   ENDIF_MPI
   DO iv = 1,obst_nbvar
     IF(l_obst_filedistri(iv)) THEN
       IF_MPI(MASTER) THEN
         WRITE(iscreenlog,*) ' '
         WRITE(iscreenlog,*) 'Obstruction variable : ',obst_varname(iv)
         WRITE(iscreenlog,*) 'file defining vertical distribution of obstuctions :'
         WRITE(iscreenlog,*) TRIM(obst_fn_distrib(iv))
       ENDIF_MPI
       INQUIRE(file=obst_fn_distrib(iv),exist=ex)
       IF(ex) THEN
         OPEN(53,file=obst_fn_distrib(iv),form='formatted')
         READ(53,*) ! Filename
         READ(53,*) ! Title
         READ(53,*) obst_nbhnorm(iv)
         IF_MPI(MASTER) THEN
           WRITE(iscreenlog,*) ' '
           WRITE(iscreenlog,*) 'Number of vertical discretization : ',obst_nbhnorm(iv)
         ENDIF_MPI
         CLOSE(53)
       ELSE
         IF_MPI(MASTER) THEN
           WRITE(ierrorlog,*) ' '
           WRITE(ierrorlog,*) ' '
           WRITE(ierrorlog,*) '************************************************************************'
           WRITE(ierrorlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_READFILE_DISTRI *****'
           WRITE(ierrorlog,*) '************************************************************************'
           WRITE(ierrorlog,*) ' ERROR : File ',TRIM(obst_fn_distrib(iv)),'does not exist'
           WRITE(ierrorlog,*) ' --> THE SIMULATION IS STOPPED !!! '
           WRITE(ierrorlog,*) '************************************************************************'
           CALL_MPI MPI_FINALIZE(IERR_MPI)
           STOP
         ENDIF_MPI
       ENDIF
     ELSE
       obst_nbhnorm(iv) = 11
       IF_MPI(MASTER) THEN
         WRITE(iscreenlog,*) ' '
         WRITE(iscreenlog,*) 'Obstruction variable : ',obst_varname(iv)
         WRITE(iscreenlog,*) 'An homogene distribution is applied'
         WRITE(iscreenlog,*) ' '
         WRITE(iscreenlog,*) 'Number of vertical discretization : ',obst_nbhnorm(iv)
       ENDIF_MPI
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
     IF(l_obst_filedistri(iv)) THEN
       OPEN(53,file=obst_fn_distrib(iv),form='formatted')
       READ(53,*) ! Filename
       READ(53,*) ! Title
       READ(53,*) !NBhnorm
       READ(53,*) ! Column titles
       kk=1
       DO WHILE (kk.LE.obst_nbhnorm(iv))
         READ(53,*) obst_height_norm(iv,kk),obst_dens_norm(iv,kk)
         WRITE(iscreenlog,*) 'kk',kk
         WRITE(iscreenlog,*) 'height : ',obst_height_norm(iv,kk), ' dens : ',obst_dens_norm(iv,kk)
         kk = kk + 1
       ENDDO
       CLOSE(53)
     ELSE
       ! Homogene distribution of obstructions distribution
       ! Homogene distribution of obstructions distribution
#if defined key_casobstflume2D_Test06_3DO || defined key_casobstflume2D_Test07_Multi
       l_obst_filedistri(iv) = .TRUE.
#endif
       obst_dens_norm(iv,:) = 100.0_rsh
       obst_height_norm(iv,:) = 0.001_rsh
       DO kk=1,obst_nbhnorm(iv)
         obst_height_norm(iv,kk) = (REAL(kk-1))*10.0_rsh
         obst_dens_norm(iv,kk) = 100.0_rsh
#ifdef key_casobstflume2D_Test06_3DO
         IF(kk.LE.8) THEN
            obst_dens_norm(iv,kk) = 0.0_rsh
         ENDIF
#endif
#ifdef key_casobstflume2D_Test07_Multi
         IF((obst_varnum(iv).EQ.3).AND.(kk.LE.9))THEN
           obst_dens_norm(iv,kk) = 0.0_rsh
         ENDIF
#endif
       ENDDO
       obst_height_norm(iv,obst_nbhnorm(iv))=100.01_rsh
     ENDIF
   ENDDO

     PRINT_DBG*, 'END OBSTRUCTIONS_READFILE_DISTRI'
END SUBROUTINE OBSTRUCTIONS_readfile_distri

   !!==========================================================================================================

#endif
END MODULE initOBSTRUCTIONS
