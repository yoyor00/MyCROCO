MODULE init_HYBIOSED

#include "cppdefs.h"

#if defined HYBIOSED
   !!===========================================================================
   !!                   ***  MODULE  init_HYBIOSED  ***
   !!
   !! ** Purpose : Interactions and feedbacks between hydrodynamics, biology and
   !!              sediment dynamics
   !!
   !! ** Description :
   !!
   !!  subroutine hbs_init             ! Initialization
   !!
   !!===========================================================================

   !! * Modules used
   USE com_HYBIOSED
   USE HYBIOSED1DV, ONLY: hbs1dv_alloc, hbs1dv_init

   IMPLICIT NONE

   !! * Accessibility
   PUBLIC hbs_init

   ! default
   PRIVATE

CONTAINS

   !!===========================================================================

   SUBROUTINE hbs_init(h0fond_in)
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE HYBIOSED_init ***
   !!
   !! ** Description : Initialize variables
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for GLOBAL_2D_ARRAY, N, stdout, mynode
      USE comMUSTANG, ONLY: hsed ! hsed init by MUSTANG before HYBIOSED init
      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: h0fond_in

      INTEGER                    :: iv, k, lstr, lenstr

      NAMELIST /hbs_main/ hbs_nbvar
      NAMELIST /hbs_input/ hbs_fn_position, hbs_fn_var
      NAMELIST /hbs_obstsedim/ l_hbs_suspsed_trapp, hbs_c_suspsed_trapp_max, &
         hbs_c_suspsed_trapp_exp, &
         l_hbs_suspsed_block, hbs_c_suspsed_block_max, &
         hbs_c_suspsed_block_exp
      NAMELIST /hbs_output/ l_hbsout_pos, l_hbsout_susp_trapp, &
         l_hbsout_susp_block, l_hbsout_rbiom, l_hbsout_zroot, &
         l_hbsout_thickroot

      hbs_kmax = N
      iscreenlog = stdout
      ierrorlog = stdout
      iwarnlog = stdout

      ! ************************
      ! * READING NAMELIST
      ! ************************
      MPI_master_only WRITE (iscreenlog, *) &
         ' '
      MPI_master_only WRITE (iscreenlog, *) &
         '************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***** module HYBIOSED, subroutine HBS_INIT *****'
      MPI_master_only WRITE (iscreenlog, *) &
         '************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         ' Reading file ', TRIM(hbsname)
      lstr = lenstr(hbsname)
      OPEN (50, file=hbsname(1:lstr), status='old', &
            form='formatted', access='sequential')
      READ (50, hbs_main); REWIND (50)
      READ (50, hbs_obstsedim); REWIND (50)
      READ (50, hbs_output); REWIND (50)

      IF (hbs_nbvar > 0) THEN
         ALLOCATE (hbs_fn_var(1:hbs_nbvar))
         READ (50, hbs_input); REWIND (50)
      ELSE ! hbs_nbvar > 0
         MPI_master_only WRITE (iscreenlog, *) &
            ' No Hybiosed, hbs_nbvar = ', hbs_nbvar
      END IF ! hbs_nbvar > 0

      CLOSE (50)

      ! ************************************
      ! * ALLOCATION OF VARIABLES PARAMETERS
      ! ************************************
      CALL hbs_alloc_nbvar

      ! ***********************
      ! * READING VARIABLES
      ! ***********************
      DO iv = 1, hbs_nbvar
         CALL hbs_readvar(iv)
      END DO

      ! **********************
      ! * TABLES ALLOCATION
      ! ***********************
      CALL hbs_alloc_xyz

      ! **********************
      ! * CHECKS AND SUMMARY
      ! ***********************
      CALL hbs_compatibility
      CALL hbs_write_summary

      ! *******************************
      ! * INITIALIZATION OF OUTPUT
      ! *******************************
      CALL hbs_vname

      ! **********************
      ! * INITIALIZE VARIABLES
      ! **********************
      IF (hbs_nbvar > 0) THEN
         ! ***********************
         ! * READING POSITION FILE
         ! ***********************
         CALL hbs_init_pos
         ! ********************************************************
         ! * READING SPATIAL INIT FILE or init from constant values
         ! ********************************************************
         CALL hbs_init_BRoot_ZDZroot
      END IF

      ! ***********************
      ! * OTHER INITIALIZATIONS
      ! ***********************
      CALL hbs1dv_alloc(hbs_nbvar, hbs_kmax)
      CALL hbs1dv_init(h0fond_in, l_hbs_suspsed_trapp, &
                       hbs_c_suspsed_trapp_max, hbs_c_suspsed_trapp_exp, &
                       l_hbs_suspsed_block, hbs_c_suspsed_block_max, &
                       hbs_c_suspsed_block_exp, &
                       hbs_c_tauceroot_x0, hbs_c_tauceroot_x1, &
                       hbs_c_E0root_x0, hbs_c_E0root_x1, &
                       hbs_c_root_accomod_vel_max, hbs_c_root_accomod_day, &
                       hbs_c_opt_root_depth, hbs_c_opt_root_thick, &
                       hbs_c_root_accomod_type, l_hbs_root_accomodation, &
                       l_hbs_root_erosion, l_hbs_bedsedstab)

      ! ***********************
      ! * OTHER INITIALIZATIONS
      ! ***********************
      hbs_hsed_prev(:, :) = hsed(:, :)

   END SUBROUTINE hbs_init
   !!===========================================================================

   SUBROUTINE hbs_readvar(iv)
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_readvar  ***
   !!
   !! ** Description : Load variables parameters
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(IN) :: iv
      !! * Local declaration
      CHARACTER(len=lchain) :: filepc
      !! hbs_var_main variables
      CHARACTER(len=lchain) :: r_hbs_varname, r_hbs_vardesc
      !! hbs_var_option variables
      LOGICAL  :: r_l_hbs_bedsedstab
      !! hbs_var_init variables
      LOGICAL               :: r_l_hbs_initfromfile
      CHARACTER(len=lchain) :: r_hbs_fn_initfile
      REAL(KIND=rsh)        :: r_hbs_i_rbiom, r_hbs_i_zroot, r_hbs_i_thickroot
      !! hbs_var_sedcoupl variables
      LOGICAL               :: r_l_hbs_root_accomodation, r_l_hbs_root_erosion
      INTEGER               :: r_hbs_c_root_accomod_type
      REAL(KIND=rsh)        :: r_hbs_c_tauceroot_x0, r_hbs_c_tauceroot_x1, &
                               r_hbs_c_E0root_x0, r_hbs_c_E0root_x1, &
                               r_hbs_c_root_accomod_vel_max, &
                               r_hbs_c_root_accomod_day, &
                               r_hbs_c_opt_root_depth, r_hbs_c_opt_root_thick
      !! * Namelists
      NAMELIST /hbs_var_main/ r_hbs_varname, r_hbs_vardesc
      NAMELIST /hbs_var_option/ r_l_hbs_bedsedstab
      NAMELIST /hbs_var_init/ r_l_hbs_initfromfile, &
         r_hbs_fn_initfile, r_hbs_i_rbiom, r_hbs_i_zroot, r_hbs_i_thickroot
      NAMELIST /hbs_var_sedcoupl/ r_l_hbs_root_accomodation, &
         r_l_hbs_root_erosion, &
         r_hbs_c_tauceroot_x0, r_hbs_c_tauceroot_x1, &
         r_hbs_c_E0root_x0, r_hbs_c_E0root_x1, &
         r_hbs_c_root_accomod_type, r_hbs_c_root_accomod_vel_max, &
         r_hbs_c_root_accomod_day, &
         r_hbs_c_opt_root_depth, r_hbs_c_opt_root_thick

      !!----------------------------------------------------------------------
      !! * Executable part

      ! save into simu.log
      !-------------------
      MPI_master_only WRITE (iscreenlog, *) &
         ' '
      MPI_master_only WRITE (iscreenlog, *) &
         '***************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***** module HYBIOSED, subroutine hbs_READVAR *****'
      MPI_master_only WRITE (iscreenlog, *) &
         '***************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         ' Reading file ', TRIM(hbs_fn_var(iv))
      MPI_master_only WRITE (iscreenlog, *) &
         ' defining hybiosed parameters'
      ! ******************************
      ! * Start reading variables file
      ! ******************************
      filepc = hbs_fn_var(iv)
      OPEN (55, file=filepc, status='old', form='formatted', &
            access='sequential')
      READ (55, hbs_var_main); REWIND (55)
      READ (55, hbs_var_option); REWIND (55)
      READ (55, hbs_var_init); REWIND (55)
      READ (55, hbs_var_sedcoupl); REWIND (55)
      CLOSE (55)
      ! ************************************
      ! * Alocate to corresponding parameter
      ! ************************************
      ! * For namelist hbs_var_main
      hbs_varname(iv) = TRIM(r_hbs_varname)
      hbs_vardesc(iv) = TRIM(r_hbs_vardesc)
      ! * For namelist hbs_var_option
      l_hbs_bedsedstab(iv) = r_l_hbs_bedsedstab
      ! * For namelist hbs_var_init
      l_hbs_initfromfile(iv) = r_l_hbs_initfromfile
      hbs_fn_initfile(iv) = TRIM(r_hbs_fn_initfile)
      hbs_i_rbiom(iv) = r_hbs_i_rbiom
      hbs_i_zroot(iv) = r_hbs_i_zroot
      hbs_i_thickroot(iv) = r_hbs_i_thickroot
      ! * For hbs_var_sedcoupl
      l_hbs_root_accomodation(iv) = r_l_hbs_root_accomodation
      l_hbs_root_erosion(iv) = r_l_hbs_root_erosion
      hbs_c_tauceroot_x0(iv) = r_hbs_c_tauceroot_x0
      hbs_c_tauceroot_x1(iv) = r_hbs_c_tauceroot_x1
      hbs_c_E0root_x0(iv) = r_hbs_c_E0root_x0
      hbs_c_E0root_x1(iv) = r_hbs_c_E0root_x1
      hbs_c_root_accomod_type(iv) = r_hbs_c_root_accomod_type
      hbs_c_root_accomod_vel_max(iv) = r_hbs_c_root_accomod_vel_max
      hbs_c_root_accomod_day(iv) = r_hbs_c_root_accomod_day
      hbs_c_opt_root_depth(iv) = r_hbs_c_opt_root_depth
      hbs_c_opt_root_thick(iv) = r_hbs_c_opt_root_thick
      !!**********************************
   END SUBROUTINE hbs_readvar
   !!===========================================================================

   SUBROUTINE hbs_alloc_nbvar
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_alloc_nbvar  ***
   !!
   !! ** Purpose : Allocation of tables depending on number of HYBIOSED
   !!              variables
   !!
   !!---------------------------------------------------------------------------
      IMPLICIT NONE

      IF (hbs_nbvar > 0) THEN
         ALLOCATE (hbs_c_tauceroot_x0(1:hbs_nbvar))
         hbs_c_tauceroot_x0(:) = 0.0_rsh
         ALLOCATE (hbs_c_tauceroot_x1(1:hbs_nbvar))
         hbs_c_tauceroot_x1(:) = 0.0_rsh
         ALLOCATE (hbs_c_E0root_x0(1:hbs_nbvar))
         hbs_c_E0root_x0(:) = 0.0_rsh
         ALLOCATE (hbs_c_E0root_x1(1:hbs_nbvar))
         hbs_c_E0root_x1(:) = 0.0_rsh
         ALLOCATE (hbs_c_root_accomod_vel_max(1:hbs_nbvar))
         hbs_c_root_accomod_vel_max(:) = 0.0_rsh
         ALLOCATE (hbs_c_root_accomod_day(1:hbs_nbvar))
         hbs_c_root_accomod_day(:) = 0.0_rsh
         ALLOCATE (hbs_c_opt_root_depth(1:hbs_nbvar))
         hbs_c_opt_root_depth(:) = 0.0_rsh
         ALLOCATE (hbs_c_opt_root_thick(1:hbs_nbvar))
         hbs_c_opt_root_thick(:) = 0.0_rsh
         ALLOCATE (hbs_c_root_accomod_type(1:hbs_nbvar))
         hbs_c_root_accomod_type(:) = 0
         ALLOCATE (l_hbs_root_accomodation(1:hbs_nbvar))
         l_hbs_root_accomodation(:) = .FALSE.
         ALLOCATE (l_hbs_root_erosion(1:hbs_nbvar))
         l_hbs_root_erosion(:) = .FALSE.
         ALLOCATE (hbs_varname(1:hbs_nbvar))
         hbs_varname(:) = "."
         ALLOCATE (hbs_vardesc(1:hbs_nbvar))
         hbs_vardesc(:) = "."
         ALLOCATE (l_hbs_bedsedstab(1:hbs_nbvar))
         l_hbs_bedsedstab(:) = .FALSE.

         ALLOCATE (l_hbs_initfromfile(1:hbs_nbvar))
         l_hbs_initfromfile(:) = .FALSE.
         ALLOCATE (hbs_fn_initfile(1:hbs_nbvar))
         hbs_fn_initfile(:) = "."
         ALLOCATE (hbs_i_rbiom(1:hbs_nbvar))
         hbs_i_rbiom(:) = 0.0_rsh
         ALLOCATE (hbs_i_zroot(1:hbs_nbvar))
         hbs_i_zroot(:) = 0.0_rsh
         ALLOCATE (hbs_i_thickroot(1:hbs_nbvar))
         hbs_i_thickroot(:) = 0.0_rsh
      END IF

   END SUBROUTINE hbs_alloc_nbvar
   !!===========================================================================

   SUBROUTINE hbs_alloc_xyz
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_alloc_nbvar  ***
   !!
   !! ** Purpose : Allocation of spatial HYBIOSED tables
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for GLOBAL_2D_ARRAY
      IMPLICIT NONE
      ALLOCATE (hbs_ws_trapp(1:hbs_kmax, GLOBAL_2D_ARRAY))
      hbs_ws_trapp(:, :, :) = 1.0_rsh
      ALLOCATE (hbs_ws_block(1:hbs_kmax, GLOBAL_2D_ARRAY))
      hbs_ws_block(:, :, :) = 1.0_rsh
      ALLOCATE (hbs_hsed_prev(GLOBAL_2D_ARRAY))
      hbs_hsed_prev(:, :) = 0.0_rsh

      IF (hbs_nbvar > 0) THEN
         ALLOCATE (hbs_zup_root0(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_zup_root0 = 0.0_rsh
         ALLOCATE (hbs_zup_root(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_zup_root = 0.0_rsh
         ALLOCATE (hbs_thick_root0(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_thick_root0 = 0.0_rsh
         ALLOCATE (hbs_thick_root(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_thick_root = 0.0_rsh
         ALLOCATE (hbs_dz_root(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_dz_root = 0.0_rsh
         ALLOCATE (hbs_dthick_root(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_dthick_root = 0.0_rsh
         ALLOCATE (hbs_position_wat(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_position_wat = 0.0_rsh
         ALLOCATE (hbs_position_bed(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_position_bed = 0.0_rsh
         ALLOCATE (hbs_root_biomass(1:hbs_nbvar, GLOBAL_2D_ARRAY))
         hbs_root_biomass = 0.0_rsh
      END IF

   END SUBROUTINE hbs_alloc_xyz
   !!===========================================================================

   SUBROUTINE hbs_write_summary
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_write_summary  ***
   !!
   !! ** Description : Write a summary of parameterization
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for mynode
      IMPLICIT NONE

      INTEGER :: iv

      MPI_master_only WRITE (iscreenlog, *) &
         ' '
      MPI_master_only WRITE (iscreenlog, *) &
         '******************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '******************** module module HYBIOSED **********************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***************** subroutine hbs_write_summary********************'
      MPI_master_only WRITE (iscreenlog, *) &
         '******************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         ' File position: ', TRIM(hbs_fn_position)
      IF (l_hbs_suspsed_trapp) THEN
         MPI_master_only WRITE (iscreenlog, *) &
            ' Suspended Particle Trapping by obstructions:   TRUE'
         MPI_master_only WRITE (iscreenlog, *) &
            ' Maximum correction factor on settling velocity due to', &
            ' particles trapping: ', hbs_c_suspsed_trapp_max
         MPI_master_only WRITE (iscreenlog, *) &
            ' Exponential coefficient for correction factor on settling', &
            ' velocity due to particles trapping: ', hbs_c_suspsed_trapp_exp
      ELSE
         MPI_master_only WRITE (iscreenlog, *) &
            ' Suspended Particle Trapping by obstructions:   FALSE'
         MPI_master_only WRITE (iscreenlog, *) &
            'Maximum correction factor on settling velocity due to particles', &
            ' trapping: NOT USED'
         MPI_master_only WRITE (iscreenlog, *) &
            ' Exponential coefficient for correction factor on settling', &
            ' velocity due to particles trapping: NOT USED'
      END IF
      IF (l_hbs_suspsed_block) THEN
         MPI_master_only WRITE (iscreenlog, *) &
            ' Suspended Particle Blockage by obstructions:   TRUE'
         MPI_master_only WRITE (iscreenlog, *) &
            ' Maximum correction factor on settling velocity due to', &
            ' particles blockage: ', hbs_c_suspsed_block_max
         MPI_master_only WRITE (iscreenlog, *) &
            ' Exponential coefficient for correction factor on settling', &
            ' velocity due to particles blockage: ', hbs_c_suspsed_block_exp
      ELSE
         MPI_master_only WRITE (iscreenlog, *) &
            ' Suspended Particle Blockage by obstructions:   FALSE'
         MPI_master_only WRITE (iscreenlog, *) &
            ' Maximum correction factor on settling velocity due to', &
            ' particles blockage: NOT USED'
         MPI_master_only WRITE (iscreenlog, *) &
            ' Exponential coefficient for correction factor on settling', &
            ' velocity due to particles blockage: NOT USED'
      END IF

      DO iv = 1, hbs_nbvar
         MPI_master_only WRITE (iscreenlog, *) &
            '***************************************************************'
         MPI_master_only WRITE (iscreenlog, *) &
            'NAMELIST : hbs_var_main'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Name (identifier) of the variable : ', TRIM(hbs_varname(iv))
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Short description of the variable : ', TRIM(hbs_vardesc(iv))
         MPI_master_only WRITE (iscreenlog, *) &
            'NAMELIST : hbs_var_option'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - If the current variable impacts bed sediment', &
            ' erosion/deposition processes : ', l_hbs_bedsedstab(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            'NAMELIST : hbs_var_init'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Filename for heterogeneous root biomass, root level ', &
            'depth and/or thickness : ', hbs_fn_initfile(iv)
         IF (.not. l_hbs_initfromfile(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial homogeneous root biomass', &
               ' : ', hbs_i_rbiom(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial homogeneous root level depth', &
               ' : ', hbs_i_zroot(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial homogeneous root level thickness', &
               ' : ', hbs_i_thickroot(iv)
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            'NAMELIST : hbs_var_sedcoupl'
         IF (l_hbs_bedsedstab(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - To activate root accomodation to erosion or deposition', &
               ' : ', l_hbs_root_accomodation(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - To activate root erosion (if .FALSE. root thickness will', &
               ' never change)   : ', l_hbs_root_erosion(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for effect of root on critical shear', &
               ' stress for erosion  : ', hbs_c_tauceroot_x0(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for effect of root on critical shear', &
               ' stress for erosion : ', hbs_c_tauceroot_x1(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for effect of roots on erosion fluxes', &
               ' : ', hbs_c_E0root_x0(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for effect of roots on erosion fluxes', &
               ' : ', hbs_c_E0root_x1(iv)
            IF (l_hbs_root_accomodation(iv)) THEN
            IF (hbs_c_root_accomod_type(iv) .EQ. 0) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Type of formulation for root accomodation velocity', &
                  ' : Constant'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Maximum (or constant) root accomodation velocity', &
                  ' (mm.month-1) : ', hbs_c_root_accomod_vel_max(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Day of the year when maximum vertical root', &
                  ' accomodation is reached : NOT USED'
            ELSEIF (hbs_c_root_accomod_type(iv) .EQ. 1) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Type of formulation for root accomodation velocity', &
                  ' : Variable (day of year)'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Maximum (or constant) root accomodation velocity', &
                  ' (mm.month-1) : ', hbs_c_root_accomod_vel_max(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Day of the year when maximum vertical root', &
                  ' accomodation is reached : ', hbs_c_root_accomod_day(iv)
            ELSEIF (hbs_c_root_accomod_type(iv) .EQ. 2) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Type of formulation for root accomodation velocity', &
                  ' : Growth model'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Maximum (or constant) root accomodation velocity', &
                  ' (mm.month-1) : ', hbs_c_root_accomod_vel_max(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Day of the year when maximum vertical root', &
                  ' accomodation is reached : NOT USED'
            END IF
            ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Type of formulation for root accomodation velocity', &
               ' : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Maximum (or constant) root accomodation velocity', &
               ' (mm.month-1) : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Day of the year when maximum vertical root', &
               ' accomodation is reached : NOT USED'
            END IF
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Optimal root level depth', &
               ' : ', hbs_c_opt_root_depth(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Optimal (maximum) root level thickness', &
               ' : ', hbs_c_opt_root_thick(iv)
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - To activate root accomodation to erosion or deposition', &
               ' : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - To activate root erosion (if .FALSE. root thickness will', &
               ' never change) : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for effect of root on critical shear', &
               ' stress for erosion : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for effect of root on critical shear', &
               ' stress for erosion : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for effect of roots on erosion fluxes', &
               ' : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for effect of roots on erosion fluxes', &
               ' : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Type of formulation for root accomodation velocity', &
               ' : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Maximum (or constant) root accomodation velocity', &
               ' (mm.month-1) : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Day of the year when maximum vertical root accomodation', &
               ' is reached       : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Optimal root level depth : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Optimal (maximum) root level thickness : NOT USED'
         END IF
      END DO
      MPI_master_only WRITE (iscreenlog, *) &
         '***************************************************************'

   END SUBROUTINE hbs_write_summary
   !!===========================================================================

   SUBROUTINE hbs_compatibility
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_compatibility  ***
   !!
   !! ** Description : Performed some compatibility tests
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for mynode
      IMPLICIT NONE

#ifndef MUSTANG
      MPI_master_only WRITE (ierrorlog, *) &
         ' '
      MPI_master_only WRITE (ierrorlog, *) &
         '****************************************************************'
      MPI_master_only WRITE (ierrorlog, *) &
         '************************ module HYBIOSED ***********************'
      MPI_master_only WRITE (ierrorlog, *) &
         '************* subroutine hbs_compatibility *********************'
      MPI_master_only WRITE (ierrorlog, *) &
         ' ERROR : HYBIOSED only available with MUSTANG'
      MPI_master_only WRITE (ierrorlog, *) &
         ' --> THE SIMULATION IS STOPPED !!!'
      STOP 1
#endif

   END SUBROUTINE hbs_compatibility
   !!===========================================================================

   SUBROUTINE hbs_vname
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_vname  ***
   !!
   !! ** Description : Performed some compatibility tests
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED
      IMPLICIT NONE

      INTEGER :: iv, indvar

      hbs_nout_pos_bed = 'pos_bed'        ! Name of in output file
      hbs_nout_pos_wat = 'pos_wat'        ! Name of in output file
      hbs_nout_susp_trapp = 'Corr_SuspTrapp' ! Name variable correction
      ! factor on settling velocity
      ! due to trapping
      hbs_nout_susp_block = 'Corr_SuspBlock' ! Name variable correction
      ! factor on settling velocity
      ! due to blockage
      hbs_nout_rbiom = 'RBiom'          ! Name variable root biomass
      hbs_nout_zroot = 'ZRoot'          ! Name variable depth of root
      ! level
      hbs_nout_thickroot = 'ThickRoot'      ! Name variable thickness of
      ! root level

      ALLOCATE (hisHbs(1:5*hbs_nbvar + 2))
      ALLOCATE (avgHbs(1:5*hbs_nbvar + 2))
      ALLOCATE (outHbs(1:5*hbs_nbvar + 2))
      ALLOCATE (out2DHbs(1:5*hbs_nbvar + 2))
      ALLOCATE (vname_hbs(20, 1:5*hbs_nbvar + 2))

      outHbs(:) = .FALSE.
      out2DHbs(:) = .FALSE.

      DO iv = 1, hbs_nbvar
         indvar = iv
         vname_hbs(1, indvar) = &
            TRIM(hbs_nout_pos_wat)//'_'//TRIM(hbs_varname(iv))
         vname_hbs(2, indvar) = &
            'Obstruction occupation rate for '//TRIM(hbs_varname(iv))
         vname_hbs(3, indvar) = '-                                    '
         vname_hbs(4, indvar) = '                                     '
         vname_hbs(5, indvar) = '                                     '
         vname_hbs(6, indvar) = 'time lat_rho lon_rho                 '
         vname_hbs(7, indvar) = '                                     '
         IF (l_hbsout_pos) THEN
            outHbs(indvar) = .TRUE.
            out2DHbs(indvar) = .TRUE.
         END IF

         indvar = hbs_nbvar + iv
         vname_hbs(1, indvar) = &
            TRIM(hbs_nout_pos_bed)//'_'//TRIM(hbs_varname(iv))
         vname_hbs(2, indvar) = &
            'Root occupation rate for '//TRIM(hbs_varname(iv))
         vname_hbs(3, indvar) = '-                                    '
         vname_hbs(4, indvar) = '                                     '
         vname_hbs(5, indvar) = '                                     '
         vname_hbs(6, indvar) = 'time lat_rho lon_rho                 '
         vname_hbs(7, indvar) = '                                     '
         IF (l_hbsout_pos) THEN
            outHbs(indvar) = .TRUE.
            out2DHbs(indvar) = .TRUE.
         END IF

         indvar = 2*hbs_nbvar + iv
         vname_hbs(1, indvar) = &
            TRIM(hbs_nout_rbiom)//'_'//TRIM(hbs_varname(iv))
         vname_hbs(2, indvar) = &
            'Root biomass for '//TRIM(hbs_varname(iv))
         vname_hbs(3, indvar) = '-                                    '
         vname_hbs(4, indvar) = '                                     '
         vname_hbs(5, indvar) = '                                     '
         vname_hbs(6, indvar) = 'time lat_rho lon_rho                 '
         vname_hbs(7, indvar) = '                                     '
         IF (l_hbsout_rbiom) THEN
            outHbs(indvar) = .TRUE.
            out2DHbs(indvar) = .TRUE.
         END IF

         indvar = 3*hbs_nbvar + iv
         vname_hbs(1, indvar) = &
            TRIM(hbs_nout_zroot)//'_'//TRIM(hbs_varname(iv))
         vname_hbs(2, indvar) = &
            'Root level depth for '//TRIM(hbs_varname(iv))
         vname_hbs(3, indvar) = 'm                                    '
         vname_hbs(4, indvar) = '                                     '
         vname_hbs(5, indvar) = '                                     '
         vname_hbs(6, indvar) = 'time lat_rho lon_rho                 '
         vname_hbs(7, indvar) = '                                     '
         IF (l_hbsout_zroot) THEN
            outHbs(indvar) = .TRUE.
            out2DHbs(indvar) = .TRUE.
         END IF

         indvar = 4*hbs_nbvar + iv
         vname_hbs(1, indvar) = &
            TRIM(hbs_nout_thickroot)//'_'//TRIM(hbs_varname(iv))
         vname_hbs(2, indvar) = &
            'Root thickness for '//TRIM(hbs_varname(iv))
         vname_hbs(3, indvar) = 'm                                    '
         vname_hbs(4, indvar) = '                                     '
         vname_hbs(5, indvar) = '                                     '
         vname_hbs(6, indvar) = 'time lat_rho lon_rho                 '
         vname_hbs(7, indvar) = '                                     '
         IF (l_hbsout_thickroot) THEN
            outHbs(indvar) = .TRUE.
            out2DHbs(indvar) = .TRUE.
         END IF
      END DO

      indvar = 5*hbs_nbvar + 1
      vname_hbs(1, indvar) = &
         TRIM(hbs_nout_susp_trapp)
      vname_hbs(2, indvar) = &
         'Correction factor for settling velocity due to trapping'
      vname_hbs(3, indvar) = '-                                    '
      vname_hbs(4, indvar) = '                                     '
      vname_hbs(5, indvar) = '                                     '
      vname_hbs(6, indvar) = 'time N  lat_rho lon_rho              '
      vname_hbs(7, indvar) = '                                     '
      IF (l_hbsout_susp_trapp) THEN
         outHbs(indvar) = .TRUE.
         out2DHbs(indvar) = .False.
      END IF

      indvar = 5*hbs_nbvar + 2
      vname_hbs(1, indvar) = &
         TRIM(hbs_nout_susp_block)
      vname_hbs(2, indvar) = &
         'Correction factor for settling velocity due to blockage'
      vname_hbs(3, indvar) = '-                                    '
      vname_hbs(4, indvar) = '                                     '
      vname_hbs(5, indvar) = '                                     '
      vname_hbs(6, indvar) = 'time N  lat_rho lon_rho              '
      vname_hbs(7, indvar) = '                                     '
      IF (l_hbsout_susp_block) THEN
         outHbs(indvar) = .TRUE.
         out2DHbs(indvar) = .False.
      END IF

   END SUBROUTINE hbs_vname
   !!===========================================================================

   SUBROUTINE hbs_nccheck(status, msg)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs_nccheck  ***
   !!
   !! ** Purpose : check netcdf function
   !!---------------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: status
      CHARACTER*(*), INTENT(IN) :: msg

      IF (status /= NF90_NOERR) THEN
         WRITE (ierrorlog, *) 'hbs_nc_check(): '
         WRITE (ierrorlog, *) msg
         WRITE (ierrorlog, *) trim(NF90_STRERROR(status))
         STOP
      END IF

   END SUBROUTINE hbs_nccheck

   !!===========================================================================

   SUBROUTINE hbs_ncget2D(ncid, varid, tmp, msg)
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs_nccheck  ***
   !!
   !! ** Purpose : get 2D from netcdf
   !!             (from nf_fread.F translate in f90)
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for GLOBAL_2D_ARRAY
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ncid, varid
      REAL, INTENT(INOUT) ::  tmp(GLOBAL_2D_ARRAY)
      CHARACTER*(*), INTENT(IN) :: msg

      INTEGER :: imin, imax, jmin, jmax, start(3), count(3)

      jmin = 0
      imin = 0
      start(1) = 1
      start(2) = 1
      start(3) = 1

#ifdef MPI
      if (ii .gt. 0) then
         start(1) = 1 - imin + iminmpi
         imin = 1
      end if
      if (ii .eq. NP_XI - 1) then
         imax = Lmmpi + 1
      else
         imax = Lmmpi
      end if
      if (jj .gt. 0) then
         start(2) = 1 - jmin + jminmpi
         jmin = 1
      end if
      if (jj .eq. NP_ETA - 1) then
         jmax = Mmmpi + 1
      else
         jmax = Mmmpi
      end if
#else
      imax = Lm + 1
      jmax = Mm + 1
#endif

      count(1) = imax - imin + 1
      count(2) = jmax - jmin + 1
      count(3) = 1

      CALL hbs_nccheck(NF90_GET_VAR(ncid, varid, tmp(imin:imax, jmin:jmax), &
                                    start, count), msg)
#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
! /* exchange needed */
      CALL exchange_r2d_tile(1, LOCALLM, 1, LOCALMM, tmp)
#endif /* MPI */

   END SUBROUTINE hbs_ncget2D
   !!===========================================================================

   SUBROUTINE hbs_init_pos
   !!---------------------------------------------------------------------------
   !!                 ***  ROUTINE hbs_init_pos  ***
   !!
   !! ** Description : read the hybiosed position file
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for GLOBAL_2D_ARRAY
      USE netcdf
      IMPLICIT NONE

      INTEGER iv, ncid, varid
      CHARACTER(30) name
      REAL tmp(GLOBAL_2D_ARRAY)

      CALL hbs_nccheck(NF90_OPEN(hbs_fn_position, NF90_NOWRITE, ncid), &
                       "hbs_init_pos nf90_open "//hbs_fn_position)

      DO iv = 1, hbs_nbvar
         name = TRIM(hbs_nout_pos_bed)//'_'//TRIM(hbs_varname(iv))
         CALL hbs_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                          "hbs_init_pos nf90_inq_varid "//name)
         CALL hbs_ncget2D(ncid, varid, tmp, &
                          "hbs_init_pos getvar "//name)
         hbs_position_bed(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

         name = TRIM(hbs_nout_pos_wat)//'_'//TRIM(hbs_varname(iv))
         CALL hbs_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                          "hbs_init_pos nf90_inq_varid "//name)
         CALL hbs_ncget2D(ncid, varid, tmp, &
                          "hbs_init_pos getvar "//name)
         hbs_position_wat(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)
      END DO

      ! Close input NetCDF file.
      CALL hbs_nccheck(NF90_CLOSE(ncid), &
                       "hbs_init_pos nf90_close "//hbs_fn_position)

      RETURN

   END SUBROUTINE hbs_init_pos
   !!===========================================================================

   SUBROUTINE hbs_init_BRoot_ZDZroot
   !!---------------------------------------------------------------------------
   !!                 *** SUBROUTINE hbs_init_BRoot_ZDZroot  ***
   !!
   !! ** Purpose : read the hybiosed spatial initialization files
   !!              or apply constant values
   !!
   !!---------------------------------------------------------------------------
      USE module_HYBIOSED ! for GLOBAL_2D_ARRAY
      USE netcdf
      IMPLICIT NONE

      INTEGER iv, ncid, varid
      CHARACTER(30) name
      REAL tmp(GLOBAL_2D_ARRAY)

      DO iv = 1, hbs_nbvar
         IF (l_hbs_initfromfile(iv)) THEN
            CALL hbs_nccheck( &
               NF90_OPEN(hbs_fn_initfile(iv), NF90_NOWRITE, ncid), &
               "hbs_init_spatial nf90_open "//hbs_fn_initfile(iv))

            name = TRIM(hbs_nout_rbiom)//'_'//TRIM(hbs_varname(iv))
            CALL hbs_nccheck( &
               NF90_INQ_VARID(ncid, name, varid), &
               "hbs_init_spatial nf90_inq_varid "//hbs_fn_initfile(iv)//" "// &
               name)
            CALL hbs_ncget2D( &
               ncid, varid, tmp, &
               "hbs_init_spatial getvar "//hbs_fn_initfile(iv)//" "//name)
            hbs_root_biomass(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            name = TRIM(hbs_nout_zroot)//'_'//TRIM(hbs_varname(iv))
            CALL hbs_nccheck( &
               NF90_INQ_VARID(ncid, name, varid), &
               "hbs_init_spatial nf90_inq_varid "//hbs_fn_initfile(iv)//" "// &
               name)
            CALL hbs_ncget2D( &
               ncid, varid, tmp, &
               "hbs_init_spatial getvar "//hbs_fn_initfile(iv)//" "//name)
            hbs_zup_root(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            name = TRIM(hbs_nout_thickroot)//'_'//TRIM(hbs_varname(iv))
            CALL hbs_nccheck( &
               NF90_INQ_VARID(ncid, name, varid), &
               "hbs_init_spatial nf90_inq_varid "//hbs_fn_initfile(iv)//" "// &
               name)
            CALL hbs_ncget2D( &
               ncid, varid, tmp, &
               "hbs_init_spatial getvar "//hbs_fn_initfile(iv)//" "//name)
            hbs_thick_root(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            ! Close input NetCDF file.
            CALL hbs_nccheck( &
               NF90_CLOSE(ncid), &
               "hbs_init_spatial nf90_close "//hbs_fn_initfile(iv))
         ELSE
            hbs_root_biomass(iv, :, :) = hbs_i_rbiom(iv)
            hbs_zup_root(iv, :, :) = hbs_i_zroot(iv)
            hbs_thick_root(iv, :, :) = hbs_i_thickroot(iv)
         END IF
      END DO
      hbs_zup_root0(:, :, :) = hbs_zup_root(:, :, :)
      hbs_thick_root0(:, :, :) = hbs_thick_root(:, :, :)

      WHERE (hbs_position_bed == 0.)
         hbs_root_biomass = 0.
         hbs_zup_root = 0.
         hbs_thick_root = 0.
      END WHERE

      RETURN

   END SUBROUTINE hbs_init_BRoot_ZDZroot

#endif /* HYBIOSED */
END MODULE init_HYBIOSED
