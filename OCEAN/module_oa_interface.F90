#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS
!------------------------------------------------------------------------------
!                               NHOMS
!                Non Hydrostatic Ocean Modeling System      
!------------------------------------------------------------------------------
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!!
!!     equipe d'oceanographie cotiere du laboratoire d'aerologie
!!     laboratoire d'aerologie
!!     cnrs - universite paul sabatier - observatoire midi pyrenees
!!     14 avenue edouard belin - 31400 toulouse - france
!
! DESCRIPTION: Module d'analyse spectrale et d'analyse en ondelettes en ligne pour Croco 
!
!> @brief OA Module : Online spectral and wavelet analysis module 
!
!> @details Procedures included in the module are:
!!
!! - init_parameter_oa       : reads the user namelist with global parameters (configuration, variable, domain, scales, frequency,..) 
!! - init_oa                 : initialization of the spatial, temporal and scale structures needed to stack/unstack analyses.
!! - rdnotebook_init_oa      : reading of the 1st namelist which defines the config-variables to analyse during the simulation.
!! - rdnotebook_oa           : reading of specific namelists each describing a perticular OA analysis to perform.
!! - init_configuration_oa   : example of initialization of "mixed" variables configuration (composite).
!! - var_grid_oa             : sets the grid point type for each defined OA variable (C grid).
!! - var_space_oa       : initializes the spatial structure of the state vector to analyse.
!! - var_per_oa         : initializes the frequency parameters of the state vector to analyse. 
!! - var_time_oa        : initializes the time parameters of the state vector to analyse.
!! - var_rec_oa         : calculates the "reconstruction factor" for Morlet wavelet, Fourier and Dirac (replaces inverse transformation).
!! - var_copy_oa        : used for "mixed" variables (composite), duplicates the variable information in another variable. 
!! - upd_init_oa        : prepares output arrays var2d_oa and var3d_oa and related parameters.
!! - main_oa            : main routine which performs the time convolution between the requested field and the OA atom (Fourier, wavelet, Dirac)
!!                        method without pointer, calls function var_oa.
!! - init_test_oa       : removed 
!! - test_oa            : OA test functions application.
!! - psi_oa             : wavelet fonction
!! - psi_p2s_oa         : function which transforms the time period to the wavelet scale (s).
!! - var_oa             : function which returns the field value at one grid point i,j,k according to the analysis code
!! - var_upd_oa         : updates output arrays var2d, var3d with the analysis when available.
!! - box_oa             : initialisation of the size of the heisenberg boxes
!! - lev_init_oa        : initilizes of the isopycne levels to analyse their positions.
!! - lev_upd_oa         : updates the wlev structured type with the curret position of the target density.
!! - update_level_oa    : searches the depth interval where the a target density is reached.
!! - subsave_oa         : calls var_upd_oa (outputs in file have been removed).
!! - allocate__*_oa     : dynamic variable allocation.
!! - deallocate_*_oa    : deallocation.
!! REMOVED :
!! - struct_oa          : REMOVED (sauvegarde de la structure spatio-temporelle de la configuration).
!! - pressure_i_oa      : REMOVED (fonction pour le calcul de l'anomalie de pression).
!! - subsave_init_oa    : REMOVED (preparation des fichiers de sauvegarde)
!! ADDED module routines
!! - count_nmsimult_oa  : pre-calculates the state vector dimension as the number of simultaneously openened convolution windows.
!! ADDED sources :
!! - output_oa          : OA updates global output array looping on tiles to finally send global array to XIOS
!! - var2d_oa, var3d_oa : called by output_oa
!! - module_parameter_oa: recycled module to allocate var3d_oa, var2d_oa to the croco GLOBAL_3D_ARRAY, 2D_ARRAY 
!! - online_spectral_diags : main routine called by Croco main.F and step.F 
!!                           supporting tile-threads and possible dual OpenMP-MPI parallelization of the horizontal domain
!
! REVISION HISTORY:
!
!> @authors 
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version (2006)
!! - B. Lemieux-Dudon
!!  - Croco-OA tile-thread compliant version (2021) 
!!    supporting tile-threads and possible dual OpenMP-MPI parallelization of the horizontal domain
!!     => online_spectral_diags : main routine called by Croco main.F and step.F
!!     => module_tile_oa : tile-compliant state vector spatial structure
!!  - Handles the analysis netCDF output with specific routines 
!!    that are tile-thread compatible
!!     => output_oa : strictly updates array tiles to finally send the global array to XIOS
!!     => var2d_oa_out, var3d_oa_out : 2D using Croco working arrays work2d, work2d2
!!     => module_parameter_oa : allocating var3d_oa, var2d_oa to the croco GLOBAL_3D_ARRAY, 2D_ARRAY 
!!  - Memory optimization with by pre-calculates the state vector dimension and reduce the allocated size
!!  - Croco-OA interface 1st version (2020)
!!  - Memory oprimization when tracking isopycnal levels (lev_init_oa,..)
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!
!> @date 2021
!> @todo BLXD Modified 2020 Croco-OnlineA module interface.
!! - The CROCO-OnlineAnalysis module interface must be modified : 
!!   Interface = mix between the Stand-alone and Original OnlineA versions 
!!   => Croco arrays are passed as arguments to the OnlineA routines with reduced 
!!   dimension range. It leads to memory duplication.
!> @todo BLXD general requirement to organize precision consistency
!!    eventhough compiling with double prec option
!!    (ifort, gfortran jobcomp -r8 + mpc.F preproc not applied to F90)
!!    double precision type is obsolete 
!!    Module with wp declaration and use of KIND(wp) ?
!> @todo remove from Makefile this useless version of module_parameter_oa
!------------------------------------------------------------------------------

      module module_interface_oa

      !use module_
      implicit none


      !! MPI subdomain dependent param
      !> sub. init_parameter_oa sets these parameters which are MPI subdomain dependent 
      integer :: mynode
      logical :: if_print_node

      !! Non MPI dependent param
      !> sub. init_parameter_oa sets the logical to .true. if user effectively requests for analysis
      logical             :: if_oa

      ! TODO SEE IF THOSE CROCO TIME PARAM DECLARED HERE ONLY
      !> First and last simulation iteration indices
      integer :: kount0, nt_max                                      

      !> External time step
      double precision :: dti

      !  BLXD BUG !!!!  TODO io_unit also declared in module_oa variables
      integer :: io_unit
 
      !> Directory in is currently hardcoded to INPUT, directory out is set in the namelist_oa 
      character(len=200)  :: directory_in_oa, directory_out_oa
      integer             :: verbose_oa
      integer, parameter  :: io_nml=90, io_hist=91
      character(len=1), parameter :: txtslash='/'

      !TODO hist_io --> history, io_unit --> alternative to stdout? 
      logical             :: if_sphe_deg2rad, if_extend_dom, if_chks ! TODO remove if_checks
      real                :: unite_oa    !< namelist parameters time units 
      real                :: pi_oa       !< TODO? double precision, real(8)

      !! CHECK IF PARAM DECLARED HERE
      integer :: nper_test
      double precision, allocatable, dimension(:) :: period_test_oa, amp_test_oa 

      ! BLXD in var_time_oa Maximum number of openned convolution windows
      ! integer :: nmsimult_oa_max

      !> Maximum number of predefined configuration/var types (TODO check public/private attributes)
      integer, parameter :: maxtyp_oa = 100 ! #BLXD set the # of possible config. maxtyp_oa   

      integer, parameter :: nmsgnv_oa = 99 ! #BLXD set the # of possible config. maxtyp_oa   

      !> Maximum number of analysis of the same type requested simultaneously
      integer, parameter :: max_idcfg_oa = 10                           
      
      !> Maximum number of mixed variables (composites) ( rm module_oa_upd )
      integer, parameter :: nmvar_oa = 5                              

      ! #BLXD brought from module_oa_space
      integer,dimension(maxtyp_oa) ::                                  &  ! changed 100 to 200 jwf 20070619
            tgv3d_oa                                                   &  !
            ,tgv_oa                                                       ! type de point de la grille c auquel est associee la variable symphonie
       
      character(len=5),dimension(maxtyp_oa) :: tgvnam_oa

      ! #BLXD brought from module_oa_variable
      integer:: nzvc_oa(maxtyp_oa)

      !> Last dimension of the allocatable array var?d_oa
      integer :: nzupd3d_oa, nzupd2d_oa
 
      !> Analysed 3D fields
      !! Dimensions are space dimensions i, j, k, and the index of the analysed variable has returned by tvar_oa
      complex(8), dimension(:,:,:,:), allocatable :: var3d_oa
      
      !> Analysed 2D fields
      !! Dimensions are space dimensions i, j, k, and the index of the analysed variable has returned by tvar_oa
      complex(8), dimension(:,:,:), allocatable   :: var2d_oa
      
      !> Returns the index of the analysed variable (required for var3d_oa and var2d_oa)
      !! dimension 1 : OA configuration-variable as requested in the namelist, tc_oa(ic), e.g., isopycne analysis has OA analysis index 20 
      !! dimension 2 : if several OA configuration-variable of the same type are requested, their rank, tvc_oa(ic) 
      !!               e.g. 1 for first OA analysis 20, 2 for the second OA analysis 20,... etc.
      !! dimension 3 : # of variables involved in the configuration (see composite config. with several variables)
      ! integer, dimension(1:maxtyp_oa,1:max_idcfg_oa,1:10)  :: tvar_oa
      integer, dimension(1:maxtyp_oa,1:max_idcfg_oa,1:nmvar_oa)  :: tvar_oa

      ! BLXD added in module_oa_variables croco2021 May 06
      !      removed from module variables declared in module_interface_oa
      !      with wrong size
      !      integer, dimension(1:6000) :: des_oa  
      !
      ! Variable useful for NHOMS energy analysis if you wish to call OA at specific point in the code
      ! for specific variable that are updated at some irregular code place or time step (e.g., time splitting)
      ! depnds on iv_m so it should have nzv_oa size has cnb_oa 

#ifdef OA_TRACES
      ! BLXD integer, dimension(1:nzv_oa) :: io_traces
#endif

      integer :: itgt_glob, jtgt_glob, ktgt


      contains

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: lecture des namelists pour initialiser les parametres
!  d'analyse (configutation, variable, frequence, echelles analysees).
!
!
!> @brief reads OA namelists and set global ananalysis parameters.
!
!> @details sets the configuration, variables, frequency of analysis
!   scales and area to analyse,... 
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - namelists, doxygen comments, stand-alone, optimization (2015)
!!  - Croco-OA interface 1st version (2020)
!!  - Croco Tile-thread compliant version (2021)
!> @date 2021
!> @todo specific module to set consitent precisions 
!!       with real(KIND=wp)... etc 
!----------------------------------------------------------------------
!
      !subroutine finalize_parameter_oa
      !
      !use module_tile_oa, only : deallocate_tile_space_str
      !implicit none 
      !return
      !end subroutine finalize_parameter_oa

!     logical function init_parameter_oa(     &
      subroutine init_parameter_oa(     &
     &  io_unit_oa                      &
     & ,if_print_node_oa                &
     & ,mynode_oa                       &
     & ,dti_oa                          & 
     & ,kount0_oa                       &
     & ,nt_max_oa                       &
     & ,ntiles )

      use module_oa_variables
      use module_oa_time
      use module_oa_space 
      use module_oa_periode
      use module_oa_stock
      use module_tile_oa, only : allocate_tile_space_str
      use module_grd_oa, only  : allocate_tile_grid_tmp_str
      !use module_isopycne_oa, only  : allocate_tile_isopycne_str

      implicit none

      !> Time integration step (s)
      double precision, intent(in) :: dti_oa

!     !> Current model integration iteration
!     integer, intent(in) :: iic_oa 

      !> First and last simulation iteration indices
      integer, intent(in) :: kount0_oa, nt_max_oa                                       

      !> Number of Croco tiles
      integer, intent(in) :: ntiles                                       

      integer, intent(in) :: io_unit_oa, mynode_oa
      logical, intent(in) :: if_print_node_oa

      integer :: ic_s

      ! Parameter initializion
      ! TODO remove from module_oa_interface 
      pi_oa = acos(-1.)

      directory_in_oa  = "INPUT"

      ! Test configuration set to .false. by default
      ifl_test_oa = 0

      ! Isopycne configuration set to .false. by default
      isopycne_analysis = .false.

      ! BLXD pointers
      io_unit = io_unit_oa
      mynode  = mynode_oa
      if_print_node = if_print_node_oa
      dti = dti_oa
      kount0 = kount0_oa 
      nt_max = nt_max_oa

      !open(unit=io_unit,file=trim(directory_out)//'/'//'trace.out')    
      if(if_print_node) write(io_unit,*) ' '
      if(if_print_node) write(io_unit,*) '..Begin Online Analysis OA LOG FILE'

!.....nombre de variables par configuration:
!     de 1 à 99 => configuration n'impliquant qu'une seule variable
!     sup a 100 => conf. a +s variables
!     BLXD TODO Specific init routine for nzvc_oa ?  
!     From 100 to parameter maxtyp_oa=200 => configuration including several variables

      do ic_s=1,nmsgnv_oa
         nzvc_oa(ic_s) = 1
      enddo
      nzvc_oa(100) = 2 

!..... BLXD removing des_oa initialization : 
!      WRONG type (integer) - WRONG array size 1:6000 instead of 1:nzv_oa
!..... No energy variable so far!
!      des_oa=0.

!.....lecture des nzc_oa configuration et calcul du nombre total de var. associe  nzv_oa
!     BLXD include composite configuration with several variables requires nzvc_oa
      call rdnotebook_init_oa

!.....test si calcul oa ou pas
      if_no_analysis_requested : if(nzv_oa.eq.0) then

        if(if_print_node) write(io_unit,*) '...OA LOG FILE : WARNING zero analysis requested (nzv_oa=0)'

!.....flag returning false to croco, either function or module flag
!        init_parameter_oa = .false.
         if_oa = .false.

        return

      else

!.....allocations dynamiques des parametres pour analyse ssi nzv_oa /= 0
        call allocate_part1_oa( )

!.....BLXD MOVING des_oa INITIALISATION HERE with correct type (integer) an array size (nzv_oa)
!     just allocated in allocate_part1
!       No energy variable so far!
        des_oa(1:nzv_oa) = 0 


!.....definition des ondelettes de morlet:

        ! BLXD appropriate place to init wavelet param, namelist can be a better place ?
        ! TODO use pi_oa or croco pi parameter
        fb_oa=2.                ! definition de l''ondelette de morlet complexe
        fc_oa=6./2./3.14159274  ! par defaut fb_oa=2. et fc_oa=6./3.14159274/2.
        fc_oa=6./2./pi_oa

!.....lecture du notebook.oa: 
!      #BLXD removing the use of 2D grid parameter. Still contains the call to var_copy_oa
        call rdnotebook_oa( )

        call var_grid_oa

!.....BLXD here we know the # of tiles ntiles from module_tile_oa
        call allocate_tile_space_str(ntiles)
        call allocate_tile_grid_tmp_str(ntiles)

        !if ( isopycne_analysis ) then
        !    call allocate_tile_isopycne_str(ntiles)
        !endif

! TODO DEALLOCATE NULLIFY AFTER init_parameter_oa terminated for all trd
! REQUIRES OMP BARRIER

!.....if function returning true : init_parameter_oa = .true.
        if_oa = .true.

      end if if_no_analysis_requested

      return 
      end subroutine init_parameter_oa
!      end function init_parameter_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Reads a first namelist which defines the number and names of all the OA analysis to 
!! simultaneously perform during the simulation.
!!
!> @details The namelist name is "namelist_oa". It provides the following parameter value:
!! - nzc_oa        \: number of simultaneous OA analysis to perform during the simulation,
!! - nzc_oa_names  \: a coma seperated list of the specific names of the requested analysis,
!! - nzc_oa_vartyp \: a coma seperated list of the specific OA configuration code corresponding to
!!   the requested analysis.
!! The number and names of OA analysis to perform are applied to read the supplementary namelists 
!! that specifically describe each type of analysis. 
!! The names of the specific namelist are constructed as follows:
!! - "namelist_" \+ nzc_oa_names 
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - namelists, doxygen comments, stand-alone version (2015)
!!  - Croco version (2020)
!> @date 2021
!> @todo BLXD 
!! - organize the io unit parameters, cleaning
!------------------------------------------------------------------------------

      subroutine rdnotebook_init_oa

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      integer             :: ic_r  !< Dummy configuration index, test function # of harmonics 
      character(len=250)  :: filin !< namelist filename

      namelist / oa_numbers / nzc_oa
      namelist / oa_names   / nzc_oa_names, nzc_oa_vartyp
      namelist / oa_test    / nper_test
      namelist / oa_test_signal / period_test_oa, amp_test_oa 
      namelist / oa_checks / if_chks, verbose_oa, directory_out_oa &
                            ,itgt_glob, jtgt_glob, ktgt
      namelist / oa_params / if_extend_dom, unite_oa, nmsimult_oa_max, if_sphe_deg2rad

      ! Namelist default values 
      ! TODO BLXD also read directory_in_oa in namelist 
      ! TODO BLXD check directory_in, io_unit in module_oa_variables

      directory_out_oa = "OUTPUT"
      if_chks    = .false.
      verbose_oa    = 6
      unite_oa      = 1
      if_extend_dom = .false.
      itgt_glob=-999 ; jtgt_glob=-999 ; ktgt=-999
      ! TODO RM nmsimult_oa_max : should't be useful anymore
      nmsimult_oa_max=500
      if_sphe_deg2rad=.false.

      ! TODO read verbose_oa, directory_out_oa in namelist, history_oa ?

      filin= trim(directory_in_oa) // txtslash // 'namelist_oa'

      open(unit=io_nml, file = filin)
      read(unit=io_nml, nml  = oa_numbers)

      oa_analysis_requested : if (nzc_oa /= 0) then
    
      call allocate_namelist_oa 

      read(unit=io_nml, nml  = oa_names)
     
      nzv_oa = 0

      do ic_r = 1 , nzc_oa
       
         if(if_print_node) write(io_unit,*) 'OA analysis names and variable types' 
         if(if_print_node) write(io_unit,*) ic_r, nzc_oa, nzc_oa_names(ic_r), nzc_oa_vartyp(ic_r) 

!........A une config ic_r peut-etre associe plusieurs variables 
!        nzvc_oa est le nombre de variables impliquees dans la config
         nzv_oa = nzv_oa + nzvc_oa( nzc_oa_vartyp(ic_r) )

        ! If one of several var-config code set to 99 => Test function
        ! A single analytical function can be tested against different atoms
        if ( (nzc_oa_vartyp(ic_r)==99) .and. (nper_test==0) ) then
          read(unit=io_nml, nml  = oa_test)
        end if

      enddo


      oa_test_requested : if (nper_test > 0) then
          
          call allocate_namelist_test_oa(nper_test) 
          read(unit=io_nml, nml  = oa_test_signal)

      else if  (nper_test == -1 ) then      
          ! BLXD nper_test = -1 hardcoded test function 
          if(if_print_node) write(io_unit,*) 'OA hardcoded test function'
          nper_test = 1 
          call allocate_namelist_test_oa(nper_test) 
          amp_test_oa(1)    = 0.D0 
          period_test_oa(1) = 0.D0

      endif oa_test_requested 

      if (nper_test > 0) then

          do ic_r = 1, nper_test
            if(if_print_node) write(io_unit,*) 'OA test function (period,amplitude) ', period_test_oa(ic_r), amp_test_oa(ic_r)
          end do

      end if

      read(unit=io_nml, nml  = oa_checks)

      read(unit=io_nml, nml  = oa_params)

      endif oa_analysis_requested 

      close(90)

      return
      end subroutine rdnotebook_init_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: allocate/deallocate namelist parameters 
!
!> @brief allocation of names and type of requested OA analysis
!
!> @details variable allocation according to nzc_oa the number of type of 
!!  requested OA analysis :
!! - nzc_oa_names  \: names of the requested analysis,
!! - nzc_oa_vartyp \: configuration code of the requested analysis.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - namelists, doxygen comments, stand-alone (2015)
!!  - Croco version (2020)
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine allocate_namelist_oa

      use module_oa_variables , only : nzc_oa_names, nzc_oa_vartyp
      use module_oa_periode , only : nzc_oa

      implicit none

          allocate( nzc_oa_names(1:nzc_oa) )
          allocate( nzc_oa_vartyp(1:nzc_oa) )

      end subroutine allocate_namelist_oa

      subroutine deallocate_namelist_oa

      use module_oa_variables , only : nzc_oa_names, nzc_oa_vartyp

      implicit none

          deallocate( nzc_oa_names )
          deallocate( nzc_oa_vartyp )

      end subroutine deallocate_namelist_oa


      subroutine allocate_namelist_test_oa(n)

      implicit none
      integer, intent(in) :: n

          allocate( period_test_oa(1:n) )
          allocate( amp_test_oa(1:n) )

      end subroutine allocate_namelist_test_oa

      subroutine deallocate_namelist_test_oa

      implicit none

          deallocate( period_test_oa, amp_test_oa )

      end subroutine deallocate_namelist_test_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: allocate/deallocate namelits parameters related to
!  a given variable-configuration analysis.
!
!> @brief allocate/deallocate namelist parameters related to
!  a given variable-configuration analysis (e.g., rho, M2)
!
!> @details parameters are array of size nzv_oa 
!!  
!
! REVISION HISTORY:
!
!> @authors
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine allocate_part1_oa( )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

!.....definition de toutes les tailles maximales des differentes structures
!     (temporelle-spatiale-frequentielle)
!
!     #BLXD old NHOMS allocation
!     allocate (rhphat_oa_t(1:imax,1:jmax,1:kmax))
!     allocate (temhat_oa_t(1:imax,1:jmax,1:kmax))
!     allocate (salhat_oa_t(1:imax,1:jmax,1:kmax))       

!.....module_var:

! BLXD bug tvc_oa must have the size of the # of configurations 
      allocate (                                                      &
         swt_d_oa(nzv_oa)                                             &   !< caracteristiques spatiales de la variable (voir notebook_oa)
        ,swt_t_oa(nzv_oa)                                             &   !< caracteristiques temporelles de la variable (voir notebook_oa)          
        ,tv_oa   (nzv_oa)                                             &   !< code variable associe 
        ,cnb_oa  (nzv_oa)                                             &   !< call number: pos. of call regarding time splitting
!        ,tvc_oa  (nzv_oa)                                             &   !< configuration associee a une variable    
        ,updv_oa (nzv_oa)                                             &   !< flag de remise a jour    
        ,save_oa (nzv_oa)                                             &   !< flag de sauvegarde
        ,tupd_oa (nzv_oa)                                             &   !< pour variables communes
        ,ltrec_oa (nzv_oa)                                            &   !< pour garder le 1er OA record de la simu
        ,if_first_rec_oa (nzv_oa)                                     &   ! ...
        !,tvar_oa (200,10,10)                                         &   !STDANDLONE known size, allocated elsewere
         )

      ! BLXD bug des_oa allocated here at the right size
       allocate ( des_oa(nzv_oa) ) 

      ! BLXD croco2021 tile-thread => spatial structure now handled in module_tile_oa.F90
      !allocate (                                                      &
      !   begvs_oa(nzv_oa+1)                                           &   ! structure 2d du vecteur d etat
      !  ,begvt_oa(nzv_oa+1)                                           &   ! structure temporelle du vecteur d etat
      !   )

      allocate (                                                      &
         begvt_oa(nzv_oa+1)                                           &   ! structure temporelle du vecteur d etat
         )

      allocate (                                                      &
         lat_oa  (2,nzv_oa)                                           &   ! latitude  min, max de la structure 2d du vecteur d etat
        ,lon_oa  (2,nzv_oa)                                           &   ! longitude  ...
        ,h_oa    (2,nzv_oa)                                           &   ! profondeur ...  
        ,k_oa    (2,nzv_oa)                                           &   ! niveaux verticaux min et max de la structure 3d du vecteur d etat
        ,ptij_oa (2,nzv_oa)                                           &   ! point particulier demande par l utilisateur
         )
  
      allocate (                                                      &
         dx_oa   (nzv_oa)                                             &   ! resolution horizontale suivant x demandee par l utilisateur
        ,dy_oa   (nzv_oa)                                             &   ! resolution horizontale suivant y demandee par l utilisateur
        ,dk_oa   (nzv_oa)                                             &   ! resolution verticale   suivant z demandee par l utilisateur
         )
 
      allocate(                                                       &
         nzpt_per_oa(nzv_oa)                                          &   ! nombre de points de discretisation par periode
         )

      allocate (                                                      &
         kount_user_oa(3,nzv_oa)                                      &   ! description des periodes choisies par l utilisateur
         )

      allocate (                                                      &
         t0_oa (nzv_oa)                                               &   ! date de la premiere sortie
                 )
 
      allocate (                                                      &
         tpsi_oa (nzv_oa)                                             &   ! type d atome utilise
         )

      allocate (                                                      &
         swt_wfpf_oa (nzv_oa)                                        &    ! calcul du coef wf ou du spectre pf (choix utilisateur)
        ,fl_rec_oa   (nzv_oa)                                         &   ! flag de reconstruction
         )

      allocate (                                                      &
         dori_oa      (nzv_oa)                                        &   ! configuration frequentielle choisie par l utilisateur
        ,delta_t_oa   (nzv_oa)                                        &   ! nombre de periodes etudiees
         )

      allocate (                                                      &
         per_oa (3,nzv_oa)                                            &   ! preriodes min, max, delta de chaque variable (ou configuration)
         )

      allocate (                                                      &
         begvp_oa (nzv_oa+1)                                          &   ! structure frequentielle du vecteur d etat
         )

      allocate (                                                      &   ! Configuration type index
         tc_oa (nzc_oa)                                               &
        ,tvc_oa(nzc_oa)                                               &   ! If several config. of the type is requested
         )   

      allocate (                                                      &
         begc_oa (nzc_oa+1)                                           &
         )   

      return
      end subroutine allocate_part1_oa

      subroutine deallocate_part1_oa( )
      implicit none
      !deallocate()  
      return
      end subroutine deallocate_part1_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Read each specific namelist parameters associated to each requested OA analysis.
!!
!> @details The namelist parameters describe the configuration-variable online analysis, 
!! with especially the horizontal and vertical part of the domain to analyse, 
!! at which simulation time steps (e.g., defines the convolution window time characteristic,...), 
!! at which time scales (e.g., wavelet scale, fourier period).
!! The specific namelist names are constructed on the basis of parameter nzc_oa_names : 
!! "namelist_" \+ nzc_oa_names \+ "_oa"
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - namelists, doxygen comments, stand-alone version (2015)
!!  - add new croco sub. tool_datetosec
!!  - Croco version (2020)
!> @date 2021
!> @todo BLXD 
!! - test new croco sub. tool_datetosec
!! - organize the io unit parameters, cleaning
!------------------------------------------------------------------------------

      subroutine rdnotebook_oa( )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock
!     use module_oa_upd
!     use module_nrj

      implicit none


      !> Local Namelist parameters
      integer :: TPSI_W,     &  !< Type of atom (0:Dirac, 1:Ondelette de Morlet, 2:Windowed Fourier, 3: Fourier)
                 SAVE_W,     &  !< Not applied (0,1:variable utilisateur,2:variable commune)
                 UPDV_W,     &  !< Variable remise à jour interactivement? Keep it to 2.
                 SWT_WFPF_W, &  !< Outputs: real part Wf (1), Energy Pf=|Wf|^2 (2) Absolute value |Wf| (4) Full complex analysis Wf
                 FL_REC_W,   &  !< Reconstruction flag 0: raw coefficient, 1: reconstruction
                 DORI_W,     &  !< Frequency analysis configuration: Discrete (1), Integration (2)
                 NZPT_PER_W, &  !< Number of points per period of analysis
                 CNB_W,      &  !< call number: position of the call in the baroclinic / barotropic time step. Not applied. Keep it to zero.
                 SWT_T_W,    &  !< 1: between two dates [T1:T2] , 2: [T1,end] , 3: single date T1, 4: all the simulation (#BLXD ONLY 4 TESTED)
                 SWT_D_W,    &  !< Domaine spatial (1: région [lat1,lat2],[lon1,lon2], 2: point (I,J), 3: Région [Hmin,Hmax] dans domaine (Lat,Lon)
                 DX_W,       &  ! Résolution en X (1 point sur...)
                 DY_W,       &  ! Résolution en Y (1 point sur...)
                 DK_W           ! Résolution en Z (1 point sur...)


      integer , dimension(1:2) :: K_W    !< Niveaux verticaux Min et Max K1,K2  (-99 pour colonne entiere)

      !> TODO test the routine to convert date YYYY,MM,DD,HH,MM,SS to the simulation time step index
      integer, dimension(1:6) :: DATE_DEB, & !< Annee, Mois, Jour, Heure, Min., Sec. --> ex : 2011,01,1,00,01,18
                                 DATE_FIN    !< Annee, Mois, Jour, Heure, Min., Sec. --> ex : 2011,01,1,00,15,15

      real  :: DT_W,  &         !< Temporal period in second at wich analysis are repeated (if SWT_T_W=1 ou 4)
               T0_W,  &         !< Time lag to launch the first analysis (SWT_T_W=4)
               DELTA_T_W        !< 1/2 Largeur de l'atome ou nombre de périodes par atome (1 for Dirac, 2 or more for Wavelet)


      real, dimension(1:2)  :: LAT_W,  & !< Latitude Min et Max
                               LON_W,  & !< Longitude Min et Max
                               HH_W,   & !< Profondeur Min et Max
                               PTIJ_W    !< Point particulier demande par l utilisateur

      real, dimension(1:3)  :: PER_W     !< Minimum, Maximum time period (s) (wavelets & fourier windows) or Length (Dirac brush), 0 for delta Dirac
                                         !! step determining # of periods to be examined (-99=optimum)

      !> Local 
      character(len=250)  :: filin !< namelist filename

      integer ::                                                      &
           iv_r                                                       &
          ,ip_r                                                       &
          ,ic_r                                                      

      real    :: dt_r

      double precision    :: time_from_croco_tref_in_sec
      integer             :: izv_oa

      namelist /oa_parameters/ TPSI_W, SAVE_W, UPDV_W, SWT_WFPF_W, FL_REC_W, PER_W,    &
                               DORI_W, DELTA_T_W, NZPT_PER_W, CNB_W, SWT_T_W, SWT_D_W, &
                               DATE_DEB, DATE_FIN, DT_W, T0_W, LAT_W, LON_W,           &
                               HH_W, K_W, DX_W, DY_W, DK_W, PTIJ_W
   
      oa_analysis_requested : if (nzc_oa /=0) then
    
!..... add unite_oa to namelist i) seconds unite_oa=1, ii) hours unite_oa=3600
      
      izv_oa = 0 ! Must be equal to nzv_oa at the end of the oa_config loop

      if(if_print_node) write(io_unit,*) '...Loop on Online Analysis configurations'
      oa_config_loop : do ic_r = 1 , nzc_oa
       
         if(if_print_node) write(io_unit,*) '...Configuration index ', ic_r

!..... namelist lue:
 
         filin = trim(directory_in_oa) // txtslash // 'namelist_' // trim( nzc_oa_names(ic_r) ) // '_oa'

         open(unit=io_nml, file=filin)
         read(unit=io_nml, nml  = oa_parameters)

!..... HYP : config. ic_r has a single variable ie, nzvc_oa( nzc_oa_vartyp(ic_r) ) = 1
         izv_oa = izv_oa + 1

!......lecture d'une configuration:
!
         tv_oa(izv_oa) = nzc_oa_vartyp(ic_r)        ! type of variable (of config. izv_oa)
         if(if_print_node) write(io_unit,*) '...Variable type ', tv_oa(izv_oa)
  
!......mise a jour de la config.:

         tc_oa (ic_r)  = tv_oa (izv_oa)             ! type of variable(of config izv_oa)
         begc_oa(ic_r) = izv_oa                     ! config pointer : several variables per config

         tpsi_oa (izv_oa)    = TPSI_W          ! type of atom
         save_oa (izv_oa)    = SAVE_W      
         updv_oa (izv_oa)    = UPDV_W      
         swt_wfpf_oa(izv_oa) = SWT_WFPF_W
         fl_rec_oa(izv_oa)   = FL_REC_W
         per_oa (1,izv_oa)   = PER_W(1)       ! minimum period (h)
         per_oa (2,izv_oa)   = PER_W(2)       ! maximum period (h)
         per_oa (3,izv_oa)   = PER_W(3)       ! step determining # of periods to be examined : evenly spaced or optimal
         dori_oa(izv_oa)     = DORI_W
         delta_t_oa (izv_oa) = DELTA_T_W      ! # periods per atom
         nzpt_per_oa(izv_oa) = NZPT_PER_W     ! # calculation points per period
         cnb_oa(izv_oa)      = CNB_W          ! call number which may be used to define where a variable needs to be called.
!        pour les bilans energetiques cette variable est mise a jour automatiquement un peu plus loin.

         if ( per_oa(1,izv_oa).gt. per_oa(2,izv_oa) ) then
            if(if_print_node) write (io_unit,*) " ERROR : per_oa (1,",izv_oa,") trop grand"
            stop
         endif  

         if (per_oa (2,izv_oa).ne.per_oa (1,izv_oa).and.                        &
              per_oa (3,izv_oa).gt.(per_oa (2,izv_oa)-per_oa (1,izv_oa))        &
              .and.per_oa (3,izv_oa).ne.-99.) then                                                         
            if(if_print_node) write (io_unit,*) " ERROR : per_oa (3,",izv_oa,") trop grand"
            stop
         endif  
         
         do ip_r = 1 , 3
            per_oa  (ip_r,izv_oa) = per_oa(ip_r,izv_oa) * unite_oa 
         enddo

         swt_t_oa(izv_oa) = SWT_T_W
         swt_d_oa(izv_oa) = SWT_D_W

! #BLXD 2020 dev new croco sub. tool_datetosec
!            convert date YYYY, MM, DD, HH, MM, SS to the model time step index
!            TODO test with swt_t_oa=4 (even though not used ) then remove

         call tool_datetosec( DATE_DEB(1), DATE_DEB(2), DATE_DEB(3), &
          &                   DATE_DEB(4), DATE_DEB(5), DATE_DEB(6), &
          &                   time_from_croco_tref_in_sec)
         kount_user_oa(1,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         
         call tool_datetosec( DATE_FIN(1), DATE_FIN(2), DATE_FIN(3), &
          &                   DATE_FIN(4), DATE_FIN(5), DATE_FIN(6), &
          &                   time_from_croco_tref_in_sec)
         kount_user_oa(2,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         
!...>    (1) [T1,T2], (2) [T1,end], (3) single date T1             
         if  ( (swt_t_oa(izv_oa) == 1) .or. ( swt_t_oa(izv_oa) == 2)  .or. ( swt_t_oa(izv_oa) == 3) ) then 
         call tool_datetosec( DATE_DEB(1), DATE_DEB(2), DATE_DEB(3), &
          &                   DATE_DEB(4), DATE_DEB(5), DATE_DEB(6), &
          &                   time_from_croco_tref_in_sec)

             kount_user_oa(1,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         end if
         if  (swt_t_oa(izv_oa) == 1) then 
         call tool_datetosec( DATE_FIN(1), DATE_FIN(2), DATE_FIN(3), &
          &                   DATE_FIN(4), DATE_FIN(5), DATE_FIN(6), &
          &                   time_from_croco_tref_in_sec)

             kount_user_oa(2,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         end if
         ! #BLXD OnlineA 2020 change : nt_max-1 -> nt_max
         ! For wavelet and Fourier ( see var_time_oa )
         ! The nzvt_oa convolution windows are defined over the model time index 
         ! intervals [ kountv_oa(1,:), kountv_oa(2,:) ]
         ! Being given the width of the convolution window 2*dkount_tot_t
         ! and its time index center at k and being given the simulation time window :
         ! - the 1st possible convolution window can start at the 1st model time index
         ! which is the 'now' time index of the restart just before the 1st 
         ! model time step from now to after, ie kountv_oa(1, izvt_oa=1 ) = kount0
         ! - the last possible convolution window can have its last window time index
         ! at the last model time index, which is to say at 'after' state of the last
         ! - the last time index of the last possible convolution window si allowed to 
         !   match the 'after' index of the last time step of the model,
         !   kountv_oa( 2, izvt_oa=nzvt_oa ) = nt_max
         ! #BLXD TODO check with Dirac atom 
         !     a) if perv_oa = 0. => dkount_tot_t = 0 (OK) 
         !     b) else               dkount_tot_t =  INT(perv_oa / 2 / dti - 0.5) + 1

         !if  ( swt_t_oa(izv_oa) == 2) kount_user_oa(2,izv_oa) = nt_max-1
         if  ( swt_t_oa(izv_oa) == 2) kount_user_oa(2,izv_oa) = nt_max
         if  ( swt_t_oa(izv_oa) == 3) kount_user_oa(2,izv_oa) = kount_user_oa(1,izv_oa)
       
         if (swt_t_oa(izv_oa).eq.1.and.kount_user_oa(2,izv_oa).lt.kount_user_oa(1,izv_oa))   then
              if(if_print_node) write (io_unit,*) "probleme dates choisies t1 > t2" 
         end if
         !if (swt_t_oa(izv_oa).eq.1.and.kount_user_oa(1,izv_oa).gt.nt_max-1)                  then 
         if (swt_t_oa(izv_oa).eq.1.and.kount_user_oa(1,izv_oa).gt.nt_max)                  then 
              if(if_print_node) write (io_unit,*) "probleme dates choisies t1 > t final"
         end if
         if ( swt_t_oa(izv_oa) == 4) then
            kount_user_oa(1,izv_oa) = kount0
            ! BLD 2020 change kount_user_oa(2,izv_oa) = nt_max-1
            kount_user_oa(2,izv_oa) = nt_max
         endif

         if(if_print_node) write(io_unit,*) '=> SWT_T_OA option is ',swt_t_oa(izv_oa)
         if(if_print_node) write(io_unit,*) '   OA analysis over Time steps ',kount_user_oa(1,izv_oa),kount_user_oa(2,izv_oa)
         if(if_print_node) write(io_unit,*) '   Simu Time step interval including initial and last',kount0,nt_max
         if  ( (swt_t_oa(izv_oa) == 1) .or. ( swt_t_oa(izv_oa) == 2)  .or. ( swt_t_oa(izv_oa) == 3) ) then 
             if(if_print_node) write(io_unit,*) 'ERRROR'
             if(if_print_node) write(io_unit,*) 'This option is not allowed (developments required), please set SWT_T_W to 4 in the OA namelist'  
             stop
         end if                     
        
         dt_r = DT_W
         kount_user_oa(3,izv_oa) = max(1,int ( dt_r * unite_oa / dti ))

         if(if_print_node) write(io_unit,*) '   OA analysis time discretization ',kount_user_oa(3,izv_oa)
        
         t0_oa(izv_oa) = T0_W
         t0_oa(izv_oa) = t0_oa(izv_oa) * unite_oa
         kount_user_oa(1,izv_oa) = kount_user_oa(1,izv_oa) + int(t0_oa(izv_oa) / dti) 

         if(if_print_node) write(io_unit,*) '   OA analysis starting earlier with T0_W (sec) ',t0_oa(izv_oa)
         if(if_print_node) write(io_unit,*) '   OA analysis over Time steps ',kount_user_oa(1,izv_oa),kount_user_oa(2,izv_oa)

         lat_oa(1,izv_oa) = LAT_W(1)
         lat_oa(2,izv_oa) = LAT_W(2)
         lon_oa(1,izv_oa) = LON_W(1)
         lon_oa(2,izv_oa) = LON_W(2)
         h_oa(1,izv_oa) = HH_W(1)
         h_oa(2,izv_oa) = HH_W(2)

         k_oa(1,izv_oa) = K_W(1)
         k_oa(2,izv_oa) = K_W(2)

!.....TODO IMPT Move this in extend_domain
! 
!         if (k_oa(1,izv_oa).eq.-99) then
!            k_oa(1,izv_oa) = 1 
!         endif
!         if (k_oa(2,izv_oa).eq.-99) then
!            k_oa(2,izv_oa) = kmax
!         endif
        
         dx_oa(izv_oa) = DX_W
         dy_oa(izv_oa) = DY_W
         dk_oa(izv_oa) = DK_W
         ptij_oa(1,izv_oa) = PTIJ_W(1)
         ptij_oa(2,izv_oa) = PTIJ_W(2)

!  attention, test uniquement sur les points z et aps sur les points,x,y,rot
!  BLXD unit conv can be moved in init_oa
#ifdef SPHERICAL
         ! #BLXD 2020 change
         if ( if_sphe_deg2rad ) then
         if(if_print_node) write(io_unit,*) 'SPHERICAL (nml) : degrees converted to rads'

         lat_oa(1,izv_oa)=LAT_W(1)*pi_oa/180.D0
         lat_oa(2,izv_oa)=LAT_W(2)*pi_oa/180.D0     
         lon_oa(1,izv_oa)=LON_W(1)*pi_oa/180.D0
         lon_oa(2,izv_oa)=LON_W(2)*pi_oa/180.D0
         else 
         if(if_print_node) write(io_unit,*) 'SPHERICAL (nml) : degrees NOT converted to rads'
         endif
#else    
         
         if(if_print_node) write(io_unit,*) 'NOT SPHERICAL (nml) : coordinates xr, xp are in meters'
         
         lat_oa(1,izv_oa) = LAT_W(1)
         lat_oa(2,izv_oa) = LAT_W(2)
         lon_oa(1,izv_oa) = LON_W(1)
         lon_oa(2,izv_oa) = LON_W(2)
#endif

!.....NHOMS configuration treatment and nrj eliminated

!.....extend_domain_oa eliminated

!.....treatment of specific variable configuration

!.....TODO ? Move this in extend_domain. Not sure
 
         if ( tv_oa (izv_oa)==20 ) then
            isopycne_analysis = .true.
            if(if_print_node) write (io_unit,*) "....Isopycne analysis requested"
         else if ( (tv_oa (izv_oa)==99) .or. (tv_oa (izv_oa)==98) ) then
            ifl_test_oa = 1
            if(if_print_node) write (io_unit,*) "....Test analysis requested"
         endif

!.....treatment of "test" variable configuration
!     BLXD REMOVED SUBROUTINE    call init_test_oa(izv_oa)
 
!.....treatment of "mixed" variable configuration
!.....HYP : config. ic_r has 1 variable ie, nzvc_oa( nzc_oa_vartyp(ic_r) ) = 1

         if ( nzvc_oa( nzc_oa_vartyp(ic_r) ) >= 2 ) then
            if(if_print_node) write(io_unit,*) 'Composite configuration : init_configuration_oa'
            call init_configuration_oa(izv_oa)
         endif


      enddo oa_config_loop

      if ( izv_oa /= nzv_oa ) then
        if(if_print_node) write (io_unit,*) "ERROR in rdnotebook_oa : izv_oa should be equal to nzv_oa"
        stop
      endif

!.....desallocation tableau de lecture des namelists OA

      ! #BLXD DO NOT DEALLOCATE AND MAKE PUBLIC 
      ! - nzc_oa_names (config name in namelist) 
      ! - nzc_oa_vartyp (possibly composite config-variable name)
      ! TODO useful to construct OA analysis variable field_def.xml ?
      ! write( vnam_oa, fmt='(a5,5,i3.3,a1,i3.3)') tgvnam_oa(tv_oa(iv)),'3d_r_',nzc_oa_names(ic),'_',tgvnam_oa(nzc_oa_vartyp(tv_oa(iv)))

      call deallocate_namelist_oa

!.....fin de la derniere configuration:

      begc_oa(nzc_oa + 1) = nzv_oa + 1

      close(90)

!.....quelques controles:

      do ic_r=1,nzc_oa
         do iv_r=begc_oa(ic_r),begc_oa(ic_r+1)-1
            ! Fourier => analysis on the full simulation
            ! #BLXD why not choosing swt_t_oa(iv_r) = 1 or 2 ????
            if (tpsi_oa(iv_r).eq.3) then 
                if(if_print_node) write (io_unit,*) "WARNING : Fourier requested => changing SWT_T_OA to 4"
                swt_t_oa(iv_r) = 4
            endif
         enddo
      enddo

      endif oa_analysis_requested

      return

      end subroutine rdnotebook_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine must be implemented to apply new configuration with mixed variables,
!! and/or when you wish to call main_oa at different location of the calling code (e.g., 
!! time-splitting case where an updated variable is available at different location).
!
! DESCRIPTION: Example pour la construction d'analyse de variables composites. 
!
!> @brief Handles the configuration code corresponding to "mixed" variables. 
!
!> @details Mixed variables are encounterd for configuration that requires several variables of the calling code
!! to construct a "mixed variable" on the basis of "simple" variables, 
!! e.g., kinetic energy calculated from zonal and meridional velocities. 
!! 
!! Use example configuration 100 to construct new configuration with mixed variables 
!! - Configuration code from 1 to 99 corresponds to "simple" variable.
!! - Configuration code from 100 to XX corresponds to "mixed" variables. 
!! Parameter tc_oa stores the configuration codes as requested in the OA namelist. 
!!
!! Default cnb_oa parameter is set to zero for all variables. 
!! cnb_oa can be used to call main_oa at different location of the calling code (e.g., 
!! time-splitting case where an updated variable is available at different location).
!! In such case, change the cnb_oa parameter for particular variable, and adapt main_oa accodingly.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco version (2020-2021), Stand alone version (2015)
!> @date 2021
!> @todo BLXD var_oa_copy moved to module_oa_variable ? 
!------------------------------------------------------------------------------
      subroutine init_configuration_oa(izv_oa)

      use module_oa_variables, only : tv_oa

      implicit none

      integer, intent(inout) :: izv_oa ! counts for the final nzv_oa
      integer :: tmp_izv_oa

!.....configuration a multiples variables ie, nzvc_oa( nzc_oa_vartyp(ic_r) ) > 1

!-------------------------------------------
!     Example of mixed variable configuration (100)
!-------------------------------------------
!......> Composite variable with type index 100
!        Adding the variables, copying the config. and incrementing the total number of variable izv_oa  
         if (tv_oa(izv_oa).eq.100) then

            tmp_izv_oa = izv_oa

!------>1ere variable: velocity vel_u 
            tv_oa(izv_oa) = 1
            izv_oa = izv_oa + 1

!     copie toutes les infos pour la nouvelle variable:
            !call var_copy_oa ( izv_oa - 1 , izv_oa )
            call var_copy_oa ( tmp_izv_oa , izv_oa )

!------>2eme variable: velocity vel_v
            tv_oa(izv_oa) = 2
            izv_oa = izv_oa + 1

!     copie toutes les infos pour la nouvelle variable:
            !call var_copy_oa ( izv_oa - 2 , izv_oa )
            call var_copy_oa ( tmp_izv_oa , izv_oa )

!------>3eme variable: velocity vel_w
            tv_oa(izv_oa) = 3
            izv_oa = izv_oa + 1

!     copie toutes les infos pour la nouvelle variable:
            !call var_copy_oa ( izv_oa - 2 , izv_oa )
            call var_copy_oa ( tmp_izv_oa , izv_oa )

         else if (tv_oa(izv_oa)>100) then
            if(if_print_node) write (io_unit,*) "ERROR: configuration above 100 is not yet hardcoded"
            stop
         endif

! Uncomment these lines and adapt them if you wish to call the OA treatment at specific
! point of the code (cf time splitting), for specific variables that are updated
! at some irregular time step or code place. 
!-------------------------------------------
!   Eventuellement point de sortie pour les variables
!      associees aux bilans energetiques
!-------------------------------------------
         ! if (tv_oa(izv_oa).eq.30) then
         !    outw_oa(5,2) =1
         !    outw_oa(5,5) =1
         !    cnb_oa(izv_oa)=-30
         !    des_oa(izv_oa)=-5005
         ! endif

      end subroutine init_configuration_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Handles "mixed" variables configuration. 
!
!> @details Configuration codes for mixed variables involves several variables,
!  and requires to copy out configuration informations to the variable parameters. 
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco version (2020), stand alone version (2015)
!> @date 2021
!> @todo BLXD declare var_oa_copy to module_oa_variable ? 
!------------------------------------------------------------------------------

      subroutine var_copy_oa (                                       &
                  nzvold_c                                           &             
                 ,nzvnew_c   )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none
      integer                                                         &
        nzvold_c                                                      &
       ,nzvnew_c

        tv_oa  (nzvnew_c)         = tv_oa  (nzvold_c)
        t0_oa  (nzvnew_c)         = t0_oa  (nzvold_c)
        tpsi_oa(nzvnew_c)         = tpsi_oa(nzvold_c)      
        save_oa(nzvnew_c)         = save_oa(nzvold_c)
        updv_oa(nzvnew_c)         = updv_oa(nzvold_c)
        swt_wfpf_oa(nzvnew_c)     = swt_wfpf_oa(nzvold_c)
        fl_rec_oa(nzvnew_c)       = fl_rec_oa(nzvold_c)
        per_oa (1,nzvnew_c)       = per_oa (1,nzvold_c)
        per_oa (2,nzvnew_c)       = per_oa (2,nzvold_c)
        per_oa (3,nzvnew_c)       = per_oa (3,nzvold_c) 
        delta_t_oa(nzvnew_c)      = delta_t_oa(nzvold_c)
        nzpt_per_oa(nzvnew_c)     = nzpt_per_oa(nzvold_c)
        cnb_oa(nzvnew_c)          = cnb_oa(nzvold_c)
        dori_oa(nzvnew_c)         = dori_oa(nzvold_c)
        ! BLXD croco2021 added very IMPT
        des_oa(nzvnew_c)          = des_oa(nzvold_c)
        ! BLXD croco2020 
        ! added but not formally needed since not init. yet (var_time_oa)
        ! useful to check first conv. wind. for each var/cfg
        ltrec_oa       (nzvnew_c) = ltrec_oa       (nzvold_c)
        if_first_rec_oa(nzvnew_c) = if_first_rec_oa(nzvold_c)
        swt_t_oa(nzvnew_c)        = swt_t_oa(nzvold_c)
        swt_d_oa(nzvnew_c)        = swt_d_oa(nzvold_c)
        kount_user_oa(1,nzvnew_c) = kount_user_oa(1,nzvold_c)
        kount_user_oa(2,nzvnew_c) = kount_user_oa(2,nzvold_c)
        kount_user_oa(3,nzvnew_c) = kount_user_oa(3,nzvold_c)
        lat_oa(1,nzvnew_c)        = lat_oa(1,nzvold_c)
        lat_oa(2,nzvnew_c)        = lat_oa(2,nzvold_c)
        lon_oa(1,nzvnew_c)        = lon_oa(1,nzvold_c)
        lon_oa(2,nzvnew_c)        = lon_oa(2,nzvold_c)
        h_oa(1,nzvnew_c)          = h_oa(1,nzvold_c) 
        h_oa(2,nzvnew_c)          = h_oa(2,nzvold_c)
        k_oa(1,nzvnew_c)          = k_oa(1,nzvold_c)
        k_oa(2,nzvnew_c)          = k_oa(2,nzvold_c)
        dx_oa(nzvnew_c)           = dx_oa(nzvold_c)
        dy_oa(nzvnew_c)           = dy_oa(nzvold_c)
        dk_oa(nzvnew_c)           = dk_oa(nzvold_c)
        ptij_oa(1,nzvnew_c)       = ptij_oa(1,nzvold_c)
        ptij_oa(2,nzvnew_c)       = ptij_oa(2,nzvold_c)

       return
       end subroutine var_copy_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: a completer en fonction des variables a traiter  
!
!> @brief Specifies the number of dimensions of the variable to analyse, 
!! and at which grid point fields are defined.
!
!> @details This routine must be modified when adding a new variable:
!! - tgv3d_oa must be set to 3 or 2 for three and two-dimensional variable respectively.
!! - tgv_oa must be set to:
!!   - "1" for t-grid point,
!!   - "2" for u-grid point,
!!   - "3" for v-grid point,
!!   - "4" for f-grid point.
! 
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco version (2020)
!> @date 2021
!> @todo use croco grid definition ?
!------------------------------------------------------------------------------
      subroutine var_grid_oa
    
      !use module_oa_space , only : tgv3d_oa, tgv_oa, tgvnam_oa 

      implicit none

!.....vitesse vel_u:
      tgv_oa(1)   = 2 
      tgv3d_oa(1) = 3
      tgvnam_oa(1) = 'del_u'

!.....vitesse vel_v:
      tgv_oa(2)   = 3 
      tgv3d_oa(2) = 3
      tgvnam_oa(2) = 'del_v'

!     BLXD croco should avoid k=0  for w-grid point
!     2D rmask is appropriate for w-grid point
!.....vitesse vel_w:
      tgv_oa(3)   = 1 
      tgv3d_oa(3) = 3
      tgvnam_oa(3) = 'vel_w'

!.....vitesse vel_u:
      tgv_oa(4)   = 2 
      tgv3d_oa(4) = 3
      tgvnam_oa(4) = 'vel_u'

!.....vitesse vel_v:
      tgv_oa(5)   = 3 
      tgv3d_oa(5) = 3
      tgvnam_oa(5) = 'vel_v'

!.....vitesse u_bar:
      tgv_oa(6)   = 2 
      tgv3d_oa(6) = 2
      tgvnam_oa(6) = 'ubar_'

!.....vitesse v_bar:
      tgv_oa(7)   = 3 
      tgv3d_oa(7) = 2
      tgvnam_oa(7) = 'vbar_'

!.....variable temp:
      tgv_oa(8)   = 1 
      tgv3d_oa(8) = 3
      tgvnam_oa(8) = 'temp_'

!.....variable salt:
      tgv_oa(9)   = 1 
      tgv3d_oa(9) = 3
      tgvnam_oa(9) = 'salt_'

!.....variable density:
      tgv_oa(11)   = 1 
      tgv3d_oa(11) = 3
      tgvnam_oa(11) = 'rho__'

!.....variable rhp_t to track isopycne movement:
      tgv_oa(20)   = 1 
      tgv3d_oa(20) = 3
      tgvnam_oa(20) = 'isplv'

!.....test variable vardp_test_oa mais 2D:
      tgv_oa(98)   = 1 
      tgv3d_oa(98) = 2
      tgvnam_oa(98) = 'test_'

!.....test variable vardp_test_oa:
      tgv_oa(99)   = 1 
      tgv3d_oa(99) = 3
      tgvnam_oa(99) = 'test_'

!.....horizontal kinetic energy at t-grid point:
      tgv_oa(100)   = 1 
      tgv3d_oa(100) = 3
      tgvnam_oa(100) = 'comp_'

      return
      end subroutine var_grid_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION:  initilisation des structures spatiales, frequentielles et 
!!  temporelles du vecteur d'etat.
!
!> @brief construction of the spatial, temporal and scale structures needed
!!  to stack and unstack the analysis in the state vector. 
!!     
!!
!> @details several array-pointer allocation to vectorize the analysis.
!!   The state vector is an derived data type array whose dimension
!!   is the # of convolution windows (time structure) and whose field
!!   is the # of grid point to analyse (spatial structure). 
!!   Spatial structures are croco tile-thread dependent.
!!  
!! 1) Croco tile-thread version : OA module variables which are related to horizontal
!!    subdomain are now embeded in a derived type array having the dimension of the
!!    number of tiles. This routine is called within a tile loop in main.F.
!!    see, init_online_spectral_diags(Istr,Iend,Jstr,Jend,tile) 
!! 2) OA stand-alone/croco version WARNING : 
!!    Memory duplication when passing arguments to the init_oa subroutine :
!!    - the calling program should set the exact array dimension ranges 
!!      using lbound and ubound fonctions.   
!!    - for croco tile-thread version => use pointers (applied for grid arrays only)
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco Tile-thread compliant version (2021) => module_tile_oa
!!  - temporary use of pointer to croco grid arrays => module_grd_oa
!!    in order to build the spatial structure of the OA state vector.
!!  - Croco-OA interface 1st version (2020)
!!  - namelists, doxygen comments, Stand-alone version, optimization (2015)
!! @todo reintroduce pointers for isopycne tracking and/or function var_oa ?
!------------------------------------------------------------------------------
!
      subroutine init_oa( tile          &
      ,iic_oa                           &
      ,imin, imax                 &
      ,jmin, jmax                 &
      ,kmin, kmax                 &
      ,lon_t_oa                         & !(-1:imax+2,-1:jmax+2)            
      ,lat_t_oa                         & !(-1:imax+2,-1:jmax+2)            
      ,lon_u_oa                         & !(0:imax+2,0:jmax+2)              
      ,lat_u_oa                         & !(0:imax+2,0:jmax+2)              
      ,lon_v_oa                         & !(0:imax+2,0:jmax+2)              
      ,lat_v_oa                         & !(0:imax+2,0:jmax+2)              
      ,lon_f_oa                         & !(0:imax+2,0:jmax+2)              
      ,lat_f_oa                         & !(0:imax+2,0:jmax+2)              
      ,mask_t_oa                        & !(-1:imax+2,-1:jmax+2,0:kmax+1)
      ,mask_f_oa                        & !(0:imax+1,0:jmax+1,0:kmax+1)
      ,mask_u_oa                        & !(0:imax+1,0:jmax+1,0:kmax+1)
      ,mask_v_oa                        & !(0:imax+1,0:jmax+1,0:kmax+1)
      ,h_w_oa                           & !(0:imax+1,0:jmax+1)              
      ,h_u_oa                           & !(0:imax+1,0:jmax+1)              
      ,h_v_oa                           & !(0:imax+1,0:jmax+1)              
      ,h_f_oa                           & !(0:imax+1,0:jmax+1)              
      ,rhp_t                            & 
      ,depth_t                          & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
      ) 

      use module_grd_oa, only  : associate_grid_oa_ptr, nullify_grid_oa_ptr 
      use module_tile_oa, only : allocate_begvs_oa, allocate_begvs3d_oa
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock
!      use scalars

      implicit none

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      !> Grid index range for the analysis of fields passed in argument to the OA module 
      !! Grid index are Roms-Croco tile(-thread) indices.
      !! If no tile(-thread) grid index are MPI subdomain index
      !! If no tile(-thread) and no MPI subdomain decomposition => full domain grid range
      integer, intent(in) ::       &
             imin, imax            & 
            ,jmin, jmax            &
            ,kmin, kmax

      !> tile parameter from croco
      integer, intent(in) :: tile

      !> 3-Dimensional density array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: rhp_t
                                                                          
      !> 3-Dimensional depth array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: depth_t

      ! 2-Dimensional array coordinate arrays passed from the calling code
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_t_oa  !< Longitude at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_u_oa  !< Longitude at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_v_oa  !< Longitude at v-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_f_oa  !< Longitude at f-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_t_oa  !< Latitude at t-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_u_oa  !< Latitude at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_v_oa  !< Latitude at v-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_f_oa  !< Latitude at f-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_w_oa    !< Ocean thickness at w-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_u_oa    !< Ocean thickness at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_v_oa    !< Ocean thickness at v-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_f_oa    !< Ocean thickness at f-grid point           


      !> 3-Dimensional mask arrays passed from the calling code
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_t_oa
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_u_oa
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_v_oa
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_f_oa

      !> Lower and upper index bounds of the 3-Dimensional density array (see if pst%imin,...)
      integer, dimension(3)  ::  rhp_t_lbound, rhp_t_ubound


!.....History file:
      if(if_print_node) write(io_unit,*) '...OA HISTORY file initialization'
      call history_oa(4,-1,-1,-1,-1, -1)

!.....initialisation de la structure spatiale:
!
!     BLXD croco2021 tile-thread : 
!     - s(tile)%begvs_oa   : can be alloc. bef. the 1st var_space_oa call (nzv_oa size)
!     - s(tile)%begvs3d_oa : cannot be allocated. bef. size calculated in var_space_oa

      if(if_print_node) write(io_unit,*) '...ALLOCATING allocate_begvs_oa size(nzv_oa+1)'

      call allocate_begvs_oa(imin,imax,jmin,jmax,kmin,kmax,tile) 

      if(if_print_node) write(io_unit,*) '...TEMPORARY Pointer association'

      call associate_grid_oa_ptr( tile                    &
                        ,imin, imax                       &
                        ,jmin, jmax                       &
                        ,kmin, kmax                       &
                        ,lon_t_oa                         &
                        ,lat_t_oa                         &
                        ,lon_u_oa                         &
                        ,lat_u_oa                         &
                        ,lon_v_oa                         &
                        ,lat_v_oa                         &
                        ,lon_f_oa                         &
                        ,lat_f_oa                         &
                        ,mask_t_oa                        &
                        ,mask_f_oa                        &
                        ,mask_u_oa                        &
                        ,mask_v_oa                        &
                        ,h_w_oa                           &
                        ,h_u_oa                           &
                        ,h_v_oa                           &
                        ,h_f_oa                           )

      
      if(if_print_node) write(io_unit,*) '...OA ENTERING var_space_oa .false. : calc. space dim'

      call var_space_oa(.false., tile ) 

      if(if_print_node) write(io_unit,*) '...OA ALLOCATING allocate_begvs3d_oa size(nzvs_oa+1)'

      call allocate_begvs3d_oa(tile) 

      if(if_print_node) write(io_unit,*) '...OA ENTERING var_space_oa .true. : set spatial variables'

      call var_space_oa(.true., tile )

!.....BLXD horizontal grid information no longer needed
!     nullify pointers 
       
      call nullify_grid_oa_ptr( tile ) 


!.....initialisation boite d'heisenberg:
!     attention reprendre calcul selon ondelette
      call box_oa

!     BLXD TODO CHECK if kount0, nt_max, dti needs to be passed as argt
!               this will be true if they are no longer module variables
!               => var_time_oa, var_par_oa, count_nmsimult_oa ?

!.....initialisation de la structure frequentielle:

      call var_per_oa( .false., dti )

      call allocate_part3_oa   

      call var_per_oa( .true., dti )

!.....initialisation de la structure temporelle:

      call var_time_oa(.false., kount0, nt_max, dti )

      call allocate_part4_oa

      call var_time_oa(.true., kount0, nt_max, dti )

!.....initialisation du facteur de reconstruction

      call var_rec_oa( dti )

!.....Calculate the size of the module_oa_stock array structure wf_oa

      call pre_allocate_part5_oa()

      call count_nmsimult_oa( 0                       & ! ichoix
                          ,1                          & ! ivar_m
                          ,kount0                     &
                          ,nt_max                     &
                          ,dti                        &
                          ,tile                       &
                          ,nzw_oa) 

      ! BLXD Pay attention that nzw_oa is now set to nmsimult_oa in count_nmsimult_oa
      ! nmsimult_oa and nwz_oa are both used for the same purpose 
      ! one should be now removed, 
      ! TODO namelist par. nmsimult_oa_max not really useful

      call allocate_part5_oa()

!.....History file:        
      call history_oa(6,-1,-1,-1,-1, tile)

   
!.....ecriture du fichier de structure:

      !STDALONE call struct_oa

!.....History file:
      call history_oa(5,-1,-1,-1,-1, tile, iic_oa, nt_max )

!.....test initialization if any
      if (ifl_test_oa==1) then
        call test_init_oa( imin=imin, imax=imax     &
                          ,jmin=jmin, jmax=jmax     &
                          ,kmin=kmin, kmax=kmax )
      endif

!.....initialisation eventuelle des "levels"
      if (isopycne_analysis) then
        rhp_t_lbound = (/imin,jmin,kmin/)
        rhp_t_ubound = (/imax,jmax,kmax/)
        !if ( present(rhp_t) .and. present(depth_t) ) then
            call lev_init_oa( tile, rhp_t, rhp_t_lbound, rhp_t_ubound ) 
        !else
        !    stop "3D density array must be passed as argument to the OA module (see subroutine initial_oa)" 
        !endif
      endif

      ! OLD symphonie initialisations
      !......initialisation de variables symphonie:
      !    call var_upd_oa (-1,-1,-1)
      !......initialisation de l'energie
      !     call nrj_init_oa

      call upd_init_oa( imin, imax    &
                       ,jmin, jmax    &
                       ,kmin, kmax    )

! BLXD IMPT barrier instruction : see init_online_spectral_diags_tile
! DEALLOCATES THE GRID STRUCTURE, POINTERS HAVE ALREADY BEEN NULLIFIED
! C$OMP BARRIER
!        call deallocate_tile_grid_tmp_str()

      return
      end subroutine init_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION:  initialise la structure spatiale du vecteur d'etat.
!
!> @brief The size of the horizontal and vertical domain over which analyses will be conducted 
!! is calculated, i.e., nzvs_oa and nzvs3d_oa respectively (flag_s equals true). 
!! All the variable analysis are stored in a single vector, called state vector, 
!! which vectorizes all the domain points requested for analysis. 
!! Variables enabling conversion from vector index to domain indices 
!! are initialized (flag_s equals false).
!
!> @details Spatial domain options are set in the OA namelist with SWT_D_W parameter.  
!! - if set to 1: analysis is requested over the region [lat1,lat2],[lon1,lon2]. 
!! - if set to 2: analysis is requested at point (I,J).
!! - if set to 3: analysis is requested over the ocean colum [Hmin,Hmax] in the horizontal 
!!   domain the (Lat,Lon)
!! The domain to analyse is more over defined by the namelist parameters:
!! - dx_oa(iv_g), dy_oa(iv_g) : steps in i,j indices
!! - k_oa(1,iv_g),k_oa(2,iv_g),dk_oa(iv_g) : first, last k indices and k index step
!!
!! Variables:
!! - nzvs_oa,   begvs_oa    : 2d struture of the state vector
!! - nzvs3d_oa, begvs3d_oa  : 3d structure of the state vector
!! - ij2l_oa, l2i_oa, l2j_oa : conversion l <--> (i,j)
!!
!! The section of the state vector defined by indices begvs3d_oa(ls_l), begvs3d_oa(ls_l+1)-1
!! corresponds to the points to analyse over one given (i,j) ocean column for a given variable
!! with index iv_g. For each variable, the number of (i,j) points to analyse is given by index ls_l
!! which ranges from begvs_oa(iv_g) to begvs_oa(iv_g+1)-1. 
!!
!! For each variables: 
!! - the grid point is needed to get the mask value, i.e., tgv_oa.
!! - the 2d/3d attribute is required to construct the vertical part of the state vector, i.e., tgv3d_oa.
!!
!! The kmin3d_oa parameter stores the lowest analysed k index.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - Benedicte Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Croco Tile-thread compliant version (2021)
!> @date 2021
!> @todo BLXD Modify Croco-OnlineA module interface.
!! - The CROCO-OnlineAnalysis module interface must be modified : 
!!   Interface = mix between the Stand-alone and Original OnlineA versions 
!!   => Croco arrays are passed as arguments to the OnlineA routines with reduced 
!!   dimension range. It leads to memory duplication.
!------------------------------------------------------------------------------
      subroutine var_space_oa( flag_s, tile )

      use module_grd_oa, only  : grd, pg, get_3D_grid_mask &
     &                          ,get_2D_subdomain_minmax
      use module_tile_oa, only : st, pst
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      logical, intent(in) :: flag_s
      integer, intent(in) :: tile
      

      integer :: i, j, k, grd_pt_code                                 &                          
           ,iv_g                                                      &
           ,mykmin                                                    &
           ,msk_g                                                     &
           ,ls_g

      real    ::                                                      &
            lat_g                                                     &
           ,lon_g                                                     &
           ,h_g                                                         

      real    ::                                                      &
            latmin_g                                                  &
           ,lonmin_g                                                  &
           ,latmax_g                                                  &
           ,lonmax_g

      integer :: io_nodoa, num !, ii_glob, jj_glob
      integer          :: imin, imax, jmin, jmax, kmin, kmax
      integer, pointer :: nzvs_oa => null(), nzvs3d_oa => null()

!.... Grid and Space structure : temporary pointer association
      pg  => grd(tile)
      pst => st(tile) 

!.... Initialisation of the size of the spatial structure
      pst%nzvs_oa   = 0
      pst%nzvs3d_oa = 0      

!.... Local pointers pointing toward 2D/3D the spatial structure size, resp.
!     just to handle the counts
      nzvs_oa => pst%nzvs_oa ; nzvs3d_oa => pst%nzvs3d_oa

!.... Checking grid size consistency
      if(if_print_node) write (io_unit,*) 'Dimension for grid and space structure ', tile
      if(if_print_node) write (io_unit,*) pg%imin,pst%imin
      if(if_print_node) write (io_unit,*) pg%imax,pst%imax
      if(if_print_node) write (io_unit,*) pg%jmin,pst%jmin
      if(if_print_node) write (io_unit,*) pg%jmax,pst%jmax
      if(if_print_node) write (io_unit,*) pg%kmin,pst%kmin
      if(if_print_node) write (io_unit,*) pg%kmax,pst%kmax

      if ( ( pg%imin /= pst%imin ) .or. &
           ( pg%imax /= pst%imax ) .or. &
           ( pg%jmin /= pst%jmin ) .or. &
           ( pg%jmax /= pst%jmax ) .or. &
           ( pg%kmin /= pst%kmin ) .or. &
           ( pg%kmax /= pst%kmax ) ) then
        if(if_print_node) write (io_unit,*) 'STOP : error with tile dim for grid and space structure ', tile
        stop
      else
!.... Local variable (just to shorten loop instructions)
        imin=pg%imin ; imax=pg%imax ; kmin=pg%kmin
        jmin=pg%jmin ; jmax=pg%jmax ; kmax=pg%kmax
      endif

      spa_var_loop : do iv_g = 1 , nzv_oa

#ifdef OA_TRACES
         if (flag_s) then
             io_nodoa = set_io_nodoa(iv_g,mynode,8,7) ! [7/8](odd/even-iv)000+mynode
         else
             io_nodoa = set_io_nodoa(iv_g,mynode,6,5) ! [5/6](odd/even-iv)000+mynode
         end if
#endif

         if (flag_s.eq..false.) then

            if(if_print_node) write (io_unit,*) '...var_space_oa .false. loop #iv_g  = ', iv_g

!.... Check user values
            if (k_oa(1,iv_g).eq.-99) then
               k_oa(1,iv_g) = 1 
            endif
            if (k_oa(2,iv_g).eq.-99) then
               k_oa(2,iv_g) = kmax
            endif

            if ( k_oa(1,iv_g) .gt. k_oa(2,iv_g) ) then
               if(if_print_node) write (io_unit,*) " k_oa(2,:) trop grand > kmax = N !"
               stop
            endif
            if ( k_oa(2,iv_g).gt.kmax ) then
               if(if_print_node) write (io_unit,*) " k_oa(2,:) trop grand > kmax = N !"
               stop
            endif
            if ( k_oa(1,iv_g).lt.kmin ) then
               if(if_print_node) write (io_unit,*) " k_oa(1,:) trop petit < kmin !"
               stop
            endif

!..... #BLXD Croco T-grid indices 1:N, W-grid indices from 0:N and kmin=1,kmax=N 
            if ( k_oa(1,iv_g).eq.0 .and.       &
                (tgv_oa(tv_oa(iv_g)).eq.1      &
             .or.tgv_oa(tv_oa(iv_g)).eq.2      &
             .or.tgv_oa(tv_oa(iv_g)).eq.3      &
! Humm p grid point is the 2D vorticity grid point
! Should also start to 1 and not zero ( not w-grid point )
!            .or.tgv_oa(tv_oa(iv_g)).eq.3      &
               ) ) then
                k_oa(1,iv_g)=1
                if (k_oa(2,iv_g).eq.0) k_oa(2,iv_g)=1
            endif

            if_domain_lonlat_1 : if (swt_d_oa(iv_g).eq.1.or.swt_d_oa(iv_g).eq.3) then

                ! BLD TODO ADD here missed checks dependent of imin,imax,.. ?
                !if ( if_chks ) then
                !
                !    call domain_checks()
                !end if

                if ( if_extend_dom ) then

                    ! USEFUL? TO TEST BLXD
                    grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)

                    call get_2D_subdomain_minmax( imin,imax,jmin, jmax &
                                          ,grd_pt_code, pg               &
                                          ,latmin_g, lonmin_g, latmax_g, lonmax_g )

                    call extend_domain_oa(iv_g,if_extend_dom,lonmin_g,latmin_g,lonmax_g,latmax_g)

                end if

            end if if_domain_lonlat_1

         else ! flag_s .true. 

            if(if_print_node) write (io_unit,*) '...var_space_oa .false. loop #iv_g  = ', iv_g

            ! BLXD nzvs_oa => pst%nzvs_oa which is properly initialized to zero
            pst%begvs_oa( iv_g ) = nzvs_oa + 1

         end if

!.....configuration spatiale n°1 ou 3: 

         spatial_config1_3 : if (swt_d_oa(iv_g).eq.1.or.swt_d_oa(iv_g).eq.3) then

#ifdef OA_TRACES
            if (verbose_oa>=10) then
                write (io_nodoa,*) 'VAR_SPACE_OA j=',jmin, jmax, dy_oa(iv_g),tile
                write (io_nodoa,*) 'VAR_SPACE_OA i=',imin, imax, dx_oa(iv_g),tile
            end if
#endif

            jloop : do j = jmin, jmax, dy_oa(iv_g)
                iloop : do i = imin, imax, dx_oa(iv_g)

       
                  ! BLXD TODO TEST
                  grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)
                  call get_3D_grid_mask( i, j , kmax, grd_pt_code, pg, &
                                         h_g, lat_g, lon_g, msk_g )
#ifdef OA_TRACES
            if (verbose_oa>=10) then
               write (io_nodoa,*) "all lon, lat, h, i, j ", lon_g, lat_g,h_g, i,j
            end if
#endif

#ifdef SPHERICAL
         if ( if_sphe_deg2rad ) then
! BLXD 2020 change
                lat_g=lat_g*pi_oa/180.D0
                lon_g=lon_g*pi_oa/180.D0
         endif
#endif


                  if_spatial_dom_ok_loopij : if (                     &
                       msk_g.eq.1              .and.                  &
                       lat_g.ge.lat_oa(1,iv_g) .and.                  &
                       lat_g.le.lat_oa(2,iv_g) .and.                  &
                       lon_g.ge.lon_oa(1,iv_g) .and.                  &
                       lon_g.le.lon_oa(2,iv_g) .and.                  &
                       ( (h_g.ge.h_oa(1,iv_g)  .and.                  &
                       h_g.le.h_oa(2,iv_g) )                          &
                       .or.swt_d_oa(iv_g).ne.3) )                     &
                       then

#ifdef OA_TRACES
            if ( verbose_oa>=10 ) then
               write (io_nodoa,*) "unmasked lon, lat, h, i, j ", lon_g, lat_g,h_g, i,j
            end if
#endif

!-------->structure 2d:

                     nzvs_oa          = nzvs_oa + 1
                     if (flag_s) then
                     pst%l2i_oa (nzvs_oa)  = i
                     pst%l2j_oa (nzvs_oa)  = j
                     pst%ij2l_oa(i,j,iv_g)= nzvs_oa
                     endif

!-------->structure 3d:
! #BLXD line_begin - line_end => SUBROUTINE                     
                     mykmin=0
                     ! BLXD nzvs_oa => pst%nzvs_oa which is properly initialized to zero
                     if (flag_s) pst%begvs3d_oa(nzvs_oa) = nzvs3d_oa + 1
                     if (tgv3d_oa(tv_oa(iv_g)).eq.2) then
                        nzvs3d_oa          = nzvs3d_oa + 1
                     else
                        do k = k_oa(1,iv_g),k_oa(2,iv_g),dk_oa(iv_g)
                           if (mod(tgv_oa(tv_oa(iv_g)),5).eq.1) then
                              msk_g = pg%mask_t(i,j,k) 
                           elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.2) then
                              msk_g = pg%mask_u(i,j,k) 
                           elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.3) then
                              msk_g = pg%mask_v(i,j,k) 
                           else
                              msk_g = pg%mask_f(i,j,k) 
                           endif
                           if ( msk_g.eq.1 ) then
                              if ( mykmin.eq.0.and.flag_s ) then
                                 pst%kmin3d_oa(nzvs_oa)=k            
                                 mykmin=1
                              endif
                              nzvs3d_oa          = nzvs3d_oa + 1
                           endif
                        enddo
                     endif
! #BLXD line_end                  
                  endif if_spatial_dom_ok_loopij

               enddo iloop
            enddo jloop

         endif spatial_config1_3
         
!.....configuration spatiale n°2:

         spatial_config2 : if (swt_d_oa(iv_g).eq.2) then
            i = ptij_oa(1,iv_g)
            j = ptij_oa(2,iv_g)
            if (mod(tgv_oa(tv_oa(iv_g)),5).eq.1) then
               msk_g = pg%mask_t(i,j,kmax) 
            elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.2) then
               msk_g = pg%mask_u(i,j,kmax) 
            elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.3) then
               msk_g = pg%mask_v(i,j,kmax) 
            else
               msk_g = pg%mask_f(i,j,kmax) 
            endif

            if_spatial_dom_ok_pointij : if ( msk_g.eq.1 ) then
               nzvs_oa          = nzvs_oa + 1

               if (flag_s) then
                pst%l2i_oa (nzvs_oa)  = ptij_oa(1,iv_g)
                pst%l2j_oa (nzvs_oa)  = ptij_oa(2,iv_g)
                pst%ij2l_oa(i,j,iv_g)= nzvs_oa
               endif

!-------->structure 3d:
! #BLXD line_begin 2 TIMES THESE LINES => subroutine
               mykmin=0
               if (flag_s) pst%begvs3d_oa(nzvs_oa) = nzvs3d_oa + 1
               if (tgv3d_oa(tv_oa(iv_g)).eq.2) then
                  nzvs3d_oa          = nzvs3d_oa + 1
               else
                  do k = k_oa(1,iv_g),k_oa(2,iv_g),dk_oa(iv_g)
                     if (mod(tgv_oa(tv_oa(iv_g)),5).eq.1) then
                        msk_g = pg%mask_t(i,j,k) 
                     elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.2) then
                        msk_g = pg%mask_u(i,j,k) 
                     elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.3) then
                        msk_g = pg%mask_v(i,j,k) 
                     else
                        msk_g = pg%mask_f(i,j,k) 
                     endif
                     if ( msk_g.eq.1 ) then
                        if ( mykmin.eq.0.and.flag_s ) then
                           pst%kmin3d_oa(nzvs_oa)=k            
                           mykmin=1
                        endif
                        nzvs3d_oa          = nzvs3d_oa + 1
                     endif
                  enddo
               endif
! #BLXD line_end

            endif if_spatial_dom_ok_pointij

         endif spatial_config2

      enddo spa_var_loop

!.....last point:
      if (flag_s) then

       pst%begvs3d_oa( nzvs_oa + 1 ) = nzvs3d_oa + 1
       pst%begvs_oa  ( nzv_oa  + 1 ) = nzvs_oa + 1
    
      endif

!.....Local pointers nullified
      nzvs_oa   => null()
      nzvs3d_oa => null()

!.....Specific module_tile_oa and module_grd_oa pointer nullified
      pg  => null()
      pst => null() 

      return
      end subroutine var_space_oa

      subroutine extend_domain_oa( iv_g, if_extend_dom,  lonmin_g,latmin_g,lonmax_g,latmax_g )

      use module_oa_space, only : lon_oa, lat_oa
! BLXD USEFUL ? TODO REMOVE
      implicit none

      integer, intent(in) :: iv_g

      real, intent(inout) :: lonmin_g,latmin_g,lonmax_g,latmax_g

      integer :: io_nodoa

      logical, intent(in) :: if_extend_dom

#ifdef SPHERICAL
        if ( if_sphe_deg2rad ) then

        latmin_g=latmin_g*pi_oa/180. !D0
        lonmin_g=lonmin_g*pi_oa/180. !D0
        latmax_g=latmax_g*pi_oa/180. !D0
        lonmax_g=lonmax_g*pi_oa/180. !D0

        end if
#endif

      end subroutine extend_domain_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Sets the frequency parameters as requested in namelist 
!! which controls the "state vector" frequency structure.
!
!> @details Variable :
!! - nzvp_oa the total number of periods to analyse in the simulation (all variables included). 
!! - begvp_oa stores the number of periods per variables. It is refered as to the "frequency" structure
!!   of the state vector.
!!
!! Namelist parameters are:
!! PER_W(1), per_oa(1) : minimum period
!! PER_W(2), per_oa(2) : maximum period
!! PER_W(3), per_oa(3) : step determining the number of periods to be examined, dp in sec. (-99=optimum).
!! NZPT_PER_W, nzpt_per_oa : number of points in the convolution window.
!! 
!! Outputs:
!! - resv_oa : integer setting the temporal resolution (i.e., number of model time step) at which
!!   the convolution product will be calculated (e.g., one point every 4 model time steps)
!!   TODO check this definition in the case when dori_oa equals 2.
!! - nzvp_oa and begvp_oa vector : parameters controling the frequential structure of the "state vector".
!! - perv_oa(1,:) : series of atom periods which will be examined for each cfg/var ( optimal steps if PER_W(3) = -99 ) 
!! - perv_oa(2,:) : reconstruction factor ( replaces inverse transformation ).
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - Benedicte Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - Croco-OnlineA module interface (2020)
!> @date 2021
!> @todo consistency with var_time_oa change
!! - remove goto syntax ? intent attribute ?
!------------------------------------------------------------------------------
      subroutine var_per_oa( flag_p, dti )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      double precision, intent(in) :: dti                             !< Integration time step

      integer                                                         &
         iv_p                                                         &  !< Variable index
        ,iv1_p

      real                                                            &
         ip_p

      real                                                            &   
         p_p                                                          &   !< Time scale
        ,t0c,dtc,w0c,dwc                                                  !< Wavelet resolution : Heisenberg parameters

      logical :: flag_p

      t0c = 0.
      w0c = 6.0
      dtc = 0.707106781
      dwc = dtc
      
      nzvp_oa = 0

!.....precalcul des periodes:
      do iv_p = 1, nzv_oa

! cas resolution reelle
        if (per_oa(3,iv_p).eq.-99.*unite_oa) then
         
         if (flag_p) begvp_oa( iv_p ) = nzvp_oa + 1
         if (per_oa(1,iv_p).eq.0) per_oa(1,iv_p)=1. ! on commence pas a 0      
         p_p=per_oa(1,iv_p) 


         do while (p_p.le.per_oa(2,iv_p)) 
          nzvp_oa = nzvp_oa + 1

          if (flag_p) perv_oa (1,nzvp_oa) = p_p
          p_p=p_p*(1.+dwc/w0c)/(1.-dwc/w0c)
         enddo
         nzvp_oa = nzvp_oa + 1
         if (flag_p) perv_oa (1,nzvp_oa) =p_p

! cas resolution reguliere 
       else
         if (flag_p) begvp_oa( iv_p ) = nzvp_oa + 1
!---------------------------------------------------------
!.......periode entiere en secondes
!       #BLXD the following loop is removed from the original sources
!       => ONLY one period is processed when per_oa(3,iv_p) = PER_W(3) /= -99 
!       do ip_p = int(per_oa(1,iv_p)),int(per_oa(2,iv_p)),int(per_oa(3,iv_p))
!       changer aussi la declaration de ip_p
!---------------------------------------------------------
        ip_p = per_oa(1,iv_p)

 100    continue
         nzvp_oa          = nzvp_oa + 1
         if (flag_p) perv_oa(1,nzvp_oa) = ip_p
         ip_p = ip_p + per_oa(3,iv_p)
        if (ip_p.le.per_oa(2,iv_p)) goto 100 

!---------------------------------------------------------
!       enddo
!---------------------------------------------------------
       endif

      enddo
      if (flag_p) begvp_oa(nzv_oa+1)     = nzvp_oa + 1

      if (.not.(flag_p)) return
      
!.....precalcul des resolutions:

      do iv_p = 1 , nzv_oa
        do iv1_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1

!.........atome= dirac:
          if (tpsi_oa(iv_p).eq.0.and.flag_p) then
            resv_oa (iv1_p)    = 1 
          endif

!.........atome= ondelette de morlet ou windowed fourier ou fourier:
          if (tpsi_oa(iv_p).ne.0.and.flag_p) then
           if (dori_oa(iv_p).eq.1) then ! Discrete
            !#BLXD 2020 change TODO to be consistent with changes in var_time_oa
            !if ( MOD( perv_oa (1,iv1_p), dti) /= 0 ) then

            resv_oa (iv1_p)   = max ( (  int(                         &
                             perv_oa (1,iv1_p)  /   dti               &
                              ) + 1 )                                 &
                           / nzpt_per_oa (iv_p)                       &
                                     , 1 )      
            !else
            !resv_oa (iv1_p)   = max ( (  int(                         &
            !                 perv_oa (1,iv1_p)  /   dti               &
            !                  ) )                                     &
            !               / nzpt_per_oa (iv_p)                       &
            !                         , 1 )      
            !endif 

           else
            resv_oa (iv1_p)   = max ( (  int(                         &
                             perv_oa (1,begvp_oa(iv_p))               & 
                           /   dti                                    &
                              ) + 1 )                                 &
                                      , 1)                            &
                           / nzpt_per_oa (iv_p)
            ! #BLXD : The above expression (from original code) is suspicious BUG ?
            ! Why the MAX function should be here taken bef. division by nzpt_per_oa

            resv_oa (iv1_p)   = max ( (  int(                         &
                             perv_oa (1,begvp_oa(iv_p))               & 
                           /   dti                                    &
                              ) + 1 )                                 &
                                      , 1)                            &
                           / nzpt_per_oa (iv_p)

           endif
          endif

        enddo
       enddo
      return

      end subroutine var_per_oa

      subroutine allocate_part3_oa

      ! BLXD change name allocate_nzvp_oa ?
      !use module_oa_variables
      use module_oa_time
      use module_oa_periode

      implicit none

      allocate (                                                      &   
       resv_oa (nzvp_oa)                                              &  ! resolution temporelle pour le calcul de la convolution
       )
      allocate (                                                      &  ! periodes associees a la structure vectorielle du vecteur d etat, 
       perv_oa (2,nzvp_oa)                                            &  ! facteurs de reconstruction associes a l'ondelette
       )

      return
      end subroutine allocate_part3_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Records the temporal parameters for all the requested analysis.
!
!> @details Variables :
!! - nzvt_oa is the total number of time windows requested for analysis (all variables included). 
!! - begvt_oa stores the number of time windows per variables. It is refered to as the "time" structure
!!   of the state vector.
!!
!! Corresponding namelist parameters :
!! - DELTA_T_W, delta_t_oa : set to 1 for Dirac analysis, set to 2 or more for wavelets or
!!   windowed Fourier. In the latter case, the convolution window is set 2 times the period 
!!   of interest (c.f., perv_oa(1,ip_t) ).
!! - SWT_T_W, swt_t_oa : defines the time domain where analysis are requested
!!   * 1 : from T1 to T2 #BLXD to test with tool_datetosec TODO 
!!   * 2 : single date T1 at the end, NOT AVAILABLE YET TO TEST 
!!   * 3 : single date T1, NOT AVAILABLE YET TO TEST
!!   * 4 : the entire simulation.
!!   Check : the only namelist parameters are DATE_DEB, DATE_END (pattern is 2011, 01, 1, 00, 01, 18).
!!   kount_user_oa(1:2,nzv_oa) is calculated from DATE_DEB, DATE_END in terms of model time steps.
!! - DT_W : time period at which the analysis is repeated (s). Enables to deduce kount_user_oa(3,nzv_oa) 
!!   in terms of model time steps. 
!!   
!! - lt_p = begvt_oa(iv), begvt_oa(iv+1)-1
!!   is a "pointer" to all the requested analysis for a given variable-configuration
!!   including the requested number of simulation steps where to perform the analysis and all the requested period of analysis
!!   For a given lt_p, one can retreive the corresponding analysis time period per_t2p_oa(lt_p)
!!   To each lt_p corresponds :
!!   1) a convolution window [kountv_oa(1,lt_p), kountv_oa(2,lt_p)] centered at the simulation time step 
!!      [kountv_oa(2) + kountv_oa(1)]/2 and set according to :
!!   - the type of atom tpsi_oa(iv_p)
!!   - the analysis time period (delta_t_w) 
!!   2) a period of analysis recovered through ip_t = per_t2p_oa(lt_p)
!!      the requested period can be :
!!      - "discrete" and for each variable-cfg retrieved with ipt = begvp_oa(iv_t), begvp_oa(ivt+1)-1
!!         if dori_oa(iv_t) is set to 1 
!!
!! Outputs: 
!! - kountv_oa(1,nzvt_oa), kountv_oa(2,nzvt_oa) are delimiting each convolution window (for each variables).
!! - per_t2p_oa (nzvt_oa) : correspondance between "temporal" and "frequency" structure ( loop on cfg/var, loop on window => begvp_oa(iv)).
!! - dkount_oa  (1:nzv_oa): duree d'echantillonnage,
!! - nztime_v_oa(1:nzv_oa): nombre reel de points d'echantillonnage par periode.    
!! - nzvt_oa, begvt_oa    : structure temporelle du vecteur d'etat.
!! - nzw_oa , begw_oa     : structure complete (spatiale et temporelle) du vecteur d'etat.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments 
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo
!------------------------------------------------------------------------------
      subroutine var_time_oa( flag_t, kount0, nt_max, dti )
      
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      double precision, intent(in) :: dti                                !< Integration time step
      integer, intent(in)          :: kount0,                         &  !< First simulation iteration
                                      nt_max                             !< Last integration iteration

      logical, intent(in) :: flag_t

      integer                                                         &
            k                                                         &
           ,dkount_tot_t                                              &
           ,iv_t                                                      &
           !,ls_t                                                      &
           !,ls1_t                                                     &
           ,ip_t                                                      &
           ,ip1_t                                                     &
           ,ip2_t                                                     !&
           !,ls_s,ls1_s                                                &
           !,i_s,j_s,k_s

      integer, dimension(1:nzv_oa) :: comptt
      !character*40                                                    &
      !     nameid

      integer :: nzw_oa_tmp

      nzvt_oa  = 0
      nzw_oa_tmp   = 0


!---- >boucle sur toutes les variables:
      loop_on_variables : do iv_t = 1, nzv_oa

          if_first_rec_oa(iv_t)=0

!*******************************************************************
!     preparation de la structure temporelle
!*******************************************************************
         
         if (flag_t) begvt_oa(iv_t)   = nzvt_oa + 1

         if (dori_oa(iv_t).eq.1) then 
            ip1_t = begvp_oa( iv_t     ) 
            ip2_t = begvp_oa( iv_t + 1 ) - 1
         else
            ip1_t = begvp_oa( iv_t + 1 ) - 1 
            ip2_t = begvp_oa( iv_t + 1 ) - 1
         endif

!----->boucle sur toutes les periodes demandees pour cette variable:
         loop_on_requested_period : do ip_t = ip1_t , ip2_t

!----------------------------------------------------------------
!     precalcul de la largeur en kount de l'atome (from time periods perv_oa):
!     dkount_tot_t for 0: Dirac, 1: Wavelette, 2: Windowed Fourier
!----------------------------------------------------------------

!...........atome = dirac:
            if (tpsi_oa(iv_t).eq.0) then
               if (perv_oa(1,ip_t).eq.0.) then
                  dkount_tot_t    = 0
               else
                  dkount_tot_t    = int(perv_oa(1,ip_t)/2./dti-0.5)+1
               endif
            endif   

!...........atome = ondelette de morlet ou windowed fourier:     
            if (tpsi_oa(iv_t).eq.1.or.tpsi_oa(iv_t).eq.2) then
               ! #BLXD 2020 WARNING : introduction of a special case when perv is divible by dti 
               ! The integer parameter dkount_tot_t defines the width of the convolution window 
               ! The convolution window centered at the k-time index is latter defined as
               ! the the time index intervals [k-dkount_tot_t,k+dkount_tot_t].
               ! Objective : setting dkount_tot_t in order to define a convolution window 
               !             at least as large as 2 times the time period requested by 
               !             the user ( ie, perv_oa )
               ! To meet this objective :
               ! dkount_tot_t is not taken as the integer part of the division of perv_oa by dti 
               ! (ie, p) but set to p+1 (ie, ceiling operation since perv_oa is positive):
               !       p * dti < perv_oa <= (p + 1) * dti
               ! However :
               ! Being given the user defined multiplicative factor delta_t_oa >= 1, 
               ! this implies "unexpected total number of wavelet analysis" 
               ! when the simulation lasts exactly p times the time period perv_oa
               ! Proposed modification : introduction of a special case when perv is divible by dti 
               ! TODO extend this tho var_per_oa (see resv_oa)
               if (dori_oa(iv_t).eq.1) then
                  if ( MOD( perv_oa (1,ip_t), dti) /= 0 ) then
                  dkount_tot_t    = (int(                             &
                       perv_oa (1,ip_t)  /   dti                     &
                       ) + 1)                                         &
                       * delta_t_oa(iv_t)
                  else 
                    dkount_tot_t    = int(                            &
                        perv_oa (1,ip_t) * delta_t_oa(iv_t) /   dti )                             
                  endif 
               else
                  dkount_tot_t    = (int(                             &
                       perv_oa (1,ip2_t)  /   dti                  &
                       ) + 1)                                         &
                       * delta_t_oa(iv_t)
               endif
            endif

!...........atome = fourier :
            if (tpsi_oa(iv_t).eq.3) then
               dkount_tot_t = 99999.
            endif   

!...........initialisation du compteur de la variable:
            comptt(iv_t)=0

!----------------------------------------------------------------
!     configuration temporelle numero 1: periode
!     configuration temporelle numero 4: simu
!     
!     atomes differents de fourier => Dirac, Morlet or windowed Fourier
!----------------------------------------------------------------
            not_fourier : if ((swt_t_oa(iv_t).eq.1.or.swt_t_oa(iv_t).eq.4).and.tpsi_oa(iv_t).ne.3 ) then

               simulation_window : do k = kount_user_oa(1,iv_t),kount_user_oa(2,iv_t),kount_user_oa(3,iv_t)
                  nzvt_oa             = nzvt_oa + 1

                  if (flag_t) then ! second pass
                   kountv_oa(1,nzvt_oa) = k - dkount_tot_t
                   kountv_oa(2,nzvt_oa) = k + dkount_tot_t
                   per_t2p_oa (nzvt_oa)  = ip_t

!----------------------------------------------------------------
                   ! #BLXD 2020 change 
                   ! The nzvt_oa convolution windows are defined over the model time index 
                   ! intervals [ kount_oa(1,:), kount_oa(2,:) ]
                   ! Being given the width of the convolution window 2*dkount_tot_t
                   ! and its time index center k and being given the simulation time window :
                   ! - the 1st possible convolution window can start at the 1st model time index
                   ! which is the 'now' time index of the restart just before the 1st 
                   ! model time step from now to after, ie kountv_oa(1, izvt_oa=1 ) = kount0
                   ! (see call init_oa followed by main_oa at k=kount0)
                   ! - the last possible convolution window can have its last window time index
                   ! at the last model time index, which is to say at 'after' state of the last
                   ! - the last time index of the last possible convolution window is allowed to 
                   !   match the 'after' index of the last time step of the model, i.e.,
                   !   kountv_oa( 2, izvt_oa=nzvt_oa ) = nt_max
                   ! => Excluded model time steps (flag -9999) should be 
                   !    therefore defined as follows (see var_rec_oa) 
                   ! ORIGINAL TESTS
                   !if (    kountv_oa(1,nzvt_oa).lt.kount0             &
                   !    .or.kountv_oa(1,nzvt_oa).gt.nt_max-1               &
                   !    .or.kountv_oa(2,nzvt_oa).lt.kount0             &
                   !    .or.kountv_oa(2,nzvt_oa).gt.nt_max-1 )             &
                   !    then
                   if (    kountv_oa(1,nzvt_oa).lt.kount0             &
                       .or.kountv_oa(1,nzvt_oa).gt.nt_max               &
                       .or.kountv_oa(2,nzvt_oa).lt.kount0             &
                       .or.kountv_oa(2,nzvt_oa).gt.nt_max )             &
                       then


                     kountv_oa(1,nzvt_oa) = -9999
                     kountv_oa(2,nzvt_oa) = -9999
                   else
!----------------------------------------------------------------
                     ! #BLXD What is the purpose of updv_oa == 1 ?
                     ! if the above test is false non of the ORIGINAL tests
                     ! after updv_oa can be found true
                     if ( if_first_rec_oa(iv_t)==0 ) then
                      ltrec_oa(iv_t) = nzvt_oa
                      if_first_rec_oa(iv_t) = 1
                     end if 
                     if (updv_oa(iv_t).eq.1) then
                      if (kountv_oa(1,nzvt_oa).lt.kount0) then
                       kountv_oa(1,nzvt_oa) = kount0
                       kountv_oa(2,nzvt_oa) = kount0 + 2*dkount_tot_t
                      endif
                      if (kountv_oa(2,nzvt_oa).gt.nt_max-1) then
                       kountv_oa(1,nzvt_oa) = nt_max     - 2*dkount_tot_t
                       kountv_oa(2,nzvt_oa) = nt_max
                      endif
                     endif

                     comptt(iv_t)=comptt(iv_t)+1
                   endif
                  endif ! second pass
               enddo simulation_window
            endif not_fourier


!----------------------------------------------------------------
!     configuration temporelle numero 4   : simu
!     
!     atome: fourier
!----------------------------------------------------------------
            ! BLXD I would prefer
            ! if ((swt_t_oa(iv_t).eq.1.or.swt_t_oa(iv_t).eq.4).and.tpsi_oa(iv_t).eq.3 ) then
            ! But ok with init_parameter_oa : tpsi_oa == 3 => swt_t_oa(iv_t).eq.4)
            fourier : if ( swt_t_oa(iv_t).eq.4 .and. tpsi_oa(iv_t).eq.3 ) then
               nzvt_oa             = nzvt_oa + 1

               if (flag_t) then
                kountv_oa(1,nzvt_oa) = nt_max-1                       &
                    - int ( int( real(nt_max-1-kount0+1) * dti        &
                           / perv_oa(1,ip_t) )                        &
                         * perv_oa(1,ip_t)/dti ) + 1             
                kountv_oa(2,nzvt_oa) = nt_max-1  
                per_t2p_oa (nzvt_oa) = ip_t
                ! One single Fourier analysis using all the simulation window for integration  
                comptt(iv_t)       = 1
               endif

            endif fourier


!----------------------------------------------------------------
!     configuration temporelle numero 2: une date a la fin de la simulation
!     0 : dirac, 1 : ondelette, 2 : windowed Fourier
!----------------------------------------------------------------
            if (swt_t_oa(iv_t).eq.2) then
               nzvt_oa             = nzvt_oa + 1
               if (flag_t) then
                kountv_oa(1,nzvt_oa) = nt_max - 1 - 2*dkount_tot_t
                kountv_oa(2,nzvt_oa) = nt_max - 1      
! BLXD TODO change nt_max - 1  to nt_max
                !kountv_oa(1,nzvt_oa) = nt_max - 2*dkount_tot_t
                !kountv_oa(2,nzvt_oa) = nt_max      
                per_t2p_oa (nzvt_oa)  = ip_t
! BLXD TODO change nt_max - 1  to nt_max in the if instruction as well
                if (    kountv_oa(1,nzvt_oa).lt.kount0                 &
                    .or.kountv_oa(1,nzvt_oa).gt.nt_max-1                   &
                    .or.kountv_oa(2,nzvt_oa).lt.kount0                 &
                    .or.kountv_oa(2,nzvt_oa).gt.nt_max-1 )                 &
                    then
                  kountv_oa(1,nzvt_oa) = -9999
                  kountv_oa(2,nzvt_oa) = -9999
                else
! BLXD one vriable a single date of conv. window at a single date
                   comptt(iv_t)       = 1
                endif
               endif
            endif

!----------------------------------------------------------------
!     configuration temporelle numero 3: date unique
!----------------------------------------------------------------
            if (swt_t_oa(iv_t).eq.3) then
               nzvt_oa             = nzvt_oa + 1

               if (flag_t) then
! BLXD TODO change nt_max - 1  to nt_max in the if instruction as well
               kountv_oa(1,nzvt_oa) = kount_user_oa(1,iv_t) - dkount_tot_t
               kountv_oa(2,nzvt_oa) = kount_user_oa(1,iv_t) + dkount_tot_t
               per_t2p_oa (nzvt_oa) = ip_t
               if (    kountv_oa(1,nzvt_oa).lt.kount0                  &
                    .or.kountv_oa(1,nzvt_oa).gt.nt_max - 1                 &
                    .or.kountv_oa(2,nzvt_oa).lt.kount0                 &
                    .or.kountv_oa(1,nzvt_oa).gt.nt_max - 1 )               &
                    then
                  kountv_oa(1,nzvt_oa) = -9999
                  kountv_oa(2,nzvt_oa) = -9999
               else
! BLXD one variable a single date of conv. window at a single date
                   comptt(iv_t)       = 1
               endif
               endif
            endif

!.....enddo ip_t

            ! #BLD TODO CHECK counting effective analysis ( <= # of requested OA nzv_oa)
            if (flag_t) nzw_oa_tmp   = nzw_oa_tmp + comptt(iv_t)

         enddo loop_on_requested_period

         
!---- >dernier point:
         if (flag_t) begvt_oa( nzv_oa + 1 )   = nzvt_oa + 1

!.....enddo: boucle sur les variables...

      enddo loop_on_variables

#ifdef OA_TRACES
         if (flag_t) then
             if(if_print_node) write (io_unit,*) '#REQUESTED OA nzvt_oa     = ', nzvt_oa-1
             if(if_print_node) write (io_unit,*) '#EFFECTIVE OA nzw_oa_tmp  = ', nzw_oa_tmp
             nzw_oa = nzw_oa_tmp
         end if
#endif

      return
      end subroutine var_time_oa

      subroutine allocate_part4_oa

      ! BLXD change name allocate_nzvt_oa ?
      !use module_oa_variables
      use module_oa_time
      !use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      ! calcul des kounts de debut et de fin pour chaque variable (ou configuration)
      allocate (                                                      &
       kountv_oa (2,nzvt_oa)                                          &
       ) 
    
      allocate (                                                      & 
       tallocated_oa(nzvt_oa)                                         &
       )

      allocate (                                                      &
       per_t2p_oa (nzvt_oa)                                           &  ! transformation structure temporelle --> structure frequentielle
       )

      return
      end subroutine allocate_part4_oa

      subroutine pre_allocate_part5_oa

      use module_oa_stock

      implicit none

      tallocated_oa(:)   = -1
      return
      end subroutine pre_allocate_part5_oa

      subroutine allocate_part5_oa

      !use module_oa_variables
      use module_oa_time
      !use module_oa_space
      !use module_oa_periode
      use module_oa_stock

      implicit none

      integer                                                         &
       l_a


      allocate (                                                      &
       wf_oa (nmsimult_oa)                                            &           ! vecteur d etat resultat (w) contenant les analyses
       )

      do l_a=1,nmsimult_oa
         nullify(wf_oa(l_a)%coef)
      enddo

      tallocated_oa(:)   = -1
      wf_oa(:)%t_indice  = -1
      wf_oa(:)%config    = -1
      wf_oa(:)%variable  = -1


      return
      end subroutine allocate_part5_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
! Optimization de la memoire pour allouer la structure temporelle
! du vecteur d'etat a sa juste taille qui correspond au nombre
! maximum de fenetres de convolution ouvertes a un meme pas de temps 
! de simulation.
! 
!> @brief memory optimization pre-calculating the size of 
!   temporal structure of the state vector (wf_oa dimension) 
!   using the exact same algo applied in main_oa.
!
!> @details nmsimult_oa based on the number of requested 
!!  config, variable, and frequency of analysis during the simulation
!   vectorize the analysis.
!!   The state vector is an derived data type array whose dimension
!!   is the # of convolution windows (time structure) and whose field
!!   is the # of grid point to analyse (spatial structure). 
!!
! REVISION HISTORY:
!
!> @authors
!! - B. Lemieux-Dudon
!!  - Croco-OA interface 1st version (2020)
!!  - Croco Tile-thread compliant version (2021)
!!  - Memory optimization
!> @date 2021
!
!> @todo nmsimult_oa_max in namelist can be kept to stop if the
!!  number of analysis ovrepasses if the number of analyses 
!!  exceeds the memory capacity
!> @todo declare psi as real(KIND=wp) with specific module handl. prec. 
!! @todo revise if convolution windows are simulatenously defined over
!!   the slow and fast modes.
!! @todo revise as the algo for main_oa subroutine evolves. 
!----------------------------------------------------------------------

      subroutine count_nmsimult_oa( ichoix                 &
                          ,ivar_m                          &
                          ,kount0                          &
                          ,nt_max                          &
                          ,dti                             &
                          ,tile                            &
                          ,nzw_oa_tmp) 


      use module_tile_oa, only : st, pst, ntiles
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock
!      use module_oa_upd
!      use scalars

      implicit none

      !> Tile dependent calculation, ichoix set to zero in "regular" cases
      integer, intent(in) :: tile, ichoix

      !> To handle specific calls to main_oa with variables updated at different position in the calling code
      integer, intent(in), optional :: ivar_m

      !> Time integration step
      double precision, intent(in) :: dti

      !> First and last simulation iteration index
      integer, intent(in) :: nt_max, kount0                                        

      !> To set the mx number of conv. windows simultaneously under calculatio,
      integer, intent(inout) :: nzw_oa_tmp

      !TODO declare flag here if only needed here
      logical :: ifl_test_composite

      integer                                                         &
         iv_m                                                         &
        ,ic_m                                                         &
        ,lt_m                                                         &
        ,ls_m                                                         &
        ,ls1_m

      integer                                                         &
         i_m                                                          &
        ,j_m                                                          &
        ,k_m                                                          &
        ,kpt_m                                                        &
        ,lp_m                                                          

      !> For the model integration iteration
      integer :: iic_oa 

      !> Count of total size of the stock structure / process IO
      integer :: imsimult_oa, io_nodoa

      !> For checking the count of spatial structure size coef
      !! might set the logical to false if to mudch time consuling at init time
      logical, parameter :: if_count_space=.true. 
      integer, dimension(1:ntiles,1:nzc_oa,1:nzv_oa,1:nzvt_oa) ::     &
         intper_count_space                                           &
        ,sglper_count_space
   
! BLXD : double precision ? 
      complex                                                         &
         psi_m

!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

!.... The size of the module_oa_stock structure wf_oa is only dependent of
!     conf ic_m, var iv_m , and specific conv. window lt_m
!     It is independent of the tile which is related to the 2D horizontal domain
!     "tile-decomposition"
      imsimult_oa = 0 
      nmsimult_oa = 0 

!.... For each conf ic_m, var iv_m , and specific conv. window lt_m
!     One can check the # of analysed spatial points i,j,k 
      sglper_count_space(tile, :,:,:) = 0
      intper_count_space(tile, :,:,:) = 0


      oa_analysis_requested : if (nzv_oa.ne.0) then

!.... Pointer to tile structure
      pst => st(tile)

!.....flag variable composite 
      ifl_test_composite = .true.      

      loop_timestepping : do iic_oa = kount0, nt_max

      loop_ic : do ic_m = 1 , nzc_oa
      loop_iv : do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1
         !STDALONE to call main_oa at specific location for composite varaible confguration
         if (present(ivar_m)) then
            ifl_test_composite = (tv_oa(iv_m).eq.ivar_m)
         endif

!........ce test est utilise pour les appels cibles lors des calculs energetiques:
!         if ( (ichoix.ge.0.or.ifl_test_composite) .and. ( tv_oa(iv_m).ne.20) ) then
!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

     if_discrete_per : if (dori_oa(iv_m).eq.1) then

      time_of_analysis_loop : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1

!.....test pour le calcul:
       tstp_match_convol_window : if ( iic_oa.ge.kountv_oa(1,lt_m) .and. iic_oa.le.kountv_oa(2,lt_m)                  &
            .and.               mod( int( (iic_oa-kountv_oa(1,lt_m)), kind=8 ), resv_oa(per_t2p_oa(lt_m)) ).eq.0      &
            .and.               ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 )                                           &
                             ) then

! #BLXD ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 ) = calling main_oa for particular variable

!........allocation if necessary (if variable isn't allocated yet, it is done.., at user-specified kount):


!........precalcul de l'indice temporel prenant en compte la translation (u):

         kpt_m =  iic_oa - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )

          if (tallocated_oa(lt_m).eq.-1) then
            imsimult_oa = imsimult_oa + 1
            nmsimult_oa = max(nmsimult_oa,imsimult_oa)
            !tallocated_oa(lt_m) = imsimult_oa
            tallocated_oa(lt_m) = lt_m 
          endif

!........precalcul de l'indice temporel prenant en compte la translation (u):
!
         kpt_m =  iic_oa - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )

!........precalcul de psi:

          psi_m = psi_oa(                                                       &
                    tpsi_oa(iv_m)                                               &
                   ,psi_p2s_oa(                                                 &
                                     tpsi_oa(iv_m)                              &         
                                    ,perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa )  &                                     
                   ,real(kpt_m)* dti                                            &
                   ,dti*resv_oa(per_t2p_oa(lt_m))                               &           
                   ,fb_oa,fc_oa  ) 
#ifdef OA_TRACES
         if ( ( lt_m==ltrec_oa(iv_m) ) ) then
         if (iv_m==1) then
         io_nodoa = 10000+mynode
         write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)),  perv_oa(2,per_t2p_oa(lt_m))
         else if (iv_m==2) then
         io_nodoa = 20000+mynode
         write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)),  perv_oa(2,per_t2p_oa(lt_m))
         end if
         end if
#endif
!........boucle spatiale:

           !horizontal_space_loop : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
           ! i_m = pst%l2i_oa(ls_m)
           ! j_m = pst%l2j_oa(ls_m)
           !
           ! vertical_space_loop : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
           ! k_m = pst%kmin3d_oa(ls_m) + (ls1_m - pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
           ! sglper_count_space(tile, ic_m,iv_m,lt_m) = sglper_count_space(tile, ic_m,iv_m,lt_m) + 1
           !
           !
           ! enddo vertical_space_loop
           !enddo  horizontal_space_loop

       else if ( (kountv_oa(2,lt_m).le.iic_oa.or.iic_oa.eq.nt_max)      &
     &        .and.(ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)               &
     &        .and.(ichoix.eq.des_oa(iv_m)                               &
     &        .or.des_oa(iv_m).eq.0)                                     &
     &                        ) then
            if ( tallocated_oa(lt_m)==lt_m) then
               if ( imsimult_oa .lt. 1 ) then
                write(io_unit,*)'ERROR : isimult_oa cannot be samller than 1 ',mynode
               endif
               imsimult_oa = imsimult_oa - 1
               tallocated_oa(lt_m)=-1
            endif
       endif tstp_match_convol_window

      enddo time_of_analysis_loop
 
!.....periodes integrees:
!---------------------------------------------------------------------------
      elseif (dori_oa(iv_m).eq.2) then if_discrete_per
!---------------------------------------------------------------------------

      time_of_analysis_loop2 : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1

!.....test pour le calcul:

       tstp_match_convol_window2 : if ( iic_oa.ge.kountv_oa(1,lt_m).and.                                     &
                               iic_oa.le.kountv_oa(2,lt_m).and.                                      &
                               mod( int((iic_oa-kountv_oa(1,lt_m)),kind=8),                          &
                                    resv_oa(per_t2p_oa(lt_m)) ).eq.0                                 &
                                     .and.                                                           &
                                     (ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)                         &
                             ) then

!........allocation si necessaire:
          if (tallocated_oa(lt_m).eq.-1) then
            imsimult_oa = imsimult_oa + 1
            nmsimult_oa = max(nmsimult_oa,imsimult_oa)
            tallocated_oa(lt_m) = lt_m
          endif


!........boucle spatiale:
           !horizontal_space_loop2 : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
           ! i_m = pst%l2i_oa(ls_m)
           ! j_m = pst%l2j_oa(ls_m)
           !
           ! vertical_space_loop2 : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
           ! k_m = pst%kmin3d_oa(ls_m) + (ls1_m-pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
           !
           ! intper_count_space(tile, ic_m,iv_m,lt_m) = intper_count_space(tile, ic_m,iv_m,lt_m) + 1
           !
           !
           ! enddo vertical_space_loop2
           !enddo horizontal_space_loop2
         !enddo period_of_analysis_loop2

       else if ( (kountv_oa(2,lt_m).le.iic_oa.or.iic_oa.eq.nt_max)      &
             .and.(ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)                            &
             .and.(ichoix.eq.des_oa(iv_m)                                            &
             .or.des_oa(iv_m).eq.0)                                                  &
                             ) then
            if ( tallocated_oa(lt_m)==lt_m) then
               if ( imsimult_oa .lt. 1 ) then
                write(io_unit,*)'ERROR : isimult_oa cannot be samller than 1 ',mynode
               endif
            imsimult_oa = imsimult_oa - 1
            tallocated_oa(lt_m)=-1
            endif
       endif tstp_match_convol_window2

      enddo time_of_analysis_loop2
      endif if_discrete_per

!.....enddo associe a iv_m et ic_m:
      enddo loop_iv
      enddo loop_ic

      enddo loop_timestepping

      if (if_print_node) then
        write(io_unit,*)'... OA precalculation of wf_oa array size (per tile) = ',tile,nmsimult_oa
      endif
      nzw_oa_tmp = nmsimult_oa !BLXD TODO clean nzw_oa or nmsimult_oa

      !do ic_m = 1 , nzc_oa
      !do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1
      !  if (dori_oa(iv_m).eq.1) then
      !      do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1
      !      if (if_print_node) write(io_unit,*)'tile,ic,tv_oa,SGLPER OUT = ',tile,ic_m, tv_oa(iv_m),sglper_count_space(tile, ic_m,iv_m,lt_m)
      !      end do
      !  elseif (dori_oa(iv_m).eq.2) then
      !      do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1
      !      if (if_print_node) write(io_unit,*)'tile,ic,tv_oa,INTPER OUT = ',tile,ic_m, tv_oa(iv_m),intper_count_space(tile, ic_m,iv_m,lt_m)
      !      end do
      !  endif
      !enddo
      !enddo

!.....Pointer to tile struct. nullified

      pst => null()

      endif oa_analysis_requested

      return
      end subroutine count_nmsimult_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note 
!
!
! DESCRIPTION: 
!
!> @brief Initialise les coefficients de reconstruction pour chaque periode. 
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments 
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine var_rec_oa( dti )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      double precision, intent(in) :: dti    !< Integration time step

      integer kpt_m
      real    signal_r
 !STDALONE      real    pi_p
   
! BLXD : DBLE PRECISION ? 
      complex ::                                                      &
              signal_rec                                              &
             ,temp_r                                                  &
             ,psi_m

      integer                                                         &
            ip_p                                                      &
           ,iv_p                                                      &
           ,lt_p                                                      &
           ,kount_p

      integer :: io_nodoa
#ifdef OA_TRACES
      integer, dimension(1:2), parameter :: io_trc=(1000,2000)
      ! BLXD integer, dimension(1:nzv_oa) :: io_traces
#endif
!      pi_p = acos(-1.)

      do iv_p = 1 , nzv_oa
!---------------------------------------------------------------------------
!.....initialisation (pour tous les atomes):
!---------------------------------------------------------------------------
         do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1
            perv_oa(2,ip_p) = 1.
         enddo

!.....pour certains atomes, reconstruction possible si demandee:
         reconstruct_requested : if (fl_rec_oa(iv_p).eq.1) then

!---------------------------------------------------------------------------
!.....dirac:
!---------------------------------------------------------------------------
!TODO replace some if-endif by if-else if-endif

            type_of_oa_atom : if (tpsi_oa(iv_p).eq.0) then

               do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1
                  if (perv_oa(1,ip_p).eq.0) then
                     perv_oa(2,ip_p) = 1.
                  else
                     perv_oa(2,ip_p)=real(2*(int(perv_oa(1,ip_p)      & 
                          /2./dti-0.5) )/1)                        &
                          +3.
                  endif
               enddo
!TODO replace some if-endif by if-else if-endif
!            endif


!---------------------------------------------------------------------------
!.....fourier:
!---------------------------------------------------------------------------
!TODO replace some if-endif by if-else if-endif
!            if (tpsi_oa(iv_p).eq.2.or.tpsi_oa(iv_p).eq.3) then 

            else if (tpsi_oa(iv_p).eq.2.or.tpsi_oa(iv_p).eq.3) then type_of_oa_atom 
               lt_p  = begvt_oa(iv_p)
               if (dori_oa(iv_p).eq.2) stop 'a programmer...'
               do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1
                  do while ((per_t2p_oa(lt_p).ne.ip_p.or.kountv_oa(1,lt_p).eq.-9999.or.kountv_oa(2,lt_p).eq.-9999).and.lt_p.lt.begvt_oa(iv_p+1) )
                     lt_p=lt_p+1
                  enddo
                  perv_oa(2,ip_p)=real((kountv_oa(2,lt_p)-kountv_oa(1,lt_p))/resv_oa(ip_p)+1)
               enddo
              
!TODO replace some if-endif by if-else if-endif
!            endif


!---------------------------------------------------------------------------
!.....ondelette de morlet:
!---------------------------------------------------------------------------
!TODO replace some if-endif by if-else if-endif
!            if (tpsi_oa(iv_p).eq.1) then

             else if (tpsi_oa(iv_p).eq.1) then type_of_oa_atom 

!.....periodes discretisees:
               if (dori_oa(iv_p).eq.1) then
                  

                  do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1

                     temp_r = 0.

!     boucle temporelle

                     lt_p = begvt_oa(iv_p)

                     do while (lt_p.lt.begvt_oa(iv_p+1).and.(per_t2p_oa(lt_p).ne.ip_p.or.kountv_oa(1,lt_p).eq.-9999.or.kountv_oa(2,lt_p).eq.-9999))
                     ! #BLXD issue with array per_t2p_oa bounds overflow with intel ifort
                     !      Even if lt_p < begvt_oa is false the following conditions are tested (per_t2p_oa(lt_p)/=ip_p) etc
                     !      => run time failure 
                     ! Adding exit condition within the loop 
                     !do while (lt_p.lt.begvt_oa(iv_p+1).and.(kountv_oa(1,lt_p).eq.-9999.or.kountv_oa(2,lt_p).eq.-9999))

                        lt_p = lt_p + 1
                        if (lt_p.eq.begvt_oa(iv_p+1)) then
                            exit
                        end if
                     enddo
                    
                     ! #BLXD 1st temporal pointer (requested by the user) for which the convolution can be performed (if any)
                     if (lt_p.lt.begvt_oa(iv_p+1)) then
                        do  kount_p =  kountv_oa(1,lt_p),kountv_oa(2,lt_p),resv_oa(per_t2p_oa(lt_p))
                           kpt_m = kount_p - ( int(( kountv_oa(1,lt_p)+kountv_oa(2,lt_p))/2) )

                           
!     precalcul de psi:
                           
                           psi_m = psi_oa(                       &
                                tpsi_oa(iv_p)                         &
                                ,psi_p2s_oa(                     &
                                tpsi_oa(iv_p)                         &
                                ,perv_oa(1,ip_p),fb_oa,fc_oa )        &
                                ,real(kpt_m)* dti                  &   !- tempschoisi*dti
                                ,dti*resv_oa(per_t2p_oa(lt_p))     &
                                ,fb_oa,fc_oa  )

!     cosinus signal test
                           signal_r=cos(2*pi_oa/perv_oa(1,ip_p) *real(kpt_m)* dti)  
!                                                              *real(kount_p)* dti)  

!     reconstruction factor
                           temp_r = temp_r + conjg(psi_m) * signal_r
#ifdef OA_TRACES
                           if ( ( lt_p==ltrec_oa(iv_p) ) ) then
                             if (iv_p==1) then 
                             io_nodoa = io_trc(1)+mynode !1000+mynode
                             write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)), signal_r
                             else if (iv_p==2) then
                             io_nodoa = io_trc(2)+mynode !2000+mynode
                             write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)), signal_r
                             end if
                            end if
#endif

                        enddo

!     fin boucle temporelle

!     signal_rec=cos(2*pi_oa/perv_oa(1,ip_p)*real(kpt_m)* dti)
!    &  +(0.,1.)*sin(2*pi_oa/perv_oa(1,ip_p)*real(kpt_m)* dti) 
!     perv_oa(2,ip_p)=abs(perv_oa(2,ip_p)/signal_rec)

                        perv_oa(2,ip_p) = real(temp_r)

#ifdef OA_TRACES
                        if(if_print_node) write (io_unit,*) '       => REC FACTOR ip_p ?  ',perv_oa(1,ip_p), perv_oa(2,ip_p)
                        if ( ( lt_p==ltrec_oa(iv_p) ) ) then
                            if (iv_p==1) then
                             io_nodoa = io_trc(1)+mynode !1000+mynode
                            write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)), perv_oa(2,ip_p)
                            else if (iv_p==2) then
                             io_nodoa = io_trc(2)+mynode !2000+mynode
                            write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)), perv_oa(2,ip_p)
                            end if
                        end if
#endif
                     else
 

                     endif
                  enddo

!.....periodes integrees:
               elseif (dori_oa(iv_p).eq.2) then

                  do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1

                     temp_r = 0.

!     boucle temporelle
                     
                     lt_p = begvt_oa(iv_p)
                     do  kount_p = kountv_oa(1,lt_p) , kountv_oa(2,lt_p),resv_oa(per_t2p_oa(lt_p))
                        
                        kpt_m = kount_p - ( int((kountv_oa(1,lt_p) +kountv_oa(2,lt_p))/2) )

!     precalcul de psi:
                        
                        psi_m = psi_oa(                          &
                             tpsi_oa(iv_p)                            &
                             ,psi_p2s_oa(                        &
                             tpsi_oa(iv_p)                            &
                             ,perv_oa(1,ip_p),fb_oa,fc_oa )           &
                             ,real(kpt_m)* dti                     &  !- tempschoisi*dti
                             ,dti*resv_oa(per_t2p_oa(lt_p))        &
                             ,fb_oa,fc_oa  )

!     cosinus signal test
                        signal_r=cos(2*pi_oa/perv_oa(1,ip_p)*real(kpt_m) * dti)  
                        
!     reconstruction factor
                        temp_r = temp_r + conjg(psi_m) * signal_r

                     enddo
!     fin boucle temporelle

!     signal_rec=cos(2*pi_oa/perv_oa(1,ip_p)*real(! tempschoisi)* dti)
!     &  +(0.,1.)*sin(2*pi_oa/perv_oa(1,ip_p)*real(! tempschoisi)* dti) 
!     perv_oa(2,ip_p)=real(perv_oa(2,ip_p)/signal_rec)

                     perv_oa(2,ip_p) = real(temp_r)

                  enddo

!.....endif associe au test sur dori_oa
               endif

            endif type_of_oa_atom

!........endif associe au test: fl_rec_w=1
         endif reconstruct_requested

!......enddo associe a la variable iv_p:
      enddo
     
      return
      end subroutine var_rec_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note 
!
! DESCRIPTION: preparation des variables communes.
!
!> @brief Allocate and initializes output variables and related parameters.
!
!> @details
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments, cleaning
!!  - tvar_oa modified alloc dim., bug with tvc_oa alloc. dim nzc_oa <- nzv_oa
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!------------------------------------------------------------------------------

      subroutine upd_init_oa(     &
       imin, imax                 &
      ,jmin, jmax                 &
      ,kmin, kmax                 )


      use module_oa_variables
      use module_parameter_oa, only : allocate_3D_glob_array_oa_cplx &
                                     ,allocate_2D_glob_array_oa_cplx
      use module_oa_periode

      !use module_oa_time
      !use module_oa_space
      !use module_oa_stock

      implicit none

      !> Grid index range for the analysis of fields passed in argument to the OA module 
      !! Fields are analysed inside grid point domain [imin,imax]x[jmin,jmax]x[kmin,kmax]:
      integer, intent(in) :: &
       imin, imax            & 
      ,jmin, jmax            &
      ,kmin, kmax

      integer                                                         &
       ic_u                                                           &
      ,iv_u                                                           &
      ,nz_u                                                           &
      ,i                                                              &
      ,j                                                              &
      ,l

      integer                                                         &
       nzc_u

      nzupd3d_oa = 0
      nzupd2d_oa = 0

      tvar_oa(1:maxtyp_oa,1:max_idcfg_oa,1:nmvar_oa) = 0.

      set_tvar_oa_loop_config : do ic_u = 1 , nzc_oa

       nzc_u = 1
!      nzc_u = 0
!100  continue
!STDALONE : search possible configuration of the same type ie, having the same config. code tc_oa
     do while (tvar_oa(tc_oa(ic_u),nzc_u,1).ne.0)
       nzc_u = nzc_u + 1
     enddo
!     if (tvar_oa(tc_oa(ic_u),nzc_u,1).ne.0) goto 100
!STDALONE tvc_oa : if several configuration of the same type are requested tvc_oa ordinates them.
!#BLXD CORRECTED BUG tvc_oa should be allocated to nzc_oa not to nzv_oa
      tvc_oa (ic_u) = nzc_u
      var_in_cfg_loop : do iv_u = begc_oa(ic_u),begc_oa(ic_u+1)-1
      if (updv_oa(iv_u).eq.2) then ! #BLXD "variable remise a jour interactivement" what if UPDV_W=1 ?
        if (tgv3d_oa(tv_oa(iv_u)).eq.3) then
          nzupd3d_oa    = nzupd3d_oa + 1
          tupd_oa(iv_u) = nzupd3d_oa 
          tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd3d_oa
        else
          nzupd2d_oa    = nzupd2d_oa + 1
          tupd_oa(iv_u) = nzupd2d_oa 
          tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd2d_oa
        endif
      ! #BLXD 3rd dimension of tvar_oa is the number of variable per configuration
      !       in simple case it is always 1, but for composite var... begc_oa(ic_u+1)-begc_oa(ic_u)-1
      ! ( tc_oa(ic_u), tvc_oa(ic_u), 1 )
      endif
      enddo var_in_cfg_loop

      enddo set_tvar_oa_loop_config

    !STDALONE !.....ajout d'une variable pour le calcul de l'energie:
    !STDALONE       if (tc_oa(ic_u).eq.112.and.updv_oa(begc_oa(ic_u)).eq.2) then
    !STDALONE           nzupd3d_oa    = nzupd3d_oa + 1
    !STDALONE           tvar_oa(tc_oa(ic_u),nzc_u,begc_oa(ic_u+1)-begc_oa(ic_u)+1) = nzupd3d_oa
    !STDALONE !         nzupd2d_oa    = nzupd2d_oa + 1
    !STDALONE !         tvar_oa(tc_oa(ic_u),nzc_u,begc_oa(ic_u+1)-begc_oa(ic_u)+2) = nzupd2d_oa
    !STDALONE       endif
    !STDALONE       if (tc_oa(ic_u).eq.110.and.updv_oa(begc_oa(ic_u)).eq.2) then
    !STDALONE           nzupd3d_oa    = nzupd3d_oa + 1
    !STDALONE           tvar_oa(tc_oa(ic_u),nzc_u,begc_oa(ic_u+1)-begc_oa(ic_u)+1) = nzupd3d_oa
    !STDALONE       endif
    !STDALONE       enddo  

!.....allocations dynamiques:

      !STDALONE change allocation limits
      !STDALONE allocate (                                                      &
      !STDALONE  var2d_oa(0:imax+1,0:jmax+1     ,nzupd2d_oa)                    &
      !STDALONE ,var3d_oa(0:imax+1,0:jmax+1,0:kmax+1,nzupd3d_oa)                &
      !STDALONE          )

      ! BLXD dans la version "tile-thread compliant" les variables var3d et var2d 
      !       ont pour dimension GLOBAL_3D_ARRAY et GLOBAL_2D_ARRAY (voir module_parameter_oa)
      !       Ces variables sont traites par le module OA  par sous-sous domaines 
      !       horizontaux compatibles avec les boucles croco en "tile-threads"
      !       (double parallelisation possible pour les sous-domaines horizontaux)
      !       c-a-d que var2/3d_oa ont les dimensions des PRIVATE_SCRATCH_ARRAY.       
      !       Actuellement, on a imin,imax=Istr-2:Iend+2, jmin,jmax=Jstr-2,Jend+2
      !       Un ajustement est possible via ip_ran et I_PRI_RANGE, J_PRI_RANGE,..etc 
      !       (point interieur seulement avec ou sans les bdy points)
      !       voir online_spectral_diags.F
      !       Le traitement fait par var_space_oa, qui definit le pointeur 
      !       de domaine spatial begvs_oa  
      !       (structure vectoriel qui empile/desempile les analyses OA)
      !       est sensible a la cle cppc MASKING. 
      !       La conversion index de pointeur begsv_oa vers indices de domaine horizontal
      !       se fait via les variables ls2i_oa/l2sj_oa qui a actuellement pour
      !       domaine d'application {Istr-2:Iend+2, Jstr-2:Jend+2}.
      !       Au final, les points (i,j) mis a jour avec les analyses OA 
      !       et restitue dans var3d_oa/var2d_oa sont sensibles a la strategie adoptee 
      !       (bdy, points interieurs, points de frontiere ouverte,...)  
      !       Si var3d_oa, var2d_oa sont uniquement des variables diagnostics
      !       on peut tout de meme souhaiter avoir les analyses aus points frontieres 
      !       (forçage)
      !       Si var3d_oa, var2d_oa deviennent des variables actives 
      !       avec par exemple des derivees spatiales par exemple, il faudra envisager
      !       les echanges MPI necessaires sur les ghost points et faire attention
      !       a la compatibilite avec la cle cpp THREE_GHOST_POINTS et les schema d'ordre eleve
      !
      !allocate (                                                              &
      ! var2d_oa( imin-1:imax+1, jmin-1:jmax+1, nzupd2d_oa )                   &
      !,var3d_oa( imin-1:imax+1, jmin-1:jmax+1, kmin:kmax, nzupd3d_oa )        &

      call allocate_3D_glob_array_oa_cplx( var3d_oa, nzupd3d_oa )
      call allocate_2D_glob_array_oa_cplx( var2d_oa, nzupd2d_oa )

      var3d_oa(:,:,:,:) = ( 0.D0, 0.D0 )
      var2d_oa(:,:,:)   = ( 0.D0, 0.D0 )

!.....History file:

      call history_oa(11,-1,-1,-1,-1, -1)

      return
      end subroutine upd_init_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note 
!
! DESCRIPTION: 
!
!> @brief
!
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors 
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon : 
!!  - Memory oprimization when tracking isopycnal levels (lev_init_oa,..)
!!      - ifl_l flag eliminated, replaced by counting analysis of type 20 in nzlevel_oa, 
!!        and testing if nzlevel_oa>0.
!!      - enables to diminish the size of the structured type array wlev, now sized 
!!        to nzlevel_oa intead of nzv_oa (the total number of OA analysis requested 
!!        in the simulation).
!!  - stand-alone version : 
!!      - reading isopycne data in a file not maintained. ut only the possibility
!!      - isopycne values can be picked along the rhp_t field profile (according to namelist parameters).
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo
!! - reading isopycne targets in a parameter file ?
!------------------------------------------------------------------------------
      subroutine  lev_init_oa( tile, rhp_t, rhp_t_lbound, rhp_t_ubound )

      use module_tile_oa, only : st, pst
      use module_oa_variables
      !use module_oa_time
      use module_oa_space
      !use module_oa_periode
      !use module_oa_stock
      use module_oa_level

      implicit none

      !> Lower and upper index bounds of the 3-Dimensional density array
      integer, dimension(3), intent(in) ::  &
        rhp_t_lbound, rhp_t_ubound

      ! 3-Dimensional density array passed in arguments to the OA module
      double precision, dimension(rhp_t_lbound(1):rhp_t_ubound(1),rhp_t_lbound(2):rhp_t_ubound(2),rhp_t_lbound(3):rhp_t_ubound(3)) :: &
        rhp_t  !< Density array at t-grid point

      ! tile number
      integer, intent(in) :: tile 

      integer                                                         &             
       i_l                                                            &
      ,j_l                                                            &
      ,k_l                                                            &
      ,ls_l                                                           &
      ,ls1_l                                                          &
      ,iv_l                                                           &
      ,dim_l                                                          &
      ,izlevel_oa

   !STDALONE   double precision                                                &
   !STDALONE    val_l(1:kmax+1)
   !STDALONE
   !STDALONE   integer ifl_l
   !STDALONE
   !STDALONE   character*30 file_l

!.....Counts of the number of requested isopycne analysis
!     Test nzlevel_oa>0 replaces ifl_l flag
!     This loop is applied to reduce the size of the wlev structure allocation to nzlevel_oa
      nzlevel_oa = 0
      do iv_l = 1 , nzv_oa
         if (tv_oa(iv_l).eq.20) then
            nzlevel_oa = nzlevel_oa + 1
            !ifl_l = 1
         endif
      enddo

!.....si variables de type 20:

   !STDALONE if (ifl_l.eq.1) then
   ! BLXD identical to flag isopycne_analysis (true if at least one variable as type tv_oa(iv_l)==20)
      
      isopycne_requested : if ( nzlevel_oa>0 ) then

!.... Pointer to tile structure
          pst => st(tile)

   !STDALONE The array of structured type wlev is sized to nzlevel_oa intead of
   !STDALONE of nzv_oa (the total number of OA analysis requested in trhe simulation)

          call allocate_lev_part1_oa

          izlevel_oa = 0

    !STDALONE To read the target density in a file, adpat the following lines to your directories
    !STDALONE !.....lecture du fichier contenant les caracteristiques des isopycnes
    !STDALONE !     a suivre:
    !STDALONE       open(unit=10,file='DATA/isopycnes.in')
    !STDALONE       file_l ='DATA/'//dom_c//'isopycnes.dat'
    !STDALONE       open(unit=11,file=file_l)
    !STDALONE         iv_l = 1
    !STDALONE         ls_l  =begvs_oa(iv_l)
    !STDALONE         do ls1_l=begvs3d_oa(ls_l),begvs3d_oa(ls_l+1)-1
    !STDALONE           read (10,*) val_l(ls1_l-begvs3d_oa(ls_l)+1)
    !STDALONE           write(11,*) val_l(ls1_l-begvs3d_oa(ls_l)+1)
    !STDALONE         enddo
    !STDALONE       close(10)
    !STDALONE       close(11)

!.....construction de la structure de variables necessaire
!     au suivi des isopycnes:


      variable_loop : do iv_l = 1 , nzv_oa

        if (tv_oa(iv_l).eq.20) then

   !STDALONE allocation is performed earlier with the derived type array wlev which is sized to nzlevel_oa
   !STDALONE instead of being sized to nzv_oa (the total number of OA analysis requested in trhe simulation)        
   !STDALONE if (izlevel_oa.eq.0) call allocate_lev_part1_oa

            izlevel_oa            = izlevel_oa + 1
            lev2v_oa (izlevel_oa) = iv_l
            v2lev_oa (iv_l)       = izlevel_oa
            dim_l                 = pst%begvs3d_oa(pst%begvs_oa(iv_l+1)) &
                                  - pst%begvs3d_oa(pst%begvs_oa(iv_l))
            call allocate_lev_part2_oa(                               &
              izlevel_oa                                              &
             ,iv_l                                                    &
             ,dim_l )
  !           ,begvs3d_oa(begvs_oa(iv_l+1))-begvs3d_oa(begvs_oa(iv_l)) )
            do ls_l = pst%begvs_oa(iv_l),pst%begvs_oa(iv_l+1)-1
             i_l = pst%l2i_oa(ls_l)
             j_l = pst%l2j_oa(ls_l)
             do ls1_l=pst%begvs3d_oa(ls_l),pst%begvs3d_oa(ls_l+1)-1
              k_l = pst%kmin3d_oa(ls_l) + (ls1_l - pst%begvs3d_oa(ls_l))* dk_oa(iv_l)
              wlev_oa(izlevel_oa)%rhp(ls1_l - pst%begvs3d_oa(pst%begvs_oa(iv_l)) + 1 ) &
                  = rhp_t(i_l,j_l,k_l)
! BLXD TODO optimize successive precision conversion from module_oa_type where rhp and Z are real precision, here real(8)
              !if(if_print_node) write(io_unit2,*) i_l, j_l, k_l, rhp_t(i_l,j_l,k_l)
   !STDALONE                                                              = val_l(ls1_l - begvs3d_oa(ls_l)+1)
             enddo
            enddo
         endif

       enddo variable_loop

!.....Specific module_tile_oa pointer nullified
      pst => null()


      endif isopycne_requested

      return
      end subroutine  lev_init_oa

      subroutine allocate_lev_part1_oa 

      use module_oa_level

      implicit none

      ! BLXD reduced allocated size to nzlevel_oa instead of nzv_oa
      !      see lev_init_oa modifications.
      !STDALONE allocate(wlev_oa (nzv_oa))
      !STDALONE allocate(lev2v_oa(nzv_oa))
      !STDALONE allocate(v2lev_oa(nzv_oa))
      allocate(wlev_oa (nzlevel_oa))
      allocate(lev2v_oa(nzlevel_oa))
      allocate(v2lev_oa(nzlevel_oa))

      return
      end subroutine allocate_lev_part1_oa 

      subroutine allocate_lev_part2_oa (l_a,lv_a,dim_a)

      !use module_oa_variables
      !use module_oa_time
      !use module_oa_space
      !use module_oa_periode
      !use module_oa_stock
      use module_oa_level

      implicit none

      integer, intent(in) :: dim_a
      integer, intent(in) :: l_a, lv_a

      !STDANLONE BLXD l_a now ranges from 1 to nzlevel_oa
      allocate( wlev_oa(l_a)%z  (dim_a) )
      allocate( wlev_oa(l_a)%rhp(dim_a) )
      allocate( wlev_oa(l_a)%k  (dim_a) )

!.....History file:
      call history_oa(1,l_a,lv_a,dim_a,-1, -1)

      wlev_oa(l_a)%k(:) = 0
      wlev_oa(l_a)%z(:) = 0.

      return
      end subroutine allocate_lev_part2_oa


!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
! subroutine appelee par le programme principal
! pour chaque variable ("configuration") main_oa teste
! si l'on se trouve dans une periode de "calcul" et si tel
! est le cas il gere le calcul du coef...
! outputs: wf_oa: variable contenant tous les coefs...
!
!> @brief performs the time convolution between numerical fields and 
!! the requested "atoms" (i.e., wavelets, Fourier, Dirac). 
!
!> @details For each OA variariable-configuration (fields, scale), at each time step 
!! cumulates the convolution product between the evolving numerical field 
!! and the requested "atom" (i.e., wavelets, Fourier, Dirac).
!! The resulting coefficient is stored in the wf_oa structure.
!!
!! In the case of Fourier or wavelets, two types of analysis are available:
!! - dori_oa set to 1 : no scale integration, analysis at a discrete given scale 
!! - dori_oa set to 2 : scale integration
!!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments, cleaning, intent and optional attributes
!!  - Croco-OA interface 1st version (2020)
!!  - Croco Tile-thread compliant version (2021)
!> @date 2021
!
!> @todo BLXD Modify Croco-OnlineA module interface.
!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!   Interface = mix between the Stand-alone and Original OnlineA versions 
!!   => Croco arrays are passed as arguments to the OnlineA routines with 
!!   reduced dimension range. Check memory duplication.
!! @todo BLXD specific module to set consitent precisions 
!!       with real(KIND=wp)... etc 
!!       eg, var_oa function defined as real => real(8)
!!       eg, complex =>
!!       double precision is an obsolete type TO REPLACE
!! @todo BLXD introduce tiny instead of 1.e-30 ?
!----------------------------------------------------------------------

      subroutine main_oa(  tile                       &
                          ,ichoix                     &
                          ,ivar_m                     &
                          ,iic_oa                     &
                          ,imin, imax                 &
                          ,jmin, jmax                 &
                          ,kmin, kmax                 &
                          ,rhp_t                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
                          ,depth_t                    ) 


      use module_tile_oa, only : pst, st, deallocate_tile_varspace_oa
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock
!      use module_oa_upd
!      use scalars

      implicit none

      integer, intent(in) :: ichoix

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      !> To handle specific calls to main_oa with variables updated at different position in the calling code
      integer, intent(in), optional :: ivar_m
     
      !> Grid index range for the analysis of fields passed in argument to the OA module 
      !! Grid index are Roms-Croco tile(-thread) indices.
      !! If no tile(-thread) grid index are MPI subdomain index
      !! If no tile(-thread) and no MPI subdomain decomposition => full domain grid range
      integer, intent(in) ::       &
             imin, imax            & 
            ,jmin, jmax            &
            ,kmin, kmax

      !> tile parameter from croco
      integer, intent(in) :: tile

      !> 3-Dimensional density array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: rhp_t
                                                                          
      !> 3-Dimensional depth array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: depth_t

      !> Lower and upper index bounds of the 3-Dimensional density array (see if pst%imin,...)
      integer, dimension(3)  ::  rhp_t_lbound, rhp_t_ubound

      !complex :: psi_oa
      !real    :: psi_p2s_oa

      real    :: var_oa
      
      !TODO declare flag here if only needed here
      logical :: ifl_test_composite

      integer                                                         &
         iv_m                                                         &
        ,ic_m                                                         &
        ,lt_m                                                         &
        ,ls_m                                                         &
        ,ls1_m

      integer                                                         &
         i_m                                                          &
        ,j_m                                                          &
        ,k_m                                                          &
        ,kpt_m                                                        &
        ,lp_m                                                         &
        ,la_m
    
! BLXD TODO : what precision ? 
      complex                                                         &
         psi_m

! BLXD var_oa function real(4) ?
      real :: tmpvar

! BLXD introduce tiny instead of 1.e-30

      integer :: io_nodoa

!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

!     if (iic_oa.eq.kount0.and.ifl_init_oa.eq.0) then

!STDALONE      if (ifl_init_oa.eq.0) then
!STDALONE         call initial_oa
!STDALONE         ifl_init_oa=1
!STDALONE      endif


      oa_analysis_requested : if (nzv_oa.ne.0) then

!.....Pointer association to specific module_tile_oa
      pst => st(tile)

!.....flag variable composite 
      ifl_test_composite = .true.      

!.....test variable
       if (ifl_test_oa==1) then
        call test_oa( iic_oa=iic_oa  & 
          ,imin=imin, imax=imax     &
          ,jmin=jmin, jmax=jmax     &
          ,kmin=kmin, kmax=kmax )
       endif

      !STDALONE if (flag_nrj_oa.eq.1) call nrj_upd_oa(1)

      if ( isopycne_analysis ) then
        rhp_t_lbound = (/imin,jmin,kmin/)
        rhp_t_ubound = (/imax,jmax,kmax/)
        call lev_upd_oa(     tile                                          &
                             ,rhp_t                                        &
                             ,depth_t                                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
                             ,rhp_t_lbound                                 &
                             ,rhp_t_ubound  )
      endif


      do ic_m = 1 , nzc_oa
      do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1
         !STDALONE to call main_oa at specific location for composite varaible confguration
         if (present(ivar_m)) then
            ifl_test_composite = tv_oa(iv_m).eq.ivar_m
         endif

!........ce test est utilise pour les appels cibles lors des calculs energetiques:
         if_energy_calc : if (ichoix.ge.0.or.ifl_test_composite) then
!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

     if_discrete_per : if (dori_oa(iv_m).eq.1) then

      time_of_analysis_loop : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1


!.....test pour le calcul:
       tstp_match_convol_window : if ( iic_oa.ge.kountv_oa(1,lt_m) .and. iic_oa.le.kountv_oa(2,lt_m)                               &
        &  .and.               mod( int( (iic_oa-kountv_oa(1,lt_m)), kind=8 ), resv_oa(per_t2p_oa(lt_m)) ).eq.0            &
        &  .and.               ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 )                                                 &
                             ) then

! #BLXD ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 ) = calling main_oa for particular variable

!........allocation if necessary (if variable isn't allocated yet, it is done.., at user-specified kount):



          if (tallocated_oa(lt_m).eq.-1) then
#ifdef OA_TRACES
            if (verbose_oa>=10) then
            if(if_print_node) write (io_unit,*) 'MAIN_OA : alloc_win_oa',ic_m,iv_m,lt_m,mynode
            endif
#endif
            call allocate_win_oa (                                                             &
                    lt_m, ic_m, iv_m,                                                          &
                    pst%begvs3d_oa(pst%begvs_oa(iv_m+1)) - pst%begvs3d_oa(pst%begvs_oa(iv_m)), &
                    iic_oa )
            !lt_m,ic_m,iv_m,begvs3d_oa(begvs_oa(iv_m+1)) - begvs3d_oa(begvs_oa(iv_m)) )
          endif



!........precalcul de l'indice temporel prenant en compte la translation (u):

         kpt_m =  iic_oa - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )

!........precalcul de psi:

          psi_m = psi_oa(                                                     &
                tpsi_oa(iv_m)                                                 &
               ,psi_p2s_oa(                                                   &
                                 tpsi_oa(iv_m)                                &         
                                ,perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa )    &                                     
               ,real(kpt_m)* dti                                              &
               ,dti*resv_oa(per_t2p_oa(lt_m))                                 &           
               ,fb_oa,fc_oa  ) 
  
!........boucle spatiale:


           horizontal_space_loop : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
            i_m = pst%l2i_oa(ls_m)
            j_m = pst%l2j_oa(ls_m)

            vertical_space_loop : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
            k_m = pst%kmin3d_oa(ls_m) + (ls1_m - pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
            tmpvar =  var_oa( tv_oa(iv_m)                                        &
                                   ,ichoix                                       &
                                   ,i_m                                          &
                                   ,j_m                                          &
                                   ,k_m                                          &
                                   ,iv_m                                         &
                                   ,ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1       ) 

             wf_oa(tallocated_oa(lt_m))%coef (                                   &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )   =           &
             wf_oa(tallocated_oa(lt_m))%coef (                                   &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )               &
                     + conjg(psi_m) * tmpvar                                     &
                     / max(perv_oa(2,per_t2p_oa(lt_m)),1.e-30)

            enddo vertical_space_loop
           enddo  horizontal_space_loop


       endif tstp_match_convol_window
      enddo time_of_analysis_loop
 
!.....periodes integrees:
!---------------------------------------------------------------------------
      elseif (dori_oa(iv_m).eq.2) then if_discrete_per
!---------------------------------------------------------------------------

      time_of_analysis_loop2 : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1

!.....test pour le calcul:

       tstp_match_convol_window2 : if ( iic_oa.ge.kountv_oa(1,lt_m).and.                                     &
                               iic_oa.le.kountv_oa(2,lt_m).and.                                      &
                               mod( int((iic_oa-kountv_oa(1,lt_m)),kind=8) ,                         &
                                                                 resv_oa(per_t2p_oa(lt_m))).eq.0          &
                                     .and.                                                                &
                                     (ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)                              &
                             ) then

!........allocation si necessaire:

         if (tallocated_oa(lt_m).eq.-1) then
            call allocate_win_oa (                                                             &
                    lt_m, ic_m, iv_m,                                                          &
                    pst%begvs3d_oa(pst%begvs_oa(iv_m+1)) - pst%begvs3d_oa(pst%begvs_oa(iv_m)), &
                    iic_oa )
         endif

!........precalcul de l'indice temporel prenant en compte la translation (u):

         kpt_m =  iic_oa  - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )
     
!........boucle sur toutes les periodes:

         period_of_analysis_loop2 : do lp_m = begvp_oa(iv_m) , begvp_oa(iv_m+1) - 1

!........precalcul de psi:
         
          psi_m = psi_oa(                                                       &
                tpsi_oa(iv_m)                                                   &
               ,psi_p2s_oa(                                                     &
                                 tpsi_oa(iv_m)                                  &
                                ,perv_oa(1,lp_m),fb_oa,fc_oa )                  &
               ,real(kpt_m)* dti                                             &
               ,dti*resv_oa(per_t2p_oa(lt_m))                                &
               ,fb_oa,fc_oa   ) 

!........boucle spatiale:
           horizontal_space_loop2 : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
            i_m = pst%l2i_oa(ls_m)
            j_m = pst%l2j_oa(ls_m)

            vertical_space_loop2 : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
            k_m = pst%kmin3d_oa(ls_m) + (ls1_m-pst%begvs3d_oa(ls_m))* dk_oa(iv_m)

            tmpvar = var_oa( tv_oa(iv_m)                                        &
                                   ,ichoix                                      &
                                   ,i_m                                         &
                                   ,j_m                                         &
                                   ,k_m                                         &
                                   ,iv_m                                        &
                                   ,ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )

             wf_oa(tallocated_oa(lt_m))%coef (                                  &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 ) =            &
             wf_oa(tallocated_oa(lt_m))%coef (                                  &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )              &
                     + conjg(psi_m) * tmpvar                                    &
                     / max(perv_oa(2,lp_m),1.e-30)

            enddo vertical_space_loop2
           enddo horizontal_space_loop2
         enddo period_of_analysis_loop2
       endif tstp_match_convol_window2

      enddo time_of_analysis_loop2
      endif if_discrete_per

!.....enddo associe a iv_m et ic_m:
      endif if_energy_calc
      enddo
      enddo

!.....eventuelle desalocation:

      do la_m = 1,nmsimult_oa
         if ( associated(wf_oa(la_m)%coef) ) then

         !if ((kountv_oa(2,wf_oa(la_m)%t_indice).le.iic_oa.or.iic_oa.eq.nt_max-1)      &
         if ((kountv_oa(2,wf_oa(la_m)%t_indice).le.iic_oa.or.iic_oa.eq.nt_max)      &

           .and.(ichoix.eq.cnb_oa(wf_oa(la_m)%variable).or.ichoix.lt.0)         &
           .and.(ichoix.eq.des_oa(wf_oa(la_m)%variable)                          &
                 .or.des_oa(wf_oa(la_m)%variable).eq.0)                          &
                 ) then                

            lt_m = wf_oa(la_m)%t_indice
            iv_m = wf_oa(la_m)%variable
            ic_m = wf_oa(la_m)%config
#ifdef OA_TRACES
            if(verbose_oa>=10) then
            if(if_print_node) write (io_unit,*) 'MAIN_OA : dealloc_win_oa',ic_m,iv_m,lt_m,mynode
            endif
#endif
            call subsave_oa ( tile, la_m )                           
            call deallocate_win_oa(lt_m)
         endif
         endif
      enddo

!.....Specific module_tile_oa pointer nullified
      pst => null()

!.....If end of simulation Deallocate deallocate_tile_varspace_oa
!     consistent with the deallocate_tile_varspace_oa in online_spectral_diags
!     after the $OMP BARRIER
      if (iic_oa.eq.nt_max) call deallocate_tile_varspace_oa(tile)


      endif oa_analysis_requested

      return
      
      end subroutine main_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: sauvegarde des coefficients qui viennent d'etre calcules.
!
!> @brief test variable is an homogeneous LC of cosine functions
!
!> @details
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - 2020-2021 Croco version. 
!!  - 2014 namelists, optimization, Stand-alone version
!> @date 2021
!> @todo
!----------------------------------------------------------------------

      subroutine subsave_oa( tile, la_s )
 
      use module_oa_variables
      use module_oa_stock

      implicit none

      integer, intent(in) :: tile, la_s                                   
      integer ::                                                      &
            iv_s                                                      &
           ,lt_s                                                      &
           ,ic_s                                                       

      ic_s = wf_oa(la_s)%config
      iv_s = wf_oa(la_s)%variable
      lt_s = wf_oa(la_s)%t_indice
      
      if (updv_oa(iv_s).ne.0) then
       call var_upd_oa ( tile, ic_s,iv_s,lt_s )
      endif

      return
      end subroutine subsave_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
! DESCRIPTION: 
!
! DESCRIPTION: 
!
!> @brief Updates output arrays var2d, var3d with the analysis when available.
!
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments,  cleaning
!!  - intent in/out specification, changing if statements since swt_wfpf_oa has a single value
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine var_upd_oa( tile                                     &
      ,ic_u                                                           &
      ,iv_u                                                           &
      ,lt_u                                                           &
                                )

      use module_tile_oa, only : st, pst
      use module_oa_variables !updv_oa, kmin3d_oa, begvs3d_oa, begvsoa, l2i_oa, l2j_oa, tupd_oa
      use module_oa_time
      use module_oa_space     !tgv3d_oa, dk_oa
   
      use module_oa_periode   !swt_wfpf_oa
      use module_oa_stock     !wf_oa, tallocated_oa
!      use module_oa_level
!      use module_oa_upd      !var3d_oa
!      use scalars, only : iminmpi, jminmpi
      use module_parameter_oa, only : iminmpi, jminmpi

      implicit none

      integer, intent(in) :: tile, ic_u, iv_u, lt_u 

      integer    ::                                                   &
       i_u                                                            &
      ,j_u                                                            &
      ,k_u                                                            &
      ,ls_u                                                           &
      ,ls1_u                                                           
      
      integer, pointer :: pla_u => null(), ps1_u => null()
   
      integer :: io_nodoa, iu_glob, ju_glob, num

 
      variable_commune : if (iv_u.ne.-1.and.updv_oa(iv_u).eq.2) then
!---------------------------------------------------------------------------
!.....mise a jour des variables communes:
!---------------------------------------------------------------------------

!.....History file:

      call history_oa(14,lt_u,iv_u,-1,-1, -1)

!.....Local pointer to module_tile_oa array structure
       
      pst => st(tile) 

!.....Loacal pointer with iv_u, lt_u args in
      pla_u => tallocated_oa(lt_u)
      ps1_u => pst%begvs3d_oa(pst%begvs_oa(iv_u))

       var_number_of_dimensions : if (tgv3d_oa(tv_oa(iv_u)).eq.3) then                            ! three-dimensional variables

! BLXD WARNING ONLY var3d_oa  swt_wfpf_oa(iv_u)=4 TESTED 
!.......sortie de la partie reelle du coeff
        if (swt_wfpf_oa(iv_u).eq.1) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           !var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = real(wf_oa(tallocated_oa(lt_u))%coef (ls1_u-begvs3d_oa(begvs_oa(iv_u))+1 ) )
           var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = real( wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
          enddo
         enddo
        ! BLXD exclusive options if => else if
        !endif

!.......sortie du coef au carre:
        else if (swt_wfpf_oa(iv_u).eq.2) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           !var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = abs(wf_oa(tallocated_oa(lt_u))%coef (ls1_u-begvs3d_oa(begvs_oa(iv_u))+1 ) )**2
           var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = abs( wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )**2
          enddo
         enddo
        !endif

!.......sortie du valeur absolue du coef:
        else if (swt_wfpf_oa(iv_u).eq.3) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           !var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = abs(wf_oa(tallocated_oa(lt_u))%coef (ls1_u-begvs3d_oa(begvs_oa(iv_u))+1 ) )
           var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = abs(wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )

          enddo
         enddo
        !endif

!.......sortie du coef complexe:
        else if (swt_wfpf_oa(iv_u).eq.4) then

!            pla_u => tallocated_oa(lt_u)
!           pla_u => null()

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1


          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)


          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1

           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = wf_oa(pla_u)%coef( ls1_u-ps1_u+1 )

#if defined OA_TRACES
             if_verb3d : if (verbose_oa>=6) then
                 iu_glob = i_u + iminmpi-1 
                 ju_glob = j_u + jminmpi-1 
                 !if ( (lt_u==ltrec_oa(iv_u)) ) then
# ifdef IGW
                if ( iu_glob ==itgt_glob .and. (k_u==ktgt) ) then
# else
#  ifdef MILES
                if ( iu_glob==itgt_glob .and. ju_glob==jtgt_glob .and. k_u==ktgt ) then
#  endif
# endif
                     io_nodoa = set_io_nodoa(iv_u,mynode,4,3)! [3/4](odd/even-iv)000+mynode
                     ! old script iv_u =1 => 30000+mynode iv_u=2 => 40000+mynode
                     write (io_nodoa,fmt='(i5,i5,i3,2(1x,ES22.15E2))')i_u,j_u,k_u &
                         ,REAL(DBLE( var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) )),REAL(DIMAG( var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) ))
                end if
             endif if_verb3d
#endif
          enddo
         enddo

        endif

       else var_number_of_dimensions                                                          ! two-dimensional variables

!.......sortie de la partie reelle du coef
        if (swt_wfpf_oa(iv_u).eq.1) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           ! BLXD useless
           !k_u = kmin3d_oa(ls_u) + (ls1_u-begvs3d_oa(ls_u))* dk_oa(iv_u)
           var2d_oa(i_u,j_u,tupd_oa(iv_u)) = real(wf_oa(pla_u)%coef(ls1_u-ps1_u+1 ) )
          enddo
         enddo
        ! BLXD exclusive options if => else if
        !endif

!.......sortie du coef au carre (energy):
        else if (swt_wfpf_oa(iv_u).eq.2) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           !k_u = kmin3d_oa(ls_u) + (ls1_u-begvs3d_oa(ls_u))* dk_oa(iv_u)
           var2d_oa(i_u,j_u,tupd_oa(iv_u)) = abs(wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )**2
          enddo
         enddo
        !endif

!.......sortie de la valeur absolue du coef: 
        else if (swt_wfpf_oa(iv_u).eq.3) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1  
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           !k_u = kmin3d_oa(ls_u) + (ls1_u-begvs3d_oa(ls_u))* dk_oa(iv_u)
           var2d_oa(i_u,j_u,tupd_oa(iv_u)) = abs(wf_oa(pla_u)%coef(ls1_u-ps1_u+1 ) )
          enddo
         enddo
        !endif

!.......sortie du coeff complexe:
        else if (swt_wfpf_oa(iv_u).eq.4) then
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1  
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           !k_u = kmin3d_oa(ls_u) + (ls1_u-begvs3d_oa(ls_u))* dk_oa(iv_u)
           var2d_oa(i_u,j_u,tupd_oa(iv_u)) = (wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )

#if defined OA_TRACES
             if_verb2d : if (verbose_oa>=6) then
                 iu_glob = i_u + iminmpi-1 
                 ju_glob = j_u + jminmpi-1 
                 !if ( (lt_u==ltrec_oa(iv_u)) ) then
                if ( iu_glob ==itgt_glob .and. (ju_glob==jtgt_glob) ) then
                     io_nodoa = set_io_nodoa(iv_u,mynode,4,3)! [3/4](odd/even-iv)000+mynode
                     ! old script iv_u =1 => 30000+mynode iv_u=2 => 40000+mynode
                     write (io_nodoa,fmt='(i5,i5,2(1x,ES22.15E2))')i_u,j_u &
                         ,REAL(DBLE( var2d_oa(i_u,j_u,tupd_oa(iv_u)) )),REAL(DIMAG( var2d_oa(i_u,j_u,tupd_oa(iv_u)) ))
                end if
             endif if_verb2d
#endif

          enddo
         enddo
        endif

        endif var_number_of_dimensions

!.....Nullify pointer to module_oa_tile struct.
       
      pst => null()

!.....Nullify local pointer
      pla_u => null()
      ps1_u => null()

      endif variable_commune

      return
      end subroutine  var_upd_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION:
!
!
!> @brief Allocates/deallocates array fields of derived type wf_oa
!
!> @details  wf_oa stores all the online analysis before dumping them
!            into var3d_oa, var2d_oa.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco version (2020) added pre-calculated size of wf_oa array
!!    i.e., max number of convolution window being simultaneously 
!!    under calculation
!> @date 2021
!> @todo BLXD general requirement to organize precision consistency
!!    eventhough compiling with double prec option
!!    Module with wp declaration and use of KIND(wp) ?
!----------------------------------------------------------------------

      subroutine allocate_win_oa (lti_a,lc_a,lv_a,dim_a,iic_oa)

      !use module_oa_time
      !use module_oa_space
      !use module_oa_periode
      use module_oa_stock , only    : wf_oa,           &
                                      tallocated_oa,   &
                                      nzw_oa

      use module_oa_time , only : nmsimult_oa   !< Pre calculated size for array wf_oa 


      implicit none

      integer, intent(in) :: iic_oa

      integer, intent(in) :: &
        lc_a                 &  !< index related to a specific configuration of analysis 
       ,lv_a                 &  !< index related to a specific variable to analyse
       ,lti_a                   !< index related to a specific time of analysis (i.e., time convolution window)

      integer, intent(in) :: dim_a !< Size of the spatial domain to analyse for variable lv_a 
      integer l_a

      ! Shall we keep track of imsimult_oa ?
      l_a=1
      do while ( (associated(wf_oa(l_a)%coef)).and.l_a.lt.nmsimult_oa)
         l_a = l_a + 1
      enddo

      ! BLXD TODO CLEAN we set nzw_oa to the pre-calc size of wf_oa array, ie nmsimult_oa
      ! ie, nzw_oa = nmsimult_oa 
      if (l_a.eq.nzw_oa ) then
         if(if_print_node) write (io_unit,*) 'Pre-calculated size for wf_oa now reached ',nzw_oa, nmsimult_oa
      else if (l_a.gt.nzw_oa) then
         if(if_print_node) write (io_unit,*) 'STOP : above pre-calculated size for wf_oa ',nzw_oa, nmsimult_oa
         stop
      endif

      if (l_a.eq.nmsimult_oa.and.associated(wf_oa(l_a)%coef) ) then
         if(if_print_node) write (io_unit,*) 'nmsimult_oa trop petit!', nmsimult_oa
         stop
      endif


      allocate( wf_oa(l_a)%coef(dim_a) )

!.....History file:
      call history_oa(2,lc_a,lv_a,l_a,-1,lti_a,iic_oa)

      wf_oa(l_a)%coef(:)               = (0.D0,0.D0) ! #BLDX double prec
      tallocated_oa(lti_a)             = l_a
      wf_oa(l_a)%t_indice              = lti_a 
      wf_oa(l_a)%config                = lc_a 
      wf_oa(l_a)%variable              = lv_a 

      return
      end subroutine allocate_win_oa


      subroutine deallocate_win_oa (lt_a)    

      use module_oa_stock , only : wf_oa, tallocated_oa
      use module_oa_time , only : nmsimult_oa

      implicit none

      integer, intent(in) :: lt_a

      
!.....History file:
      call history_oa(3,lt_a,-1,-1,-1, -1)

      deallocate (wf_oa(tallocated_oa(lt_a))%coef )

      nullify( wf_oa(tallocated_oa(lt_a))%coef )

      ! BLXD I put back the line where it was originally  
      ! changing declaration of dummy arg lt_a with option intent(in)
      ! Indeed deallocate_win_oa subroutine is called with effecive arg 
      ! being wf_oa(tallocated_oa(lt_a))%t_indice
      ! if %t_indice is set to -1 before other instructions
      ! lt_a is also put to -1 and it is the arg of tallocated array! 
      
      wf_oa(tallocated_oa(lt_a))%t_indice   = -1          ! 20070608 moved down a couple of lines..
      wf_oa(tallocated_oa(lt_a))%config     = -1
      wf_oa(tallocated_oa(lt_a))%variable   = -1 
      tallocated_oa(lt_a)                   = -1

!      wf_oa(tallocated_oa(lt_a))%t_indice   = -1

! for some reason:
! wf_oa(tallocated(lt_a))%t_indice was used in the call to this subroutine
! as a consequence, if that's changed, the subroutine can't use it any longer??
! anyhow, by moving this command to the end of the subroutine,
! the problem seems solved
      return
      end subroutine deallocate_win_oa 

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief returns the complex wavelet value.
!
!> @details 
!
!
! REVISION HISTORY:
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments,  cleaning
!!  - intent in/out specification, changing if statements since tpsi_p has a single value
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo BLXD
!   check double precision 
!------------------------------------------------------------------------------
      complex function psi_oa(                                        & 
                 tpsi_p                                               & 
                ,scale_p                                              & 
                ,t_p                                                  & 
                ,dti_p                                             & 
                ,fb_p                                                 & 
                ,fc_p             )


      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

!....attention il faut encore tenir compte des normalisations...

      integer, intent(in) :: tpsi_p

      double precision, intent(in) :: dti_p, t_p

      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p

      real ::                                                         & 
        w0_p                                                          & 
       ,rec_p 

!STDALONE      pi_p = acos(-1.)

! attention exemple de psi incomplet pour windowed fourier...
! les coefs pour la reconstruction n'ont
! pas ete ajoutes.

!.....dirac:

      if (tpsi_p .eq. 0 ) then      
       psi_oa = 1.

!.....ondelette de morlet:
  
      else if (tpsi_p .eq. 1 ) then      
       psi_oa =                                                        & 
          exp( - ( t_p /scale_p ) ** 2 / fb_p )                        & 
        * exp( t_p / scale_p * fc_p*(2.*pi_oa) * (0.,1.) )             & 
!       / (2*pi_oa)**0.5                                               & 
        / (pi_oa)**0.25                                                & 
!       * sqrt( dti_p/ scale_p )    
        / sqrt( scale_p )    

!....."windowed fourier":
      else if (tpsi_p .eq. 2.or.tpsi_p.eq.3) then
       psi_oa =  exp(  2.* pi_oa* ( t_p - dti_p ) / scale_p * (0.,1.) )
      endif 

      end function psi_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
! DESCRIPTION: 
!
!
! DESCRIPTION: 
!
!> @brief Transforms the requested time period of the analysis into the wavelet scale (s).
!
!> @details transformation periode --> echelle (s) pour l'ondelette choisie.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments,  cleaning
!!  - intent in/out specification, changing if statements since tpsi_p has a single value
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo 
!
!------------------------------------------------------------------------------

      real function psi_p2s_oa(                                       &
                 tpsi_p                                               &
                ,per_p                                                &
                ,fb_p                                                 &
                ,fc_p                )

      use module_oa_variables
      implicit none

      integer, intent(in) :: tpsi_p

      real, intent(in)    ::                                          &
       per_p                                                          &
      ,fb_p                                                           &
      ,fc_p

!STDALONE      ,pi_p 

!STDALONE       pi_p = acos(-1.)

!.....dirac:

      if (tpsi_p .eq. 0 ) then
       psi_p2s_oa = 0.

!.....ondelette de morlet:

      else if (tpsi_p .eq. 1 ) then
       psi_p2s_oa = per_p / ( 4. * pi_oa ) * (fc_p*2.*pi_oa + sqrt(2.+(fc_p*2.*pi_oa)**2))  
                                                 ! t/torrence central frequency

!....."windowed fourier":

      else if (tpsi_p .eq. 2.or.tpsi_p.eq.3 ) then
       psi_p2s_oa = per_p
      endif 

      return

      end function psi_p2s_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note Calculus must be adapted according to the wavelet shape (Morlet,...)
!
! DESCRIPTION: 
!
!> @brief In the case of analysis with "integrated time periods", calculates the 
!! wavelet scale resolution or Heisenberg box.
!
!> @details calcul de la taille des boites d'heisenberg et sauvegarde
!! attention reprendre calcul selon ondelette
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments,  cleaning
!!  - intent in/out specification, changing if statements since tpsi_p has a single value
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo BLXD
!! - exponential overflow underflow treatment!!!!!
!------------------------------------------------------------------------------

      subroutine box_oa

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock

      implicit none

      integer ib

      double precision tb(20001),wb(10001),w0b,npsi2,dtb,dwb
      double precision w0bb,dtbb,dwbb,t0b 

      w0b=fc_oa*2*pi_oa
      dtb=0.1
      dwb=0.1


      tb(1)=-1000.
      do ib=1,20000
       tb(ib+1)=tb(ib)+dtb
      enddo

      wb(1)=0.
      do ib=1,10000
       wb(ib+1)=wb(ib)+dwb
      enddo

      npsi2=0.
      do ib=1,20001
       npsi2=npsi2+exp(-tb(ib)*tb(ib)/fb_oa*2)
      enddo
      npsi2=npsi2*dtb/(pi_oa**0.5) !torrence

!---------------------------------------------------------------
!      npsi2=npsi2*dtb/(2.*pi_oa)  !matlab
!---------------------------------------------------------------

      t0b=0.
!---------------------------------------------------------------
! parametres (ivane)
!
!      do ib=1,20001
!      t0b=t0b+tb(ib)*exp(-tb(ib)*tb(ib)/fb_oa*2.)
!      enddo
!      t0b=sqrt(t0b*dtb/npsi2/(2.*pi_oa)) !matlab
!!      t0b=sqrt(t0b*dtb/npsi2/(pi_oa**0.5)) !torrence
!---------------------------------------------------------------

      w0bb=0.
      do ib=1,10001
      w0bb=w0bb+wb(ib)*exp(-(wb(ib)-w0b)*(wb(ib)-w0b)*fb_oa/2.)*fb_oa/2. ! torrence
      enddo
      w0bb=w0bb*dwb/npsi2/(pi_oa**0.5) !torrence

!---------------------------------------------------------------
!      w0bb=w0bb*dwb/npsi2/(2.*pi_oa) !matlab
!---------------------------------------------------------------

      dtbb=0.
      do ib=1,20001
      dtbb=dtbb+(tb(ib)-t0b)*(tb(ib)-t0b)*exp(-tb(ib)*tb(ib)/fb_oa*2.)
      enddo
      dtbb=sqrt(dtbb*dtb/npsi2/(pi_oa**0.5)) !torrence

!---------------------------------------------------------------
!      dtbb=sqrt(dtbb*dtb/npsi2/(2.*pi_oa)) !matlab
!---------------------------------------------------------------

      dwbb=0.
      do ib=1,10001
      dwbb=dwbb+(wb(ib)-w0bb)*(wb(ib)-w0bb)*exp(-(wb(ib)-w0b)*(wb(ib)-w0b)*fb_oa/2.)*fb_oa/2.   ! torrence
      enddo

      dwbb=sqrt(dwbb*dwb/npsi2/(pi_oa**0.5)) !torrence

!---------------------------------------------------------------
!      dwbb=sqrt(dwbb*dwb/npsi2/(2.*pi_oa)) !matlab
!---------------------------------------------------------------
 
      return
      end subroutine box_oa      

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: fonctions tests pour le developpement 
!
!> @brief test variable is an homogeneous LC of cosine functions
!
!> @details  the test variable is an homogeneous 3D field 
!! varying as a LC of cosine functions with periods and amplitude
!! set in the OA namelist 
!!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - 2020-2021 Croco version. 
!!  - 2014 namelists, doxygen comments, Stand-alone version, optimization
!> @date 2021
!> @todo
!----------------------------------------------------------------------

      subroutine test_init_oa(    &   
       imin, imax                 &
      ,jmin, jmax                 &
      ,kmin, kmax )

      use module_oa_variables, only : vardp_test_oa

      implicit none

      ! BLXD TODO check if using tile st(tile)%imin,... indices
      !> Grid index range for the analysis of fields passed in argument to the OA module 
      !! Croco tile(-thread) dimensions or MPI subdomain or full domain
      integer, intent(in) :: &
       imin, imax            & 
      ,jmin, jmax            &
      ,kmin, kmax
      
#ifdef OA_TRACES
      if(if_print_node) write (io_unit,*) '...Allocating vardp_test_oa'
      if(if_print_node) write (io_unit,*) imin,imax
      if(if_print_node) write (io_unit,*) jmin,jmax
      if(if_print_node) write (io_unit,*) kmin,kmax
#endif

      if ( .not. allocated(vardp_test_oa) ) then
        allocate (vardp_test_oa(imin:imax,jmin:jmax,kmin:kmax))
      end if

      return
      end subroutine test_init_oa

      subroutine test_oa(         &   
!       ichoix                     & 
       iic_oa                     & 
!      ,dti                        & 
      ,imin, imax                 &
      ,jmin, jmax                 &
      ,kmin, kmax )


      use module_oa_variables, only : vardp_test_oa

      implicit none

      ! BLXD TODO clean
      !!> Time integration step
      !!double precision, intent(in) :: dti

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      !> Grid index range might be retreived with tile structure st(tile)%imin
      integer, intent(in) :: &
       imin, imax            & 
      ,jmin, jmax            &
      ,kmin, kmax

      double precision  :: time_t, var_t
      integer           :: i, j, k, iper

       time_t = dti * real(iic_oa)
     
       var_t=0.D0 
       do iper = 1, nper_test
#ifdef OA_TRACES
         if(if_print_node) write(io_unit,*) 'OA test function (period,amplitude) ', iic_oa, period_test_oa(iper), amp_test_oa(iper)
#endif
         var_t = var_t + amp_test_oa(iper) * cos(2.d0*pi_oa/period_test_oa(iper)*time_t)
       end do

       ! BLXD removing temporary hardcoded option 
       !if (nper_test /= -1) then 
       do k=kmin,kmax
         do j=jmin,jmax
           do i=imin,imax !0,imax+1
            vardp_test_oa(i,j,k) = var_t
           enddo
         enddo
       enddo
       !end if


      return

      end subroutine test_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!> @brief Updates the isopycne localization (closest level and depth).
!
!> @details lev_upd_oa calls update_level_oa which searches for current
!! depth position of the target density (configuration code 20 requested in the OA namelist)
!! and stores the closest level and depth in the wlev_oa structured type array.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - 2020-2021 Croco version. 
!!  - 2014 namelists, doxygen comments, Stand-alone version, optimization
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine  lev_upd_oa( tile                  & 
      ,rhp_t                                        &
      ,depth_t                                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
      ,rhp_t_lbound                                 &
      ,rhp_t_ubound  )

      use module_tile_oa, only : pst, st
      use module_oa_space
      use module_oa_level

      implicit none

      !> Tile index
      integer, intent(in)  :: tile
 
      !> Lower and upper index bounds of the 3-Dimensional density array
      integer, dimension(3), intent(in) :: rhp_t_lbound, rhp_t_ubound

      !> 3D depth and density array at t-grid point
      double precision,                                                                                                  &
       dimension(rhp_t_lbound(1):rhp_t_ubound(1),rhp_t_lbound(2):rhp_t_ubound(2),rhp_t_lbound(3):rhp_t_ubound(3)), &
        intent(in) :: rhp_t, depth_t

      ! Local indices 
      integer ::                                                      &
       i_l                                                            &
      ,j_l                                                            &
      ,k_l                                                            &
      ,ls_l                                                           &
      ,ls1_l                                                          &
      ,iv_l                                                           &
      ,il_l

      integer, pointer :: ps1_l => null()

      pst => st(tile)

      do il_l = 1 , nzlevel_oa
         iv_l = lev2v_oa(il_l)
           do ls_l = pst%begvs_oa(iv_l),pst%begvs_oa(iv_l+1)-1
            i_l = pst%l2i_oa(ls_l)
            j_l = pst%l2j_oa(ls_l)
            ps1_l=> pst%begvs3d_oa(pst%begvs_oa(iv_l))
            do ls1_l = pst%begvs3d_oa(ls_l),pst%begvs3d_oa(ls_l+1)-1
              k_l = pst%kmin3d_oa(ls_l) + (ls1_l - pst%begvs3d_oa(ls_l))* dk_oa(iv_l)
              call update_level_oa(                                      &
                i_l                                                   &
               ,j_l                                                   &
               ,wlev_oa(il_l)%rhp(ls1_l-ps1_l+1) &
               ,wlev_oa(il_l)%k  (ls1_l-ps1_l+1) &
               ,wlev_oa(il_l)%z  (ls1_l-ps1_l+1) &
               ,pst%kmin, pst%kmax                           &
               ,rhp_t                                        &
               ,depth_t                                      &
               ,rhp_t_lbound                               &
               ,rhp_t_ubound  )

            enddo
           enddo

      enddo
      ps1_l => null()

      return
      end subroutine  lev_upd_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief
!
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments, cleaning
!!  - update_level_oa :
!!      - lower/upper k levels set to kmin,kmax instead of 1,kmax+1.
!!      - when the search fails, extrapolation treatment is removed (requires supplementary arguments).
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo update_level_oa
!! - extrapolation treatment?
!------------------------------------------------------------------------------

      subroutine update_level_oa (                  &
       i_c                                          &
      ,j_c                                          &
      ,rh_c                                         &
      ,k_c                                          &
      ,z_c                                          &
      ,kmin, kmax                                   &
      ,rhp_t                                        &
      ,depth_t                                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
      ,rhp_t_lbound                                 &
      ,rhp_t_ubound  )


      implicit none

      !> Vertical grid index range for the analysis of density field passed in argument to the OA module 
      integer, intent(in) :: kmin, kmax

      !> Lower and upper index bounds of the 3-Dimensional density array
      integer, dimension(3), intent(in) :: rhp_t_lbound, rhp_t_ubound

      ! BLXD TODO optimize successive precision conversion from module_oa_type where rhp and Z are real precision, here real(8)
      !> 3-Dimensional density array at t-grid point
      double precision, dimension(rhp_t_lbound(1):rhp_t_ubound(1),rhp_t_lbound(2):rhp_t_ubound(2),rhp_t_lbound(3):rhp_t_ubound(3)) &
                      , intent(in) :: rhp_t, depth_t


      integer                                                         &
       i_c                                                            &
      ,j_c                                                            &
      ,k_c                                                            &
      ,iflag

      real                                                            &
       z_c                                                            &
      ,rh_c                                                           &
      ,z2_c

      integer                                                         &
       l_c                                                            &
      ,dk1_c                                                          &
      ,dk2_c                                                          &
      ,nk_c                                                           &
      ,k2_c


!.....diverses initialisations:
      nk_c = kmax
      !STDALONE k_c   = min(nk_c,max(1,k_c)) 
      k_c   = min(nk_c,max(kmin,k_c)) 
      dk1_c = 0
      dk2_c = 0
      l_c   = 0
      iflag = 0
      k2_c  = 0

!.....boucle sur l'indice l_c:
!                              l_c impair: on cherche au dessus,
!                              l_c pair  : on cherche en dessous.
!
!
!
! la sortie de la boucle est le point delicat. elle est actuellement
! basee sur 3 tests:
!
!
! test 1: la boucle est executee au moins une fois,
! test 2: 1ere partie: la recherche a ete concluante (une profondeur a ete calculee)
!         2eme partie: on sort des bornes au dessus et au dessous (iflag=2)... abandon
! test 3: la premiere recherche "en dessous" n'a rien donne et sort des bornes (<1), on cherche au dessus...
!

      do while (                                                      &
        ( dk1_c.eq.0 .and. dk2_c.eq.0 ).or.                           &
!STDALONE        (  (rh_c-rhp_t(i_c,j_c,min(nk_c,max(1,k_c+dk1_c))))*(rh_c-rhp_t(i_c,j_c,min(nk_c,max(1,k_c+dk2_c)))).gt.0   &
        (  (rh_c-rhp_t(i_c,j_c,min(nk_c,max(kmin,k_c+dk1_c))))*(rh_c-rhp_t(i_c,j_c,min(nk_c,max(kmin,k_c+dk2_c)))).gt.0   &
         .and.iflag.ne.2                                              &
         .and.k_c+dk1_c.ge.1                                          &
         .and.k_c+dk2_c.ge.1                                          &
         .and.k_c+dk1_c.le.nk_c                                       &
         .and.k_c+dk1_c.le.nk_c )                                     &
         .or.                                                         &
         k_c+dk1_c.lt.1                                               &
         )

!......incrementation de l_c l'indice de boucle:
       l_c = l_c + 1

!......calcul des increments verticaux rã©els:
       if (mod(l_c,2).ne.0) then
          dk1_c = -l_c/2-1
          dk2_c = -l_c/2
       else
          dk1_c = l_c/2-1
          dk2_c = l_c/2
       endif

!......bornes verticales:
!STDALONE      if (k_c+dk1_c.ge.1.and.k_c+dk2_c.le.nk_c) then
     if (k_c+dk1_c.ge.kmin.and.k_c+dk2_c.le.nk_c) then

!.........remise a zero de iflag (utilise pour les sorties de domaines):
          iflag = 0

!........calcul de la profondeur dans le cas ou la densite n'est pas constante.
     if (rhp_t(i_c,j_c,k_c+dk2_c)-rhp_t(i_c,j_c,k_c+dk1_c).ne.0.) then
          z2_c = depth_t(i_c,j_c,k_c+dk1_c) +                        &
          ( depth_t(i_c,j_c,k_c+dk2_c)-depth_t(i_c,j_c,k_c+dk1_c) )  &
          * ( rh_c                    -rhp_t(i_c,j_c,k_c+dk1_c) )    &  
          / ( rhp_t(i_c,j_c,k_c+dk2_c)-rhp_t(i_c,j_c,k_c+dk1_c) )     
     else
      z2_c=0.
         endif
!........calcul du niveau le plus proche (utilise le cas echeant pour demarrer
!        la recherche au prochain pas de temps).
         if ( abs(rh_c-rhp_t(i_c,j_c,k_c+dk1_c)).le. abs(rh_c-rhp_t(i_c,j_c,k_c+dk2_c))   ) then
            k2_c = k_c + dk1_c
         else 
                k2_c = k_c + dk2_c
         endif
       else
!.......dans le cas ou l'on sort du domaine (par le haut ou par le bas), 
!       on increment iflag. lorsque iflag=2, i.e. l'on sort successivement
!       par le haut et par le bas, on arrete la recherche... pas de solution! 
        iflag = iflag + 1
       endif
!.....sortie de la boucle principale:
      enddo 
!.....si l'on a quitte la boucle precedente normalement (i.e. sans etre sorti
!     du domaine...) on valide la solution z (et l'on conserve le k le plus proche
!     pour gagner du temps lors de la prochaine recherche). 
!STDALONE      if (k_c+dk1_c.ge.1.and.k_c+dk2_c.le.nk_c) then
      if (k_c+dk1_c.ge.kmin.and.k_c+dk2_c.le.nk_c) then
        z_c = z2_c
        k_c = k2_c
      else
!.......extrapolation eventuelle: attention il s'agit d'une solution tres discutable...
!       mais qui permet de fournir une solution tres approximative.
        !STDALONE
        !TODO extrapolation?

        !if (rh_c.le.rhp_t(i_c,j_c,nk_c).and.rhp_t(i_c,j_c,nk_c-1) -rhp_t(i_c,j_c,nk_c).ne.0.) then
        ! z_c = min((hssh_w(i_c,j_c,1)-h_w(i_c,j_c))*0. ,                            &
        !   depth_t(i_c,j_c,nk_c) +                                  &
        ! ( depth_t(i_c,j_c,nk_c-1)-depth_t(i_c,j_c,nk_c) )         &
        ! * ( rh_c                  -rhp_t(i_c,j_c,nk_c) )            & 
        ! / ( rhp_t(i_c,j_c,nk_c-1) -rhp_t(i_c,j_c,nk_c) )            &
        !          )
        ! k_c = nk_c
        !elseif (rhp_t(i_c,j_c,1) -rhp_t(i_c,j_c,2).ne.0.) then
        ! z_c = max(-h_w(i_c,j_c),                                    &
        !   depth_t(i_c,j_c,2) +                                     &
        ! ( depth_t(i_c,j_c,1)-depth_t(i_c,j_c,2) )                 &
        ! * ( rh_c                  -rhp_t(i_c,j_c,2) )               &
        ! / ( rhp_t(i_c,j_c,1) -rhp_t(i_c,j_c,2) )                    &
        !          )
        ! k_c = 1
        !else
          z_c = 0.
          k_c = 0
        !endif
      endif
            
      return
      end subroutine update_level_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Handling OA module history output.
!
!> @details 
!
!
! REVISION HISTORY:
!!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo
!
!------------------------------------------------------------------------------

      subroutine history_oa ( ichoix, i1_h, i2_h, i3_h, i4_h, tile &
                              ,iic_oa, nt_max )

      use module_tile_oa, only : pst, st
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      use module_oa_stock
!      use module_oa_upd
! BLXD croco 2021 why scalars ?
!      use scalars
      implicit none
 
      integer, intent(in) :: tile
      integer, intent(in) :: ichoix
      integer, intent(in) :: i1_h, i2_h, i3_h, i4_h

      !> Current model integration iteration
      integer, intent(in), optional :: iic_oa 

      !> Last simulation iteration index
      integer, intent(in), optional :: nt_max                                        

      integer ::                                                     &
          k_m                                                        &
         ,ir_o                                                       &
         ,l_a                                                        &
         ,lv_a                                                       &
         ,lc_a                                                       &
         ,lti_a                                                      &
         ,ic_o                                                       &
         ,ic_u                                                       &
         ,lt_u                                                       &
         ,lt_o                                                       &
         ,lp_o                                                       &
         ,iv_o                                                       &
         ,la_s                                                       &
         ,lt_a                                                       &
         ,dim_a                                                      &
         ,iv_u                                                       &
         ,iv_s

       character(len=250) :: file_hist

#ifdef MPI
      if (.not. if_print_node) return
#endif

!---->Fichier de sortie:

      file_hist = trim(directory_out_oa) // txtslash // 'history_oa.dat'

      if (ichoix.eq.1) then
!*******************************************************************************
! allocate_lev_part2_oa
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
       !if (ifl_oa_out.eq.1) call sequential_begin()
       open(unit=io_hist,file=trim(file_hist),position='append')   
       write (io_hist,*) 'allocation: level n°',i1_h,',var n°',i2_h,',dim=',i3_h
       close(io_hist)
       !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.2) then
!*******************************************************************************
! allocate_win_oa  
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
       !if (ifl_oa_out.eq.1) call sequential_begin()
       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write (io_hist,*) 'allocation: fenetre n°',l_a,',iic_oa=',iic_oa       &
!                    ,',lc=',i1_h                                       &
!                    ,',lv=',i2_h                                       &
!                    ,',dim=',i3_h                                      &
!                    ,',lti=',i4_h
       write (io_hist,*) 'allocation: fenetre n°',i3_h,',iic_oa=',iic_oa       &
                    ,',lc=',i1_h                                       &
                    ,',lv=',i2_h                                       &
                    ,',dim=?'                                      &
                    ,',lti=',i4_h
       close(io_hist)
       !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.3) then
!*******************************************************************************
! deallocate_win_oa  
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write (io_hist,*) 'desallocation: fenetre n°',tallocated_oa(i1_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif


      if (ichoix.eq.4) then
!*******************************************************************************
! initial_oa
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*)
        write(io_hist,*) '*************************'
        write(io_hist,*) 'subroutine: initial_oa'
        write(io_hist,*) '*************************'
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      chx5 : if (ichoix.eq.5) then
!*******************************************************************************
! initial_oa   
!*******************************************************************************
!******************************************************************
!        caracteristiques d'une variable:      
!******************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()

         pst => st(tile)

         open(unit=io_hist,file=trim(file_hist),position='append')   
         loop_config : do ic_o=1,nzc_oa
            write(io_hist,*)
            write(io_hist,*)
            write(io_hist,*) '*************************************'
            write(io_hist,*) ' configuration de variables n°',ic_o
            write(io_hist,*) '*************************************'
            if (tc_oa(ic_o).ge.100) then
               write (io_hist,*) 'type de configuration:',tc_oa(ic_o)
               write (io_hist,*) 'WARNING ! composite configuration are hardcoded'
               write (io_hist,*) 'Composite configuration never tested with croco'
               write (io_hist,*) 'ERROR : check/change the configuration code, remove stop and recompile'
               stop 
            else
               write (io_hist,*) 'type de cfg-variable simple:',tc_oa(ic_o)
            endif

            loop_variable : do iv_o=begc_oa(ic_o),begc_oa(ic_o+1)-1
               write(io_hist,*)
               write(io_hist,*) '*************************************'
               write(io_hist,*) '        variable n°',iv_o
               write(io_hist,*) '*************************************'

               write(io_hist,*) '----------------------'
               write(io_hist,*) '---caracteristiques---'
               write(io_hist,*) '----------------------'
               write(io_hist,*) 'type            : ',tv_oa(iv_o) 
               if (swt_wfpf_oa(iv_o).eq.1) then
                  write(io_hist,*) 'calcul du coef. reel'
               endif 
               if (swt_wfpf_oa(iv_o).eq.2) then
                  write(io_hist,*) 'calcul du module au carre'
               endif 
               if (swt_wfpf_oa(iv_o).eq.3) then
                  write(io_hist,*) 'calcul du module'
               endif 
               if (swt_wfpf_oa(iv_o).eq.4) then
                  write(io_hist,*) 'calcul du coef. complexe'
               endif 
               if (save_oa(iv_o).eq.0) then
                  write(io_hist,*) 'variable non sauvegardee dans un fichier'
               endif 
               if (save_oa(iv_o).eq.1) then
                  write(io_hist,*) 'variable sauvegardee dans un fichier'
               endif 
               if (updv_oa(iv_o).eq.1) then
                  write(io_hist,*) 'variable utilisee pour une mise a jour',iv_s
               endif 

               write(io_hist,*) 
               write(io_hist,*) '----------------------------------'
               write(io_hist,*) '---caracteristiques de l''atome---'
               write(io_hist,*) '----------------------------------'

               write(io_hist,*) 'type de l''atome                : ',tpsi_oa(iv_o)
               
               if (tpsi_oa(iv_o).eq.0) then
                  write(io_hist,*) 'dirac'
               !endif
               else if (tpsi_oa(iv_o).eq.1) then
                  write(io_hist,*) 'ondelette morlet complexe fb_oa,fc_oa: ',fb_oa,fc_oa
               !endif
               else if (tpsi_oa(iv_o).eq.2) then
                  write(io_hist,*) 'windows fourier'
               !endif
               else if (tpsi_oa(iv_o).eq.3) then
                  write(io_hist,*) 'transformee de fourier classique'
               endif

               if (tpsi_oa(iv_o).ne.3) then
                  write(io_hist,*) 'largeur (h ou nbre de periodes): ',delta_t_oa(iv_o)
               endif
               write(io_hist,*) 'nombre de points par periode   : ',nzpt_per_oa(iv_o)

               write(io_hist,*) 
               write(io_hist,*) '------------------------------'
               write(io_hist,*) '---configuration temporelle---'
               write(io_hist,*) '------------------------------'
               if (unite_oa.eq.1.) then
                write(io_hist,*) 'unites                     : ','secondes'
               else
               if (unite_oa.eq.3600.) then
                 write(io_hist,*) 'unites                     : ','heures'
                else
                 write(io_hist,*) 'unites                     : ','autre'
                 endif
               endif
               write(io_hist,*) 'type echantillonnage       : ',swt_t_oa(iv_o)
               write(io_hist,*) 'premiere sortie (h/s)      : ',t0_oa(iv_o)
               write(io_hist,*) 'discretisation/integration : ',dori_oa(iv_o)

               if (swt_t_oa(iv_o).eq.1.or.swt_t_oa(iv_o).eq.4) then
                  write(io_hist,*) 'extraction: kount initial  : ',kount_user_oa(1,iv_o)
                  write(io_hist,*) '            kount final    : ',kount_user_oa(2,iv_o)
                  write(io_hist,*) '            delta(kount)   : ',kount_user_oa(3,iv_o)
               else if (swt_t_oa(iv_o).eq.3) then
                  write(io_hist,*) 'extraction: kount initial: ',kount_user_oa(1,iv_o)
               !else if (swt_t_oa(iv_o).eq.3) then #BLDX croco nt_max ?
                  write(io_hist,*) 'extraction: kount final  : ',nt_max-1
               endif 

               do lt_o  = begvt_oa(iv_o) , begvt_oa(iv_o+1) - 1
                  ! #BLXD Corrected Log. 
                  ! see call var_oa if (   kountv_oa(1,lt_o).ne.-9999.or.kountv_oa(2,lt_o).ne.-9999) then
                  if (   kountv_oa(1,lt_o).ne.-9999.and.kountv_oa(2,lt_o).ne.-9999) then
                     write(io_hist,*) 'BLXD localisation de l''atome n°',lt_o-begvt_oa(iv_o)+1,' (iic_oa) :',kountv_oa(1,lt_o), kountv_oa(2,lt_o)
                  else
                     write(io_hist,*) 'pas d''analyse pour la localisation n°',lt_o-begvt_oa(iv_o)+1
                  endif
               enddo

               write(io_hist,*) 
               write(io_hist,*) '----------------------------'
               write(io_hist,*) '---configuration spatiale---'
               write(io_hist,*) '----------------------------'
               write(io_hist,*) 'type d''echantillonnage        : ',swt_d_oa(iv_o)  

               if (swt_d_oa(iv_o).eq.1.or.swt_d_oa(iv_o).eq.3) then
#ifdef SPHERICAL
                  write(io_hist,*) ' BLXD CHECK lon, lat min max'
         if ( if_sphe_deg2rad ) then
                  write(io_hist,*) 'extraction: lat min,lat max   : ',lat_oa(1,iv_o)*180/pi_oa,lat_oa(2,iv_o)*180/pi_oa
                  write(io_hist,*) '            lon min,lon max   : ',lon_oa(1,iv_o)*180/pi_oa,lon_oa(2,iv_o)*180/pi_oa
         else
                  write(io_hist,*) 'extraction: lat min,lat max   : ',lat_oa(1,iv_o),lat_oa(2,iv_o)
                  write(io_hist,*) '            lon min,lon max   : ',lon_oa(1,iv_o),lon_oa(2,iv_o)
         endif
#else
                  write(io_hist,*) 'extraction: lat min,lat max   : ',lat_oa(1,iv_o),lat_oa(2,iv_o)
                  write(io_hist,*) '            lon min,lon max   : ',lon_oa(1,iv_o),lon_oa(2,iv_o)
#endif
               else if (swt_d_oa(iv_o).eq.2) then
                  write(io_hist,*) 'extraction: point (i,j)       : ',ptij_oa(1,iv_o),ptij_oa(2,iv_o)
               endif

               if (swt_d_oa(iv_o).eq.3)                               &
                  write(io_hist,*) '            prof min,prof max : ',h_oa(1,iv_o),h_oa(2,iv_o)
               write(io_hist,*)    '            kmin,kmax         : ',k_oa(1,iv_o),k_oa(2,iv_o)
               if (swt_d_oa(iv_o).ne.2) then
                  write(io_hist,*) '            di,dj,dk          : ',dx_oa(iv_o),dy_oa(iv_o),dk_oa(iv_o)
               else
                  write(io_hist,*) '            dk                : ',dk_oa(iv_o)
               endif

               write(io_hist,*) 
               write(io_hist,*) '---------------------------------'
               write(io_hist,*) '---configuration frequentielle---'
               write(io_hist,*) '---------------------------------'
               if (dori_oa(iv_o).eq.1) then
                  write(io_hist,*) 'periodes (t) discretes'
               else
                  write(io_hist,*) 'periodes (t) integrees'
               endif
               write(io_hist,*)    'tmin,dt,tmax (h ou s)      : ',per_oa (1,iv_o)/unite_oa,per_oa (3,iv_o)/unite_oa,per_oa (2,iv_o)/unite_oa
               do lp_o = begvp_oa(iv_o) , begvp_oa(iv_o+1)-1
                  write (io_hist,*) 'periode n°',lp_o-begvp_oa(iv_o)+1,'en h ou s     :',perv_oa(1,lp_o)/unite_oa
                  write (io_hist,*) 'coef. reconstruction  :',perv_oa(2,lp_o)
                  write( io_hist,*) 'resolution de l''atome :',resv_oa(lp_o)
               enddo

            enddo loop_variable
         enddo loop_config

         write(io_hist,*)
         write(io_hist,*)
         write(io_hist,*) '*************************************'
         write(io_hist,*) '          tailles en memoire'
         write(io_hist,*) '       des differentes structures'
         write(io_hist,*) '*************************************'
         write(io_hist,*)
! #BLXD
         write(io_hist,*) 'allocations simultanees (max) : ',nmsimult_oa_max
         write(io_hist,*) 'allocations simultanees (calc): ',nmsimult_oa
         write(io_hist,*) 'structure temporelle          : ',nzvt_oa
         !write(io_hist,*) 'structure spatiale (2d)       : ',nzvs_oa
         write(io_hist,*) 'BLXD CHECK tile structure spatiale (2d)       : ',pst%nzvs_oa
         !write(io_hist,*) 'structure spatiale (3d)       : ',nzvs3d_oa
         write(io_hist,*) 'BLXD CHECK tile structure spatiale (3d)       : ',pst%nzvs3d_oa
         write(io_hist,*) 'structure frequentielle       : ',nzvp_oa
         close(io_hist)
         !if (ifl_oa_out.eq.1) call sequential_end()
       !endif

         pst => null()

      endif chx5

      chx6 : if (ichoix.eq.6) then
!*******************************************************************************
! initial_oa   
!*******************************************************************************
!******************************************************************
!     alocations dynamiques:
!******************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
      
        pst => st(tile)

        open(unit=io_hist,file=trim(file_hist),position='append')   
        write (io_hist,*)
        write(io_hist,*) '*************************'
        write (io_hist,*) 'subroutine: allocate::'
        write(io_hist,*) '*************************'
        write (io_hist,*)
! #BLXD
        write (io_hist,*) 'nmsimult_oa_max =',nmsimult_oa_max
        write (io_hist,*) 'nmsimult_oa =',nmsimult_oa

! #BLD
        write (io_hist,*) 'nzw_oa      =',nzw_oa
        write (io_hist,*) 'nzv_oa      =',nzv_oa
        write (io_hist,*) 'BLXD nzvs_oa     =',pst%nzvs_oa
        write (io_hist,*) 'nzvt_oa     =',nzvt_oa
        write (io_hist,*) 'nzvp_oa     =',nzvp_oa
! #BLD typo 
!        write (io_hist,*) 'nzvc_oa     =',nzc_oa
         write (io_hist,*) 'nzvc_oa     =',nzvc_oa
        write (io_hist,*) 'BLXD nzvs3d_oa   =',pst%nzvs3d_oa
        write (io_hist,*)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
        pst => null()

      endif chx6

!STDALONE struct_oa eliminated      if (ichoix.eq.7) then
!*******************************************************************************
! struct_oa   
!*******************************************************************************
!      !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       !if (ifl_oa_out.eq.1) call sequential_begin()
!       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write(io_hist,*)
!       write(io_hist,*) '*************************'
!       write(io_hist,*) 'subroutine: struct_oa'
!       write(io_hist,*) '*************************'
!       close(io_hist)
!       !if (ifl_oa_out.eq.1) call sequential_end()
!      !endif
!     endif

!STDALONE struct_oa eliminated     if (ichoix.eq.8) then
!*******************************************************************************
! struct_oa   
!*******************************************************************************
!      !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       !if (ifl_oa_out.eq.1) call sequential_begin()
!       !open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write(3,*) '...fichier structure initialise.'
!       close(3)
!       !if (ifl_oa_out.eq.1) call sequential_end()
!      !endif
!     endif

!STDALONE subsave_init_oa eliminated     if (ichoix.eq.9) then
!*******************************************************************************
! subsave_init_oa   
!*******************************************************************************
!      !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       !if (ifl_oa_out.eq.1) call sequential_begin()
!       !open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write(io_hist,*)
!       write(io_hist,*) '*****************************'
!       write(io_hist,*) 'subroutine: save_init_oa'
!       write(io_hist,*) '*****************************'
!       close(io_hist)
!       !if (ifl_oa_out.eq.1) call sequential_end()
!      !endif
!     endif

      if (ichoix.eq.10) then
!*******************************************************************************
! subsave_oa   
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) '*** sauvegarde de la fenetre n°',i1_h,'***'
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
        !endif
      endif

      if (ichoix.eq.11) then
!*******************************************************************************
! upd_init_oa   
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*)
        write(io_hist,*) '*************************************'
        write(io_hist,*) ' allocation dynamique des variables 2d/3d'
        write(io_hist,*) 'nzupd2d_oa,nzupd3d_oa:',nzupd2d_oa,nzupd3d_oa
        write(io_hist,*) '*************************************'
        write(io_hist,*)
        close (io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif


      if (ichoix.eq.12) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
!STDALONE ichoix==12 NOT FOUND
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
       !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'initialisation des variables mises a jour'
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.13) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
!STDALONE ichoix==13 NOT FOUND
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'mise a jour des variables specifiques it=',i1_h,tallocated_oa(i1_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.14) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'mise a jour des variables specifiques it=',i1_h,tallocated_oa(i1_h)        
        write(io_hist,*) ' - tupd_oa    ',tupd_oa(i2_h)
        write(io_hist,*) ' - 2d-3d var  ',tgv3d_oa(tv_oa(i2_h))
        write(io_hist,*) ' - coeff type ',swt_wfpf_oa(i2_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.15) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'mise a jour des variables specifiques it=',i1_h,tallocated_oa(i1_h)        
        write(io_hist,*) ' - tupd_oa    ',tupd_oa(i2_h)
        write(io_hist,*) ' - 2d-3d var  ',tgv3d_oa(tv_oa(i2_h))
        write(io_hist,*) ' - coeff type ',swt_wfpf_oa(i2_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

!STDALONE energie_oaavelet removed      if (ichoix.eq.15) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
!      if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       if (ifl_oa_out.eq.1) call sequential_begin()
!       open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       write (io_hist,*) 'appel de la routine energie_oaavelet: ','ic=',i1_h
!       close (io_hist)
!       if (ifl_oa_out.eq.1) call sequential_end()
!      endif
!      endif

!STDALONE energie_oaavelet removed      if (ichoix.eq.16) then
!*******************************************************************************
! var_upd_oa   
!*******************************************************************************
!      if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       if (ifl_oa_out.eq.1) call sequential_begin()
!       open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       write (io_hist,*) 'appel de la routine energie_oaavelet: ' ,'ic=',i1_h
!       close (io_hist)
!       if (ifl_oa_out.eq.1) call sequential_end()
!      endif
!      endif

       return

       end subroutine history_oa

       function set_io_nodoa(iv, mynod, inteven, intodd)
        use module_oa_variables, only : nzv_oa
        implicit none
        integer, intent(in) :: mynod, inteven, intodd, iv
        integer :: ionodoa, set_io_nodoa
        if (nzv_oa >= 100) then
           ionodoa=1000
        else if (nzv_oa >= 10) then
           ionodoa=100
        else
           ionodoa=10
        endif 
        if ( MOD(iv,2)==0 ) then
         set_io_nodoa = (inteven*ionodoa+iv)*1000+mynod
        else
         set_io_nodoa = (intodd*ionodoa+iv)*1000+mynod
        end if
        return
       end function set_io_nodoa

      end module module_interface_oa
 
#else /* ONLINE_ANALYSIS */
      module module_interface_oa_empty
      end module
#endif /* ONLINE_ANALYSIS */

