! Copyright (C) 2022-2026 IFREMER
! License: CeCILL-C
! See LICENSES/LICENSE_MUSTANG.txt

#include "cppdefs.h"

MODULE comMUSTANG

#ifdef MUSTANG

    !!============================================================================
    !! ***  MODULE  comMUSTANG  ***
    !! Purpose : declare all common variables related to sediment dynamics
    !!============================================================================

    !! * Modules used
    USE comsubstance ! for lchain, rsh, rlg, riosh

    implicit none

    ! default
    public
  
    !! * Shared or public variables for MUSTANG 

    ! parameters
    REAL(kind=rlg), PARAMETER :: epsi30_MUSTANG = 1.e-30
    REAL(kind=rlg), PARAMETER :: epsilon_MUSTANG = 1.e-09 
    REAL(kind=rsh), PARAMETER :: valmanq = 999.0
    REAL(kind=riosh), PARAMETER :: rg_valmanq_io = 999.0

    ! namelists

    ! namsedim_init
    CHARACTER(len=19) :: date_start_dyninsed = '0000/01/01 00:00:00' ! starting
        ! date for dynamic processes in sediment; format '01/01/0000 00:00:00'
    CHARACTER(len=19) :: date_start_morpho = '0000/01/01 00:00:00' ! starting
        ! date for morphodynamic; format '01/01/0000 00:00:00'
    LOGICAL  :: l_initsed_vardiss = .false. !set to .true. if initialization
        ! of dissolved variables, temperature and salinity in sediment
        ! (will be done with concentrations in water at bottom (k=1))
    LOGICAL  :: l_unised = .true. !set to .true. for a uniform bottom initialization
    LOGICAL  :: l_unised_adjust_hsed = .false. ! set to .true. if we want to
        ! adjust the sediment thickness in order to be coherent with sediment
        ! parameters (calculation of a new hseduni based on
        ! cseduni, cvolmax values, and csed_ini of each sediment)
    CHARACTER(len=lchain) :: filrepsed = './' ! file path from which the model
        ! is initialized for the continuation of a previous run or non uniform init
    REAL(KIND=rsh) :: cseduni = 1500.0_rsh ! initial sediment concentration (kg/m3)
    REAL(KIND=rsh) :: hseduni = 1._rsh ! initial uniform sediment thickness (m)
    REAL(KIND=rsh) :: csed_mud_ini = 300.0_rsh ! real, mud concentration into initial
        ! sediment (kg/m3) (if = 0. ==> csed_mud_ini = cfreshmud)
    INTEGER        :: ksmiuni = 1 ! lower grid cell index in the sediment
    INTEGER        :: ksmauni = 1 ! upper grid cell index in the sediment
    REAL(KIND=rsh) :: sini_sed = 35.5_rsh ! initial interstitial water uniform salinity
    REAL(KIND=rsh) :: tini_sed = 10._rsh ! initial interstitial water uniform temperature
    REAL(KIND=rsh) :: poro_mud_ini = 0._rsh !only if key_MUSTANG_V2, initial porosity of
        ! mud fraction


    ! namsedim_layer
    LOGICAL  :: l_dzsmaxuni = .true. ! set to .true. dzsmax = dzsmaxuni ,
        ! if set to .false. then linearly computed in MUSTANG_sedinit
        ! from dzsmaxuni to dzsmaxuni/100 depending on water depth
    LOGICAL  :: l_dzsminuni = .false. !only if key_MUSTANG_V2, set to .false.
        ! if dzsmin vary with sediment bed composition, else dzsmin =  dzsminuni
    REAL(KIND=rsh) :: dzsminuni = 1.0e-3_rsh !only if key_MUSTANG_V2, minimum sediment
        ! layer thickness (m)
    REAL(KIND=rsh) :: dzsmin = 0.5e-2_rsh ! minimum sediment layer thickness (m)
    REAL(KIND=rsh) :: dzsmaxuni = 2.0e-2_rsh ! uniform maximum thickness for the superficial
        ! sediment layer (m), must be >0
    REAL(KIND=rsh) :: dzsmax_bottom = 2.0_rsh ! maximum thickness of bottom layers
        ! which result from the fusion when ksdmax is exceeded (m)
    REAL(KIND=rsh) :: k1HW97 = 0.07_rsh !only if key_MUSTANG_V2,
        ! ref value k1HW97 = 0.07, parameter to compute active layer
        ! thickness (Harris and Wiberg, 1997)
    REAL(KIND=rsh) :: k2HW97 = 6.0_rsh !only if key_MUSTANG_V2
        ! rref value k2HW97 = 6.0,parameter to compute active layer
        ! thickness (Harris and Wiberg, 1997)
    REAL(KIND=rsh) :: fusion_para_activlayer = 1._rsh !only if key_MUSTANG_V2
        ! criterion cohesiveness for fusion in active layer
        ! 0 : no fusion,
        ! = 1 : frmudcr1,
        ! > 1 : between frmudcr1 & frmudcr2
    INTEGER :: nlayer_surf_sed = 5 ! number of layers below the sediment surface
        ! that can not be melted (max thickness = dzsmax)
    LOGICAL :: l_splitlayersurf = .false. ! set to .true. to split surface sediment
        ! layers for a regular, precise discretization at the surface when
        ! too thick (over nlayer_surf_sed layers below the surface) ; if
        ! .false., the excess is simply moved into one new layer above


    ! namsedim_bottomstress
    LOGICAL  :: l_z0seduni = .true. ! boolean, set to .false. for z0sed computation from
        ! sediment diameter (if true, z0seduni is used)
    REAL(KIND=rsh) :: z0seduni = 0.00002_rsh ! uniform bed roughness (m)
    REAL(KIND=rsh) :: z0sedmud = 0.0001_rsh ! mud (i.e.minimum) bed roughness (m)
        ! (used only if l_unised is false)
    REAL(KIND=rsh) :: z0sedbedrock = 0.005_rsh ! bed roughness for bedrock (no sediment) (m)
        ! (used only if l_unised is false)
    LOGICAL :: l_tauskin_center = .false. ! boolean, to compute tauskin_c at rho point directly
    LOGICAL :: l_tauskin_ubar = .false. ! boolean, to use depth-averaged velocity instead of u(k=1)
    LOGICAL :: l_tauskin_upwind = .false. ! boolean, to use upwind interpolation for tauskin_x/y
    LOGICAL :: l_fricwave = .true. ! boolean, set to .true. if using wave related friction
        ! factor for bottom shear stress (from wave orbital velocity and period)
        ! if .false. then fricwav namelist value is used
    REAL(KIND=rsh) :: fricwav = 0.06_rsh ! default value is 0.06, wave related friction
        !factor (used for bottom shear stress computation)
    LOGICAL :: l_z0hydro_coupl_init = .false. ! boolean, set to .true. if evaluation of
        ! z0 hydro depends on sediment composition at the beginning
        ! of the simulation
    LOGICAL :: l_z0hydro_coupl = .false. ! boolean, set to .true. if evaluation of
        ! z0 hydro depends on sediment composition along the run
    REAL(KIND=rsh) :: coef_z0_coupl = 1._rsh ! parameter to compute z0hydro in the
        ! first centimeter : z0hydro = coef_z0_coupl * sand diameter
    REAL(KIND=rsh) :: z0_hydro_mud = 0.0001_rsh ! z0hydro if pure mud (m)
    REAL(KIND=rsh) :: z0_hydro_bed = 0.005_rsh ! z0hydro if no sediment (m)


    ! namsedim_deposition
    REAL(KIND=rsh) :: cfreshmud = 550.0_rsh ! fresh deposit concentration (kg/m3)
        ! (must be around 100 if consolidation
        ! or higher (300-500 if no consolidation)
    REAL(KIND=rsh) :: csedmin = 100._rsh ! concentration of the upper layer under
        ! which there is fusion with the underlying sediment cell (kg/m3)
    REAL(KIND=rsh) :: cmudcr = 2000._rsh ! critical relative concentration of the surface
        ! layer above which no mixing is allowed with the underlying
        ! sediment (kg/m3)
    REAL(KIND=rsh) :: aref_sand = 0.02_rsh ! parameter used in sandconcextrap,
        ! reference height above sediment, used for computing of
        ! sand deposit. Parameter used for sand extrapolation on water
        ! column and correct sand transport, value by default = 0.02
        ! correspond to Van Rijn experiments
        ! DO NOT CHANGED IF NOT EXPERT
    LOGICAL :: l_corflux = .true. ! boolean to activate correction on horizontal sand
        ! fluxes
    REAL(KIND=rsh) :: cvolmaxsort = 0.58_rsh ! max volumic concentration of sorted sand
    REAL(KIND=rsh) :: cvolmaxmel = 0.67_rsh ! maxvolumic concentration of mixed sediments
    LOGICAL :: l_slipdeposit = .false. ! boolean to activate sliding (avalanching) of
        ! deposited sediment on steep slopes
    REAL(KIND=rsh) :: slopefac = 0.01_rsh ! slope effect multiplicative on deposit
        ! (used only if l_slipdeposit)


    ! namsedim_erosion
    REAL(KIND=rsh) :: activlayer = 0.02_rsh ! active layer thickness (m)
    REAL(KIND=rsh) :: frmudcr2 = 0.7_rsh ! critical mud fraction under which the
        ! behaviour is intermediate between sand and mud and above which the
        ! behavior is purely muddy
    REAL(KIND=rsh) :: coef_frmudcr1 = 1000._rsh ! to compute critical mud fract. frmudcr1
        ! underwhich the behaviour is purely sandy
        ! (frmudcr1=min(coef_frmudcr1*d50 sand,frmudcr2))
    REAL(KIND=rsh) :: x1toce_mud = 0.1_rsh ! coef. for the formulation of the critical
        ! erosion stress in mud behavior toce=x1toce*csed**x2toce
    REAL(KIND=rsh) :: x2toce_mud = 0._rsh ! coef. for the formulation of the critical
    ! erosion stress in mud behavior toce=x1toce*csed**x2toce
    REAL(KIND=rsh) :: E0_sand_para = 1._rsh ! coefficient used to modulate erosion
        ! flux for sand (=1 if no correction )
    REAL(KIND=rsh) :: n_eros_sand = 1.5_rsh ! parameter for erosion flux for sand
        ! (E0_sand*(tenfo/toce-1.)**n_eros_sand )
        ! WARNING : choose parameters compatible with E0_sand_option
        ! (example : n_eros_sand=1.6 for E0_sand_option=1)
    REAL(KIND=rsh) :: E0_mud = 0.00001_rsh ! erosion flux for mud
    REAL(KIND=rsh) :: n_eros_mud = 1._rsh ! E0_mud*(tenfo/toce-1.)**n_eros_mud
    INTEGER        :: ero_option = 3 ! choice of erosion formulation for mixing
        ! sand-mud
        ! ero_option= 0 : pure mud behavior
        ! ero_option= 1 : linear interpolation between sand and mud behavior,
        !   depend on proportions of the mixture
        ! ero_option= 2 : formulation derived from that of J.Vareilles (2013)
        ! ero_option= 3 : formulations proposed by B. Mengual (2015) with
        !   exponential coefficients depend on proportions of the mixture
    INTEGER        :: E0_sand_option = 0 ! integer, choice of formulation for
        ! E0_sand evaluation :
        ! E0_sand_option = 0 E0_sand = E0_sand_Cst
        ! E0_sand_option = 1 E0_sand evaluated with Van Rijn (1984)
        ! E0_sand_option = 2 E0_sand evaluated with erodimetry
        !    (min(0.27,1000*d50-0.01)*toce**n_eros_sand)
        ! E0_sand_option = 3 E0_sand evaluated with Wu and Lin (2014)
    REAL(KIND=rsh) :: xexp_ero = 40.0_rsh !used only if ero_option=3 : adjustment on
        ! exponential variation  (more brutal when xexp_ero high)
    REAL(KIND=rsh) :: E0_sand_Cst = 0.00594_rsh ! constant erosion flux for sand
        ! (used if E0_sand_option= 0)
    REAL(KIND=rsh) :: E0_mud_para_indep = 1._rsh !only if key_MUSTANG_V2,
        ! parameter to correct E0_mud in case of erosion
        ! class by class in non cohesive regime
    LOGICAL        :: l_peph_suspension = .false. !only if key_MUSTANG_V2,
        ! set to .true. if hindering / exposure processes in critical
        ! shear stress estimate for suspension
    LOGICAL        :: l_eroindep_noncoh = .true. !only if key_MUSTANG_V2,
        ! set to .true. in order to activate independant erosion for
        ! the different sediment classes sands and muds
        ! set to .false. to have the mixture mud/sand eroded as in V1
    LOGICAL        :: l_eroindep_mud = .false. !only if key_MUSTANG_V2,
        ! set to .true. if mud erosion independant for sands erosion
        ! set to .false. if mud erosion proportionnal to total sand erosion
    LOGICAL        :: l_xexp_ero_cst = .false. !only if key_MUSTANG_V2, set to .true.
        ! if xexp_ero estimated from empirical formulation, depending on
        ! frmudcr1
    INTEGER        :: tau_cri_option = 0 !only if key_MUSTANG_V2,
        ! choice of critical stress formulation ,
        ! 0: Shields 1: Wu and Lin (2014)
    INTEGER        :: tau_cri_mud_option_eroindep = 1 !only if key_MUSTANG_V2
        ! choice of mud critical stress formulation
        ! 0: x1toce_mud*cmudr**x2toce_mud
        ! 1: toce_meansan if somsan>eps (else->case0)
        ! 2: minval(toce_sand*cvsed/cvsed+eps) if >0 (else->case0)
        ! 3: min( case 0; toce(isand2) )


    ! namsedim_poro
    INTEGER :: poro_option = 2 ! choice of porosity formulation
        ! 1: Wu and Li (2017) (incompatible with consolidation))
        ! 2: mix ideal coarse/fine packing
    REAL(KIND=rsh) :: Awooster = 0.42_rsh ! parameter of the formulation of
        ! Wooster et al. (2008) for estimating porosity associated to the
        ! non-cohesive sediment see Cui et al. (1996) ref value = 0.42
    REAL(KIND=rsh) :: Bwooster = -0.458_rsh ! parameter of the formulation of
    ! Wooster et al. (2008) for estimating porosity associated to the
    ! non-cohesive sediment see Cui et al. (1996) ref value = -0,458
    REAL(KIND=rsh) :: Bmax_wu = 0.65_rsh ! maximum portion of the coarse sediment class
        ! participating in filling , ref value = 0.65
    REAL(KIND=rsh) :: poro_min = 0.2_rsh ! minimum porosity below which consolidation
        ! is stopped


#if defined key_MUSTANG_V2
    ! namsedim_bedload
    LOGICAL :: l_peph_bedload = .false. ! set to .true. if hindering / exposure processes
        ! in critical shear stress estimate for bedload
    LOGICAL :: l_slope_effect_bedload = .true. ! set to .true. if accounting for slope
        ! effects in bedload fluxes (Lesser formulation)
    LOGICAL :: l_fsusp = .false. ! limitation erosion fluxes of non-coh sediment in case
        ! of simultaneous bedload transport, according to Wu & Lin formulations
        ! set to .true. if erosion flux is fitted to total transport
        ! should be set to .false. if E0_sand_option=3 (Wu & Lin)
    REAL(KIND=rsh) :: alphabs = 1.0_rsh ! coefficient for slope effects (default
        ! coefficients Lesser et al. (2004), alphabs = 1.)
    REAL(KIND=rsh) :: alphabn = 1.5_rsh ! coefficient for slope effects (default
        ! coefficients Lesser et al. (2004), default alphabn is 1.5 but
        ! can be higher, until 5-10 (Gerald Herling experience))
    REAL(KIND=rsh) :: hmin_bedload = 0.1_rsh ! no bedload in u/v directions if
        ! h0+ssh <= hmin_bedload in neighbouring cells
#endif

    ! namsedim_lateral_erosion

    LOGICAL        :: l_erolat = .false. ! set to .true to activate lateral erosion
    REAL(KIND=rsh) :: htncrit_eros = 0.0_rsh ! critical water height so as to prevent
        ! erosion under a given threshold (the threshold value is different for
        ! flooding or ebbing, cf. Hibma's PhD, 2004, page 78)
    REAL(KIND=rsh) :: coef_erolat = 0.002_rsh ! slope effect multiplicative factor
    REAL(KIND=rsh) :: coef_tauskin_lat = 5.0_rsh ! parameter to evaluate the lateral
        ! stress as a function of the average tangential velocity on the
        ! vertical
    LOGICAL        :: l_erolat_wet_cell = .false. ! set to .true in order to take
        ! into account wet cells lateral erosion

    ! namsedim_consolidation
    LOGICAL        :: l_consolid = .false. ! set to .true. if sediment consolidation is
        ! accounted for
    REAL(KIND=rsh) :: dt_consolid = 600.0_rsh ! time step for consolidation processes
    REAL(KIND=rlg) :: subdt_consol = 30.0_rlg ! sub time step for consolidation and
                                   ! particulate bioturbation  in sediment
    REAL(KIND=rsh) :: csegreg = 250.0_rsh ! NOT CHANGE VALUE if not expert, default 250.0
    REAL(KIND=rsh) :: csandseg = 1250.0_rsh ! NOT CHANGE VALUE if not expert, default 1250.0
    REAL(KIND=rsh) :: xperm1 = 4.0e-12_rsh ! permeability=xperm1*d50*d50*voidratio**xperm2
    REAL(KIND=rsh) :: xperm2 = -6.0_rsh ! permeability=xperm1*d50*d50*voidratio**xperm2
    REAL(KIND=rsh) :: xsigma1 = 6.0e+05_rsh ! parameter used in Merckelback & Kranenburg s
        ! (2004) formulation NOT CHANGE VALUE if not expert, default 6.0e+05
    REAL(KIND=rsh) :: xsigma2 = 6._rsh ! parameter used in Merckelback & Kranenburg s
        ! (2004) formulation NOT CHANGE VALUE if not expert, default 6


    ! namsedim_diffusion
    LOGICAL        :: l_diffused = .false. ! set to .true. if taking into account
        ! dissolved diffusion in sediment and at the water/sediment interface
    REAL(KIND=rsh) :: dt_diffused = 500.0_rsh ! time step for diffusion in sediment
    INTEGER        :: choice_flxdiss_diffsed = 3 ! choice for expression of
        ! dissolved fluxes at sediment-water interface
        ! 1 : Fick law : gradient between Cv_wat at dz(1)/2
        ! 2 : Fick law : gradient between Cv_wat at distance epdifi
    REAL(KIND=rsh) :: xdifs1 = 1.e-8_rsh ! diffusion coefficients within the sediment
    REAL(KIND=rsh) :: xdifsi1 = 1.e-6_rsh ! diffusion coefficients at the water-sediment
        ! interface
    REAL(KIND=rsh) :: epdifi = 0.01_rsh ! diffusion thickness in the water at the
        ! sediment-water interface
    REAL(KIND=rsh) :: fexcs = 0.5_rsh ! factor of eccentricity of concentrations in
        ! vertical fluxes evaluation (.5 a 1)


    ! namsedim_bioturb
    LOGICAL        :: l_bioturb = .false. ! set to .true. if taking into account
        ! particulate bioturbation (diffusive mixing) in sediment
    LOGICAL        :: l_biodiffs = .false. ! set to .true. if taking into account
        ! dissolved bioturbation diffusion in sediment
    REAL(KIND=rsh) :: dt_bioturb = 600.0_rsh ! time step for bioturbation in sediment
    REAL(KIND=rsh) :: subdt_bioturb = 30.0_rsh ! sub time step for bioturbation
    REAL(KIND=rsh) :: xbioturbmax_part = 8.2e-12_rsh ! max particular bioturbation
        ! coefficient by bioturbation Db (in surface)
    REAL(KIND=rsh) :: xbioturbk_part = 6.0_rsh ! for part. bioturbation coefficient
        ! between max Db at sediment surface and 0 at bottom
    REAL(KIND=rsh) :: dbiotu0_part = 0.1_rsh ! max depth beneath the sediment
        ! surface below which there is no bioturbation
    REAL(KIND=rsh) :: dbiotum_part = 0.1_rsh ! sediment thickness where the
        ! part-bioturbation coefficient Db is constant (max)
    REAL(KIND=rsh) :: xbioturbmax_diss = 1.157e-09_rsh ! max diffusion coeffient by
        ! biodiffusion Db (in surface)
    REAL(KIND=rsh) :: xbioturbk_diss = 6.0_rsh ! coef (slope) for biodiffusion
        ! coefficient between max Db at sediment surface and 0 at bottom
    REAL(KIND=rsh) :: dbiotu0_diss = 0.1_rsh ! max depth beneath the sediment
        ! surface below which there is no bioturbation
    REAL(KIND=rsh) :: dbiotum_diss = 0.005_rsh ! sediment thickness where the
        ! diffsolved-bioturbation coefficient Db is constant (max)
    REAL(KIND=rsh) :: frmud_db_min = 0.6_rsh ! mud fraction limit (min) below which
        ! there is no Biodiffusion
    REAL(KIND=rsh) :: frmud_db_max = 0.8_rsh ! mud fraction limit (max)above which
        ! the biodiffusion coefficient Db is maximum (muddy sediment)


    ! namsedim_morpho
    LOGICAL :: l_morphocoupl = .false. ! set to .true if coupling module morphodynamic
    LOGICAL :: l_MF_dhsed = .false. ! set to .true. if morphodynamic applied with
        ! sediment height variation amplification
        ! (MF_dhsed = MF; then MF will be = 0)
        ! set to .false. if morphodynamic is applied with
        ! erosion/deposit fluxes amplification (MF_dhsed not used)
    REAL(KIND=rsh) :: MF = 1.0_rsh ! morphological factor : multiplication factor for
        ! morphologicalevolutions, equivalent to a "time acceleration"
        ! (morphological evolutions over a MF*T duration are assumed to be
        ! equal to MF * the morphological evolutions over T).
    REAL(KIND=rlg) :: dt_morpho = 0.1_rlg ! time step for morphodynamic (s)


#if  ! defined key_noTSdiss_insed
    ! namtempsed
    REAL(KIND=rsh) :: mu_tempsed1 = 8.e-7_rsh ! parameters used to estimate thermic
        ! diffusitiyfunction of mud fraction
    REAL(KIND=rsh) :: mu_tempsed2 = -1.4e-6_rsh ! parameters used to estimate thermic
        ! diffusitiyfunction of mud fraction
    REAL(KIND=rsh) :: mu_tempsed3 = 9.e-7_rsh ! parameters used to estimate thermic
        ! diffusitiyfunction of mud fraction
    REAL(KIND=rsh) :: epsedmin_tempsed = 0.2_rsh ! sediment thickness limits for
         ! estimation heat loss at bottom, if hsed < epsedmin_tempsed :
        ! heat loss at sediment bottom = heat flux a sediment surface
    REAL(KIND=rsh) :: epsedmax_tempsed = 2._rsh ! sediment thickness limits for
        ! estimation heat loss at bottom, if hsed > epsedmax_tempsed :
        ! heat loss at sediment bottom = 0.
#endif


    ! namsedoutput
    LOGICAL :: l_outsed_nb_lay_sed = .true. ! To output the current number of layer
    LOGICAL :: l_outsed_hsed = .true. ! To the sediment thickness
    LOGICAL :: l_outsed_tauskin = .true. ! To output the total skin stress
    LOGICAL :: l_outsed_tauskin_cw = .false. ! To output the current skin stress and the wave skin stress
    LOGICAL :: l_outsed_poro = .false. ! To output the porosity
    LOGICAL :: l_outsed_dzs = .true. ! To output the sediment thickness of each layer
    LOGICAL :: l_outsed_temp_sed = .true. ! To output the temperature in sediment
    LOGICAL :: l_outsed_salt_sed = .true. ! To output the salinity in sediment
    LOGICAL :: l_outsed_cv_sed = .true. ! To output each sediment class concentration
    LOGICAL :: l_outsed_ws = .false. ! To output each MUD class settling velocities
    LOGICAL :: l_outsed_toce = .false. ! To output the critical erosion stress
    LOGICAL :: l_outsed_flx_s2w_w2s = .false. ! To output the sediment to water and the water to sediment fluxes
    LOGICAL :: l_outsed_pephm_fcor = .false. ! To output the hindering exposure factor
    LOGICAL :: l_outsed_bedload = .false. ! To output the bedload flux along x/y-axis  andthe divergence of bedload flux
    LOGICAL :: l_outsed_fsusp = .false. ! To output the fraction of transport in susp
    LOGICAL :: l_outsed_frmudsup = .false. ! To output mud fraction in the ksmax layer
    LOGICAL :: l_outsed_dzs_ksmax = .false. ! To output layer thickness at sed. surface
    LOGICAL :: l_outsed_theoric_active_layer = .false. ! To output theoric act. lay.
    LOGICAL :: l_outsed_ero_details = .false. ! To output iterations in sed_erosion and part in coh and noncoh during time step
    LOGICAL :: l_outsed_z0sed = .false. ! Skin roughness length
    LOGICAL :: l_outsed_z0hydro = .false. ! Hydrodynamic roughness length
    LOGICAL :: l_outsed_consolidation = .false. ! To output consolidation variables
    INTEGER :: nk_nivsed_out = 5 ! number of saved sediment layers
        ! =ksdmax if choice_nivsed_out = 1
        ! <=ksdmax if choice_nivsed_out = 2,
        ! unused if choice_nivsed_out = 3
        !  <6 if choice_nivsed_out = 4,
    INTEGER :: choice_nivsed_out = 1 ! choice of saving output  (1 to 4)
    REAL(KIND=rsh), DIMENSION(5) :: ep_nivsed_out = (/0._rsh, 0._rsh, 0._rsh, 0._rsh, 0._rsh/) ! 5 values of sediment
        ! layer thickness (mm), beginning with surface layer
        ! (used if choice_nivsed_out=4)
    REAL(KIND=rsh) :: epmax_nivsed_out = 0._rsh ! maximum thickness (mm) for
        ! output each layers of sediment (used if choice_nivsed_out=3).
        ! Below the layer which bottom level exceed this thickness,
        ! an addition layer is an integrative layer till bottom


    ! namflocmod
    LOGICAL :: l_flocmod = .false. ! set to .true. to activate the FLOCMOD flocculation
        ! module for mud settling velocity
    LOGICAL :: l_ASH = .true. ! set to .true. if aggregation by shear
    LOGICAL :: l_ADS = .false. ! set to .true. if aggregation by differential settling
    LOGICAL :: l_COLLFRAG = .false. ! set to .true. if fragmentation by collision
    INTEGER :: f_ero_iv = 1 ! fragment class (mud variable index corresponding to
        ! the eroded particle size - typically 1)
    REAL(KIND=rsh) :: f_ater = 0.0_rsh ! ternary fragmentation factor : proportion of
        ! flocs fragmented as half the size of the initial binary fragments
        ! (0.0 if full binary fragmentation, 0.5 if ternary fragmentation)
    REAL(KIND=rsh) :: f_dmin_frag = 0.00001_rsh ! minimum diameter for fragmentation
        ! (default 10e-6 microns)
    REAL(KIND=rsh) :: f_ero_frac = 0.0_rsh ! floc erosion (% of floc mass eroded)
        ! (default 0.05)
    REAL(KIND=rsh) :: f_ero_nbfrag = 2.0_rsh ! nb of fragments produced by erosion
        ! (default 2.0)
    REAL(KIND=rsh) :: f_mneg_param = 0.001_rsh ! negative mass after
        ! flocculation/fragmentation allowed before redistribution
        ! (default 0.001 g/l)
    REAL(KIND=rsh) :: f_collfragparam = 0.01_rsh ! fraction of shear aggregation leading
        ! to fragmentation by collision (default 0.0, must be less than 1.0)
    REAL(KIND=rsh) :: f_cfcst = 0.1875_rsh ! fraction of mass lost by flocs if fragmentation
        ! by collision .. (default : =3._rsh/16._rsh)
    REAL(KIND=rsh) :: f_fp = 0.1_rsh ! relative depth of inter particle penetration
        ! (default =0.1) (McAnally, 1999)
    REAL(KIND=rsh) :: f_fy = 1.0e-10_rsh ! floc yield strength  (default= 1.0e-10)
        ! (Winterwerp, 2002)
    REAL(KIND=rsh) :: f_dp0 = 4.e-6_rsh ! primary particle size (default 4.e-6 m)
    REAL(KIND=rsh) :: f_alpha = 0.15_rsh ! flocculation efficiency parameter
        ! (default 0.15)
    REAL(KIND=rsh) :: f_beta = 0.150_rsh ! floc break up parameter (default 0.1)
    REAL(KIND=rsh) :: f_nb_frag = 2.0_rsh ! nb of fragments of equal size by shear
        ! fragmentation (default 2.0 as binary fragmentation)
    REAL(KIND=rsh) :: f_nf = 2.0_rsh ! fractal dimension (default 2.0, usual range from
        ! 1.6 to 2.8)
    REAL(KIND=rsh) :: f_clim = 0.001_rsh ! min concentration below which flocculation
        !processes are not calculated

    CHARACTER(len=lchain) :: dredging_location_file = '' ! TODO DREDGING
    CHARACTER(len=lchain) :: dredging_settings_file = '' ! TODO DREDGING
    CHARACTER(len=lchain) :: dredging_out_file = '' ! TODO DREDGING
    INTEGER :: dredging_dumping_layer = 1 ! TODO DREDGING
    REAL(KIND=rsh) :: dredging_dt = 3600._rsh ! TODO DREDGING
    REAL(KIND=rsh) :: dredging_dt_out = 3600._rsh ! TODO DREDGING

! end namelist variables


    REAL(KIND=rsh) :: h0fond  ! residual thickness in water  (m)

    ! fwet =1 if not used  
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: fwet

    ! Initialization 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: cini_sed
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: cv_sedini
    REAL(KIND=rsh) :: hsed_new

    ! Fluxes at the interface water-sediment
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_sum
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_CROCO
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_CROCO
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_sum_CROCO

    ! Sediment parameters
    REAL(KIND=rsh)            :: ros_sand_homogen 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: typart
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: diamstar
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: ws_sand
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: rosmrowsros
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: stresscri0
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: tetacri0


    INTEGER :: nv_use
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ksmi
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ksma
    REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE   :: cv_sed
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: c_sedtot
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: poro
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: dzs       
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: dzsmax
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: gradvit       

    ! Sediment height
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed_previous

    ! Bottom stress variables
    REAL(KIND=rsh) :: fws2  ! fricwav/2   
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE     :: z0sed ! roughness (m)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin ! max bottom stress due to the combinaison current/wave (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_c ! bottom stress due to current (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_w ! bottom stress due to wave (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_x ! bottom stress - component on x axis
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_y ! bottom stress - component on y axis
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: ustarbot ! (m/s)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: raphbx ! adim.
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: raphby ! adim.

    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w
    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w_consol
    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w_drycell

    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: htot
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: alt_cw1

    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: sal_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: temp_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: epn_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: cw_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: ws3_bottom_MUSTANG ! settling velocities in  bottom cell (m/s)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: roswat_bot

    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE   :: corflux
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE   :: corfluy

    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: rouse2D ! Rouse2D number
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: rouse2D_integral ! see integrate_rouse_profile


    ! Dynamic in sediment (consolidation/diffusion/bioturbation)
    REAL(KIND=rlg)   :: tstart_dyninsed ! time beginning dynamic in sediment
    REAL(KIND=rlg)   :: t_dyninsed      ! time of next dynamic in sediment step
    REAL(KIND=rlg)   :: dt_dyninsed     ! time step for dynamic in sediment (min of dt for each process)
    LOGICAL :: l_dyn_insed ! true if (l_consolid .OR. l_bioturb .OR. l_diffused .OR. l_biodiffs)
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fludif
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fluconsol
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fluconsol_drycell
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flu_dyninsed
   
    ! Diffusion
    REAL(KIND=rsh) :: cexcs

    ! Morphodynamic
    REAL(KIND=rlg) :: tstart_morpho   ! time beginning morphodynamic
    REAL(KIND=rlg) :: t_morpho        ! time of next morphodynamic step
    REAL(KIND=rsh) :: MF_dhsed

#ifdef key_MUSTANG_V2
    REAL(KIND=rsh) :: coeff_dzsmin
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: l_isitcohesive
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: psi_sed
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: poro_mud
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: crel_mud
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: sigmapsg
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: stateconsol
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: permeab
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: E0_sand
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE  :: flx_bx
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE  :: flx_by
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: slope_dhdx
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: slope_dhdy
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: sedimask_h0plusxe
#if defined MORPHODYN
            INTEGER :: it_morphoYes
#endif
#endif

    ! Sedim output
    INTEGER                            :: rstMust_nbvar
    INTEGER, DIMENSION(:), ALLOCATABLE :: rstMust                    ! Output identifier
    INTEGER, DIMENSION(:), ALLOCATABLE :: hisMust                    ! Output identifier
    INTEGER, DIMENSION(:), ALLOCATABLE :: avgMust                    ! Output identifier
    LOGICAL, DIMENSION(:), ALLOCATABLE :: rstoutintegerMust          ! To indicate if the output variable is integer
    LOGICAL, DIMENSION(:), ALLOCATABLE :: rstout2DMust               ! To indicate if the output variable is 2D or 3D
    LOGICAL, DIMENSION(:), ALLOCATABLE :: rstout3DsedMust            ! To indicate if the output variable is 3D with sediment vertical axis
    INTEGER                            :: outMust_nbvar              ! Number of variables available in output
    LOGICAL, DIMENSION(:), ALLOCATABLE :: outMust                    ! To choose which variable is outputed
    LOGICAL, DIMENSION(:), ALLOCATABLE :: out2DMust                  ! To indicate if the output variable is 2D or 3D
    LOGICAL, DIMENSION(:), ALLOCATABLE :: out3DsedMust            ! To indicate if the output variable is 3D with sediment vertical axis
    CHARACTER(LEN=75), DIMENSION(:, :), ALLOCATABLE :: vname_Must    ! Vector of characteristics of each outputed variables
    CHARACTER(LEN=75), DIMENSION(:, :), ALLOCATABLE :: vname_rstMust    ! Vector of characteristics of each outputed variables

    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: ep_nivsed_outp1
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE :: var2D_hsed
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_cvsed
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_poro
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_dzs
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_TEMP
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_SAL
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_flx_s2w
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_flx_w2s
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_pephm_fcor
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_flx_bx
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_flx_by
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_bil_bedload
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_fsusp
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_toce
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_frmudsup
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_dzs_ksmax
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_theoric_active_layer
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_tero_noncoh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_tero_coh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_pct_iter_noncoh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_pct_iter_coh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_niter_ero
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_flx_s2w_coh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_flx_w2s_coh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_flx_s2w_noncoh 
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_flx_w2s_noncoh
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_flx_bx_int
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_flx_by_int
    REAL(KIND=riosh), DIMENSION(:,:), ALLOCATABLE  :: var2D_bil_bedload_int

    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_loadograv
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_permeab
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_sigmapsg
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_dtsdzs
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_hinder
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_sed_rate
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_sigmadjge
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3Dksed_stateconsol

#if defined key_BLOOM_insed
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_diagsed
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_diagsed
#endif

!  used in lateral_erosion only 
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corim1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corip1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corjm1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corjp1
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corim1
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corip1
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corjm1
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corjp1

   ! slipdeposit fluxes (used only if l_slipdeposit)
   !  used in accretion (settling) only bud exchange and dimensions could depend on grid model
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corin
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corjp1


#if ! defined key_noTSdiss_insed
    ! Temperature in sediment 
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phitemp_s
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ivdiss
    INTEGER       , DIMENSION(:), ALLOCATABLE :: D0_funcT_opt
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: D0_m0
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: D0_m1
#endif

#if ! defined key_nofluxwat_IWS && ! defined key_noTSdiss_insed
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: WATER_FLUX_INPUTS ! not operationnal, stil to code **TODO**
#endif

#ifdef key_BLOOM_insed
    LOGICAL :: l_out_subs_diag_sed
#endif


    CONTAINS
 
#endif /* ifdef MUSTANG */

END MODULE comMUSTANG
