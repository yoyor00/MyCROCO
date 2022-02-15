
#if defined BIOLink

   !&E==========================================================================
   !&E                   ***  coupleur_define_BIOLink  ***
   !&E
   !&E                 *********** modele CROCO *******************
   !&E
   !&E ** Purpose : definitions of dimensions, variables and parameters 
   !&E               for module BIOLink         
   !&E 
   !&E ** Description : must be completed by the user
   !&E          when coupling with a hydrodynamic model and a module BIO
   !&E
   !&E ** History :
   !&E     ! 2018-11  (B.Thouvenin )    : creation for portability and 
   !&E                                    coupling with any hydro model 
   !&E
   !&E     ! 2022-02 (G. Koenig) : Ordering of variables and updates of comments
   !&E==========================================================================

/*************************************************************************/
/*************************************************************************/
/*********************General variables and size**************************/
/*********************of the computationnal domain************************/
/*************************************************************************/
/*************************************************************************/


/*!====================================================================
  ! Definition of general and math variables
  !====================================================================*/


# define NUMBER_PI pi /* The number pi */

/*!====================================================================
  ! Variables related to the namelists and other external files
  !====================================================================*/
# define REPFICNAMELIST 'FIC_NAMELIST' /* Directory containing the namelists */
# define SUFFIX_RUNFILE 1 /* I do not really know */
# define TEST_NOT_RESTART_SUBS .NOT.l_restart_subs /* Boolean to check 
                                                      if a restart was used */
# define TEST_NOT_INITFROMFILE .NOT.l_initfromfile /* Boolean to check 
                                                      if a init file was used */

/*!====================================================================
  ! MPI and system commands
  !====================================================================*/

# if ! defined key_MARS
!/* MPI */
#   define CALL_MPI !call_MPI
#   define IF_MPI !if
#   define USE_MPI !use
#   define ENDIF_MPI !endif
#   define PRINT_DBG !print
#   define IF_AGRIF !if
#   define OMPMPI  !OMPMPI
!/*CROCO */
#   define iscreenlog stdout
#   define ierrorlog stdout
#   define iwarnlog stdout
# endif

/*!====================================================================
  ! Definition of the hydrodynamical model used
  !====================================================================*/

# define BIOLink_HOST BIOLink_CROCO /* The physical model */

/*=====================================================================
! Spatial Grid limits definition  of loops
! inside the domain - except meshes at open boundaries
!======================================================================*/

# define IMIN_BOUCL 1
# define IMAX_BOUCL Lm /* Zonal dimension*/
# define JMIN_BOUCL 1
# define JMAX_BOUCL Mm /* Meridional dimension */

/*=====================================================================
!  Dimension of the different kind of tables
!======================================================================*/

# define COMPLETE_ARRAY      1   :Lm    , 1   :Mm /* Size of a 2D array */
# define PROC_IN_ARRAY       GLOBAL_2D_ARRAY   
# define PROC_IN_ARRAY_m1p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m1p1  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m2p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_0p1   GLOBAL_2D_ARRAY /* Size of an array divided by MPI. 
                                                Here BIOLink uses different
                                                Size of arrays while the physical model
                                                Only has GLOBAL_2D_ARRAY. Thus different
                                                Names are attributed to GLOBAL_2D_ARRAY */

/*=====================================================================
! Geometry of the computational grid
!======================================================================*/

# define NB_LAYER_WAT N /* Number of vertical levels */
# define CELL_SURF surf_cell /* Surface a of computationnal cell */

/*=====================================================================
!  Position variables
!======================================================================*/

# define LONGITUDE lonr /* Longitude */
# define LATITUDE latr  /* Latitude */

/*************************************************************************/
/*************************************************************************/
/****************BIOLink tracer and time variables************************/
/*************************************************************************/
/*************************************************************************/

/*=====================================================================
!  Time-related variables: Timesteps and dates
!======================================================================*/

# define TRANSPORT_TIME_STEP dt /* Hydro model time step */
# define ECO_TIME_STEP dt_bio_update /* Time step of the biological model */
# define CURRENT_TIME time /* Current time/date for the biological model */
# define TIME_BEGIN tdeb  /* Starting time/date for the biological model */

/*!====================================================================
!  BIOLink time variables
!======================================================================*/

# if defined ECO3M
#   define BIO_TIME_STEP dt_bio /* Time step for the biological model */
#   define DT_CONSERV_BIOLINK dt_bio_conserv /* Time step for BIOLink conservativity routine */
# endif /* ECO3M */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define BIO_TIME_STEP dtbio /* Time step for the biological model */
#   define DT_CONSERV_BIOLINK dt_conserv_BIO /* Time step for BIOLink conservatity routine */
# endif /* PEPTIC/BLOOM/METeOR */


/*=====================================================================
! Variables related to the names and characteristics of tracers
!======================================================================*/

# if ! defined key_MARS
#   define NAME_SUBS vname(1,indxT+ntrc_salt+isubs) /* Array of name of substance variables */
# endif /* key_MARS */

/*=====================================================================
! Number of tracer variables and table of the hydro model where 
! they are stored
!======================================================================*/

# define NBVARADV_TOT  ntrc_subs /* Total number of tracer variables added by BIOLink */
# define WATER_CONCENTRATION t  /*  hydro model water concentration, 
                                    in croco it is stored in an expanded 
                                    temperature (t) array */

/*=====================================================================
!  Counter variables for the internal loops of BIOLink
!======================================================================*/

# if defined BIOLink_UPDATE_CONCBIO 
#   define IRANGE1 1,nv_adv /* Index for the tracer variable,
                               first loop of BIOLink update_concbio */
#   define IRANGE2 1,N /* Index for the vertical levels,
                          second loop of BIOLink_update_concbio !*/
#   define IRANGE3 jfirst,jlast /* Index for the meridional direction, 
                                   third loop of BIOLink_update_concbio */
#   define IRANGE4 ifirst,ilast /* Index for the zonal direction,
                                   Fourth loop of the BIOLink_update_concbio */
# endif /* BIOLink_UPDATE_CONCBIO */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define IVERIF_BIOLINK i_BIOLink_verif /* Counter for the zonal direction in BIOLink */
#   define JVERIF_BIOLINK j_BIOLink_verif /* Counter for the meridional direction in BIOLink */
# endif /* PEPTIC/BLOOM/METeOR */

/*=====================================================================
!  Indexes for navigating the arrays of tracer variables
!======================================================================*/

# if defined BIOLink_UPDATE_CONCBIO 
#   define WATCONC_INDEX_EQ i4,i3,i2,nnew,itemp+ntrc_salt+i1 /* index for the concentration
                                                                with the BIOLink counter 
                                                                variables i4,i3,i2 */
#   define BIOSKSC_INDEX_EQ i4,i3,i2,i1                     /* index for sink and source terms
                                                               with the BIOLink counter 
                                                               variables i4,i3,i2 */
# endif /* BIOLink_UPDATE_CONCBIO */

# define STATE_VAR_INDEXkij 1:nv_state,k,i,j /* Index of the state variable */

/*=====================================================================
! Vertical indexes to navigate the water column
!======================================================================*/

# define LOOPK_SURF_TO_BOTTOM_WAT kmaxmod,1,-1 /* From bottom to surface */
# define LOOPK_SUBSURF_TO_BOTTOM_WAT kmaxmod-1,1,-1 /* From subsurface to bottom */
# define ABOVE_K k+1 /* Above the k level */

/*=====================================================================
! Tables related to the advected variables
! There does not appear to be a cvadv_wat table
!======================================================================*/

# define ARRAY_WATER_CONC GLOBAL_2D_ARRAY,N,3,NT /* Size of the array to store 
                                                    the concentration of tracer Variables, 
                                                    including salinity and temperature.
                                                    I guess the '3' is a time dimensionnality */

# define ARRAY_WATER_CONC0 GLOBAL_2D_ARRAY,N               /* Size of array to store a 
                                                              concentration, but of what ? */

# define WAT_CONCADV_ivkij t(i,j,k,nnew,itemp+ntrc_salt+iv)  /* Concentration for the index 
                                                                i,j,k,nnew and the variable
                                                                itemp+ntrc_salt+iv */
# define ADV_VAR_INDEXkij i,j,k,1:nv_adv /* Index of tracer variables */

/*=====================================================================
! Tables related to the advected variables corrected to have positive
! concentrations.
!======================================================================*/

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define WATCONCPOS cvadv_wat_pos /* Positive concentration of tracer variables */

# elif defined ECO3M

#   define WATCONCPOS  VAR(:)%conc(:,:,:)       /* Table  of positive tracer 
                                                   concentrations, the dimensions 
                                                   are (iv,i,j,k) */
#   define WATCONCPOS_ij  VAR(:)%conc(i,j,:)    /* Table of positive tracer 
                                                   concentrations, the dimension 
                                                   are (iv,k), or the tracer
                                                   indexes and the depth */
#   define WATCONCPOS_ivijk  VAR(iv)%conc(i,j,k) /* Positive concentration of the
                                                    tracer iv at index i,j,k*/
# endif  /* PEPTIC/BLOOM/METeOR/ECO3M */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define WATCONCPOS_ij  cvadv_wat_pos(:,:,i,j)     /* Table of positive concentration 
                                                        of the tracer variables at index 
                                                        i,j. The dimensions are (iv,k) */
#   define WATCONCPOS_ivijk  cvadv_wat_pos(iv,k,i,j) /* Positive concentration of 
                                                        tracer variables at index k,i,j */
# elif defined ECO3M
#   define WATCONCPOS_ij  VAR(:)%conc(i,j,:)    /* Table of positive tracer concentrations, 
                                                   the dimension are (iv,k), or the tracer
                                                   indexes and the depth */
#   define WATCONCPOS_ivijk  VAR(iv)%conc(i,j,k) /* Positive concentration of the tracer 
                                                    iv at index i,j,k*/
# endif  /* PEPTIC/BLOOM/METeOR/ECO3M */

/*=====================================================================
! Tables related to the sink and source terms of the advected variables 
!======================================================================*/

# define ARRAY_SINKSOURCES PROC_IN_ARRAY,NB_LAYER_WAT,nv_adv    /*  Size of array to store the sink
                                                                    and source terms computed by
                                                                    the biological model*/

# define BIO_SINKSOURCES BIO_SKSC_ADV /* Source and sink terms of the tracers 
                                         provided by the biological model */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   if defined key_MARS
#     define BIO_SKSC_ADV dcdt    /* Sink and source terms for the concentration 
                                     tracer variables */
#   else
#     define BIO_SKSC_ADV dcdt_adv    /* Sink and source terms for the concentration 
                                        of tracer variables */
#   endif/* key_MARS */

# elif defined ECO3M

#   define BIO_SKSC_ADV tend    /* Source and sink terms computed by BIOLink at index i,j */
# endif  /* PEPTIC/BLOOM/METeOR/ECO3M */

/*=====================================================================
! Tables related to the fixed variables
!======================================================================*/

# define FIXED_VAR_CONC cvfix_wat /* Hydro model concentration of fix variables. 
                                     For this a new array cvfix_wat had had to be 
                                     created */

# define FIXED_VAR_INDEXkij i,j,k,1:nv_fix /* Index of the fixed variables */
# define FIXED_VAR_INDEX_CUMULPROD i,j,k,ivfix_cumulprod_first-nv_adv:ivfix_cumulprod_last-nv_adv 
                                                   /* Index for the cumulative productions */

# define WAT_CONCFIX_ivkij cvfix_wat(i,j,k,iv)  /* Concentration for the index i,j,k of the fixed
                                                   variable of index iv */
# define WAT_CONCFIX_ifixkij cvfix_wat(i,j,k,iv)  /* Concentrtation for the index i,j,k of the fixed
                                                     variable of index iv with a different order */

/*=====================================================================
! Tables related to the positive concentration of fixed variables
!======================================================================*/

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define FIXCONCPOS cvfix_wat_pos /* Concentration of fixed varibles */
# endif /* PEPTIC/BLOOM/METeOR */

/*=====================================================================
! Tables related to the sink and source terms of fixed variables
!======================================================================*/

# define ARRAY_FIXED_SKSC nv_fix,NB_LAYER_WAT,GLOBAL_2D_ARRAY   /* Size of the array for the source
                                                                   and sink terms of the fixed 
                                                                   variables */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   if defined key_MARS
#     define BIO_SKSC_FIX dcdt   /*  Sink and source terms for the concentration 
                                     of fixed variables */
#   else
#     define BIO_SKSC_FIX dcdt_fix   /* Sink and source terms for the concentration
                                        of fixed variables */
#   endif/* key_MARS */
# endif /* PEPTIC/BLOOM/METeOR */


# define FIXED_SKSC_INDEXkij 1:nv_fix,k,i,j /* Index for the source and sink terms */

/*=====================================================================
! Tables related to the benthic variables
!======================================================================*/

# if defined BLOOM
#   define BENTH_INDEX_RANGE :,:,2:NB_LAYER_WAT,iv-nv_adv /* Index range of benthic variables */
#   define BENTH_INDEXij i,j,1,iv-nv_adv /* Index of benthic variables in the cvfix table */
# endif /* BLOOM */

# if defined BLOOM
#   define BENTHIC_CONCENTRATION cvfix_wat /* Concentration of fixed variable in the benthos. 
                                              Currently it is the same table as the one for 
                                              fixed variables */
# endif /* BLOOM */

# define WAT_CONCBENT_ivij c_bent(iv,i,j)  /* Concentration for the index i,j of the benthic 
                                              variable iv */

/*=====================================================================
! Tables related to the positive concentration of benthic variables
!======================================================================*/

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define BENTCONCPOS cv_bent_pos  /* Concentration of benthic variables */
# endif /* PEPTIC/BLOOM/METeOR */

/*************************************************************************/
/*************************************************************************/
/****************Variables for the internal routines**********************/
/**********************Or sediment model**********************************/
/*************************************************************************/
/*************************************************************************/

/*=====================================================================
!  Sinking velocity of tracers
!======================================================================*/

# define ARRAY_WAT_SETTL GLOBAL_2D_ARRAY,N,itsubs1:itsubs2 /* Size of the array to store the 
                                                              settling velocity of the BIOLink 
                                                              tracers */

# define WAT_SETTL ws_part /* Sinking velocity of particulate variables */
# define WAT_SETTL_ivkij ws_part(i,j,k,itemp+ntrc_salt+iv) /* Water settling velocities 
                                                              For the index i,j,k and the variable
                                                              itemp+ntrc_salt+iv */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define  WS_BIOLink ws3             /* Settling velocity computed by BIOLink at index i,j,k */
# endif /* PEPTIC/BLOOM/METeOR */

/*=====================================================================
!  Radiation variables
!======================================================================*/

# define ROOF_LAYER_RAD alumplafond /* Radiation at the top layer */
# define SOLAR_RAD srflx /* Incoming solar radiation */
# define RAD_SRFSCALE srf_scale /* Scaling factor = 1/(rho0*Cp) */
# define EXTINCTION_RAD extinction /* Coefficient of light 
                                      Extinction in the water column */

/*=====================================================================
!  BIOLink variables related to the photosynthetic available radiation
!======================================================================*/

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define SPMTOT_MGL cmes_3dmgl /* Concentration of suspended matter */
# endif /* PEPTIC/BLOOM/METeOR */

# ifdef PEPTIC
#   define PARAM_WAT_EXTINCT bd_fp%extincwat /* Extinction parameter */
#   define PARAM_CHLORO1_EXTINCT bd_fp%extincChl1 /* Parameter of extinction 
                                                     due to chlorophyll */
#   define PARAM_CHLORO2_EXTINCT bd_fp%extincChl2 /* Second parameter of extinction
                                                     due to chlorophyll */
#   define INIT_PEPTIC_PLCT_INDEX :,:,:,:,plct(i_plkt)%num_mod_mars(iquo) 
                                                  /* Index of plankton, I am not sure */
# endif /* PEPTIC */

# ifdef BLOOM
#   define PARAM_WAT_EXTINCT p_extincwat /* Extinction parameter */
#   define PARAM_CHLORO1_EXTINCT p_extincChl1 /* Extinction parameter due 
                                                 to chlorophyll */
#   define PARAM_CHLORO2_EXTINCT p_extincChl2 /* Second extinction parameter 
                                                 due to chlorophyll */
# endif /* BLOOM */

# if defined METeOR 
#   if ! defined PEPTIC && ! defined BLOOM
#     define PARAM_WAT_EXTINCT p_extincwat /* Extinction parameter */
#   endif /* PEPTIC/BLOOM */
# endif /* METeOR */

/*=====================================================================
!  Wind speed variables
!======================================================================*/

# ifdef BULK_FLUX
#   define WIND_SPEED wspd  /*  m.s-1 */
# else
#   define WIND_SPEED win10  /*  m.s-1 */
# endif /* BULK_FLUX */

# define ARRAY_WINDSPEED GLOBAL_2D_ARRAY /* The table to store the windspeed */

/*************************************************************************/
/*************************************************************************/
/*********************Variables from the hydro model**********************/
/*************************************************************************/
/*************************************************************************/

/*=====================================================================
! Height of the water column variables
!======================================================================*/

# define BATHY_H0 h /* Depth */
# ifdef WET_DRY
#   define RESIDUAL_THICKNESS_WAT D_wetdry
# else
#   define RESIDUAL_THICKNESS_WAT 0.
# endif /* Residual thickness if the wetting and drying is activated or not */
# define WATER_ELEVATION zeta /* Water elevation variations */
# define WATER_ELEV_ij zeta(i,j,1) /* Water elevation at the index i,j */

/*=====================================================================
!  Variables related to the vertical dimension of the water column
!======================================================================*/

# if defined ECO3M
#   define TOTAL_WATER_HEIGHT prof     /* Total height of the water column computed by BIOLink at
                                          index i,j */
#   define THICKLAYERWC dz /* Thickness of the vertical grid cells at index i,j */
# endif /* ECO3M */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define TOTAL_WATER_HEIGHT htot    /* Total height of the water column computed by BIOLink at
                                         index i,j */
#   define THICKLAYERWC thicklayerW_C  /* Thickness of the vertical grid cells measured from their
                                          center at index i,j */
#   define THICKLAYERWW thicklayerW_W  /* Thickness of the vertical grid cells measured from their
                                          upper faces at index i,j */
# endif /* PEPTIC/BLOOM/BIOLink */

/*=====================================================================
!  Current variables
!======================================================================*/

# define ZONAL_VELOCITY_UZ u /* Zonal velocity from the hydro model */
# define MERIDIONAL_VELOCITY_VZ v /* Meridional velocity from the hydro model */

# define INDEX_UVZBOTT_ij (1,i,j,1) /* Index for the bottom current */

/*=====================================================================
!  Bottom current variables
!======================================================================*/

/* If METeOR is used we compute a center approximation for the magnitude of bottom current */
# if defined METeOR
#   define BOTTOM_CURRENT_kmimjm 0.5_rsh*sqrt((ZONAL_VELOCITY_UZ(im-1,jm,km,nnew)+ZONAL_VELOCITY_UZ(im,jm,km,nnew))**2.+(MERIDIONAL_VELOCITY_VZ(im,jm,km,nnew)+MERIDIONAL_VELOCITY_VZ(im,jm-1,km,nnew))**2.)
# endif /* METeOR */

# if defined BIOLink_UPDATE_CONCBIO 

#   if defined BLOOM
#     if defined key_benthos
#       define BOTTOM_CURRENT_ij 0.5_rsh*sqrt((ZONAL_VELOCITY_UZ(i-1,j,1,nnew)+ZONAL_VELOCITY_UZ(i,j,1,nnew))**2.+(MERIDIONAL_VELOCITY_VZ(i,j,1,nnew)+MERIDIONAL_VELOCITY_VZ(i,j-1,1,nnew))**2.)
                              /* Interpolated velocity of the bottom current at index i,j */
#     endif /* key_benthos */
#     define BOTTCURRENTBIO bottom_current /* Velocity of the bottom current */
#   endif /*BLOOM*/

# endif /* BIOLink_UPDATE_CONCBIO */



/*=====================================================================
!  Kinetic energy variables
!======================================================================*/

# define CIN_TURBULENT_ENERGY  trb /* turbulent kinetic  energy  evaluated by  hydro model  (in m2/s2) */
# define ARRAY_CIN_TURB_ENERGY GLOBAL_2D_ARRAY,0:NB_LAYER_WAT,2,NGLS /* The array to store the 
                                                                        kinetic energy, and also
                                                                        some mixing variables I 
                                                                        guess */
# if defined GLS_MIXING
#   define ECT_kij  trb(i,j,k,1,1) /* Turbulent kinetic energy at index i,j,k  
                                      evaluated by the hydro model */
# else
#   define ECT_kij  0. /* Turbulent kinetic energy at index i,j,k evaluated 
                          by the hydro model  */
# endif /*GLS_MIXING */

/*=====================================================================
!  Density of seawater
!======================================================================*/

# define RHOREF rho0 /* Reference density */
# define WATER_DENSITY_kp1ij 0. /* water density at mesh i,j,k+1 if already known and evaluated by  hydro model  */
# define WATER_DENSITY_kij   0. /* water density at mesh i,j,k if already known and evaluated by  hydro model  */

/*=====================================================================
!  Variables related to the temperature and salinity
!======================================================================*/

# define SALHYDRO_ijk t(i,j,k,1,itemp+1) /* Salinity in the hydro model at the index i,j,k */
# define TEMPHYDRO_ijk t(i,j,k,1,itemp) /* Temperature in the hydro model at the index i,j,k */

# if defined ECO3M
#   define TEMP_BIOLink temp_bio       /* Temperature computed by BIOLink at index i,j,k  */
#   define SAL_BIOLink sal_bio         /* Salinity computed by BIOLink at index i,j,k  */
# endif /* ECO3M */

# if defined PEPTIC || defined BLOOM || defined METeOR 
#   define TEMP_BIOLink temp           /* Temperature computed by BIOLink at index k,i,j */
#   define SAL_BIOLink sal             /* Salinity computed by BIOLink at index i,j,k  */
# endif /* PEPTIC/BLOOM/METeOR */

#endif /* BIOLink */
