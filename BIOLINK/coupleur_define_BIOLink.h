
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
   !&E==========================================================================

# define BIOLink_HOST BIOLink_CROCO

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some specific commands used in MARS are to be commented out automatically 
!   for another hydro model (we use  define)
!   + There may be other definitions to add in order to replace compatible 
!      terms with the hydro coupled model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

!!#if ! defined MUSTANG
!/* not defined these variables if already defined with module MUSTANG */

#if ! defined key_MARS
!/* MPI */
# define CALL_MPI !call_MPI
# define IF_MPI !if
# define USE_MPI !use
# define ENDIF_MPI !endif
# define PRINT_DBG !print
# define IF_AGRIF !if
# define OMPMPI  !OMPMPI
!/*CROCO */
# define iscreenlog stdout
# define ierrorlog stdout
# define iwarnlog stdout
# define NAME_SUBS vname(1,indxT+ntrc_salt+isubs)
! MARS : #define NAME_SUBS name_var(irk_fil(isubs))

!/* end not MARS */
#endif

!/* end not module MUSTANG */
!!#endif

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*         Directory where are namelists files  
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
# define REPFICNAMELIST 'FIC_NAMELIST'

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         VOCABULARY     Modele HYDRO   TO   BIOLink
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

!/* Spatial Grid limits definition  of loops
!*   inside the domain - except meshes at open boundaries
!*/
# define IMIN_BOUCL 1
# define IMAX_BOUCL Lm
# define JMIN_BOUCL 1
# define JMAX_BOUCL Mm

!/* dimensions table definition 
!*/
# define COMPLETE_ARRAY      1   :Lm    , 1   :Mm
# define PROC_IN_ARRAY       GLOBAL_2D_ARRAY   
# define PROC_IN_ARRAY_m1p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m1p1  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m2p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_0p1   GLOBAL_2D_ARRAY

!/* total number of state variables (need for BLOOM and PEPTIC when adding variables)    
!*/
# define NBVARADV_TOT  ntrc_subs  

!/* general variable hydro , bathy, time ... defined in hydro model and used  by BIOLink 
!*/
# define NB_LAYER_WAT N
# define BATHY_H0 h
# ifdef WET_DRY
# define RESIDUAL_THICKNESS_WAT D_wetdry
# else
# define RESIDUAL_THICKNESS_WAT 0.
# endif
# define WATER_ELEVATION zeta
# define CELL_SURF surf_cell
# define TRANSPORT_TIME_STEP dt 
# define CURRENT_TIME time
# define TIME_BEGIN tdeb
!# define SALINITY_MOD sal
!# define TEMPERATURE_MOD temp
# define WAT_SETTL ws_part
# define WATER_CONCENTRATION t  /* water concentration in hydro model (=cv_wat in MARS)*/
# define BIO_SINKSOURCES BIO_SKSC_ADV
# define ECO_TIME_STEP dt_bio_update
# define ROOF_LAYER_RAD alumplafond
# define SOLAR_RAD srflx
# define RAD_SRFSCALE srf_scale /* facteur d echelle= 1/(rho0*Cp)in CROCO available in forces.h (=1.0_rsh in MARS)*/
# define EXTINCTION_RAD extinction
# define FIXED_VAR_CONC cvfix_wat
# define SUFFIX_RUNFILE 1
# define NUMBER_PI pi
#ifdef BULK_FLUX
# define WIND_SPEED wspd  /*  m.s-1 */
#else
# define WIND_SPEED win10  /*  m.s-1 */
#endif
# define ZONAL_VELOCITY_UZ u
# define MERIDIONAL_VELOCITY_VZ v
# define LONGITUDE lonr
# define LATITUDE latr
# define CIN_TURBULENT_ENERGY  trb /* turbulent kinetic  energy  evaluated by  hydro model  (in m2/s2) */
# define RHOREF rho0 

!/* dimensions of variables in hydro modele used by BIOLink 
!*/
# define ARRAY_WAT_SETTL GLOBAL_2D_ARRAY,N,itsubs1:itsubs2
# define ARRAY_WATER_CONC GLOBAL_2D_ARRAY,N,3,NT
# define ARRAY_WATER_CONC0 GLOBAL_2D_ARRAY,N
# define ARRAY_SINKSOURCES PROC_IN_ARRAY,NB_LAYER_WAT,nv_adv    /*  CROCO - resolution dans BIOLink ou non !*/
# define ARRAY_WINDSPEED GLOBAL_2D_ARRAY
# define ARRAY_CIN_TURB_ENERGY GLOBAL_2D_ARRAY,0:NB_LAYER_WAT,2,NGLS
# define ARRAY_FIXED_SKSC nv_fix,NB_LAYER_WAT,GLOBAL_2D_ARRAY
!

!/* expressions or indexes
!*/

# define STATE_VAR_INDEXkij 1:nv_state,k,i,j
# define ADV_VAR_INDEXkij i,j,k,1:nv_adv
# define WAT_SETTL_ivkij ws_part(i,j,k,itemp+ntrc_salt+iv)
# define WAT_CONCADV_ivkij t(i,j,k,nnew,itemp+ntrc_salt+iv)  /* variables d etat advectees in hydro model into the mesh i,j,k ATTENTION order of indexes */
# define WATER_ELEV_ij zeta(i,j,1)
# define WAT_CONCFIX_ivkij cvfix_wat(i,j,k,iv)  /* variables d etat fixees in hydro model into the mesh i,j,k ATTENTION order of indexes */
# define WAT_CONCFIX_ifixkij cvfix_wat(i,j,k,iv)  /* variables d etat fixees in hydro model into the mesh i,j,k ATTENTION order of indexes */
# define WAT_CONCBENT_ivij c_bent(iv,i,j)  /* variables d etat fixees in hydro model into the mesh i,j,k ATTENTION order of indexes */
# define SALHYDRO_ijk t(i,j,k,1,itemp+1) /* salinity in hydro model into the mesh i,j,k ATTENTION order of indexes */
# define TEMPHYDRO_ijk t(i,j,k,1,itemp) /* temperature in hydro model into the mesh i,j,k ATTENTION order of indexes */
# define INDEX_UVZBOTT_ij (1,i,j,1) /* indexes order to define courrant in bottom layer in hydro model into the mesh i,j ATTENTION order of indexes */
# define LOOPK_SURF_TO_BOTTOM_WAT kmaxmod,1,-1
# define LOOPK_SUBSURF_TO_BOTTOM_WAT kmaxmod-1,1,-1
# define ABOVE_K k+1
# define TEST_NOT_RESTART_SUBS .NOT.l_restart_subs
# define TEST_NOT_INITFROMFILE .NOT.l_initfromfile
# define WATER_DENSITY_kp1ij 0. /* water density at mesh i,j,k+1 if already known and evaluated by  hydro model  */
# define WATER_DENSITY_kij   0. /* water density at mesh i,j,k if already known and evaluated by  hydro model  */
#if defined GLS_MIXING
# define ECT_kij  trb(i,j,k,1,1)
 /* turbulent kinetic  energy at mesh i,j,k evaluated by  hydro model  */
#else
# define ECT_kij  0. /* turbulent kinetic  energy at mesh i,j,k evaluated by  hydro model  */
#endif
# define FIXED_VAR_INDEXkij i,j,k,1:nv_fix /* could be not the same index order or not the same numbre varibale index  if not MARS */
# define FIXED_VAR_INDEX_CUMULPROD i,j,k,ivfix_cumulprod_first-nv_adv:ivfix_cumulprod_last-nv_adv /* could be not the same index order or not the same numbre varibale index  if not MARS */
# define FIXED_SKSC_INDEXkij 1:nv_fix,k,i,j

#if defined BIOLink_UPDATE_CONCBIO 
!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   info for resolution equation in BIOLink (see coupl_BIOLink2hydro in coupleur_BIOLink.F90) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
# define IRANGE1 1,nv_adv /* first loop for resolution bio equation (index 4 )  :  variable index for CROCO!*/
# define IRANGE2 1,N                                      /* second loop for resolution bio equation (index 3 ) : number of water layer for CROCO !*/
# define IRANGE3 jfirst,jlast                    /* third loop for resolution bio equation (index 2 )  : j grid mesh !*/
# define IRANGE4 ifirst,ilast                    /* last loop for resolution bio equation (index 1 )  : i grid mesh !*/
# define WATCONC_INDEX_EQ i4,i3,i2,nnew,itemp+ntrc_salt+i1                /* index for water concentration  !*/
# define BIOSKSC_INDEX_EQ i4,i3,i2,i1                     /* index for bio sink sources  !*/
#endif

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         VOCABULARY     BIOLink  to module BIO 
!    (arrays which are stored in comBIOLink after conversion because not the same index order or number)
!     ( or other variable which module BIO need and evaluated bt BIOLink)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#if defined ECO3M
# define BIO_TIME_STEP dt_bio
# define DT_CONSERV_BIOLINK dt_bio_conserv
# define TEMP_BIOLink temp_bio       /* temperature array in BIOLink module index order if k,i,j or same as hydro model  */
# define SAL_BIOLink sal_bio         /* salinity array in BIOLink module index order if k,i,j   */

# define TOTAL_WATER_HEIGHT prof     /* hauteur d eau totale calculee par BIOLink dans le temps i,j */
# define BIO_SKSC_ADV tend    /* puits et sources bio calcules par BIOLink dans le temps i,j (pour toutes les variables d etat advectees ou non) = BIO_SINKSOURCES si meme tableau */

! verifier que dz est bien centre en c ? oui
#define THICKLAYERWC dz
! est ce que dz en w est necessaire?
!#define THICKLAYERWW thicklayerW_W
!#define SPMTOT_MGL cmes_3dmgl
# define WATCONCPOS_ivijk  VAR(iv)%conc(i,j,k)
# define WATCONCPOS_ij  VAR(:)%conc(i,j,:)
# define WATCONCPOS  VAR(:)%conc(:,:,:)

#endif  /* ECO3M */

#if defined PEPTIC || defined BLOOM || defined METeOR 
# define BIO_TIME_STEP dtbio
# define DT_CONSERV_BIOLINK dt_conserv_BIO
# define IVERIF_BIOLINK i_BIOLink_verif 
# define JVERIF_BIOLINK j_BIOLink_verif 
# define TEMP_BIOLink temp           /* temperature array in BIOLink module index order if k,i,j or same as hydro model  */
# define SAL_BIOLink sal             /* salinity array in BIOLink module index order if k,i,j   */
# define  WS_BIOLink ws3             /* settling rate array in BIOLink module index order iv,k,i,j  */
# define TOTAL_WATER_HEIGHT htot    /* hauteur d eau totale calculee par BIOLink dans le temps i,j */

#ifdef key_MARS

# define BIO_SKSC_FIX dcdt   /*  sinksources for fixed variables (=dcdt for MARS)  index : FIXED_SKSC_ALLINDEX*/
# define BIO_SKSC_ADV dcdt    /* puits et sources bio calcules par BIOLink dans le temps i,j (pour toutes les variables d etat advectees ou non) = BIO_SINKSOURCES si meme tableau */

#else
/* CROCO */
# define BIO_SKSC_FIX dcdt_fix   /*  sinksources for fixed variables (= dcdt_fixedvar : BIOLink array) index : FIXED_SKSC_ALLINDEX*/
# define BIO_SKSC_ADV dcdt_adv    /* puits et sources bio calcules par BIOLink dans le temps i,j (pour toutes les variables d etat advectees ou non) = BIO_SINKSOURCES si meme tableau */

#endif

#define THICKLAYERWC thicklayerW_C
#define THICKLAYERWW thicklayerW_W
#define SPMTOT_MGL cmes_3dmgl
#define WATCONCPOS cvadv_wat_pos 
#define FIXCONCPOS cvfix_wat_pos 
#define BENTCONCPOS cv_bent_pos
# define WATCONCPOS_ivijk  cvadv_wat_pos(iv,k,i,j)
# define WATCONCPOS_ij  cvadv_wat_pos(:,:,i,j)

#endif

#ifdef PEPTIC
# define PARAM_WAT_EXTINCT bd_fp%extincwat
# define PARAM_CHLORO1_EXTINCT bd_fp%extincChl1
# define PARAM_CHLORO2_EXTINCT bd_fp%extincChl2
# define INIT_PEPTIC_PLCT_INDEX :,:,:,:,plct(i_plkt)%num_mod_mars(iquo)
!(MARS)# define INIT_PEPTIC_PLCT_INDEX :,:,:, plct(i_plkt)%num_mod_mars(iquo)
#endif

#ifdef BLOOM
# define BENTHIC_CONCENTRATION cvfix_wat
# define BENTH_INDEX_RANGE :,:,2:NB_LAYER_WAT,iv-nv_adv
# define BENTH_INDEXij i,j,1,iv-nv_adv
# define PARAM_WAT_EXTINCT p_extincwat
# define PARAM_CHLORO1_EXTINCT p_extincChl1
# define PARAM_CHLORO2_EXTINCT p_extincChl2
#ifdef key_benthos
# define BOTTOM_CURRENT_ij 0.5_rsh*sqrt((ZONAL_VELOCITY_UZ(i-1,j,1,nnew)+ZONAL_VELOCITY_UZ(i,j,1,nnew))**2.+(MERIDIONAL_VELOCITY_VZ(i,j,1,nnew)+MERIDIONAL_VELOCITY_VZ(i,j-1,1,nnew))**2.)
#endif
#define BOTTCURRENTBIO bottom_current
#endif

#if defined METeOR 
#if ! defined PEPTIC && ! defined BLOOM
# define PARAM_WAT_EXTINCT p_extincwat
#endif
# define BOTTOM_CURRENT_kmimjm 0.5_rsh*sqrt((ZONAL_VELOCITY_UZ(im-1,jm,km,nnew)+ZONAL_VELOCITY_UZ(im,jm,km,nnew))**2.+(MERIDIONAL_VELOCITY_VZ(im,jm,km,nnew)+MERIDIONAL_VELOCITY_VZ(im,jm-1,km,nnew))**2.)
#endif


#endif
