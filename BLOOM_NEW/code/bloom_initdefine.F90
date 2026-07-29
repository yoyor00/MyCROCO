
!---------------------------------------------------------------------------
!
                     MODULE bloom_initdefine
!
!---------------------------------------------------------------------------
   

#include "cppdefs.h"

#if defined SUBSTANCE && defined BIOLink && defined BLOOM

  !!======================================================================
  !!                   ***  MODULE  bloom_initdefine  ***
  !! Ocean dynamics Bio :  Initializations, reading of files *.dat (rivers, discharges, bio variables...)
  !!
  !!   History :
  !!    !  2019-08 (B. Thouvenin) issued from bioloinit and pepticinit for portability adaptation
  !!======================================================================

#include "coupleur_define_BIOLink.h"

  USE comBIOLink
  USE comBIOLink_physics
  USE comBIOLink_helping
  USE comsubstance, ONLY : lchain,rsh,rlg
  USE module_BIOLink
  USE COMBLOOM

  IMPLICIT NONE

  !! * Accessibility
  PUBLIC bloom_param                   ! routine called by init.F90
  PUBLIC bloom_init_iv,bloom_init_id   ! routine called by subreaddat.F90 or initBIOLink.F90


  !! * Private variables
   INTEGER                  :: idmin=99999,idmax=0, &  ! to get min and max indexes of diagnostic var.
                               ivmin=99999,ivmax=0     ! to get min and max indexes of substance var.
 
  CONTAINS
  !!======================================================================


  SUBROUTINE bloom_param(rw)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE bloom_param  ***
  !&E
  !&E ** Purpose : Read and write namelist relative to the module of bloom
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : init
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E       !  2009-10  (V. Garnier)  Original code
  !&E       !  2015-09  (M. Arancio)  Adapation for Quota (darwin model)
  !&E
  !&E---------------------------------------------------------------------
  !! * Modules used

  !! * Arguments
  CHARACTER(LEN=1), INTENT( in ) :: rw

  !! * Local declarations
  CHARACTER(LEN=lchain) :: filepc


  INTEGER            :: iso,IERR_MPI
  INTEGER           :: lstr,lenstr

#ifdef key_N_tracer
   CHARACTER(LEN=19) :: date_start_tracerN
   REAL(KIND=rlg)        :: tool_datosec
#elif defined key_P_tracer
   CHARACTER(LEN=19) :: date_start_tracerP
   REAL(KIND=rlg)        :: tool_datosec
#endif

   NAMELIST/namBIOLink/l_bioretro_extinct,dt_bio_update,l_waterdensity_known,  &
                     i_BIOLink_verif,j_BIOLink_verif,DT_CONSERV_BIOLINK
   NAMELIST/namoptions/ l_physadaptation,l_filtbenthsinus,l_filtbenthmes,l_SNeffect_settle,l_phyzoodeteffect_settle,  &
                        l_ChlNratio_var,p_sali_thhold_bio,l_bioretro_extinct
#ifdef key_BLOOM_opt2
   NAMELIST/namorgmat/ p_N_remin,p_P_remin,p_BSi_dissEau,p_T_effect,p_det_fragm,&
                       p_diss_regmod,p_micz_diss
   NAMELIST/namphosphor/ p_P_adsor,p_P_desor,p_P_adsormaxspim,p_P_adsormaxsed
#else
   !NAMELIST/namorgmat/ p_N_remin,p_P_remin,p_Si_diss,p_T_effect,p_reminbenth
   NAMELIST/namorgmat/ p_N_remin,p_P_remin,p_T_effect,p_kO2_reminO2,p_kO2_nit,p_nitrif,  &
                       p_DNO3_denit,p_kNO3_reminssO2,p_kiO2_remin0O2,p_kiO2_denit,p_kiNO3_remin0O2,     &
                       p_kO2_reoxyd,p_ODU_oxy,p_ODU_precip,p_N_reminR,p_P_reminR,          &
                       p_GO2_Norg,p_GODU_Norg,p_GNO3_Norg,p_GO2_NorgR,p_GODU_NorgR,p_GNO3_NorgR,p_GO2_NH4,   &
                       p_aging_MO,p_burried,      &
!                       p_Si_diss,p_kiO2_precSi,p_Si_precip
                       p_Si_Eq,p_Si_EqPrec,p_BSi_dissEau,p_BSi_dissSurfSed,p_BSi_dissFondSed,   &
                       p_T_effectSi,p_k_remin,p_xflimz,p_kSi,p_Si_precip
   NAMELIST/namphosphor/ p_P_speedup_reminanaer,p_P_adsor,p_P_desor,p_P_adsormaxspim, &
                         p_P_adsormaxsed,p_P_precFeO2,p_kO2_precPFe,p_kNO3_precPFe,p_kiO2_dissPFe,p_kiO2_desorP,  &
                         p_kiNO3_dissPFe,p_P_precFeNO3,p_P_dissFe,p_GO2_PFe,p_GNO3_PFe
#endif
   NAMELIST/namgenphyto/ p_phyto_ChlNratio,p_phyto_SiNratio,p_phyto_NPratio,p_phyto_CNratio, &
                         p_phyto_ChlNratiomax,p_phyto_ChlN_ksmithextinct
   NAMELIST/namnanophyto/ p_nano_mumax,p_nano_kNO3,p_nano_kNH4,p_nano_kPO4,p_nano_mort, &
                          p_nano_iksmith,p_nano_thhold_mort
   NAMELIST/namdiatom/ p_diat_mumax,p_diat_kNO3,p_diat_kNH4,p_diat_kSi,p_diat_kPO4,     &
                       p_diat_iksmith,p_diat_mort,p_diat_thhold_mort
   NAMELIST/namdino/ p_dino_mumax,p_dino_thhold_ect,p_dino_kNO3,p_dino_kNH4,p_dino_kPO4,p_dino_mort,      &
                     p_dino_iksmith,p_dino_thhold_mort
   NAMELIST/nammesozoo/ p_mesz_thrN,p_mesz_kivlev,p_mesz_mumax,p_mesz_assim,      &
                        p_mesz_excret,p_mesz_mort1,p_mesz_mort2,p_zoo_CDWratio,    &
                        p_zoo_CNratio,p_mesz_thhold_mes_kivlev,p_mesz_thhold_mort
   NAMELIST/nammicrozoo/ p_micz_mumax,p_micz_kgraz,p_micz_assim,p_micz_excret,p_micz_mort,   &
                         p_micz_thrnano,p_micz_thhold_mort

#ifdef key_BLOOM_opt2
   NAMELIST/namoxygen/ p_phyto_photoratio,p_kO2_reminO2,p_nitrif,p_phyto_resp,p_zoo_resp
#else
   NAMELIST/namoxygen/ p_phyto_photoratio,p_phyto_resp,p_zoo_resp,p_KO2sed_aeration,p_Kzsed_aeration
#endif
   NAMELIST/namoptics/ p_extincwat,p_extincspim,p_extincChl1,p_extincChl2,p_parradratio
#if defined MUSTANG
   NAMELIST/namsediment/ p_diat_mort_sed, p_detzoo_wsed,p_detphy_wsed
#ifdef key_sedim_MPB
   NAMELIST/namMPB/p_Tmax_MPB,p_Topt_MPB,p_beta_temp,p_K_PARsed,p_KuptC_MPBN,p_gamma_P,  &
                   p_EPS_leaching,p_uptN_alpha2,p_KuptN_DIN,p_KuptN_MPBC,p_migexu_alpha1,   &
                   p_mig_delta_NC,p_prop_resp_MPBC,p_Tmax_Bact,p_Topt_Bact,p_uptN_Bact_alpha3, &
                   p_KuptBactN_DON,p_ratio_uptakeCsN,p_uptBactN_EPS_beta1,p_morta_MPB,  &
                   p_morta_Bact,p_qN2C_migupMPB,p_qC2N_migdwnMPB,p_Kext_sed,p_ratioC_Chla, &
                   p_ratio_nu_ud,p_nu_up,p_GO2_respMPB,attenu_w,xK_migd_N,ratioC_Chla,kd_chla,  &
                   seuil_NH4_uptakNH4only,seuil_NH4_uptakNH4NO3,Chlamax_m2microm,          &
                   ratio_uptakeNO3sNH4,QR_diat
#endif
#endif
   NAMELIST/namgrazing/ p_mesz_captdiat,p_mesz_captdino,p_mesz_captmicz,p_micz_captdiat,      &
                        p_micz_captnano,p_micz_captdet,p_micz_captdino,p_txfiltbenthmax



  !!----------------------------------------------------------------------
  !! * Executable part

IF(rw == 'r')THEN

! Lecture des namelists
!==========================

#ifdef key_BLOOM_opt2
     filepc=REPFICNAMELIST2//'/parabloom_opt2.txt'
     IF_AGRIF (.NOT. Agrif_Root()) filepc='./parabloom_opt2'//TRIM(Agrif_Cfixed())//'.txt'
#else
     lstr = lenstr(parafilename)
     filepc = parafilename(1:lstr)
     !IF_AGRIF (.NOT. Agrif_Root()) filepc='./para_BLOOM'//TRIM(Agrif_Cfixed())//'.txt'
     IF_AGRIF (.NOT. Agrif_Root()) filepc=filepc//TRIM(Agrif_Cfixed())//'.txt'
#endif


   OPEN(50,file=filepc,status='old',form='formatted', access='sequential')
   READ(50,namBIOLink)
   READ(50,namoptions)
   READ(50,namorgmat)
   READ(50,namphosphor)
   READ(50,namgenphyto)
   ! facteurs de conversion
   rapsiaz=p_phyto_SiNratio
   rappaz=1.0_rsh/p_phyto_NPratio

   READ(50,namnanophyto)
   READ(50,namdiatom)
   READ(50,namdino)

   READ(50,nammesozoo)
   READ(50,nammicrozoo)
   READ(50,namoxygen)
   READ(50,namoptics)

#if defined MUSTANG
   READ(50,namsediment)
#ifdef key_sedim_MPB
   READ(50,namMPB)
   invMPBCmaxm3=12.e-3/(45._rsh*Chlamax_m2microm) *1.e-3  ! divise par 1000 passer de mole en mmole
   MPBCmaxm3=1._rsh/invMPBCmaxm3
   QR_diat_inv=1.0_rsh/QR_diat
#endif
#endif

   READ(50,namgrazing)


  CLOSE(50)                                  ! end reading of namelists

ELSE   !IF(rw == 'w')


! writing of namelists and different simulated plankton species
!================================================================


   !:::::::::: module BLOOM ::::::::::::

   !IF_MPI (MASTER) THEN
     MPI_master_only  WRITE(iscreenlog,namBIOLink)
     MPI_master_only  WRITE(iscreenlog,namoptions)
     MPI_master_only  WRITE(iscreenlog,namorgmat)
     MPI_master_only  WRITE(iscreenlog,namphosphor)
     MPI_master_only  WRITE(iscreenlog,namgenphyto)
     MPI_master_only  WRITE(iscreenlog,namnanophyto)
     MPI_master_only  WRITE(iscreenlog,namdiatom)
     MPI_master_only  WRITE(iscreenlog,namdino)
     MPI_master_only  WRITE(iscreenlog,nammesozoo)
     MPI_master_only  WRITE(iscreenlog,nammicrozoo)
     MPI_master_only  WRITE(iscreenlog,namoxygen)
     MPI_master_only  WRITE(iscreenlog,namoptics)

#if defined MUSTANG
     MPI_master_only WRITE(iscreenlog,namsediment)
#ifdef key_sedim_MPB
     MPI_master_only WRITE(iscreenlog,namMPB)
#endif
#endif

     MPI_master_only WRITE(iscreenlog,namgrazing)

  ! ENDIF_MPI

ENDIF


  END SUBROUTINE bloom_param

  !!======================================================================
  SUBROUTINE bloom_init_iv(iv,standname,icall)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_init_iv  ***
   !&E
   !&E ** Purpose : Initialize indexes of state variables (substances only)
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : sub_read_var
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2009-10  (V. Garnier)  Original code
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   INTEGER,          INTENT( in ) :: iv
   INTEGER,INTENT( in ),OPTIONAL  :: icall
   CHARACTER(LEN=*), INTENT( in ) :: standname

   !! * Local declarations
   INTEGER               :: sumverif,IERR_MPI

   !!----------------------------------------------------------------------
   !! * Executable part

   SELECT CASE( TRIM(ADJUSTL(ADJUSTR(standname))) )

     CASE('mole_concentration_of_nitrate_in_sea_water')
       iv_nutr_NO3 = irk_mod(iv)
     CASE('mole_concentration_of_silicate_in_sea_water')
       iv_nutr_SiOH = irk_mod(iv)
     CASE('mole_concentration_of_phosphate_in_sea_water')
       iv_nutr_PO4 = irk_mod(iv)
     CASE('mole_concentration_of_sorbed_phosphate_in_sea_water')
       iv_nutr_Pads = irk_mod(iv)
     CASE('mole_concentration_of_diatoms_expressed_as_nitrogen_in_sea_water')
       iv_phyto_diat_N = irk_mod(iv)
     CASE('mole_concentration_of_dinoflagellates_expressed_as_nitrogen_in_sea_water')
       iv_phyto_dino_N = irk_mod(iv)
     CASE('mole_concentration_of_organic_detritus_expressed_as_nitrogen_in_sea_water')
       iv_detr_N = irk_mod(iv)
     CASE('mole_concentration_of_organic_detritus_expressed_as_silicon_in_sea_water')
       iv_detr_Si = irk_mod(iv)
     CASE('mole_concentration_of_organic_detritus_expressed_as_phosphorus_in_sea_water')
       iv_detr_P = irk_mod(iv)
     CASE('mole_concentration_of_ammonium_in_sea_water')
       iv_nutr_NH4 = irk_mod(iv)
#if defined key_BLOOM_opt2
     CASE('mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water')
#else
     CASE('mole_concentration_of_mesozooplankton_expressed_as_nitrogen_in_sea_water')
#endif
       iv_zoo_meso_N = irk_mod(iv)
     CASE('mole_concentration_of_nanoplankton_expressed_as_nitrogen_in_sea_water')
       iv_phyto_nano_N = irk_mod(iv)
#if defined key_BLOOM_opt2
     CASE('mole_concentration_of_microzooplankton_expressed_as_carbon_in_sea_water')
#else
     CASE('mole_concentration_of_microzooplankton_expressed_as_nitrogen_in_sea_water')
#endif
       iv_zoo_micr_N = irk_mod(iv)
     CASE('mass_concentration_of_suspended_matter_in_sea_water')
       iv_spim = irk_mod(iv)
#if defined key_BLOOM_opt2
     CASE('cumulative_diatom_production_expressed_as_nitrogen_in_sea_water')
#else
     CASE('cumulative_diatom_production_expressed_as_carbon_in_sea_water')
     MPI_master_only PRINT*, 'cumulative_diatom_production_expressed_as_carbon_in_sea_water'
#endif
       iv_phyto_diat_pp = irk_mod(iv)
       MPI_master_only PRINT*,'Jai passé le case'
#if defined key_BLOOM_opt2
     CASE('cumulative_dinoflagellate_production_expressed_as_nitrogen_in_sea_water')
#else
     CASE('cumulative_dinoflagellate_production_expressed_as_carbon_in_sea_water')
#endif
       iv_phyto_dino_pp = irk_mod(iv)
#if defined key_BLOOM_opt2
     CASE('cumulative_nanoflagellate_production_expressed_as_nitrogen_in_sea_water')
#else
     CASE('cumulative_nanoflagellate_production_expressed_as_carbon_in_sea_water')
#endif
       iv_phyto_nano_pp = irk_mod(iv)


#ifdef key_BLOOM_insed
     CASE('oxygen_demand_anaerobic_unit')
       iv_ODU = irk_mod(iv)
     CASE('iron_bound_phosphorus')
       iv_PFe = irk_mod(iv)
     CASE('mole_concentration_of_refractory_organic_detritus_expressed_as_nitrogen_in_sea_water')
       iv_detrR_N = irk_mod(iv)
     CASE('mole_concentration_of_refractory_organic_detritus_expressed_as_phosphorus_in_sea_water')
       iv_detrR_P = irk_mod(iv)
#endif


     CASE('dissolved_oxygen_in_water_column')
       iv_oxygen = irk_mod(iv)

#ifdef key_BLOOM_opt2
    !insertionVSjuin2010 3 nvelles variables dissoutes
     CASE('mole_concentration_of_dissolved_nitrogen')
       iv_diss_N = irk_mod(iv)
     CASE('mole_concentration_of_dissolved_phosphate')
       iv_diss_P = irk_mod(iv)
     CASE('mole_concentration_of_dissolved_silicate')
       iv_diss_Si = irk_mod(iv)
     !fin insertion
     !InsertionVSnov2010 2 nvelles variables pr equilibre au fond
     CASE('mole_concentration_of_dissolved_N2')
       iv_diss_fond_Nitr = irk_mod(iv)
      !fin insertion Nov 2010
#endif     

    CASE DEFAULT

   END SELECT



   ! smallest and largest indexes
   IF (iv >= ivmax) ivmax=iv
   IF (iv <= ivmin) ivmin=iv


   IF ( iv==nv_state) THEN
     IF_MPI (MASTER) THEN
       MPI_master_only WRITE(iscreenlog,*) ' '
       MPI_master_only WRITE(iscreenlog,*) '********************************'
       MPI_master_only WRITE(iscreenlog,*) '          bloom_INIT_IV         '
       MPI_master_only WRITE(iscreenlog,*) ' Initialization of iv_variables '
       MPI_master_only WRITE(iscreenlog,*) '********************************'
       MPI_master_only WRITE(iscreenlog,*) ' '
       MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_NH4 = ',iv_nutr_NH4
       MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_NO3 = ',iv_nutr_NO3
       MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_SiOH = ',iv_nutr_SiOH
       MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_PO4 = ',iv_nutr_PO4
       MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_Pads = ',iv_nutr_Pads
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_nano_N = ',iv_phyto_nano_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_diat_N = ',iv_phyto_diat_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_dino_N = ',iv_phyto_dino_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_meso_N = ',iv_zoo_meso_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_micr_N = ',iv_zoo_micr_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_detr_N = ',iv_detr_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_detr_Si = ',iv_detr_Si
       MPI_master_only WRITE(iscreenlog,*) 'iv_detr_P = ',iv_detr_P
       IF(iv_spim == 0 .and. nv_mud .ne. 0) then
         iv_spim=imud1
       ENDIF
       MPI_master_only WRITE(iscreenlog,*) 'iv_spim = ',iv_spim
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_nano_pp = ',iv_phyto_nano_pp
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_diat_pp = ',iv_phyto_diat_pp
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_dino_pp = ',iv_phyto_dino_pp
#if defined key_BLOOM_opt2
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_N = ',iv_diss_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_P = ',iv_diss_P
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_Si = ',iv_diss_Si
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_fond_Nitr = ',iv_diss_fond_Nitr
#endif


#ifdef key_BLOOM_insed
       MPI_master_only WRITE(iscreenlog,*) 'iv_ODU = ',iv_ODU
       MPI_master_only WRITE(iscreenlog,*) 'iv_PFe = ',iv_PFe
       MPI_master_only WRITE(iscreenlog,*) 'iv_detrR_N = ',iv_detrR_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_detrR_P = ',iv_detrR_P
#endif

       MPI_master_only WRITE(iscreenlog,*) 'iv_oxygen = ',iv_oxygen

       MPI_master_only WRITE(iscreenlog,*) ' '
       MPI_master_only WRITE(iscreenlog,*) '********************************'
       MPI_master_only WRITE(iscreenlog,*) ' '

#if ! defined MUSTANG && ! defined key_conta
     sumverif = iv_nutr_NH4 + iv_nutr_NO3 + iv_nutr_SiOH + iv_nutr_PO4 + iv_nutr_Pads + &
                iv_phyto_nano_N + iv_phyto_diat_N + iv_phyto_dino_N +                   &
                iv_zoo_meso_N + iv_zoo_micr_N + iv_detr_N + iv_detr_Si + iv_detr_P +    &
#if defined key_BLOOM_opt2
                iv_diss_N+ iv_diss_P + iv_diss_Si + iv_diss_fond_Nitr+                  &
#endif

#if defined key_BLOOM_insed
                iv_ODU + iv_PFe + iv_detrR_N + iv_detrR_P            &
#endif

                iv_oxygen + &

                iv_spim + iv_phyto_nano_pp + iv_phyto_diat_pp + iv_phyto_dino_pp


       IF ( sumverif /= sumindex(ivmin,ivmax) ) THEN
         MPI_master_only WRITE(ierrorlog,*)' '
         MPI_master_only WRITE(ierrorlog,*)'bloom_INIT_IV'
         MPI_master_only WRITE(ierrorlog,*)'the initialization of indexes of substances'
         MPI_master_only WRITE(ierrorlog,*)'is not correct'
         MPI_master_only WRITE(ierrorlog,*)'Rank min=',ivmin,'Rank max=',ivmax,'sumindex(ivmin,ivmax)=',sumindex(ivmin,ivmax)
         MPI_master_only WRITE(ierrorlog,*)'sumverif=',sumverif
         MPI_master_only WRITE(ierrorlog,*)'Have a look at simu.log and variable.dat files'
         MPI_master_only WRITE(ierrorlog,*)'and inquire into the mistakes'
         MPI_master_only WRITE(ierrorlog,*)' '
         MPI_master_only WRITE(ierrorlog,*)'The simulation is stopped'
         STOP
       END IF
#endif
     ENDIF_MPI


     ivfix_cumulprod_first=min(iv_phyto_diat_pp,iv_phyto_dino_pp,iv_phyto_nano_pp)
     ivfix_cumulprod_last=max(iv_phyto_diat_pp,iv_phyto_dino_pp,iv_phyto_nano_pp)

     MPI_master_only WRITE(iscreenlog,*)'*******'
     MPI_master_only WRITE(iscreenlog,*)' first fixed variable production =',ivfix_cumulprod_first
     MPI_master_only WRITE(iscreenlog,*)' last fixed variable production =',ivfix_cumulprod_last
     MPI_master_only WRITE(iscreenlog,*)'*******'

  ENDIF   ! ends test on nv_state

  END SUBROUTINE bloom_init_iv

   !!======================================================================

  SUBROUTINE bloom_init_id(id,standname)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_init_id  ***
   !&E
   !&E ** Purpose : Initialize indexes of diagnostic variables
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : sub_read_vardiag
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2009-10  (V. Garnier)  Original code
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER,          INTENT( in ) :: id
   CHARACTER(LEN=*), INTENT( in ) :: standname

   !! * Local declarations
   INTEGER               :: sumverif

   !!----------------------------------------------------------------------
   !! * Executable part

   SELECT CASE( TRIM(ADJUSTL(ADJUSTR(standname))) )

     CASE('maximum_diatom_mass_concentration_in_sea_water')
       id_diat_max = id
     CASE('date_of_maximum_diatom_mass_concentration_in_sea_water')
       id_diat_datemax = id
     CASE('depth_of_maximum_diatom_mass_concentration_in_sea_water')
       id_diat_depthmax = id
     CASE('maximum_dinoflagellate_mass_concentration_in_sea_water')
       id_dino_max = id
     CASE('date_of_maximum_dinoflagellate_mass_concentration_in_sea_water')
       id_dino_datemax = id
     CASE('depth_of_maximum_dinoflagellate_mass_concentration_in_sea_water')
       id_dino_depthmax = id
     CASE('maximum_nanoflagellate_mass_concentration_in_sea_water')
       id_nano_max = id
     CASE('date_of_maximum_nanoflagellate_mass_concentration_in_sea_water')
       id_nano_datemax = id
     CASE('depth_of_maximum_nanoflagellate_mass_concentration_in_sea_water')
       id_nano_depthmax = id
     CASE('maximum_vertical_gradient_of_sea_water_salinity')
       id_gradsali_max = id
     CASE('depth_of_maximum_vertical_gradient_of_sea_water_salinity')
       id_gradsali_depthmax = id
     CASE('maximum_vertical_gradient_of_sea_water_temperature' )
       id_gradtemp_max = id
     CASE('maximum_vertical_gradient_of_sea_water_density')
       id_graddens_max = id
     CASE('depth_of_maximum_vertical_gradient_of_sea_water_temperature')
       id_gradtemp_depthmax = id
     CASE('depth_of_maximum_vertical_gradient_of_sea_water_density')
       id_graddens_depthmax = id
#ifdef key_BLOOM_insed
     CASE('diffusive_flux_of_ammonium_through_water_sediment_interface')
       id_diffuflux_NH4 = id
     CASE('diffusive_flux_of_nitrate_through_water_sediment_interface')
       id_diffuflux_NO3 = id
     CASE('diffusive_flux_of_phosphate_through_water_sediment_interface')
       id_diffuflux_PO4 = id
     CASE('diffusive_flux_of_oxygen_through_water_sediment_interface')
       id_diffuflux_O2D = id
     CASE('flux_aerobic_miner_Norg')
       id_remin_aerN = id
     CASE('flux_aerobic_miner_Porg')
       id_remin_aerP = id
     CASE('flux_aerobic_miner_biogenic_Si')
       id_remin_aerSi = id
     CASE('flux_anaerobic_miner_Norg')
       id_remin_anaerN = id
     CASE('flux_anaerobic_miner_Porg')
       id_remin_anaerP = id
     CASE('flux_nitrate_miner_Norg')
       id_remin_nitrateN = id
     CASE('flux_dnra_miner_Norg')
       id_remin_dnraN = id
     CASE('flux_denitrification_miner_Norg')
       id_remin_denitN = id
     CASE('flux_nitrate_miner_Porg')
       id_remin_nitrateP = id
     CASE('flux_nitrification')
       id_nitrif = id
     CASE('flux_ODU_oxyd_solid')
       id_oxyd_solid_ODU = id
     CASE('flux_adsorb_desorb_P')
       id_adsor_desorb_P = id
     CASE('flux_dissolution_PFe')
       id_dissol_PFe = id
     CASE('flux_precipitation_P')
       id_precipit_P = id
     CASE('flux_precipitation_Si')
       id_precipit_Si = id
     CASE('cum_flux_mortality_phyto_sed')
       id_morta_phyto = id
     CASE('cum_flux_benthic_filter_grazing')
       id_filtr_benth = id
     CASE('cum_flux_aeration_sediment')
       id_fluxsed_aeration = id
     CASE('sediment_porosity')
        id_porosite_sed= id
#endif
#if defined key_BLOOM_opt2
     CASE('corrected_concentration_of_suspended_matter_in_sea_water')
       id_spim_satused = id
#endif
     CASE('light_limitation_of_diatom_growth')
       id_diat_limlight = id
     CASE('nitrogen_limitation_of_diatom_growth')
       id_diat_limN = id
     CASE('silicon_limitation_of_diatom_growth')
       id_diat_limSi = id
     CASE('phosphate_limitation_of_diatom_growth')
       id_diat_limP = id
     CASE('light_limitation_of_dinoflagellate_growth')
       id_dino_limlight = id
     CASE('nitrogen_limitation_of_dinoflagellate_growth')
       id_dino_limN = id
     CASE('phosphate_limitation_of_dinoflagellate_growth')
       id_dino_limP = id
     CASE('light_limitation_of_nanopicoplankton_growth')
       id_nano_limlight = id
     CASE('nitrogen_limitation_of_nanopicoplankton_growth')
       id_nano_limN = id
     CASE('phosphate_limitation_of_nanopicoplankton_growth')
       id_nano_limP = id
     CASE('light_extinction_in_sea_water')
       id_extinctioncoeff=id
#if defined key_BLOOM_opt2
     CASE('cumulated_production_of_diatoms_in_sea_water_column_expressed_in_nitrogen')
#else
     CASE('cumulated_production_of_diatoms_in_sea_water_column_expressed_in_carbon')
#endif
       id_diat_columnprod=id
#if defined key_BLOOM_opt2
     CASE('cumulated_production_of_dinoflagellates_in_sea_water_column_expressed_in_nitrogen')
#else
     CASE('cumulated_production_of_dinoflagellates_in_sea_water_column_expressed_in_carbon')
#endif
       id_dino_columnprod=id
#if defined key_BLOOM_opt2
     CASE('cumulated_production_of_nanoflagellates_in_sea_water_column_expressed_in_nitrogen')
#else
     CASE('cumulated_production_of_nanoflagellates_in_sea_water_column_expressed_in_carbon')
#endif
       id_nano_columnprod=id
     CASE('chlorophyll_mass_concentration_in_sea_water')
       id_totalchl=id
#if defined key_BLOOM_opt2
     CASE('cumulated_total_production_in_sea_water_column_expressed_in_nitrogen')
#else
     CASE('cumulated_total_production_in_sea_water_column_expressed_in_carbon')
#endif
       id_columnprodtotal=id
     CASE('diatom_sinking_velocity_in_sea_water')
       id_diatsettling=id
     CASE('detritus_sinking_velocity_in_sea_water')
       id_detsettling=id
     CASE('total_suspended_matter_in_sea_water')
       id_spm_total=id
     CASE('txfiltbenth')
        id_benthos_txf=id     
     CASE('oxygen_saturation')
       id_oxy_sat=id  

     CASE DEFAULT
             MPI_master_only WRITE(iscreenlog,*) 'Invalid standard name of diagnostic variable : ', &
                                             TRIM(ADJUSTL(ADJUSTR(standname)))
   END SELECT


   ! smallest and largest indexes
   IF (id >= idmax) idmax=id
   IF (id <= idmin) idmin=id


   IF ( id==ndiag_tot) THEN
   IF_MPI (MASTER) THEN
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) '********************************'
     MPI_master_only WRITE(iscreenlog,*) '          bloom_INIT_ID         '
     MPI_master_only WRITE(iscreenlog,*) ' Initialization of id_variables '
     MPI_master_only WRITE(iscreenlog,*) '********************************'
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_max = ',id_diat_max
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_datemax = ',id_diat_datemax
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_depthmax = ',id_diat_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_max = ',id_dino_max
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_datemax = ',id_dino_datemax
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_depthmax = ',id_dino_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_max = ',id_nano_max
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_datemax = ',id_nano_datemax
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_depthmax = ',id_nano_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_gradsali_max = ',id_gradsali_max
     MPI_master_only WRITE(iscreenlog,*) 'id_gradsali_depthmax = ',id_gradsali_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_gradtemp_max = ',id_gradtemp_max
     MPI_master_only WRITE(iscreenlog,*) 'id_gradtemp_depthmax = ',id_gradtemp_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_graddens_max = ',id_graddens_max
     MPI_master_only WRITE(iscreenlog,*) 'id_graddens_depthmax = ',id_graddens_depthmax
#if defined key_BLOOM_opt2
     MPI_master_only WRITE(iscreenlog,*) 'id_spim_satused = ',id_spim_satused
#endif
#ifdef key_BLOOM_insed
     MPI_master_only WRITE(iscreenlog,*) 'id_diffuflux_NO3 = ',id_diffuflux_NO3
     MPI_master_only WRITE(iscreenlog,*) 'id_diffuflux_NH4 = ',id_diffuflux_NH4
     MPI_master_only WRITE(iscreenlog,*) 'id_diffuflux_PO4 = ',id_diffuflux_PO4
     MPI_master_only WRITE(iscreenlog,*) 'id_diffuflux_O2D = ',id_diffuflux_O2D
#endif
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_limlight = ',id_diat_limlight
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_limN = ',id_diat_limN
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_limSi = ',id_diat_limSi
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_limP = ',id_diat_limP
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_limlight = ',id_dino_limlight
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_limN = ',id_dino_limN
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_limP = ',id_dino_limP
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_limlight = ',id_nano_limlight
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_limN = ',id_nano_limN
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_limP = ',id_nano_limP
     MPI_master_only WRITE(iscreenlog,*) 'id_extinctioncoeff = ',id_extinctioncoeff
     MPI_master_only WRITE(iscreenlog,*) 'id_diat_columnprod = ',id_diat_columnprod
     MPI_master_only WRITE(iscreenlog,*) 'id_dino_columnprod = ',id_dino_columnprod
     MPI_master_only WRITE(iscreenlog,*) 'id_nano_columnprod = ',id_nano_columnprod
     MPI_master_only WRITE(iscreenlog,*) 'id_totalchl = ',id_totalchl
     MPI_master_only WRITE(iscreenlog,*) 'id_columnprodtotal = ',id_columnprodtotal
     MPI_master_only WRITE(iscreenlog,*) 'id_diatsettling = ',id_diatsettling
     MPI_master_only WRITE(iscreenlog,*) 'id_detsettling = ',id_detsettling
     MPI_master_only WRITE(iscreenlog,*) 'id_spm_total = ',id_spm_total
     MPI_master_only WRITE(iscreenlog,*) 'id_benthos_txf=' ,id_benthos_txf

    MPI_master_only WRITE(iscreenlog,*) 'id_oxy_sat=',id_oxy_sat

#ifdef key_BLOOM_insed
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_aerN=',id_remin_aerN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_aerP=',id_remin_aerP
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_aerSi=',id_remin_aerSi
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_anaerN=',id_remin_anaerN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_anaerP=',id_remin_anaerP
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_nitrateN=',id_remin_nitrateN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_dnraN=',id_remin_dnraN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_denitN=',id_remin_denitN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_nitrateP=',id_remin_nitrateP
    MPI_master_only WRITE(iscreenlog,*) 'id_dissol_PFe=',id_dissol_PFe
    MPI_master_only WRITE(iscreenlog,*) 'id_nitrif =',id_nitrif
    MPI_master_only WRITE(iscreenlog,*) 'id_oxyd_solid_ODU=',id_oxyd_solid_ODU
    MPI_master_only WRITE(iscreenlog,*) 'id_adsor_desorb_P=',id_adsor_desorb_P
    MPI_master_only WRITE(iscreenlog,*) 'id_precipit_P=',id_precipit_P
    MPI_master_only WRITE(iscreenlog,*) 'id_precipit_Si=',id_precipit_Si
    MPI_master_only WRITE(iscreenlog,*) 'id_morta_phyto=',id_morta_phyto
    MPI_master_only WRITE(iscreenlog,*) 'id_filtr_benth =',id_filtr_benth
    MPI_master_only WRITE(iscreenlog,*) 'id_fluxsed_aeration =',id_fluxsed_aeration
    MPI_master_only WRITE(iscreenlog,*) 'id_porosite_sed =',id_porosite_sed
#endif


     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) '********************************'
     MPI_master_only WRITE(iscreenlog,*) ' '

 
#if ! defined MUSTANG && ! defined key_conta
     sumverif = id_diat_max + id_diat_datemax + id_diat_depthmax + &
              id_dino_max + id_dino_datemax + id_dino_depthmax + &
              id_nano_max + id_nano_datemax + id_nano_depthmax + &
              id_gradsali_max + id_gradsali_depthmax +           &
              id_gradtemp_max + id_gradtemp_depthmax  + &
              id_graddens_max + id_graddens_depthmax  + &
#if defined key_BLOOM_opt2
              id_spim_satused + &
#endif
              id_diat_limlight + id_diat_limN + id_diat_limSi + id_diat_limP + &
              id_dino_limlight + id_dino_limN + id_dino_limP +   &
              id_nano_limlight + id_nano_limN + id_nano_limP +   &
              id_extinctioncoeff + id_diat_columnprod +          &
              id_dino_columnprod + id_nano_columnprod +          &
              id_benthos_txf + &

              id_oxy_sat  + &              
              id_totalchl + id_columnprodtotal + id_spm_total +  &
              id_diatsettling + id_detsettling

     IF ( sumverif /= sumindex(idmin,idmax) ) THEN
       write(ierrorlog,*)' '
       write(ierrorlog,*)'bloom_INIT_ID'
       write(ierrorlog,*)'the initialization of indexes of diagnostic variables'
       write(ierrorlog,*)'is not correct'
       write(ierrorlog,*)'Rank min=',idmin,'Rank max=',idmax,'sumindex(idmin,idmax)=',sumindex(idmin,idmax)
       write(ierrorlog,*)'sumverif=',sumverif
       write(ierrorlog,*)'Have a look at simu.log and vardiag.dat files'
       write(ierrorlog,*)'and inquire into the mistakes'
       write(ierrorlog,*)' '
       write(ierrorlog,*)'The simulation is stopped'
       STOP
     END IF
#endif
   ENDIF_MPI
   END IF   ! ends test on ndiag_tot

  END SUBROUTINE bloom_init_id



   !!======================================================================
  FUNCTION sumindex(nstart,nend)

   !&E---------------------------------------------------------------------
   !&E                 ***  FUNCTION sumindex  ***
   !&E
   !&E ** Purpose : sum of indexes associated to biological variables
   !&E              to make sure all indexes are initialized
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
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
   !&E       !  2009-10  (V. Garnier)  Original code
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT( in ) :: nstart,nend          ! index of the begining/end of the series
   INTEGER               :: sumindex             ! result of the function

   !! * Local declarations
   INTEGER               :: i                    ! loop index

   !!----------------------------------------------------------------------
   !! * Executable part

   sumindex = nstart
   i = nstart + 1

   DO WHILE (i <= nend)
     sumindex = sumindex + i
     i = i + 1
   END DO

  END FUNCTION sumindex

   !!======================================================================

#endif

END MODULE
