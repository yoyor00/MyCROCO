
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

#ifdef MUSTANG
  USE comMUSTANG, ONLY : valmanq
# endif

  IMPLICIT NONE

  !! * Accessibility
  PUBLIC bloom_param                   ! routine called by init.F90
  PUBLIC bloom_init_iv,bloom_init_id   ! routine called by subreaddat.F90 or initBIOLink.F90
#if defined key_N_tracer || defined key_P_tracer
  PUBLIC bloom_create_vartracer,bloom_create_vardiagtracer   ! routine called by subreaddat.F90
#endif
  PUBLIC bloom_userinit                
 
   !! * Interface
  INTERFACE bloom_userinit
!     MODULE PROCEDURE bloom_init_dinodiatzoo
     MODULE PROCEDURE bloom_init_nut4phy3zoo2
  END INTERFACE

#if defined key_N_tracer || defined key_P_tracer
  INTERFACE bloom_create_vartracer
#ifdef key_MARS
     MODULE PROCEDURE bloom_create_vartracer_MARS
#else
     MODULE PROCEDURE bloom_create_vartracer_other
#endif
  END INTERFACE
#endif

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
#ifdef key_CROCO
   INTEGER           :: lstr,lenstr
#endif
#ifdef key_N_tracer
   CHARACTER(LEN=19) :: date_start_tracerN
   REAL(KIND=rlg)        :: tool_datosec
#elif defined key_P_tracer
   CHARACTER(LEN=19) :: date_start_tracerP
   REAL(KIND=rlg)        :: tool_datosec
#endif

   NAMELIST/namBIOLink/l_bioretro_extinct,dt_bio_update,l_waterdensity_known,  &
                     i_BIOLink_verif,j_BIOLink_verif,DT_CONSERV_BIOLINK
   NAMELIST/namoptions/ l_filtbenthsinus,l_filtbenthmes,l_SNeffect_settle,l_phyzoodeteffect_settle,  &
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
#if defined GAMELAG
                       p_kO2_nit_wat, p_nitrif_wat,p_labi,p_K_lP_max,                   &             
                       p_K_rP_max,p_Temp_opt_OM,p_Temp_width_OM,p_Vmax_lD,p_Vmax_rD,    &
                       p_K_l_N,p_K_r_N,p_K_l_P,p_K_r_P,p_gamma,p_epsilon,p_K_min,       &
#endif
                       p_T_effectSi,p_k_remin,p_xflimz,p_kSi,p_Si_precip
   NAMELIST/namphosphor/ p_P_speedup_reminanaer,p_P_adsor,p_P_desor,p_P_adsormaxspim, &
                         p_P_adsormaxsed,p_P_precFeO2,p_kO2_precPFe,p_kNO3_precPFe,p_kiO2_dissPFe,p_kiO2_desorP,  &
                         p_kiNO3_dissPFe,p_P_precFeNO3,p_P_dissFe,p_GO2_PFe,p_GNO3_PFe
#endif
   NAMELIST/namgenphyto/ p_phyto_ChlNratio,p_phyto_SiNratio,p_phyto_NPratio,p_phyto_CNratio, &
                         p_phyto_ChlNratiomax,p_phyto_ChlN_ksmithextinct
   NAMELIST/namnanophyto/ p_nano_mumax,p_nano_kNO3,p_nano_kNH4,p_nano_kPO4,p_nano_mort, &
#if defined GAMELAG
                          p_nano_phy_NH4,p_nano_excret,p_nano_aO2,p_nano_DOcrit, &
#endif
                          p_nano_iksmith,p_nano_thhold_mort
   NAMELIST/namdiatom/ p_diat_mumax,p_diat_kNO3,p_diat_kNH4,p_diat_kSi,p_diat_kPO4,     &
#if defined GAMELAG
                       p_diat_excret,p_diat_aO2,p_diat_DOcrit, &
#endif
                       p_diat_iksmith,p_diat_mort,p_diat_thhold_mort
   NAMELIST/namdino/ p_dino_mumax,p_dino_thhold_ect,p_dino_kNO3,p_dino_kNH4,p_dino_kPO4,p_dino_mort,      &
                     p_dino_iksmith,p_dino_thhold_mort
#if defined key_psnz
   NAMELIST/nampsnz/ p_psnz_mumax,p_psnz_tempopt,p_psnz_kNO3,p_psnz_kNH4,p_psnz_kSi,p_psnz_kPO4,     &
                     p_psnz_vmaxSi,p_psnz_qminSi,p_psnz_kqSi,p_psnz_iksmith,p_psnz_mort,p_mesz_captpsnz,  &
                     p_psnz_thhold_Si,p_psnz_prod_domoic,p_psnz_decay_domoic,p_psnz_templethal,p_psnz_beta
#endif
#if defined key_karenia
   NAMELIST/namkarenia/ p_karenia_mumax,p_karenia_thhold_ect,p_karenia_tempopt,p_karenia_vmaxN,p_karenia_kqN,    &
                        p_karenia_qminN,p_karenia_kNO3,p_karenia_kNH4,p_karenia_vmaxP,p_karenia_kqP,             &
                        p_karenia_qminP,p_karenia_kPO4,p_karenia_iksmith,p_karenia_mort,p_karenia_wmax,  &
                        p_karenia_templethal,p_karenia_beta,p_mesz_captkarenia,p_micz_captkarenia
#endif
#if defined key_phaeocystis
   NAMELIST/namphaeocystis/ p_phaeocystis_mumaxcolo,p_phaeocystis_mumaxcell,p_phaeocystis_tempopt,    &
                            p_phaeocystis_kNO3colo,p_phaeocystis_kNO3cell,    &
                            p_phaeocystis_kNH4,p_phaeocystis_kPO4colo,p_phaeocystis_kPO4cell,  &
                            p_phaeocystis_iksmith,p_phaeocystis_mort,p_phaeocystis_lyse,p_phaeocystis_coloinit,  &
                            p_phaeocystis_coloseuil,p_phaeocystis_mucus_decay,p_mesz_captphaeocolo,   &
                            p_phaeocystis_templethal, p_phaeocystis_beta
#endif
   NAMELIST/nammesozoo/ p_mesz_thrN,p_mesz_kivlev,p_mesz_mumax,p_mesz_assim,      &
                        p_mesz_excret,p_mesz_mort1,p_mesz_mort2,p_zoo_CDWratio,    &
#if defined GAMELAG
                        p_mesz_reg,p_zoo_NPratio,p_mesz_DOcrit,p_mesz_aO2,           &
#endif                        
                        p_zoo_CNratio,p_mesz_thhold_mes_kivlev,p_mesz_thhold_mort
   NAMELIST/nammicrozoo/ p_micz_mumax,p_micz_kgraz,p_micz_assim,p_micz_excret,p_micz_mort,   &
#if defined GAMELAG
                         p_micz_reg,p_micz_DOcrit,p_micz_aO2,  &
#endif   
                         p_micz_thrnano,p_micz_thhold_mort

#ifdef key_BLOOM_opt2
   NAMELIST/namoxygen/ p_phyto_photoratio,p_kO2_reminO2,p_nitrif,p_phyto_resp,p_zoo_resp
#else
   NAMELIST/namoxygen/ p_phyto_photoratio,p_phyto_resp,p_zoo_resp,                              &
#if defined GAMELAG
                       p_O2_Threshold,p_O2SED_Threshold,p_R_photo,p_Q_photo,                                      &
#endif
                       p_KO2sed_aeration,p_Kzsed_aeration
#endif
   NAMELIST/namoptics/ p_extincwat,p_extincspim,p_extincChl1,p_extincChl2,p_parradratio
#if defined key_benthos
   NAMELIST/namsedbenthos/ p_erodflux,p_erodvitcrit,p_depovitcrit,p_detzoo_wsed,p_detphy_wsed,p_erodcrittau,  &
                           p_burial,p_reminbenth
#elif defined MUSTANG
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
#if defined GAMELAG
                        p_mesz_captnano,                                                      &
#endif                        
                        p_micz_captnano,p_micz_captdet,p_micz_captdino,p_txfiltbenthmax
#ifdef key_ulvas
   NAMELIST/namulva/ p_ulv_mumax,p_ulv_iksmith,p_ulv_kNO3,p_ulv_kPO4,p_ulv_maxabsN,   &
                     p_ulv_maxabsP,p_ulv_minNWratio,p_ulv_maxNWratio,p_ulv_minPWratio,  &
                     p_ulv_maxPWratio,p_ulv_dryfreshratio,p_ulv_susmort,p_ulv_depmort,   &
                     p_ulv_extincwat,p_ulv_erodflux,p_ulv_erodcritv,p_ulv_settlcrittv
#endif

#ifdef key_zostera
  NAMELIST/namzostera/ p_zost_klai,p_zost_leafabscoef,p_zost_Ik,p_zost_maxprod0,p_zost_T_prod, &
                        p_zost_respf0,p_zost_respr0,p_zost_T_respf,p_zost_T_respr, p_zost_qphotos, &
                        p_zost_LNquotamin,p_zost_LNquotamax,p_zost_RNquotamin,p_zost_RNquotamax, &
                        p_zost_KNH4L,p_zost_KNH4R,p_zost_KNO3L,p_zost_delta1,p_zost_delta2, &
                        p_zost_LNnh4maxabs, p_zost_LNno3maxabs, &
                        p_zost_RNmaxabs,p_zost_LPquotamin,p_zost_LPquotamax,p_zost_RPquotamin, &
                        p_zost_RPquotamax,p_zost_KPL,p_zost_KPR,p_zost_LPmaxabs,p_zost_RPmaxabs, &
                        p_zost_K,p_zost_transfrate,p_zost_reclamax,p_zost_mort,p_zost_T_mort,p_zost_RECRmax, &
                        p_zost_T_recr,p_zost_KREC1,p_zost_KREC2,p_zost_SB1,p_zost_KGERM,p_zost_SB0, &
                        p_zost_ERSmax,p_zost_seedprod,p_zost_GERmax,p_zost_Smort
#endif
#ifdef key_oyster_benthos
   NAMELIST/namoysterfiltbent/ p_huitre_DW,p_huitre_thr_colmat,p_huitre_temp1,p_huitre_temp_opt,   &
                               p_huitre_thr_mes,p_huitre_loi_mes1,p_huitre_loi_mes2,p_huitre_filt_std, &
                               p_huitre_exp_allom,p_huitre_loi_colmat
#endif
#ifdef key_oyster_SFG
   NAMELIST/namoysterSFG/ p_oys_absdet,p_oys_absmaxmop,p_oys_kreten,p_oys_allfiltr,p_oys_allres,p_oys_yreten,p_oys_CtoDWratio, &
                       p_oys_stdfiltr,p_oys_paramfiltr,p_oys_tempopt_filtr,p_oys_sestfiltr_thr,p_oys_yfilt,     &
                       p_oys_kfilt,p_oys_colmparam,p_oys_pseufeci_exp,p_oys_pseufeco_exp,p_oys_prodpseufeci_lev,   &
                       p_oys_prodpseufeco_lev,p_oys_retmin,p_oys_retmin_thr,p_oys_prodpseufec_thr,        &
                       p_oys_subratio,p_oys_colmat_thr
#endif
#ifdef key_oyster_DEB
   NAMELIST/namoysterdeb/ p_oysDEB_TA,p_oysDEB_TAL,p_oysDEB_TAH,p_oysDEB_TL, &
                          p_oysDEB_THing,p_oysDEB_THresp,p_oysDEB_Tseuilponte, &
                          p_oysDEB_Eg,p_oysDEB_Em,p_oysDEB_Kappa,p_oysDEB_Vp, &
                          p_oysDEB_kR,p_oysDEB_shape,p_oysDEB_pxm0,p_oysDEB_pammax, &
                          p_oysDEB_ae,p_oysDEB_pm0,p_oysDEB_ERlim,p_oysDEB_kchl, &
                          p_oysDEB_muE
#endif
#ifdef key_oyster_DEB_GAMELAG
   NAMELIST/namoysterdeb_gamelag/Ta,Tl,Th,Tal,Tah,T1 , &
                                 Jxm,Pm_deb,XkN_SPAT,XkP_SPAT,XkN_OYST,XkP_OYST, &
                                 KappaX,Em,Eg,Kappa,deltam,Lp,T_spawn, GSR_spawn, &
                                 muE, d_v, T_im, dgo,Eggo,Ygo, &
                                 Yv, KappaGo, muE_N,muV,muGo,alpha_PS,alpha_PL,alpha_ZS,alpha_ZL, &
                                 alpha_PON,alpha_POP,NP_oyster,C_oyster, &
                                 BioDpo,DOcrit,DOcrit_MR,bDO,OMR_DO,OMR_Dpo, &
                                 epsOyst,oyster_mortality,DOcrit
#endif
#if defined key_microtracers
   NAMELIST/nammicrotrace/ p_trace_debitinject,p_trace_depth1inject,p_trace_depth2inject,p_trace_depth3inject
#endif
#if defined key_N_tracer
   NAMELIST/namtracerN/ nb_source_river_tracerN,nb_source_marin_tracerN,date_start_tracerN,   &
                        p_source_river1_tracerN,p_source_river2_tracerN,p_source_river3_tracerN,  &
                        p_source_river4_tracerN,p_marqueNH4,p_marqueNO3,p_marquenanoN,p_marquedinoN,p_marquediatN,  &
                        p_marquemicrN,p_marquemesoN,p_marquedetN
#endif
#if defined key_P_tracer
   NAMELIST/namtracerP/ nb_source_river_tracerP,nb_source_marin_tracerP,date_start_tracerP,   &
                        p_source_river1_tracerP,p_source_river2_tracerP,p_source_river3_tracerP,  &
                        p_source_river4_tracerP,p_marquePO4,p_marquePads,p_marquedetP,p_marquenanoP,p_marquedinoP,  &
                        p_marquediatP,p_marquemicrP,p_marquemesoP
#endif
   !!----------------------------------------------------------------------
#ifdef key_benthos_gener
   NAMELIST/nammeio/    p_meio_satur_proie,p_meio_capdet,p_meio_thrdet,p_meio_capdiat,p_meio_thrdiat,  &
                        p_meio_mort,p_meio_pphysio,p_meio_egest_det,p_meio_thrmax_inhib,  &
                        p_meio_thrmin_inhib,p_meio_inges
   NAMELIST/namdeposi/  p_deposi_satur_proie,p_deposi_capdet,p_deposi_thrdet,p_deposi_capdiat,  &
                        p_deposi_capmeio,p_deposi_thrmeio,p_deposi_mort,p_deposi_pphysio,  &
                        p_deposi_egest_meio,p_deposi_egest_det,p_deposi_egest_diat,  &
                        p_deposi_thrmax_inhib,p_deposi_thrmin_inhib,p_deposi_inges
   NAMELIST/namsuspens/ p_susp_satur_proie,p_susp_capdet,p_susp_thrdet,p_susp_capphy,p_susp_thrphy,      &
                        p_susp_mort,p_susp_pphysio,p_susp_egest_det,p_susp_egest_phy,p_susp_thrmax_inhib,   &
                        p_susp_thrmin_inhib,p_susp_inges
   NAMELIST/namcarn/    p_carn_satur_proie,p_carn_capmeio,p_carn_thrmeio,p_carn_capdeposi,  &
                        p_carn_thrdeposi,p_carn_satur_proie_susp,p_carn_capsusp,p_carn_thrsusp, &
                        p_carn_mort,p_carn_pphysio,p_carn_egest,p_carn_thrmax_inhib,p_carn_thrmin_inhib, &
                        p_carn_inges
#endif



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
#if defined key_psnz
   READ(50,nampsnz)
#endif
#if defined key_karenia
   READ(50,namkarenia)
#endif
#if defined key_phaeocystis
   READ(50,namphaeocystis)
#endif
   READ(50,nammesozoo)
   READ(50,nammicrozoo)
   READ(50,namoxygen)
   READ(50,namoptics)
#if defined key_benthos
   READ(50,namsedbenthos)
#elif defined MUSTANG
   READ(50,namsediment)
#ifdef key_sedim_MPB
   READ(50,namMPB)
   invMPBCmaxm3=12.e-3/(45._rsh*Chlamax_m2microm) *1.e-3  ! divise par 1000 passer de mole en mmole
   MPBCmaxm3=1._rsh/invMPBCmaxm3
   QR_diat_inv=1.0_rsh/QR_diat

#endif
#endif
   READ(50,namgrazing)
#ifdef key_ulvas
   READ(50,namulva)
#endif
#ifdef key_zostera
   READ(50,namzostera)
#endif
#ifdef key_oyster_benthos
   READ(50,namoysterfiltbent)
#endif
#ifdef key_oyster_SFG
   READ(50,namoysterSFG)
#endif
#ifdef key_oyster_DEB
   READ(50,namoysterDEB)
#endif
#ifdef key_oyster_DEB_GAMELAG
   READ(50,namoysterDEB_GAMELAG)
#endif
#if defined key_microtracers
   READ(50,nammicrotrace)
#endif
#if defined key_N_tracer
   READ(50,namtracerN)
   tdeb_tracerN=tool_datosec(date_start_tracerN)
   nb_source_tracerN=nb_source_river_tracerN+nb_source_marin_tracerN
   ALLOCATE(name_source_tracerN(nb_source_tracerN))
   DO iso=1,nb_source_river_tracerN
      SELECT CASE(iso)
            CASE (1)
             name_source_tracerN(iso)=p_source_river1_tracerN
            CASE (2)
             name_source_tracerN(iso)=p_source_river2_tracerN
            CASE (3)
             name_source_tracerN(iso)=p_source_river3_tracerN
            CASE (4)
             name_source_tracerN(iso)=p_source_river4_tracerN
            CASE DEFAULT
               MPI_master_only  WRITE(ierrorlog,*) ' '
               MPI_master_only  WRITE(ierrorlog,*) 'ERROR: NOT MORE THAN 4 RIVER SOURCES for TRACER N'
               MPI_master_only  WRITE(ierrorlog,*) '       CHECK in parabloom.txt : nb_source_river_tracerN must be < 5' 
               MPI_master_only  WRITE(ierrorlog,*) '  Simulation stopped'
               CALL_MPI MPI_FINALIZE(IERR_MPI)
               STOP
      END SELECT
   ENDDO
   DO iso=nb_source_river_tracerN+1,nb_source_tracerN
     name_source_tracerN(iso)='MARIN'
   ENDDO
#endif
#if defined key_P_tracer
   READ(50,namtracerP)
   tdeb_tracerP=tool_datosec(date_start_tracerP)
   nb_source_tracerP=nb_source_river_tracerP+nb_source_marin_tracerP
   ALLOCATE(name_source_tracerP(nb_source_tracerP))
   DO iso=1,nb_source_river_tracerP
      SELECT CASE(iso)
            CASE (1)
             name_source_tracerP(iso)=p_source_river1_tracerP
            CASE (2)
             name_source_tracerP(iso)=p_source_river2_tracerP
            CASE (3)
             name_source_tracerP(iso)=p_source_river3_tracerP
            CASE (4)
             name_source_tracerP(iso)=p_source_river4_tracerP
            CASE DEFAULT
               MPI_master_only  WRITE(ierrorlog,*) ' '
               MPI_master_only  WRITE(ierrorlog,*) 'ERROR: NOT MORE THAN 4 RIVER SOURCES for TRACER N'
               MPI_master_only  WRITE(ierrorlog,*) '       CHECK in parabloom.txt : nb_source_river_tracerN must be < 5' 
               MPI_master_only  WRITE(ierrorlog,*) '  Simulation stopped'
               CALL_MPI MPI_FINALIZE(IERR_MPI)
               STOP
      END SELECT
   ENDDO
   DO iso=nb_source_river_tracerP+1,nb_source_tracerP
     name_source_tracerP(iso)='MARIN'
   ENDDO
   IF(nb_source_tracerP /= nb_sourceP)THEN
      MPI_master_only  WRITE(ierrorlog,*)'ERROR : total number source of tracer in parabloom.txt is different from nb_source in parameters'
      MPI_master_only  WRITE(ierrorlog,*)'The simulation is stopped'
      CALL_MPI MPI_FINALIZE(IERR_MPI)
      STOP
   ENDIF
#endif
#ifdef key_benthos_gener
   READ(50,nammeio)
   READ(50,namdeposi)
   READ(50,namsuspens)
   READ(50,namcarn)
#endif

  CLOSE(50)                                  ! end reading of namelists

ELSE

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
#if defined key_benthos
     MPI_master_only WRITE(iscreenlog,namsedbenthos)
#elif defined MUSTANG
     MPI_master_only WRITE(iscreenlog,namsediment)
#ifdef key_sedim_MPB
     MPI_master_only WRITE(iscreenlog,namMPB)
#endif
#endif
     MPI_master_only WRITE(iscreenlog,namgrazing)

#if defined key_microtracers
     MPI_master_only  WRITE(iscreenlog,nammicrotrace)
#endif
#ifdef key_ulvas
     MPI_master_only  WRITE(iscreenlog,namulva)
#endif
#ifdef key_oyster_benthos
     MPI_master_only  WRITE(iscreenlog,namoysterfiltbent)
#endif
#ifdef key_oyster_SFG
     MPI_master_only  WRITE(iscreenlog,namoysterSFG)
#endif
#ifdef key_oyster_DEB
     MPI_master_only  WRITE(iscreenlog,namoysterDEB)
#endif
#ifdef key_oyster_DEB_GAMELAG
     MPI_master_only  WRITE(iscreenlog,namoysterDEB_GAMELAG)
#endif
#if defined key_psnz
     MPI_master_only  WRITE(iscreenlog,nampsnz)
#endif	     
#if defined key_karenia
     MPI_master_only  WRITE(iscreenlog,namkarenia)
#endif	
#if defined key_phaeocystis
     MPI_master_only  WRITE(iscreenlog,namphaeocystis)
#endif     
#if defined key_N_tracer
     MPI_master_only  WRITE(iscreenlog,namtracerN)
#endif
#if defined key_P_tracer
     MPI_master_only  WRITE(iscreenlog,namtracerP)
#endif
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
#if defined key_N_tracer || defined key_P_tracer
   INTEGER               :: is,ivtra
#endif

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
#if defined GAMELAG
     CASE('mole_concentration_of_dissolved_organic_detritus_expressed_as_nitrogen_in_sea_water')
       iv_diss_detr_N = irk_mod(iv)
     CASE('mole_concentration_of_dissolved_refractory_organic_detritus_expressed_as_nitrogen_in_sea_water')
       iv_diss_detrR_N = irk_mod(iv)
     CASE('mole_concentration_of_dissolved_organic_detritus_expressed_as_phosphorus_in_sea_water')
       iv_diss_detr_P = irk_mod(iv)
     CASE('mole_concentration_of_dissolved_refractory_organic_detritus_expressed_as_phosphorus_in_sea_water')
       iv_diss_detrR_P = irk_mod(iv)       
#endif
#if defined key_psnz
     CASE('mole_concentration_of_pseudonitzschia_expressed_as_nitrogen_in_sea_water')
       iv_phyto_psnz_N = irk_mod(iv)
     CASE('mole_concentration_of_pseudonitzschia_expressed_as_silicon_in_sea_water')
       iv_phyto_psnz_Si = irk_mod(iv)
     CASE('mass_concentration_of_domoic_acid_in_sea_water')
       iv_phyto_psnz_ad = irk_mod(iv)
     CASE('cumulative_pseudonitzschia_production_expressed_as_carbon_in_sea_water')
       iv_phyto_psnz_pp = irk_mod(iv)
#endif
#if defined key_karenia
     CASE('mole_concentration_of_karenia_expressed_as_carbon_in_sea_water')
       iv_phyto_karenia_C = irk_mod(iv)
     CASE('mole_concentration_of_karenia_expressed_as_nitrogen_in_sea_water')
       iv_phyto_karenia_N = irk_mod(iv)
     CASE('mole_concentration_of_karenia_expressed_as_phosphorus_in_sea_water')
       iv_phyto_karenia_P = irk_mod(iv)
     CASE('cumulative_karenia_production_expressed_as_carbon_in_sea_water')
       iv_phyto_karenia_pp = irk_mod(iv)
#endif
#if defined key_phaeocystis
     CASE('mole_concentration_of_colonial_phaeocystis_expressed_as_nitrogen_in_sea_water')
       iv_phyto_phaeocystis_colo_N = irk_mod(iv)
     CASE('mole_concentration_of_phaeocystis_cells_expressed_as_nitrogen_in_sea_water')
       iv_phyto_phaeocystis_cell_N = irk_mod(iv)
     CASE('concentration_of_phaeocystis_mucus_expressed_as_mass_in_sea_water')
       iv_phyto_phaeocystis_mucus = irk_mod(iv)
     CASE('cumulative_phaeocystis_production_expressed_as_carbon_in_sea_water')
       iv_phyto_phaeocystis_pp = irk_mod(iv)
#endif
#ifdef key_zoo_prod
     CASE('cumulative_microzooplankton_production_expressed_as_carbon_in_sea_water')
       iv_zoo_micr_ps = irk_mod(iv)
     CASE('cumulative_mesozooplankton_production_expressed_as_carbon_in_sea_water')
       iv_zoo_meso_ps = irk_mod(iv)
#endif
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
#ifdef key_ulvas
     CASE('suspended_ulvas_nitrogen')
       iv_ulv_sus_N = irk_mod(iv)
     CASE('benthic_ulvas_nitrogen')
       iv_ulv_benth_N = irk_mod(iv)
     CASE('suspended_ulvas_dry_weight')
       iv_ulv_susdrywght = irk_mod(iv)
     CASE('benthic_ulvas_dry_weight')
       iv_ulv_benthdrywght = irk_mod(iv)
     CASE('suspended_ulvas_phosphorus')
       iv_ulv_sus_P = irk_mod(iv)
     CASE('benthic_ulvas_phosphorus')
       iv_ulv_benth_P = irk_mod(iv) 
     CASE('cumulative_suspended_ulvas_production')
       iv_ulv_suspp = irk_mod(iv)
     CASE('cumulative_benthic_ulvas_production')
       iv_ulv_benthpp = irk_mod(iv) 
#endif
#ifdef key_zostera
    CASE('zostera_leaf_biomass_per_m2')
       iv_zost_LB = irk_mod(iv)
     CASE('zostera_root_biomass_per_m2')
       iv_zost_RB = irk_mod(iv)
     CASE('zostera_shoot_density_per_m2')
       iv_zost_D = irk_mod(iv)
     CASE('zostera_leaf_nitrogen_pool_per_m2')
       iv_zost_LN = irk_mod(iv)
     CASE('zostera_root_nitrogen_pool_per_m2')
       iv_zost_RN = irk_mod(iv)
     CASE('zostera_leaf_phosphorus_pool_per_m2')
       iv_zost_LP = irk_mod(iv)
     CASE('zostera_root_phosphorus_pool_per_m2')
       iv_zost_RP = irk_mod(iv) 
     CASE('mole_concentration_of_zostera_floating_detritus_expressed_as_nitrogen')
       iv_detr_zost_N = irk_mod(iv)
     CASE('mole_concentration_of_zostera_floating_detritus_expressed_as_phosphorus')
       iv_detr_zost_P = irk_mod(iv)         
     CASE('mole_concentration_of_zostera_detritus_in_benthos_expressed_as_nitrogen')
       iv_zost_benth_N = irk_mod(iv)
     CASE('mole_concentration_of_zostera_detritus_in_benthos_expressed_as_phosphorus')
       iv_zost_benth_P = irk_mod(iv)         
     CASE('mole_concentration_of_zostera_seeds_in_water_column')
       iv_zost_seed = irk_mod(iv)
     CASE('mole_concentration_of_zostera_seeds_in_benthos')
       iv_zost_benth_seed = irk_mod(iv)         
     CASE('mole_concentration_of_ammonium_in_benthos')
       iv_benth_NH4 = irk_mod(iv)         
     CASE('mole_concentration_of_dissolved_phosphorus_in_benthos')
       iv_benth_PO4 = irk_mod(iv) 
     CASE('mole_concentration_of_benthic_diatoms_expressed_as_nitrogen_in_benthos')
       iv_benth_phyto_diat_N = irk_mod(iv)        
!     CASE('cumulative_zostera_production')
!       iv_zost_pp = irk_mod(iv)         
!     CASE('cumulative_zostera_leaf_N_uptake')
!       iv_zost_LNuptake=irk_mod(iv)
!     CASE('cumulative_zostera_leaf_P_uptake')
!       iv_zost_LPuptake=irk_mod(iv)
!     CASE('cumulative_zostera_root_N_uptake')
!       iv_zost_RNuptake=irk_mod(iv)
!     CASE('cumulative_zostera_root_P_uptake')
!       iv_zost_RPuptake=irk_mod(iv)
#endif 
#ifdef key_oyster_SFG
     CASE('oyster_dry_soma_weight')
       iv_oys_so = irk_mod(iv)
     CASE('oyster_dry_resgon_weight')
       iv_oys_go = irk_mod(iv)
     CASE('oyster_coq_weight')
       iv_oys_co = irk_mod(iv)
#endif      
#ifdef key_oxygen
     CASE('dissolved_oxygen_in_water_column')
       iv_oxygen = irk_mod(iv)
#endif
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
#ifdef key_oyster_DEB
     CASE('oysdeb_res')
       iv_oysdeb_res = irk_mod(iv)
     CASE('oysdeb_str')
       iv_oysdeb_str = irk_mod(iv)
     CASE('oysdeb_gon')
       iv_oysdeb_gon = irk_mod(iv)
     CASE('oysdeb_res2')
       iv_oysdeb_res2 = irk_mod(iv)
     CASE('oysdeb_str2')
       iv_oysdeb_str2 = irk_mod(iv)
     CASE('oysdeb_gon2')
       iv_oysdeb_gon2 = irk_mod(iv)
!     CASE('oysdeb_res3')
!       iv_oysdeb_res3 = irk_mod(iv)
!     CASE('oysdeb_str3')
!       iv_oysdeb_str3 = irk_mod(iv)
!     CASE('oysdeb_gon3')
!       iv_oysdeb_gon3 = irk_mod(iv)
#endif
#ifdef key_oyster_DEB_GAMELAG
     CASE('oysdeb')
       iv_oysdeb = irk_mod(iv)
     CASE('oysdeb_E_V')
       iv_oysdeb_E_V = irk_mod(iv)
     CASE('oysdeb_E_GO')
       iv_oysdeb_E_GO = irk_mod(iv)
     CASE('oysdeb_E')
       iv_oysdeb_E = irk_mod(iv)
     CASE('oysdeb_E_R')
       iv_oysdeb_E_R = irk_mod(iv)
#endif
#ifdef key_benthos
     CASE('mole_concentration_of_organic_detritus_expressed_as_nitrogen_in_benthos')
       iv_benth_N = irk_mod(iv)
     CASE('mole_concentration_of_organic_detritus_expressed_as_silicon_in_benthos')
       iv_benth_Si = irk_mod(iv)
     CASE('mole_concentration_of_organic_detritus_expressed_as_phosphorus_in_benthos')
       iv_benth_P = irk_mod(iv)
#if defined key_diatbenth 
     CASE('mole_concentration_of_benthic_diatoms_expressed_as_nitrogen_in_benthos')
       iv_benth_phyto_diat_N = irk_mod(iv)
#endif
#if defined key_NPbenth 
     CASE('mole_concentration_of_ammonium_in_benthos')
       iv_benth_NH4 = irk_mod(iv)         
     CASE('mole_concentration_of_dissolved_phosphorus_in_benthos')
       iv_benth_PO4 = irk_mod(iv)         
#endif
#if defined key_zostera
     CASE('mole_concentration_of_benthic_diatoms_expressed_as_nitrogen_in_benthos')
       iv_benth_phyto_diat_N = irk_mod(iv)
     CASE('mole_concentration_of_ammonium_in_benthos')
       iv_benth_NH4 = irk_mod(iv)         
     CASE('mole_concentration_of_dissolved_phosphorus_in_benthos')
       iv_benth_PO4 = irk_mod(iv)         
#endif

#endif
#ifdef key_suspensivores
     CASE('mole_concentration_of_suspensivores_expressed_as_nitrogen_in_benthos')
       iv_benth_susp_N = irk_mod(iv)
#endif
#ifdef key_benthos_gener
     CASE('mole_concentration_of_bacterie_meiofaune_expressed_as_nitrogen_in_benthos')
       iv_benth_meio_N = irk_mod(iv)
     CASE('mole_concentration_of_deposivore_herbivore_expressed_as_nitrogen_in_benthos')
       iv_benth_deposi_N = irk_mod(iv)
     CASE('mole_concentration_of_suspensivores_expressed_as_nitrogen_in_benthos')
       iv_benth_susp_N = irk_mod(iv)
     CASE('mole_concentration_of_carnivore_expressed_as_nitrogen_in_benthos')
       iv_benth_carn_N = irk_mod(iv)
#endif

#ifdef key_larve
    CASE('concentration_of_larvae_zone_1')
       iv_larve_1 = irk_mod(iv)
    CASE('concentration_of_larvae_zone_2')   
        iv_larve_2 = irk_mod(iv)
    CASE('concentration_of_larvae_zone_3')   
        iv_larve_3 = irk_mod(iv)
    CASE('concentration_of_larvae_zone_4')   
        iv_larve_4 = irk_mod(iv)
    CASE('concentration_of_larvae_zone_5')   
        iv_larve_5 = irk_mod(iv)
    CASE('concentration_of_larvae_zone_6')   
        iv_larve_6 = irk_mod(iv)
#endif
    CASE DEFAULT
#if defined key_N_tracer || defined key_P_tracer
     IF (icall ==1) THEN   ! advected variables in CROCO, all variable in MARS
#if defined key_N_tracer 
       DO is=1,nb_source_tracerN
         DO ivtra=1,nb_vara_tracerN ! not MARS : variables advected only
          IF(iv==isubs_tracer_N(ivtra,is)) THEN
            iv_tracer_N(ivtra,is)=irk_mod(iv)
            iv_signed_N(ivtra,is)=irk_mod(isubs_signed_N(ivtra,is))
#if defined key_age_tracer
            iv_age_N(ivtra,is)=irk_mod(iv+1)
#endif
            SELECT CASE( TRIM(ADJUSTL(ADJUSTR(name_var(isubs_signed_N(ivtra,is)))) ))
                CASE('mesozooplankton_nitrogen')
                   iv_zoo_meso_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_zoo_meso_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('microzooplankton_nitrogen')
                   iv_zoo_micr_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_zoo_micr_age_tra_N(is) = irk_mod(iv+1)
#endif
#if defined key_benthos
                CASE('organic_nitrogen_benth')
                   iv_benth_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_benth_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif
                CASE('ammonium')
                   iv_nutr_NH4_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_nutr_NH4_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('nitrate')
                   iv_nutr_NO3_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_nutr_NO3_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('nanopicoplankton_nitrogen')
                   iv_phyto_nano_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_nano_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('diatom_nitrogen')
                   iv_phyto_diat_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_diat_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('dinoflagellate_nitrogen')
                   iv_phyto_dino_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_dino_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('detrital_nitrogen')
                   iv_detr_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_detr_age_tra_N(is) = irk_mod(iv+1)
#endif
#if defined key_karenia
                CASE('karenia_nitrogen')
                   iv_phyto_karenia_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_karenia_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif
#if defined key_psnz
                CASE('pseudonitzschia_nitrogen')
                   iv_phyto_psnz_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_psnz_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif
#if defined key_phaeocystis
                CASE('colonial_phaeocystis_nitrogen')
                   iv_phyto_phaeocystis_colo_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_phaeocystis_colo_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('phaeocystis_cells_nitrogen')
                   iv_phyto_phaeocystis_cell_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_phaeocystis_cell_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif
#if defined key_ulvas
!!! ATTENTION 	: NOMS DES VARIABLES A VERIFIER
                CASE('ulves_nitrogen_benth')
                   iv_ulv_benth_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_ulv_benth_age_tra_N(is) = irk_mod(iv+1)
#endif
                CASE('ulves_nitrogen_susp')
                   iv_ulv_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_ulv_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif

            END SELECT
          ENDIF
         ENDDO
       ENDDO
#else 
!! tracer P
       DO is=1,nb_source_tracerP
         DO ivtra=1,nb_vara_tracerP ! not MARS : variables advected only
          IF(iv==isubs_tracer_P(ivtra,is)) THEN
            iv_tracer_P(ivtra,is)=irk_mod(iv)
            iv_signed_P(ivtra,is)=irk_mod(isubs_signed_P(ivtra,is))
#if defined key_age_tracer
            iv_age_P(ivtra,is)=irk_mod(iv+1)
#endif
            SELECT CASE( TRIM(ADJUSTL(ADJUSTR(name_var(isubs_signed_P(ivtra,is)))) ))
                CASE('mesozooplankton_nitrogen')
                   iv_zoo_meso_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_zoo_meso_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('microzooplankton_nitrogen')
                   iv_zoo_micr_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_zoo_micr_age_tra_P(is) = irk_mod(iv+1)
#endif
#if defined key_benthos
                CASE('detrital_phosphorus_benth')
                   iv_benth_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_benth_age_tra_P(is) = irk_mod(iv+1)
#endif
#endif
                CASE('dissolved_phosphate')
                   iv_nutr_PO4_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_nutr_PO4_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('sorbed_phosphate')
                   iv_nutr_Pads_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_nutr_Pads_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('detrital_phosphorus')
                   iv_detr_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_detr_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('nanopicoplankton_nitrogen')
                   iv_phyto_nano_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_nano_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('diatom_nitrogen')
                   iv_phyto_diat_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_diat_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('dinoflagellate_nitrogen')
                   iv_phyto_dino_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_dino_age_tra_P(is) = irk_mod(iv+1)
#endif
#if defined key_karenia
                CASE('karenia_phosphorus')
                   iv_phyto_karenia_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_karenia_age_tra_P(is) = irk_mod(iv+1)
#endif
#endif
#if defined key_psnz
                CASE('pseudonitzschia_nitrogen')
                   iv_phyto_psnz_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_psnz_age_tra_P(is) = irk_mod(iv+1)
#endif
#endif
#if defined key_phaeocystis
                CASE('colonial_phaeocystis_nitrogen')
                   iv_phyto_phaeocystis_colo_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_phaeocystis_colo_age_tra_P(is) = irk_mod(iv+1)
#endif
                CASE('phaeocystis_cells_nitrogen')
                   iv_phyto_phaeocystis_cell_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_phyto_phaeocystis_cell_age_tra_P(is) = irk_mod(iv+1)
#endif
#endif
            END SELECT
          ENDIF
         ENDDO
       ENDDO
#endif

    ELSE IF (icall == 2) THEN    ! fixed variables (not in MARS becuse only one array for all variables in MARS)

    ! in incellwat_bloom_nitrogentracer : working with c() which bringing together advected and fixed variables as in MARS
#if defined key_N_tracer 
       DO is=1,nb_source_tracerN
         DO ivtra=nb_vara_tracerN+1,nb_vara_tracerN+nb_fix_tracerN
          IF(iv==nv_adv+isubs_tracer_N(ivtra,is)) THEN
            iv_tracer_N(ivtra,is)=irk_mod(iv)
            iv_signed_N(ivtra,is)=irk_mod(nv_adv+isubs_signed_N(ivtra,is))
#if defined key_age_tracer
            iv_age_N(ivtra,is)=irk_mod(iv+1)
#endif
            SELECT CASE( TRIM(ADJUSTL(ADJUSTR(name_var_fix(isubs_signed_N(ivtra,is)))) ))
#if defined key_benthos
                CASE('organic_nitrogen_benth')
                   iv_benth_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_benth_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif
#if defined key_ulvas
!!! ATTENTION 	: NOMS DES VARIABLES A VERIFIER
                CASE('ulves_nitrogen_benth')
                   iv_ulv_benth_tra_N(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_ulv_benth_age_tra_N(is) = irk_mod(iv+1)
#endif
#endif

            END SELECT
          ENDIF
         ENDDO
       ENDDO
#else 
!! tracer P
       DO is=1,nb_source_tracerP
         DO ivtra=nb_vara_tracerP+1,nb_vara_tracerP+nb_fix_tracerP  
          IF(iv==nv_adv+isubs_tracer_P(ivtra,is)) THEN
            iv_tracer_P(ivtra,is)=irk_mod(iv)
            iv_signed_P(ivtra,is)=irk_mod(nv_adv+isubs_signed_P(ivtra,is))
#if defined key_age_tracer
            iv_age_P(ivtra,is)=irk_mod(iv+1)
#endif
            SELECT CASE( TRIM(ADJUSTL(ADJUSTR(name_var_fix(isubs_signed_P(ivtra,is)))) ))
#if defined key_benthos
                CASE('detrital_phosphorus_benth')
                   iv_benth_tra_P(is) = irk_mod(iv)
#if defined key_age_tracer
                   iv_benth_age_tra_P(is) = irk_mod(iv+1)
#endif
#endif
            END SELECT
          ENDIF
         ENDDO
       ENDDO
#endif
#ifdef key_benthic
    ELSE IF (icall == 3) THEN    ! benthic variables (not in MARS)
#endif
    ENDIF

#else
       MPI_master_only write(iscreenlog,*) 'Invalid standard name of variable : ',TRIM(ADJUSTL(ADJUSTR(standname)))
#endif
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
#if defined GAMELAG
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_detr_N = ',iv_diss_detr_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_detrR_N = ',iv_diss_detrR_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_detr_P = ',iv_diss_detr_P
       MPI_master_only WRITE(iscreenlog,*) 'iv_diss_detrR_P = ',iv_diss_detrR_P
#endif
#if defined key_psnz
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_N = ',iv_phyto_psnz_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_Si = ',iv_phyto_psnz_Si
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_ad = ',iv_phyto_psnz_ad
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_pp = ',iv_phyto_psnz_pp          
#endif
#if defined key_karenia
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_C = ',iv_phyto_karenia_C
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_N = ',iv_phyto_karenia_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_P = ',iv_phyto_karenia_P
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_pp = ',iv_phyto_karenia_pp          
#endif
#if defined key_phaeocystis
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_colo_N = ',iv_phyto_phaeocystis_colo_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_cell_N = ',iv_phyto_phaeocystis_cell_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_mucus = ',iv_phyto_phaeocystis_mucus
       MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_pp = ',iv_phyto_phaeocystis_pp          
#endif
#ifdef key_zoo_prod
       MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_micr_ps = ',iv_zoo_micr_ps
       MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_meso_ps = ',iv_zoo_meso_ps
#endif
#ifdef key_BLOOM_insed
       MPI_master_only WRITE(iscreenlog,*) 'iv_ODU = ',iv_ODU
       MPI_master_only WRITE(iscreenlog,*) 'iv_PFe = ',iv_PFe
       MPI_master_only WRITE(iscreenlog,*) 'iv_detrR_N = ',iv_detrR_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_detrR_P = ',iv_detrR_P
#endif

#ifdef key_ulvas
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_sus_N = ',iv_ulv_sus_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_benth_N = ',iv_ulv_benth_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_susdrywght = ',iv_ulv_susdrywght
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_benthdrywght = ',iv_ulv_benthdrywght 
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_sus_P = ',iv_ulv_sus_P
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_benth_P = ',iv_ulv_benth_P 
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_suspp = ',iv_ulv_suspp
       MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_benthpp = ',iv_ulv_benthpp  
#endif
#ifdef key_zostera
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_LB = ',iv_zost_LB
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_RB = ',iv_zost_RB
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_D = ',iv_zost_D
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_LN = ',iv_zost_LN
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_RN = ',iv_zost_RN
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_LP = ',iv_zost_LP
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_RP = ',iv_zost_RP
       MPI_master_only WRITE(iscreenlog,*) 'iv_detr_zost_N = ',iv_detr_zost_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_detr_zost_P = ',iv_detr_zost_P
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_benth_N = ',iv_zost_benth_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_benth_P = ',iv_zost_benth_P
       MPI_master_only WRITE(iscreenlog,*) 'iv_zost_seed = ',iv_zost_seed
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_seed = ',iv_zost_benth_seed
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_NH4 = ',iv_benth_NH4
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_PO4 = ',iv_benth_PO4
!      MPI_master_only WRITE(iscreenlog,*) 'iv_zost_pp = ',iv_zost_pp
!      MPI_master_only WRITE(iscreenlog,*) 'iv_zost_LNuptake = ',iv_zost_LNuptake
!      MPI_master_only WRITE(iscreenlog,*) 'iv_zost_LPuptake = ',iv_zost_LPuptake
!      MPI_master_only WRITE(iscreenlog,*) 'iv_zost_RNuptake = ',iv_zost_RNuptake
!      MPI_master_only WRITE(iscreenlog,*) 'iv_zost_RPuptake = ',iv_zost_RPuptake
#endif  
#ifdef key_oyster_SFG
       MPI_master_only WRITE(iscreenlog,*) 'iv_oys_so = ',iv_oys_so
       MPI_master_only WRITE(iscreenlog,*) 'iv_oys_go = ',iv_oys_go
       MPI_master_only WRITE(iscreenlog,*) 'iv_oys_co = ',iv_oys_co
#endif
#ifdef key_oxygen
       MPI_master_only WRITE(iscreenlog,*) 'iv_oxygen = ',iv_oxygen
#endif
#ifdef key_oyster_DEB
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_res = ',iv_oysdeb_res
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_str = ',iv_oysdeb_str
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_gon = ',iv_oysdeb_gon
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_res2 = ',iv_oysdeb_res2
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_str2 = ',iv_oysdeb_str2
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_gon2 = ',iv_oysdeb_gon2
!     WRITE(iscreenlog,*) 'iv_oysdeb_res3 = ',iv_oysdeb_res3
!     WRITE(iscreenlog,*) 'iv_oysdeb_str3 = ',iv_oysdeb_str3
!     WRITE(iscreenlog,*) 'iv_oysdeb_gon3 = ',iv_oysdeb_gon3
#endif
#ifdef key_oyster_DEB
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb = ',iv_oysdeb
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_E = ',iv_oysdeb_E
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_E_V = ',iv_oysdeb_E_V
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_E_R = ',iv_oysdeb_E_R
     MPI_master_only WRITE(iscreenlog,*) 'iv_oysdeb_E_GO = ',iv_oysdeb_E_GO
#endif
#ifdef key_benthos
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_N = ',iv_benth_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_Si = ',iv_benth_Si
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_P = ',iv_benth_P
#if defined key_diatbenth 
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_phyto_diat_N = ',iv_benth_phyto_diat_N
#endif
#if defined key_NPbenth 
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_NH4 = ',iv_benth_NH4
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_PO4 = ',iv_benth_PO4
#endif
#if defined key_zostera
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_phyto_diat_N = ',iv_benth_phyto_diat_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_NH4 = ',iv_benth_NH4
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_PO4 = ',iv_benth_PO4
#endif
#endif
#ifdef key_suspensivores
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_susp_N = ',iv_benth_susp_N
#endif
#ifdef key_benthos_gener
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_meio_N = ',iv_benth_meio_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_deposi_N = ',iv_benth_deposi_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_susp_N = ',iv_benth_susp_N
       MPI_master_only WRITE(iscreenlog,*) 'iv_benth_carn_N = ',iv_benth_carn_N
#endif
#ifdef key_larve
       MPI_master_only WRITE(iscreenlog,*) 'iv_larve_1 = ',iv_larve_1
       MPI_master_only WRITE(iscreenlog,*) 'iv_larve_2 = ',iv_larve_2
       MPI_master_only WRITE(iscreenlog,*) 'iv_larve_3 = ',iv_larve_3
       MPI_master_only WRITE(iscreenlog,*) 'iv_larve_4 = ',iv_larve_4
       MPI_master_only WRITE(iscreenlog,*) 'iv_larve_5 = ',iv_larve_5
       MPI_master_only WRITE(iscreenlog,*) 'iv_larve_6 = ',iv_larve_6
#endif    
#if defined key_N_tracer
#if defined key_BLOOM_opt2
       MPI_master_only WRITE(ierrorlog,*)'ERROR : option 2 of key_BLOOM (-key_BLOOM_opt2) is incompatible with key_N_tracer'
       MPI_master_only WRITE(ierrorlog,*)'The simulation is stopped'
       CALL_MPI MPI_FINALIZE(IERR_MPI)
       STOP
#endif
       DO is=1,nb_source_tracerN
         MPI_master_only WRITE(iscreenlog,*) 'tracer_N source ',is
         DO ivtra=1,nb_vara_tracerN 

          MPI_master_only WRITE(iscreenlog,*) iv_tracer_N(ivtra,is) ,'variable : ',  &
                                             TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_N(ivtra,is)))))
#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) iv_age_N(ivtra,is) ,'variable : ',  &
                                             TRIM(ADJUSTL(ADJUSTR(name_var(isubs_age_N(ivtra,is)))))
#endif
         ENDDO
         DO ivtra=nb_vara_tracerN+1,nb_vara_tracerN+nb_fix_tracerN 
          MPI_master_only WRITE(iscreenlog,*) iv_tracer_N(ivtra,is) ,'variable : ',  &
                                             TRIM(ADJUSTL(ADJUSTR(name_var_fix(isubs_tracer_N(ivtra,is)))))
#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) iv_age_N(ivtra,is) ,'variable : ',  &
                                             TRIM(ADJUSTL(ADJUSTR(name_var_fix(isubs_age_N(ivtra,is)))))
#endif
         ENDDO
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_NH4_tra_N(',is,' = ',iv_nutr_NH4_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_NO3_tra_N(',is,' = ',iv_nutr_NO3_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_nano_tra_N(',is,' = ',iv_phyto_nano_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_diat_tra_N(',is,' = ',iv_phyto_diat_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_dino_tra_N(',is,' = ',iv_phyto_dino_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_micr_tra_N(',is,' = ',iv_zoo_micr_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_meso_tra_N(',is,' = ',iv_zoo_meso_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_detr_tra_N(',is,' = ',iv_detr_tra_N(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_NH4_age_tra_N(',is,' = ',iv_nutr_NH4_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_NO3_age_tra_N(',is,' = ',iv_nutr_NO3_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_nano_age_tra_N(',is,' = ',iv_phyto_nano_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_diat_age_tra_N(',is,' = ',iv_phyto_diat_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_dino_age_tra_N(',is,' = ',iv_phyto_dino_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_micr_age_tra_N(',is,' = ',iv_zoo_micr_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_meso_age_tra_N(',is,' = ',iv_zoo_meso_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_detr_age_tra_N(',is,' = ',iv_detr_age_tra_N(is)          
#endif
#if defined key_benthos
         MPI_master_only WRITE(iscreenlog,*) 'iv_benth_tra_N(',is,' = ',iv_benth_tra_N(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_benth_age_tra_N(',is,' = ',iv_benth_age_tra_N(is)
#endif
#endif
#if defined key_psnz
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_tra_N(',is,' = ',iv_phyto_psnz_tra_N(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_age_tra_N(',is,' = ',iv_phyto_psnz_age_tra_N(is)
#endif
#endif
#if defined key_phaeocystis
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_colo_tra_N(',is,' = ',iv_phyto_phaeocystis_colo_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_cell_tra_N(',is,' = ',iv_phyto_phaeocystis_cell_tra_N(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_colo_age_tra_N(',is,' = ',iv_phyto_phaeocystis_colo_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_cell_age_tra_N(',is,' = ',iv_phyto_phaeocystis_cell_age_tra_N(is)
#endif
#endif
#if defined key_karenia
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_tra_N(',is,' = ',iv_phyto_karenia_tra_N(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_age_tra_N(',is,' = ',iv_phyto_karenia_age_tra_N(is)
#endif
#endif
#if defined key_ulvas
         MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_benth_tra_N(',is,' = ',iv_ulv_benth_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_tra_N(',is,' = ',iv_ulv_tra_N(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_benth_age_tra_N(',is,' = ',iv_ulv_benth_age_tra_N(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_ulv_age_tra_N(',is,' = ',iv_ulv_age_tra_N(is)

#endif
#endif
       ENDDO
#endif
#if defined key_P_tracer
#if defined key_BLOOM_opt2
       MPI_master_only WRITE(ierrorlog,*)'ERROR : option 2 of key_BLOOM (-key_BLOOM_opt2) is incompatible with key_P_tracer'
       MPI_master_only WRITE(ierrorlog,*)'The simulation is stopped'
       CALL_MPI MPI_FINALIZE(IERR_MPI)
       STOP
#endif
       DO is=1,nb_source_tracerP
         MPI_master_only WRITE(iscreenlog,*) 'tracer_P source ',is
         DO ivtra=1,nb_vara_tracerP 

          MPI_master_only WRITE(iscreenlog,*) iv_tracer_P(ivtra,is) ,'variable : ',  &
                              TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_P(ivtra,is)))))
#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) iv_age_P(ivtra,is) ,'variable : ',  &
                              TRIM(ADJUSTL(ADJUSTR(name_var(isubs_age_P(ivtra,is)))))
#endif
         ENDDO
         DO ivtra=nb_vara_tracerP+1,nb_vara_tracerP+nb_fix_tracerP 
          MPI_master_only WRITE(iscreenlog,*) iv_tracer_P(ivtra,is) ,'variable : ',  &
                                             TRIM(ADJUSTL(ADJUSTR(name_var_fix(isubs_tracer_P(ivtra,is)))))
#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) iv_age_P(ivtra,is) ,'variable : ',  &
                                             TRIM(ADJUSTL(ADJUSTR(name_var_fix(isubs_age_P(ivtra,is)))))
#endif
         ENDDO
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_PO4_tra_P(',is,' = ',iv_nutr_PO4_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_Pads_tra_P(',is,' = ',iv_nutr_Pads_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_nano_tra_P(',is,' = ',iv_phyto_nano_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_diat_tra_P(',is,' = ',iv_phyto_diat_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_dino_tra_P(',is,' = ',iv_phyto_dino_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_micr_tra_P(',is,' = ',iv_zoo_micr_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_meso_tra_P(',is,' = ',iv_zoo_meso_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_detr_tra_P(',is,' = ',iv_detr_tra_P(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_PO4_age_tra_P(',is,' = ',iv_nutr_PO4_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_nutr_Pads_age_tra_P(',is,' = ',iv_nutr_Pads_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_nano_age_tra_P(',is,' = ',iv_phyto_nano_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_diat_age_tra_P(',is,' = ',iv_phyto_diat_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_dino_age_tra_P(',is,' = ',iv_phyto_dino_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_micr_age_tra_P(',is,' = ',iv_zoo_micr_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_zoo_meso_age_tra_P(',is,' = ',iv_zoo_meso_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_detr_age_tra_P(',is,' = ',iv_detr_age_tra_P(is)          
#endif
#if defined key_benthos
         MPI_master_only WRITE(iscreenlog,*) 'iv_benth_tra_P(',is,' = ',iv_benth_tra_P(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_benth_age_tra_P(',is,' = ',iv_benth_age_tra_P(is)
#endif
#endif
#if defined key_psnz
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_tra_P(',is,' = ',iv_phyto_psnz_tra_P(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_psnz_age_tra_P(',is,' = ',iv_phyto_psnz_age_tra_P(is)
#endif
#endif
#if defined key_phaeocystis
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_colo_tra_P(',is,' = ',iv_phyto_phaeocystis_colo_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_cell_tra_P(',is,' = ',iv_phyto_phaeocystis_cell_tra_P(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_colo_age_tra_P(',is,' = ',iv_phyto_phaeocystis_colo_age_tra_P(is)
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_phaeocystis_cell_age_tra_P(',is,' = ',iv_phyto_phaeocystis_cell_age_tra_P(is)
#endif
#endif
#if defined key_karenia
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_tra_P(',is,' = ',iv_phyto_karenia_tra_P(is)
#if defined key_age_tracer
         MPI_master_only WRITE(iscreenlog,*) 'iv_phyto_karenia_age_tra_P(',is,' = ',iv_phyto_karenia_age_tra_P(is)
#endif
#endif
       ENDDO
#endif
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
#if defined key_psnz
                iv_phyto_psnz_N + iv_phyto_psnz_Si + iv_phyto_psnz_ad + iv_phyto_psnz_pp + &
#endif
#if defined key_karenia
                iv_phyto_karenia_C + iv_phyto_karenia_N + iv_phyto_karenia_P + iv_phyto_karenia_pp + &
#endif
#if defined key_phaeocystis
                iv_phyto_phaeocystis_colo_N + iv_phyto_phaeocystis_cell_N +      &
                iv_phyto_phaeocystis_mucus + iv_phyto_phaeocystis_pp + &
#endif
#if defined key_zoo_prod
                iv_zoo_micr_ps + iv_zoo_meso_ps + &
#endif
#if defined key_BLOOM_insed
                iv_ODU + iv_PFe + iv_detrR_N + iv_detrR_P            &
#endif
#ifdef key_ulvas
                iv_ulv_sus_N + iv_ulv_benth_N +iv_ulv_susdrywght + &
                iv_ulv_benthdrywght + iv_ulv_sus_P + iv_ulv_benth_P + &
                iv_ulv_suspp + iv_ulv_benthpp + & 
#endif
#ifdef key_zostera
               iv_zost_LB+iv_zost_RB+iv_zost_D+iv_zost_LN+iv_zost_RN+iv_zost_LP+iv_zost_RP+ &
               iv_detr_zost_N+iv_detr_zost_P+iv_zost_benth_N+iv_zost_benth_P+iv_zost_seed+iv_zost_benth_seed +   &
!               iv_zost_pp+iv_zost_LNuptake+iv_zost_LPuptake+iv_zost_RNuptake+iv_zost_RPuptake
#endif
#ifdef key_oyster_SFG
                iv_oys_so+iv_oys_go+iv_oys_co + &
#endif  
#ifdef key_oxygen
                iv_oxygen + &
#endif
#ifdef key_oyster_DEB
                iv_oysdeb_res+iv_oysdeb_str+iv_oysdeb_gon + &
                iv_oysdeb_res2+iv_oysdeb_str2+iv_oysdeb_gon2 + &
!                iv_oysdeb_res3+iv_oysdeb_str3+iv_oysdeb_gon3 + &
#endif
#ifdef key_oyster_DEB_GAMELAG
                iv_oysdeb+iv_oysdeb_E+iv_oysdeb_E_V+iv_oysdeb_E_R+iv_oysdeb_E_GO + &
#endif
#ifdef key_benthos
                iv_benth_N+iv_benth_Si+iv_benth_P + &
#if defined key_diatbenth
                iv_benth_phyto_diat_N + &
#if defined key_NPbenth
                iv_benth_NH4+iv_benth_PO4 + &
#endif
#elif defined  key_zostera
                iv_benth_phyto_diat_N +iv_benth_NH4+iv_benth_PO4 + &
#endif
#endif
#if defined key_microtracers
                iv_tracer_depth1 + iv_tracer_depth2 + iv_tracer_depth3 + &
                iv_tracer_agedepth1 + iv_tracer_agedepth2 + iv_tracer_agedepth3 + &
                iv_tracer_cumuldepth1 + iv_tracer_cumuldepth2 + iv_tracer_cumuldepth3 + &
#endif
                iv_spim + iv_phyto_nano_pp + iv_phyto_diat_pp + iv_phyto_dino_pp

#if defined key_N_tracer
       DO is=1,nb_source_tracerN
         DO ivtra=1,nb_var_tracerN 
          sumverif = sumverif + iv_tracer_N(ivtra,is)
#if defined key_age_tracer
          sumverif = sumverif + iv_age_N(ivtra,is) 
#endif
         ENDDO
       ENDDO
#endif
#if defined key_P_tracer
       DO is=1,nb_source_tracerP
         DO ivtra=1,nb_var_tracerP 
          sumverif = sumverif + iv_tracer_P(ivtra,is)
#if defined key_age_tracer
          sumverif = sumverif + iv_age_P(ivtra,is) 
#endif
         ENDDO
       ENDDO
#endif 
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
#if defined key_psnz
  ivfix_cumulprod_first=min(ivfix_cumulprod_first,iv_phyto_psnz_pp)
  ivfix_cumulprod_last=max(ivfix_cumulprod_last,iv_phyto_psnz_pp)
#endif
#if defined key_karenia
  ivfix_cumulprod_first=min(ivfix_cumulprod_first,iv_phyto_karenia_pp)
  ivfix_cumulprod_last=max(ivfix_cumulprod_last,iv_phyto_karenia_pp)
#endif
#if defined key_phaeocystis
  ivfix_cumulprod_first=min(ivfix_cumulprod_first,iv_phyto_phaeocystis_pp)
  ivfix_cumulprod_last=max(ivfix_cumulprod_last,iv_phyto_phaeocystis_pp)
#endif
#if defined key_zoo_prod
  ivfix_cumulprod_first=min(ivfix_cumulprod_first,iv_zoo_micr_ps)
  ivfix_cumulprod_first=min(ivfix_cumulprod_first,iv_zoo_meso_ps)
  ivfix_cumulprod_last=max(ivfix_cumulprod_last,iv_zoo_micr_ps)
  ivfix_cumulprod_last=max(ivfix_cumulprod_last,iv_zoo_meso_ps)
#endif
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
#if defined key_N_tracer || defined key_P_tracer
   INTEGER               :: is,ivtra,idiag
#endif

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
     CASE('Fluxcum_aerobic_miner_Norg')
       id_remin_aerN = id
     CASE('Fluxcum_aerobic_miner_Porg')
       id_remin_aerP = id
     CASE('Fluxcum_aerobic_miner_Si')
       id_remin_aerSi = id
     CASE('Fluxcum_anaerobic_miner_Norg')
       id_remin_anaerN = id
     CASE('Fluxcum_anaerobic_miner_Porg')
       id_remin_anaerP = id
     CASE('Fluxcum_nitrate_miner_Norg')
       id_remin_nitrateN = id
     CASE('Fluxcum_drna_miner_Norg')
       id_remin_drnaN = id
     CASE('Fluxcum_denit_miner_Norg')
       id_remin_denitN = id
     CASE('Fluxcum_nitrate_miner_Porg')
       id_remin_nitrateP = id
     CASE('Fluxcum_nitrification')
       id_nitrif = id
     CASE('Fluxcum_ODU_oxyd_solid')
       id_oxyd_solid_ODU = id
     CASE('Fluxcum_adsorb_desorb_P')
       id_adsor_desorb_P = id
     CASE('Fluxcum_dissolution_PFe')
       id_dissol_PFe = id
     CASE('Fluxcum_precipitation_P')
       id_precipit_P = id
     CASE('Fluxcum_precipitation_Si')
       id_precipit_Si = id
     CASE('Fluxcum_mortality_phyto_sed')
       id_morta_phyto = id
     CASE('Fluxcum_benthic_filtre_grazing')
       id_filtr_benth = id
     CASE('Fluxcum_aeration_sediment')
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
     CASE('light_limitation_of_nanoflagellate_growth')
       id_nano_limlight = id
     CASE('nitrogen_limitation_of_nanoflagellate_growth')
       id_nano_limN = id
     CASE('phosphate_limitation_of_nanoflagellate_growth')
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
#if defined key_psnz
     CASE('maximum_pseudonitzschia_mass_concentration_in_sea_water')
       id_psnz_max = id
     CASE('date_of_maximum_pseudonitzschia_mass_concentration_in_sea_water')
       id_psnz_datemax = id
     CASE('depth_of_maximum_pseudonitzschia_mass_concentration_in_sea_water')
       id_psnz_depthmax = id
     CASE('light_limitation_of_pseudonitzschia_growth')
       id_psnz_limlight = id
     CASE('nitrogen_limitation_of_pseudonitzschia_growth')
       id_psnz_limN = id
     CASE('silicon_limitation_of_pseudonitzschia_growth')
       id_psnz_limSi = id
     CASE('phosphate_limitation_of_pseudonitzschia_growth')
       id_psnz_limP = id
     CASE('cumulated_production_of_pseudonitzschia_in_sea_water_column_expressed_in_carbon')
       id_psnz_columnprod=id
#endif
#if defined key_karenia
     CASE('maximum_karenia_mass_concentration_in_sea_water')
       id_karenia_max = id
     CASE('date_of_maximum_karenia_mass_concentration_in_sea_water')
       id_karenia_datemax = id
     CASE('depth_of_maximum_karenia_mass_concentration_in_sea_water')
       id_karenia_depthmax = id
     CASE('light_limitation_of_karenia_growth')
       id_karenia_limlight = id
     CASE('nitrogen_limitation_of_karenia_growth')
       id_karenia_limN = id
     CASE('phosphate_limitation_of_karenia_growth')
       id_karenia_limP = id
     CASE('cumulated_production_of_karenia_in_sea_water_column_expressed_in_carbon')
       id_karenia_columnprod=id
#endif
#if defined key_phaeocystis
     CASE('maximum_phaeocystis_mass_concentration_in_sea_water')
       id_phaeocystis_max = id
     CASE('date_of_maximum_phaeocystis_mass_concentration_in_sea_water')
       id_phaeocystis_datemax = id
     CASE('depth_of_maximum_phaeocystis_mass_concentration_in_sea_water')
       id_phaeocystis_depthmax = id
     CASE('light_limitation_of_phaeocystis_growth')
       id_phaeocystis_limlight = id
     CASE('nitrogen_limitation_of_phaeocystis_growth')
       id_phaeocystis_limN = id
     CASE('phosphate_limitation_of_phaeocystis_growth')
       id_phaeocystis_limP = id
     CASE('cumulated_production_of_phaeocystis_in_sea_water_column_expressed_in_carbon')
       id_phaeocystis_columnprod=id
#endif
#if defined key_zoo_prod
     CASE('cumulated_production_of_microzooplankton_in_sea_water_column_expressed_in_carbon')
       id_zoo_micr_columnprod=id
     CASE('cumulated_production_of_mesozooplankton_in_sea_water_column_expressed_in_carbon')
       id_zoo_meso_columnprod=id
     CASE('cumulated_total_zooplancton_production_in_sea_water_column_expressed_in_carbon')
       id_zoo_columnprodtotal=id
#endif
#ifdef key_ulvas
     CASE('suspended_ulva_growth')
       id_ulv_ratiosus=id
     CASE('suspended_ulva_mortality')
       id_ulv_susmort=id
     CASE('light_limitation_of_suspended_ulva_growth')
       id_ulv_suslimlight=id
     CASE('nitrogen_limitation_of_suspended_ulva_growth')
       id_ulv_suslimN=id
     CASE('phosphate_limitation_of_suspended_ulva_growth')
       id_ulv_suslimP=id
     CASE('nitrogen_pumping_by_suspended_ulva')
       id_ulv_suspumpN=id
     CASE('phosphate_pumping_by_suspended_ulva')
       id_ulv_suspumpP=id
     CASE('ulva_sinking_velocity_in_sea_water')
       id_ulv_settling=id
     CASE('settled_ulva_growth')
       id_ulv_ratiobenth=id
     CASE('settled_ulva_mortality')
       id_ulv_benthmort=id
     CASE('light_limitation_of_settled_ulva_growth')
       id_ulv_benthlimlight=id
     CASE('nitrogen_limitation_of_settled_ulva_growth')
      id_ulv_benthlimN=id
     CASE('phosphate_limitation_of_settled_ulva_growth')
      id_ulv_benthlimP=id
     CASE('nitrogen_pumping_by_settled_ulva')
      id_ulv_benthpumpN=id
     CASE('phosphate_pumping_by_settled_ulva')
      id_ulv_benthpumpP=id
     CASE('settl_resusp_bilan')
      id_ulv_resusbil=id
     CASE('nitrogen_resusp')
      id_ulv_resusN=id
     CASE('phosphate_resusp')
      id_ulv_resusP=id
#endif
#ifdef key_zostera
!     CASE('zostera_leaf_biomass_expressed_as_gDW_per_m2')
!       id_zost_LBDW=id
!     CASE('zostera_root_biomass_expressed_as_gDW_per_m2')
!       id_zost_RBDW=id
     CASE('light_at_the_top_of_the_canopy')
       id_zost_Qcan=id
     CASE('zostera_leaf_area_index')
       id_zost_LAI=id
     CASE('light_limitation_for_zostera_production')
       id_zost_limlum=id
     CASE('nitrogen_limitation_for_zostera_production')
       id_zost_limprodLN=id
     CASE('phosphorous_limitation_for_zostera_production')
       id_zost_limprodLP=id
     CASE('nitrogen_limitation_for_zostera_recruitment')
       id_zost_limRN=id
     CASE('phosphorous_limitation_for_zostera_recruitment')
       id_zost_limRP=id
!     CASE('temperature_effect_on_zostera_production')
!       id_zost_effetchaleurzprod=id
!     CASE('temperature_effect_on_leaf_respiration')
!       id_zost_effetchaleurzrespf=id
!     CASE('temperature_effect_on_root_respiration')
!       id_zost_effetchaleurzrespr=id
!     CASE('ammonium_limitation_for_leaf_absorption')
!       id_zost_limabsLnh4=id
!     CASE('nitrate_limitation_for_leaf_absorption')
!       id_zost_limabsLno3=id
!     CASE('phosphate_limitation_for_leaf_absorption')
!       id_zost_limabsLpo4=id
!     CASE('ammonium_limitation_for_root_absorption')
!       id_zost_limabsRnh4=id
!     CASE('phosphate_limitation_for_root_absorption')
!       id_zost_limabsRpo4=id
     CASE('selfshading_limitation_for_zostera_recruitment')
       id_zost_limselfshad=id
     CASE('root_biomass_limitation_for_zostera_recruitment')
       id_zost_limRB=id
!     CASE('wind_effect_on_leaf_mortality')
!       id_zost_effetvent=id
!     CASE('zostera_seed_germination_rate')
!       id_zost_Sgerm=id
!     CASE('zostera_sexual_reproduction_effort')
!       id_zost_ERS=id
     CASE('daily_zostera_production')
       id_zost_prod=id
#endif
#ifdef key_oyster_SFG
     CASE('occurence_of_spawn_coq')
       id_oys_tcoq=id
     CASE('occurence_of_spawn_gam')
       id_oys_tgam=id
     CASE('gain_coq')
       id_oys_gaincoq=id
     CASE('gain_soma') 
       id_oys_gainsoma=id
     CASE('gain_reprod')
       id_oys_gainrepr=id
     CASE('organic_ingered_fraction')
       id_oys_orgafring=id
     CASE('bilansoma')
       id_oys_bilansoma=id
     CASE('bilancoq')
       id_oys_bilancoq=id
     CASE('oyster_filtration')
       id_oys_filt=id
     CASE('oyster_respiration')
       id_oys_resp=id
     CASE('oyster_gameto')
       id_oys_gameto=id
     CASE('oyster_thresh')
       id_oys_seuilponte=id
     CASE('reposcoq')
       id_oys_reposcoq=id
     CASE('reposgam')
       id_oys_reposgam=id
     CASE('ponte')
       id_oys_ponte=id
     CASE('absorg')
       id_oys_absorg=id   
#endif
#ifdef key_oxygen
     CASE('oxygen_saturation')
       id_oxy_sat=id  
#endif
#if defined key_oyster_DEB  
     CASE('oyster_filtration')
       id_oysDEB_FRpop=id
     CASE('oyster_pseudofeces')
       id_oysDEB_PFpop=id
     CASE('oyster_feces')
       id_oysDEB_Fpop=id
     CASE('oyster_dry_weigth')
        id_oysDEB_DWtot=id
     CASE('variable_chlorophylle')
        id_oysDEB_kchlvble=id
     CASE('oyster_gonad-soma_ratio')
        id_oysDEB_ER=id
     CASE('oyster_filtration_rate')
        id_oysDEB_filtrate=id
     CASE('total_weigth')
        id_oysDEB_Wtot=id
     CASE('oyster_stress_indice')
        id_oysDEB_Istress=id
     CASE('oyster_total_length')
        id_oysDEB_Lgtot=id
     CASE('chloro_carbone_ratio')
        id_oysDEB_ChlC=id
 !    CASE('oyster_dry_weigth2')
 !       id_oysDEB_DWtot2=id
 !    CASE('total_weigth2')
!        id_oysDEB_Wtot2=id
!     CASE('oyster_total_length2')
!        id_oysDEB_Lgtot2=id
!     CASE('oyster_dry_weigth3')
!        id_oysDEB_DWtot3=id
!     CASE('total_weigth3')
!        id_oysDEB_Wtot3=id
!     CASE('oyster_total_length3')
!        id_oysDEB_Lgtot3=id
     CASE('number_oyster_co1')
        id_oysDEB_nbhuit=id
!     CASE('number_oyster_co2')
!        id_oysDEB_nbhuit2=id
#endif
#if defined key_microtracers
     CASE('maximum_microtraceur_1_mass_concentration_in_sea_water')
       id_tracer_maxdepth1 = id
     CASE('mean_microtraceur_1_mass_concentration_in_sea_water')
       id_tracer_meandepth1 = id
     CASE('age_of_microtraceur_1_mass_concentration_in_sea_water')
       id_tracer_agedepth1 = id       
     CASE('maximum_microtraceur_2_mass_concentration_in_sea_water')
       id_tracer_maxdepth2 = id
     CASE('mean_microtraceur_2_mass_concentration_in_sea_water')
       id_tracer_meandepth2 = id
     CASE('age_of_microtraceur_2_mass_concentration_in_sea_water')
       id_tracer_agedepth2 = id              
     CASE('maximum_microtraceur_3_mass_concentration_in_sea_water')
       id_tracer_maxdepth3 = id
     CASE('mean_microtraceur_3_mass_concentration_in_sea_water')
       id_tracer_meandepth3 = id
     CASE('age_of_microtraceur_3_mass_concentration_in_sea_water')
       id_tracer_agedepth3 = id              
#endif     
     CASE DEFAULT
#if ! defined key_N_tracer && ! defined key_P_tracer
       !IF_MPI (MASTER) THEN
             MPI_master_only WRITE(iscreenlog,*) 'Invalid standard name of diagnostic variable : ', &
                                             TRIM(ADJUSTL(ADJUSTR(standname)))
       !ENDIF_MPI
#endif     
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
#if defined key_psnz
     WRITE(iscreenlog,*) 'id_psnz_max = ',id_psnz_max
     WRITE(iscreenlog,*) 'id_psnz_datemax = ',id_psnz_datemax
     WRITE(iscreenlog,*) 'id_psnz_depthmax = ',id_psnz_depthmax
     WRITE(iscreenlog,*) 'id_psnz_limlight = ',id_psnz_limlight
     WRITE(iscreenlog,*) 'id_psnz_limN = ',id_psnz_limN
     WRITE(iscreenlog,*) 'id_psnz_limP = ',id_psnz_limP
     WRITE(iscreenlog,*) 'id_psnz_limSi = ',id_psnz_limSi
     WRITE(iscreenlog,*) 'id_psnz_columnprod = ',id_psnz_columnprod
#endif
#if defined key_karenia
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_max = ',id_karenia_max
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_datemax = ',id_karenia_datemax
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_depthmax = ',id_karenia_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_limlight = ',id_karenia_limlight
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_limN = ',id_karenia_limN
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_limP = ',id_karenia_limP
     MPI_master_only WRITE(iscreenlog,*) 'id_karenia_columnprod = ',id_karenia_columnprod
#endif
#if defined key_phaeocystis
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_max = ',id_phaeocystis_max
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_datemax = ',id_phaeocystis_datemax
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_depthmax = ',id_phaeocystis_depthmax
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_limlight = ',id_phaeocystis_limlight
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_limN = ',id_phaeocystis_limN
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_limP = ',id_phaeocystis_limP
     MPI_master_only WRITE(iscreenlog,*) 'id_phaeocystis_columnprod = ',id_phaeocystis_columnprod
#endif
#if defined key_zoo_prod
     MPI_master_only WRITE(iscreenlog,*) 'id_zoo_micr_columnprod = ',id_zoo_micr_columnprod
     MPI_master_only WRITE(iscreenlog,*) 'id_zoo_meso_columnprod = ',id_zoo_meso_columnprod
     MPI_master_only WRITE(iscreenlog,*) 'id_zoo_columnprodtotal = ',id_zoo_columnprodtotal
#endif
#ifdef key_ulvas
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_ratiosus=',id_ulv_ratiosus
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_susmort=',id_ulv_susmort
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_suslimlight=',id_ulv_suslimlight
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_suslimN=',id_ulv_suslimN
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_suslimP=',id_ulv_suslimP
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_suspumpN=',id_ulv_suspumpN
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_suspumpP=',id_ulv_suspumpP
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_settling=',id_ulv_settling
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_ratiobenth=',id_ulv_ratiobenth
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_benthmort=',id_ulv_benthmort
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_benthlimlight=',id_ulv_benthlimlight
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_benthlimN=',id_ulv_benthlimN
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_benthlimP=',id_ulv_benthlimP
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_benthpumpN=',id_ulv_benthpumpN
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_benthpumpP=',id_ulv_benthpumpP
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_resusbil=',id_ulv_resusbil
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_resusN=', id_ulv_resusN   
     MPI_master_only WRITE(iscreenlog,*) 'id_ulv_resusP=',id_ulv_resusP  
#endif
#ifdef key_zostera
!     WRITE(iscreenlog,*) 'id_zost_LBDW  = ',id_zost_LBDW
!     WRITE(iscreenlog,*) 'id_zost_RBDW  = ',id_zost_RBDW
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_Qcan  = ',id_zost_Qcan
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_LAI  = ',id_zost_LAI
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limlum  = ',id_zost_limlum
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limprodLN  = ',id_zost_limprodLN
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limprodLP  = ',id_zost_limprodLP
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limRN  = ',id_zost_limRN
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limRP  = ',id_zost_limRP
!     WRITE(iscreenlog,*) 'id_zost_effetchaleurzprod  = ',id_zost_effetchaleurzprod
!     WRITE(iscreenlog,*) 'id_zost_effetchaleurzrespf  = ',id_zost_effetchaleurzrespf
!     WRITE(iscreenlog,*) 'id_zost_effetchaleurzrespr  = ',id_zost_effetchaleurzrespr
!     WRITE(iscreenlog,*) 'id_zost_limabsLnh4  = ',id_zost_limabsLnh4
!     WRITE(iscreenlog,*) 'id_zost_limabsLno3  = ',id_zost_limabsLno3
!     WRITE(iscreenlog,*) 'id_zost_limabsLpo4  = ',id_zost_limabsLpo4
!     WRITE(iscreenlog,*) 'id_zost_limabsRnh4  = ',id_zost_limabsRnh4
!     WRITE(iscreenlog,*) 'id_zost_limabsRpo4  = ',id_zost_limabsRpo4
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limselfshad  = ',id_zost_limselfshad
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_limRB  = ',id_zost_limRB
!     WRITE(iscreenlog,*) 'id_zost_effetvent  = ',id_zost_effetvent
!     WRITE(iscreenlog,*) 'id_zost_Sgerm  = ',id_zost_Sgerm
!     WRITE(iscreenlog,*) 'id_zost_ERS  = ',id_zost_ERS
     MPI_master_only WRITE(iscreenlog,*) 'id_zost_prod  = ',id_zost_prod
#endif
#ifdef key_oyster_SFG
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_tcoq=',id_oys_tcoq
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_tgam=',id_oys_tgam
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_gaincoq=',id_oys_gaincoq
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_gainsoma=',id_oys_gainsoma
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_gainrepr=',id_oys_gainrepr
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_orgafring=',id_oys_orgafring
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_bilansoma=',id_oys_bilansoma
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_bilancoq=',id_oys_bilancoq
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_filt=',id_oys_filt
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_resp=',id_oys_resp
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_gameto=',id_oys_gameto
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_seuilponte=',id_oys_seuilponte
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_reposcoq=',id_oys_reposcoq
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_reposgam=',id_oys_reposgam
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_ponte=',id_oys_ponte
     MPI_master_only WRITE(iscreenlog,*) 'id_oys_absorg=',id_oys_absorg                
#endif
#if defined key_oyster_DEB
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_FRpop=' ,id_oysDEB_FRpop
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_PFpop=', id_oysDEB_PFpop
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_Fpop=',id_oysDEB_Fpop
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_DWtot=',id_oysDEB_DWtot
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_kchlvble=',id_oysDEB_kchlvble
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_ER=',id_oysDEB_ER
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_filtrate=',id_oysDEB_filtrate
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_Wtot=',id_oysDEB_Wtot
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_Istress=',id_oysDEB_Istress
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_Lgtot=',id_oysDEB_Lgtot
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_ChlC=',id_oysDEB_ChlC
!     WRITE(iscreenlog,*) 'id_oysDEB_DWtot2=',id_oysDEB_DWtot2
!     WRITE(iscreenlog,*) 'id_oysDEB_Wtot2=',id_oysDEB_Wtot2
!     WRITE(iscreenlog,*) 'id_oysDEB_Lgtot2=',id_oysDEB_Lgtot2
!     WRITE(iscreenlog,*) 'id_oysDEB_DWtot3=',id_oysDEB_DWtot3
!     WRITE(iscreenlog,*) 'id_oysDEB_Wtot3=',id_oysDEB_Wtot3
!     WRITE(iscreenlog,*) 'id_oysDEB_Lgtot3=',id_oysDEB_Lgtot3
     MPI_master_only WRITE(iscreenlog,*) 'id_oysDEB_nbhuit=',id_oysDEB_nbhuit
!     WRITE(iscreenlog,*) 'id_oysDEB_nbhuit2=',id_oysDEB_nbhuit2
#endif
#ifdef key_oxygen
    MPI_master_only WRITE(iscreenlog,*) 'id_oxy_sat=',id_oxy_sat
#endif
#ifdef key_BLOOM_insed
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_aerN=',id_remin_aerN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_aerP=',id_remin_aerP
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_aerSi=',id_remin_aerSi
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_anaerN=',id_remin_anaerN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_anaerP=',id_remin_anaerP
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_nitrateN=',id_remin_nitrateN
    MPI_master_only WRITE(iscreenlog,*) 'id_remin_drnaN=',id_remin_drnaN
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
#if defined key_microtracers
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_maxdepth1 = ',id_tracer_maxdepth1
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_meandepth1 = ',id_tracer_meandepth1
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_agedepth1 = ',id_tracer_agedepth1
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_maxdepth2 = ',id_tracer_maxdepth2
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_meandepth2 = ',id_tracer_meandepth2
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_agedepth2 = ',id_tracer_agedepth2
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_maxdepth3 = ',id_tracer_maxdepth3
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_meandepth3 = ',id_tracer_meandepth3
     MPI_master_only WRITE(iscreenlog,*) 'id_tracer_agedepth3 = ',id_tracer_agedepth3
#endif
#if defined key_N_tracer
       DO is=1,nb_source_tracerN
         DO ivtra=1,nb_var_tracerN 
          MPI_master_only WRITE(iscreenlog,*) 'tracer_N source ',is,' : ',id_tracer_signN(ivtra,is) ,' : ',  &
                               TRIM(ADJUSTL(ADJUSTR(name_vardiag(id_tracer_signN(ivtra,is)))))
#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) 'age_N source ',is,' : ',id_tracer_ageN(ivtra,is) ,' : ',  &
                               TRIM(ADJUSTL(ADJUSTR(name_vardiag(id_tracer_ageN(ivtra,is)))))
#endif
         ENDDO
#if defined key_age_tracer
         idiag=id_tracer_signN(nb_var_tracerN,is)+2
#else
         idiag=id_tracer_signN(nb_var_tracerN,is)+1
#endif

#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) 'tracer_N source ',is,' : ',idiag ,' : ',  &
                         TRIM(ADJUSTL(ADJUSTR(name_vardiag(id-1))))
          idiag=idiag+1
          MPI_master_only WRITE(iscreenlog,*) 'age_N source ',is,' : ',idiag ,' : ',   &
                         TRIM(ADJUSTL(ADJUSTR(name_vardiag(id))))
#else
         MPI_master_only WRITE(iscreenlog,*) 'tracer_N source ',is,' : ',idiag ,' : ',  &
                          TRIM(ADJUSTL(ADJUSTR(name_vardiag(id))))
#endif      
      ENDDO
#endif
#if defined key_P_tracer
       DO is=1,nb_source_tracerP
         MPI_master_only WRITE(iscreenlog,*) 'tracer_P source ',is
         DO ivtra=1,nb_var_tracerP 
          MPI_master_only WRITE(iscreenlog,*) id_tracer_signP(ivtra,is) ,' : ', &
                        TRIM(ADJUSTL(ADJUSTR(name_vardiag(id_tracer_signP(ivtra,is)))))
#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) id_tracer_ageP(ivtra,is) ,' : ',  &
                      TRIM(ADJUSTL(ADJUSTR(name_vardiag(id_tracer_ageP(ivtra,is)))))
#endif
         ENDDO
#if defined key_age_tracer
         idiag=id_tracer_signP(nb_var_tracerP,is)+2
#else
         idiag=id_tracer_signP(nb_var_tracerP,is)+1
#endif

#if defined key_age_tracer
          MPI_master_only WRITE(iscreenlog,*) idiag ,' : ',TRIM(ADJUSTL(ADJUSTR(name_vardiag(id-1))))
          idiag=idiag+1
          MPI_master_only WRITE(iscreenlog,*) idiag ,' : ',TRIM(ADJUSTL(ADJUSTR(name_vardiag(id))))
#else
         MPI_master_only WRITE(iscreenlog,*) idiag ,' : ',TRIM(ADJUSTL(ADJUSTR(name_vardiag(id))))
#endif      
      ENDDO
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
#if defined key_oyster_DEB
              id_oysDEB_FRpop+ id_oysDEB_PFpop+id_oysDEB_Fpop + &
              id_oysDEB_DWtot+id_oysDEB_kchlvble+id_oysDEB_ER+ &
              id_oysDEB_filtrate+id_oysDEB_Wtot+id_oysDEB_Istress+ &
              id_oysDEB_Lgtot + id_oysDEB_ChlC + &
!              id_oysDEB_DWtot2+id_oysDEB_Wtot2+id_oysDEB_Lgtot2 + &
!              id_oysDEB_DWtot3+id_oysDEB_Wtot3+id_oysDEB_Lgtot3 + &
!              id_oysDEB_nbhuit2 + &
              id_oysDEB_nbhuit  + &
#endif

#if defined key_psnz
              id_psnz_max + id_psnz_datemax + id_psnz_depthmax + &
              id_psnz_limlight + id_psnz_limN + id_psnz_limP + id_psnz_limSi + id_psnz_columnprod + id_psnz_proddomoic + &
#endif
#if defined key_karenia
              id_karenia_max + id_karenia_datemax + id_karenia_depthmax + &
              id_karenia_limlight + id_karenia_limN + id_karenia_limP + id_karenia_columnprod + &
#endif
#if defined key_phaeocystis
              id_phaeocystis_max + id_phaeocystis_datemax + id_phaeocystis_depthmax + &
              id_phaeocystis_limlight + id_phaeocystis_limN + id_phaeocystis_limP + id_phaeocystis_columnprod + &
#endif
#ifdef key_zoo_prod
              id_zoo_micr_columnprod + id_zoo_meso_columnprod + id_zoo_columnprodtotal + &      
#endif
#ifdef key_zostera
              id_zost_Qcan + id_zost_LAI + id_zost_limlum + id_zost_limprodLN + id_zost_limprodLP + &
              id_zost_limRN + id_zost_limRP + id_zost_limselfshad + id_zost_limRB + id_zost_prod + &
!              id_zost_LBDW + id_zost_RBDW + id_zost_effetchaleurzprod + id_zost_effetchaleurzrespf + &
!              id_zost_effetchaleurzrespr + id_zost_limabsLnh4 + id_zost_limabsLno3 + id_zost_limabsLpo4 + &
!              id_zost_limabsRnh4 + id_zost_limabsRpo4 + id_zost_effetvent + id_zost_Sgerm + id_zost_ERS + &
#endif
#ifdef key_ulvas
              id_ulv_ratiosus+id_ulv_susmort+id_ulv_suslimlight+  &
              id_ulv_suslimN+id_ulv_suslimP+id_ulv_suspumpN+id_ulv_suspumpP+  &
              id_ulv_settling+id_ulv_ratiobenth+id_ulv_benthmort+  &
              id_ulv_benthlimlight+id_ulv_benthlimN+id_ulv_benthlimP+   &
              id_ulv_benthpumpN+id_ulv_benthpumpP+id_ulv_resusbil+   &
              id_ulv_resusN+id_ulv_resusP + &
#endif
#ifdef key_oyster_SFG 
              id_oys_tcoq+id_oys_tgam+id_oys_gaincoq+id_oys_gainsoma + &
              id_oys_gainrepr+id_oys_orgafring+id_oys_bilansoma+id_oys_bilancoq + &
              id_oys_filt+id_oys_resp+id_oys_gameto+id_oys_seuilponte+id_oys_reposcoq + &
              id_oys_reposgam+id_oys_ponte+id_oys_absorg + &               
#endif
#ifdef key_oxygen 
              id_oxy_sat  + &              
#endif
#if defined key_microtracers
              id_tracer_maxdepth1 + id_tracer_meandepth1 + id_tracer_agedepth1 + &
              id_tracer_maxdepth2 + id_tracer_meandepth2 + id_tracer_agedepth2 + &
              id_tracer_maxdepth3 + id_tracer_meandepth3 + id_tracer_agedepth3 + &
#endif
              id_totalchl + id_columnprodtotal + id_spm_total +  &
              id_diatsettling + id_detsettling

#if defined key_N_tracer
       DO is=1,nb_source_tracerN
         DO ivtra=1,nb_var_tracerN 
          sumverif = sumverif + id_tracer_signN(ivtra,is)
#if defined key_age_tracer
          sumverif = sumverif + id_tracer_ageN(ivtra,is) 
#endif
         ENDDO
#if defined key_age_tracer
         sumverif=sumverif+id_tracer_signN(nb_var_tracerN,is)+2+id_tracer_signN(nb_var_tracerN,is)+3
#else
         sumverif=sumverif+id_tracer_signN(nb_var_tracerN,is)+1
#endif
       ENDDO
#endif
#if defined key_P_tracer
       DO is=1,nb_source_tracerP
         DO ivtra=1,nb_var_tracerP 
          sumverif = sumverif + id_tracer_signP(ivtra,is)
#if defined key_age_tracer
          sumverif = sumverif + id_tracer_ageP(ivtra,is) 
#endif
         ENDDO
#if defined key_age_tracer
         sumverif=sumverif+id_tracer_signP(nb_var_tracerP,is)+2+id_tracer_signP(nb_var_tracerP,is)+3
#else
         sumverif=sumverif+id_tracer_signP(nb_var_tracerP,is)+1
#endif
       ENDDO
#endif
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
#if defined key_N_tracer || defined key_P_tracer
  SUBROUTINE bloom_create_vartracer_other(flx_atm_r,cv_rain_r,cini_wat_r,cini_air_r,    &
                                  l_out_subs_r,init_cv_name_r,obc_cv_name_r,            &
#if defined MUSTANG
                                  ws_free_opt_r,ws_free_para_r,ws_hind_opt_r,           &
                                  ws_hind_para_r,tocd_r,diam_r,ros_r,cini_sed_r         &
#endif
                                  !itypv_r,name_var,long_name_var,standard_name_var,     &
                                  itypv_r,long_name_var,                                &
                                  name_varpc_assoc,unit_var,l_out_subs_fix,             &
                                  long_name_var_fix,unit_var_fix,          &
                                  ws_free_min_r,ws_free_max_r)


   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_create_vartracer_other ***
   !&E
   !&E ** Purpose : create variables tracer 
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
   REAL(KIND=rsh), DIMENSION(NBVARADV_TOT),INTENT(inout)         :: flx_atm_r
   REAL(KIND=rsh), DIMENSION(NBVARADV_TOT),INTENT(inout)         :: cini_wat_r,cv_rain_r,cini_air_r
   CHARACTER(LEN=lchain), DIMENSION(NBVARADV_TOT),INTENT(inout)  :: obc_cv_name_r,init_cv_name_r
   CHARACTER(LEN=lchain),DIMENSION(ntfix),INTENT(inout)          :: long_name_var_fix,unit_var_fix
   LOGICAL, DIMENSION(NBVARADV_TOT),INTENT(inout)                :: l_out_subs_r
   LOGICAL,DIMENSION(ntfix),INTENT(inout)                        :: l_out_subs_fix
   REAL(KIND=rsh), DIMENSION(NBVARADV_TOT),intent(inout)         :: ws_free_min_r,ws_free_max_r
#if defined MUSTANG
   REAL(KIND=rsh), DIMENSION(2,NBVARADV_TOT),intent(inout)       :: ws_hind_para_r
   INTEGER       , DIMENSION(NBVARADV_TOT),intent(inout)         :: ws_hind_opt_r,ws_free_opt_r
   REAL(KIND=rsh), DIMENSION(4,NBVARADV_TOT),intent(inout)       :: ws_free_para_r
   REAL(KIND=rsh), DIMENSION(NBVARADV_TOT)                       :: diam_r,ros_r,tocd_r,cini_sed_r
#endif
   INTEGER,DIMENSION(NBVARADV_TOT),intent(inout)                 :: itypv_r
   CHARACTER(LEN=lchain),DIMENSION(NBVARADV_TOT),intent(inout)   :: name_varpc_assoc,unit_var
   !CHARACTER(LEN=lchain),DIMENSION(NBVARADV_TOT),intent(inout)   :: name_var,long_name_var,standard_name_var
   CHARACTER(LEN=lchain),DIMENSION(NBVARADV_TOT),intent(inout)   :: long_name_var 



   !! * Local declarations
    INTEGER                  :: nn,eof,ind_white,is,js,ks,ivtra,ii,nb_var_tracer,IERR_MPI
    INTEGER                  :: iss,isubs,isubs_age,iv,ivv,nvfix0
    CHARACTER(LEN=5)         :: comment
    CHARACTER(LEN=lchain), DIMENSION(NBVARADV_TOT) :: name_var_tracer
    CHARACTER(LEN=lchain)    :: name_ss
    LOGICAL                  :: ex
#if defined key_P_tracer
    CHARACTER(LEN=lchain)    :: name_varP,name_varP2
    INTEGER                  :: lenname,pos
#endif
   !!----------------------------------------------------------------------
   !! * Executable part
     nb_var_tracer=0
#if defined key_N_tracer
! definiton des variables a tracer pour l azote (8 dans le modele de base)
     nb_var_tracer=nb_var_tracer+8
     name_var_tracer(1)='ammonium'
     name_var_tracer(2)='nitrate'
     name_var_tracer(3)='nanopicoplankton_nitrogen'
     name_var_tracer(4)='diatom_nitrogen'
     name_var_tracer(5)='dinoflagellate_nitrogen'
     name_var_tracer(6)='microzooplankton_nitrogen'
     name_var_tracer(7)='mesozooplankton_nitrogen'
     name_var_tracer(8)='detrital_nitrogen'
     
#if defined key_psnz
     nb_var_tracer=nb_var_tracer+1
     name_var_tracer(nb_var_tracer)='pseudonitzschia_nitrogen'
#endif
#if defined key_karenia
     nb_var_tracer=nb_var_tracer+1
     name_var_tracer(nb_var_tracer)='karenia_nitrogen'
#endif
#if defined key_phaeocystis
     name_var_tracer(nb_var_tracer+1)='colonial_phaeocystis_nitrogen'
     name_var_tracer(nb_var_tracer+2)='phaeocystis_cells_nitrogen'
     nb_var_tracer=nb_var_tracer+2
#endif
#if defined key_benthos
     nb_var_tracer=nb_var_tracer+1
     name_var_tracer(nb_var_tracer)='organic_nitrogen_benth'
#endif
#if defined key_ulvas
!!!!!!! ATTENTION VERIFIER LES NOMS !!!!!!!!!!!!!!!!!
     name_var_tracer(nb_var_tracer+1)='ulves_susp'
     name_var_tracer(nb_var_tracer+2)='ulves_benth'
     nb_var_tracer=nb_var_tracer+2
#endif
#endif
#if defined key_P_tracer
! definiton des variables a tracer pour le phosphore (8 dans le modele de base)

     name_var_tracer(nb_var_tracer+1)='dissolved_phosphate'
     name_var_tracer(nb_var_tracer+2)='sorbed_phosphate'
     name_var_tracer(nb_var_tracer+3)='nanopicoplankton_nitrogen'
     name_var_tracer(nb_var_tracer+4)='diatom_nitrogen'
     name_var_tracer(nb_var_tracer+5)='dinoflagellate_nitrogen'
     name_var_tracer(nb_var_tracer+6)='microzooplankton_nitrogen'
     name_var_tracer(nb_var_tracer+7)='mesozooplankton_nitrogen'
     name_var_tracer(nb_var_tracer+8)='detrital_phosphorus'
     nb_var_tracer=nb_var_tracer+8
#if defined key_psnz
     nb_var_tracer=nb_var_tracer+1
     name_var_tracer(nb_var_tracer)='pseudonitzschia_nitrogen'
#endif
#if defined key_karenia
     nb_var_tracer=nb_var_tracer+1
     name_var_tracer(nb_var_tracer)='karenia_phosphorus'
#endif
#if defined key_phaeocystis
     name_var_tracer(nb_var_tracer+1)='colonial_phaeocystis_nitrogen'
     name_var_tracer(nb_var_tracer+2)='phaeocystis_cells_nitrogen'
     nb_var_tracer=nb_var_tracer+2
#endif
#if defined key_benthos
     nb_var_tracer=nb_var_tracer+1
     name_var_tracer(nb_var_tracer)='detrital_phosphorus_benth'
#endif
#endif

     call bloom_alloc_vartracer
     
! advected variables first
     isubs=nv_adv   ! tableaux separes

#if defined key_N_tracer
     DO is=1,nb_source_tracerN
       name_ss=name_source_tracerN(is)
       name_ss=name_ss(1:5)
       DO ivtra=1,nb_var_tracerN
        DO iv=1,nv_adv
         IF(TRIM(ADJUSTL(ADJUSTR(name_var(iv)))) == name_var_tracer(ivtra)) THEN
            isubs_signed_N(ivtra,is)=iv
            isubs=isubs+1
            isubs_tracer_N(ivtra,is)=isubs
#if defined key_age_tracer
            isubs_age=isubs+1
            isubs_age_N(ivtra,is)=isubs_age
            name_var(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_var(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            long_name_var(isubs_age)=TRIM(ADJUSTL(ADJUSTR(long_name_var(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            standard_name_var(isubs_age)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            unit_var(isubs_age)='days'
#else
            isubs_age=isubs
#endif
            name_var(isubs)=TRIM(ADJUSTL(ADJUSTR(name_var(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            long_name_var(isubs)=TRIM(ADJUSTL(ADJUSTR(long_name_var(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            standard_name_var(isubs)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            unit_var(isubs)=unit_var(iv)
            itypv_r(isubs:isubs_age)=itypv_r(iv)
            IF (itypv_r(iv)< 5) THEN
               ws_free_min_r(isubs:isubs_age)=ws_free_min_r(iv)
               ws_free_max_r(isubs:isubs_age)=ws_free_max_r(iv)
#ifdef MUSTANG
               DO ivv=isubs,isubs_age
                 ws_free_para_r(:,ivv)=ws_free_para_r(:,iv)
                 ws_hind_para_r(:,ivv)=ws_hind_para_r(:,iv)
               ENDDO
               ws_free_opt_r(isubs:isubs_age)=ws_free_opt_r(iv)
               ws_hind_opt_r(isubs:isubs_age)=ws_hind_opt_r(iv)
               tocd_r(isubs:isubs_age)= tocd_r(iv)     
               ros_r(isubs:isubs_age)= ros_r(iv)      
               diam_r(isubs:isubs_age)=diam_r(iv)
#endif
            ELSE IF(itypv_r(iv)== 5) THEN
               name_varpc_assoc(isubs:isubs_age)= name_varpc_assoc(iv)
            END IF
            flx_atm_r(isubs:isubs_age)=0.0_rsh
            cv_rain_r(isubs:isubs_age)=0.0_rsh
            cini_wat_r(isubs:isubs_age)=0.0_rsh
#ifdef MUSTANG
            cini_sed_r(isubs:isubs_age)=0.0_rsh
#endif
            cini_air_r(isubs:isubs_age)=0.0_rsh
            l_out_subs_r(isubs:isubs_age)=.false.
            init_cv_name_r(isubs:isubs_age)='none'
            obc_cv_name_r(isubs:isubs_age)='none'
            IF (itypv_r(isubs)==1) THEN
             nv_grav=nv_grav+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==2) THEN
             nv_sand=nv_sand+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==3) THEN
             nv_mud=nv_mud+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==4 ) THEN
             nv_ncp=nv_ncp+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==5) THEN
             nv_sorb=nv_sorb+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==6) THEN
             nv_dis=nv_dis+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==7) THEN
             nv_fix=nv_fix+1+(isubs_age-isubs)

            END IF
#if defined key_age_tracer
            isubs=isubs+1
#endif

            exit
         ENDIF
       ENDDO 
       IF(isubs_signed_N(ivtra,is) .NE. 0) THEN
! definition des p_marque (seulement parmi les 8 variables nitrogen du modele de base)
        SELECT CASE( TRIM(ADJUSTL(ADJUSTR(name_var(isubs_signed_N(ivtra,is))))) )
          CASE ('ammonium') 
           p_marque_tracerN(ivtra,is)=p_marqueNH4
          CASE('nitrate') 
           p_marque_tracerN(ivtra,is)=p_marqueNO3 
          CASE('nanopicoplankton_nitrogen') 
           p_marque_tracerN(ivtra,is)=p_marquenanoN 
          CASE('diatom_nitrogen') 
           p_marque_tracerN(ivtra,is)=p_marquediatN 
          CASE('dinoflagellate_nitrogen') 
           p_marque_tracerN(ivtra,is)=p_marquedinoN 
          CASE('microzooplankton_nitrogen') 
           p_marque_tracerN(ivtra,is)=p_marquemicrN
          CASE('mesozooplankton_nitrogen') 
           p_marque_tracerN(ivtra,is)=p_marquemesoN 
          CASE('detrital_nitrogen') 
           p_marque_tracerN(ivtra,is)=p_marquedetN 
          CASE DEFAULT
           p_marque_tracerN(ivtra,is)=0
          END SELECT
         ENDIF
      ENDDO
     ENDDO
#endif
#if defined key_P_tracer
     DO is=1,nb_source_tracerP
      name_ss=name_source_tracerP(is)
      name_ss=name_ss(1:5)
       DO ivtra=1,nb_var_tracerP
        DO iv=1,nv_adv
         IF(TRIM(ADJUSTL(ADJUSTR(name_var(iv)))) == name_var_tracer(ivtra)) THEN
            isubs_signed_P(ivtra,is)=iv
            isubs=isubs+1
            isubs_tracer_P(ivtra,is)=isubs
            name_varP=name_var(iv)
            pos=index(name_varP,"nitrogen")
            lenname=len(name_varP)
            IF (pos==0) then
              name_varP2=name_varP
            ELSE
              name_varP2=name_varP(1:pos-1)//'phosphorus'//name_varP(pos+8:lenname)
            ENDIF
#if defined key_age_tracer
            isubs_age=isubs+1
            isubs_age_P(ivtra,is)=isubs_age
            name_var(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            long_name_var(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            standard_name_var(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            unit_var(isubs_age)='days'
#else
            isubs_age=isubs
#endif
            name_var(isubs)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            long_name_var(isubs)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            standard_name_var(isubs)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            unit_var(isubs)=unit_var(iv)
            itypv_r(isubs:isubs_age)=itypv_r(iv)
            IF (itypv_r(iv)< 5) THEN
               ws_free_min_r(isubs:isubs_age)=ws_free_min_r(iv)
               ws_free_max_r(isubs:isubs_age)=ws_free_max_r(iv)
#ifdef MUSTANG
               DO ivv=isubs,isubs_age
                 ws_free_para_r(:,ivv)=ws_free_para_r(:,iv)
                 ws_hind_para_r(:,ivv)=ws_hind_para_r(:,iv)
               ENDDO
               ws_free_opt_r(isubs:isubs_age)=ws_free_opt_r(iv)
               ws_hind_opt_r(isubs:isubs_age)=ws_hind_opt_r(iv)
               tocd_r(isubs:isubs_age)= tocd_r(iv)     
               ros_r(isubs:isubs_age)= ros_r(iv)      
               diam_r(isubs:isubs_age)=diam_r(iv)
#endif
            ELSE IF(itypv_r(iv)== 5) THEN
               name_varpc_assoc(isubs:isubs_age)= name_varpc_assoc(iv)
            END IF
            flx_atm_r(isubs:isubs_age)=0.0_rsh
            cv_rain_r(isubs:isubs_age)=0.0_rsh
            cini_wat_r(isubs:isubs_age)=0.0_rsh
#ifdef MUSTANG
            cini_sed_r(isubs:isubs_age)=0.0_rsh
#endif
            cini_air_r(isubs:isubs_age)=0.0_rsh
            l_out_subs_r(isubs:isubs_age)=.false.
            init_cv_name_r(isubs:isubs_age)='none'
            obc_cv_name_r(isubs:isubs_age)='none'
            IF (itypv_r(isubs)==1) THEN
             nv_grav=nv_grav+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==2) THEN
             nv_sand=nv_sand+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==3) THEN
             nv_mud=nv_mud+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==4 ) THEN
             nv_ncp=nv_ncp+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==5) THEN
             nv_sorb=nv_sorb+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==6) THEN
             nv_dis=nv_dis+1+(isubs_age-isubs)
            ELSE IF (itypv_r(isubs)==7) THEN
             nv_fix=nv_fix+1+(isubs_age-isubs)

            END IF
#if defined key_age_tracer
            isubs=isubs+1
#endif

            exit
         ENDIF
       ENDDO 
! definition des p_marque (seulement parmi les 3 variables phosphore du modele de base)
       IF(isubs_signed_P(ivtra,is) .NE. 0) THEN
        SELECT CASE( TRIM(ADJUSTL(ADJUSTR(name_var(isubs_signed_P(ivtra,is))))) )
         CASE ('dissolved_phosphate') 
           p_marque_tracerP(ivtra,is)=p_marquePO4
         CASE('sorbed_phosphate') 
           p_marque_tracerP(ivtra,is)=p_marquePads 
         CASE('nanopicoplankton_nitrogen') 
           p_marque_tracerP(ivtra,is)=p_marquenanoP 
         CASE('diatom_nitrogen') 
           p_marque_tracerP(ivtra,is)=p_marquediatP 
         CASE('dinoflagellate_nitrogen') 
           p_marque_tracerP(ivtra,is)=p_marquedinoP 
         CASE('microzooplankton_nitrogen') 
           p_marque_tracerP(ivtra,is)=p_marquemicrP
         CASE('mesozooplankton_nitrogen') 
           p_marque_tracerP(ivtra,is)=p_marquemesoP 
         CASE('detrital_phosphorus') 
           p_marque_tracerP(ivtra,is)=p_marquedetP 
         CASE DEFAULT
           p_marque_tracerP(ivtra,is)=0
         END SELECT
       ENDIF
      ENDDO
     ENDDO
#endif

!! fixed variables
!!!!!!!!!!!!!!!!!!!!!!!
  isubs=nv_fix
  nvfix0=nv_fix
#if defined key_N_tracer
     DO is=1,nb_source_tracerN
      name_ss=name_source_tracerN(is)
      name_ss=name_ss(1:5)
       DO ivtra=1,nb_var_tracerN
        DO iv=1,nvfix0
         IF(TRIM(ADJUSTL(ADJUSTR(name_var_fix(iv)))) == name_var_tracer(ivtra)) THEN
            isubs_signed_N(ivtra,is)=iv
            isubs=isubs+1
            isubs_tracer_N(ivtra,is)=isubs
#if defined key_age_tracer
            isubs_age=isubs+1
            isubs_age_N(ivtra,is)=isubs_age
            name_var_fix(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_var_fix(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            long_name_var_fix(isubs_age)=TRIM(ADJUSTL(ADJUSTR(long_name_var_fix(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            standard_name_var_fix(isubs_age)=TRIM(ADJUSTL(ADJUSTR(standard_name_var_fix(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            unit_var_fix(isubs_age)='days'
#else
            isubs_age=isubs
#endif
            name_var_fix(isubs)=TRIM(ADJUSTL(ADJUSTR(name_var_fix(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            long_name_var_fix(isubs)=TRIM(ADJUSTL(ADJUSTR(long_name_var_fix(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            standard_name_var_fix(isubs)=TRIM(ADJUSTL(ADJUSTR(standard_name_var_fix(iv))))//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            unit_var_fix(isubs)=unit_var_fix(iv)
            !itypv_r(isubs:isubs_age)=7 ! fixed variables
            l_out_subs_fix(isubs:isubs_age)=.false.
            nv_fix=nv_fix+1+(isubs_age-isubs)
#if defined key_age_tracer
            isubs=isubs+1
#endif
            exit
         ENDIF
       ENDDO 
       p_marque_tracerN(ivtra,is)=0
      ENDDO
     ENDDO
#endif
#if defined key_P_tracer
     DO is=1,nb_source_tracerP
      name_ss=name_source_tracerP(is)
      name_ss=name_ss(1:5)
       DO ivtra=1,nb_var_tracerP
        DO iv=1,nvfix0
         IF(TRIM(ADJUSTL(ADJUSTR(name_var_fix(iv)))) == name_var_tracer(ivtra)) THEN
            isubs_signed_P(ivtra,is)=iv
            isubs=isubs+1
            isubs_tracer_P(ivtra,is)=isubs
            name_varP=name_var_fix(iv)
            pos=index(name_varP,"nitrogen")
            lenname=len(name_varP)
            IF (pos==0) then
              name_varP2=name_varP
            ELSE
              name_varP2=name_varP(1:pos-1)//'phosphorus'//name_varP(pos+8:lenname)
            ENDIF
#if defined key_age_tracer
            isubs_age=isubs+1
            isubs_age_P(ivtra,is)=isubs_age
            name_var_fix(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            long_name_var_fix(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            standard_name_var_fix(isubs_age)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_age'
            unit_var_fix(isubs_age)='days'
#else
            isubs_age=isubs
#endif
            name_var_fix(isubs)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            long_name_var_fix(isubs)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            standard_name_var_fix(isubs)=TRIM(ADJUSTL(ADJUSTR(name_varP2)))//'_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))//'_tracer'
            unit_var_fix(isubs)=unit_var_fix(iv)
            !itypv_r(isubs:isubs_age)=7
            l_out_subs_fix(isubs:isubs_age)=.false.
            cini_wat_fix(isubs:isubs_age)=cini_wat_fix(iv)
            nv_fix=nv_fix+1+(isubs_age-isubs)
#if defined key_age_tracer
            isubs=isubs+1
#endif
            exit
         ENDIF
       ENDDO 
! definition des p_marque (seulement parmi les 3 variables phosphore du modele de base)
           p_marque_tracerP(ivtra,is)=0
      ENDDO
     ENDDO
#endif

     nvpc=nv_mud+nv_sand+nv_grav
     nvp=nvpc+nv_ncp+nv_sorb
     nv_adv=nvp+nv_dis
     nv_state=nv_adv+nv_fix
     nv_tot=nv_state


     IF(nv_tot /= ntrc_substot)THEN
       IF_MPI (MASTER) THEN
         MPI_master_only WRITE(ierrorlog,*)'WARNING with the total number of variables with tracer and age variables :',nv_tot
         MPI_master_only WRITE(ierrorlog,*)'nv_tot new=',nv_tot,' is DIFFERENT from ntrc_substot',ntrc_substot
         do iv=1,nv_adv
          MPI_master_only WRITE(ierrorlog,*)'advected var : ',iv,name_var(iv)
         enddo
         do iv=1,nv_fix
          MPI_master_only WRITE(ierrorlog,*)'fixed var :', iv,name_var_fix(iv)
         enddo
         MPI_master_only WRITE(ierrorlog,*)'The simulation is stopped'
         CALL_MPI MPI_FINALIZE(IERR_MPI)
         STOP
       ENDIF_MPI
     END IF
 
  END SUBROUTINE bloom_create_vartracer_other
#endif

   !!======================================================================
#if defined key_N_tracer || defined key_P_tracer
  SUBROUTINE bloom_alloc_vartracer

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_alloc_vartracer ***
   !&E
   !&E ** Purpose : allocate variables tracer 
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : bloom_create_vartracer
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
#if defined key_N_tracer
     ALLOCATE(p_marque_tracerN(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(isubs_tracer_N(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(iv_tracer_N(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(iv_signed_N(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(isubs_signed_N(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(id_tracer_signN(nb_var_tracerN,nb_source_tracerN))
     p_marque_tracerN(:,:)=0
     isubs_signed_N(:,:)=0
     isubs_tracer_N(:,:)=0
#if defined key_age_tracer
     ALLOCATE(iv_age_N(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(isubs_age_N(nb_var_tracerN,nb_source_tracerN))
     ALLOCATE(id_tracer_ageN(nb_var_tracerN,nb_source_tracerN))
     isubs_age_N(:,:)=0
     iv_age_N(:,:)=0
#endif

       ALLOCATE(iv_zoo_meso_tra_N(nb_source_tracerN))
       ALLOCATE(iv_zoo_micr_tra_N(nb_source_tracerN))
       ALLOCATE(iv_nutr_NH4_tra_N(nb_source_tracerN))
       ALLOCATE(iv_nutr_NO3_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_nano_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_diat_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_dino_tra_N(nb_source_tracerN))
       ALLOCATE(iv_detr_tra_N(nb_source_tracerN))
#if defined key_age_tracer
       ALLOCATE(iv_zoo_meso_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_zoo_micr_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_nutr_NH4_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_nutr_NO3_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_nano_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_diat_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_dino_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_detr_age_tra_N(nb_source_tracerN))
#endif
#if defined key_benthos
       ALLOCATE(iv_benth_tra_N(nb_source_tracerN))
#if defined key_age_tracer
       ALLOCATE(iv_benth_age_tra_N(nb_source_tracerN))
#endif
#endif
#if defined key_karenia
       ALLOCATE(iv_phyto_karenia_tra_N(nb_source_tracerN))
#if defined key_age_tracer
       ALLOCATE(iv_phyto_karenia_age_tra_N(nb_source_tracerN))
#endif
#endif
#if defined key_psnz
       ALLOCATE(iv_phyto_psnz_tra_N(nb_source_tracerN))
#if defined key_age_tracer
       ALLOCATE(iv_phyto_psnz_age_tra_N(nb_source_tracerN))
#endif
#endif
#if defined key_phaeocystis
       ALLOCATE(iv_phyto_phaeocystis_colo_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_phaeocystis_cell_tra_N(nb_source_tracerN))
#if defined key_age_tracer
       ALLOCATE(iv_phyto_phaeocystis_colo_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_phyto_phaeocystis_cell_age_tra_N(nb_source_tracerN))
#endif
#endif
#if defined key_ulvas
       ALLOCATE(iv_ulv_benth_tra_N(nb_source_tracerN))
       ALLOCATE(iv_ulv_tra_N(nb_source_tracerN))
#if defined key_age_tracer
       ALLOCATE(iv_ulv_benth_age_tra_N(nb_source_tracerN))
       ALLOCATE(iv_ulv_age_tra_N(nb_source_tracerN))
#endif
#endif

! P tracer
!!!!!!!!!!!!!
#else
     ALLOCATE(p_marque_tracerP(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(isubs_tracer_P(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(iv_tracer_P(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(iv_signed_P(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(isubs_signed_P(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(id_tracer_signP(nb_var_tracerP,nb_source_tracerP))
     p_marque_tracerP(:,:)=0
     isubs_signed_P(:,:)=0
     isubs_tracer_P(:,:)=0
#if defined key_age_tracer
     ALLOCATE(iv_age_P(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(isubs_age_P(nb_var_tracerP,nb_source_tracerP))
     ALLOCATE(id_tracer_ageP(nb_var_tracerP,nb_source_tracerP))
     isubs_age_P(:,:)=0
     iv_age_P(:,:)=0
#endif
       ALLOCATE(iv_nutr_PO4_tra_P(nb_source_tracerP))
       ALLOCATE(iv_nutr_Pads_tra_P(nb_source_tracerP))
       ALLOCATE(iv_detr_tra_P(nb_source_tracerP))
       ALLOCATE(iv_zoo_meso_tra_P(nb_source_tracerP))
       ALLOCATE(iv_zoo_micr_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_nano_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_diat_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_dino_tra_P(nb_source_tracerP))
#if defined key_age_tracer
       ALLOCATE(iv_zoo_meso_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_zoo_micr_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_nutr_PO4_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_nutr_Pads_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_detr_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_nano_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_diat_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_dino_age_tra_P(nb_source_tracerP))
#endif
#if defined key_benthos
       ALLOCATE(iv_benth_tra_P(nb_source_tracerP))
#if defined key_age_tracer
       ALLOCATE(iv_benth_age_tra_P(nb_source_tracerP))
#endif
#endif
#if defined key_karenia
       ALLOCATE(iv_phyto_karenia_tra_P(nb_source_tracerP))
#if defined key_age_tracer
       ALLOCATE(iv_phyto_karenia_age_tra_P(nb_source_tracerP))
#endif
#endif
#if defined key_psnz
       ALLOCATE(iv_phyto_psnz_tra_P(nb_source_tracerP))
#if defined key_age_tracer
       ALLOCATE(iv_phyto_psnz_age_tra_P(nb_source_tracerP))
#endif
#endif
#if defined key_phaeocystis
       ALLOCATE(iv_phyto_phaeocystis_colo_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_phaeocystis_cell_tra_P(nb_source_tracerP))
#if defined key_age_tracer
       ALLOCATE(iv_phyto_phaeocystis_colo_age_tra_P(nb_source_tracerP))
       ALLOCATE(iv_phyto_phaeocystis_cell_age_tra_P(nb_source_tracerP))
#endif
#endif

#endif

 END SUBROUTINE  bloom_alloc_vartracer
#endif

   !!======================================================================

#if defined key_N_tracer || defined key_P_tracer
  SUBROUTINE bloom_create_vardiagtracer

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_create_vardiagtracer ***
   !&E
   !&E ** Purpose : create diagnostics variables tracer 
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

   !! * Local declarations
    INTEGER                  :: id,is,ivtra
    CHARACTER(LEN=lchain)    :: name_ss
   !!----------------------------------------------------------------------
   !! * Executable part
      id=ndiag_1d+ndiag_2d+ndiag_3d
#if defined key_N_tracer
       DO is=1,nb_source_tracerN
        name_ss=name_source_tracerN(is)
        name_ss=name_ss(1:5)
         DO ivtra=1,nb_var_tracerN 
            id=id+1 
            unit_vardiag(id)='unitless'
            l_diagBIOLink_out(id)=.TRUE.
            name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_N(ivtra,is)))))//'_sign'
            long_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_N(ivtra,is)))))//'_sign'
            standard_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs_tracer_N(ivtra,is)))))//'_sign'
            id_tracer_signN(ivtra,is)=id
            idimv_r(id)=3
            ndiag_3d=ndiag_3d+1
            ndiag_3d_wat=ndiag_3d_wat+1
#if defined key_age_tracer
            id=id+1
            unit_vardiag(id)='days'
            l_diagBIOLink_out(id)=.TRUE.
            name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_N(ivtra,is)))))//'_age'
            long_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_N(ivtra,is)))))//'_age'
            standard_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs_tracer_N(ivtra,is)))))//'_age'
            idimv_r(id)=3
            ndiag_3d=ndiag_3d+1
            ndiag_3d_wat=ndiag_3d_wat+1
            id_tracer_ageN(ivtra,is)=id
#endif
         ENDDO
 !!une (ou 2 pour age) variable diagnostique en plus pour phyto total
         id=id+1 
         name_vardiag(id)='phytoplankton_sign_N_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         long_name_vardiag(id)='phytoplankton_sign_N'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         standard_name_vardiag(id)='nitrogen_fraction_in_phytoplankton_from_source_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         unit_vardiag(id)='unitless'
         l_diagBIOLink_out(id)=.TRUE.
         idimv_r(id)=3
         ndiag_3d=ndiag_3d+1
         ndiag_3d_wat=ndiag_3d_wat+1
#if defined key_age_tracer
         id=id+1
         name_vardiag(id)='phytoplankton_age_N_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         long_name_vardiag(id)='phytoplankton_age_N_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         standard_name_vardiag(id)='age_of_nitrogen_fraction_in_phytoplankton_from_source_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         unit_vardiag(id)='days'
         l_diagBIOLink_out(id)=.TRUE.
         idimv_r(id)=3
         ndiag_3d=ndiag_3d+1
         ndiag_3d_wat=ndiag_3d_wat+1
#endif         
       ENDDO
! P tracer
!!!!!!!!!!
#else
       DO is=1,nb_source_tracerP
        name_ss=name_source_tracerP(is)
        name_ss=name_ss(1:5)
         DO ivtra=1,nb_var_tracerP 
            id=id+1 
            unit_vardiag(id)='unitless'
            l_diagBIOLink_out(id)=.TRUE.
            name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_P(ivtra,is)))))//'_sign'
            long_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_P(ivtra,is)))))//'_sign'
            standard_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs_tracer_P(ivtra,is)))))//'_sign'
            id_tracer_signP(ivtra,is)=id
            idimv_r(id)=3
            ndiag_3d=ndiag_3d+1
            ndiag_3d_wat=ndiag_3d_wat+1
#if defined key_age_tracer
            id=id+1
            unit_vardiag(id)='days'
            l_diagBIOLink_out(id)=.TRUE.
            name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_P(ivtra,is)))))//'_age'
            long_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(name_var(isubs_tracer_P(ivtra,is)))))//'_age'
            standard_name_vardiag(id)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs_tracer_P(ivtra,is)))))//'_age'
            idimv_r(id)=3
            ndiag_3d=ndiag_3d+1
            ndiag_3d_wat=ndiag_3d_wat+1
            id_tracer_ageP(ivtra,is)=id
#endif
         ENDDO
 !!une (ou 2 pour age) variable diagnostique en plus pour phyto total
         id=id+1 
         name_vardiag(id)='phytoplankton_sign_P_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         long_name_vardiag(id)='phytoplankton_sign_P'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         standard_name_vardiag(id)='phosphorus_fraction_in_phytoplankton_from_source_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         unit_vardiag(id)='unitless'
         l_diagBIOLink_out(id)=.TRUE.
         idimv_r(id)=3
         ndiag_3d=ndiag_3d+1
         ndiag_3d_wat=ndiag_3d_wat+1
#if defined key_age_tracer
         id=id+1
         name_vardiag(id)='phytoplankton_age_P_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         long_name_vardiag(id)='phytoplankton_age_P_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         standard_name_vardiag(id)='age_of_phorphorus_fraction_in_phytoplankton_from_source_'//TRIM(ADJUSTL(ADJUSTR(name_ss)))
         unit_vardiag(id)='days'
         l_diagBIOLink_out(id)=.TRUE.
         idimv_r(id)=3
         ndiag_3d=ndiag_3d+1
         ndiag_3d_wat=ndiag_3d_wat+1
#endif         
       ENDDO
#endif
 END SUBROUTINE  bloom_create_vardiagtracer
#endif
   !!======================================================================
  SUBROUTINE bloom_init_nut4phy3zoo2(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE bloom_init_nut4phy3zoo2  ***
   !&E
   !&E ** Purpose : Local initialisations, reading of special files
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : BIOLink_init
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2006-10  (M. Sourisseau)  Original code
   !&E       !  2008-01  (M. Sourisseau)  Introduction MES satellitale (climato)
   !&E       !  2008-01  (P. Karleskind)  Initialisation idateinf=1
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   
   IMPLICIT NONE

   !! * Arguments
   INTEGER, INTENT(IN)                                        :: ifirst,ilast,jfirst,jlast

   !! * Local declarations
   INTEGER  :: i,j,iv,k
#if defined key_oyster_benthos || defined key_oyster_DEB
   !REAL(kind=rsh)                :: surfcad,nbfilt,nbmail
   REAL(kind=rsh)                :: biom_huitre,poids_indv
   REAL(kind=rsh)                :: surface_cadastre,densite_huitre
   INTEGER                       :: ih,jh
   !INTEGER                       :: itm
   !CHARACTER(LEN=lchain)         :: fic_benthos
   !REAL(kind=rsh), DIMENSION(COMPLETE_ARRAY)        :: zlit1
#endif
#if defined key_oyster_DEB
   !CHARACTER(LEN=lchain)         :: fic_chloro
#elif defined key_oyster_DEB_GAMELAG
   REAL(kind=rsh)                :: V
#elif defined key_oyster_SFG
    INTEGER                       :: compt,compt2   ! compteurs qui servent pour tester la profondeur de la maille
#endif

   !!----------------------------------------------------------------------
   !! * Executable part

   IF_MPI (MASTER) THEN
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) '**************************************************'
     MPI_master_only WRITE(iscreenlog,*) '**************** BLOOM_USERINIT.F90 ****************'
     MPI_master_only WRITE(iscreenlog,*) '**************************************************'
   ENDIF_MPI

   ! Initializations
   !----------------
   cmes_3dmgl(:,:,:)  =0.0_rsh

   ! Modifications of global initializations
   !----------------------------------------
   !!! ATENTION , ordre des index du tableau des concentrations peut etre different
   !!!       selon le code hydro
   
   !IF(TEST_NOT_INITFROMFILE) THEN

   !END IF
   

   !!!!!!!!!!!!!!!!!!!!
   !! module BLOOM   !!
   !!!!!!!!!!!!!!!!!!!!
   
#ifdef key_benthos
   !remise a zero des vars "benthiques" quand k<>1
     iv=iv_benth_N
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_benth_Si
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_benth_P
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
#ifdef key_NPbenth
     iv=iv_benth_PO4
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_benth_NH4
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
#endif
#ifdef key_zostera
     iv=iv_benth_PO4
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_benth_NH4
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_zost_benth_N
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_zost_benth_P
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
     iv=iv_zost_benth_seed
     BENTHIC_CONCENTRATION(BENTH_INDEX_RANGE) = 0.0_rsh
#endif
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   initialisation de la couverture benthos
!   a programmer par l utilisateur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! des exemples :
! Baie de Quiberon
!   open(119,file='../../inputs/huitres_vilo.txt')
!121    read(119,*,end=120)ih,jh 
!       IF(ih.ge.limin.and.ih.le.limax.and.jh.ge.ljmin.and.jh.le.ljmax) THEN  
!         nbhuitre(ih,jh)=1.
!         hautable(ih,jh)=0.
!       ENDIF
!      goto 121
!120 close(119) 
!    surface_cadastre =  2.7781778E+07
!    DO i=limin,limax
!      DO j=ljmin,ljmax
!        IF (nbhuitre(i,j).ne.0) THEN
!              nbhuitre(i,j) = (20.e9/35.)/ surface_cadastre*CELL_SURF(i,j)
!        ENDIF
!      ENDDO
!    ENDDO

!  exemple Baie de Bougneuf
!   fic_benthos='../../inputs/BENTHOS/couvbenthos_bnuf_vilo.nc'
!   call ionc4_openr(fic_benthos)
!   WRITE(iscreenlog,*)'------------------------------------------------'
!   WRITE(iscreenlog,*)'lecture couverture huitre cultivees'
!   WRITE(iscreenlog,*)'     donnees en gPF/m2', trim(fic_bgnf)
!   call ionc4_read_xy(fic_benthos,'couverture_huitres',zlit1,imin,imax,jmin,jmax)!
!   DO j=ljmin,ljmax
!       DO i=MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
!         IF(zlit1(i,j).gt.0.) THEN
!            nbhuitre(i,j)=11.5*CELL_SURF(i,j)
!            nbhuitre2(i,j)=60.*CELL_SURF(i,j)
! !            nbhuitre3(i,j)=30.*CELL_SURF(i,j)
!            hautable(i,j)=1.
!         ENDIF
!       ENDDO
!   ENDDO
!   DO j=ljmin,ljmax
!     DO i=MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
!       IF (nbhuitre(i,j).ne.0) then
!         WRITE(*,*) 'nombre huitre init par maille i,j =' , i,j,nbhuitre(i,j)
!         WRITE(*,*) 'densitehuitre init par maille i,j =' , i,j,nbhuitre(i,j)/CELL_SURF(i,j)
!         WRITE(*,*) 'nombre huitre2 init par maille i,j =' , i,j,nbhuitre2(i,j)
!         WRITE(*,*) 'densite huitre2 init par maille i,j =' , i,j,nbhuitre2(i,j)/CELL_SURF(i,j)
!!         WRITE(*,*) 'nombre huitre3 init par maille i,j =' , i,j,nbhuitre3(i,j)
!!         WRITE(*,*) 'densite huitre3 init par maille i,j =' , i,j,nbhuitre3(i,j)/CELL_SURF(i,j)
!         WRITE(*,*) 'surface maille i,j =' , CELL_SURF(i,j)
!         WRITE(*,*) 'hauteur table conchyl maille i,j =' , i,j,hautable(i,j)
!       ENDIF
!     ENDDO
!   ENDDO

#ifdef key_oyster_DEB
    !ALLOCATE(chlorofondmoy(limin:limax,ljmin:ljmax))
    !chlorofondmoy=0.
    !ALLOCATE(tempfondmax(limin:limax,ljmin:ljmax))
    !tempfondmax=0.
!lecture de la chlorophylle moyenne pour calcul du Xk variable
!   fic_chloro='../../inputs/BENTHOS/chloro_fond_moyenne.nc'
    !fic_chloro='../../inputs/BENTHOS/chloro_temp_fond.nc'
    !call ionc4_openr(fic_chloro)

    !WRITE(iscreenlog,*)'------------------------------------------------'
    !WRITE(iscreenlog,*)'lecture de la chloro fond moyenne pour calcul de Xk'
    !WRITE(iscreenlog,*)'     donnees en microg.l-1', trim(fic_chloro)
    !call ionc4_read_xy(fic_chloro,'Chloro_moy_fond',zlit1,imin,imax,jmin,jmax)!
    !DO j=ljmin,ljmax
    !   DO i=MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
    !     IF(zlit1(i,j).ge.0.) THEN
    !       chlorofondmoy(i,j)=zlit1(i,j)
    !     ENDIF
    !  ENDDO
    !ENDDO
    !write(*,*)'chloro_fond_moy=',chlorofondmoy(78,99)
    !zlit1(:,:)=0.
    !WRITE(iscreenlog,*)'------------------------------------------------'
    !WRITE(iscreenlog,*)'lecture de la temperature fond max pour calcul de Xk'
    !WRITE(iscreenlog,*)'     donnees en degres Celsius', trim(fic_chloro)
    !call ionc4_read_xy(fic_chloro,'temp_max_fond',zlit1,imin,imax,jmin,jmax)!
    !DO j=ljmin,ljmax
    !   DO i=MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)
    !     IF(zlit1(i,j).ge.0.) then
    !       tempfondmax(i,j)=zlit1(i,j)
    !     ENDIF
    !  ENDDO
    !ENDDO
    !write(*,*)'tempfondmax=',tempfondmax(78,99)
    !call ionc4_close(fic_chloro)
#endif    
#ifdef key_oyster_SFG
    nbhuitre=0.0_rsh
    DO compt2=jfirst,jlast
       DO compt=ifirst,ilast
                IF (BATHY_H0(compt,compt2) < 30.0_rsh) THEN
                        nbhuitre(compt,compt2)=1.0_rsh
                        write(iscreenlog,*)'1 huitre ajoutee en ',compt,compt2
                ENDIF
        ENDDO
    ENDDO
#endif
#if defined key_oyster_benthos || defined key_oyster_DEB
  ! autre version
   ! surface_cadastre=0.
   ! open(110,file='../../inputs/cadastre.txt')
!112    read(110,*,end=111)ih,jh,nbhuitre(ih,jh)
   !    surface_cadastre=surface_cadastre+CELL_SURF(ih,jh)
   !    goto 112
 !! 111 close(110)
    surface_cadastre =  2.7781778E+10

!     biomasse totale huitre
    biom_huitre=3000e3
!     poids moyen d une huitre
    poids_indv=0.1
!   densite d huitre sur le secteur des parcs ind/m2
    densite_huitre=biom_huitre/poids_indv/surface_cadastre
    IF_MPI (MASTER) THEN
       MPI_master_only WRITE(iscreenlog,*)'densite_huitre=',densite_huitre
    ENDIF_MPI
!nombre d huitres par maille
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        IF(BATHY_H0(i,j).lt.25.0_rsh .AND. BATHY_H0(i,j) > -valmanq) THEN
            nbhuitre(i,j)=30
            hautable(i,j)=1.0_rsh
        ENDIF
      ENDDO
    ENDDO
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        IF(nbhuitre(i,j).gt.0.) THEN
         IF_MPI (MASTER) THEN
            MPI_master_only WRITE(iscreenlog,*)'i,j,nbhuitre=',i,j,nbhuitre(i,j)
         ENDIF_MPI
        ENDIF
      ENDDO
    ENDDO
#endif

! if WATER_CONCENTRATION has been assigned in this routine
!   CALL_MPI ex_i_rsh(-1,2,NBVARADV_TOT*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,WATER_CONCENTRATION(:,:,liminm1:limaxp2,ljminm1:ljmaxp2))
!   CALL_MPI ex_j_rsh(-1,2,NBVARADV_TOT*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,WATER_CONCENTRATION(:,:,liminm1:limaxp2,ljminm1:ljmaxp2))

#ifdef key_oyster_DEB
      write(*,*)'debut init huitre deb'
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        !IF ((nbhuitre(i,j).ne.0.0_rsh) .and. (.NOT.l_init_rtime)) THEN
        IF ((nbhuitre(i,j).ne.0.0_rsh)) THEN
         BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_res-nv_adv)=50.0_rsh
         BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_gon-nv_adv)=500.0_rsh
         BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_str-nv_adv)=310.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_gon,1,i,j)=2656.0_rsh
        ENDIF
        !IF ((nbhuitre2(i,j).ne.0.0_rsh) .and. (.NOT.l_init_rtime)) THEN
!         WATER_CONCENTRATION(iv_oysdeb_res,1,i,j)=5000.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_gon,1,i,j)=500.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_str,1,i,j)=1900.0_rsh
!       Poids initial référence
!         WATER_CONCENTRATION(iv_oysdeb_res2,1,i,j)=880.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_gon2,1,i,j)=891.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_str2,1,i,j)=5449.0_rsh
!       Poids initial faible
!         WATER_CONCENTRATION(iv_oysdeb_res2,1,i,j)=30.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_gon2,1,i,j)=0.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_str2,1,i,j)=1170.0_rsh
!       Poids initial fort
        ! WATER_CONCENTRATION(iv_oysdeb_res2,1,i,j)=465.0_rsh
        ! WATER_CONCENTRATION(iv_oysdeb_gon2,1,i,j)=271.0_rsh
        ! WATER_CONCENTRATION(iv_oysdeb_str2,1,i,j)=8214.0_rsh
        !ENDIF
!	if ((nbhuitre3(i,j).ne.0.0_rsh).and.(.NOT.l_init_rtime)) then
!         WATER_CONCENTRATION(iv_oysdeb_res3,1,i,j)=1641.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_gon3,1,i,j)=8339.0_rsh
!         WATER_CONCENTRATION(iv_oysdeb_str3,1,i,j)=10160.0_rsh
!        endif
      ENDDO
      ENDDO 
      MPI_master_only WRITE(*,*)'fin init. huitre deb'
!   Mortalites journaliere par mois et par cohortes
!       do itm=1,12
!         txmorthuitco1(itm)=2.7e-4*2.
!         txmorthuitco2(itm)=2.7e-4
!       enddo 

!    Cas faibles mortalites
!    Cohorte 1
!      txmorthuitco1(1)=1.33e-4
!      txmorthuitco1(2)=1.33e-4
!      txmorthuitco1(3)=4.1e-5
!      txmorthuitco1(4)=0.
!      txmorthuitco1(5)=2.76e-4
!      txmorthuitco1(6)=6.0e-4
!      txmorthuitco1(7)=9.07e-4
!      txmorthuitco1(8)=1.21e-3
!      txmorthuitco1(9)=1.52e-3
!      txmorthuitco1(10)=1.12e-3
!      txmorthuitco1(11)=7.17e-4
!      txmorthuitco1(12)=5.0e-4
!    Cohorte 2
!      txmorthuitco2(1)=0.
!      txmorthuitco2(2)=2.67e-4
!      txmorthuitco2(3)=2.86e-4
!      txmorthuitco2(4)=2.81e-4
!      txmorthuitco2(5)=3.5e-4
!      txmorthuitco2(6)=8.37e-4
!      txmorthuitco2(7)=1.63e-4
!      txmorthuitco2(8)=6.07e-5
!      txmorthuitco2(9)=5.43e-4
!      txmorthuitco2(10)=8.97e-5
!      txmorthuitco2(11)=2.99e-4
!      txmorthuitco2(12)=2.35e-4
!
!    Cas fortes mortalites
!    Cohorte 1
      txmorthuitco1(1)=0.
      txmorthuitco1(2)=0.
      txmorthuitco1(3)=2.36e-4
      txmorthuitco1(4)=1.11e-3
      txmorthuitco1(5)=3.63e-3
      txmorthuitco1(6)=9.17e-3
      txmorthuitco1(7)=1.39e-3
      txmorthuitco1(8)=9.83e-4
      txmorthuitco1(9)=2.06e-3
      txmorthuitco1(10)=7.53e-4
      txmorthuitco1(11)=5.27e-4
      txmorthuitco1(12)=3.50e-4
!    Cohorte 2
      !txmorthuitco2(1)=1.63e-4
      !txmorthuitco2(2)=0.
      !txmorthuitco2(3)=2.01e-4
      !txmorthuitco2(4)=1.20e-3
      !txmorthuitco2(5)=1.39e-3
      !txmorthuitco2(6)=2.46e-3
      !txmorthuitco2(7)=4.83e-4
      !txmorthuitco2(8)=1.66e-4
      !txmorthuitco2(9)=1.25e-3
      !txmorthuitco2(10)=3.73e-4
      !txmorthuitco2(11)=1.86e-4
      !txmorthuitco2(12)=6.63e-4
      
! if cv_wat has been assigned in this routine
!   CALL_MPI ex_i_rsh(-1,2,NBVARADV_TOT*kmax,liminm1,limaxp2,ljminm1,ljmaxp2,cv_wat(:,:,liminm1:limaxp2,ljminm1:ljmaxp2))

#endif 

#ifdef key_oyster_DEB_GAMELAG
write(*,*)'debut init huitre deb'
  DO j=jfirst,jlast
    DO i=ifirst,ilast
      !oyster number
      BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb-nv_adv)=30.0_rsh
    ENDDO
  ENDDO 
  DO j=jfirst,jlast
    DO i=ifirst,ilast
      IF ((BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb-nv_adv).ne.0.0_rsh)) THEN
        V=4.73
        BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_E_V-nv_adv) = V * (muV*d_v)
        BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_E-nv_adv)=14700.0_rsh
        BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_E_R-nv_adv)=4704.0_rsh
        BENTHIC_CONCENTRATION(i,j,2,iv_oysdeb_E_GO-nv_adv)=0.0_rsh      
      ENDIF
    ENDDO
  ENDDO 
#endif 
    
!#ifdef key_larve
!		i_ponte=0
!#endif

  END SUBROUTINE bloom_init_nut4phy3zoo2

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
