!---------------------------------------------------------------------------
!
!                 include    combloom
!
!---------------------------------------------------------------------------

 LOGICAL, PUBLIC     :: l_filtbenthsinus,l_SNeffect_settle, &
                        l_filtbenthmes,l_phyzoodeteffect_settle,l_ChlNratio_var


  !!!!  ---------------------------------------
   !!!!    state variable indexes
   !!!!  ---------------------------------------
   
   INTEGER, PUBLIC   :: ihour_previous
   INTEGER, PUBLIC   :: &
            iv_nutr_NO3,    &   ! index for nitrate
            iv_nutr_SiOH,   &   ! index for silicium
            iv_nutr_PO4,    &   ! index for phosphate
            iv_nutr_Pads,   &   ! index for adsorbed phosphate
            iv_phyto_diat_N, &   ! index for diatoms as nitrogen
            iv_phyto_dino_N, &   ! index for dinoflagellates as nitrogen
            iv_detr_N,      &   ! index for organic detritus as nitrogen
            iv_detr_Si,     &   ! index for organic detritus as silicon
            iv_detr_P,      &   ! index for organic detritus as phosphate
            iv_nutr_NH4,    &   ! index for ammonium
            iv_zoo_meso_N,  &   ! index for mesozooplankton
            iv_phyto_nano_N,      &   ! index for nanopicoplankton
            iv_zoo_micr_N,  &   ! index for microzooplankton
            iv_spim,        &   ! index for suspended particulate inorganic matter
            iv_phyto_diat_pp,     &   ! index for diatoms production
            iv_phyto_dino_pp,     &   ! index for dinoflagellates production
            iv_phyto_nano_pp,     &    ! index for nanoflagellates production
            ivfix_cumulprod_first,&   ! index of first fixed variable for cumulated production
            ivfix_cumulprod_last      ! index of last fixed variable for cumulated production
            
#if defined key_zoo_prod
   INTEGER, PUBLIC   :: &
            iv_zoo_micr_ps,       &   ! index for microzooplankton production
            iv_zoo_meso_ps            ! index for mesozooplankton production
#endif	    
#if defined key_BLOOM_insed
   INTEGER, PUBLIC   :: &
            iv_detrR_N       ,       &   ! index for organic detritus more refractory as nitrogen
            iv_detrR_P       ,       &   ! index for organic detritus more refractory as phosphate
            iv_ODU       ,       &   ! index for oxygen demand unit for compounds that mineralizes organic matter under anaerobic conditions (FeOH3, SO4, CH4..)
            iv_PFe                     ! index for Fer precipite
#if defined key_sedim_MPB
   INTEGER, PUBLIC   :: &
            iv_MPB_C,  &
            iv_MPB_N,  &
            iv_Bact,  &
            iv_EPS
#endif
#endif

#if defined key_BLOOM_opt2
   INTEGER, PUBLIC   :: &
!! Insertion VS juin2010
            iv_diss_N,           & ! index for dissolved nitrogen
            iv_diss_P,           & ! index for dissolved phosphate
            iv_diss_Si,          & ! index for dissolved silicate
!! fin insertion nvelles variables
! Insertion VS nov 2010
            iv_detr_fond_N,      & ! index for detritus cumulated in bottom of last wat layer
            iv_diss_fond_Nitr      ! index for diss nitr (N2) issued from det in last wat layer
#endif
#if defined key_psnz
   INTEGER, PUBLIC   :: &
            iv_phyto_psnz_N, &   ! index for PseudoNitzschia as nitrogen
            iv_phyto_psnz_Si,&   ! index for PseudoNitzschia as silicon
            iv_phyto_psnz_ad,&   ! index for PseudoNitzschia as domoic acid
            iv_phyto_psnz_pp     ! index for PseudoNitzschia primary production 
#endif
#if defined key_karenia
   INTEGER, PUBLIC   :: &
            iv_phyto_karenia_C, &   ! index for karenia as carbon
            iv_phyto_karenia_N, &   ! index for karenia as nitrogen
            iv_phyto_karenia_P, &   ! index for karenia as phosphorus
            iv_phyto_karenia_pp     ! index for karenia primary production 
#endif
#if defined key_phaeocystis
   INTEGER, PUBLIC   :: &
            iv_phyto_phaeocystis_colo_N, &   ! index for phaeocystis colonies as nitrogen
            iv_phyto_phaeocystis_cell_N, &   ! index for phaeocystis cells as nitrogen
            iv_phyto_phaeocystis_mucus,  &   ! index for phaeocystis mucus as mass
            iv_phyto_phaeocystis_pp          ! index for phaeocystis primary production 
#endif
#ifdef key_ulvas
   INTEGER, PUBLIC   :: &
            iv_ulv_sus_N,     &     ! index for suspended ulvas as nitrogen
            iv_ulv_benth_N,     &   ! index for benthic ulvas as nitrogen
            iv_ulv_susdrywght, &   ! index for suspended ulvas as dry weight
            iv_ulv_benthdrywght, & ! index for benthic ulvas as dry weight
            iv_ulv_sus_P,      &    ! index for suspended ulvas as phosphorus
            iv_ulv_benth_P,     &   ! index for benthic ulvas as phosphorus
            iv_ulv_suspp,      &   ! index for suspended ulvas production
            iv_ulv_benthpp         ! index for benthic ulvas production	    
#endif
#ifdef key_zostera
   INTEGER, PUBLIC   ::               &
                   iv_zost_LB,        &   ! index for zostera leaf biomass
                   iv_zost_RB,        &   ! index for zostera rhizomes and roots biomass
                   iv_zost_D,         &   ! index for zostera shoot density
                   iv_zost_LN,        &   ! index for zostera leaf nitrogen pool
                   iv_zost_RN,        &   ! index for zostera rhizomes and roots nitrogen pool
                   iv_zost_LP,        &   ! index for zostera leaf phophorus pool
                   iv_zost_RP,        &   ! index for zostera rhizomes and roots phosphorus pool
                   iv_detr_zost_N,    &   ! index for floating detrital N
                   iv_detr_zost_P,    &   ! index for floating detrital P
                   iv_zost_benth_N,   &   ! index for floating detrital N
                   iv_zost_benth_P,   &   ! index for floating detrital P
                   iv_zost_seed,      &   ! index for zostera seeds in the water column
                   iv_zost_benth_seed     ! index for zostera seeds in the sediment
!                   iv_zost_pp,        &   ! index for zostera cumulative production
!                   iv_zost_LNuptake,  &   ! index for zostera cumulative leaf N uptake
!                   iv_zost_LPuptake,  &   ! index for zostera cumulative leaf P uptake
!                   iv_zost_RNuptake,  &   ! index for zostera cumulative root N uptake
!                   iv_zost_RPuptake       ! index for zostera cumulative root P uptake
#endif
#ifdef key_oyster_SFG
   INTEGER, PUBLIC   :: &
            iv_oys_so,     &   ! index for oyster dry weight
            iv_oys_go,     &   ! index for oyster resgon weight
            iv_oys_co          ! index for oyster coq weight
#endif
#ifdef key_oxygen
   INTEGER, PUBLIC   :: &
            iv_oxygen          ! index for oxygen
#endif
   INTEGER, PUBLIC :: id_benthos_txf
#ifdef key_oyster_DEB
   INTEGER, PUBLIC   :: &
            iv_oysdeb_res,     &   ! index for oyster supplies
            iv_oysdeb_str,     &   ! index for oyster structures
            iv_oysdeb_gon          ! index for oyster gonades
        INTEGER, PUBLIC   :: &
            iv_oysdeb_res2,     &   ! index for oyster supplies
            iv_oysdeb_str2,     &   ! index for oyster structures
            iv_oysdeb_gon2          ! index for oyster gonades
        INTEGER, PUBLIC   :: &
            iv_oysdeb_res3,     &   ! index for oyster supplies
            iv_oysdeb_str3,     &   ! index for oyster structures
            iv_oysdeb_gon3          ! index for oyster gonades
#endif
#ifdef key_benthos
   INTEGER, PUBLIC   :: &
            iv_benth_N,     &   ! index for organic nitrogen in benthos
            iv_benth_Si,     &  ! index for organic silicon in benthos
            iv_benth_P          ! index for organic phosphorus in benthos
#ifdef key_diatbenth
   INTEGER, PUBLIC   :: iv_benth_phyto_diat_N   !index for benthic diatoms as nitrogen
#endif
#ifdef key_NPbenth
   INTEGER, PUBLIC   :: iv_benth_PO4,  iv_benth_NH4
#endif
#ifdef key_zostera
   INTEGER, PUBLIC   :: iv_benth_phyto_diat_N,  &  !index for benthic diatoms as nitrogen
                             iv_benth_PO4,           &
                             iv_benth_NH4
#endif
#endif
#ifdef key_microtracers
   INTEGER, PUBLIC   :: &
            iv_tracer_depth1, &       ! index for micro-tracer from depth 1
            iv_tracer_depth2,&        ! index for micro-tracer from depth 2
            iv_tracer_depth3,&        ! index for micro-tracer from depth 3
            iv_tracer_agedepth1,&     ! index for age of the micro-tracer from depth 1
            iv_tracer_agedepth2,&     ! index for age of the micro-tracer from depth 2
            iv_tracer_agedepth3,&     ! index for age of the micro-tracer from depth 3
            iv_tracer_cumuldepth1,&   ! index for cumulatiuve micro-tracer from depth 1
            iv_tracer_cumuldepth2,&   ! index for cumulatiuve micro-tracer from depth 2
            iv_tracer_cumuldepth3     ! index for cumulatiuve micro-tracer from depth 3
#endif
#if defined key_N_tracer
! index (rank model et rank file) for tracer variable coming from source
   INTEGER, dimension(:,:),ALLOCATABLE, PUBLIC   ::  iv_tracer_N,isubs_tracer_N 
! index for signed variable which is traced by corresponding tracer_variable (rank model et rank file)    
   INTEGER, dimension(:,:),ALLOCATABLE, PUBLIC   ::  iv_signed_N,isubs_signed_N    
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_detr_tra_N,iv_nutr_NH4_tra_N,iv_nutr_NO3_tra_N,    &
                                                 iv_phyto_nano_tra_N,iv_phyto_diat_tra_N,iv_phyto_dino_tra_N,  &
                                                 iv_zoo_micr_tra_N,iv_zoo_meso_tra_N
#if defined key_benthos
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_benth_tra_N
#endif
#if defined key_psnz
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_psnz_tra_N
#endif
#if defined key_phaeocystis
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_phaeocystis_colo_tra_N,iv_phyto_phaeocystis_cell_tra_N
#endif
#if defined key_karenia
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_karenia_tra_N
#endif
#if defined key_ulvas
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_ulv_benth_tra_N,iv_ulv_tra_N
#endif

#ifdef key_age_tracer
! index for the age (rank model et rank file)
   INTEGER, dimension(:,:),ALLOCATABLE, PUBLIC   ::  isubs_age_N,iv_age_N 
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: &
                         iv_detr_age_tra_N,iv_nutr_NH4_age_tra_N,iv_nutr_NO3_age_tra_N,iv_phyto_nano_age_tra_N, &
                         iv_phyto_diat_age_tra_N,iv_phyto_dino_age_tra_N,iv_zoo_micr_age_tra_N,iv_zoo_meso_age_tra_N
#if defined key_benthos
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_benth_age_tra_N
#endif
#if defined key_psnz
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_psnz_age_tra_N
#endif
#if defined key_phaeocystis
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_phaeocystis_colo_age_tra_N,iv_phyto_phaeocystis_cell_age_tra_N
#endif
#if defined key_karenia
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_karenia_age_tra_N
#endif
#if defined key_ulvas
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_ulv_benth_age_tra_N,iv_ulv_age_tra_N
#endif
#endif
#endif

#if defined key_P_tracer
! index (rank model et rank file) for tracer variable coming from source 
   INTEGER, dimension(:,:),ALLOCATABLE, PUBLIC   ::  iv_tracer_P,isubs_tracer_P  
! index for signed variable which is traced by corresponding tracer_variable (rank model et rank file)
   INTEGER, dimension(:,:),ALLOCATABLE, PUBLIC   ::  iv_signed_P,isubs_signed_P  
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_detr_tra_P,iv_nutr_PO4_tra_P,iv_nutr_Pads_tra_P,    &
                                                 iv_phyto_nano_tra_P,iv_phyto_diat_tra_P,iv_phyto_dino_tra_P, &
                                                 iv_zoo_micr_tra_P,iv_zoo_meso_tra_P
#if defined key_benthos
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_benth_tra_P
#endif
#if defined key_psnz
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_psnz_tra_P
#endif
#if defined key_phaeocystis
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_phaeocystis_colo_tra_P,iv_phyto_phaeocystis_cell_tra_P
#endif
#if defined key_karenia
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_karenia_tra_P
#endif
#if defined key_ulvas
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_ulv_benth_tra_P,iv_ulv_tra_P
#endif

#ifdef key_age_tracer
   INTEGER, dimension(:,:),ALLOCATABLE, PUBLIC   ::  isubs_age_P,iv_age_P ! index for the age (rank model et rank file)
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_detr_age_tra_P,iv_nutr_PO4_age_tra_P,iv_nutr_Pads_age_tra_P, &
                                                 iv_phyto_nano_age_tra_P,iv_phyto_diat_age_tra_P, &
                                                 iv_phyto_dino_age_tra_P,iv_zoo_micr_age_tra_P,iv_zoo_meso_age_tra_P
#if defined key_benthos
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_benth_age_tra_P
#endif
#if defined key_psnz
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_psnz_age_tra_P
#endif
#if defined key_phaeocystis
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_phaeocystis_colo_age_tra_P,iv_phyto_phaeocystis_cell_age_tra_P
#endif
#if defined key_karenia
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_phyto_karenia_age_tra_P
#endif
#if defined key_ulvas
   INTEGER, dimension(:),ALLOCATABLE,PUBLIC   :: iv_ulv_benth_age_tra_P,iv_ulv_age_tra_P
#endif
#endif
#endif

#ifdef key_suspensivores
   INTEGER, PUBLIC   :: &
            iv_benth_susp_N     ! index for suspensivores as nitrogen in benthos
#endif
#ifdef key_benthos_gener
   INTEGER, PUBLIC   :: &
            iv_benth_meio_N, &    ! index for bacterie+meiofaune as nitrogen in benthos
            iv_benth_deposi_N, &    ! index for deposivore+herbivore as nitrogen in benthos
            iv_benth_susp_N,   &  ! index for suspensivores as nitrogen in benthos
            iv_benth_carn_N     ! index for carnivore as nitrogen in benthos
#endif

#ifdef key_larve
         INTEGER, PUBLIC   :: &  
          iv_larve_1, &          !index for larve zone 1
          iv_larve_2, &          !index for larve zone 2
          iv_larve_3, &          !index for larve zone 3
          iv_larve_4, &          !index for larve zone 4
          iv_larve_5, &         !index for larve zone 5
          iv_larve_6, &
          i_ponte

#endif

   !!!!  ---------------------------------------
   !!!!    bio parameters
   !!!!  ---------------------------------------
   
 REAL(KIND=rsh), PUBLIC       ::     &
                   p_sali_thhold_bio,  &
                   p_N_remin,          &
                   p_P_remin,          &
                   p_N_reminR,         &
                   p_P_reminR,         &
                   p_BSi_dissEau,      &
                   p_BSi_dissSurfSed,  &
                   p_BSi_dissFondSed,  &
                   p_kSi,              &
                   p_k_remin,          &
                   p_xflimz,           &
                   p_T_effect,         &
                   p_T_effectSi,       &
                   p_burial,           &
                   p_P_speedup_reminanaer, &
                   p_P_adsor,          &
                   p_P_desor,          &
                   p_P_adsormaxspim,   &
                   p_P_adsormaxsed,    &
                   p_phyto_ChlNratio,  &
                   rappaz,rapsiaz,     &
                   p_phyto_SiNratio,   &
                   p_phyto_NPratio,    &
                   p_phyto_CNratio,    &
                   p_phyto_ChlNratiomax, &
                   p_phyto_ChlN_ksmithextinct, &
                   p_nano_mumax,       &
                   p_nano_kNO3,        &
                   p_nano_kNH4,        &
                   p_nano_kPO4,        &
                   p_nano_thhold_mort, &
                   p_nano_mort,        &
                   p_nano_iksmith,     &
                   p_diat_mumax,       &
                   p_diat_kNO3,        &
                   p_diat_kNH4,        &
                   p_diat_kSi,         &
                   p_diat_kPO4,        &
                   p_diat_iksmith,     &
!#if defined key_BLOOM_opt2
                   p_diat_mort,        &
!#else
!                   p_diat_mortmax,     &
!#endif
                   p_diat_mort_sed,    &
                   p_diat_thhold_mort, &
                   p_dino_mumax,       &
                   p_dino_thhold_ect,  &
                   p_dino_kNO3,        &
                   p_dino_kNH4,        &
                   p_dino_kPO4,        &
                   p_dino_mort,        &
                   p_dino_thhold_mort, &
                   p_dino_iksmith,     &
                   p_mesz_thrN,      &
                   p_mesz_kivlev,      &
                   p_mesz_thhold_mes_kivlev, &
                   p_mesz_thhold_mort, &
                   p_mesz_mumax,       &
                   p_mesz_assim,       &
                   p_mesz_excret,      &
                   p_mesz_mort1,       &
                   p_mesz_mort2,       &
                   p_zoo_CDWratio,     &
                   p_zoo_CNratio,      &
                   p_micz_mumax,       &
                   p_micz_kgraz,       &
                   p_micz_assim,       &
                   p_micz_excret,      &
                   p_micz_mort,        &
                   p_micz_thhold_mort, &
                   p_micz_thrnano,     &
                   p_phyto_photoratio, &
                   !p_Nbent_remin,      &
                   !p_Pbent_remin,      &
                   !p_Sibent_dissol,    &
                   p_kO2_reminO2,       &
                   p_nitrif,        &
                   p_phyto_resp,       &
                   p_zoo_resp,         &
                   p_extincwat,        &
                   p_extincspim,       &
                   p_extincChl1,       &
                   p_extincChl2,       &
                   p_parradratio,      &
                   p_erodflux,         &
                   p_erodvitcrit,      &
                   p_depovitcrit, &
                   p_erodcrittau,      &
                   p_detzoo_wsed,      &
                   p_detphy_wsed,      &
                   p_mesz_captdiat,    &
                   p_mesz_captpsnz,    & 
                   p_mesz_captdino,    &
                   p_mesz_captkarenia, &  
                   p_mesz_captmicz,    &
                   p_micz_captdiat,    &
                   p_micz_captdino,    &
                   p_micz_captnano,    &
                   p_micz_captdet,     &
                   p_micz_captkarenia, &
                   p_txfiltbenthmax,   &
                   p_kO2_nit,          &
                   p_kO2_reoxyd,       &
                   p_kNO3_reminssO2,   &
                   p_kiO2_remin0O2,    &
                   p_kiO2_denit,       &
                   p_kiNO3_remin0O2,   &
                   p_ODU_oxy,         &
                   p_ODU_precip,      &
                   p_DNO3_denit,       &
                   p_GO2_Norg,         &
                   p_GODU_Norg,        &
                   p_GNO3_Norg,        &
                   p_GO2_NorgR,         &
                   p_GODU_NorgR,        &
                   p_GNO3_NorgR,        &
                   p_GO2_NH4,          &    
                   p_kO2_precPFe,      &
                   p_kNO3_precPFe,     &
                   p_kiO2_dissPFe,     &
                   p_kiO2_desorP,      &
                   p_kiNO3_dissPFe,    &
                   p_P_precFeO2,        &
                   p_P_precFeNO3,       &
                   p_P_dissFe,        &
                   p_GO2_PFe,          &
                   p_GNO3_PFe,         &
                   p_Si_precip,         &
                   p_aging_MO,         &
                   p_burried,         &
                   p_KO2sed_aeration,    &
                   p_Kzsed_aeration,     &
!                   p_kiO2_precSi
                   p_Si_Eq,              &
                   p_Si_EqPrec


#if defined key_BLOOM_opt2
!! Insertion nveaux parametres en rapport avec ajout variables dissoutes VS juin 2010
 REAL(KIND=rsh), PUBLIC       ::       &
                   p_det_fragm,        &
                   p_diss_regmod,      & 
                   p_micz_diss
#else
 REAL(KIND=rsh), PUBLIC       ::   p_reminbenth
#endif

#if defined key_BLOOM_insed && defined key_sedim_MPB
!    parametres des reactions dans le sediment
 REAL(KIND=rsh), PUBLIC   :: p_Tmax_MPB,p_Topt_MPB,p_beta_temp,p_K_PARsed,p_KuptC_MPBN,p_gamma_P,  &
                             p_EPS_leaching,p_uptN_alpha2,p_KuptN_DIN,p_KuptN_MPBC,p_migexu_alpha1,   &
                             p_mig_delta_NC,p_prop_resp_MPBC,p_Tmax_Bact,p_Topt_Bact,p_uptN_Bact_alpha3, &
                             p_KuptBactN_DON,p_ratio_uptakeCsN,p_uptBactN_EPS_beta1,p_morta_MPB,ratio_uptakeNO3sNH4,  &
                             p_morta_Bact,p_qN2C_migupMPB,p_qC2N_migdwnMPB,p_Kext_sed,p_ratioC_Chla, &
                             p_ratio_nu_ud,p_nu_up,p_GO2_respMPB,attenu_w,xK_migd_N,ratioC_Chla,kd_chla,QR_diat_inv,  &
                             invMPBCmaxm3,MPBCmaxm3,seuil_NH4_uptakNH4only,seuil_NH4_uptakNH4NO3,Chlamax_m2microm,QR_diat
#endif

#if defined key_psnz
   REAL(KIND=rsh), PUBLIC       ::     &
                   p_psnz_mumax,       &
                   p_psnz_tempopt,     &
                   p_psnz_templethal,  &   
                   p_psnz_tempsigma,   &
                   p_psnz_beta,        &
                   p_psnz_kNO3,        &
                   p_psnz_kNH4,        &
                   p_psnz_kSi,         &
                   p_psnz_kPO4,        &
                   p_psnz_vmaxSi,      &
                   p_psnz_qminSi,      &
                   p_psnz_kqSi,        &
                   p_psnz_iksmith,     &
                  ! p_psnz_mortmax,     &
                   p_psnz_mort,        &
                   p_psnz_thhold_Si,   &
                   p_psnz_prod_domoic, &
                   p_psnz_decay_domoic
#endif
#if defined key_karenia
   REAL(KIND=rsh), PUBLIC       ::     &
                   p_karenia_mumax,       &
                   p_karenia_thhold_ect,  &
                   p_karenia_tempopt,     &
                   p_karenia_tempsigma,   &
                   p_karenia_templethal,  &
                   p_karenia_beta,        &
                   p_karenia_vmaxN,       &
                   p_karenia_kqN,         &
                   p_karenia_qminN,       &
                   p_karenia_vmaxP,       &
                   p_karenia_kqP,         &
                   p_karenia_qminP,       &
                   p_karenia_kNO3,        &
                   p_karenia_kNH4,        &
                   p_karenia_kPO4,        &
                   p_karenia_iksmith,     &
                   p_karenia_mort,        &
                   p_karenia_wmax
#endif
#if defined key_phaeocystis
   REAL(KIND=rsh), PUBLIC       ::     &
                   p_phaeocystis_mumaxcolo,   &
                   p_phaeocystis_mumaxcell,   &
                   p_phaeocystis_tempopt,     &
                   p_phaeocystis_templethal,  &   
                   p_phaeocystis_tempsigma,   &
                   p_phaeocystis_beta,        &
                   p_phaeocystis_kNO3colo,    &
                   p_phaeocystis_kNO3cell,    &
                   p_phaeocystis_kNH4,        &
                   p_phaeocystis_kPO4colo,    &
                   p_phaeocystis_kPO4cell,    &
                   p_phaeocystis_iksmith,     &
                   p_phaeocystis_mort,        &
                   p_phaeocystis_lyse,        &
                   p_phaeocystis_coloinit,    &
                   p_phaeocystis_coloseuil,   &
                   p_phaeocystis_mucus_decay, &
                   p_mesz_captphaeocolo
#endif
#ifdef key_ulvas
   REAL(KIND=rsh), PUBLIC       ::     &
                   p_ulv_mumax,        &                      
                   p_ulv_iksmith,      &     
                   p_ulv_kNO3,         &        
                   p_ulv_kPO4,         &       
                   p_ulv_maxabsN,      &               
                   p_ulv_maxabsP,      &              
                   p_ulv_minNWratio,   &     
                   p_ulv_maxNWratio,   &      
                   p_ulv_minPWratio,   &      
                   p_ulv_maxPWratio,   &      
                   p_ulv_dryfreshratio, &      
                   p_ulv_susmort,      &      
                   p_ulv_depmort,      &       
                   p_ulv_extincwat,    &     
                   p_ulv_erodflux,     &    
                   p_ulv_erodcritv,    &        
                   p_ulv_settlcrittv 
#endif
#ifdef key_zostera
  REAL(KIND=rsh), PUBLIC       ::      &
                   p_zost_klai,        &
                   p_zost_leafabscoef, &
                   p_zost_Ik,          &
                   p_zost_maxprod0,    &
                   p_zost_T_prod,      &
                   p_zost_respf0,      &
                   p_zost_respr0,      &
                   p_zost_T_respf,     &
                   p_zost_T_respr,     &
                   p_zost_qphotos,     &
                   p_zost_LNquotamin,  &
                   p_zost_LNquotamax,  &
                   p_zost_RNquotamin,  &
                   p_zost_RNquotamax,  &
                   p_zost_KNH4L,       &
                   p_zost_KNH4R,       &
                   p_zost_KNO3L,       &
                   p_zost_delta1,      &
                   p_zost_delta2,      &
                   p_zost_LNnh4maxabs, &
                   p_zost_LNno3maxabs, &
                   p_zost_RNmaxabs,    &
                   p_zost_LPquotamin,  &
                   p_zost_LPquotamax,  &
                   p_zost_RPquotamin,  &
                   p_zost_RPquotamax,  &
                   p_zost_KPL,         &
                   p_zost_KPR,         &
                   p_zost_LPmaxabs,    &
                   p_zost_RPmaxabs,    &
                   p_zost_K,           &
                   p_zost_transfrate,  &
                   p_zost_reclamax,    &
                   p_zost_mort,        &
                   p_zost_T_mort,      &
                   p_zost_RECRmax,     &
                   p_zost_T_recr,      &
                   p_zost_KREC1,       &
                   p_zost_KREC2,       &
                   p_zost_SB1,         & 
                   p_zost_KGERM,       &
                   p_zost_SB0,         &
                   p_zost_ERSmax,      &
                   p_zost_seedprod,    &
                   p_zost_GERmax,      &          
                   p_zost_Smort
#endif
#ifdef key_oyster_benthos
           REAL(KIND=rsh), PUBLIC ::     &
                   p_huitre_DW,          &
                   p_huitre_thr_colmat,  &
                   p_huitre_temp1,       &
                   p_huitre_temp_opt,    &
                   p_huitre_thr_mes,     &
                   p_huitre_loi_mes1,    &
                   p_huitre_loi_mes2,    &
                   p_huitre_filt_std,    &
                   p_huitre_exp_allom,   &
                   p_huitre_loi_colmat
#endif
#ifdef key_oyster_SFG
           REAL(KIND=rsh), PUBLIC       ::     &
                    p_oys_absdet,              &
                    p_oys_absmaxmop,           &
                    p_oys_kreten,              &
                    p_oys_allfiltr,            &
                    p_oys_allres,              &
                    p_oys_yreten,              &
                    p_oys_CtoDWratio,          &
                    p_oys_stdfiltr,            & 
                    p_oys_paramfiltr,          &
                    p_oys_tempopt_filtr,       &
                    p_oys_sestfiltr_thr,       &
                    p_oys_yfilt,               &
                    p_oys_kfilt,               &
                    p_oys_colmparam,           &
                    p_oys_pseufeci_exp,        &
                    p_oys_pseufeco_exp,        &
                    p_oys_prodpseufeci_lev,    &
                    p_oys_prodpseufeco_lev,    &
                    p_oys_retmin,              &
                    p_oys_retmin_thr,          &
                    p_oys_prodpseufec_thr,     &
                    p_oys_subratio,            &
                    p_oys_colmat_thr
#endif
#ifdef key_oyster_DEB
           REAL(KIND=rsh), PUBLIC       ::     &
                   p_oysDEB_TA,                &
                   p_oysDEB_TAL,               &
                   p_oysDEB_TAH,               &
                   p_oysDEB_TL,                &
                   p_oysDEB_THing,             &
                   p_oysDEB_THresp,            &
                   p_oysDEB_Tseuilponte,       &                
                   p_oysDEB_Eg,                &
                   p_oysDEB_Em,                &
                   p_oysDEB_Kappa,             &    
                   p_oysDEB_Vp,                &
                   p_oysDEB_kR,                &
                   p_oysDEB_shape,             &
                   p_oysDEB_pxm0,              &
                   p_oysDEB_pammax,            &
                   p_oysDEB_ae,                &
                   p_oysDEB_pm0,               &
                   p_oysDEB_ERlim,             &
                   p_oysDEB_kchl,              &
                   p_oysDEB_muE 
#endif	  
#ifdef key_microtracers
   REAL(KIND=rsh), PUBLIC       ::      &
                   p_trace_debitinject, &
                   p_trace_depth1inject,&
                   p_trace_depth2inject,&
                   p_trace_depth3inject
#endif
#if defined key_N_tracer
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE, PUBLIC       :: p_marque_tracerN
   INTEGER,PUBLIC                                          :: nb_source_tracerN,nb_source_river_tracerN,nb_source_marin_tracerN
   CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE, PUBLIC  :: name_source_tracerN
   REAL(KIND=rsh)                                          :: p_marqueNH4,p_marqueNO3,p_marquenanoN,p_marquedinoN,p_marquediatN,  &
                                                              p_marquemicrN,p_marquemesoN,p_marquedetN
   REAL(KIND=rsh),PUBLIC                                   :: tdeb_tracerN
   CHARACTER(LEN=lchain)                                   :: p_source_river1_tracerN,p_source_river2_tracerN,   &
                                                              p_source_river3_tracerN,p_source_river4_tracerN
#endif
#if defined key_P_tracer
   REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE, PUBLIC       :: p_marque_tracerP
   INTEGER,PUBLIC                                          :: nb_source_tracerP,nb_source_river_tracerP,nb_source_marin_tracerP
   CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE, PUBLIC  :: name_source_tracerP
   REAL(KIND=rsh)                                          :: p_marquePO4,p_marquePads,p_marquedetP,p_marquenanoP,  &
                                                              p_marquedinoP,p_marquediatP,p_marquemicrP,p_marquemesoP
   REAL(KIND=rsh),PUBLIC                                   :: tdeb_tracerP
   CHARACTER(LEN=lchain)                                   :: p_source_river1_tracerP,p_source_river2_tracerP,  &
                                                              p_source_river3_tracerP,p_source_river4_tracerP
#endif

REAL(KIND=rsh), PUBLIC       ::     &
                   p_susp_satur_proie, &
                   p_susp_capdet,      &
                   p_susp_thrdet,      &
                   p_susp_thrphy,      &
                   p_susp_capphy,      &
                   p_susp_mort,        &
                   p_susp_pphysio,     &
                   p_susp_egest_det,   &
                   p_susp_egest_phy,   &
                   p_susp_thrmax_inhib,&
                   p_susp_inges,       &
                   p_susp_thrmin_inhib,&
                   p_meio_satur_proie, &
                   p_meio_capdet,      &
                   p_meio_thrdet,      &
                   p_meio_capdiat,     &
                   p_meio_thrdiat,     &
                   p_meio_mort,        &
                   p_meio_pphysio,     &
                   p_meio_egest_det,   &
                   p_meio_thrmax_inhib,&
                   p_meio_thrmin_inhib,&
                   p_meio_inges,       &
                   p_deposi_satur_proie,&
                   p_deposi_capdet,    &
                   p_deposi_thrdet,    &
                   p_deposi_capdiat,   &
                   p_deposi_capmeio,   &
                   p_deposi_thrmeio,   &
                   p_deposi_mort,      &
                   p_deposi_pphysio,   &
                   p_deposi_egest_meio,&
                   p_deposi_egest_det, &
                   p_deposi_egest_diat,&
                   p_deposi_thrmax_inhib,&
                   p_deposi_thrmin_inhib,&
                   p_deposi_inges,     &
                   p_carn_satur_proie, &
                   p_carn_capmeio,     &
                   p_carn_thrmeio,     &
                   p_carn_capdeposi,   &
                   p_carn_thrdeposi,   &
                   p_carn_satur_proie_susp,&
                   p_carn_capsusp,     &
                   p_carn_thrsusp,     &
                   p_carn_mort,        &
                   p_carn_pphysio,     &
                   p_carn_thrmax_inhib,&
                   p_carn_thrmin_inhib,&
                   p_carn_inges,       &
                   p_carn_egest

   ! indexes of diagnostic variables
   INTEGER, PUBLIC :: id_diat_max,          &
                      id_diat_datemax,      &
                      id_diat_depthmax,     &
                      id_dino_max,          &
                      id_dino_datemax,      &
                      id_dino_depthmax,     &
                      id_nano_max,          &
                      id_nano_datemax,      &
                      id_nano_depthmax,     &
                      id_gradsali_max,      &
                      id_gradsali_depthmax, &
                      id_gradtemp_max,      &
                      id_graddens_max,      &
                      id_gradtemp_depthmax, &
                      id_graddens_depthmax, &
#ifdef key_BLOOM_insed
                      id_diffuflux_NO3,     &
                      id_diffuflux_PO4,     &
                      id_diffuflux_NH4,     &
                      id_diffuflux_O2D,     &
#endif
#ifdef key_BLOOM_opt2
                      id_spim_satused,      &
#endif
                      id_diat_limlight,     &
                      id_diat_limN,         &
                      id_diat_limSi,        &
                      id_diat_limP,         &
                      id_dino_limlight,     &
                      id_dino_limN,         &
                      id_dino_limP,         &
                      id_nano_limlight,     &
                      id_nano_limN,         &
                      id_nano_limP,         &
                      id_extinctioncoeff,   &
                      id_diat_columnprod,   &
                      id_dino_columnprod,   &
                      id_nano_columnprod,   &
                      id_totalchl,          &
                      id_columnprodtotal,   &
                      id_diatsettling,      &
                      id_detsettling,       &
                      id_spm_total
#if defined key_zoo_prod
   INTEGER, PUBLIC :: id_zoo_micr_columnprod,   &
                      id_zoo_meso_columnprod,   &
                      id_zoo_columnprodtotal
#endif	  		      
#if defined key_oyster_DEB
   INTEGER, PUBLIC :: id_oysDEB_FRpop,     &
                      id_oysDEB_PFpop,  &
                      id_oysDEB_Fpop,  &
                      id_oysDEB_DWtot,  &
                      id_oysDEB_kchlvble,  &
                      id_oysDEB_ER,  &
                      id_oysDEB_filtrate,  &
                      id_oysDEB_Wtot,  &
                      id_oysDEB_Istress,  &
                      id_oysDEB_Lgtot,  &
                      id_oysDEB_ChlC
 !  INTEGER, PUBLIC :: id_oysDEB_DWtot2,  &
 !                     id_oysDEB_Wtot2,  &
 !                     id_oysDEB_Lgtot2
 !  INTEGER, PUBLIC :: id_oysDEB_DWtot3,  &
 !                     id_oysDEB_Wtot3,  &
 !                     id_oysDEB_Lgtot3
   INTEGER, PUBLIC :: id_oysDEB_nbhuit
 !  INTEGER, PUBLIC :: id_oysDEB_nbhuit2
#endif
#if defined key_psnz
   INTEGER, PUBLIC :: id_psnz_limlight,     &
                      id_psnz_limN,         &
                      id_psnz_limSi,        &
                      id_psnz_limP,         &
                      id_psnz_proddomoic,   &
                      id_psnz_max,          &
                      id_psnz_datemax,      &
                      id_psnz_depthmax,     &
                      id_psnz_columnprod
#endif
#if defined key_karenia
   INTEGER, PUBLIC :: id_karenia_limlight,     &
                      id_karenia_limN,         &
                      id_karenia_limP,         &
                      id_karenia_max,          &
                      id_karenia_datemax,      &
                      id_karenia_depthmax,     &
                      id_karenia_columnprod
#endif
#if defined key_phaeocystis
   INTEGER, PUBLIC :: id_phaeocystis_limlight,     &
                      id_phaeocystis_limN,         &
                      id_phaeocystis_limP,         &
                      id_phaeocystis_max,          &
                      id_phaeocystis_datemax,      &
                      id_phaeocystis_depthmax,     &
                      id_phaeocystis_columnprod
#endif
#ifdef key_ulvas
   INTEGER, PUBLIC :: id_ulv_ratiosus,     &
                      id_ulv_susmort,      &
                      id_ulv_suslimlight,  &
                      id_ulv_suslimN,      &
                      id_ulv_suslimP,      &
                      id_ulv_suspumpN,     &
                      id_ulv_suspumpP,     &
                      id_ulv_settling,     &
                      id_ulv_ratiobenth,   &
                      id_ulv_benthmort,    &
                      id_ulv_benthlimlight,& 
                      id_ulv_benthlimN,    &
                      id_ulv_benthlimP,    &
                      id_ulv_benthpumpN,   &
                      id_ulv_benthpumpP,   &
                      id_ulv_resusbil,     &
                      id_ulv_resusN,       &
                      id_ulv_resusP
#endif
#ifdef key_zostera
   INTEGER, PUBLIC :: id_zost_LBDW,                &
                      id_zost_RBDW,                &
                      id_zost_Qcan,                &
                      id_zost_LAI,                 &
                      id_zost_limlum,              &
                      id_zost_limprodLN,           &
                      id_zost_limprodLP,           &
                      id_zost_limRN,               &
                      id_zost_limRP,               &
                      id_zost_effetchaleurzprod,   &
                      id_zost_effetchaleurzrespf,  &
                      id_zost_effetchaleurzrespr,  &
                      id_zost_limabsLnh4,          &
                      id_zost_limabsLno3,          &
                      id_zost_limabsLpo4,          &
                      id_zost_limabsRnh4,          &
                      id_zost_limabsRpo4,          &
                      id_zost_limselfshad,         &
                      id_zost_limRB,               &
                      id_zost_effetvent,           &
                      id_zost_Sgerm,               &
                      id_zost_ERS,                 &
                      id_zost_prod
#endif
#ifdef key_oyster_SFG
   INTEGER, PUBLIC :: id_oys_tcoq,          &
                      id_oys_tgam,          &
                      id_oys_gaincoq,       &
                      id_oys_gainsoma,      &
                      id_oys_gainrepr,      &
                      id_oys_orgafring,     &
                      id_oys_bilansoma,     &
                      id_oys_bilancoq,      &
                      id_oys_filt,          &
                      id_oys_resp,          &
                      id_oys_gameto,        &
                      id_oys_seuilponte,    &
                      id_oys_reposcoq,      &
                      id_oys_reposgam,      &
                      id_oys_ponte,         &
                      id_oys_absorg         
#endif
#ifdef key_oxygen
   INTEGER, PUBLIC :: id_oxy_sat                  
#endif
#ifdef key_BLOOM_insed
   INTEGER, PUBLIC :: id_remin_aerN ,      &
                      id_remin_aerP,       &
                      id_remin_aerSi,       &
                      id_remin_anaerN,       &
                      id_remin_anaerP,       &
                      id_remin_nitrateN,       &
                      id_remin_drnaN,       &
                      id_remin_denitN,       &
                      id_remin_nitrateP,       &
                      id_nitrif,       &
                      id_oxyd_solid_ODU,       &
                      id_adsor_desorb_P,       &
                      id_dissol_PFe,       &
                      id_precipit_P,       &
                      id_precipit_Si,       &
                      id_morta_phyto,       &
                      id_fluxsed_aeration,  &
                      id_porosite_sed,      &
                      id_filtr_benth
#ifdef key_sedim_MPB
   INTEGER, PUBLIC :: id_fluxsed_migC ,      &
                      id_fluxsed_migN,       &
                      id_fluxsed_exuN,       &
                      id_fluxsed_mortBAC,       &
                      id_fluxsed_mortMPBN,       &
                      id_fluxsed_mortMPBC,       &
                      id_fluxsed_proEPS,       &
                      id_fluxsed_respC,       &
                      id_fluxsed_resuspEPS,   &
                      id_fluxsed_uptBAC,       &
                      id_fluxsed_uptC,       &
                      id_fluxsedz_uptC,       &
                      id_fluxsed_uptN,       &
                      id_fluxsed_photoMPB_o2,   &
                      id_fluxsed_respMPB_o2,   &
                      id_lightlim,       &
                      id_light_z,        &
                      id_ratioCN ,    &
                      id_zactiv,    &
                      id_nkactiv,    &
                      id_rad
#endif
#endif

#ifdef key_microtracers
   INTEGER, PUBLIC   :: &
            id_tracer_maxdepth1, &   ! index for maximum value of micro-tracer from depth 1
            id_tracer_agedepth1, &   ! index for age of micro-tracer from depth 1
            id_tracer_meandepth1,&   ! index for mean value of micro-tracer from depth 1
            id_tracer_maxdepth2,&    ! index for maximum value of micro-tracer from depth 2
            id_tracer_agedepth2, &   ! index for age of micro-tracer from depth 2
            id_tracer_meandepth2,&   ! index for mean value of micro-tracer from depth 2
            id_tracer_maxdepth3,&    ! index for maximum value of micro-tracer from depth 3
            id_tracer_agedepth3, &   ! index for age of micro-tracer from depth 3
            id_tracer_meandepth3     ! index for mean value of micro-tracer from depth 3
#endif

#if defined key_N_tracer
   INTEGER, PUBLIC   :: ndiag_tracerN
! index for N-fraction from source  in signed variable
   INTEGER,DIMENSION(:,:),ALLOCATABLE, PUBLIC   :: id_tracer_signN 
#ifdef key_age_tracer
 ! index for age of N-fraction from source  in signed variable
   INTEGER,DIMENSION(:,:),ALLOCATABLE, PUBLIC   :: id_tracer_ageN 
#endif
#endif

#if defined key_P_tracer
   INTEGER, PUBLIC   :: ndiag_tracerP
! index for P-fraction from source  in signed variable
   INTEGER,DIMENSION(:,:),ALLOCATABLE, PUBLIC   :: id_tracer_signP 
#ifdef key_age_tracer
! index for age of P-fraction from source  in signed variable
   INTEGER,DIMENSION(:,:),ALLOCATABLE, PUBLIC   :: id_tracer_ageP  
#endif
#endif

#ifdef key_oyster_benthos
   REAL(KIND=rlg),DIMENSION(:,:),ALLOCATABLE :: nbhuitre,couvmicrophyto,hautable
   REAL(KIND=rlg),DIMENSION(:,:),ALLOCATABLE :: nbhuitre2,nbhuitre3
#elif defined key_oyster_SFG
   REAL(kind=rlg), DIMENSION(:,:),ALLOCATABLE :: nbhuitre,hautable,tpostpontecoq,tpostpontegam
#elif defined key_oyster_DEB
   REAL(KIND=rlg),DIMENSION(1:12) :: txmorthuitco1,txmorthuitco2,txmorthuitco3
   REAL(KIND=rlg),DIMENSION(:,:),ALLOCATABLE :: nbhuitre,hautable,chlorofondmoy,tempfondmax
#endif


#if defined key_MANGAbiovague
   ! configuration MANGAbio (reference) lecture d un fichier de vagues climato pour mieux prevoir MES
   ! nom du fichier lu dans namelist nammessat parasubs comme pour key_messat
   CHARACTER(LEN=lchain)         :: filevague_ubr !,filevague_hs
#endif

