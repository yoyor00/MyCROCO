
MODULE COMBLOOM

#include "cppdefs.h"
   
   USE comsubstance, only: rsh
!---------------------------------------------------------------------------
!
!                 include    combloom
!
!---------------------------------------------------------------------------

   IMPLICIT NONE
   
   LOGICAL, PUBLIC     :: l_filtbenthsinus,l_SNeffect_settle, &
                        l_physadaptation, l_filtbenthmes,l_phyzoodeteffect_settle,l_ChlNratio_var


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


   INTEGER, PUBLIC   :: &
            iv_oxygen          ! index for oxygen
   INTEGER, PUBLIC :: id_benthos_txf


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
                   p_diat_mort,        &
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
   INTEGER, PUBLIC :: id_oxy_sat                  
#ifdef key_BLOOM_insed
   INTEGER, PUBLIC :: id_remin_aerN ,      &
                      id_remin_aerP,       &
                      id_remin_aerSi,       &
                      id_remin_anaerN,       &
                      id_remin_anaerP,       &
                      id_remin_nitrateN,       &
                      id_remin_dnraN,       &
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



END MODULE 
