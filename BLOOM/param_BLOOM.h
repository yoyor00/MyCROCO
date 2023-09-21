! $Id: param_BLOOM_tracer.h  $

!!!***********************************************************************************
!  number of variables for module BLOOM depends on CPPkeys (modules bio add variables)
!!!***********************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  number of Bio advected variables without additional modules
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_BLOOM_opt2
      INTEGER,PARAMETER :: nb_var_bio=16
      INTEGER,PARAMETER :: nb_fix_bio=4
#else
#ifdef key_BLOOM_insed
      ! add PFe & ODU & NdetR & PdetR
      INTEGER,PARAMETER :: nb_var_bio=17
#else
      INTEGER,PARAMETER :: nb_var_bio=13
#endif
      INTEGER,PARAMETER :: nb_fix_bio=3
#endif

  !!***************************************
  !! optionnal modules
  !! number of additional variables
  !!========================================
#ifdef key_psnz
      INTEGER,PARAMETER :: nb_var_psnz=3
      INTEGER,PARAMETER :: nb_fix_psnz=1
#else
      INTEGER,PARAMETER :: nb_var_psnz=0
      INTEGER,PARAMETER :: nb_fix_psnz=0
#endif

#ifdef key_karenia
      INTEGER,PARAMETER :: nb_var_karenia=3
      INTEGER,PARAMETER :: nb_fix_karenia=1
#else
      INTEGER,PARAMETER :: nb_var_karenia=0
      INTEGER,PARAMETER :: nb_fix_karenia=0
#endif

#ifdef key_phaeocystis
      INTEGER,PARAMETER :: nb_var_phaeocystis=3
      INTEGER,PARAMETER :: nb_fix_phaeocystis=1
#else
      INTEGER,PARAMETER :: nb_var_phaeocystis=0
      INTEGER,PARAMETER :: nb_fix_phaeocystis=0
#endif

#ifdef key_zoo_prod
      INTEGER,PARAMETER :: nb_fix_zoo_ps=2
#else
      INTEGER,PARAMETER :: nb_fix_zoo_ps=0
#endif

#if defined key_oyster_SFG || defined key_oyster_DEB
      INTEGER,PARAMETER :: nb_fix_oys=3
#else
      INTEGER,PARAMETER :: nb_fix_oys=0
#endif

#ifdef key_oxygen
      INTEGER,PARAMETER :: nb_var_oxy=1
#else
      INTEGER,PARAMETER :: nb_var_oxy=0
#endif

#ifdef key_benthos
      INTEGER,PARAMETER :: nb_fix_benth_det=3
#ifdef key_diatbenth
#ifdef key_NPbenth
      INTEGER,PARAMETER :: nb_fix_benth_plus=3
#else
      INTEGER,PARAMETER :: nb_fix_benth_plus=1
#endif
#elif defined key_zostera
      INTEGER,PARAMETER :: nb_fix_benth_plus=3
#else
      INTEGER,PARAMETER :: nb_fix_benth_plus=0
#endif
#else
      INTEGER,PARAMETER :: nb_fix_benth_det=0
      INTEGER,PARAMETER :: nb_fix_benth_plus=0
#endif
#ifdef key_zostera
      INTEGER,PARAMETER :: nb_var_zost=3
      INTEGER,PARAMETER :: nb_fix_zost=10
#else
      INTEGER,PARAMETER :: nb_var_zost=0
      INTEGER,PARAMETER :: nb_fix_zost=0
#endif

#ifdef key_N_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! module N_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER,PARAMETER :: nb_sourceN=1
      INTEGER,PARAMETER :: nb_var_tracerN_defaut=8

!  modele de base=8, ajouter 1 pour Karenia ou pour Pseudo-Nitzschia
#ifdef key_psnz
      INTEGER,PARAMETER :: nb_var_tracerN_psnz=1
#else
      INTEGER,PARAMETER :: nb_var_tracerN_psnz=0
#endif
#ifdef key_karenia
      INTEGER,PARAMETER :: nb_var_tracerN_kar=1
#else
      INTEGER,PARAMETER :: nb_var_tracerN_kar=0
#endif
#ifdef key_phaeocystis
      INTEGER,PARAMETER :: nb_var_tracerN_phae=2
#else
      INTEGER,PARAMETER :: nb_var_tracerN_phae=0
#endif
#ifdef key_benthos
      INTEGER,PARAMETER :: nb_fix_tracerN=1
#else
      INTEGER,PARAMETER :: nb_fix_tracerN=0
#endif

      INTEGER,PARAMETER :: nb_vara_tracerN= nb_var_tracerN_defaut+
     &                                     nb_var_tracerN_phae+
     &                                     nb_var_tracerN_kar+
     &                                     nb_var_tracerN_psnz

#ifdef key_age_tracer
      INTEGER,PARAMETER :: nb_vara_N_age_tracer=nb_vara_tracerN
      INTEGER,PARAMETER :: nb_fix_N_age_tracer=nb_fix_tracerN
#else
      INTEGER,PARAMETER :: nb_vara_N_age_tracer=0
      INTEGER,PARAMETER :: nb_fix_N_age_tracer=0
#endif
      INTEGER,PARAMETER :: nb_var_tracerN= nb_vara_tracerN+nb_fix_tracerN
#else
      INTEGER,PARAMETER :: nb_sourceN=0
      INTEGER,PARAMETER :: nb_vara_tracerN=0
      INTEGER,PARAMETER :: nb_vara_N_age_tracer=0
      INTEGER,PARAMETER :: nb_fix_tracerN=0
      INTEGER,PARAMETER :: nb_fix_N_age_tracer=0
#endif


#ifdef key_P_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! module P_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER,PARAMETER :: nb_sourceP=2
      INTEGER,PARAMETER :: nb_var_tracerP_defaut=8
!       modele de base=8
#ifdef key_psnz
      INTEGER,PARAMETER :: nb_var_tracerP_psnz=1
#else
      INTEGER,PARAMETER :: nb_var_tracerP_psnz=0
#endif
#ifdef key_karenia
      INTEGER,PARAMETER :: nb_var_tracerP_kar=1
#else
      INTEGER,PARAMETER :: nb_var_tracerP_kar=0
#endif
#ifdef key_phaeocystis
      INTEGER,PARAMETER :: nb_var_tracerP_phae=2
#else
      INTEGER,PARAMETER :: nb_var_tracerP_phae=0
#endif
#ifdef key_benthos
      INTEGER,PARAMETER :: nb_fix_tracerP=1
#else
      INTEGER,PARAMETER :: nb_fix_tracerP=0
#endif
      INTEGER,PARAMETER :: nb_vara_tracerP= nb_var_tracerP_defaut+
     &                                     nb_var_tracerP_phae+
     &                                     nb_var_tracerP_kar+
     &                                     nb_var_tracerP_psnz
#ifdef key_age_tracer
      INTEGER,PARAMETER :: nb_vara_P_age_tracer=nb_vara_tracerP
      INTEGER,PARAMETER :: nb_fix_P_age_tracer=nb_fix_tracerP
#else
      INTEGER,PARAMETER :: nb_vara_P_age_tracer=0
      INTEGER,PARAMETER :: nb_fix_P_age_tracer=0
#endif
      INTEGER,PARAMETER :: nb_var_tracerP= nb_vara_tracerP+nb_fix_tracerP
#else
      INTEGER,PARAMETER :: nb_sourceP=0
      INTEGER,PARAMETER :: nb_vara_tracerP=0
      INTEGER,PARAMETER :: nb_vara_P_age_tracer=0
      INTEGER,PARAMETER :: nb_fix_tracerP=0
      INTEGER,PARAMETER :: nb_fix_P_age_tracer=0
#endif

      INTEGER, PARAMETER :: ntrc_biol=nb_var_bio+
     &                                  nb_var_psnz+
     &                                  nb_var_karenia+
     &                                  nb_var_phaeocystis+
     &                                  nb_var_oxy + 
     &                                  nb_vara_tracerN*nb_sourceN+
     &                                  nb_vara_N_age_tracer*nb_sourceN+
     &                                  nb_vara_tracerP*nb_sourceP+ 
     &                                  nb_vara_P_age_tracer*nb_sourceP+
     &                                  nb_var_zost
      INTEGER, PARAMETER :: ntfix     = nb_fix_bio+
     &                                  nb_fix_psnz+
     &                                  nb_fix_karenia+
     &                                  nb_fix_phaeocystis+
     &                                  nb_fix_oys+
     &                                  nb_fix_zoo_ps+
     &                                  nb_fix_benth_det+
     &                                  nb_fix_benth_plus+
     &                                  nb_fix_tracerP*nb_sourceP+
     &                                  nb_fix_P_age_tracer*nb_sourceP+
     &                                  nb_fix_tracerN*nb_sourceN+
     &                                  nb_fix_N_age_tracer*nb_sourceN+
     &                                  nb_fix_zost

