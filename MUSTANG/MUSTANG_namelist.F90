! Copyright (C) 2022-2026 IFREMER
! License: CeCILL-C
! See LICENSES/LICENSE_MUSTANG.txt

MODULE MUSTANG_namelist

#include "cppdefs.h"

#ifdef MUSTANG

!&E============================================================================
!&E                   ***  MODULE  MUSTANG_namelist  ***
!&E
!&E ** Purpose : declare the MUSTANG namelist groups (default values are set
!&E              directly on the corresponding variables in comMUSTANG) and
!&E              concerns reading and logging of the user namelist
!&E
!&E ** Description :
!&E     subroutine MUSTANG_readnml        ! read namelist
!&E     subroutine MUSTANG_param_log      ! write parameters for information
!&E
!&E============================================================================
    !! * Modules used
    USE comMUSTANG
    USE flocmod, ONLY : f_ws, f_diam, f_vol, f_rho, f_mass

    IMPLICIT NONE

    !! * Accessibility
    PUBLIC MUSTANG_readnml, MUSTANG_param_log

    PRIVATE

    !! * Local declarations
    ! declaration of namelists, variables described and defaulted in comMUSTANG
    namelist /namsedim_init/ filrepsed, l_unised,                   &
                             date_start_morpho, date_start_dyninsed,          &
                             hseduni, cseduni, ksmiuni, ksmauni,              &
                             sini_sed, tini_sed,                              &
                             l_unised_adjust_hsed, csed_mud_ini,              &
                             l_initsed_vardiss, poro_mud_ini

    namelist /namsedim_layer/ l_dzsminuni, dzsminuni,                         &
                              l_dzsmaxuni, dzsmaxuni,                         &
                              dzsmax_bottom, dzsmin,                          &
                              nlayer_surf_sed, l_splitlayersurf,              &
                              k1HW97, k2HW97,                                 &
                              fusion_para_activlayer

    namelist /namsedim_erosion/ activlayer, frmudcr2, coef_frmudcr1,          &
                                x1toce_mud, x2toce_mud,                       &
                                E0_sand_option, E0_sand_para, n_eros_sand,    &
                                E0_mud, n_eros_mud,                           &
                                ero_option, xexp_ero,                        &
                                E0_sand_Cst,                    &
                                tau_cri_option,                               &
                                tau_cri_mud_option_eroindep,                  &
                                l_peph_suspension, l_xexp_ero_cst,            &
                                l_eroindep_mud, l_eroindep_noncoh,            &
                                E0_mud_para_indep

    namelist /namsedim_bottomstress/ l_z0seduni,                              &
                                     z0seduni, z0sedmud, z0sedbedrock,        &
                                     l_tauskin_center,l_tauskin_ubar,         &
                                     l_tauskin_upwind, l_fricwave, fricwav,   &
                                     l_z0hydro_coupl_init,                    &
                                     l_z0hydro_coupl,                         &
                                     coef_z0_coupl,                           &
                                     z0_hydro_mud,                            &
                                     z0_hydro_bed

    namelist /namsedim_deposition/ cfreshmud, csedmin, cmudcr, aref_sand,     &
                                   l_corflux, cvolmaxsort, cvolmaxmel,        &
                                   l_slipdeposit, slopefac

    namelist /namsedim_lateral_erosion/ l_erolat, coef_erolat,&
                                        coef_tauskin_lat, l_erolat_wet_cell, &
                                        htncrit_eros

    namelist /namsedim_consolidation/ l_consolid, xperm1, xperm2, xsigma1,    &
                                      xsigma2, csegreg, csandseg,             &
                                      dt_consolid, subdt_consol

    namelist /namsedim_diffusion/ l_diffused, choice_flxdiss_diffsed,         &
                                  xdifs1, xdifsi1,           &
                                  epdifi, fexcs, dt_diffused

    namelist /namsedim_bioturb/ l_bioturb, l_biodiffs,                        &
                                xbioturbmax_part, xbioturbk_part,             &
                                dbiotu0_part, dbiotum_part,                   &
                                xbioturbmax_diss, xbioturbk_diss,             &
                                dbiotu0_diss, dbiotum_diss,                   &
                                frmud_db_min, frmud_db_max,                   &
                                dt_bioturb, subdt_bioturb

    namelist /namsedim_morpho/ l_morphocoupl, MF, l_MF_dhsed, dt_morpho

    namelist /namsedoutput/ choice_nivsed_out,                                &
                            nk_nivsed_out, ep_nivsed_out, epmax_nivsed_out,   &
                            l_outsed_nb_lay_sed, &
                            l_outsed_hsed, &
                            l_outsed_tauskin, &
                            l_outsed_tauskin_cw, &
                            l_outsed_poro, &
                            l_outsed_dzs, &
                            l_outsed_temp_sed, &
                            l_outsed_salt_sed, &
                            l_outsed_cv_sed, &
                            l_outsed_ws, &
                            l_outsed_toce, &
                            l_outsed_flx_s2w_w2s, &
                            l_outsed_pephm_fcor, &
                            l_outsed_bedload, &
                            l_outsed_fsusp, &
                            l_outsed_frmudsup, &
                            l_outsed_dzs_ksmax, &
                            l_outsed_theoric_active_layer, &
                            l_outsed_ero_details, &
                            l_outsed_z0sed, &
                            l_outsed_z0hydro, &
                            l_outsed_consolidation

    namelist /namdredging/ dredging_location_file, dredging_settings_file,    &
            dredging_out_file, dredging_dumping_layer, &
            dredging_dt, dredging_dt_out

    namelist /namsedim_poro/ poro_option, poro_min,                           &
                             Awooster, Bwooster, Bmax_wu
#ifdef key_MUSTANG_V2
    namelist /namsedim_bedload/ l_peph_bedload, l_slope_effect_bedload,       &
                                alphabs, alphabn, hmin_bedload, l_fsusp
#endif

    namelist /namflocmod/ l_flocmod, l_ADS, l_ASH, l_COLLFRAG,                &
                          f_dp0, f_nf, f_nb_frag, f_alpha, f_beta, f_ater,    &
                          f_ero_frac, f_ero_nbfrag, f_ero_iv, f_mneg_param,   &
                          f_collfragparam, f_dmin_frag, f_cfcst, f_fp, f_fy,  &
                          f_clim
#if !defined key_noTSdiss_insed
    namelist /namtempsed/ mu_tempsed1, mu_tempsed2, mu_tempsed3,              &
                          epsedmin_tempsed,                                   &
                          epsedmax_tempsed
#endif

CONTAINS

!!=============================================================================
    SUBROUTINE MUSTANG_readnml(filein)
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_readnml  ***
    !&E
    !&E ** Purpose : reads namelist file from filein
    !&E
    !&E ** Description : read namelist file paraMUSTANG
    !&E
    !&E ** Called by :  MUSTANG_init
    !&E
    !&E--------------------------------------------------------------------------
    !! * Arguments
    CHARACTER(len=lchain), INTENT(IN) :: filein

    !! * Executable part
    OPEN(unit = 50, file = filein, status = 'old', action = 'read')

    MPI_master_only write(*,*) '*****************************************************'
    MPI_master_only write(*,*) 'READING MUSTANG input file'
    MPI_master_only write(*,*) TRIM(filein)
    MPI_master_only write(*,*) '*****************************************************'

    READ(50, namsedim_init); rewind(50)
    READ(50, namsedim_layer); rewind(50)
    READ(50, namsedim_bottomstress); rewind(50)
    READ(50, namsedim_deposition); rewind(50)
    READ(50, namsedim_erosion); rewind(50)
    READ(50, namsedim_poro); rewind(50)
#ifdef key_MUSTANG_V2
    READ(50, namsedim_bedload); rewind(50)
#endif
    READ(50, namsedim_lateral_erosion); rewind(50)
    READ(50, namsedim_consolidation); rewind(50)
    READ(50, namsedim_diffusion); rewind(50)
    READ(50, namsedim_bioturb); rewind(50)
    READ(50, namsedim_morpho); rewind(50)
#if !defined key_noTSdiss_insed
    READ(50, namtempsed); rewind(50)
#endif
    READ(50, namsedoutput); rewind(50)
    ! module FLOCULATION
    READ(50, namflocmod); rewind(50)
    READ(50, namdredging); rewind(50)

    CLOSE(50)

    END SUBROUTINE MUSTANG_readnml
!!===========================================================================

    SUBROUTINE MUSTANG_param_log(iscreenlog)
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_param_log  ***
    !&E
    !&E ** Purpose : write parameters in log file
    !&E
    !&E ** Description : namelists writing for information
    !&E
    !&E ** Called by :  MUSTANG_init
    !&E
    !&E--------------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT(IN) :: iscreenlog

    !! * Local declarations
    INTEGER   :: iv

    MPI_master_only WRITE(iscreenlog, namsedim_init)
    MPI_master_only WRITE(iscreenlog, namsedim_layer)
    MPI_master_only WRITE(iscreenlog, namsedim_bottomstress)
    MPI_master_only WRITE(iscreenlog, namsedim_deposition)
    MPI_master_only WRITE(iscreenlog, namsedim_poro)
#ifdef key_MUSTANG_V2
    MPI_master_only WRITE(iscreenlog, namsedim_bedload)
#else
    MPI_master_only WRITE(iscreenlog, namsedim_erosion)
#endif
    MPI_master_only WRITE(iscreenlog, namsedim_lateral_erosion)
    MPI_master_only WRITE(iscreenlog, namsedim_consolidation)
    MPI_master_only write(iscreenlog, *) 'MUSTANG_param_log dt_consolid=', dt_consolid
    MPI_master_only WRITE(iscreenlog, namsedim_diffusion)
    MPI_master_only WRITE(iscreenlog, namsedim_bioturb)
    MPI_master_only WRITE(iscreenlog, namsedim_morpho)
#if !defined key_noTSdiss_insed
    MPI_master_only WRITE(iscreenlog, namtempsed)
#endif

    IF (l_flocmod) THEN
    !! module floculation
    MPI_master_only WRITE(iscreenlog, *) ' '
    MPI_master_only WRITE(iscreenlog, *) '    FLOCMOD'
    MPI_master_only WRITE(iscreenlog, *) '***********************'
    MPI_master_only WRITE(iscreenlog, *) 'class  diameter  volume  density  mass Ws'
    DO iv = 1, nv_mud
        MPI_master_only WRITE(iscreenlog, *) iv, f_diam(iv), f_vol(iv), f_rho(iv), f_mass(iv), f_ws(iv)
    ENDDO
    MPI_master_only WRITE(iscreenlog, *) ' '
    MPI_master_only WRITE(iscreenlog, *) ' *** PARAMETERS ***'
    MPI_master_only WRITE(iscreenlog, *) &
        'Primary particle size (f_dp0)                                : ', f_dp0
    MPI_master_only WRITE(iscreenlog, *) &
        'Fractal dimension (f_nf)                                     : ', f_nf
    MPI_master_only WRITE(iscreenlog, *) &
        'Flocculation efficiency (f_alpha)                            : ', f_alpha
    MPI_master_only WRITE(iscreenlog, *) &
        'Floc break up parameter (f_beta)                             : ', f_beta
    MPI_master_only WRITE(iscreenlog, *) &
        'Nb of fragments (f_nb_frag)                                  : ', f_nb_frag
    MPI_master_only WRITE(iscreenlog, *) &
        'Ternary fragmentation (f_ater)                               : ', f_ater
    MPI_master_only WRITE(iscreenlog, *) &
        'Floc erosion (% of mass) (f_ero_frac)                        : ', f_ero_frac
    MPI_master_only WRITE(iscreenlog, *) &
        'Nb of fragments by erosion (f_ero_nbfrag)                    : ', f_ero_nbfrag
    MPI_master_only WRITE(iscreenlog, *) &
        'fragment class (f_ero_iv)                                    : ', f_ero_iv
    MPI_master_only WRITE(iscreenlog, *) &
        'negative mass tolerated before redistribution (f_mneg_param) : ', f_mneg_param
    MPI_master_only WRITE(iscreenlog, *) &
        'Boolean for differential settling aggregation (L_ADS)        : ', l_ADS
    MPI_master_only WRITE(iscreenlog, *) &
        'Boolean for shear aggregation (L_ASH)                        : ', l_ASH
    MPI_master_only WRITE(iscreenlog, *) &
        'Boolean for collision fragmenation (L_COLLFRAG)              : ', l_COLLFRAG
    MPI_master_only WRITE(iscreenlog, *) &
        'Collision fragmentation parameter (f_collfragparam)          : ', f_collfragparam
    MPI_master_only WRITE(iscreenlog, *) &
        'Min concentration below which flocculation is not calculated : ', f_clim
    MPI_master_only WRITE(iscreenlog, *) ' '
    MPI_master_only WRITE(iscreenlog, *) '*** END FLOCMOD INIT *** '
    ENDIF

    END SUBROUTINE MUSTANG_param_log
!!===========================================================================

#endif /* ifdef MUSTANG */

END MODULE MUSTANG_namelist
