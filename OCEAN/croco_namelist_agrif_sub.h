!=======================================================================
!  croco_namelist_agrif_sub.h: Multi-grid save/restore subroutines for AGRIF
!=======================================================================
   subroutine croco_namelist_save_grid(gid_in)
      integer, intent(in) :: gid_in
      integer :: gid
      if (gid_in <= 0) then
         gid = 0
      else
         gid = gid_in
      endif
      if (gid > MAX_AGRIF_GRIDS) return
      agrif_nml_storage(gid)%fname_nml = fname_nml
      agrif_nml_storage(gid)%testcase_name = testcase_name
      agrif_nml_storage(gid)%lts_hadv = lts_hadv
      agrif_nml_storage(gid)%lts_rot_or_diag = lts_rot_or_diag
      agrif_nml_storage(gid)%lts_per = lts_per
      agrif_nml_storage(gid)%lekman_apg = lekman_apg
      agrif_nml_storage(gid)%loscnonrot = loscnonrot
      agrif_nml_storage(gid)%lnonrot = lnonrot
      agrif_nml_storage(gid)%UG_forced = UG_forced
      agrif_nml_storage(gid)%wave_ramp_duration = wave_ramp_duration
      agrif_nml_storage(gid)%title = title
      agrif_nml_storage(gid)%logname = logname
      agrif_nml_storage(gid)%dt = dt
      agrif_nml_storage(gid)%ntimes = ntimes
      agrif_nml_storage(gid)%ndtfast = ndtfast
      agrif_nml_storage(gid)%ninfo = ninfo
      agrif_nml_storage(gid)%m2filter_alpha = m2filter_alpha
      agrif_nml_storage(gid)%ldefhis = ldefhis
      agrif_nml_storage(gid)%nwrt = nwrt
      agrif_nml_storage(gid)%nrpfhis = nrpfhis
      agrif_nml_storage(gid)%hisname = hisname
      agrif_nml_storage(gid)%nrrec = nrrec
      agrif_nml_storage(gid)%ininame = ininame
      agrif_nml_storage(gid)%nrst = nrst
      agrif_nml_storage(gid)%nrpfrst = nrpfrst
      agrif_nml_storage(gid)%rstname = rstname
#if defined NBQ || defined K3FAST   /* H1-fix 2026-08-26: directive invalide (#ifdef NBQ || ...) */
      agrif_nml_storage(gid)%ndtnbq_in = ndtnbq_in
      agrif_nml_storage(gid)%csound_nbq = csound_nbq
      agrif_nml_storage(gid)%visc2read_nbq = visc2read_nbq
#endif
#ifdef SOLVE3D
      agrif_nml_storage(gid)%theta_s = theta_s
      agrif_nml_storage(gid)%theta_b = theta_b
      agrif_nml_storage(gid)%hc = hc
      agrif_nml_storage(gid)%Tcline = Tcline
#endif
#ifdef USE_CALENDAR
      agrif_nml_storage(gid)%start_date = start_date
      agrif_nml_storage(gid)%end_date = end_date
      agrif_nml_storage(gid)%dt_his = dt_his
      agrif_nml_storage(gid)%dt_avg = dt_avg
      agrif_nml_storage(gid)%dt_rst = dt_rst
#endif
#ifndef ANA_GRID
      agrif_nml_storage(gid)%grdname = grdname
#endif
      agrif_nml_storage(gid)%frcname = frcname
      agrif_nml_storage(gid)%use_frcname = use_frcname
      agrif_nml_storage(gid)%bulkname = bulkname
#if (  defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || (defined ZCLIMATOLOGY&& !defined ANA_SSH) || (defined M2CLIMATOLOGY&& !defined ANA_M2CLIMA) || (defined M3CLIMATOLOGY&& !defined ANA_M3CLIMA)
      agrif_nml_storage(gid)%clmname = clmname
#endif
#if !defined ANA_BRY && defined FRC_BRY
      agrif_nml_storage(gid)%bry_file = bry_file
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      agrif_nml_storage(gid)%brywkb_file = brywkb_file
#endif
#ifdef WKB_WWAVE
      agrif_nml_storage(gid)%wkb_amp = wkb_amp
      agrif_nml_storage(gid)%wkb_ang = wkb_ang
      agrif_nml_storage(gid)%wkb_prd = wkb_prd
      agrif_nml_storage(gid)%wkb_tide = wkb_tide
      agrif_nml_storage(gid)%wkb_btg = wkb_btg
      agrif_nml_storage(gid)%wkb_gam = wkb_gam
#  ifdef WAVE_ROLLER
      agrif_nml_storage(gid)%wkb_rsb = wkb_rsb
      agrif_nml_storage(gid)%wkb_roller = wkb_roller
#  endif
#endif
#ifdef WAVE_MAKER
      agrif_nml_storage(gid)%wmaker_amp = wmaker_amp
      agrif_nml_storage(gid)%wmaker_prd = wmaker_prd
      agrif_nml_storage(gid)%wmaker_dir = wmaker_dir
      agrif_nml_storage(gid)%wmaker_dsp = wmaker_dsp
      agrif_nml_storage(gid)%wmaker_fsp = wmaker_fsp
      agrif_nml_storage(gid)%wmaker_wf1 = wmaker_wf1
      agrif_nml_storage(gid)%wmaker_wf2 = wmaker_wf2
      agrif_nml_storage(gid)%wmaker_wa1 = wmaker_wa1
      agrif_nml_storage(gid)%wmaker_wa2 = wmaker_wa2
#endif
#ifdef AVERAGES
      agrif_nml_storage(gid)%ntsavg = ntsavg
      agrif_nml_storage(gid)%navg = navg
      agrif_nml_storage(gid)%nrpfavg = nrpfavg
      agrif_nml_storage(gid)%avgname = avgname
#endif
#if defined OUTPUTS_SURFACE && !defined XIOS
      agrif_nml_storage(gid)%ldefsurf = ldefsurf
      agrif_nml_storage(gid)%nwrtsurf = nwrtsurf
      agrif_nml_storage(gid)%nrpfsurf = nrpfsurf
      agrif_nml_storage(gid)%surfname = surfname
      agrif_nml_storage(gid)%out_his_surf = out_his_surf
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefsurf_avg = ldefsurf_avg
      agrif_nml_storage(gid)%ntssurf_avg = ntssurf_avg
      agrif_nml_storage(gid)%nwrtsurf_avg = nwrtsurf_avg
      agrif_nml_storage(gid)%nrpfsurf_avg = nrpfsurf_avg
      agrif_nml_storage(gid)%surfname_avg = surfname_avg
      agrif_nml_storage(gid)%out_avg_surf = out_avg_surf
#  endif
#endif
#if defined ABL1D && !defined XIOS
      agrif_nml_storage(gid)%ldefablhis = ldefablhis
      agrif_nml_storage(gid)%nwrtablhis = nwrtablhis
      agrif_nml_storage(gid)%nrpfablhis = nrpfablhis
      agrif_nml_storage(gid)%ablname = ablname
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefablavg = ldefablavg
      agrif_nml_storage(gid)%ntsablavg = ntsablavg
      agrif_nml_storage(gid)%nwrtablavg = nwrtablavg
      agrif_nml_storage(gid)%nrpfablavg = nrpfablavg
      agrif_nml_storage(gid)%ablname_avg = ablname_avg
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      agrif_nml_storage(gid)%wave_file = wave_file
#endif
#if defined BIOLOGY && defined PISCES
      agrif_nml_storage(gid)%bioname = bioname
#endif
#ifdef BODYFORCE
      agrif_nml_storage(gid)%levsfrc = levsfrc
      agrif_nml_storage(gid)%levbfrc = levbfrc
#endif
      agrif_nml_storage(gid)%R0 = R0
      agrif_nml_storage(gid)%T0 = T0
      agrif_nml_storage(gid)%S0 = S0
      agrif_nml_storage(gid)%Tcoef = Tcoef
      agrif_nml_storage(gid)%Scoef = Scoef
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
      agrif_nml_storage(gid)%ltra_min = ltra_min
      agrif_nml_storage(gid)%ltra_max = ltra_max
#endif
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
      agrif_nml_storage(gid)%ldyn_min = ldyn_min
      agrif_nml_storage(gid)%ldyn_max = ldyn_max
#endif
#ifdef SEDIMENT
      agrif_nml_storage(gid)%sedname = sedname
#endif
#ifdef MUSTANG
      agrif_nml_storage(gid)%sedname_must = sedname_must
#endif
#ifdef SUBSTANCE
      agrif_nml_storage(gid)%subsfilename = subsfilename
#endif
#ifdef OBSTRUCTION
      agrif_nml_storage(gid)%obstname = obstname
#endif
#ifdef XIOS
      agrif_nml_storage(gid)%xios_origin_date = xios_origin_date
#endif
#ifdef ASSIMILATION
      agrif_nml_storage(gid)%aparnam = aparnam
      agrif_nml_storage(gid)%assname = assname
#endif
      agrif_nml_storage(gid)%rho0 = rho0
      agrif_nml_storage(gid)%rdrg = rdrg
      agrif_nml_storage(gid)%rdrg2 = rdrg2
      agrif_nml_storage(gid)%Zobt = Zobt
      agrif_nml_storage(gid)%Cdb_min = Cdb_min
      agrif_nml_storage(gid)%Cdb_max = Cdb_max
      agrif_nml_storage(gid)%gamma2 = gamma2
      agrif_nml_storage(gid)%visc2 = visc2
      agrif_nml_storage(gid)%visc4 = visc4
#ifdef SOLVE3D
#  ifdef TRACERS
      if (allocated(tnu2)) then
         if (allocated(agrif_nml_storage(gid)%tnu2)) deallocate(agrif_nml_storage(gid)%tnu2)
         allocate(agrif_nml_storage(gid)%tnu2(size(tnu2)))
         agrif_nml_storage(gid)%tnu2 = tnu2
      endif
      if (allocated(tnu4)) then
         if (allocated(agrif_nml_storage(gid)%tnu4)) deallocate(agrif_nml_storage(gid)%tnu4)
         allocate(agrif_nml_storage(gid)%tnu4(size(tnu4)))
         agrif_nml_storage(gid)%tnu4 = tnu4
      endif
#  endif
#  if !defined LMD_MIXING
      agrif_nml_storage(gid)%Akv_bak = Akv_bak
#    ifdef TRACERS
      if (allocated(Akt_bak)) then
         if (allocated(agrif_nml_storage(gid)%Akt_bak)) deallocate(agrif_nml_storage(gid)%Akt_bak)
         allocate(agrif_nml_storage(gid)%Akt_bak(size(Akt_bak)))
         agrif_nml_storage(gid)%Akt_bak = Akt_bak
      endif
#    endif
#  endif
#endif
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
      agrif_nml_storage(gid)%psource_Nsrc = psource_Nsrc
      if (allocated(psource_Isrc)) then
         if (allocated(agrif_nml_storage(gid)%psource_Isrc)) deallocate(agrif_nml_storage(gid)%psource_Isrc)
         allocate(agrif_nml_storage(gid)%psource_Isrc(size(psource_Isrc)))
         agrif_nml_storage(gid)%psource_Isrc = psource_Isrc
      endif
      if (allocated(psource_Jsrc)) then
         if (allocated(agrif_nml_storage(gid)%psource_Jsrc)) deallocate(agrif_nml_storage(gid)%psource_Jsrc)
         allocate(agrif_nml_storage(gid)%psource_Jsrc(size(psource_Jsrc)))
         agrif_nml_storage(gid)%psource_Jsrc = psource_Jsrc
      endif
      if (allocated(psource_Dsrc)) then
         if (allocated(agrif_nml_storage(gid)%psource_Dsrc)) deallocate(agrif_nml_storage(gid)%psource_Dsrc)
         allocate(agrif_nml_storage(gid)%psource_Dsrc(size(psource_Dsrc)))
         agrif_nml_storage(gid)%psource_Dsrc = psource_Dsrc
      endif
      agrif_nml_storage(gid)%psource_qbarname = psource_qbarname
      if (allocated(psource_qbardir)) then
         if (allocated(agrif_nml_storage(gid)%psource_qbardir)) deallocate(agrif_nml_storage(gid)%psource_qbardir)
         allocate(agrif_nml_storage(gid)%psource_qbardir(size(psource_qbardir)))
         agrif_nml_storage(gid)%psource_qbardir = psource_qbardir
      endif
      if (allocated(psource_Qbar)) then
         if (allocated(agrif_nml_storage(gid)%psource_Qbar)) deallocate(agrif_nml_storage(gid)%psource_Qbar)
         allocate(agrif_nml_storage(gid)%psource_Qbar(size(psource_Qbar)))
         agrif_nml_storage(gid)%psource_Qbar = psource_Qbar
      endif
# if defined TRACERS
      if (allocated(psource_Lsrc)) then
         if (allocated(agrif_nml_storage(gid)%psource_Lsrc)) deallocate(agrif_nml_storage(gid)%psource_Lsrc)
         allocate(agrif_nml_storage(gid)%psource_Lsrc(size(psource_Lsrc,1), size(psource_Lsrc,2)))
         agrif_nml_storage(gid)%psource_Lsrc = psource_Lsrc
      endif
      if (allocated(psource_Tsrc0)) then
         if (allocated(agrif_nml_storage(gid)%psource_Tsrc0)) deallocate(agrif_nml_storage(gid)%psource_Tsrc0)
         allocate(agrif_nml_storage(gid)%psource_Tsrc0(size(psource_Tsrc0,1), size(psource_Tsrc0,2)))
         agrif_nml_storage(gid)%psource_Tsrc0 = psource_Tsrc0
      endif
# endif
#endif
#if defined DIAGNOSTICS_TS
      agrif_nml_storage(gid)%ldefdia = ldefdia
      agrif_nml_storage(gid)%nwrtdia = nwrtdia
      agrif_nml_storage(gid)%nrpfdia = nrpfdia
      agrif_nml_storage(gid)%dianame = dianame
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdia_avg = ldefdia_avg
      agrif_nml_storage(gid)%ntsdia_avg = ntsdia_avg
      agrif_nml_storage(gid)%nwrtdia_avg = nwrtdia_avg
      agrif_nml_storage(gid)%nrpfdia_avg = nrpfdia_avg
      agrif_nml_storage(gid)%dianame_avg = dianame_avg
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
      agrif_nml_storage(gid)%mld_crit_D = mld_crit_D
      agrif_nml_storage(gid)%mld_crit_T = mld_crit_T
#  endif
#endif
#if defined DIAGNOSTICS_UV
      agrif_nml_storage(gid)%ldefdiaM = ldefdiaM
      agrif_nml_storage(gid)%nwrtdiaM = nwrtdiaM
      agrif_nml_storage(gid)%nrpfdiaM = nrpfdiaM
      agrif_nml_storage(gid)%dianameM = dianameM
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdiaM_avg = ldefdiaM_avg
      agrif_nml_storage(gid)%ntsdiaM_avg = ntsdiaM_avg
      agrif_nml_storage(gid)%nwrtdiaM_avg = nwrtdiaM_avg
      agrif_nml_storage(gid)%nrpfdiaM_avg = nrpfdiaM_avg
      agrif_nml_storage(gid)%dianameM_avg = dianameM_avg
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      agrif_nml_storage(gid)%ldefdiags_vrt = ldefdiags_vrt
      agrif_nml_storage(gid)%nwrtdiags_vrt = nwrtdiags_vrt
      agrif_nml_storage(gid)%nrpfdiags_vrt = nrpfdiags_vrt
      agrif_nml_storage(gid)%diags_vrtname = diags_vrtname
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdiags_vrt_avg = ldefdiags_vrt_avg
      agrif_nml_storage(gid)%ntsdiags_vrt_avg = ntsdiags_vrt_avg
      agrif_nml_storage(gid)%nwrtdiags_vrt_avg = nwrtdiags_vrt_avg
      agrif_nml_storage(gid)%nrpfdiags_vrt_avg = nrpfdiags_vrt_avg
      agrif_nml_storage(gid)%diags_vrtname_avg = diags_vrtname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_KE
      agrif_nml_storage(gid)%ldefdiags_ek = ldefdiags_ek
      agrif_nml_storage(gid)%nwrtdiags_ek = nwrtdiags_ek
      agrif_nml_storage(gid)%nrpfdiags_ek = nrpfdiags_ek
      agrif_nml_storage(gid)%diags_ekname = diags_ekname
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdiags_ek_avg = ldefdiags_ek_avg
      agrif_nml_storage(gid)%ntsdiags_ek_avg = ntsdiags_ek_avg
      agrif_nml_storage(gid)%nwrtdiags_ek_avg = nwrtdiags_ek_avg
      agrif_nml_storage(gid)%nrpfdiags_ek_avg = nrpfdiags_ek_avg
      agrif_nml_storage(gid)%diags_ekname_avg = diags_ekname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      agrif_nml_storage(gid)%ldefdiags_pv = ldefdiags_pv
      agrif_nml_storage(gid)%nwrtdiags_pv = nwrtdiags_pv
      agrif_nml_storage(gid)%nrpfdiags_pv = nrpfdiags_pv
      agrif_nml_storage(gid)%diags_pvname = diags_pvname
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdiags_pv_avg = ldefdiags_pv_avg
      agrif_nml_storage(gid)%ntsdiags_pv_avg = ntsdiags_pv_avg
      agrif_nml_storage(gid)%nwrtdiags_pv_avg = nwrtdiags_pv_avg
      agrif_nml_storage(gid)%nrpfdiags_pv_avg = nrpfdiags_pv_avg
      agrif_nml_storage(gid)%diags_pvname_avg = diags_pvname_avg
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdiags_eddy_avg = ldefdiags_eddy_avg
      agrif_nml_storage(gid)%ntsdiags_eddy_avg = ntsdiags_eddy_avg
      agrif_nml_storage(gid)%nwrtdiags_eddy_avg = nwrtdiags_eddy_avg
      agrif_nml_storage(gid)%nrpfdiags_eddy_avg = nrpfdiags_eddy_avg
      agrif_nml_storage(gid)%diags_eddyname_avg = diags_eddyname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      agrif_nml_storage(gid)%ldefdiabio = ldefdiabio
      agrif_nml_storage(gid)%nwrtdiabio = nwrtdiabio
      agrif_nml_storage(gid)%nrpfdiabio = nrpfdiabio
      agrif_nml_storage(gid)%dianamebio = dianamebio
#  ifdef AVERAGES
      agrif_nml_storage(gid)%ldefdiabio_avg = ldefdiabio_avg
      agrif_nml_storage(gid)%ntsdiabio_avg = ntsdiabio_avg
      agrif_nml_storage(gid)%nwrtdiabio_avg = nwrtdiabio_avg
      agrif_nml_storage(gid)%nrpfdiabio_avg = nrpfdiabio_avg
      agrif_nml_storage(gid)%dianamebio_avg = dianamebio_avg
#  endif
#endif
#ifdef STATIONS
      agrif_nml_storage(gid)%ldefsta = ldefsta
      agrif_nml_storage(gid)%nsta = nsta
      agrif_nml_storage(gid)%nrpfsta = nrpfsta
      agrif_nml_storage(gid)%staposname = staposname
      agrif_nml_storage(gid)%staname = staname
      agrif_nml_storage(gid)%sta_grd = sta_grd
      agrif_nml_storage(gid)%sta_temp = sta_temp
      agrif_nml_storage(gid)%sta_salt = sta_salt
      agrif_nml_storage(gid)%sta_rho = sta_rho
      agrif_nml_storage(gid)%sta_vel = sta_vel
#endif
#if (defined SPONGE && !defined SPONGE_GRID)
      agrif_nml_storage(gid)%x_sponge = x_sponge
#endif
#ifdef WET_DRY
      agrif_nml_storage(gid)%D_wetdry = D_wetdry
#endif
#if defined MRL_WCI && defined WAVE_DRY
      agrif_nml_storage(gid)%D_wavedry = D_wavedry
#endif
#if defined T_FRC_BRY     || defined M2_FRC_BRY    || defined M3_FRC_BRY||defined Z_FRC_BRY|| defined W_FRC_BRY||defined NBQ_FRC_BRY|| defined TCLIMATOLOGY||defined M2CLIMATOLOGY|| defined M3CLIMATOLOGY||defined ZCLIMATOLOGY|| defined WCLIMATOLOGY||defined NBQCLIMATOLOGY|| defined TNUDGING||defined ZNUDGING|| defined M2NUDGING||defined M3NUDGING|| defined NBQ_NUDGING||defined SPONGE
      agrif_nml_storage(gid)%tauT_in = tauT_in
      agrif_nml_storage(gid)%tauT_out = tauT_out
      agrif_nml_storage(gid)%tauM_in = tauM_in
      agrif_nml_storage(gid)%tauM_out = tauM_out
#endif
#if defined BHFLUX || (defined BWFLUX && defined SALINITY)
      agrif_nml_storage(gid)%btfname = btfname
#endif
      agrif_nml_storage(gid)%out_his_zeta = out_his_zeta
      agrif_nml_storage(gid)%out_his_ubar = out_his_ubar
      agrif_nml_storage(gid)%out_his_vbar = out_his_vbar
#ifdef SOLVE3D
      agrif_nml_storage(gid)%out_his_u = out_his_u
      agrif_nml_storage(gid)%out_his_v = out_his_v
# ifdef TRACERS
      if (allocated(out_his_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_his_tracer)) deallocate(agrif_nml_storage(gid)%out_his_tracer)
         allocate(agrif_nml_storage(gid)%out_his_tracer(size(out_his_tracer)))
         agrif_nml_storage(gid)%out_his_tracer = out_his_tracer
      endif
# endif
#endif
#ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_zeta = out_avg_zeta
      agrif_nml_storage(gid)%out_avg_ubar = out_avg_ubar
      agrif_nml_storage(gid)%out_avg_vbar = out_avg_vbar
# ifdef SOLVE3D
      agrif_nml_storage(gid)%out_avg_u = out_avg_u
      agrif_nml_storage(gid)%out_avg_v = out_avg_v
#  ifdef TRACERS
      if (allocated(out_avg_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_tracer)) deallocate(agrif_nml_storage(gid)%out_avg_tracer)
         allocate(agrif_nml_storage(gid)%out_avg_tracer(size(out_avg_tracer)))
         agrif_nml_storage(gid)%out_avg_tracer = out_avg_tracer
      endif
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_TS
# ifdef TRACERS
      if (allocated(out_his_dia3D_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_his_dia3D_tracer)) deallocate(agrif_nml_storage(gid)%out_his_dia3D_tracer)
         allocate(agrif_nml_storage(gid)%out_his_dia3D_tracer(size(out_his_dia3D_tracer)))
         agrif_nml_storage(gid)%out_his_dia3D_tracer = out_his_dia3D_tracer
      endif
#  ifdef DIAGNOSTICS_TS_MLD
      if (allocated(out_his_dia2D_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_his_dia2D_tracer)) deallocate(agrif_nml_storage(gid)%out_his_dia2D_tracer)
         allocate(agrif_nml_storage(gid)%out_his_dia2D_tracer(size(out_his_dia2D_tracer)))
         agrif_nml_storage(gid)%out_his_dia2D_tracer = out_his_dia2D_tracer
      endif
#  endif
# endif
# if defined AVERAGES && defined TRACERS
      if (allocated(out_avg_dia3D_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_dia3D_tracer)) deallocate(agrif_nml_storage(gid)%out_avg_dia3D_tracer)
         allocate(agrif_nml_storage(gid)%out_avg_dia3D_tracer(size(out_avg_dia3D_tracer)))
         agrif_nml_storage(gid)%out_avg_dia3D_tracer = out_avg_dia3D_tracer
      endif
#  ifdef DIAGNOSTICS_TS_MLD
      if (allocated(out_avg_dia2D_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_dia2D_tracer)) deallocate(agrif_nml_storage(gid)%out_avg_dia2D_tracer)
         allocate(agrif_nml_storage(gid)%out_avg_dia2D_tracer(size(out_avg_dia2D_tracer)))
         agrif_nml_storage(gid)%out_avg_dia2D_tracer = out_avg_dia2D_tracer
      endif
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_UV
      agrif_nml_storage(gid)%out_his_diagM_u = out_his_diagM_u
      agrif_nml_storage(gid)%out_his_diagM_v = out_his_diagM_v
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_diagM_u = out_avg_diagM_u
      agrif_nml_storage(gid)%out_avg_diagM_v = out_avg_diagM_v
# endif
#endif
#ifdef DIAGNOSTICS_VRT
      agrif_nml_storage(gid)%out_his_diags_vrt = out_his_diags_vrt
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_diags_vrt = out_avg_diags_vrt
# endif
#endif
#ifdef DIAGNOSTICS_KE
      agrif_nml_storage(gid)%out_his_diags_ek = out_his_diags_ek
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_diags_ek = out_avg_diags_ek
# endif
#endif
#if defined DIAGNOSTICS_PV && defined TRACERS
      if (allocated(out_his_diags_pv_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_his_diags_pv_tracer)) deallocate(agrif_nml_storage(gid)%out_his_diags_pv_tracer)
         allocate(agrif_nml_storage(gid)%out_his_diags_pv_tracer(size(out_his_diags_pv_tracer)))
         agrif_nml_storage(gid)%out_his_diags_pv_tracer = out_his_diags_pv_tracer
      endif
# ifdef AVERAGES
      if (allocated(out_avg_diags_pv_tracer)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_diags_pv_tracer)) deallocate(agrif_nml_storage(gid)%out_avg_diags_pv_tracer)
         allocate(agrif_nml_storage(gid)%out_avg_diags_pv_tracer(size(out_avg_diags_pv_tracer)))
         agrif_nml_storage(gid)%out_avg_diags_pv_tracer = out_avg_diags_pv_tracer
      endif
# endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_diags_eddy = out_avg_diags_eddy
# endif
#endif
#ifdef DIAGNOSTICS_BIO
      if (allocated(out_his_diagbioFlux)) then
         if (allocated(agrif_nml_storage(gid)%out_his_diagbioFlux)) deallocate(agrif_nml_storage(gid)%out_his_diagbioFlux)
         allocate(agrif_nml_storage(gid)%out_his_diagbioFlux(size(out_his_diagbioFlux)))
         agrif_nml_storage(gid)%out_his_diagbioFlux = out_his_diagbioFlux
      endif
      if (allocated(out_his_diagbioVSink)) then
         if (allocated(agrif_nml_storage(gid)%out_his_diagbioVSink)) deallocate(agrif_nml_storage(gid)%out_his_diagbioVSink)
         allocate(agrif_nml_storage(gid)%out_his_diagbioVSink(size(out_his_diagbioVSink)))
         agrif_nml_storage(gid)%out_his_diagbioVSink = out_his_diagbioVSink
      endif
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      if (allocated(out_his_diagbioGasExc)) then
         if (allocated(agrif_nml_storage(gid)%out_his_diagbioGasExc)) deallocate(agrif_nml_storage(gid)%out_his_diagbioGasExc)
         allocate(agrif_nml_storage(gid)%out_his_diagbioGasExc(size(out_his_diagbioGasExc)))
         agrif_nml_storage(gid)%out_his_diagbioGasExc = out_his_diagbioGasExc
      endif
# endif
# ifdef AVERAGES
      if (allocated(out_avg_diagbioFlux)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_diagbioFlux)) deallocate(agrif_nml_storage(gid)%out_avg_diagbioFlux)
         allocate(agrif_nml_storage(gid)%out_avg_diagbioFlux(size(out_avg_diagbioFlux)))
         agrif_nml_storage(gid)%out_avg_diagbioFlux = out_avg_diagbioFlux
      endif
      if (allocated(out_avg_diagbioVSink)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_diagbioVSink)) deallocate(agrif_nml_storage(gid)%out_avg_diagbioVSink)
         allocate(agrif_nml_storage(gid)%out_avg_diagbioVSink(size(out_avg_diagbioVSink)))
         agrif_nml_storage(gid)%out_avg_diagbioVSink = out_avg_diagbioVSink
      endif
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      if (allocated(out_avg_diagbioGasExc)) then
         if (allocated(agrif_nml_storage(gid)%out_avg_diagbioGasExc)) deallocate(agrif_nml_storage(gid)%out_avg_diagbioGasExc)
         allocate(agrif_nml_storage(gid)%out_avg_diagbioGasExc(size(out_avg_diagbioGasExc)))
         agrif_nml_storage(gid)%out_avg_diagbioGasExc = out_avg_diagbioGasExc
      endif
#  endif
# endif
#endif
#ifdef ONLINE
      agrif_nml_storage(gid)%yearnum = yearnum
      agrif_nml_storage(gid)%monthnum = monthnum
      agrif_nml_storage(gid)%recordsperday = recordsperday
      agrif_nml_storage(gid)%yearend = yearend
      agrif_nml_storage(gid)%monthend = monthend
      agrif_nml_storage(gid)%pathbulk = pathbulk
#endif
#ifdef STOGEN
      agrif_nml_storage(gid)%out_his_xi2d = out_his_xi2d
      agrif_nml_storage(gid)%out_his_xi3d = out_his_xi3d
#endif
#if defined SOLVE3D && defined GLS_MIXING
      agrif_nml_storage(gid)%out_his_tke = out_his_tke
      agrif_nml_storage(gid)%out_his_gls = out_his_gls
      agrif_nml_storage(gid)%out_his_lscale = out_his_lscale
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_tke = out_avg_tke
      agrif_nml_storage(gid)%out_avg_gls = out_avg_gls
      agrif_nml_storage(gid)%out_avg_lscale = out_avg_lscale
# endif
#endif
#if defined SOLVE3D && defined SEDIMENT
      agrif_nml_storage(gid)%out_his_sed_athk = out_his_sed_athk
      agrif_nml_storage(gid)%out_his_sed_bthk = out_his_sed_bthk
      agrif_nml_storage(gid)%out_his_sed_bpor = out_his_sed_bpor
      if (allocated(out_his_sed_bfra)) then
         if (allocated(agrif_nml_storage(gid)%out_his_sed_bfra)) deallocate(agrif_nml_storage(gid)%out_his_sed_bfra)
         allocate(agrif_nml_storage(gid)%out_his_sed_bfra(size(out_his_sed_bfra)))
         agrif_nml_storage(gid)%out_his_sed_bfra = out_his_sed_bfra
      endif
# ifdef SUSPLOAD
      if (allocated(out_his_sed_dflx)) then
         if (allocated(agrif_nml_storage(gid)%out_his_sed_dflx)) deallocate(agrif_nml_storage(gid)%out_his_sed_dflx)
         allocate(agrif_nml_storage(gid)%out_his_sed_dflx(size(out_his_sed_dflx)))
         agrif_nml_storage(gid)%out_his_sed_dflx = out_his_sed_dflx
      endif
      if (allocated(out_his_sed_eflx)) then
         if (allocated(agrif_nml_storage(gid)%out_his_sed_eflx)) deallocate(agrif_nml_storage(gid)%out_his_sed_eflx)
         allocate(agrif_nml_storage(gid)%out_his_sed_eflx(size(out_his_sed_eflx)))
         agrif_nml_storage(gid)%out_his_sed_eflx = out_his_sed_eflx
      endif
# endif
# ifdef BEDLOAD
      if (allocated(out_his_sed_bdlu)) then
         if (allocated(agrif_nml_storage(gid)%out_his_sed_bdlu)) deallocate(agrif_nml_storage(gid)%out_his_sed_bdlu)
         allocate(agrif_nml_storage(gid)%out_his_sed_bdlu(size(out_his_sed_bdlu)))
         agrif_nml_storage(gid)%out_his_sed_bdlu = out_his_sed_bdlu
      endif
      if (allocated(out_his_sed_bdlv)) then
         if (allocated(agrif_nml_storage(gid)%out_his_sed_bdlv)) deallocate(agrif_nml_storage(gid)%out_his_sed_bdlv)
         allocate(agrif_nml_storage(gid)%out_his_sed_bdlv(size(out_his_sed_bdlv)))
         agrif_nml_storage(gid)%out_his_sed_bdlv = out_his_sed_bdlv
      endif
# endif
# if defined MIXED_BED || defined COHESIVE_BED
      agrif_nml_storage(gid)%out_his_sed_btcr = out_his_sed_btcr
# endif
#endif
#ifdef BBL
      agrif_nml_storage(gid)%out_his_abed = out_his_abed
      agrif_nml_storage(gid)%out_his_hripple = out_his_hripple
      agrif_nml_storage(gid)%out_his_lripple = out_his_lripple
      agrif_nml_storage(gid)%out_his_zbnot = out_his_zbnot
      agrif_nml_storage(gid)%out_his_zbapp = out_his_zbapp
      agrif_nml_storage(gid)%out_his_bostrw = out_his_bostrw
#endif
#ifdef MRL_WCI
      agrif_nml_storage(gid)%out_his_sup = out_his_sup
      agrif_nml_storage(gid)%out_his_ust2d = out_his_ust2d
      agrif_nml_storage(gid)%out_his_vst2d = out_his_vst2d
# ifdef SOLVE3D
      agrif_nml_storage(gid)%out_his_ust = out_his_ust
      agrif_nml_storage(gid)%out_his_vst = out_his_vst
      agrif_nml_storage(gid)%out_his_wst = out_his_wst
      agrif_nml_storage(gid)%out_his_akb = out_his_akb
      agrif_nml_storage(gid)%out_his_akw = out_his_akw
      agrif_nml_storage(gid)%out_his_kvf = out_his_kvf
      agrif_nml_storage(gid)%out_his_calp = out_his_calp
      agrif_nml_storage(gid)%out_his_kaps = out_his_kaps
# endif
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_sup = out_avg_sup
      agrif_nml_storage(gid)%out_avg_ust2d = out_avg_ust2d
      agrif_nml_storage(gid)%out_avg_vst2d = out_avg_vst2d
#  ifdef SOLVE3D
      agrif_nml_storage(gid)%out_avg_ust = out_avg_ust
      agrif_nml_storage(gid)%out_avg_vst = out_avg_vst
      agrif_nml_storage(gid)%out_avg_wst = out_avg_wst
      agrif_nml_storage(gid)%out_avg_akb = out_avg_akb
      agrif_nml_storage(gid)%out_avg_akw = out_avg_akw
      agrif_nml_storage(gid)%out_avg_kvf = out_avg_kvf
      agrif_nml_storage(gid)%out_avg_calp = out_avg_calp
      agrif_nml_storage(gid)%out_avg_kaps = out_avg_kaps
#  endif
# endif
#endif
#if defined MRL_WCI || defined OW_COUPLING
      agrif_nml_storage(gid)%out_his_hrm = out_his_hrm
      agrif_nml_storage(gid)%out_his_frq = out_his_frq
      agrif_nml_storage(gid)%out_his_action = out_his_action
      agrif_nml_storage(gid)%out_his_k_xi = out_his_k_xi
      agrif_nml_storage(gid)%out_his_k_eta = out_his_k_eta
      agrif_nml_storage(gid)%out_his_eps_b = out_his_eps_b
      agrif_nml_storage(gid)%out_his_eps_d = out_his_eps_d
      agrif_nml_storage(gid)%out_his_erol = out_his_erol
      agrif_nml_storage(gid)%out_his_eps_r = out_his_eps_r
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_hrm = out_avg_hrm
      agrif_nml_storage(gid)%out_avg_frq = out_avg_frq
      agrif_nml_storage(gid)%out_avg_action = out_avg_action
      agrif_nml_storage(gid)%out_avg_k_xi = out_avg_k_xi
      agrif_nml_storage(gid)%out_avg_k_eta = out_avg_k_eta
      agrif_nml_storage(gid)%out_avg_eps_b = out_avg_eps_b
      agrif_nml_storage(gid)%out_avg_eps_d = out_avg_eps_d
      agrif_nml_storage(gid)%out_avg_erol = out_avg_erol
      agrif_nml_storage(gid)%out_avg_eps_r = out_avg_eps_r
# endif
#endif
#if defined ABL1D && !defined XIOS
      agrif_nml_storage(gid)%out_his_abl_pu_dta = out_his_abl_pu_dta
      agrif_nml_storage(gid)%out_his_abl_pv_dta = out_his_abl_pv_dta
      agrif_nml_storage(gid)%out_his_abl_pt_dta = out_his_abl_pt_dta
      agrif_nml_storage(gid)%out_his_abl_pq_dta = out_his_abl_pq_dta
      agrif_nml_storage(gid)%out_his_abl_pgu_dta = out_his_abl_pgu_dta
      agrif_nml_storage(gid)%out_his_abl_pgv_dta = out_his_abl_pgv_dta
      agrif_nml_storage(gid)%out_his_abl_u_abl = out_his_abl_u_abl
      agrif_nml_storage(gid)%out_his_abl_v_abl = out_his_abl_v_abl
      agrif_nml_storage(gid)%out_his_abl_t_abl = out_his_abl_t_abl
      agrif_nml_storage(gid)%out_his_abl_q_abl = out_his_abl_q_abl
      agrif_nml_storage(gid)%out_his_abl_tke_abl = out_his_abl_tke_abl
      agrif_nml_storage(gid)%out_his_abl_mxlm_abl = out_his_abl_mxlm_abl
      agrif_nml_storage(gid)%out_his_abl_mxld_abl = out_his_abl_mxld_abl
      agrif_nml_storage(gid)%out_his_abl_avm_abl = out_his_abl_avm_abl
      agrif_nml_storage(gid)%out_his_abl_avt_abl = out_his_abl_avt_abl
      agrif_nml_storage(gid)%out_his_abl_ablh_abl = out_his_abl_ablh_abl
      agrif_nml_storage(gid)%out_his_abl_zr_abl = out_his_abl_zr_abl
      agrif_nml_storage(gid)%out_his_abl_zw_abl = out_his_abl_zw_abl
      agrif_nml_storage(gid)%out_his_abl_Hzr_abl = out_his_abl_Hzr_abl
      agrif_nml_storage(gid)%out_his_abl_Hzw_abl = out_his_abl_Hzw_abl
# ifdef AVERAGES
      agrif_nml_storage(gid)%out_avg_abl_pu_dta = out_avg_abl_pu_dta
      agrif_nml_storage(gid)%out_avg_abl_pv_dta = out_avg_abl_pv_dta
      agrif_nml_storage(gid)%out_avg_abl_pt_dta = out_avg_abl_pt_dta
      agrif_nml_storage(gid)%out_avg_abl_pq_dta = out_avg_abl_pq_dta
      agrif_nml_storage(gid)%out_avg_abl_pgu_dta = out_avg_abl_pgu_dta
      agrif_nml_storage(gid)%out_avg_abl_pgv_dta = out_avg_abl_pgv_dta
      agrif_nml_storage(gid)%out_avg_abl_u_abl = out_avg_abl_u_abl
      agrif_nml_storage(gid)%out_avg_abl_v_abl = out_avg_abl_v_abl
      agrif_nml_storage(gid)%out_avg_abl_t_abl = out_avg_abl_t_abl
      agrif_nml_storage(gid)%out_avg_abl_q_abl = out_avg_abl_q_abl
      agrif_nml_storage(gid)%out_avg_abl_tke_abl = out_avg_abl_tke_abl
      agrif_nml_storage(gid)%out_avg_abl_mxlm_abl = out_avg_abl_mxlm_abl
      agrif_nml_storage(gid)%out_avg_abl_mxld_abl = out_avg_abl_mxld_abl
      agrif_nml_storage(gid)%out_avg_abl_avm_abl = out_avg_abl_avm_abl
      agrif_nml_storage(gid)%out_avg_abl_avt_abl = out_avg_abl_avt_abl
      agrif_nml_storage(gid)%out_avg_abl_ablh_abl = out_avg_abl_ablh_abl
      agrif_nml_storage(gid)%out_avg_abl_zr_abl = out_avg_abl_zr_abl
      agrif_nml_storage(gid)%out_avg_abl_zw_abl = out_avg_abl_zw_abl
      agrif_nml_storage(gid)%out_avg_abl_Hzr_abl = out_avg_abl_Hzr_abl
      agrif_nml_storage(gid)%out_avg_abl_Hzw_abl = out_avg_abl_Hzw_abl
# endif
#endif
# if defined SOLVE3D || defined RIP
      agrif_nml_storage(gid)%out_his_rho = out_his_rho
      agrif_nml_storage(gid)%out_his_omega = out_his_omega
      agrif_nml_storage(gid)%out_his_w = out_his_w
      agrif_nml_storage(gid)%out_his_akv = out_his_akv
      agrif_nml_storage(gid)%out_his_bostr = out_his_bostr
      agrif_nml_storage(gid)%out_his_bustr = out_his_bustr
      agrif_nml_storage(gid)%out_his_bvstr = out_his_bvstr
      agrif_nml_storage(gid)%out_his_wstr = out_his_wstr
      agrif_nml_storage(gid)%out_his_ustr = out_his_ustr
      agrif_nml_storage(gid)%out_his_vstr = out_his_vstr
# endif
# if defined SOLVE3D && defined TEMPERATURE
      agrif_nml_storage(gid)%out_his_akt = out_his_akt
      agrif_nml_storage(gid)%out_his_shflx = out_his_shflx
      agrif_nml_storage(gid)%out_his_shflx_rsw = out_his_shflx_rsw
# endif
# if defined SOLVE3D && defined SALINITY
      agrif_nml_storage(gid)%out_his_aks = out_his_aks
      agrif_nml_storage(gid)%out_his_swflx = out_his_swflx
# endif
# if defined SOLVE3D && defined BULK_FLUX
      agrif_nml_storage(gid)%out_his_shflx_rlw = out_his_shflx_rlw
      agrif_nml_storage(gid)%out_his_shflx_lat = out_his_shflx_lat
      agrif_nml_storage(gid)%out_his_shflx_sen = out_his_shflx_sen
# endif
# if defined SOLVE3D && defined BHFLUX
      agrif_nml_storage(gid)%out_his_bhflx = out_his_bhflx
# endif
# if defined SOLVE3D && defined BWFLUX && defined SALINITY
      agrif_nml_storage(gid)%out_his_bwflx = out_his_bwflx
# endif
# if defined SOLVE3D && (defined ANA_VMIX||defined LMD_MIXING|| defined LMD_SKPP||defined LMD_BKPP|| defined GLS_MIXING)
      agrif_nml_storage(gid)%out_his_bvf = out_his_bvf
# endif
# if defined SOLVE3D && (defined LMD_SKPP || defined GLS_MIXING)
      agrif_nml_storage(gid)%out_his_hbl = out_his_hbl
# endif
# if defined SOLVE3D && defined LMD_BKPP
      agrif_nml_storage(gid)%out_his_hbbl = out_his_hbbl
# endif
# if defined SOLVE3D && defined VIS_COEF_3D
      agrif_nml_storage(gid)%out_his_visc3d = out_his_visc3d
# endif
# if defined SOLVE3D && defined DIF_COEF_3D
      agrif_nml_storage(gid)%out_his_diff3d = out_his_diff3d
# endif
# if defined SOLVE3D && defined BIOLOGY && !defined PISCES
      agrif_nml_storage(gid)%out_his_hel = out_his_hel
# endif
# if defined SOLVE3D && defined BIO_NChlPZD
      agrif_nml_storage(gid)%out_his_chc = out_his_chc
      agrif_nml_storage(gid)%out_his_u10 = out_his_u10
      agrif_nml_storage(gid)%out_his_kvo2 = out_his_kvo2
      agrif_nml_storage(gid)%out_his_o2sat = out_his_o2sat
# endif
# if defined SOLVE3D && defined BIO_BioEBUS
      agrif_nml_storage(gid)%out_his_aou = out_his_aou
      agrif_nml_storage(gid)%out_his_wind10 = out_his_wind10
# endif
# if defined SOLVE3D && defined MORPHODYN
      agrif_nml_storage(gid)%out_his_hm = out_his_hm
# endif
# if defined AVERAGES && defined SOLVE3D
      agrif_nml_storage(gid)%out_avg_rho = out_avg_rho
      agrif_nml_storage(gid)%out_avg_omega = out_avg_omega
      agrif_nml_storage(gid)%out_avg_w = out_avg_w
      agrif_nml_storage(gid)%out_avg_akv = out_avg_akv
      agrif_nml_storage(gid)%out_avg_bostr = out_avg_bostr
      agrif_nml_storage(gid)%out_avg_bustr = out_avg_bustr
      agrif_nml_storage(gid)%out_avg_bvstr = out_avg_bvstr
      agrif_nml_storage(gid)%out_avg_wstr = out_avg_wstr
      agrif_nml_storage(gid)%out_avg_ustr = out_avg_ustr
      agrif_nml_storage(gid)%out_avg_vstr = out_avg_vstr
# endif
# if defined AVERAGES && defined SOLVE3D && defined TEMPERATURE
      agrif_nml_storage(gid)%out_avg_akt = out_avg_akt
      agrif_nml_storage(gid)%out_avg_shflx = out_avg_shflx
      agrif_nml_storage(gid)%out_avg_shflx_rsw = out_avg_shflx_rsw
# endif
# if defined AVERAGES && defined SOLVE3D && defined SALINITY
      agrif_nml_storage(gid)%out_avg_aks = out_avg_aks
      agrif_nml_storage(gid)%out_avg_swflx = out_avg_swflx
# endif
# if defined AVERAGES && defined SOLVE3D && defined BULK_FLUX
      agrif_nml_storage(gid)%out_avg_shflx_rlw = out_avg_shflx_rlw
      agrif_nml_storage(gid)%out_avg_shflx_lat = out_avg_shflx_lat
      agrif_nml_storage(gid)%out_avg_shflx_sen = out_avg_shflx_sen
# endif
# if defined AVERAGES && defined SOLVE3D && defined BHFLUX
      agrif_nml_storage(gid)%out_avg_bhflx = out_avg_bhflx
# endif
# if defined AVERAGES && defined SOLVE3D && defined BWFLUX && defined SALINITY
      agrif_nml_storage(gid)%out_avg_bwflx = out_avg_bwflx
# endif
# if defined AVERAGES && defined SOLVE3D && (defined ANA_VMIX||defined LMD_MIXING|| defined LMD_SKPP||defined LMD_BKPP|| defined GLS_MIXING)
      agrif_nml_storage(gid)%out_avg_bvf = out_avg_bvf
# endif
# if defined AVERAGES && defined SOLVE3D && (defined LMD_SKPP||defined GLS_MIXING)
      agrif_nml_storage(gid)%out_avg_hbl = out_avg_hbl
# endif
# if defined AVERAGES && defined SOLVE3D && defined LMD_BKPP
      agrif_nml_storage(gid)%out_avg_hbbl = out_avg_hbbl
# endif
# if defined AVERAGES && defined SOLVE3D && defined VIS_COEF_3D
      agrif_nml_storage(gid)%out_avg_visc3d = out_avg_visc3d
# endif
# if defined AVERAGES && defined SOLVE3D && defined DIF_COEF_3D
      agrif_nml_storage(gid)%out_avg_diff3d = out_avg_diff3d
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIOLOGY && !defined PISCES
      agrif_nml_storage(gid)%out_avg_hel = out_avg_hel
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_NChlPZD
      agrif_nml_storage(gid)%out_avg_chc = out_avg_chc
      agrif_nml_storage(gid)%out_avg_u10 = out_avg_u10
      agrif_nml_storage(gid)%out_avg_kvo2 = out_avg_kvo2
      agrif_nml_storage(gid)%out_avg_o2sat = out_avg_o2sat
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_BioEBUS
      agrif_nml_storage(gid)%out_avg_aou = out_avg_aou
      agrif_nml_storage(gid)%out_avg_wind10 = out_avg_wind10
# endif
# if defined AVERAGES && defined SOLVE3D && defined MORPHODYN
      agrif_nml_storage(gid)%out_avg_hm = out_avg_hm
# endif
      agrif_nml_stored(gid) = .true.
      agrif_nml_current_grid = gid
   end subroutine croco_namelist_save_grid

   subroutine croco_namelist_restore_grid(gid_in)
      integer, intent(in) :: gid_in
      integer :: gid
      if (gid_in <= 0) then
         gid = 0
      else
         gid = gid_in
      endif
      if (gid > MAX_AGRIF_GRIDS) return
      if (.not. agrif_nml_stored(gid)) return
      fname_nml = agrif_nml_storage(gid)%fname_nml
      testcase_name = agrif_nml_storage(gid)%testcase_name
      lts_hadv = agrif_nml_storage(gid)%lts_hadv
      lts_rot_or_diag = agrif_nml_storage(gid)%lts_rot_or_diag
      lts_per = agrif_nml_storage(gid)%lts_per
      lekman_apg = agrif_nml_storage(gid)%lekman_apg
      loscnonrot = agrif_nml_storage(gid)%loscnonrot
      lnonrot = agrif_nml_storage(gid)%lnonrot
      UG_forced = agrif_nml_storage(gid)%UG_forced
      wave_ramp_duration = agrif_nml_storage(gid)%wave_ramp_duration
      title = agrif_nml_storage(gid)%title
      logname = agrif_nml_storage(gid)%logname
      dt = agrif_nml_storage(gid)%dt
      ntimes = agrif_nml_storage(gid)%ntimes
      ndtfast = agrif_nml_storage(gid)%ndtfast
      ninfo = agrif_nml_storage(gid)%ninfo
      m2filter_alpha = agrif_nml_storage(gid)%m2filter_alpha
      ldefhis = agrif_nml_storage(gid)%ldefhis
      nwrt = agrif_nml_storage(gid)%nwrt
      nrpfhis = agrif_nml_storage(gid)%nrpfhis
      hisname = agrif_nml_storage(gid)%hisname
      nrrec = agrif_nml_storage(gid)%nrrec
      ininame = agrif_nml_storage(gid)%ininame
      nrst = agrif_nml_storage(gid)%nrst
      nrpfrst = agrif_nml_storage(gid)%nrpfrst
      rstname = agrif_nml_storage(gid)%rstname
#if defined NBQ || defined K3FAST   /* H1-fix 2026-08-26: directive invalide (#ifdef NBQ || ...) */
      ndtnbq_in = agrif_nml_storage(gid)%ndtnbq_in
      csound_nbq = agrif_nml_storage(gid)%csound_nbq
      visc2read_nbq = agrif_nml_storage(gid)%visc2read_nbq
#endif
#ifdef SOLVE3D
      theta_s = agrif_nml_storage(gid)%theta_s
      theta_b = agrif_nml_storage(gid)%theta_b
      hc = agrif_nml_storage(gid)%hc
      Tcline = agrif_nml_storage(gid)%Tcline
#endif
#ifdef USE_CALENDAR
      start_date = agrif_nml_storage(gid)%start_date
      end_date = agrif_nml_storage(gid)%end_date
      dt_his = agrif_nml_storage(gid)%dt_his
      dt_avg = agrif_nml_storage(gid)%dt_avg
      dt_rst = agrif_nml_storage(gid)%dt_rst
#endif
#ifndef ANA_GRID
      grdname = agrif_nml_storage(gid)%grdname
#endif
      frcname = agrif_nml_storage(gid)%frcname
      use_frcname = agrif_nml_storage(gid)%use_frcname
      bulkname = agrif_nml_storage(gid)%bulkname
#if (  defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || (defined ZCLIMATOLOGY&& !defined ANA_SSH) || (defined M2CLIMATOLOGY&& !defined ANA_M2CLIMA) || (defined M3CLIMATOLOGY&& !defined ANA_M3CLIMA)
      clmname = agrif_nml_storage(gid)%clmname
#endif
#if !defined ANA_BRY && defined FRC_BRY
      bry_file = agrif_nml_storage(gid)%bry_file
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      brywkb_file = agrif_nml_storage(gid)%brywkb_file
#endif
#ifdef WKB_WWAVE
      wkb_amp = agrif_nml_storage(gid)%wkb_amp
      wkb_ang = agrif_nml_storage(gid)%wkb_ang
      wkb_prd = agrif_nml_storage(gid)%wkb_prd
      wkb_tide = agrif_nml_storage(gid)%wkb_tide
      wkb_btg = agrif_nml_storage(gid)%wkb_btg
      wkb_gam = agrif_nml_storage(gid)%wkb_gam
#  ifdef WAVE_ROLLER
      wkb_rsb = agrif_nml_storage(gid)%wkb_rsb
      wkb_roller = agrif_nml_storage(gid)%wkb_roller
#  endif
#endif
#ifdef WAVE_MAKER
      wmaker_amp = agrif_nml_storage(gid)%wmaker_amp
      wmaker_prd = agrif_nml_storage(gid)%wmaker_prd
      wmaker_dir = agrif_nml_storage(gid)%wmaker_dir
      wmaker_dsp = agrif_nml_storage(gid)%wmaker_dsp
      wmaker_fsp = agrif_nml_storage(gid)%wmaker_fsp
      wmaker_wf1 = agrif_nml_storage(gid)%wmaker_wf1
      wmaker_wf2 = agrif_nml_storage(gid)%wmaker_wf2
      wmaker_wa1 = agrif_nml_storage(gid)%wmaker_wa1
      wmaker_wa2 = agrif_nml_storage(gid)%wmaker_wa2
#endif
#ifdef AVERAGES
      ntsavg = agrif_nml_storage(gid)%ntsavg
      navg = agrif_nml_storage(gid)%navg
      nrpfavg = agrif_nml_storage(gid)%nrpfavg
      avgname = agrif_nml_storage(gid)%avgname
#endif
#if defined OUTPUTS_SURFACE && !defined XIOS
      ldefsurf = agrif_nml_storage(gid)%ldefsurf
      nwrtsurf = agrif_nml_storage(gid)%nwrtsurf
      nrpfsurf = agrif_nml_storage(gid)%nrpfsurf
      surfname = agrif_nml_storage(gid)%surfname
      out_his_surf = agrif_nml_storage(gid)%out_his_surf
#  ifdef AVERAGES
      ldefsurf_avg = agrif_nml_storage(gid)%ldefsurf_avg
      ntssurf_avg = agrif_nml_storage(gid)%ntssurf_avg
      nwrtsurf_avg = agrif_nml_storage(gid)%nwrtsurf_avg
      nrpfsurf_avg = agrif_nml_storage(gid)%nrpfsurf_avg
      surfname_avg = agrif_nml_storage(gid)%surfname_avg
      out_avg_surf = agrif_nml_storage(gid)%out_avg_surf
#  endif
#endif
#if defined ABL1D && !defined XIOS
      ldefablhis = agrif_nml_storage(gid)%ldefablhis
      nwrtablhis = agrif_nml_storage(gid)%nwrtablhis
      nrpfablhis = agrif_nml_storage(gid)%nrpfablhis
      ablname = agrif_nml_storage(gid)%ablname
#  ifdef AVERAGES
      ldefablavg = agrif_nml_storage(gid)%ldefablavg
      ntsablavg = agrif_nml_storage(gid)%ntsablavg
      nwrtablavg = agrif_nml_storage(gid)%nwrtablavg
      nrpfablavg = agrif_nml_storage(gid)%nrpfablavg
      ablname_avg = agrif_nml_storage(gid)%ablname_avg
#  endif
#endif
#if defined WAVE_OFFLINE && defined MUSTANG
      wave_file = agrif_nml_storage(gid)%wave_file
#endif
#if defined BIOLOGY && defined PISCES
      bioname = agrif_nml_storage(gid)%bioname
#endif
#ifdef BODYFORCE
      levsfrc = agrif_nml_storage(gid)%levsfrc
      levbfrc = agrif_nml_storage(gid)%levbfrc
#endif
      R0 = agrif_nml_storage(gid)%R0
      T0 = agrif_nml_storage(gid)%T0
      S0 = agrif_nml_storage(gid)%S0
      Tcoef = agrif_nml_storage(gid)%Tcoef
      Scoef = agrif_nml_storage(gid)%Scoef
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
      ltra_min = agrif_nml_storage(gid)%ltra_min
      ltra_max = agrif_nml_storage(gid)%ltra_max
#endif
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
      ldyn_min = agrif_nml_storage(gid)%ldyn_min
      ldyn_max = agrif_nml_storage(gid)%ldyn_max
#endif
#ifdef SEDIMENT
      sedname = agrif_nml_storage(gid)%sedname
#endif
#ifdef MUSTANG
      sedname_must = agrif_nml_storage(gid)%sedname_must
#endif
#ifdef SUBSTANCE
      subsfilename = agrif_nml_storage(gid)%subsfilename
#endif
#ifdef OBSTRUCTION
      obstname = agrif_nml_storage(gid)%obstname
#endif
#ifdef XIOS
      xios_origin_date = agrif_nml_storage(gid)%xios_origin_date
#endif
#ifdef ASSIMILATION
      aparnam = agrif_nml_storage(gid)%aparnam
      assname = agrif_nml_storage(gid)%assname
#endif
      rho0 = agrif_nml_storage(gid)%rho0
      rdrg = agrif_nml_storage(gid)%rdrg
      rdrg2 = agrif_nml_storage(gid)%rdrg2
      Zobt = agrif_nml_storage(gid)%Zobt
      Cdb_min = agrif_nml_storage(gid)%Cdb_min
      Cdb_max = agrif_nml_storage(gid)%Cdb_max
      gamma2 = agrif_nml_storage(gid)%gamma2
      visc2 = agrif_nml_storage(gid)%visc2
      visc4 = agrif_nml_storage(gid)%visc4
#ifdef SOLVE3D
#  ifdef TRACERS
      if (allocated(agrif_nml_storage(gid)%tnu2)) then
         if (allocated(tnu2)) deallocate(tnu2)
         allocate(tnu2(size(agrif_nml_storage(gid)%tnu2)))
         tnu2 = agrif_nml_storage(gid)%tnu2
      endif
      if (allocated(agrif_nml_storage(gid)%tnu4)) then
         if (allocated(tnu4)) deallocate(tnu4)
         allocate(tnu4(size(agrif_nml_storage(gid)%tnu4)))
         tnu4 = agrif_nml_storage(gid)%tnu4
      endif
#  endif
#  if !defined LMD_MIXING
      Akv_bak = agrif_nml_storage(gid)%Akv_bak
#    ifdef TRACERS
      if (allocated(agrif_nml_storage(gid)%Akt_bak)) then
         if (allocated(Akt_bak)) deallocate(Akt_bak)
         allocate(Akt_bak(size(agrif_nml_storage(gid)%Akt_bak)))
         Akt_bak = agrif_nml_storage(gid)%Akt_bak
      endif
#    endif
#  endif
#endif
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
      psource_Nsrc = agrif_nml_storage(gid)%psource_Nsrc
      if (allocated(agrif_nml_storage(gid)%psource_Isrc)) then
         if (allocated(psource_Isrc)) deallocate(psource_Isrc)
         allocate(psource_Isrc(size(agrif_nml_storage(gid)%psource_Isrc)))
         psource_Isrc = agrif_nml_storage(gid)%psource_Isrc
      endif
      if (allocated(agrif_nml_storage(gid)%psource_Jsrc)) then
         if (allocated(psource_Jsrc)) deallocate(psource_Jsrc)
         allocate(psource_Jsrc(size(agrif_nml_storage(gid)%psource_Jsrc)))
         psource_Jsrc = agrif_nml_storage(gid)%psource_Jsrc
      endif
      if (allocated(agrif_nml_storage(gid)%psource_Dsrc)) then
         if (allocated(psource_Dsrc)) deallocate(psource_Dsrc)
         allocate(psource_Dsrc(size(agrif_nml_storage(gid)%psource_Dsrc)))
         psource_Dsrc = agrif_nml_storage(gid)%psource_Dsrc
      endif
      psource_qbarname = agrif_nml_storage(gid)%psource_qbarname
      if (allocated(agrif_nml_storage(gid)%psource_qbardir)) then
         if (allocated(psource_qbardir)) deallocate(psource_qbardir)
         allocate(psource_qbardir(size(agrif_nml_storage(gid)%psource_qbardir)))
         psource_qbardir = agrif_nml_storage(gid)%psource_qbardir
      endif
      if (allocated(agrif_nml_storage(gid)%psource_Qbar)) then
         if (allocated(psource_Qbar)) deallocate(psource_Qbar)
         allocate(psource_Qbar(size(agrif_nml_storage(gid)%psource_Qbar)))
         psource_Qbar = agrif_nml_storage(gid)%psource_Qbar
      endif
# if defined TRACERS
      if (allocated(agrif_nml_storage(gid)%psource_Lsrc)) then
         if (allocated(psource_Lsrc)) deallocate(psource_Lsrc)
         allocate(psource_Lsrc(size(agrif_nml_storage(gid)%psource_Lsrc,1), size(agrif_nml_storage(gid)%psource_Lsrc,2)))
         psource_Lsrc = agrif_nml_storage(gid)%psource_Lsrc
      endif
      if (allocated(agrif_nml_storage(gid)%psource_Tsrc0)) then
         if (allocated(psource_Tsrc0)) deallocate(psource_Tsrc0)
         allocate(psource_Tsrc0(size(agrif_nml_storage(gid)%psource_Tsrc0,1), size(agrif_nml_storage(gid)%psource_Tsrc0,2)))
         psource_Tsrc0 = agrif_nml_storage(gid)%psource_Tsrc0
      endif
# endif
#endif
#if defined DIAGNOSTICS_TS
      ldefdia = agrif_nml_storage(gid)%ldefdia
      nwrtdia = agrif_nml_storage(gid)%nwrtdia
      nrpfdia = agrif_nml_storage(gid)%nrpfdia
      dianame = agrif_nml_storage(gid)%dianame
#  ifdef AVERAGES
      ldefdia_avg = agrif_nml_storage(gid)%ldefdia_avg
      ntsdia_avg = agrif_nml_storage(gid)%ntsdia_avg
      nwrtdia_avg = agrif_nml_storage(gid)%nwrtdia_avg
      nrpfdia_avg = agrif_nml_storage(gid)%nrpfdia_avg
      dianame_avg = agrif_nml_storage(gid)%dianame_avg
#  endif
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_DENS
      mld_crit_D = agrif_nml_storage(gid)%mld_crit_D
      mld_crit_T = agrif_nml_storage(gid)%mld_crit_T
#  endif
#endif
#if defined DIAGNOSTICS_UV
      ldefdiaM = agrif_nml_storage(gid)%ldefdiaM
      nwrtdiaM = agrif_nml_storage(gid)%nwrtdiaM
      nrpfdiaM = agrif_nml_storage(gid)%nrpfdiaM
      dianameM = agrif_nml_storage(gid)%dianameM
#  ifdef AVERAGES
      ldefdiaM_avg = agrif_nml_storage(gid)%ldefdiaM_avg
      ntsdiaM_avg = agrif_nml_storage(gid)%ntsdiaM_avg
      nwrtdiaM_avg = agrif_nml_storage(gid)%nwrtdiaM_avg
      nrpfdiaM_avg = agrif_nml_storage(gid)%nrpfdiaM_avg
      dianameM_avg = agrif_nml_storage(gid)%dianameM_avg
#  endif
#endif
#ifdef DIAGNOSTICS_VRT
      ldefdiags_vrt = agrif_nml_storage(gid)%ldefdiags_vrt
      nwrtdiags_vrt = agrif_nml_storage(gid)%nwrtdiags_vrt
      nrpfdiags_vrt = agrif_nml_storage(gid)%nrpfdiags_vrt
      diags_vrtname = agrif_nml_storage(gid)%diags_vrtname
#  ifdef AVERAGES
      ldefdiags_vrt_avg = agrif_nml_storage(gid)%ldefdiags_vrt_avg
      ntsdiags_vrt_avg = agrif_nml_storage(gid)%ntsdiags_vrt_avg
      nwrtdiags_vrt_avg = agrif_nml_storage(gid)%nwrtdiags_vrt_avg
      nrpfdiags_vrt_avg = agrif_nml_storage(gid)%nrpfdiags_vrt_avg
      diags_vrtname_avg = agrif_nml_storage(gid)%diags_vrtname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_KE
      ldefdiags_ek = agrif_nml_storage(gid)%ldefdiags_ek
      nwrtdiags_ek = agrif_nml_storage(gid)%nwrtdiags_ek
      nrpfdiags_ek = agrif_nml_storage(gid)%nrpfdiags_ek
      diags_ekname = agrif_nml_storage(gid)%diags_ekname
#  ifdef AVERAGES
      ldefdiags_ek_avg = agrif_nml_storage(gid)%ldefdiags_ek_avg
      ntsdiags_ek_avg = agrif_nml_storage(gid)%ntsdiags_ek_avg
      nwrtdiags_ek_avg = agrif_nml_storage(gid)%nwrtdiags_ek_avg
      nrpfdiags_ek_avg = agrif_nml_storage(gid)%nrpfdiags_ek_avg
      diags_ekname_avg = agrif_nml_storage(gid)%diags_ekname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_PV
      ldefdiags_pv = agrif_nml_storage(gid)%ldefdiags_pv
      nwrtdiags_pv = agrif_nml_storage(gid)%nwrtdiags_pv
      nrpfdiags_pv = agrif_nml_storage(gid)%nrpfdiags_pv
      diags_pvname = agrif_nml_storage(gid)%diags_pvname
#  ifdef AVERAGES
      ldefdiags_pv_avg = agrif_nml_storage(gid)%ldefdiags_pv_avg
      ntsdiags_pv_avg = agrif_nml_storage(gid)%ntsdiags_pv_avg
      nwrtdiags_pv_avg = agrif_nml_storage(gid)%nwrtdiags_pv_avg
      nrpfdiags_pv_avg = agrif_nml_storage(gid)%nrpfdiags_pv_avg
      diags_pvname_avg = agrif_nml_storage(gid)%diags_pvname_avg
#  endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
#  ifdef AVERAGES
      ldefdiags_eddy_avg = agrif_nml_storage(gid)%ldefdiags_eddy_avg
      ntsdiags_eddy_avg = agrif_nml_storage(gid)%ntsdiags_eddy_avg
      nwrtdiags_eddy_avg = agrif_nml_storage(gid)%nwrtdiags_eddy_avg
      nrpfdiags_eddy_avg = agrif_nml_storage(gid)%nrpfdiags_eddy_avg
      diags_eddyname_avg = agrif_nml_storage(gid)%diags_eddyname_avg
#  endif
#endif
#ifdef DIAGNOSTICS_BIO
      ldefdiabio = agrif_nml_storage(gid)%ldefdiabio
      nwrtdiabio = agrif_nml_storage(gid)%nwrtdiabio
      nrpfdiabio = agrif_nml_storage(gid)%nrpfdiabio
      dianamebio = agrif_nml_storage(gid)%dianamebio
#  ifdef AVERAGES
      ldefdiabio_avg = agrif_nml_storage(gid)%ldefdiabio_avg
      ntsdiabio_avg = agrif_nml_storage(gid)%ntsdiabio_avg
      nwrtdiabio_avg = agrif_nml_storage(gid)%nwrtdiabio_avg
      nrpfdiabio_avg = agrif_nml_storage(gid)%nrpfdiabio_avg
      dianamebio_avg = agrif_nml_storage(gid)%dianamebio_avg
#  endif
#endif
#ifdef STATIONS
      ldefsta = agrif_nml_storage(gid)%ldefsta
      nsta = agrif_nml_storage(gid)%nsta
      nrpfsta = agrif_nml_storage(gid)%nrpfsta
      staposname = agrif_nml_storage(gid)%staposname
      staname = agrif_nml_storage(gid)%staname
      sta_grd = agrif_nml_storage(gid)%sta_grd
      sta_temp = agrif_nml_storage(gid)%sta_temp
      sta_salt = agrif_nml_storage(gid)%sta_salt
      sta_rho = agrif_nml_storage(gid)%sta_rho
      sta_vel = agrif_nml_storage(gid)%sta_vel
#endif
#if (defined SPONGE && !defined SPONGE_GRID)
      x_sponge = agrif_nml_storage(gid)%x_sponge
#endif
#ifdef WET_DRY
      D_wetdry = agrif_nml_storage(gid)%D_wetdry
#endif
#if defined MRL_WCI && defined WAVE_DRY
      D_wavedry = agrif_nml_storage(gid)%D_wavedry
#endif
#if defined T_FRC_BRY     || defined M2_FRC_BRY    || defined M3_FRC_BRY||defined Z_FRC_BRY|| defined W_FRC_BRY||defined NBQ_FRC_BRY|| defined TCLIMATOLOGY||defined M2CLIMATOLOGY|| defined M3CLIMATOLOGY||defined ZCLIMATOLOGY|| defined WCLIMATOLOGY||defined NBQCLIMATOLOGY|| defined TNUDGING||defined ZNUDGING|| defined M2NUDGING||defined M3NUDGING|| defined NBQ_NUDGING||defined SPONGE
      tauT_in = agrif_nml_storage(gid)%tauT_in
      tauT_out = agrif_nml_storage(gid)%tauT_out
      tauM_in = agrif_nml_storage(gid)%tauM_in
      tauM_out = agrif_nml_storage(gid)%tauM_out
#endif
#if defined BHFLUX || (defined BWFLUX && defined SALINITY)
      btfname = agrif_nml_storage(gid)%btfname
#endif
      out_his_zeta = agrif_nml_storage(gid)%out_his_zeta
      out_his_ubar = agrif_nml_storage(gid)%out_his_ubar
      out_his_vbar = agrif_nml_storage(gid)%out_his_vbar
#ifdef SOLVE3D
      out_his_u = agrif_nml_storage(gid)%out_his_u
      out_his_v = agrif_nml_storage(gid)%out_his_v
# ifdef TRACERS
      if (allocated(agrif_nml_storage(gid)%out_his_tracer)) then
         if (allocated(out_his_tracer)) deallocate(out_his_tracer)
         allocate(out_his_tracer(size(agrif_nml_storage(gid)%out_his_tracer)))
         out_his_tracer = agrif_nml_storage(gid)%out_his_tracer
      endif
# endif
#endif
#ifdef AVERAGES
      out_avg_zeta = agrif_nml_storage(gid)%out_avg_zeta
      out_avg_ubar = agrif_nml_storage(gid)%out_avg_ubar
      out_avg_vbar = agrif_nml_storage(gid)%out_avg_vbar
# ifdef SOLVE3D
      out_avg_u = agrif_nml_storage(gid)%out_avg_u
      out_avg_v = agrif_nml_storage(gid)%out_avg_v
#  ifdef TRACERS
      if (allocated(agrif_nml_storage(gid)%out_avg_tracer)) then
         if (allocated(out_avg_tracer)) deallocate(out_avg_tracer)
         allocate(out_avg_tracer(size(agrif_nml_storage(gid)%out_avg_tracer)))
         out_avg_tracer = agrif_nml_storage(gid)%out_avg_tracer
      endif
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_TS
# ifdef TRACERS
      if (allocated(agrif_nml_storage(gid)%out_his_dia3D_tracer)) then
         if (allocated(out_his_dia3D_tracer)) deallocate(out_his_dia3D_tracer)
         allocate(out_his_dia3D_tracer(size(agrif_nml_storage(gid)%out_his_dia3D_tracer)))
         out_his_dia3D_tracer = agrif_nml_storage(gid)%out_his_dia3D_tracer
      endif
#  ifdef DIAGNOSTICS_TS_MLD
      if (allocated(agrif_nml_storage(gid)%out_his_dia2D_tracer)) then
         if (allocated(out_his_dia2D_tracer)) deallocate(out_his_dia2D_tracer)
         allocate(out_his_dia2D_tracer(size(agrif_nml_storage(gid)%out_his_dia2D_tracer)))
         out_his_dia2D_tracer = agrif_nml_storage(gid)%out_his_dia2D_tracer
      endif
#  endif
# endif
# if defined AVERAGES && defined TRACERS
      if (allocated(agrif_nml_storage(gid)%out_avg_dia3D_tracer)) then
         if (allocated(out_avg_dia3D_tracer)) deallocate(out_avg_dia3D_tracer)
         allocate(out_avg_dia3D_tracer(size(agrif_nml_storage(gid)%out_avg_dia3D_tracer)))
         out_avg_dia3D_tracer = agrif_nml_storage(gid)%out_avg_dia3D_tracer
      endif
#  ifdef DIAGNOSTICS_TS_MLD
      if (allocated(agrif_nml_storage(gid)%out_avg_dia2D_tracer)) then
         if (allocated(out_avg_dia2D_tracer)) deallocate(out_avg_dia2D_tracer)
         allocate(out_avg_dia2D_tracer(size(agrif_nml_storage(gid)%out_avg_dia2D_tracer)))
         out_avg_dia2D_tracer = agrif_nml_storage(gid)%out_avg_dia2D_tracer
      endif
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_UV
      out_his_diagM_u = agrif_nml_storage(gid)%out_his_diagM_u
      out_his_diagM_v = agrif_nml_storage(gid)%out_his_diagM_v
# ifdef AVERAGES
      out_avg_diagM_u = agrif_nml_storage(gid)%out_avg_diagM_u
      out_avg_diagM_v = agrif_nml_storage(gid)%out_avg_diagM_v
# endif
#endif
#ifdef DIAGNOSTICS_VRT
      out_his_diags_vrt = agrif_nml_storage(gid)%out_his_diags_vrt
# ifdef AVERAGES
      out_avg_diags_vrt = agrif_nml_storage(gid)%out_avg_diags_vrt
# endif
#endif
#ifdef DIAGNOSTICS_KE
      out_his_diags_ek = agrif_nml_storage(gid)%out_his_diags_ek
# ifdef AVERAGES
      out_avg_diags_ek = agrif_nml_storage(gid)%out_avg_diags_ek
# endif
#endif
#if defined DIAGNOSTICS_PV && defined TRACERS
      if (allocated(agrif_nml_storage(gid)%out_his_diags_pv_tracer)) then
         if (allocated(out_his_diags_pv_tracer)) deallocate(out_his_diags_pv_tracer)
         allocate(out_his_diags_pv_tracer(size(agrif_nml_storage(gid)%out_his_diags_pv_tracer)))
         out_his_diags_pv_tracer = agrif_nml_storage(gid)%out_his_diags_pv_tracer
      endif
# ifdef AVERAGES
      if (allocated(agrif_nml_storage(gid)%out_avg_diags_pv_tracer)) then
         if (allocated(out_avg_diags_pv_tracer)) deallocate(out_avg_diags_pv_tracer)
         allocate(out_avg_diags_pv_tracer(size(agrif_nml_storage(gid)%out_avg_diags_pv_tracer)))
         out_avg_diags_pv_tracer = agrif_nml_storage(gid)%out_avg_diags_pv_tracer
      endif
# endif
#endif
#if defined DIAGNOSTICS_EDDY && !defined XIOS
# ifdef AVERAGES
      out_avg_diags_eddy = agrif_nml_storage(gid)%out_avg_diags_eddy
# endif
#endif
#ifdef DIAGNOSTICS_BIO
      if (allocated(agrif_nml_storage(gid)%out_his_diagbioFlux)) then
         if (allocated(out_his_diagbioFlux)) deallocate(out_his_diagbioFlux)
         allocate(out_his_diagbioFlux(size(agrif_nml_storage(gid)%out_his_diagbioFlux)))
         out_his_diagbioFlux = agrif_nml_storage(gid)%out_his_diagbioFlux
      endif
      if (allocated(agrif_nml_storage(gid)%out_his_diagbioVSink)) then
         if (allocated(out_his_diagbioVSink)) deallocate(out_his_diagbioVSink)
         allocate(out_his_diagbioVSink(size(agrif_nml_storage(gid)%out_his_diagbioVSink)))
         out_his_diagbioVSink = agrif_nml_storage(gid)%out_his_diagbioVSink
      endif
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      if (allocated(agrif_nml_storage(gid)%out_his_diagbioGasExc)) then
         if (allocated(out_his_diagbioGasExc)) deallocate(out_his_diagbioGasExc)
         allocate(out_his_diagbioGasExc(size(agrif_nml_storage(gid)%out_his_diagbioGasExc)))
         out_his_diagbioGasExc = agrif_nml_storage(gid)%out_his_diagbioGasExc
      endif
# endif
# ifdef AVERAGES
      if (allocated(agrif_nml_storage(gid)%out_avg_diagbioFlux)) then
         if (allocated(out_avg_diagbioFlux)) deallocate(out_avg_diagbioFlux)
         allocate(out_avg_diagbioFlux(size(agrif_nml_storage(gid)%out_avg_diagbioFlux)))
         out_avg_diagbioFlux = agrif_nml_storage(gid)%out_avg_diagbioFlux
      endif
      if (allocated(agrif_nml_storage(gid)%out_avg_diagbioVSink)) then
         if (allocated(out_avg_diagbioVSink)) deallocate(out_avg_diagbioVSink)
         allocate(out_avg_diagbioVSink(size(agrif_nml_storage(gid)%out_avg_diagbioVSink)))
         out_avg_diagbioVSink = agrif_nml_storage(gid)%out_avg_diagbioVSink
      endif
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      if (allocated(agrif_nml_storage(gid)%out_avg_diagbioGasExc)) then
         if (allocated(out_avg_diagbioGasExc)) deallocate(out_avg_diagbioGasExc)
         allocate(out_avg_diagbioGasExc(size(agrif_nml_storage(gid)%out_avg_diagbioGasExc)))
         out_avg_diagbioGasExc = agrif_nml_storage(gid)%out_avg_diagbioGasExc
      endif
#  endif
# endif
#endif
#ifdef ONLINE
      yearnum = agrif_nml_storage(gid)%yearnum
      monthnum = agrif_nml_storage(gid)%monthnum
      recordsperday = agrif_nml_storage(gid)%recordsperday
      yearend = agrif_nml_storage(gid)%yearend
      monthend = agrif_nml_storage(gid)%monthend
      pathbulk = agrif_nml_storage(gid)%pathbulk
#endif
#ifdef STOGEN
      out_his_xi2d = agrif_nml_storage(gid)%out_his_xi2d
      out_his_xi3d = agrif_nml_storage(gid)%out_his_xi3d
#endif
#if defined SOLVE3D && defined GLS_MIXING
      out_his_tke = agrif_nml_storage(gid)%out_his_tke
      out_his_gls = agrif_nml_storage(gid)%out_his_gls
      out_his_lscale = agrif_nml_storage(gid)%out_his_lscale
# ifdef AVERAGES
      out_avg_tke = agrif_nml_storage(gid)%out_avg_tke
      out_avg_gls = agrif_nml_storage(gid)%out_avg_gls
      out_avg_lscale = agrif_nml_storage(gid)%out_avg_lscale
# endif
#endif
#if defined SOLVE3D && defined SEDIMENT
      out_his_sed_athk = agrif_nml_storage(gid)%out_his_sed_athk
      out_his_sed_bthk = agrif_nml_storage(gid)%out_his_sed_bthk
      out_his_sed_bpor = agrif_nml_storage(gid)%out_his_sed_bpor
      if (allocated(agrif_nml_storage(gid)%out_his_sed_bfra)) then
         if (allocated(out_his_sed_bfra)) deallocate(out_his_sed_bfra)
         allocate(out_his_sed_bfra(size(agrif_nml_storage(gid)%out_his_sed_bfra)))
         out_his_sed_bfra = agrif_nml_storage(gid)%out_his_sed_bfra
      endif
# ifdef SUSPLOAD
      if (allocated(agrif_nml_storage(gid)%out_his_sed_dflx)) then
         if (allocated(out_his_sed_dflx)) deallocate(out_his_sed_dflx)
         allocate(out_his_sed_dflx(size(agrif_nml_storage(gid)%out_his_sed_dflx)))
         out_his_sed_dflx = agrif_nml_storage(gid)%out_his_sed_dflx
      endif
      if (allocated(agrif_nml_storage(gid)%out_his_sed_eflx)) then
         if (allocated(out_his_sed_eflx)) deallocate(out_his_sed_eflx)
         allocate(out_his_sed_eflx(size(agrif_nml_storage(gid)%out_his_sed_eflx)))
         out_his_sed_eflx = agrif_nml_storage(gid)%out_his_sed_eflx
      endif
# endif
# ifdef BEDLOAD
      if (allocated(agrif_nml_storage(gid)%out_his_sed_bdlu)) then
         if (allocated(out_his_sed_bdlu)) deallocate(out_his_sed_bdlu)
         allocate(out_his_sed_bdlu(size(agrif_nml_storage(gid)%out_his_sed_bdlu)))
         out_his_sed_bdlu = agrif_nml_storage(gid)%out_his_sed_bdlu
      endif
      if (allocated(agrif_nml_storage(gid)%out_his_sed_bdlv)) then
         if (allocated(out_his_sed_bdlv)) deallocate(out_his_sed_bdlv)
         allocate(out_his_sed_bdlv(size(agrif_nml_storage(gid)%out_his_sed_bdlv)))
         out_his_sed_bdlv = agrif_nml_storage(gid)%out_his_sed_bdlv
      endif
# endif
# if defined MIXED_BED || defined COHESIVE_BED
      out_his_sed_btcr = agrif_nml_storage(gid)%out_his_sed_btcr
# endif
#endif
#ifdef BBL
      out_his_abed = agrif_nml_storage(gid)%out_his_abed
      out_his_hripple = agrif_nml_storage(gid)%out_his_hripple
      out_his_lripple = agrif_nml_storage(gid)%out_his_lripple
      out_his_zbnot = agrif_nml_storage(gid)%out_his_zbnot
      out_his_zbapp = agrif_nml_storage(gid)%out_his_zbapp
      out_his_bostrw = agrif_nml_storage(gid)%out_his_bostrw
#endif
#ifdef MRL_WCI
      out_his_sup = agrif_nml_storage(gid)%out_his_sup
      out_his_ust2d = agrif_nml_storage(gid)%out_his_ust2d
      out_his_vst2d = agrif_nml_storage(gid)%out_his_vst2d
# ifdef SOLVE3D
      out_his_ust = agrif_nml_storage(gid)%out_his_ust
      out_his_vst = agrif_nml_storage(gid)%out_his_vst
      out_his_wst = agrif_nml_storage(gid)%out_his_wst
      out_his_akb = agrif_nml_storage(gid)%out_his_akb
      out_his_akw = agrif_nml_storage(gid)%out_his_akw
      out_his_kvf = agrif_nml_storage(gid)%out_his_kvf
      out_his_calp = agrif_nml_storage(gid)%out_his_calp
      out_his_kaps = agrif_nml_storage(gid)%out_his_kaps
# endif
# ifdef AVERAGES
      out_avg_sup = agrif_nml_storage(gid)%out_avg_sup
      out_avg_ust2d = agrif_nml_storage(gid)%out_avg_ust2d
      out_avg_vst2d = agrif_nml_storage(gid)%out_avg_vst2d
#  ifdef SOLVE3D
      out_avg_ust = agrif_nml_storage(gid)%out_avg_ust
      out_avg_vst = agrif_nml_storage(gid)%out_avg_vst
      out_avg_wst = agrif_nml_storage(gid)%out_avg_wst
      out_avg_akb = agrif_nml_storage(gid)%out_avg_akb
      out_avg_akw = agrif_nml_storage(gid)%out_avg_akw
      out_avg_kvf = agrif_nml_storage(gid)%out_avg_kvf
      out_avg_calp = agrif_nml_storage(gid)%out_avg_calp
      out_avg_kaps = agrif_nml_storage(gid)%out_avg_kaps
#  endif
# endif
#endif
#if defined MRL_WCI || defined OW_COUPLING
      out_his_hrm = agrif_nml_storage(gid)%out_his_hrm
      out_his_frq = agrif_nml_storage(gid)%out_his_frq
      out_his_action = agrif_nml_storage(gid)%out_his_action
      out_his_k_xi = agrif_nml_storage(gid)%out_his_k_xi
      out_his_k_eta = agrif_nml_storage(gid)%out_his_k_eta
      out_his_eps_b = agrif_nml_storage(gid)%out_his_eps_b
      out_his_eps_d = agrif_nml_storage(gid)%out_his_eps_d
      out_his_erol = agrif_nml_storage(gid)%out_his_erol
      out_his_eps_r = agrif_nml_storage(gid)%out_his_eps_r
# ifdef AVERAGES
      out_avg_hrm = agrif_nml_storage(gid)%out_avg_hrm
      out_avg_frq = agrif_nml_storage(gid)%out_avg_frq
      out_avg_action = agrif_nml_storage(gid)%out_avg_action
      out_avg_k_xi = agrif_nml_storage(gid)%out_avg_k_xi
      out_avg_k_eta = agrif_nml_storage(gid)%out_avg_k_eta
      out_avg_eps_b = agrif_nml_storage(gid)%out_avg_eps_b
      out_avg_eps_d = agrif_nml_storage(gid)%out_avg_eps_d
      out_avg_erol = agrif_nml_storage(gid)%out_avg_erol
      out_avg_eps_r = agrif_nml_storage(gid)%out_avg_eps_r
# endif
#endif
#if defined ABL1D && !defined XIOS
      out_his_abl_pu_dta = agrif_nml_storage(gid)%out_his_abl_pu_dta
      out_his_abl_pv_dta = agrif_nml_storage(gid)%out_his_abl_pv_dta
      out_his_abl_pt_dta = agrif_nml_storage(gid)%out_his_abl_pt_dta
      out_his_abl_pq_dta = agrif_nml_storage(gid)%out_his_abl_pq_dta
      out_his_abl_pgu_dta = agrif_nml_storage(gid)%out_his_abl_pgu_dta
      out_his_abl_pgv_dta = agrif_nml_storage(gid)%out_his_abl_pgv_dta
      out_his_abl_u_abl = agrif_nml_storage(gid)%out_his_abl_u_abl
      out_his_abl_v_abl = agrif_nml_storage(gid)%out_his_abl_v_abl
      out_his_abl_t_abl = agrif_nml_storage(gid)%out_his_abl_t_abl
      out_his_abl_q_abl = agrif_nml_storage(gid)%out_his_abl_q_abl
      out_his_abl_tke_abl = agrif_nml_storage(gid)%out_his_abl_tke_abl
      out_his_abl_mxlm_abl = agrif_nml_storage(gid)%out_his_abl_mxlm_abl
      out_his_abl_mxld_abl = agrif_nml_storage(gid)%out_his_abl_mxld_abl
      out_his_abl_avm_abl = agrif_nml_storage(gid)%out_his_abl_avm_abl
      out_his_abl_avt_abl = agrif_nml_storage(gid)%out_his_abl_avt_abl
      out_his_abl_ablh_abl = agrif_nml_storage(gid)%out_his_abl_ablh_abl
      out_his_abl_zr_abl = agrif_nml_storage(gid)%out_his_abl_zr_abl
      out_his_abl_zw_abl = agrif_nml_storage(gid)%out_his_abl_zw_abl
      out_his_abl_Hzr_abl = agrif_nml_storage(gid)%out_his_abl_Hzr_abl
      out_his_abl_Hzw_abl = agrif_nml_storage(gid)%out_his_abl_Hzw_abl
# ifdef AVERAGES
      out_avg_abl_pu_dta = agrif_nml_storage(gid)%out_avg_abl_pu_dta
      out_avg_abl_pv_dta = agrif_nml_storage(gid)%out_avg_abl_pv_dta
      out_avg_abl_pt_dta = agrif_nml_storage(gid)%out_avg_abl_pt_dta
      out_avg_abl_pq_dta = agrif_nml_storage(gid)%out_avg_abl_pq_dta
      out_avg_abl_pgu_dta = agrif_nml_storage(gid)%out_avg_abl_pgu_dta
      out_avg_abl_pgv_dta = agrif_nml_storage(gid)%out_avg_abl_pgv_dta
      out_avg_abl_u_abl = agrif_nml_storage(gid)%out_avg_abl_u_abl
      out_avg_abl_v_abl = agrif_nml_storage(gid)%out_avg_abl_v_abl
      out_avg_abl_t_abl = agrif_nml_storage(gid)%out_avg_abl_t_abl
      out_avg_abl_q_abl = agrif_nml_storage(gid)%out_avg_abl_q_abl
      out_avg_abl_tke_abl = agrif_nml_storage(gid)%out_avg_abl_tke_abl
      out_avg_abl_mxlm_abl = agrif_nml_storage(gid)%out_avg_abl_mxlm_abl
      out_avg_abl_mxld_abl = agrif_nml_storage(gid)%out_avg_abl_mxld_abl
      out_avg_abl_avm_abl = agrif_nml_storage(gid)%out_avg_abl_avm_abl
      out_avg_abl_avt_abl = agrif_nml_storage(gid)%out_avg_abl_avt_abl
      out_avg_abl_ablh_abl = agrif_nml_storage(gid)%out_avg_abl_ablh_abl
      out_avg_abl_zr_abl = agrif_nml_storage(gid)%out_avg_abl_zr_abl
      out_avg_abl_zw_abl = agrif_nml_storage(gid)%out_avg_abl_zw_abl
      out_avg_abl_Hzr_abl = agrif_nml_storage(gid)%out_avg_abl_Hzr_abl
      out_avg_abl_Hzw_abl = agrif_nml_storage(gid)%out_avg_abl_Hzw_abl
# endif
#endif
# if defined SOLVE3D || defined RIP
      out_his_rho = agrif_nml_storage(gid)%out_his_rho
      out_his_omega = agrif_nml_storage(gid)%out_his_omega
      out_his_w = agrif_nml_storage(gid)%out_his_w
      out_his_akv = agrif_nml_storage(gid)%out_his_akv
      out_his_bostr = agrif_nml_storage(gid)%out_his_bostr
      out_his_bustr = agrif_nml_storage(gid)%out_his_bustr
      out_his_bvstr = agrif_nml_storage(gid)%out_his_bvstr
      out_his_wstr = agrif_nml_storage(gid)%out_his_wstr
      out_his_ustr = agrif_nml_storage(gid)%out_his_ustr
      out_his_vstr = agrif_nml_storage(gid)%out_his_vstr
# endif
# if defined SOLVE3D && defined TEMPERATURE
      out_his_akt = agrif_nml_storage(gid)%out_his_akt
      out_his_shflx = agrif_nml_storage(gid)%out_his_shflx
      out_his_shflx_rsw = agrif_nml_storage(gid)%out_his_shflx_rsw
# endif
# if defined SOLVE3D && defined SALINITY
      out_his_aks = agrif_nml_storage(gid)%out_his_aks
      out_his_swflx = agrif_nml_storage(gid)%out_his_swflx
# endif
# if defined SOLVE3D && defined BULK_FLUX
      out_his_shflx_rlw = agrif_nml_storage(gid)%out_his_shflx_rlw
      out_his_shflx_lat = agrif_nml_storage(gid)%out_his_shflx_lat
      out_his_shflx_sen = agrif_nml_storage(gid)%out_his_shflx_sen
# endif
# if defined SOLVE3D && defined BHFLUX
      out_his_bhflx = agrif_nml_storage(gid)%out_his_bhflx
# endif
# if defined SOLVE3D && defined BWFLUX && defined SALINITY
      out_his_bwflx = agrif_nml_storage(gid)%out_his_bwflx
# endif
# if defined SOLVE3D && (defined ANA_VMIX||defined LMD_MIXING|| defined LMD_SKPP||defined LMD_BKPP|| defined GLS_MIXING)
      out_his_bvf = agrif_nml_storage(gid)%out_his_bvf
# endif
# if defined SOLVE3D && (defined LMD_SKPP || defined GLS_MIXING)
      out_his_hbl = agrif_nml_storage(gid)%out_his_hbl
# endif
# if defined SOLVE3D && defined LMD_BKPP
      out_his_hbbl = agrif_nml_storage(gid)%out_his_hbbl
# endif
# if defined SOLVE3D && defined VIS_COEF_3D
      out_his_visc3d = agrif_nml_storage(gid)%out_his_visc3d
# endif
# if defined SOLVE3D && defined DIF_COEF_3D
      out_his_diff3d = agrif_nml_storage(gid)%out_his_diff3d
# endif
# if defined SOLVE3D && defined BIOLOGY && !defined PISCES
      out_his_hel = agrif_nml_storage(gid)%out_his_hel
# endif
# if defined SOLVE3D && defined BIO_NChlPZD
      out_his_chc = agrif_nml_storage(gid)%out_his_chc
      out_his_u10 = agrif_nml_storage(gid)%out_his_u10
      out_his_kvo2 = agrif_nml_storage(gid)%out_his_kvo2
      out_his_o2sat = agrif_nml_storage(gid)%out_his_o2sat
# endif
# if defined SOLVE3D && defined BIO_BioEBUS
      out_his_aou = agrif_nml_storage(gid)%out_his_aou
      out_his_wind10 = agrif_nml_storage(gid)%out_his_wind10
# endif
# if defined SOLVE3D && defined MORPHODYN
      out_his_hm = agrif_nml_storage(gid)%out_his_hm
# endif
# if defined AVERAGES && defined SOLVE3D
      out_avg_rho = agrif_nml_storage(gid)%out_avg_rho
      out_avg_omega = agrif_nml_storage(gid)%out_avg_omega
      out_avg_w = agrif_nml_storage(gid)%out_avg_w
      out_avg_akv = agrif_nml_storage(gid)%out_avg_akv
      out_avg_bostr = agrif_nml_storage(gid)%out_avg_bostr
      out_avg_bustr = agrif_nml_storage(gid)%out_avg_bustr
      out_avg_bvstr = agrif_nml_storage(gid)%out_avg_bvstr
      out_avg_wstr = agrif_nml_storage(gid)%out_avg_wstr
      out_avg_ustr = agrif_nml_storage(gid)%out_avg_ustr
      out_avg_vstr = agrif_nml_storage(gid)%out_avg_vstr
# endif
# if defined AVERAGES && defined SOLVE3D && defined TEMPERATURE
      out_avg_akt = agrif_nml_storage(gid)%out_avg_akt
      out_avg_shflx = agrif_nml_storage(gid)%out_avg_shflx
      out_avg_shflx_rsw = agrif_nml_storage(gid)%out_avg_shflx_rsw
# endif
# if defined AVERAGES && defined SOLVE3D && defined SALINITY
      out_avg_aks = agrif_nml_storage(gid)%out_avg_aks
      out_avg_swflx = agrif_nml_storage(gid)%out_avg_swflx
# endif
# if defined AVERAGES && defined SOLVE3D && defined BULK_FLUX
      out_avg_shflx_rlw = agrif_nml_storage(gid)%out_avg_shflx_rlw
      out_avg_shflx_lat = agrif_nml_storage(gid)%out_avg_shflx_lat
      out_avg_shflx_sen = agrif_nml_storage(gid)%out_avg_shflx_sen
# endif
# if defined AVERAGES && defined SOLVE3D && defined BHFLUX
      out_avg_bhflx = agrif_nml_storage(gid)%out_avg_bhflx
# endif
# if defined AVERAGES && defined SOLVE3D && defined BWFLUX && defined SALINITY
      out_avg_bwflx = agrif_nml_storage(gid)%out_avg_bwflx
# endif
# if defined AVERAGES && defined SOLVE3D && (defined ANA_VMIX||defined LMD_MIXING|| defined LMD_SKPP||defined LMD_BKPP|| defined GLS_MIXING)
      out_avg_bvf = agrif_nml_storage(gid)%out_avg_bvf
# endif
# if defined AVERAGES && defined SOLVE3D && (defined LMD_SKPP||defined GLS_MIXING)
      out_avg_hbl = agrif_nml_storage(gid)%out_avg_hbl
# endif
# if defined AVERAGES && defined SOLVE3D && defined LMD_BKPP
      out_avg_hbbl = agrif_nml_storage(gid)%out_avg_hbbl
# endif
# if defined AVERAGES && defined SOLVE3D && defined VIS_COEF_3D
      out_avg_visc3d = agrif_nml_storage(gid)%out_avg_visc3d
# endif
# if defined AVERAGES && defined SOLVE3D && defined DIF_COEF_3D
      out_avg_diff3d = agrif_nml_storage(gid)%out_avg_diff3d
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIOLOGY && !defined PISCES
      out_avg_hel = agrif_nml_storage(gid)%out_avg_hel
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_NChlPZD
      out_avg_chc = agrif_nml_storage(gid)%out_avg_chc
      out_avg_u10 = agrif_nml_storage(gid)%out_avg_u10
      out_avg_kvo2 = agrif_nml_storage(gid)%out_avg_kvo2
      out_avg_o2sat = agrif_nml_storage(gid)%out_avg_o2sat
# endif
# if defined AVERAGES && defined SOLVE3D && defined BIO_BioEBUS
      out_avg_aou = agrif_nml_storage(gid)%out_avg_aou
      out_avg_wind10 = agrif_nml_storage(gid)%out_avg_wind10
# endif
# if defined AVERAGES && defined SOLVE3D && defined MORPHODYN
      out_avg_hm = agrif_nml_storage(gid)%out_avg_hm
# endif
   end subroutine croco_namelist_restore_grid

   subroutine croco_namelist_switch_grid(gid)
      integer, intent(in) :: gid
      integer :: g
      if (gid <= 0) then
         g = 0
      else
         g = gid
      endif
      if (g > MAX_AGRIF_GRIDS) return
      if (g == agrif_nml_current_grid) return
      call croco_namelist_restore_grid(g)
      agrif_nml_current_grid = g
   end subroutine croco_namelist_switch_grid

