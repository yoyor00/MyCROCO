#!/usr/bin/env python3
"""
convert_in_to_nml.py
====================
Convert a CROCO `croco.in` + `cppdefs.h` + `param.h` set into
`croco.nml` + `cppdefs_config.h` + `param_config.h`.

Takes these 3 inputs and an optional one to specified CROCO source path
if the script is used elsewhere.

Usage
-----
    python convert_in_to_nml.py croco.in cppdefs.h param.h
    python convert_in_to_nml.py croco.in cppdefs.h param.h -o /path/to/rundir

Outputs (written to --outdir, default: current directory)
-----------------------------------------------------------
    croco.nml          Fortran namelist.  Array variables (tnu2, Akt_bak, …)
                        use compact Fortran repeat notation  NT*value  when
                        all values are equal, where NT is resolved from
                        param.h.

    param_config.h      {--src}/param.h (e.g. OCEAN/param.h) used as a
                        structural template - same line count, comments,
                        #ifdef sections, and #include lines - with each
                        simple literal value (grid dimensions, MPI, Tides,
                        point sources, tracer counts, ...) filled in from
                        the input param.h wherever cppdefs.h makes that
                        line active. Values absent from the input param.h
                        keep the template's own value; expression-valued
                        lines and #include lines are left untouched.

    cppdefs_config.h     The input cppdefs.h, copied through as-is.

Extending the mapping
---------------------
Add a row to MAPPINGS:

    # card_name       line  pos  type      nml_name                  nml_var
    ("time_stepping",   0,   0,  "int",   "&croco_time_stepping",   "ntimes"),

  card_name : key that starts the card in the .in file  (before the colon)
  line      : 0-based index of the value line inside that card block
  pos       : 0-based token index on that line (ignored for float_array)
  type      : "int" | "float" | "float_array" | "bool" | "bool_array[:SIZE]" | "str" | "str_line"
             | "labeled_bool"
             bool_array reads tokens[pos:] and trims/pads to SIZE elements, where
             SIZE is a param name (NT, NST, NumFluxTerms, …).  Plain "bool_array"
             defaults to NT.
             labeled_bool: pos is a label name (case-insensitive); the converter
             looks it up in the header comment of that card to find the column.
  nml_name  : namelist block name (with leading &)
  nml_var   : variable name inside that block

Multiple rows sharing the same nml_name are merged into one block.
Cards absent from the .in file are silently skipped (no block written).
"""

import re
import sys
import os
import shutil
import argparse
import subprocess

# ---------------------------------------------------------------------------
# Mapping table
# ---------------------------------------------------------------------------
MAPPINGS = [
    ("title", 0, 0, "str_line", "&croco_title", "title"),
    ("logfile", 0, 0, "str", "&croco_logfile", "logname"),
    ("time_stepping", 0, 0, "int", "&croco_time_stepping", "ntimes"),
    ("time_stepping", 0, 1, "float", "&croco_time_stepping", "dt"),
    ("time_stepping", 0, 2, "int", "&croco_time_stepping", "ndtfast"),
    ("time_stepping", 0, 3, "int", "&croco_time_stepping", "ninfo"),
    ("time_stepping_nbq", 0, 1, "float", "&croco_time_stepping_nbq", "csound_nbq"),
    ("time_stepping_nbq", 0, 2, "float", "&croco_time_stepping_nbq", "visc2_nbq"),
    ("S-coord", 0, 0, "float", "&croco_s_coord", "theta_s"),
    ("S-coord", 0, 1, "float", "&croco_s_coord", "theta_b"),
    ("S-coord", 0, 2, "float", "&croco_s_coord", "hc"),
    ("start_date", 0, 0, "str", "&croco_use_calendar", "start_date"),
    ("end_date", 0, 0, "str", "&croco_use_calendar", "end_date"),
    ("output_time_steps", 0, 0, "float", "&croco_use_calendar", "dt_his"),
    ("output_time_steps", 0, 1, "float", "&croco_use_calendar", "dt_avg"),
    ("output_time_steps", 0, 2, "float", "&croco_use_calendar", "dt_rst"),
    ("history", 0, 0, "bool", "&croco_history", "ldefhis"),
    ("history", 0, 1, "int", "&croco_history", "nwrt"),
    ("history", 0, 2, "int", "&croco_history", "nrpfhis"),
    ("history", 1, 0, "str", "&croco_history", "hisname"),
    ("initial", 0, 0, "int", "&croco_initial", "nrrec"),
    ("initial", 1, 0, "str", "&croco_initial", "ininame"),
    ("restart", 0, 0, "int", "&croco_restart", "nrst"),
    ("restart", 0, 1, "int", "&croco_restart", "nrpfrst"),
    ("restart", 1, 0, "str", "&croco_restart", "rstname"),
    ("grid", 0, 0, "str", "&croco_grid", "grdname"),
    ("forcing", 0, 0, "str", "&croco_forcing", "frcname"),
    ("sponge", 0, 0, "float", "&croco_sponge", "x_sponge"),
    ("nudg_cof", 0, 0, "float", "&croco_nudging", "tauT_in"),
    ("nudg_cof", 0, 1, "float", "&croco_nudging", "tauT_out"),
    ("nudg_cof", 0, 2, "float", "&croco_nudging", "tauM_in"),
    ("nudg_cof", 0, 3, "float", "&croco_nudging", "tauM_out"),
    ("bottom_forcing", 0, 0, "str", "&croco_bottom_forcing", "btfname"),
    ("bulk_forcing", 0, 0, "str", "&croco_bulk_forcing", "bulkname"),
    ("climatology", 0, 0, "str", "&croco_climatology", "clmname"),
    ("wave_offline", 0, 0, "str", "&croco_wave_offline", "wave_file"),
    ("biology", 0, 0, "str", "&croco_biology", "bioname"),
    ("boundary", 0, 0, "str", "&croco_boundary", "bry_file"),
    ("wkb_boundary", 0, 0, "str", "&croco_wkb_boundary", "brywkb_file"),
    ("bodyforce", 0, 0, "int", "&croco_bodyforce", "levsfrc"),
    ("bodyforce", 0, 1, "int", "&croco_bodyforce", "levbfrc"),
    ("lin_EOS_cff", 0, 0, "float", "&croco_lin_eos", "R0"),
    ("lin_EOS_cff", 0, 1, "float", "&croco_lin_eos", "T0"),
    ("lin_EOS_cff", 0, 2, "float", "&croco_lin_eos", "S0"),
    ("lin_EOS_cff", 0, 3, "float", "&croco_lin_eos", "Tcoef"),
    ("lin_EOS_cff", 0, 4, "float", "&croco_lin_eos", "Scoef"),
    ("abl_nudg_tra_cof", 0, 0, "float", "&croco_abl_nudg_tra", "ltra_min"),
    ("abl_nudg_tra_cof", 0, 1, "float", "&croco_abl_nudg_tra", "ltra_max"),
    ("abl_nudg_dyn_cof", 0, 0, "float", "&croco_abl_nudg_dyn", "ldyn_min"),
    ("abl_nudg_dyn_cof", 0, 1, "float", "&croco_abl_nudg_dyn", "ldyn_max"),
    ("sediments", 0, 0, "str", "&croco_sediments", "sedname"),
    ("sediments_mustang", 0, 0, "str", "&croco_sediments_mustang", "sedname_must"),
    ("substance", 0, 0, "str", "&croco_substance", "subsfilename"),
    ("obstruction", 0, 0, "str", "&croco_obstruction", "obstname"),
    ("xios_origin_date", 0, 0, "str_line", "&croco_xios_origin_date", "xios_origin_date"),
    ("assimilation", 0, 0, "str", "&croco_assimilation", "aparnam"),
    ("assimilation", 1, 0, "str", "&croco_assimilation", "assname"),
    ("rho0", 0, 0, "float", "&croco_rho0", "rho0"),
    ("bottom_drag", 0, 0, "float", "&croco_bottom_drag", "rdrg"),
    ("bottom_drag", 0, 1, "float", "&croco_bottom_drag", "rdrg2"),
    ("bottom_drag", 0, 2, "float", "&croco_bottom_drag", "Zobt"),
    ("bottom_drag", 0, 3, "float", "&croco_bottom_drag", "Cdb_min"),
    ("bottom_drag", 0, 4, "float", "&croco_bottom_drag", "Cdb_max"),
    ("gamma2", 0, 0, "float", "&croco_gamma2", "gamma2"),
    ("lateral_visc", 0, 0, "float", "&croco_lateral_visc", "visc2"),
    ("lateral_visc", 0, 1, "float", "&croco_lateral_visc", "visc4"),
    ("tracer_diff2", 0, 0, "float_array", "&croco_tracer_diff2", "tnu2"),
    ("tracer_diff4", 0, 0, "float_array", "&croco_tracer_diff4", "tnu4"),
    ("vertical_mixing", 0, 0, "float", "&croco_vertical_mixing", "Akv_bak"),
    ("vertical_mixing", 0, 1, "float_array", "&croco_vertical_mixing", "Akt_bak"),
    ("wkb_wwave", 0, 0, "float", "&croco_wkb_wwave", "wkb_amp"),
    ("wkb_wwave", 0, 1, "float", "&croco_wkb_wwave", "wkb_ang"),
    ("wkb_wwave", 0, 2, "float", "&croco_wkb_wwave", "wkb_prd"),
    ("wkb_wwave", 0, 3, "float", "&croco_wkb_wwave", "wkb_tide"),
    ("wkb_wwave", 0, 4, "float", "&croco_wkb_wwave", "wkb_btg"),
    ("wkb_wwave", 0, 5, "float", "&croco_wkb_wwave", "wkb_gam"),
    ("wkb_roller", 0, 0, "float", "&croco_wkb_roller", "wkb_rsb"),
    ("wkb_roller", 0, 1, "float", "&croco_wkb_roller", "wkb_roller"),
    ("wave_maker", 0, 0, "float", "&croco_wave_maker", "wmaker_amp"),
    ("wave_maker", 0, 1, "float", "&croco_wave_maker", "wmaker_prd"),
    ("wave_maker", 0, 2, "float", "&croco_wave_maker", "wmaker_dir"),
    ("wave_maker", 0, 3, "float", "&croco_wave_maker", "wmaker_dsp"),
    ("wave_maker", 0, 4, "float", "&croco_wave_maker", "wmaker_fsp"),
    ("averages", 0, 0, "int", "&croco_averages", "ntsavg"),
    ("averages", 0, 1, "int", "&croco_averages", "navg"),
    ("averages", 0, 2, "int", "&croco_averages", "nrpfavg"),
    ("averages", 1, 0, "str", "&croco_averages", "avgname"),
    ("surf", 0, 0, "bool", "&croco_surf", "ldefsurf"),
    ("surf", 0, 1, "int", "&croco_surf", "nwrtsurf"),
    ("surf", 0, 2, "int", "&croco_surf", "nrpfsurf"),
    ("surf", 1, 0, "str", "&croco_surf", "surfname"),
    ("surf_avg", 0, 0, "bool", "&croco_surf_avg", "ldefsurf_avg"),
    ("surf_avg", 0, 1, "int", "&croco_surf_avg", "ntssurf_avg"),
    ("surf_avg", 0, 2, "int", "&croco_surf_avg", "nwrtsurf_avg"),
    ("surf_avg", 0, 3, "int", "&croco_surf_avg", "nrpfsurf_avg"),
    ("surf_avg", 1, 0, "str", "&croco_surf_avg", "surfname_avg"),
    ("stations", 0, 0, "bool", "&croco_stations", "ldefsta"),
    ("stations", 0, 1, "int", "&croco_stations", "nsta"),
    ("stations", 0, 2, "int", "&croco_stations", "nrpfsta"),
    ("stations", 1, 0, "str", "&croco_stations", "staposname"),
    ("stations", 2, 0, "str", "&croco_stations", "staname"),
    ("online", 0, 0, "int", "&croco_online", "yearnum"),
    ("online", 0, 1, "int", "&croco_online", "monthnum"),
    ("online", 0, 2, "int", "&croco_online", "recordsperday"),
    ("online", 0, 3, "int", "&croco_online", "yearend"),
    ("online", 0, 4, "int", "&croco_online", "monthend"),
    ("online", 1, 0, "str", "&croco_online", "pathbulk"),
    ("diagnostics", 0, 0, "bool", "&croco_diagnostics_ts", "ldefdia"),
    ("diagnostics", 0, 1, "int", "&croco_diagnostics_ts", "nwrtdia"),
    ("diagnostics", 0, 2, "int", "&croco_diagnostics_ts", "nrpfdia"),
    ("diagnostics", 1, 0, "str", "&croco_diagnostics_ts", "dianame"),
    ("diag_avg", 0, 0, "bool", "&croco_diag_avg", "ldefdia_avg"),
    ("diag_avg", 0, 1, "int", "&croco_diag_avg", "ntsdia_avg"),
    ("diag_avg", 0, 2, "int", "&croco_diag_avg", "nwrtdia_avg"),
    ("diag_avg", 0, 3, "int", "&croco_diag_avg", "nrpfdia_avg"),
    ("diag_avg", 1, 0, "str", "&croco_diag_avg", "dianame_avg"),
    ("diag_mld_dens", 0, 0, "float", "&croco_diag_mld_dens", "mld_crit_D"),
    ("diag_mld_dens", 0, 1, "float", "&croco_diag_mld_dens", "mld_crit_T"),
    ("diagnosticsM", 0, 0, "bool", "&croco_diagnosticsM", "ldefdiaM"),
    ("diagnosticsM", 0, 1, "int", "&croco_diagnosticsM", "nwrtdiaM"),
    ("diagnosticsM", 0, 2, "int", "&croco_diagnosticsM", "nrpfdiaM"),
    ("diagnosticsM", 1, 0, "str", "&croco_diagnosticsM", "dianameM"),
    ("diagM_avg", 0, 0, "bool", "&croco_diagM_avg", "ldefdiaM_avg"),
    ("diagM_avg", 0, 1, "int", "&croco_diagM_avg", "ntsdiaM_avg"),
    ("diagM_avg", 0, 2, "int", "&croco_diagM_avg", "nwrtdiaM_avg"),
    ("diagM_avg", 0, 3, "int", "&croco_diagM_avg", "nrpfdiaM_avg"),
    ("diagM_avg", 1, 0, "str", "&croco_diagM_avg", "dianameM_avg"),
    ("diags_vrt", 0, 0, "bool", "&croco_diags_vrt", "ldefdiags_vrt"),
    ("diags_vrt", 0, 1, "int", "&croco_diags_vrt", "nwrtdiags_vrt"),
    ("diags_vrt", 0, 2, "int", "&croco_diags_vrt", "nrpfdiags_vrt"),
    ("diags_vrt", 1, 0, "str", "&croco_diags_vrt", "diags_vrtname"),
    ("diags_vrt_avg", 0, 0, "bool", "&croco_diags_vrt_avg", "ldefdiags_vrt_avg"),
    ("diags_vrt_avg", 0, 1, "int", "&croco_diags_vrt_avg", "ntsdiags_vrt_avg"),
    ("diags_vrt_avg", 0, 2, "int", "&croco_diags_vrt_avg", "nwrtdiags_vrt_avg"),
    ("diags_vrt_avg", 0, 3, "int", "&croco_diags_vrt_avg", "nrpfdiags_vrt_avg"),
    ("diags_vrt_avg", 1, 0, "str", "&croco_diags_vrt_avg", "diags_vrtname_avg"),
    ("diags_ek", 0, 0, "bool", "&croco_diags_ek", "ldefdiags_ek"),
    ("diags_ek", 0, 1, "int", "&croco_diags_ek", "nwrtdiags_ek"),
    ("diags_ek", 0, 2, "int", "&croco_diags_ek", "nrpfdiags_ek"),
    ("diags_ek", 1, 0, "str", "&croco_diags_ek", "diags_ekname"),
    ("diags_ek_avg", 0, 0, "bool", "&croco_diags_ek_avg", "ldefdiags_ek_avg"),
    ("diags_ek_avg", 0, 1, "int", "&croco_diags_ek_avg", "ntsdiags_ek_avg"),
    ("diags_ek_avg", 0, 2, "int", "&croco_diags_ek_avg", "nwrtdiags_ek_avg"),
    ("diags_ek_avg", 0, 3, "int", "&croco_diags_ek_avg", "nrpfdiags_ek_avg"),
    ("diags_ek_avg", 1, 0, "str", "&croco_diags_ek_avg", "diags_ekname_avg"),
    ("diags_pv", 0, 0, "bool", "&croco_diags_pv", "ldefdiags_pv"),
    ("diags_pv", 0, 1, "int", "&croco_diags_pv", "nwrtdiags_pv"),
    ("diags_pv", 0, 2, "int", "&croco_diags_pv", "nrpfdiags_pv"),
    ("diags_pv", 1, 0, "str", "&croco_diags_pv", "diags_pvname"),
    ("diags_pv_avg", 0, 0, "bool", "&croco_diags_pv_avg", "ldefdiags_pv_avg"),
    ("diags_pv_avg", 0, 1, "int", "&croco_diags_pv_avg", "ntsdiags_pv_avg"),
    ("diags_pv_avg", 0, 2, "int", "&croco_diags_pv_avg", "nwrtdiags_pv_avg"),
    ("diags_pv_avg", 0, 3, "int", "&croco_diags_pv_avg", "nrpfdiags_pv_avg"),
    ("diags_pv_avg", 1, 0, "str", "&croco_diags_pv_avg", "diags_pvname_avg"),
    ("diags_eddy_avg", 0, 0, "bool", "&croco_diags_eddy_avg", "ldefdiags_eddy_avg"),
    ("diags_eddy_avg", 0, 1, "int", "&croco_diags_eddy_avg", "ntsdiags_eddy_avg"),
    ("diags_eddy_avg", 0, 2, "int", "&croco_diags_eddy_avg", "nwrtdiags_eddy_avg"),
    ("diags_eddy_avg", 0, 3, "int", "&croco_diags_eddy_avg", "nrpfdiags_eddy_avg"),
    ("diags_eddy_avg", 1, 0, "str", "&croco_diags_eddy_avg", "diags_eddyname_avg"),
    ("diagnostics_bio", 0, 0, "bool", "&croco_diagnostics_bio", "ldefdiabio"),
    ("diagnostics_bio", 0, 1, "int", "&croco_diagnostics_bio", "nwrtdiabio"),
    ("diagnostics_bio", 0, 2, "int", "&croco_diagnostics_bio", "nrpfdiabio"),
    ("diagnostics_bio", 1, 0, "str", "&croco_diagnostics_bio", "dianamebio"),
    ("diagbio_avg", 0, 0, "bool", "&croco_diagbio_avg", "ldefdiabio_avg"),
    ("diagbio_avg", 0, 1, "int", "&croco_diagbio_avg", "ntsdiabio_avg"),
    ("diagbio_avg", 0, 2, "int", "&croco_diagbio_avg", "nwrtdiabio_avg"),
    ("diagbio_avg", 0, 3, "int", "&croco_diagbio_avg", "nrpfdiabio_avg"),
    ("diagbio_avg", 1, 0, "str", "&croco_diagbio_avg", "dianamebio_avg"),
    ("abl", 0, 0, "bool", "&croco_abl", "ldefablhis"),
    ("abl", 0, 1, "int", "&croco_abl", "nwrtablhis"),
    ("abl", 0, 2, "int", "&croco_abl", "nrpfablhis"),
    ("abl", 1, 0, "str", "&croco_abl", "ablname"),
    ("abl_averages", 0, 0, "bool", "&croco_abl_averages", "ldefablavg"),
    ("abl_averages", 0, 1, "int", "&croco_abl_averages", "ntsablavg"),
    ("abl_averages", 0, 2, "int", "&croco_abl_averages", "nwrtablavg"),
    ("abl_averages", 0, 3, "int", "&croco_abl_averages", "nrpfablavg"),
    ("abl_averages", 1, 0, "str", "&croco_abl_averages", "ablname_avg"),
    # stochastic_history_fields: xi2d xi3d
    ("stochastic_history_fields", 0, 0, "bool", "&croco_stochastic_history_fields", "out_his_xi2d"),
    ("stochastic_history_fields", 0, 1, "bool", "&croco_stochastic_history_fields", "out_his_xi3d"),
    # gls_history_fields: tke gls lscale
    ("gls_history_fields", 0, 0, "bool", "&croco_gls_history_fields", "out_his_tke"),
    ("gls_history_fields", 0, 1, "bool", "&croco_gls_history_fields", "out_his_gls"),
    ("gls_history_fields", 0, 2, "bool", "&croco_gls_history_fields", "out_his_lscale"),
    # gls_averages: tke gls lscale  (renamed to croco_gls_averages_fields)
    ("gls_averages", 0, 0, "bool", "&croco_gls_averages_fields", "out_avg_tke"),
    ("gls_averages", 0, 1, "bool", "&croco_gls_averages_fields", "out_avg_gls"),
    ("gls_averages", 0, 2, "bool", "&croco_gls_averages_fields", "out_avg_lscale"),
    # abl_history_fields: 20 bools
    ("abl_history_fields", 0, 0, "bool", "&croco_abl_history_fields", "out_his_abl_pu_dta"),
    ("abl_history_fields", 0, 1, "bool", "&croco_abl_history_fields", "out_his_abl_pv_dta"),
    ("abl_history_fields", 0, 2, "bool", "&croco_abl_history_fields", "out_his_abl_pt_dta"),
    ("abl_history_fields", 0, 3, "bool", "&croco_abl_history_fields", "out_his_abl_pq_dta"),
    ("abl_history_fields", 0, 4, "bool", "&croco_abl_history_fields", "out_his_abl_pgu_dta"),
    ("abl_history_fields", 0, 5, "bool", "&croco_abl_history_fields", "out_his_abl_pgv_dta"),
    ("abl_history_fields", 0, 6, "bool", "&croco_abl_history_fields", "out_his_abl_u_abl"),
    ("abl_history_fields", 0, 7, "bool", "&croco_abl_history_fields", "out_his_abl_v_abl"),
    ("abl_history_fields", 0, 8, "bool", "&croco_abl_history_fields", "out_his_abl_t_abl"),
    ("abl_history_fields", 0, 9, "bool", "&croco_abl_history_fields", "out_his_abl_q_abl"),
    ("abl_history_fields", 0, 10, "bool", "&croco_abl_history_fields", "out_his_abl_tke_abl"),
    ("abl_history_fields", 0, 11, "bool", "&croco_abl_history_fields", "out_his_abl_mxlm_abl"),
    ("abl_history_fields", 0, 12, "bool", "&croco_abl_history_fields", "out_his_abl_mxld_abl"),
    ("abl_history_fields", 0, 13, "bool", "&croco_abl_history_fields", "out_his_abl_avm_abl"),
    ("abl_history_fields", 0, 14, "bool", "&croco_abl_history_fields", "out_his_abl_avt_abl"),
    ("abl_history_fields", 0, 15, "bool", "&croco_abl_history_fields", "out_his_abl_ablh_abl"),
    ("abl_history_fields", 0, 16, "bool", "&croco_abl_history_fields", "out_his_abl_zr_abl"),
    ("abl_history_fields", 0, 17, "bool", "&croco_abl_history_fields", "out_his_abl_zw_abl"),
    ("abl_history_fields", 0, 18, "bool", "&croco_abl_history_fields", "out_his_abl_Hzr_abl"),
    ("abl_history_fields", 0, 19, "bool", "&croco_abl_history_fields", "out_his_abl_Hzw_abl"),
    # abl_averages_fields: 20 bools
    ("abl_averages_fields", 0, 0, "bool", "&croco_abl_averages_fields", "out_avg_abl_pu_dta"),
    ("abl_averages_fields", 0, 1, "bool", "&croco_abl_averages_fields", "out_avg_abl_pv_dta"),
    ("abl_averages_fields", 0, 2, "bool", "&croco_abl_averages_fields", "out_avg_abl_pt_dta"),
    ("abl_averages_fields", 0, 3, "bool", "&croco_abl_averages_fields", "out_avg_abl_pq_dta"),
    ("abl_averages_fields", 0, 4, "bool", "&croco_abl_averages_fields", "out_avg_abl_pgu_dta"),
    ("abl_averages_fields", 0, 5, "bool", "&croco_abl_averages_fields", "out_avg_abl_pgv_dta"),
    ("abl_averages_fields", 0, 6, "bool", "&croco_abl_averages_fields", "out_avg_abl_u_abl"),
    ("abl_averages_fields", 0, 7, "bool", "&croco_abl_averages_fields", "out_avg_abl_v_abl"),
    ("abl_averages_fields", 0, 8, "bool", "&croco_abl_averages_fields", "out_avg_abl_t_abl"),
    ("abl_averages_fields", 0, 9, "bool", "&croco_abl_averages_fields", "out_avg_abl_q_abl"),
    ("abl_averages_fields", 0, 10, "bool", "&croco_abl_averages_fields", "out_avg_abl_tke_abl"),
    ("abl_averages_fields", 0, 11, "bool", "&croco_abl_averages_fields", "out_avg_abl_mxlm_abl"),
    ("abl_averages_fields", 0, 12, "bool", "&croco_abl_averages_fields", "out_avg_abl_mxld_abl"),
    ("abl_averages_fields", 0, 13, "bool", "&croco_abl_averages_fields", "out_avg_abl_avm_abl"),
    ("abl_averages_fields", 0, 14, "bool", "&croco_abl_averages_fields", "out_avg_abl_avt_abl"),
    ("abl_averages_fields", 0, 15, "bool", "&croco_abl_averages_fields", "out_avg_abl_ablh_abl"),
    ("abl_averages_fields", 0, 16, "bool", "&croco_abl_averages_fields", "out_avg_abl_zr_abl"),
    ("abl_averages_fields", 0, 17, "bool", "&croco_abl_averages_fields", "out_avg_abl_zw_abl"),
    ("abl_averages_fields", 0, 18, "bool", "&croco_abl_averages_fields", "out_avg_abl_Hzr_abl"),
    ("abl_averages_fields", 0, 19, "bool", "&croco_abl_averages_fields", "out_avg_abl_Hzw_abl"),
    # surf_history_fields: 1 bool
    ("surf_history_fields", 0, 0, "bool", "&croco_surf_history_fields", "out_his_surf"),
    # surf_average_fields: 1 bool
    ("surf_average_fields", 0, 0, "bool", "&croco_surf_average_fields", "out_avg_surf"),
    # station_fields: grd temp salt rho vel
    ("station_fields", 0, 0, "bool", "&croco_station_fields", "sta_grd"),
    ("station_fields", 0, 1, "bool", "&croco_station_fields", "sta_temp"),
    ("station_fields", 0, 2, "bool", "&croco_station_fields", "sta_salt"),
    ("station_fields", 0, 3, "bool", "&croco_station_fields", "sta_rho"),
    ("station_fields", 0, 4, "bool", "&croco_station_fields", "sta_vel"),
    # sediment_history_fields: athk bthk bpor bfra(1:NST) [dflx eflx](NST) [bdlu bdlv](NST) [btcr]
    # CPP-conditional fields (dflx/eflx/bdlu/bdlv/btcr) are handled by
    # _build_sediment_extra_history_blocks() because their offset depends on NST and CPP keys.
    ("sediment_history_fields", 0, 0, "bool", "&croco_sediment_history_fields", "out_his_sed_athk"),
    ("sediment_history_fields", 0, 1, "bool", "&croco_sediment_history_fields", "out_his_sed_bthk"),
    ("sediment_history_fields", 0, 2, "bool", "&croco_sediment_history_fields", "out_his_sed_bpor"),
    ("sediment_history_fields", 0, 3, "bool_array:NST", "&croco_sediment_bfra_history_fields", "out_his_sed_bfra"),
    # bbl_history_fields: abed hripple lripple zbnot zbapp bostrw
    ("bbl_history_fields", 0, 0, "bool", "&croco_bbl_history_fields", "out_his_abed"),
    ("bbl_history_fields", 0, 1, "bool", "&croco_bbl_history_fields", "out_his_hripple"),
    ("bbl_history_fields", 0, 2, "bool", "&croco_bbl_history_fields", "out_his_lripple"),
    ("bbl_history_fields", 0, 3, "bool", "&croco_bbl_history_fields", "out_his_zbnot"),
    ("bbl_history_fields", 0, 4, "bool", "&croco_bbl_history_fields", "out_his_zbapp"),
    ("bbl_history_fields", 0, 5, "bool", "&croco_bbl_history_fields", "out_his_bostrw"),
    # wci_history_fields: sup ust2d vst2d [ust vst wst akb akw kvf calp kaps]
    ("wci_history_fields", 0, 0, "bool", "&croco_wci_history_fields", "out_his_sup"),
    ("wci_history_fields", 0, 1, "bool", "&croco_wci_history_fields", "out_his_ust2d"),
    ("wci_history_fields", 0, 2, "bool", "&croco_wci_history_fields", "out_his_vst2d"),
    ("wci_history_fields", 0, 3, "bool", "&croco_wci_history_3d_fields", "out_his_ust"),
    ("wci_history_fields", 0, 4, "bool", "&croco_wci_history_3d_fields", "out_his_vst"),
    ("wci_history_fields", 0, 5, "bool", "&croco_wci_history_3d_fields", "out_his_wst"),
    ("wci_history_fields", 0, 6, "bool", "&croco_wci_history_3d_fields", "out_his_akb"),
    ("wci_history_fields", 0, 7, "bool", "&croco_wci_history_3d_fields", "out_his_akw"),
    ("wci_history_fields", 0, 8, "bool", "&croco_wci_history_3d_fields", "out_his_kvf"),
    ("wci_history_fields", 0, 9, "bool", "&croco_wci_history_3d_fields", "out_his_calp"),
    ("wci_history_fields", 0, 10, "bool", "&croco_wci_history_3d_fields", "out_his_kaps"),
    # wci_average_fields: sup ust2d vst2d [ust vst wst akb akw kvf calp kaps]
    ("wci_average_fields", 0, 0, "bool", "&croco_wci_average_fields", "out_avg_sup"),
    ("wci_average_fields", 0, 1, "bool", "&croco_wci_average_fields", "out_avg_ust2d"),
    ("wci_average_fields", 0, 2, "bool", "&croco_wci_average_fields", "out_avg_vst2d"),
    ("wci_average_fields", 0, 3, "bool", "&croco_wci_average_3d_fields", "out_avg_ust"),
    ("wci_average_fields", 0, 4, "bool", "&croco_wci_average_3d_fields", "out_avg_vst"),
    ("wci_average_fields", 0, 5, "bool", "&croco_wci_average_3d_fields", "out_avg_wst"),
    ("wci_average_fields", 0, 6, "bool", "&croco_wci_average_3d_fields", "out_avg_akb"),
    ("wci_average_fields", 0, 7, "bool", "&croco_wci_average_3d_fields", "out_avg_akw"),
    ("wci_average_fields", 0, 8, "bool", "&croco_wci_average_3d_fields", "out_avg_kvf"),
    ("wci_average_fields", 0, 9, "bool", "&croco_wci_average_3d_fields", "out_avg_calp"),
    ("wci_average_fields", 0, 10, "bool", "&croco_wci_average_3d_fields", "out_avg_kaps"),
    # wave_history_fields: hrm frq action k_xi k_eta eps_b eps_d erol eps_r
    ("wave_history_fields", 0, 0, "bool", "&croco_wave_history_fields", "out_his_hrm"),
    ("wave_history_fields", 0, 1, "bool", "&croco_wave_history_fields", "out_his_frq"),
    ("wave_history_fields", 0, 2, "bool", "&croco_wave_history_fields", "out_his_action"),
    ("wave_history_fields", 0, 3, "bool", "&croco_wave_history_fields", "out_his_k_xi"),
    ("wave_history_fields", 0, 4, "bool", "&croco_wave_history_fields", "out_his_k_eta"),
    ("wave_history_fields", 0, 5, "bool", "&croco_wave_history_fields", "out_his_eps_b"),
    ("wave_history_fields", 0, 6, "bool", "&croco_wave_history_fields", "out_his_eps_d"),
    ("wave_history_fields", 0, 7, "bool", "&croco_wave_history_fields", "out_his_erol"),
    ("wave_history_fields", 0, 8, "bool", "&croco_wave_history_fields", "out_his_eps_r"),
    # wave_average_fields: hrm frq action k_xi k_eta eps_b eps_d erol eps_r
    ("wave_average_fields", 0, 0, "bool", "&croco_wave_average_fields", "out_avg_hrm"),
    ("wave_average_fields", 0, 1, "bool", "&croco_wave_average_fields", "out_avg_frq"),
    ("wave_average_fields", 0, 2, "bool", "&croco_wave_average_fields", "out_avg_action"),
    ("wave_average_fields", 0, 3, "bool", "&croco_wave_average_fields", "out_avg_k_xi"),
    ("wave_average_fields", 0, 4, "bool", "&croco_wave_average_fields", "out_avg_k_eta"),
    ("wave_average_fields", 0, 5, "bool", "&croco_wave_average_fields", "out_avg_eps_b"),
    ("wave_average_fields", 0, 6, "bool", "&croco_wave_average_fields", "out_avg_eps_d"),
    ("wave_average_fields", 0, 7, "bool", "&croco_wave_average_fields", "out_avg_erol"),
    ("wave_average_fields", 0, 8, "bool", "&croco_wave_average_fields", "out_avg_eps_r"),
    # primary_history_fields: zeta ubar vbar [u v [tracer(1:NT)]]
    ("primary_history_fields", 0, 0, "bool", "&croco_primary_history_fields", "out_his_zeta"),
    ("primary_history_fields", 0, 1, "bool", "&croco_primary_history_fields", "out_his_ubar"),
    ("primary_history_fields", 0, 2, "bool", "&croco_primary_history_fields", "out_his_vbar"),
    ("primary_history_fields", 0, 3, "bool", "&croco_primary_history_3d_fields", "out_his_u"),
    ("primary_history_fields", 0, 4, "bool", "&croco_primary_history_3d_fields", "out_his_v"),
    ("primary_history_fields", 0, 5, "bool_array", "&croco_primary_history_tracer_fields", "out_his_tracer"),
    # primary_averages: zeta ubar vbar [u v [tracer(1:NT)]]
    ("primary_averages", 0, 0, "bool", "&croco_primary_average_fields", "out_avg_zeta"),
    ("primary_averages", 0, 1, "bool", "&croco_primary_average_fields", "out_avg_ubar"),
    ("primary_averages", 0, 2, "bool", "&croco_primary_average_fields", "out_avg_vbar"),
    ("primary_averages", 0, 3, "bool", "&croco_primary_3d_average_fields", "out_avg_u"),
    ("primary_averages", 0, 4, "bool", "&croco_primary_3d_average_fields", "out_avg_v"),
    ("primary_averages", 0, 5, "bool_array", "&croco_primary_tracer_average_fields", "out_avg_tracer"),
    # auxiliary_history_fields: label-based (order varies per case)
    ("auxiliary_history_fields", 0, "rho", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_rho"),
    ("auxiliary_history_fields", 0, "Omega", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_omega"),
    ("auxiliary_history_fields", 0, "W", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_w"),
    ("auxiliary_history_fields", 0, "Akv", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_akv"),
    ("auxiliary_history_fields", 0, "Bostr", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_bostr"),
    ("auxiliary_history_fields", 0, "Bustr", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_bustr"),
    ("auxiliary_history_fields", 0, "Bvstr", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_bvstr"),
    ("auxiliary_history_fields", 0, "Wstr", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_wstr"),
    ("auxiliary_history_fields", 0, "Ustr", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_ustr"),
    ("auxiliary_history_fields", 0, "Vstr", "labeled_bool", "&croco_auxiliary_history_fields", "out_his_vstr"),
    ("auxiliary_history_fields", 0, "Akt", "labeled_bool", "&croco_temperature_history_fields", "out_his_akt"),
    ("auxiliary_history_fields", 0, "Shfl", "labeled_bool", "&croco_temperature_history_fields", "out_his_shflx"),
    ("auxiliary_history_fields", 0, "rsw", "labeled_bool", "&croco_temperature_history_fields", "out_his_shflx_rsw"),
    ("auxiliary_history_fields", 0, "Aks", "labeled_bool", "&croco_salinity_history_fields", "out_his_aks"),
    ("auxiliary_history_fields", 0, "Swfl", "labeled_bool", "&croco_salinity_history_fields", "out_his_swflx"),
    ("auxiliary_history_fields", 0, "rlw", "labeled_bool", "&croco_bulk_flux_history_fields", "out_his_shflx_rlw"),
    ("auxiliary_history_fields", 0, "lat", "labeled_bool", "&croco_bulk_flux_history_fields", "out_his_shflx_lat"),
    ("auxiliary_history_fields", 0, "sen", "labeled_bool", "&croco_bulk_flux_history_fields", "out_his_shflx_sen"),
    ("auxiliary_history_fields", 0, "Bvf", "labeled_bool", "&croco_bvf_history_fields", "out_his_bvf"),
    ("auxiliary_history_fields", 0, "HBL", "labeled_bool", "&croco_hbl_history_fields", "out_his_hbl"),
    ("auxiliary_history_fields", 0, "HBBL", "labeled_bool", "&croco_lmd_bkpp_history_fields", "out_his_hbbl"),
    ("auxiliary_history_fields", 0, "Visc3d", "labeled_bool", "&croco_vis_coef_history_fields", "out_his_visc3d"),
    ("auxiliary_history_fields", 0, "Diff3d", "labeled_bool", "&croco_dif_coef_history_fields", "out_his_diff3d"),
    ("auxiliary_history_fields", 0, "HEL", "labeled_bool", "&croco_biology_history_fields", "out_his_hel"),
    ("auxiliary_history_fields", 0, "Hm", "labeled_bool", "&croco_morphodyn_history_fields", "out_his_hm", "_morphodyn", ".true."),
    # auxiliary_averages: label-based (order varies per case)
    ("auxiliary_averages", 0, "rho", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_rho"),
    ("auxiliary_averages", 0, "Omega", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_omega"),
    ("auxiliary_averages", 0, "W", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_w"),
    ("auxiliary_averages", 0, "Akv", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_akv"),
    ("auxiliary_averages", 0, "Bostr", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_bostr"),
    ("auxiliary_averages", 0, "Bustr", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_bustr"),
    ("auxiliary_averages", 0, "Bvstr", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_bvstr"),
    ("auxiliary_averages", 0, "Wstr", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_wstr"),
    ("auxiliary_averages", 0, "Ustr", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_ustr"),
    ("auxiliary_averages", 0, "Vstr", "labeled_bool", "&croco_auxiliary_averages_fields", "out_avg_vstr"),
    ("auxiliary_averages", 0, "Akt", "labeled_bool", "&croco_temperature_averages_fields", "out_avg_akt"),
    ("auxiliary_averages", 0, "Shfl", "labeled_bool", "&croco_temperature_averages_fields", "out_avg_shflx"),
    ("auxiliary_averages", 0, "rsw", "labeled_bool", "&croco_temperature_averages_fields", "out_avg_shflx_rsw"),
    ("auxiliary_averages", 0, "Aks", "labeled_bool", "&croco_salinity_averages_fields", "out_avg_aks"),
    ("auxiliary_averages", 0, "Swfl", "labeled_bool", "&croco_salinity_averages_fields", "out_avg_swflx"),
    ("auxiliary_averages", 0, "rlw", "labeled_bool", "&croco_bulk_flux_averages_fields", "out_avg_shflx_rlw"),
    ("auxiliary_averages", 0, "lat", "labeled_bool", "&croco_bulk_flux_averages_fields", "out_avg_shflx_lat"),
    ("auxiliary_averages", 0, "sen", "labeled_bool", "&croco_bulk_flux_averages_fields", "out_avg_shflx_sen"),
    ("auxiliary_averages", 0, "Bvf", "labeled_bool", "&croco_bvf_averages_fields", "out_avg_bvf"),
    ("auxiliary_averages", 0, "HBL", "labeled_bool", "&croco_hbl_averages_fields", "out_avg_hbl"),
    ("auxiliary_averages", 0, "HBBL", "labeled_bool", "&croco_lmd_bkpp_averages_fields", "out_avg_hbbl"),
    ("auxiliary_averages", 0, "Visc3d", "labeled_bool", "&croco_vis_coef_averages_fields", "out_avg_visc3d"),
    ("auxiliary_averages", 0, "Diff3d", "labeled_bool", "&croco_dif_coef_averages_fields", "out_avg_diff3d"),
    ("auxiliary_averages", 0, "HEL", "labeled_bool", "&croco_biology_averages_fields", "out_avg_hel"),
    ("auxiliary_averages", 0, "Hm", "labeled_bool", "&croco_morphodyn_averages_fields", "out_avg_hm", "_morphodyn", ".false."),
    ("diag3D_history_fields", 0, 0, "bool_array", "&croco_diag3D_history_fields", "out_his_dia3D_tracer"),
    ("diag2D_history_fields", 0, 0, "bool_array", "&croco_diag2D_history_fields", "out_his_dia2D_tracer"),
    ("diag3D_average_fields", 0, 0, "bool_array", "&croco_diag3D_average_fields", "out_avg_dia3D_tracer"),
    ("diag2D_average_fields", 0, 0, "bool_array", "&croco_diag2D_average_fields", "out_avg_dia2D_tracer"),
    ("diagM_history_fields", 0, 0, "bool", "&croco_diagM_history_fields", "out_his_diagM_u"),
    ("diagM_history_fields", 0, 1, "bool", "&croco_diagM_history_fields", "out_his_diagM_v"),
    ("diagM_average_fields", 0, 0, "bool", "&croco_diagM_average_fields", "out_avg_diagM_u"),
    ("diagM_average_fields", 0, 1, "bool", "&croco_diagM_average_fields", "out_avg_diagM_v"),
    ("diags_vrt_history_fields", 0, 0, "bool", "&croco_diags_vrt_history_fields", "out_his_diags_vrt"),
    ("diags_vrt_average_fields", 0, 0, "bool", "&croco_diags_vrt_average_fields", "out_avg_diags_vrt"),
    ("diags_ek_history_fields", 0, 0, "bool", "&croco_diags_ek_history_fields", "out_his_diags_ek"),
    ("diags_ek_average_fields", 0, 0, "bool", "&croco_diags_ek_average_fields", "out_avg_diags_ek"),
    ("diags_pv_history_fields", 0, 0, "bool_array", "&croco_diags_pv_history_fields", "out_his_diags_pv_tracer"),
    ("diags_pv_average_fields", 0, 0, "bool_array", "&croco_diags_pv_average_fields", "out_avg_diags_pv_tracer"),
    ("diags_eddy_average_fields", 0, 0, "bool", "&croco_diags_eddy_average_fields", "out_avg_diags_eddy"),
    ("diagbioFlux_history_fields", 0, 0, "bool_array:NumFluxTerms", "&croco_diagbioFlux_history_fields", "out_his_diagbioFlux"),
    ("diagbioFlux_history_fields", 0, 0, "bool_array:NumVSinkTerms", "&croco_diagbioVSink_history_fields", "out_his_diagbioVSink"),
    ("diagbioGasExc_history_fields", 0, 0, "bool_array:NumGasExcTerms", "&croco_diagbioGasExc_history_fields", "out_his_diagbioGasExc"),
    ("diagbioFlux_average_fields", 0, 0, "bool_array:NumFluxTerms", "&croco_diagbioFlux_average_fields", "out_avg_diagbioFlux"),
    ("diagbioFlux_average_fields", 0, 0, "bool_array:NumVSinkTerms", "&croco_diagbioVSink_average_fields", "out_avg_diagbioVSink"),
    ("diagbioGasExc_average_fields", 0, 0, "bool_array:NumGasExcTerms", "&croco_diagbioGasExc_average_fields", "out_avg_diagbioGasExc"),
]

# Cards handled by dedicated code rather than the MAPPINGS table.
SPECIAL_CARDS = {"psource", "psource_ncfile"}

# Canonical namelist block order, taken from croco_full.nml. Output blocks are
# sorted to match this order; any block not listed here (should not happen —
# every nml_name produced by MAPPINGS/special handlers is covered) is appended
# at the end in encounter order.
CANONICAL_BLOCK_ORDER = (
    "&croco_title",
    "&croco_logfile",
    "&croco_testcase",
    "&croco_time_stepping",
    "&croco_time_stepping_nbq",
    "&croco_use_calendar",
    "&croco_s_coord",
    "&croco_initial",
    "&croco_grid",
    "&croco_forcing",
    "&croco_bulk_forcing",
    "&croco_climatology",
    "&croco_boundary",
    "&croco_bottom_forcing",
    "&croco_sponge",
    "&croco_wetdry",
    "&croco_wavedry",
    "&croco_nudging",
    "&croco_history",
    "&croco_averages",
    "&croco_restart",
    "&croco_surf",
    "&croco_surf_avg",
    "&croco_rho0",
    "&croco_bottom_drag",
    "&croco_gamma2",
    "&croco_lateral_visc",
    "&croco_tracer_diff2",
    "&croco_tracer_diff4",
    "&croco_vertical_mixing",
    "&croco_lin_eos",
    "&croco_diagnostics_ts",
    "&croco_diag_avg",
    "&croco_diag_mld_dens",
    "&croco_diagnosticsM",
    "&croco_diagM_avg",
    "&croco_diags_vrt",
    "&croco_diags_vrt_avg",
    "&croco_diags_ek",
    "&croco_diags_ek_avg",
    "&croco_diags_pv",
    "&croco_diags_pv_avg",
    "&croco_diags_eddy_avg",
    "&croco_diagnostics_bio",
    "&croco_diagbio_avg",
    "&croco_stations",
    "&croco_abl",
    "&croco_abl_averages",
    "&croco_abl_nudg_tra",
    "&croco_abl_nudg_dyn",
    "&croco_online",
    "&croco_wkb_boundary",
    "&croco_wkb_wwave",
    "&croco_wkb_roller",
    "&croco_wave_maker",
    "&croco_wave_offline",
    "&croco_biology",
    "&croco_bodyforce",
    "&croco_psource",
    "&croco_psource_ncfile_data",
    "&croco_psource_data",
    "&croco_psource_tracer",
    "&croco_sediments",
    "&croco_sediments_mustang",
    "&croco_substance",
    "&croco_obstruction",
    "&croco_xios_origin_date",
    "&croco_assimilation",
    "&croco_primary_history_fields",
    "&croco_primary_history_3d_fields",
    "&croco_primary_history_tracer_fields",
    "&croco_auxiliary_history_fields",
    "&croco_temperature_history_fields",
    "&croco_salinity_history_fields",
    "&croco_bulk_flux_history_fields",
    "&croco_bhflux_history_fields",
    "&croco_bwflux_history_fields",
    "&croco_bvf_history_fields",
    "&croco_hbl_history_fields",
    "&croco_lmd_bkpp_history_fields",
    "&croco_vis_coef_history_fields",
    "&croco_dif_coef_history_fields",
    "&croco_gls_history_fields",
    "&croco_biology_history_fields",
    "&croco_bio_nchlpzd_history_fields",
    "&croco_bio_bioebus_history_fields",
    "&croco_morphodyn_history_fields",
    "&croco_sediment_history_fields",
    "&croco_sediment_bfra_history_fields",
    "&croco_sediment_suspload_history_fields",
    "&croco_sediment_bedload_history_fields",
    "&croco_sediment_cohesive_history_fields",
    "&croco_bbl_history_fields",
    "&croco_wci_history_fields",
    "&croco_wci_history_3d_fields",
    "&croco_wave_history_fields",
    "&croco_abl_history_fields",
    "&croco_surf_history_fields",
    "&croco_stochastic_history_fields",
    "&croco_diag3D_history_fields",
    "&croco_diag2D_history_fields",
    "&croco_diagM_history_fields",
    "&croco_diags_vrt_history_fields",
    "&croco_diags_ek_history_fields",
    "&croco_diags_pv_history_fields",
    "&croco_diagbioFlux_history_fields",
    "&croco_diagbioVSink_history_fields",
    "&croco_diagbioGasExc_history_fields",
    "&croco_primary_average_fields",
    "&croco_primary_3d_average_fields",
    "&croco_primary_tracer_average_fields",
    "&croco_auxiliary_averages_fields",
    "&croco_temperature_averages_fields",
    "&croco_salinity_averages_fields",
    "&croco_bulk_flux_averages_fields",
    "&croco_bhflux_averages_fields",
    "&croco_bwflux_averages_fields",
    "&croco_bvf_averages_fields",
    "&croco_hbl_averages_fields",
    "&croco_lmd_bkpp_averages_fields",
    "&croco_vis_coef_averages_fields",
    "&croco_dif_coef_averages_fields",
    "&croco_gls_averages_fields",
    "&croco_biology_averages_fields",
    "&croco_bio_nchlpzd_averages_fields",
    "&croco_bio_bioebus_averages_fields",
    "&croco_morphodyn_averages_fields",
    "&croco_abl_averages_fields",
    "&croco_surf_average_fields",
    "&croco_wci_average_fields",
    "&croco_wci_average_3d_fields",
    "&croco_wave_average_fields",
    "&croco_diag3D_average_fields",
    "&croco_diag2D_average_fields",
    "&croco_diagM_average_fields",
    "&croco_diags_vrt_average_fields",
    "&croco_diags_ek_average_fields",
    "&croco_diags_pv_average_fields",
    "&croco_diags_eddy_average_fields",
    "&croco_diagbioFlux_average_fields",
    "&croco_diagbioVSink_average_fields",
    "&croco_diagbioGasExc_average_fields",
    "&croco_station_fields",
)


# ---------------------------------------------------------------------------
# CPP helpers
# ---------------------------------------------------------------------------


def _run_cpp(wrapper_src, include_dirs):
    """Run cpp -P -traditional on wrapper_src (passed via stdin). Return stdout."""
    cmd = ["cpp", "-P", "-traditional"] + [f"-I{d}" for d in include_dirs] + ["-"]
    r = subprocess.run(cmd, input=wrapper_src, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  cpp error:\n{r.stderr}", file=sys.stderr)
    return r.stdout


def cpp_define_active(cppdefs_path, define_name, extra_dirs=()):
    """Return True if define_name is set in cppdefs_path (via cpp).

    extra_dirs: additional -I directories (e.g. OCEAN/) so cppdefs.h's own
    includes like cppdefs_dev.h / set_global_definitions.h can be resolved
    even when cppdefs_path lives in an unrelated run directory.
    """
    cppdefs_abs = os.path.realpath(cppdefs_path)
    cppdefs_dir = os.path.dirname(cppdefs_abs)
    sentinel = f"__{define_name}_IS_ACTIVE__"
    wrapper = f'#include "{cppdefs_abs}"\n#ifdef {define_name}\n{sentinel}\n#endif\n'
    out = _run_cpp(wrapper, [cppdefs_dir, *extra_dirs])
    return sentinel in out


def preprocess_param(cppdefs_path, param_path, extra_dirs=()):
    """
    Run cpp on param.h using cppdefs.h for defines.
    Returns the preprocessed text (no CPP directives, resolved conditionals).

    extra_dirs: additional -I directories (e.g. OCEAN/) so cppdefs.h's own
    includes like cppdefs_dev.h / set_global_definitions.h can be resolved
    even when cppdefs_path/param_path live in an unrelated run directory.
    """
    cppdefs_abs = os.path.realpath(cppdefs_path)
    param_abs = os.path.realpath(param_path)
    wrapper = f'#include "{cppdefs_abs}"\n#include "{param_abs}"\n'
    return _run_cpp(wrapper, [os.path.dirname(param_abs), *extra_dirs])


def _candidate_value_lines(line):
    """Yield (name, old_value) for each simple `NAME = <literal>` assignment
    found on `line` (modern "integer/real, parameter :: NAME = VALUE" or
    legacy "parameter (NAME=VALUE, ...)" syntax). Lines whose value is an
    expression/macro reference (e.g. NNODES=NP_XI*NP_ETA) are not yielded —
    only plain numeric literals are considered fillable.
    """
    if line.strip().startswith("!"):
        return
    m = re.match(r"^\s*(?:integer|real)\s*,\s*parameter\s*::\s*(\w+)\s*=\s*([+-]?[\d.]+)\s*(?:!.*)?$", line)
    if m:
        yield m.group(1), m.group(2)
        return
    pm = re.search(r"parameter\s*\(([^)]*)\)", line, re.IGNORECASE)
    if pm:
        for assign in pm.group(1).split(","):
            am = re.match(r"\s*(\w+)\s*=\s*([+-]?[\d.]+)\s*$", assign.strip())
            if am:
                yield am.group(1), am.group(2)


def fill_param_template(template_path, value_source_path, cppdefs_path, extra_dirs=()):
    """
    Take template_path (e.g. {--src}/param.h) as the structural template —
    same line count, comments, #ifdef sections, and #include lines — and
    fill in each simple `NAME = <literal>` value from value_source_path
    (e.g. param_old.h), resolved against cppdefs_path. Only lines that are
    currently active per cppdefs_path are touched; if a name has no
    resolved value in value_source_path, the template's existing value is
    left as-is. Expression-valued lines (e.g. NNODES=NP_XI*NP_ETA) and
    #include lines are never modified.
    """
    template_text = open(template_path).read()

    preprocessed_values = preprocess_param(cppdefs_path, value_source_path, extra_dirs=extra_dirs)
    value_dict = read_param_text(preprocessed_values)

    # Determine which template lines are active for this cppdefs case: swap
    # the template's own #include lines for placeholders (so cpp doesn't try
    # to resolve/inline them), run cpp, and check which original lines
    # survive verbatim in the output.
    include_lines = re.findall(r'^[ \t]*#[ \t]*include\b.*$', template_text, re.MULTILINE)
    protected_text = template_text
    for i, line in enumerate(include_lines):
        protected_text = protected_text.replace(line, f"__TPL_INCLUDE_{i}__", 1)

    cppdefs_abs = os.path.realpath(cppdefs_path)
    template_dir = os.path.dirname(os.path.realpath(template_path))
    wrapper = f'#include "{cppdefs_abs}"\n{protected_text}\n'
    cpp_out = _run_cpp(wrapper, [template_dir, os.path.dirname(cppdefs_abs), *extra_dirs])
    active_lines = {l.strip() for l in cpp_out.splitlines() if l.strip()}

    out_lines = []
    for line in template_text.splitlines():
        new_line = line
        if line.strip() in active_lines:
            for name, old_val in _candidate_value_lines(line):
                if name in value_dict and str(value_dict[name]) != old_val:
                    new_line = re.sub(
                        rf"(\b{re.escape(name)}\s*=\s*){re.escape(old_val)}\b",
                        rf"\g<1>{value_dict[name]}",
                        new_line,
                        count=1,
                    )
        out_lines.append(new_line)

    return "\n".join(out_lines) + "\n"


def is_full_cppdefs(path, min_cases=10):
    """
    Return True if the file is a full, multi-case cppdefs.h master template
    rather than an already-extracted, single-case file.

    Detected structurally (name-agnostic, independent of which case if any
    is currently #define'd vs #undef'd):
      1. A case-selector block: a run of bare "#undef NAME" / "#define NAME"
         lines (no value) before any #if/#ifdef/#ifndef block.
      2. A matching top-level "#if defined X" / "#elif defined Y" / ...
         chain right after, whose branch names overlap with the selector
         block.
    Both must list at least `min_cases` names to count as "full" — a
    genuinely single-case file won't have anywhere near that many.
    """
    lines = open(path).read().splitlines()

    # 1. Case-selector names: bare #undef/#define lines before the first
    #    #if/#ifdef/#ifndef directive.
    selector_names = set()
    chain_start = None
    for i, line in enumerate(lines):
        s = line.strip()
        if re.match(r"^#\s*(?:if|ifdef|ifndef)\b", s):
            chain_start = i
            break
        m = re.match(r"^#\s*(?:undef|define)\s+(\w+)\s*(?:/\*.*)?$", s)
        if m:
            selector_names.add(m.group(1))
    if chain_start is None:
        return False

    # 2. Branches of that first top-level #if/#elif chain only (depth-tracked
    #    so nested conditionals inside a branch don't get counted).
    chain_branches = set()
    depth = 0
    for line in lines[chain_start:]:
        s = line.strip()
        if re.match(r"^#\s*(?:if|ifdef|ifndef)\b", s):
            if depth == 0:
                m = re.match(r"^#\s*if\s+defined\s+(\w+)\b", s)
                if m:
                    chain_branches.add(m.group(1))
            depth += 1
        elif re.match(r"^#\s*elif\b", s):
            if depth == 1:
                m = re.match(r"^#\s*elif\s+defined\s+(\w+)\b", s)
                if m:
                    chain_branches.add(m.group(1))
        elif re.match(r"^#\s*endif\b", s):
            depth -= 1
            if depth == 0:
                break

    shared = selector_names & chain_branches
    return len(selector_names) >= min_cases and len(chain_branches) >= min_cases and len(shared) >= min_cases


def extract_case_section(cppdefs_path):
    """
    Extract the case-specific block from cppdefs_path verbatim, preserving
    the original order and nested #ifdef structure.

    Returns (case_key, text) where text is ready to write as cppdefs_config.h:
        # define <CASE_KEY>           <- synthetic first line
        <block content, as-is>        <- inside #elif defined <CASE>
        <trailing #include lines>     <- e.g. #include "cppdefs_dev.h"
    """
    lines = open(cppdefs_path).readlines()

    # 1. Detect the active case key: the single bare "#define NAME" line in
    #    the case-selector header, i.e. before the #if/#elif chain starts.
    #    Two conventions are supported:
    #      - direct edit: the case is selected in place among the
    #        "#undef CASE / ... / #define CASE" selector list
    #        (e.g. "#define REGIONAL" or "#define UPWELLING").
    #      - patch_select_case style: "#undef REGIONAL" followed by a
    #        synthetic "#define <CASE>" line appended after it.
    #    In both cases the active key is simply the last bare #define found
    #    before the chain start (there must be exactly one, per cppdefs.h's
    #    "exactly one activated test case" rule).
    chain_start = None
    for i, line in enumerate(lines):
        if re.match(r"^#\s*if\s+defined\s+\w+\b", line.strip()):
            chain_start = i
            break
    header_lines = lines[:chain_start] if chain_start is not None else lines

    case_key = None
    for line in header_lines:
        m = re.match(r"^#\s*define\s+(\w+)\s*(?:/\*.*)?$", line.strip())
        if m:
            case_key = m.group(1)
    if case_key is None:
        return None, ""

    # 2. Find the start of the case block (#elif defined CASE or #if for REGIONAL)
    if case_key == "REGIONAL":
        start_re = re.compile(r"^#if\s+defined\s+REGIONAL\b")
    else:
        start_re = re.compile(rf"^#elif\s+defined\s+{re.escape(case_key)}\b")

    start_idx = None
    for i, line in enumerate(lines):
        if start_re.match(line.rstrip()):
            start_idx = i
            break
    if start_idx is None:
        return case_key, ""

    # 3. Extract block until next #elif/#else/#endif at depth 0,
    #    skipping the opening /* title comment */ if present.
    #    Note: single-line /* ... */ comments on the SAME line must not
    #    leave in_c_comment=True (the closing */ is on the same line).
    block = []
    depth = 0
    skip_opening_comment = True
    in_c_comment = False

    for line in lines[start_idx + 1 :]:
        s = line.rstrip()
        stripped = s.strip()

        if skip_opening_comment:
            if not in_c_comment:
                if stripped.startswith("/*"):
                    # Only enter multi-line mode when */ is NOT on the same line
                    if "*/" not in stripped:
                        in_c_comment = True
                    # Either way, this comment line belongs to the opening header
                    continue
                else:
                    # First non-comment line: end of opening-skip phase
                    skip_opening_comment = False
                    # fall through to process this line normally
            else:
                # Inside a multi-line opening comment
                if "*/" in stripped:
                    in_c_comment = False
                continue

        if depth == 0:
            if re.match(r"^#\s*(elif|else)\b", s):
                break
            if re.match(r"^#\s*endif\b", s):
                break

        if re.match(r"^#\s*if(def|ndef)?\b", s):
            depth += 1
        elif re.match(r"^#\s*endif\b", s):
            depth -= 1

        block.append(line)

    # 4. Collect trailing #include lines that follow the main closing #endif.
    # In reverse: includes/blanks appear BEFORE the #endif, so collect them
    # until we hit the #endif, then stop.
    trailing = []
    for line in reversed(lines):
        s = line.strip()
        if re.match(r"^#\s*endif\b", s):
            break  # reached the closing #endif — done
        if s.startswith("#include") or s == "":
            trailing.insert(0, line)

    # 5. Assemble: synthetic case-key line + block + trailing includes
    out = [f"# define {case_key}\n"] + block + trailing
    return case_key, "".join(out)


# ---------------------------------------------------------------------------
# param file reader  (handles both simple and preprocessed Fortran formats)
# ---------------------------------------------------------------------------


def read_param_text(text):
    """
    Parse Fortran parameter statements from preprocessed text and return a
    dict of resolved integer values.  Evaluates arithmetic expressions
    (e.g. NT = itemp + ntrc_salt + ntrc_pas + …). Handles both the modern
    "integer, parameter :: NAME = EXPR" format and the legacy
    "parameter (NAME=EXPR, …)" format — including modern declarations whose
    value is an expression rather than a literal (e.g. param_dev.h's
    "integer, parameter :: NT = itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed").
    """
    # Strip Fortran "!" comments first, so commented-out example/alternative
    # lines (e.g. "!     parameter (N_sl=40)") aren't mistaken for live code.
    text = re.sub(r"!.*", "", text)

    assignments = []  # ordered (name, expr) pairs from both syntaxes

    for m in re.finditer(r"integer\s*,\s*parameter\s*::\s*(\w+)\s*=\s*(\S+)", text):
        assignments.append((m.group(1), m.group(2).strip()))

    for stmt in re.finditer(r"parameter\s*\(([^)]+)\)", text, re.IGNORECASE):
        for assign in stmt.group(1).split(","):
            m = re.match(r"\s*(\w+)\s*=\s*(.+)", assign.strip())
            if m:
                assignments.append((m.group(1), m.group(2).strip()))

    # Resolve with retries: some expressions reference names assigned later
    # in file order (e.g. across separate #ifdef branches).
    params = {}
    changed = True
    while changed:
        changed = False
        for name, expr in assignments:
            if name in params:
                continue
            try:
                params[name] = _eval_fortran_expr(expr, params)
                changed = True
            except Exception:
                pass
    return params


def _eval_fortran_expr(expr, known):
    def replace_name(m):
        n = m.group(0)
        if n in known:
            return str(known[n])
        raise ValueError(f"Unknown: {n}")

    resolved = re.sub(r"[A-Za-z_]\w*", replace_name, expr)
    if not re.match(r"^[\d\s+\-*/()]+$", resolved):
        raise ValueError(f"Unsafe: {resolved}")
    return int(eval(resolved))  # noqa: S307 – input sanitised above


def read_param_file(path):
    return read_param_text(open(path).read())


# ---------------------------------------------------------------------------
# Output generators
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .in file parser
# ---------------------------------------------------------------------------


def parse_in_file(path):
    """Return (cards, headers).
    cards[key]   = list of value lines (stripped).
    headers[key] = list of label tokens from the header comment on the card
                   key line (the text after "key:" on that same line).
    """
    cards = {}
    headers = {}
    current_key = None
    current_lines = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if stripped == "":
                continue
            m = re.match(r"^([A-Za-z][A-Za-z0-9_\-]*):", line)
            if m:
                if current_key is not None:
                    cards[current_key] = current_lines
                current_key = m.group(1)
                headers[current_key] = line[m.end() :].split()
                current_lines = []
                continue
            if current_key is not None:
                current_lines.append(stripped)
    if current_key is not None:
        cards[current_key] = current_lines
    return cards, headers


# ---------------------------------------------------------------------------
# Value formatter
# ---------------------------------------------------------------------------


def format_value(raw, typ):
    try:
        if typ == "bool":
            if raw.upper() in ("T", ".TRUE."):
                return ".true."
            if raw.upper() in ("F", ".FALSE."):
                return ".false."
            return raw
        if typ == "int":
            return str(int(float(raw.lower().replace("d", "e"))))
        if typ == "float":
            return raw
        if typ in ("str", "str_line"):
            return f'"{raw}"'
        return raw
    except Exception:
        return raw


def expand_repeat(tokens):
    expanded = []
    for t in tokens:
        if "*" in t:
            count_str, val_str = t.split("*", 1)
            try:
                expanded.extend([val_str] * int(count_str))
            except ValueError:
                expanded.append(t)
        else:
            expanded.append(t)
    return expanded


def format_float_array(expanded, NT):
    if not expanded:
        return ""
    if NT is not None and len(expanded) < NT:
        expanded = expanded + [expanded[-1]] * (NT - len(expanded))
    if NT is not None and len(expanded) > NT:
        expanded = expanded[:NT]
    unique = set(v.lower() for v in expanded)
    if len(unique) == 1:
        count = NT if NT is not None else len(expanded)
        return f"{count}*{expanded[0]}"
    return ", ".join(expanded)


def format_bool_array(expanded, NT):
    if not expanded:
        return ""
    normalized = [".true." if v.upper() in ("T", ".TRUE.") else ".false." for v in expanded]
    if NT is not None and len(normalized) < NT:
        normalized = normalized + [normalized[-1]] * (NT - len(normalized))
    if NT is not None and len(normalized) > NT:
        normalized = normalized[:NT]
    # NT unknown: skip uniform-true arrays (module default covers them);
    # for mixed values we cannot safely emit without knowing the real size.
    if NT is None:
        return ""
    unique = set(normalized)
    if len(unique) == 1:
        return f"{NT}*{normalized[0]}"
    return ", ".join(normalized)


# ---------------------------------------------------------------------------
# psource / psource_ncfile special handlers
# ---------------------------------------------------------------------------


def _bool_token(t):
    return t.upper() in ("T", "F", ".TRUE.", ".FALSE.")


def _fmt_bool(t):
    return ".true." if t.upper() in ("T", ".TRUE.") else ".false."


def _is_numeric(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def _split_rest(tokens, nt):
    """Split tokens[4:] into (Lsrc_tokens, Tsrc0_tokens) for one source line.

    Lsrc tokens are T/F booleans; Tsrc0 tokens are numeric.  If NT is known
    both lists are trimmed to NT elements.  Trailing non-numeric labels
    (e.g. river names appended after Tsrc) are discarded.
    """
    rest = tokens[4:]
    # find where booleans end
    split = len(rest)
    for i, t in enumerate(rest):
        if not _bool_token(t):
            split = i
            break
    lsrc = rest[:split]
    tsrc = [t for t in rest[split:] if _is_numeric(t)]
    if nt is not None:
        lsrc = (lsrc + [".false."] * nt)[:nt]
        tsrc = (tsrc + ["0."] * nt)[:nt]
    return lsrc, tsrc


def _read_real_param(text, name):
    """Find a real-valued `name = VALUE` parameter declaration in preprocessed
    Fortran text, supporting both the modern (`real, parameter :: name = value`)
    and legacy (`parameter (name=value)`) syntax. Returns the value as a
    string (original formatting preserved) or None if not declared.
    """
    text = re.sub(r"!.*", "", text)  # strip comments; see read_param_text
    m = re.search(rf"real\s*,\s*parameter\s*::\s*{name}\s*=\s*([\d.eEdD+-]+)", text)
    if m:
        return m.group(1)
    m = re.search(rf"parameter\s*\(\s*{name}\s*=\s*([\d.eEdD+-]+)\s*\)", text)
    if m:
        return m.group(1)
    return None


def _build_wetdry_blocks(cards, params):
    """Return {nml_block_name: [entry_lines]} for &croco_wetdry / &croco_wavedry.

    D_wetdry/D_wavedry used to be compile-time CPP constants hard-coded per
    test case in OCEAN/param.h, now migrated to runtime namelist defaults
    (&croco_wetdry / &croco_wavedry in croco_namelist.F90). There is no .in
    card for them, so instead of guessing a value we read whatever param.h
    (preprocessed against the active cppdefs.h) actually declares for them
    via params["D_wetdry"]/params["D_wavedry"] (set in main()). If param.h
    doesn't declare them (e.g. WET_DRY inactive, or already fully migrated),
    no block is emitted and the namelist default in croco_namelist.F90 applies.
    """
    blocks = {}
    d_wetdry = params.get("D_wetdry")
    if d_wetdry is not None:
        blocks["&croco_wetdry"] = [f"  D_wetdry = {d_wetdry}"]
    d_wavedry = params.get("D_wavedry")
    if d_wavedry is not None:
        blocks["&croco_wavedry"] = [f"  D_wavedry = {d_wavedry}"]
    return blocks


def _build_psource_blocks(cards, params):
    """Return {nml_block_name: [entry_lines]} for psource/psource_ncfile cards."""
    # NT is undeclared (None) for SOLVE3D-less (2D-only) configs; treat as 0 tracers.
    NT = params.get("NT") or 0
    blocks = {}

    def add(name, line):
        blocks.setdefault(name, []).append(line)

    # Use the cpp-probed flag when available; fall back to card presence otherwise.
    use_ncfile = params.get("_psource_ncfile", "psource_ncfile" in cards)

    if not use_ncfile and "psource" in cards:
        vlines = [line for line in cards["psource"] if line.strip()]
        if not vlines:
            return blocks
        Nsrc = int(expand_repeat(vlines[0].split())[0])

        Isrc, Jsrc, Dsrc, Qbar = [], [], [], []
        lsrc_per_src, tsrc0_per_src = [], []
        for i in range(Nsrc):
            if i + 1 >= len(vlines):
                break
            toks = expand_repeat(vlines[i + 1].split())
            Isrc.append(toks[0])
            Jsrc.append(toks[1])
            Dsrc.append(toks[2])
            Qbar.append(toks[3])
            lsrc, tsrc0 = _split_rest(toks, NT)
            lsrc_per_src.append(lsrc)
            tsrc0_per_src.append(tsrc0)

        add("&croco_psource", f"  psource_Nsrc = {Nsrc}")
        for i in range(len(Isrc)):
            n = i + 1
            add(
                "&croco_psource_data",
                f"  psource_Isrc({n:>4}) = {Isrc[i]:>5},"
                + f" psource_Jsrc({n:>4}) = {Jsrc[i]:>5},"
                + f" psource_Dsrc({n:>4}) = {Dsrc[i]:>5},"
                + f" psource_Qbar({n:>4}) = {Qbar[i]:>5}",
            )
        if NT > 0:
            for isrc, lsrc in enumerate(lsrc_per_src, start=1):
                for itr, val in enumerate(lsrc, start=1):
                    add("&croco_psource_tracer", f"  psource_Lsrc({isrc},{itr}) = {val}")
            for isrc, tsrc in enumerate(tsrc0_per_src, start=1):
                for itr, val in enumerate(tsrc, start=1):
                    add("&croco_psource_tracer", f"  psource_Tsrc0({isrc},{itr}) = {val}")

    elif use_ncfile and "psource_ncfile" in cards:
        vlines = [line for line in cards["psource_ncfile"] if line.strip()]
        if len(vlines) < 2:
            return blocks
        qbarname = vlines[0].strip()
        Nsrc = int(expand_repeat(vlines[1].split())[0])

        Isrc, Jsrc, Dsrc, Qbardir = [], [], [], []
        lsrc_per_src, tsrc0_per_src = [], []
        for i in range(Nsrc):
            if i + 2 >= len(vlines):
                break
            toks = expand_repeat(vlines[i + 2].split())
            Isrc.append(toks[0])
            Jsrc.append(toks[1])
            Dsrc.append(toks[2])
            Qbardir.append(toks[3])
            lsrc, tsrc0 = _split_rest(toks, NT)
            lsrc_per_src.append(lsrc)
            tsrc0_per_src.append(tsrc0)

        add("&croco_psource", f"  psource_Nsrc = {Nsrc}")
        add("&croco_psource_ncfile_data", f'  psource_qbarname = "{qbarname}"')
        for i in range(len(Isrc)):
            n = i + 1
            add(
                "&croco_psource_ncfile_data",
                f"  psource_Isrc({n:>4}) = {Isrc[i]:>5},"
                + f" psource_Jsrc({n:>4}) = {Jsrc[i]:>5},"
                + f" psource_Dsrc({n:>4}) = {Dsrc[i]:>5},"
                + f" psource_qbardir({n:>4}) = {Qbardir[i]:>5}",
            )

        if NT > 0:
            for isrc, lsrc in enumerate(lsrc_per_src, start=1):
                for itr, val in enumerate(lsrc, start=1):
                    add("&croco_psource_tracer", f"  psource_Lsrc({isrc},{itr}) = {val}")
            for isrc, tsrc in enumerate(tsrc0_per_src, start=1):
                for itr, val in enumerate(tsrc, start=1):
                    add("&croco_psource_tracer", f"  psource_Tsrc0({isrc},{itr}) = {val}")

    return blocks


# ---------------------------------------------------------------------------
# Labeled_bool fallback handler
# ---------------------------------------------------------------------------


def _build_labeled_bool_fallback_blocks(cards, headers, params, mappings):
    """Return {nml_block_name: [entry_lines]} for labeled_bool fields that are
    absent from the .in header but declared in the namelist.

    Operates on MAPPINGS rows that carry two extra fields (cpp_key, default):
      (card, line, label, "labeled_bool", nml_block, nml_var, cpp_key, default)
    - If the label IS in the .in header, MAPPINGS labeled_bool already handled it → skip.
    - If NOT in the header and the CPP guard is active → emit the namelist default.
    """
    blocks = {}
    for row in mappings:
        if len(row) < 8:
            continue
        card, _, label, typ, nml_block, nml_var, cpp_key, default = row
        if typ != "labeled_bool":
            continue
        if card not in cards:
            continue
        hdr = [t.lower() for t in headers.get(card, [])]
        if label.lower() in hdr:
            continue  # already handled by MAPPINGS labeled_bool
        if cpp_key and not params.get(cpp_key, False):
            continue
        blocks.setdefault(nml_block, []).append(f"  {nml_var} = {default}")
    return blocks


# ---------------------------------------------------------------------------
# NML builder
# ---------------------------------------------------------------------------
# Sediment CPP-conditional output field handler
# ---------------------------------------------------------------------------


def _build_sediment_extra_history_blocks(cards, params):
    """Convert CPP-conditional sediment history fields: dflx/eflx (SUSPLOAD),
    bdlu/bdlv (BEDLOAD), btcr (MIXED_BED or COHESIVE_BED).
    Their position in the boolean array shifts with NST and active CPP keys,
    so they cannot be expressed as fixed-offset MAPPINGS entries.

    Returns {nml_block_name: [entry_lines]}.
    """
    blocks = {}
    if "sediment_history_fields" not in cards:
        return blocks
    NST = params.get("NST", None)
    if NST is None:
        return blocks
    vlines = [line for line in cards["sediment_history_fields"] if line.strip()]
    if not vlines:
        return blocks
    tokens = expand_repeat(vlines[0].split())

    offset = 3 + NST  # past: athk bthk bpor bfra(NST)

    if params.get("_suspload", False):
        dflx = tokens[offset : offset + NST]
        eflx = tokens[offset + NST : offset + 2 * NST]
        v1 = format_bool_array(dflx, NST)
        v2 = format_bool_array(eflx, NST)
        entries = []
        if v1:
            entries.append(f"  out_his_sed_dflx = {v1}")
        if v2:
            entries.append(f"  out_his_sed_eflx = {v2}")
        if entries:
            blocks["&croco_sediment_suspload_history_fields"] = entries
        offset += 2 * NST

    if params.get("_bedload", False):
        bdlu = tokens[offset : offset + NST]
        bdlv = tokens[offset + NST : offset + 2 * NST]
        v1 = format_bool_array(bdlu, NST)
        v2 = format_bool_array(bdlv, NST)
        entries = []
        if v1:
            entries.append(f"  out_his_sed_bdlu = {v1}")
        if v2:
            entries.append(f"  out_his_sed_bdlv = {v2}")
        if entries:
            blocks["&croco_sediment_bedload_history_fields"] = entries
        offset += 2 * NST

    if params.get("_mixed_or_cohesive", False) and offset < len(tokens):
        val = format_value(tokens[offset], "bool")
        blocks["&croco_sediment_cohesive_history_fields"] = [f"  out_his_sed_btcr = {val}"]

    return blocks


# ---------------------------------------------------------------------------

# Test-case keys that are always CPP-only (no ANA_GRID/ANA_VMIX dispatch).
# For these, testcase_name is still set so croco_namelist_check can validate it.
_REALISTIC_CASES = frozenset({"REGIONAL", "COASTAL"})


def _detect_case_key_from_simplified(cppdefs_path):
    """Return the first bare #define key from a simplified single-case cppdefs.h."""
    try:
        with open(cppdefs_path) as fh:
            for line in fh:
                m = re.match(r"^#\s*define\s+(\w+)", line)
                if m:
                    return m.group(1)
    except OSError:
        pass
    return None


def _build_testcase_block(params):
    """Return {nml_name: [entry_lines]} for &croco_testcase when testcase_name is known."""
    name = params.get("_testcase_name", "")
    if not name:
        return {}
    return {"&croco_testcase": [f"  testcase_name = '{name}'"]}


def build_nml(cards, headers, mappings, params):
    NT = params.get("NT", None)
    nml_entries = {}
    nml_order = []  # encounter order, used as fallback for blocks not in CANONICAL_BLOCK_ORDER

    def add_entry(nml_name, entry_line):
        if nml_name not in nml_entries:
            nml_entries[nml_name] = []
            nml_order.append(nml_name)
        nml_entries[nml_name].append(entry_line)

    for row in mappings:
        # Some labeled_bool rows carry 2 extra fields (cpp_key, default) used
        # only by _build_labeled_bool_fallback_blocks(); ignore them here.
        card_name, line_idx, pos_idx, typ, nml_name, nml_var = row[:6]
        if card_name not in cards:
            continue
        vlines = [line for line in cards[card_name] if line.strip()]
        if line_idx >= len(vlines):
            continue
        tokens = expand_repeat(vlines[line_idx].split())

        if typ == "labeled_bool":
            hdr = [lbl.lower() for lbl in headers.get(card_name, [])]
            target = str(pos_idx).lower()
            if target not in hdr:
                continue
            actual_pos = hdr.index(target)
            if actual_pos >= len(tokens):
                continue
            val = format_value(tokens[actual_pos], "bool")
            add_entry(nml_name, f"  {nml_var} = {val}")
            continue

        if typ == "float_array":
            raw_tokens = tokens[pos_idx:]
            if not raw_tokens:
                continue
            val = format_float_array(expand_repeat(raw_tokens), NT)
            if not val:
                continue
            add_entry(nml_name, f"  {nml_var} = {val}")
            continue

        if typ == "bool_array" or typ.startswith("bool_array:"):
            size_key = typ.split(":")[1] if ":" in typ else "NT"
            N = params.get(size_key, None)
            raw_tokens = tokens[pos_idx:]
            if not raw_tokens:
                continue
            val = format_bool_array(expand_repeat(raw_tokens), N)
            if not val:
                continue
            add_entry(nml_name, f"  {nml_var} = {val}")
            continue

        if typ == "str_line":
            raw = vlines[line_idx].strip()
        else:
            if pos_idx >= len(tokens):
                continue
            raw = tokens[pos_idx]

        if typ == "str" and re.match(r"^\d{4}[-/]\d{2}[-/]\d{2}$", raw) and pos_idx + 1 < len(tokens):
            raw = raw.replace("/", "-") + " " + tokens[pos_idx + 1]

        val = format_value(raw, typ)
        add_entry(nml_name, f"  {nml_var} = {val}")

    for nml_name, entry_lines in _build_wetdry_blocks(cards, params).items():
        for entry_line in entry_lines:
            add_entry(nml_name, entry_line)
    for nml_name, entry_lines in _build_psource_blocks(cards, params).items():
        for entry_line in entry_lines:
            add_entry(nml_name, entry_line)
    for nml_name, entry_lines in _build_sediment_extra_history_blocks(cards, params).items():
        for entry_line in entry_lines:
            add_entry(nml_name, entry_line)
    for nml_name, entry_lines in _build_labeled_bool_fallback_blocks(cards, headers, params, mappings).items():
        for entry_line in entry_lines:
            add_entry(nml_name, entry_line)
    for nml_name, entry_lines in _build_testcase_block(params).items():
        for entry_line in entry_lines:
            add_entry(nml_name, entry_line)

    # Order blocks to match croco_full.nml; anything unlisted falls back to
    # the order it was first encountered (should not normally happen).
    ordered_names = [n for n in CANONICAL_BLOCK_ORDER if n in nml_entries]
    leftover = [n for n in nml_order if n not in CANONICAL_BLOCK_ORDER]
    if leftover:
        print(f"Note: blocks not in canonical order, appended at end: {', '.join(leftover)}", file=sys.stderr)
    ordered_names += leftover

    lines = []
    for nml_name in ordered_names:
        lines.append(nml_name)
        lines.extend(nml_entries[nml_name])
        lines.append("/")
        lines.append("")
    lines += ["", "! End of namelist", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    default_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "OCEAN")

    parser = argparse.ArgumentParser(
        description="Convert a CROCO croco.in + cppdefs.h + param.h into croco.nml + cppdefs_config.h + param_config.h.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Takes these 3 inputs and an optional one to specified CROCO source path 
if the script is used elsewhere.
""",
    )
    parser.add_argument("croco_in", help="Path to the croco.in file")
    parser.add_argument("cppdefs", help="Path to cppdefs.h")
    parser.add_argument("param", help="Path to param.h")
    parser.add_argument("--output_namelist", default="croco.nml", help="Output .nml path (default: croco.nml)")
    parser.add_argument("--output_cppdef", default="cppdefs_config.h", help="Output cppdefs path (default: cppdefs_config.h)")
    parser.add_argument("--output_param", default="param_config.h", help="Output param path (default: param_config.h)")
    parser.add_argument(
        "--src",
        default=default_src,
        help="Path to the CROCO OCEAN/ source dir, used as an extra "
        f"include path for cpp (default: {default_src}). Needed "
        "to resolve cppdefs.h's own includes (cppdefs_dev.h, "
        "set_global_definitions.h) when cppdefs.h/param.h live "
        "in a separate run directory.",
    )
    args = parser.parse_args()

    for label, path in (("croco.in", args.croco_in), ("cppdefs.h", args.cppdefs), ("param.h", args.param)):
        if not os.path.isfile(path):
            sys.exit(f"Error: {label} file '{path}' not found.")

    src_dirs = [args.src] if os.path.isdir(args.src) else []

    out_nml = args.output_namelist
    out_cppdefs = args.output_cppdef
    out_param = args.output_param
    for out_path in (out_nml, out_cppdefs, out_param):
        out_dir = os.path.dirname(os.path.abspath(out_path))
        os.makedirs(out_dir, exist_ok=True)

    print(f"Input croco.in : {args.croco_in}")
    print(f"Input cppdefs  : {args.cppdefs}")
    print(f"Input param.h  : {args.param}")
    print()

    # ---- param.h: preprocess with cpp, resolved against cppdefs.h -----------
    # `preprocessed`/`params` fully resolve every #ifdef (including whatever
    # param.h's own #include's, e.g. param_dev.h, declare) so that NT/NST/
    # D_wetdry/etc. can be extracted below. The written param_config.h is
    # different: it's {--src}/param.h used as a structural template (same
    # line count, comments, #ifdef sections, #include lines), with each
    # simple literal value filled in from args.param wherever cppdefs.h
    # makes that line active; values absent from args.param keep the
    # template's own value.
    preprocessed = preprocess_param(args.cppdefs, args.param, extra_dirs=src_dirs)
    params = read_param_text(preprocessed)
    template_path = os.path.join(args.src, "param.h")
    with open(out_param, "w") as f:
        f.write(fill_param_template(template_path, args.param, args.cppdefs, extra_dirs=src_dirs))
    print(f"Output param   : {out_param}  (NT={params.get('NT', '?')})")

    # ---- cppdefs.h: always written out as a simplified, single-case file ----
    # If the input is a "full" cppdefs.h (multi-case #if/#elif chain), extract
    # just the active case's block. If it's already simplified, copy through.
    if is_full_cppdefs(args.cppdefs):
        case_key, section_text = extract_case_section(args.cppdefs)
        if section_text:
            with open(out_cppdefs, "w") as f:
                f.write(section_text)
            print(f"Output cppdefs : {out_cppdefs}  (case={case_key}, extracted from full cppdefs)")
        else:
            sys.exit(f"Error: could not extract a case section from full cppdefs '{args.cppdefs}'.")
    else:
        case_key = _detect_case_key_from_simplified(args.cppdefs)
        shutil.copyfile(args.cppdefs, out_cppdefs)
        print(f"Output cppdefs : {out_cppdefs}  (already simplified, case={case_key})")

    # Inject testcase_name for non-realistic configurations.
    if case_key and case_key not in _REALISTIC_CASES:
        params["_testcase_name"] = case_key
    else:
        params["_testcase_name"] = ""

    # ---- Check CPP flags that affect special-card handling -------------------
    params["_psource_ncfile"] = cpp_define_active(args.cppdefs, "PSOURCE_NCFILE", extra_dirs=src_dirs)
    params["_suspload"] = cpp_define_active(args.cppdefs, "SUSPLOAD", extra_dirs=src_dirs)
    params["_bedload"] = cpp_define_active(args.cppdefs, "BEDLOAD", extra_dirs=src_dirs)
    params["_mixed_or_cohesive"] = cpp_define_active(args.cppdefs, "MIXED_BED", extra_dirs=src_dirs) or cpp_define_active(
        args.cppdefs, "COHESIVE_BED", extra_dirs=src_dirs
    )
    params["_morphodyn"] = cpp_define_active(args.cppdefs, "MORPHODYN", extra_dirs=src_dirs)

    params["D_wetdry"] = _read_real_param(preprocessed, "D_wetdry")
    params["D_wavedry"] = _read_real_param(preprocessed, "D_wavedry")

    # ---- Convert .in → .nml --------------------------------------------------
    cards, headers = parse_in_file(args.croco_in)

    mapped_cards = {row[0] for row in MAPPINGS} | SPECIAL_CARDS
    present_cards = set(cards.keys())
    missing = mapped_cards - present_cards
    unmapped = present_cards - mapped_cards
    if missing:
        print(f"Skipped (no .in card)  : {', '.join(sorted(missing))}")
    if unmapped:
        print(f"Ignored (no mapping)   : {', '.join(sorted(unmapped))}")

    nml_text = build_nml(cards, headers, MAPPINGS, params)
    with open(out_nml, "w") as f:
        f.write(nml_text)
    print(f"Output .nml    : {out_nml}")


if __name__ == "__main__":
    main()
