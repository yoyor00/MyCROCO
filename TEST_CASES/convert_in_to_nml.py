#!/usr/bin/env python3
"""
convert_in_to_nml.py
====================
Convert a CROCO `.in` file into a Fortran namelist `.nml` file.
Optionally also produces a clean `param_<Name>.h` and `cppdefs_<Name>.h`.

Usage
-----
    # Basic: auto-discovers param_<Name>.h next to the .in file
    python convert_in_to_nml.py croco.in.Estuary

    # Full: provide cppdefs + param, get all three outputs
    python convert_in_to_nml.py croco.in.Estuary \\
        --cppdefs cppdefs.h --param ../OCEAN/param.h

    # Explicit output paths
    python convert_in_to_nml.py croco.in.Estuary \\
        --cppdefs cppdefs.h --param ../OCEAN/param.h \\
        -o croco_Estuary.nml \\
        --out-param param_Estuary.h \\
        --out-cppdefs cppdefs_Estuary.h

Outputs
-------
    .nml         Fortran namelist.  Array variables (tnu2, Akt_bak, …) use
                 compact Fortran repeat notation  NT*value  when all values
                 are equal, where NT is resolved from the param file.

    param_<N>.h  Preprocessed param.h: only the lines active for this case,
                 blank lines and Fortran comment lines removed.

    cppdefs_<N>.h  Active CPP defines for this case, one #define per line

Extending the mapping
---------------------
Add a row to MAPPINGS:

    # card_name       line  pos  type      nml_name                  nml_var
    ("time_stepping",   0,   0,  "int",   "&croco_time_stepping",   "ntimes"),

  card_name : key that starts the card in the .in file  (before the colon)
  line      : 0-based index of the value line inside that card block
  pos       : 0-based token index on that line (ignored for float_array)
  type      : "int" | "float" | "float_array" | "bool" | "str" | "str_line"
  nml_name  : namelist block name (with leading &)
  nml_var   : variable name inside that block

Multiple rows sharing the same nml_name are merged into one block.
Cards absent from the .in file are silently skipped (no block written).
"""

import re
import sys
import os
import argparse
import subprocess

# ---------------------------------------------------------------------------
# Mapping table
# ---------------------------------------------------------------------------
MAPPINGS = [
    ("title",              0,   0,  "str_line", "&croco_title",               "title"),
    ("logfile",            0,   0,  "str",   "&croco_logfile",             "logname"),
    ("time_stepping",      0,   0,  "int",   "&croco_time_stepping",       "ntimes"),
    ("time_stepping",      0,   1,  "float", "&croco_time_stepping",       "dt"),
    ("time_stepping",      0,   2,  "int",   "&croco_time_stepping",       "ndtfast"),
    ("time_stepping",      0,   3,  "int",   "&croco_time_stepping",       "ninfo"),
    ("time_stepping_nbq",  0,   0,  "int",   "&croco_time_stepping_nbq",   "ndtnbq"),
    ("time_stepping_nbq",  0,   1,  "float", "&croco_time_stepping_nbq",   "csound_nbq"),
    ("time_stepping_nbq",  0,   2,  "float", "&croco_time_stepping_nbq",   "visc2_nbq"),
    ("S-coord",            0,   0,  "float", "&croco_s_coord",              "theta_s"),
    ("S-coord",            0,   1,  "float", "&croco_s_coord",              "theta_b"),
    ("S-coord",            0,   2,  "float", "&croco_s_coord",              "Tcline"),
    ("start_date",         0,   0,  "str",   "&croco_use_calendar",        "start_date"),
    ("end_date",           0,   0,  "str",   "&croco_use_calendar",        "end_date"),
    ("output_time_steps",  0,   0,  "float", "&croco_use_calendar",        "dt_his"),
    ("output_time_steps",  0,   1,  "float", "&croco_use_calendar",        "dt_avg"),
    ("output_time_steps",  0,   2,  "float", "&croco_use_calendar",        "dt_rst"),
    ("history",            0,   0,  "bool",  "&croco_history",             "ldefhis"),
    ("history",            0,   1,  "int",   "&croco_history",             "nwrt"),
    ("history",            0,   2,  "int",   "&croco_history",             "nrpfhis"),
    ("history",            1,   0,  "str",   "&croco_history",             "hisname"),
    ("initial",            0,   0,  "int",   "&croco_initial",             "nrrec"),
    ("initial",            1,   0,  "str",   "&croco_initial",             "ininame"),
    ("restart",            0,   0,  "int",   "&croco_restart",             "nrst"),
    ("restart",            0,   1,  "int",   "&croco_restart",             "nrpfrst"),
    ("restart",            1,   0,  "str",   "&croco_restart",             "rstname"),
    ("grid",               0,   0,  "str",   "&croco_grid",                "grdname"),
    ("forcing",            0,   0,  "str",   "&croco_forcing",             "frcname"),
    ("bulk_forcing",       0,   0,  "str",   "&croco_bulk_forcing",        "bulkname"),
    ("climatology",        0,   0,  "str",   "&croco_climatology",         "clmname"),
    ("wave_offline",       0,   0,  "str",   "&croco_wave_offline",        "wave_file"),
    ("biology",            0,   0,  "str",   "&croco_biology",             "bioname"),
    ("boundary",           0,   0,  "str",   "&croco_boundary",            "bry_file"),
    ("wkb_boundary",       0,   0,  "str",   "&croco_wkb_boundary",        "brywkb_file"),
    ("bodyforce",          0,   0,  "int",   "&croco_bodyforce",           "levsfrc"),
    ("bodyforce",          0,   1,  "int",   "&croco_bodyforce",           "levbfrc"),
    ("lin_EOS_cff",        0,   0,  "float", "&croco_lin_eos",             "R0"),
    ("lin_EOS_cff",        0,   1,  "float", "&croco_lin_eos",             "T0"),
    ("lin_EOS_cff",        0,   2,  "float", "&croco_lin_eos",             "S0"),
    ("lin_EOS_cff",        0,   3,  "float", "&croco_lin_eos",             "Tcoef"),
    ("lin_EOS_cff",        0,   4,  "float", "&croco_lin_eos",             "Scoef"),
    ("abl_nudg_tra_cof",   0,   0,  "float", "&croco_abl_nudg_tra",        "ltra_min"),
    ("abl_nudg_tra_cof",   0,   1,  "float", "&croco_abl_nudg_tra",        "ltra_max"),
    ("abl_nudg_dyn_cof",   0,   0,  "float", "&croco_abl_nudg_dyn",        "ldyn_min"),
    ("abl_nudg_dyn_cof",   0,   1,  "float", "&croco_abl_nudg_dyn",        "ldyn_max"),
    ("sediments",          0,   0,  "str",   "&croco_sediments",           "sedname"),
    ("sediments_mustang",  0,   0,  "str",   "&croco_sediments_mustang",   "sedname_must"),
    ("substance",          0,   0,  "str",   "&croco_substance",           "subsfilename"),
    ("obstruction",        0,   0,  "str",   "&croco_obstruction",         "obstname"),
    ("xios_origin_date",   0,   0,  "str_line", "&croco_xios_origin_date", "xios_origin_date"),
    ("assimilation",       0,   0,  "str",   "&croco_assimilation",        "aparnam"),
    ("assimilation",       1,   0,  "str",   "&croco_assimilation",        "assname"),
    ("rho0",               0,   0,  "float", "&croco_rho0",                "rho0"),
    ("bottom_drag",        0,   0,  "float", "&croco_bottom_drag",         "rdrg"),
    ("bottom_drag",        0,   1,  "float", "&croco_bottom_drag",         "rdrg2"),
    ("bottom_drag",        0,   2,  "float", "&croco_bottom_drag",         "Zobt"),
    ("bottom_drag",        0,   3,  "float", "&croco_bottom_drag",         "Cdb_min"),
    ("bottom_drag",        0,   4,  "float", "&croco_bottom_drag",         "Cdb_max"),
    ("gamma2",             0,   0,  "float", "&croco_gamma2",              "gamma2"),
    ("lateral_visc",       0,   0,  "float", "&croco_lateral_visc",        "visc2"),
    ("lateral_visc",       0,   1,  "float", "&croco_lateral_visc",        "visc4"),
    ("tracer_diff2",       0,   0,  "float_array", "&croco_tracer_diff2",   "tnu2"),
    ("tracer_diff4",       0,   0,  "float_array", "&croco_tracer_diff4",   "tnu4"),
    ("vertical_mixing",    0,   0,  "float",       "&croco_vertical_mixing", "Akv_bak"),
    ("vertical_mixing",    0,   1,  "float_array", "&croco_vertical_mixing", "Akt_bak"),
    ("wkb_wwave",          0,   0,  "float", "&croco_wkb_wwave",           "wkb_amp"),
    ("wkb_wwave",          0,   1,  "float", "&croco_wkb_wwave",           "wkb_ang"),
    ("wkb_wwave",          0,   2,  "float", "&croco_wkb_wwave",           "wkb_prd"),
    ("wkb_wwave",          0,   3,  "float", "&croco_wkb_wwave",           "wkb_tide"),
    ("wkb_wwave",          0,   4,  "float", "&croco_wkb_wwave",           "wkb_btg"),
    ("wkb_wwave",          0,   5,  "float", "&croco_wkb_wwave",           "wkb_gam"),
    ("wkb_roller",         0,   0,  "float", "&croco_wkb_roller",          "wkb_rsb"),
    ("wkb_roller",         0,   1,  "float", "&croco_wkb_roller",          "wkb_roller"),
    ("wave_maker",         0,   0,  "float", "&croco_wave_maker",          "wmaker_amp"),
    ("wave_maker",         0,   1,  "float", "&croco_wave_maker",          "wmaker_prd"),
    ("wave_maker",         0,   2,  "float", "&croco_wave_maker",          "wmaker_dir"),
    ("wave_maker",         0,   3,  "float", "&croco_wave_maker",          "wmaker_dsp"),
    ("wave_maker",         0,   4,  "float", "&croco_wave_maker",          "wmaker_fsp"),
    ("averages",           0,   0,  "int",   "&croco_averages",            "ntsavg"),
    ("averages",           0,   1,  "int",   "&croco_averages",            "navg"),
    ("averages",           0,   2,  "int",   "&croco_averages",            "nrpfavg"),
    ("averages",           1,   0,  "str",   "&croco_averages",            "avgname"),
    ("surf",               0,   0,  "bool",  "&croco_surf",                "ldefsurf"),
    ("surf",               0,   1,  "int",   "&croco_surf",                "nwrtsurf"),
    ("surf",               0,   2,  "int",   "&croco_surf",                "nrpfsurf"),
    ("surf",               1,   0,  "str",   "&croco_surf",                "surfname"),

    ("surf_avg",           0,   0,  "bool",  "&croco_surf_avg",            "ldefsurf_avg"),
    ("surf_avg",           0,   1,  "int",   "&croco_surf_avg",            "ntssurf_avg"),
    ("surf_avg",           0,   2,  "int",   "&croco_surf_avg",            "nwrtsurf_avg"),
    ("surf_avg",           0,   3,  "int",   "&croco_surf_avg",            "nrpfsurf_avg"),
    ("surf_avg",           1,   0,  "str",   "&croco_surf_avg",            "surfname_avg"),

    ("diagnostics",        0,   0,  "bool",  "&croco_diagnostics_ts",      "ldefdia"),
    ("diagnostics",        0,   1,  "int",   "&croco_diagnostics_ts",      "nwrtdia"),
    ("diagnostics",        0,   2,  "int",   "&croco_diagnostics_ts",      "nrpfdia"),
    ("diagnostics",        1,   0,  "str",   "&croco_diagnostics_ts",      "dianame"),
    ("diag_avg",           0,   0,  "bool",  "&croco_diag_avg",            "ldefdia_avg"),
    ("diag_avg",           0,   1,  "int",   "&croco_diag_avg",            "ntsdia_avg"),
    ("diag_avg",           0,   2,  "int",   "&croco_diag_avg",            "nwrtdia_avg"),
    ("diag_avg",           0,   3,  "int",   "&croco_diag_avg",            "nrpfdia_avg"),
    ("diag_avg",           1,   0,  "str",   "&croco_diag_avg",            "dianame_avg"),
    ("diag_mld_dens",      0,   0,  "float", "&croco_diag_mld_dens",       "mld_crit_D"),
    ("diag_mld_dens",      0,   1,  "float", "&croco_diag_mld_dens",       "mld_crit_T"),
    ("diagnosticsM",       0,   0,  "bool",  "&croco_diagnosticsM",        "ldefdiaM"),
    ("diagnosticsM",       0,   1,  "int",   "&croco_diagnosticsM",        "nwrtdiaM"),
    ("diagnosticsM",       0,   2,  "int",   "&croco_diagnosticsM",        "nrpfdiaM"),
    ("diagnosticsM",       1,   0,  "str",   "&croco_diagnosticsM",        "dianameM"),
    ("diagM_avg",          0,   0,  "bool",  "&croco_diagM_avg",           "ldefdiaM_avg"),
    ("diagM_avg",          0,   1,  "int",   "&croco_diagM_avg",           "ntsdiaM_avg"),
    ("diagM_avg",          0,   2,  "int",   "&croco_diagM_avg",           "nwrtdiaM_avg"),
    ("diagM_avg",          0,   3,  "int",   "&croco_diagM_avg",           "nrpfdiaM_avg"),
    ("diagM_avg",          1,   0,  "str",   "&croco_diagM_avg",           "dianameM_avg"),
    ("diags_vrt",          0,   0,  "bool",  "&croco_diags_vrt",           "ldefdiags_vrt"),
    ("diags_vrt",          0,   1,  "int",   "&croco_diags_vrt",           "nwrtdiags_vrt"),
    ("diags_vrt",          0,   2,  "int",   "&croco_diags_vrt",           "nrpfdiags_vrt"),
    ("diags_vrt",          1,   0,  "str",   "&croco_diags_vrt",           "diags_vrtname"),
    ("diags_vrt_avg",      0,   0,  "bool",  "&croco_diags_vrt_avg",       "ldefdiags_vrt_avg"),
    ("diags_vrt_avg",      0,   1,  "int",   "&croco_diags_vrt_avg",       "ntsdiags_vrt_avg"),
    ("diags_vrt_avg",      0,   2,  "int",   "&croco_diags_vrt_avg",       "nwrtdiags_vrt_avg"),
    ("diags_vrt_avg",      0,   3,  "int",   "&croco_diags_vrt_avg",       "nrpfdiags_vrt_avg"),
    ("diags_vrt_avg",      1,   0,  "str",   "&croco_diags_vrt_avg",       "diags_vrtname_avg"),
    ("diags_ek",           0,   0,  "bool",  "&croco_diags_ek",            "ldefdiags_ek"),
    ("diags_ek",           0,   1,  "int",   "&croco_diags_ek",            "nwrtdiags_ek"),
    ("diags_ek",           0,   2,  "int",   "&croco_diags_ek",            "nrpfdiags_ek"),
    ("diags_ek",           1,   0,  "str",   "&croco_diags_ek",            "diags_ekname"),
    ("diags_ek_avg",       0,   0,  "bool",  "&croco_diags_ek_avg",        "ldefdiags_ek_avg"),
    ("diags_ek_avg",       0,   1,  "int",   "&croco_diags_ek_avg",        "ntsdiags_ek_avg"),
    ("diags_ek_avg",       0,   2,  "int",   "&croco_diags_ek_avg",        "nwrtdiags_ek_avg"),
    ("diags_ek_avg",       0,   3,  "int",   "&croco_diags_ek_avg",        "nrpfdiags_ek_avg"),
    ("diags_ek_avg",       1,   0,  "str",   "&croco_diags_ek_avg",        "diags_ekname_avg"),
    ("diags_pv",           0,   0,  "bool",  "&croco_diags_pv",            "ldefdiags_pv"),
    ("diags_pv",           0,   1,  "int",   "&croco_diags_pv",            "nwrtdiags_pv"),
    ("diags_pv",           0,   2,  "int",   "&croco_diags_pv",            "nrpfdiags_pv"),
    ("diags_pv",           1,   0,  "str",   "&croco_diags_pv",            "diags_pvname"),
    ("diags_pv_avg",       0,   0,  "bool",  "&croco_diags_pv_avg",        "ldefdiags_pv_avg"),
    ("diags_pv_avg",       0,   1,  "int",   "&croco_diags_pv_avg",        "ntsdiags_pv_avg"),
    ("diags_pv_avg",       0,   2,  "int",   "&croco_diags_pv_avg",        "nwrtdiags_pv_avg"),
    ("diags_pv_avg",       0,   3,  "int",   "&croco_diags_pv_avg",        "nrpfdiags_pv_avg"),
    ("diags_pv_avg",       1,   0,  "str",   "&croco_diags_pv_avg",        "diags_pvname_avg"),
    ("diags_eddy",         0,   0,  "bool",  "&croco_diags_eddy",          "ldefdiags_eddy"),
    ("diags_eddy",         0,   1,  "int",   "&croco_diags_eddy",          "nwrtdiags_eddy"),
    ("diags_eddy",         0,   2,  "int",   "&croco_diags_eddy",          "nrpfdiags_eddy"),
    ("diags_eddy",         1,   0,  "str",   "&croco_diags_eddy",          "diags_eddyname"),
    ("diags_eddy_avg",     0,   0,  "bool",  "&croco_diags_eddy_avg",      "ldefdiags_eddy_avg"),
    ("diags_eddy_avg",     0,   1,  "int",   "&croco_diags_eddy_avg",      "ntsdiags_eddy_avg"),
    ("diags_eddy_avg",     0,   2,  "int",   "&croco_diags_eddy_avg",      "nwrtdiags_eddy_avg"),
    ("diags_eddy_avg",     0,   3,  "int",   "&croco_diags_eddy_avg",      "nrpfdiags_eddy_avg"),
    ("diags_eddy_avg",     1,   0,  "str",   "&croco_diags_eddy_avg",      "diags_eddyname_avg"),
    ("diagnostics_bio",    0,   0,  "bool",  "&croco_diagnostics_bio",     "ldefdiabio"),
    ("diagnostics_bio",    0,   1,  "int",   "&croco_diagnostics_bio",     "nwrtdiabio"),
    ("diagnostics_bio",    0,   2,  "int",   "&croco_diagnostics_bio",     "nrpfdiabio"),
    ("diagnostics_bio",    1,   0,  "str",   "&croco_diagnostics_bio",     "dianamebio"),
    ("diagbio_avg",        0,   0,  "bool",  "&croco_diagbio_avg",         "ldefdiabio_avg"),
    ("diagbio_avg",        0,   1,  "int",   "&croco_diagbio_avg",         "ntsdiabio_avg"),
    ("diagbio_avg",        0,   2,  "int",   "&croco_diagbio_avg",         "nwrtdiabio_avg"),
    ("diagbio_avg",        0,   3,  "int",   "&croco_diagbio_avg",         "nrpfdiabio_avg"),
    ("diagbio_avg",        1,   0,  "str",   "&croco_diagbio_avg",         "dianamebio_avg"),

    ("abl",                0,   0,  "bool",  "&croco_abl",                 "ldefablhis"),
    ("abl",                0,   1,  "int",   "&croco_abl",                 "nwrtablhis"),
    ("abl",                0,   2,  "int",   "&croco_abl",                 "nrpfablhis"),
    ("abl",                1,   0,  "str",   "&croco_abl",                 "ablname"),
    ("abl_averages",       0,   0,  "bool",  "&croco_abl_averages",        "ldefablavg"),
    ("abl_averages",       0,   1,  "int",   "&croco_abl_averages",        "ntsablavg"),
    ("abl_averages",       0,   2,  "int",   "&croco_abl_averages",        "nwrtablavg"),
    ("abl_averages",       0,   3,  "int",   "&croco_abl_averages",        "nrpfablavg"),
    ("abl_averages",       1,   0,  "str",   "&croco_abl_averages",        "ablname_avg"),
]


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


def preprocess_param(cppdefs_path, param_path):
    """
    Run cpp on param.h using cppdefs.h for defines.
    Returns the preprocessed text (no CPP directives, resolved conditionals).
    """
    cppdefs_abs = os.path.realpath(cppdefs_path)
    param_dir   = os.path.dirname(os.path.realpath(param_path))
    wrapper = f'#include "{cppdefs_abs}"\n#include "param.h"\n'
    return _run_cpp(wrapper, [param_dir])


def extract_case_section(cppdefs_path):
    """
    Extract the case-specific block from cppdefs_path verbatim, preserving
    the original order and nested #ifdef structure.

    Returns (case_key, text) where text is ready to write as cppdefs_<Name>.h:
        # define <CASE_KEY>           <- synthetic first line
        <block content, as-is>        <- inside #elif defined <CASE>
        <trailing #include lines>     <- e.g. #include "cppdefs_dev.h"
    """
    lines = open(cppdefs_path).readlines()

    # 1. Detect the active case key from the top-level #define inserted by
    #    patch_select_case (the line after "#undef REGIONAL")
    case_key = None
    for i, line in enumerate(lines):
        if re.match(r'#undef\s+REGIONAL\b', line.strip()):
            if i + 1 < len(lines):
                m = re.match(r'#define\s+(\w+)', lines[i + 1].strip())
                if m:
                    case_key = m.group(1)
            break
    if case_key is None:
        for line in lines:
            if re.match(r'#define\s+REGIONAL\b', line.strip()):
                case_key = "REGIONAL"
                break
    if case_key is None:
        return None, ""

    # 2. Find the start of the case block (#elif defined CASE or #if for REGIONAL)
    if case_key == "REGIONAL":
        start_re = re.compile(r'^#if\s+defined\s+REGIONAL\b')
    else:
        start_re = re.compile(rf'^#elif\s+defined\s+{re.escape(case_key)}\b')

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

    for line in lines[start_idx + 1:]:
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
            if re.match(r'^#\s*(elif|else)\b', s):
                break
            if re.match(r'^#\s*endif\b', s):
                break

        if re.match(r'^#\s*if(def|ndef)?\b', s):
            depth += 1
        elif re.match(r'^#\s*endif\b', s):
            depth -= 1

        block.append(line)

    # 4. Collect trailing #include lines that follow the main closing #endif.
    # In reverse: includes/blanks appear BEFORE the #endif, so collect them
    # until we hit the #endif, then stop.
    trailing = []
    for line in reversed(lines):
        s = line.strip()
        if re.match(r'^#\s*endif\b', s):
            break   # reached the closing #endif — done
        if s.startswith('#include') or s == '':
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
    (e.g. NT = itemp + ntrc_salt + ntrc_pas + …).
    """
    params = {}

    # Simple "integer, parameter :: NAME = VALUE" (hand-written format)
    for m in re.finditer(r'integer\s*,\s*parameter\s*::\s*(\w+)\s*=\s*(\d+)', text):
        params[m.group(1)] = int(m.group(2))

    # Fortran "parameter (NAME=EXPR, …)" statements
    for stmt in re.finditer(r'parameter\s*\(([^)]+)\)', text, re.IGNORECASE):
        for assign in stmt.group(1).split(','):
            m = re.match(r'\s*(\w+)\s*=\s*(.+)', assign.strip())
            if not m:
                continue
            name, expr = m.group(1), m.group(2).strip()
            try:
                params[name] = _eval_fortran_expr(expr, params)
            except Exception:
                pass

    # Second pass: resolve any that depended on later definitions
    changed = True
    while changed:
        changed = False
        for stmt in re.finditer(r'parameter\s*\(([^)]+)\)', text, re.IGNORECASE):
            for assign in stmt.group(1).split(','):
                m = re.match(r'\s*(\w+)\s*=\s*(.+)', assign.strip())
                if not m:
                    continue
                name, expr = m.group(1), m.group(2).strip()
                if name not in params:
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
    resolved = re.sub(r'[A-Za-z_]\w*', replace_name, expr)
    if not re.match(r'^[\d\s+\-*/()]+$', resolved):
        raise ValueError(f"Unsafe: {resolved}")
    return int(eval(resolved))  # noqa: S307 – input sanitised above


def read_param_file(path):
    return read_param_text(open(path).read())


# ---------------------------------------------------------------------------
# Output generators
# ---------------------------------------------------------------------------

def make_clean_param(preprocessed_text):
    """
    Strip blank lines and Fortran comment lines from preprocessed param text.
    Keeps all active Fortran statements.
    """
    lines = []
    for line in preprocessed_text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith('!') or s.upper().startswith('C ') or s.upper() == 'C':
            continue
        lines.append(line.rstrip())
    return "\n".join(lines) + "\n"

# ---------------------------------------------------------------------------
# Auto-discovery helpers
# ---------------------------------------------------------------------------

def derive_case_name(input_path):
    base = os.path.basename(input_path)
    m = re.match(r'^.+\.in\.([^.]+)', base)
    return m.group(1) if m else "case"


def derive_output_name(input_path):
    basename = os.path.basename(input_path)
    m = re.match(r'^(.+)\.in\.([^.]+)(\.\d+)?$', basename)
    if m:
        return f"{m.group(1)}_{m.group(2)}.nml{m.group(3) or ''}"
    m = re.match(r'^(.+)\.in$', basename)
    if m:
        return f"{m.group(1)}.nml"
    return basename + ".nml"


def find_cppdefs(input_path, case_name):
    """
    Auto-discover the cppdefs file for this case.
    Tries cppdefs_<Name>.h first (simplified), then cppdefs_full_<Name>.h.
    Returns (path, is_full) or (None, False).
    """
    here = os.path.dirname(os.path.abspath(input_path))
    simplified = os.path.join(here, f"cppdefs_{case_name}.h")
    full       = os.path.join(here, f"cppdefs_full_{case_name}.h")
    if os.path.isfile(simplified):
        return simplified, False
    if os.path.isfile(full):
        return full, True
    return None, False


def is_full_cppdefs(path):
    """Return True if the file contains the full #if/#elif case chain."""
    text = open(path).read()
    return bool(re.search(r'^#undef\s+REGIONAL\b', text, re.MULTILINE))


# ---------------------------------------------------------------------------
# .in file parser
# ---------------------------------------------------------------------------

def parse_in_file(path):
    cards = {}
    current_key = None
    current_lines = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if stripped == "":
                continue
            m = re.match(r'^([A-Za-z][A-Za-z0-9_\-]*):', line)
            if m:
                if current_key is not None:
                    cards[current_key] = current_lines
                current_key = m.group(1)
                current_lines = []
                continue
            if current_key is not None:
                current_lines.append(stripped)
    if current_key is not None:
        cards[current_key] = current_lines
    return cards


# ---------------------------------------------------------------------------
# Value formatter
# ---------------------------------------------------------------------------

def format_value(raw, typ):
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


def expand_repeat(tokens):
    expanded = []
    for t in tokens:
        if '*' in t:
            count_str, val_str = t.split('*', 1)
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


# ---------------------------------------------------------------------------
# NML builder
# ---------------------------------------------------------------------------

def build_nml(cards, mappings, params):
    NT = params.get("NT", None)
    nml_entries = {}
    nml_order   = []

    for (card_name, line_idx, pos_idx, typ, nml_name, nml_var) in mappings:
        if card_name not in cards:
            continue
        vlines = [l for l in cards[card_name] if l.strip()]
        if line_idx >= len(vlines):
            continue
        tokens = vlines[line_idx].split()

        if typ == "float_array":
            raw_tokens = tokens[pos_idx:]
            if not raw_tokens:
                continue
            val = format_float_array(expand_repeat(raw_tokens), NT)
            if not val:
                continue
            if nml_name not in nml_entries:
                nml_entries[nml_name] = []
                nml_order.append(nml_name)
            nml_entries[nml_name].append(f"  {nml_var} = {val}")
            continue

        if typ == "str_line":
            raw = vlines[line_idx].strip()
        else:
            if pos_idx >= len(tokens):
                continue
            raw = tokens[pos_idx]

        if (typ == "str"
                and re.match(r'^\d{4}[-/]\d{2}[-/]\d{2}$', raw)
                and pos_idx + 1 < len(tokens)):
            raw = raw.replace("/", "-") + " " + tokens[pos_idx + 1]

        val = format_value(raw, typ)
        if nml_name not in nml_entries:
            nml_entries[nml_name] = []
            nml_order.append(nml_name)
        nml_entries[nml_name].append(f"  {nml_var} = {val}")

    lines = []
    for nml_name in nml_order:
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
    here = os.path.dirname(os.path.abspath(__file__))
    default_param = os.path.join(here, "..", "OCEAN", "param.h")

    parser = argparse.ArgumentParser(
        description="Convert a CROCO .in file to a Fortran namelist .nml file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Defaults (all auto-derived from croco.in.<Name>):
  cppdefs  : cppdefs_<Name>.h  (simplified)  or  cppdefs_full_<Name>.h  (full)
  param    : ../OCEAN/param.h
  output   : croco_<Name>.nml
  param out: param_<Name>.h
  cpp out  : cppdefs_<Name>.h  (only written when input was a full cppdefs)
""")
    parser.add_argument("input",
                        help="Path to the .in file (e.g. croco.in.Estuary)")
    parser.add_argument("-o", "--output",
                        help="Output .nml path")
    parser.add_argument("-c", "--cppdefs",
                        help="cppdefs file (simplified or full); auto-discovered if omitted")
    parser.add_argument("-p", "--param",
                        help=f"param.h path (default: {default_param})")
    parser.add_argument("--out-param",
                        help="Output param_<Name>.h path")
    parser.add_argument("--out-cppdefs",
                        help="Output cppdefs_<Name>.h path")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file '{args.input}' not found.")

    case_name   = derive_case_name(args.input)
    input_dir   = os.path.dirname(os.path.abspath(args.input))
    out_nml     = args.output     or os.path.join(input_dir, derive_output_name(args.input))
    out_dir     = os.path.dirname(os.path.abspath(out_nml))

    # ---- Resolve cppdefs file ------------------------------------------------
    if args.cppdefs:
        cppdefs_path = args.cppdefs
        if not os.path.isfile(cppdefs_path):
            sys.exit(f"Error: cppdefs file '{cppdefs_path}' not found.")
        full_cppdefs = is_full_cppdefs(cppdefs_path)
    else:
        cppdefs_path, full_cppdefs = find_cppdefs(args.input, case_name)

    # ---- Resolve param.h -----------------------------------------------------
    param_path = args.param or (default_param if os.path.isfile(default_param) else None)

    # ---- Log inputs ----------------------------------------------------------
    cpp_mode = "full" if full_cppdefs else "simplified"
    print(f"Input .in      : {args.input}")
    print(f"Input cppdefs  : {cppdefs_path or '(not found)'}  [{cpp_mode}]")
    print(f"Input param.h  : {param_path or '(not found)'}")
    print()

    # ---- Preprocess param.h with cpp -----------------------------------------
    params = {}
    if cppdefs_path and param_path and os.path.isfile(param_path):
        preprocessed = preprocess_param(cppdefs_path, param_path)
        params = read_param_text(preprocessed)

        out_param = args.out_param or os.path.join(out_dir, f"param_{case_name}.h")
        with open(out_param, "w") as f:
            f.write(make_clean_param(preprocessed))
        print(f"Output param   : {out_param}  (NT={params.get('NT', '?')})")

        # Extract and write clean cppdefs only when input was a full cppdefs
        if full_cppdefs:
            out_cppdefs = args.out_cppdefs or os.path.join(out_dir, f"cppdefs_{case_name}.h")
            _ck, section_text = extract_case_section(cppdefs_path)
            if section_text:
                with open(out_cppdefs, "w") as f:
                    f.write(section_text)
                print(f"Output cppdefs : {out_cppdefs}  (case={_ck})")
        else:
            print(f"Output cppdefs : (input already simplified — not rewritten)")
    else:
        if not cppdefs_path:
            print(f"Warning        : no cppdefs found for '{case_name}' — NT unknown.",
                  file=sys.stderr)
        elif not param_path or not os.path.isfile(param_path):
            print(f"Warning        : param.h not found at '{param_path}' — NT unknown.",
                  file=sys.stderr)

    # ---- Convert .in → .nml --------------------------------------------------
    cards = parse_in_file(args.input)

    mapped_cards  = {row[0] for row in MAPPINGS}
    present_cards = set(cards.keys())
    missing   = mapped_cards - present_cards
    unmapped  = present_cards - mapped_cards
    if missing:
        print(f"Skipped (no .in card)  : {', '.join(sorted(missing))}")
    if unmapped:
        print(f"Ignored (no mapping)   : {', '.join(sorted(unmapped))}")

    nml_text = build_nml(cards, MAPPINGS, params)
    with open(out_nml, "w") as f:
        f.write(nml_text)
    print(f"Output .nml    : {out_nml}")


if __name__ == "__main__":
    main()
