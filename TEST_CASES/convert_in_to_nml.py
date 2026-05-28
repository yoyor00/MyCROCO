#!/usr/bin/env python3
"""
croco_in_to_nml.py
==================
Convert a CROCO `.in` file into a Fortran namelist `.nml` file.

Usage
-----
    python croco_in_to_nml.py croco.in.Name            # -> croco_Name.nml
    python croco_in_to_nml.py croco.in                 # -> croco.nml
    python croco_in_to_nml.py croco.in -o my_out.nml   # explicit output

Extending the mapping
---------------------
Add a row to MAPPINGS:

    # card_name       line  pos  type      nml_name                  nml_var
    ("time_stepping",   0,   0,  "int",   "&croco_time_stepping",   "ntimes"),

  card_name : key that starts the card in the .in file  (before the colon)
  line      : 0-based index of the value line inside that card block
  pos       : 0-based token index on that line
  type      : "int" | "float" | "bool" | "str" | "str_line"
                int      -> written as a bare integer
                float    -> kept as-is (preserves Fortran d0 notation)
                bool     -> T/F -> .true./.false.
                str      -> double-quoted single token; if the value looks like
                            a date (YYYY-MM-DD) the time token is appended
                str_line -> double-quoted entire line (for multi-word values
                            like titles; pos is ignored)
  nml_name  : namelist block name (with leading &)
  nml_var   : variable name inside that block

Multiple rows sharing the same nml_name are merged into one block.
Cards absent from the .in file are silently skipped (no block written).
"""

import re
import sys
import os
import argparse

# ---------------------------------------------------------------------------
# Mapping table  –  add / remove rows freely
# ---------------------------------------------------------------------------
#  card_name            line  pos  type      nml_name                      nml_var
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

    # All calendar/date cards feed the same namelist block
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


    # ("xios_origin_date",   0,   0,  "str",   "&croco_use_calendar",        "xios_origin_date"),
    # ("boundary",           0,   0,  "str",   "&croco_boundary",            "fname_boundary"),



    # ("averages",           0,   0,  "int",   "&croco_averages",            "ntsavg"),
    # ("averages",           0,   1,  "int",   "&croco_averages",            "navg"),
    # ("averages",           0,   2,  "int",   "&croco_averages",            "nrpfavg"),
    # ("averages",           1,   0,  "str",   "&croco_averages",            "fname_averages"),

    # ("rho0",               0,   0,  "float", "&croco_rho0",                "rho0"),

    # ("lateral_visc",       0,   0,  "float", "&croco_lateral_visc",        "visc2"),
    # ("lateral_visc",       0,   1,  "float", "&croco_lateral_visc",        "visc4"),

    # ("tracer_diff2",       0,   0,  "float", "&croco_tracer_diff2",        "tnu2"),
    # ("tracer_diff4",       0,   0,  "float", "&croco_tracer_diff4",        "tnu4"),

    # ("vertical_mixing",    0,   0,  "float", "&croco_vertical_mixing",     "akv_bak"),
    # ("vertical_mixing",    0,   1,  "float", "&croco_vertical_mixing",     "akt_bak"),

    # ("bottom_drag",        0,   0,  "float", "&croco_bottom_drag",         "rdrg"),
    # ("bottom_drag",        0,   1,  "float", "&croco_bottom_drag",         "rdrg2"),
    # ("bottom_drag",        0,   2,  "float", "&croco_bottom_drag",         "zob"),
    # ("bottom_drag",        0,   3,  "float", "&croco_bottom_drag",         "cdb_min"),
    # ("bottom_drag",        0,   4,  "float", "&croco_bottom_drag",         "cdb_max"),

    # ("gamma2",             0,   0,  "float", "&croco_gamma2",              "gamma2"),

    # ("sponge",             0,   0,  "float", "&croco_sponge",              "x_sponge"),
    # ("sponge",             0,   1,  "float", "&croco_sponge",              "v_sponge"),

    # ("nudg_cof",           0,   0,  "float", "&croco_nudg_cof",            "taut_in"),
    # ("nudg_cof",           0,   1,  "float", "&croco_nudg_cof",            "taut_out"),
    # ("nudg_cof",           0,   2,  "float", "&croco_nudg_cof",            "taum_in"),
    # ("nudg_cof",           0,   3,  "float", "&croco_nudg_cof",            "taum_out"),

    # ("stations",           0,   0,  "bool",  "&croco_stations",            "ldefsta"),
    # ("stations",           0,   1,  "int",   "&croco_stations",            "nsta"),
    # ("stations",           0,   2,  "int",   "&croco_stations",            "nrpfsta"),
    # ("stations",           1,   0,  "str",   "&croco_stations",            "fname_stations_in"),
    # ("stations",           2,   0,  "str",   "&croco_stations",            "fname_stations"),
]


# ---------------------------------------------------------------------------
# .in file parser
# ---------------------------------------------------------------------------

def parse_in_file(path):
    """
    Return a dict  card_name -> [value_line_0, value_line_1, ...]
    The header line (the one with the colon) is excluded.
    A new CARDNAME: line always starts a new card (blank lines are ignored).
    """
    cards = {}
    current_key = None
    current_lines = []

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            stripped = line.strip()

            # Skip blank lines — they do NOT delimit cards
            if stripped == "":
                continue

            m = re.match(r'^([A-Za-z][A-Za-z0-9_\-]*):', line)
            if m:
                # Flush previous card
                if current_key is not None:
                    cards[current_key] = current_lines
                current_key = m.group(1)
                current_lines = []
                continue

            if current_key is not None:
                current_lines.append(stripped)

    if current_key is not None:          # flush last card
        cards[current_key] = current_lines

    return cards


# ---------------------------------------------------------------------------
# Value formatter
# ---------------------------------------------------------------------------

def format_value(raw, typ):
    """
    Format a single token according to its declared type.

      int   -> bare integer string  (strip trailing .0 if present)
      float -> kept as-is           (preserves Fortran d0 / scientific notation)
      bool  -> .true. / .false.
      str   -> double-quoted; if the value looks like a date (YYYY-MM-DD) the
               next token on the same line (the time part) is appended
    """
    if typ == "bool":
        if raw.upper() in ("T", ".TRUE."):
            return ".true."
        if raw.upper() in ("F", ".FALSE."):
            return ".false."
        return raw  # pass through if already .true./.false.

    if typ == "int":
        # Strip Fortran suffixes just in case (e.g. 720d0 -> 720)
        return str(int(float(raw.lower().replace("d", "e"))))

    if typ == "float":
        return raw   # keep Fortran notation verbatim

    if typ == "str_line":
        return f'"{raw}"'   # raw is already the full line, passed by build_nml

    if typ == "str":
        return f'"{raw}"'

    return raw   # unknown type: pass through


# ---------------------------------------------------------------------------
# NML builder
# ---------------------------------------------------------------------------

def build_nml(cards, mappings):
    """
    Walk MAPPINGS in order.  For each row whose card_name exists in `cards`,
    extract the requested token and format it, then accumulate into the
    appropriate namelist block.

    nml_name blocks are written in first-seen order.
    Missing cards produce no output at all.
    """
    nml_entries = {}   # nml_name -> list of "  var = val" strings
    nml_order   = []   # preserves first-seen insertion order

    for (card_name, line_idx, pos_idx, typ, nml_name, nml_var) in mappings:

        if card_name not in cards:
            continue   # card absent from .in -> skip silently

        # Filter to non-empty value lines
        vlines = [l for l in cards[card_name] if l.strip()]

        if line_idx >= len(vlines):
            continue   # line doesn't exist in this card

        tokens = vlines[line_idx].split()

        if typ == "str_line":
            raw = vlines[line_idx].strip()
        else:
            if pos_idx >= len(tokens):
                continue   # token doesn't exist on this line
            raw = tokens[pos_idx]

        # Assemble date + heure before formatting
        if (
            typ == "str"
            and re.match(r'^\d{4}[-/]\d{2}[-/]\d{2}$', raw)
            and pos_idx + 1 < len(tokens)
        ):
            raw = raw.replace("/", "-") + " " + tokens[pos_idx + 1]

        val = format_value(raw, typ)

        if nml_name not in nml_entries:
            nml_entries[nml_name] = []
            nml_order.append(nml_name)

        nml_entries[nml_name].append(f"  {nml_var} = {val}")

    # Render to text
    lines = []
    for nml_name in nml_order:
        lines.append(nml_name)
        lines.extend(nml_entries[nml_name])
        lines.append("/")
        lines.append("")

    lines.append("")
    lines.append("! End of namelist")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Output filename helper
# ---------------------------------------------------------------------------

def derive_output_name(input_path):
    """
    Examples
    --------
    croco.in.Name      -> croco_Name.nml
    croco.in.Name.1    -> croco_Name.nml.1
    croco.in           -> croco.nml
    """
    basename = os.path.basename(input_path)

    # Match: prefix.in.suffix[.number]
    m = re.match(r'^(.+)\.in\.([^.]+)(\.\d+)?$', basename)
    if m:
        prefix = m.group(1)
        suffix = m.group(2)
        trailing = m.group(3) or ""
        return f"{prefix}_{suffix}.nml{trailing}"

    m = re.match(r'^(.+)\.in$', basename)
    if m:
        return f"{m.group(1)}.nml"

    return basename + ".nml"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert a CROCO .in file to a Fortran namelist .nml file."
    )
    parser.add_argument("input", help="Path to the .in file (e.g. croco.in, croco.in.Name)")
    parser.add_argument("-o", "--output",
                        help="Output .nml path (default: derived from input name)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    out_path = args.output or derive_output_name(args.input)

    cards = parse_in_file(args.input)

    mapped_cards  = {row[0] for row in MAPPINGS}
    present_cards = set(cards.keys())
    missing       = mapped_cards - present_cards
    unmapped      = present_cards - mapped_cards

    print(f"Reading  : {args.input}")
    if missing:
        print(f"Skipped  (mapped but absent from .in) : {', '.join(sorted(missing))}")
    if unmapped:
        print(f"Ignored  (present in .in, no mapping) : {', '.join(sorted(unmapped))}")

    nml_text = build_nml(cards, MAPPINGS)

    with open(out_path, "w") as fh:
        fh.write(nml_text)

    print(f"Written  : {out_path}")


if __name__ == "__main__":
    main()
