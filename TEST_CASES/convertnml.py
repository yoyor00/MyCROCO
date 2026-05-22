#!/usr/bin/env python3

import re
import argparse
from pathlib import Path


def parse_time_stepping(text):
    # find the line with values after header
    match = re.search(
        r"time_stepping:.*?\n\s*([0-9.eEdD+-]+)\s+([0-9.eEdD+-]+)\s+([0-9.eEdD+-]+)\s+([0-9.eEdD+-]+)",
        text,
        re.IGNORECASE | re.DOTALL,
    )

    if not match:
        return 0, 0.0, 20, 1

    ntimes = int(match.group(1))
    dt = float(match.group(2))
    ndtfast = int(match.group(3))
    ninfo = int(match.group(4))

    return ntimes, dt, ndtfast, ninfo


def parse_title(text):
    match = re.search(r"title:\s*\n\s*(.+)", text, re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return "CROCO simulation"


def generate_nml_from_in(infilename, outfile=None):
    infile = Path(infilename)

    if outfile is None:
        print(infile.name[9:])
        outfile = Path("croco_" + infile.name[9:] + ".nml")

    text = infile.read_text()

    title = parse_title(text)
    ntimes, dt, ndtfast, ninfo = parse_time_stepping(text)

    # --- Write .nml ---
    with open(outfile, "w") as f:
        f.write("&croco_title\n")
        f.write(f'  title = "{title}"\n')
        f.write("/\n\n")

        f.write("&croco_time_stepping\n")
        f.write(f"  ntimes = {ntimes}\n")
        f.write(f"  dt = {dt:.3f}\n")
        f.write(f"  ndtfast= {ndtfast}\n")
        f.write(f"  ninfo = {ninfo}\n")
        f.write("/\n")

    print(f"Created: {outfile}")


def batch_convert(folder):
    folder = Path(folder)

    files = list(folder.glob("*.in.*"))
    if not files:
        print(f"No matching files found in {folder}")
        return

    for infile in files:
        generate_nml_from_in(infile)


def main():
    parser = argparse.ArgumentParser(
        description="Convert CROCO .in.* files to .nml format"
    )
    parser.add_argument("folder", help="Folder containing .in.* files")

    args = parser.parse_args()
    batch_convert(args.folder)


if __name__ == "__main__":
    main()
