#!/usr/bin/env python3
#  CROCO — Coastal and Regional Ocean COmmunity model
#  Copyright (C) 2005-2026 CROCO Development Team
#
#  SPDX-License-Identifier: CECILL-2.1

"""
Plot results from the TS_HADV_TEST test case.

Displays the surface temperature field at initial, mid, and final time,
plus the advection error T(final) - T(initial).  Prints quantitative
diagnostics: peak erosion, L2 error norm, and salinity constancy
preservation.

Works for all three variants: SOLID_BODY_ROT, SOLID_BODY_PER,
DIAGONAL_ADV.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse

# ── CLI ──────────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Plot results from the TS_HADV_TEST test case.",
    epilog="Example usage:\n"
    + "  python plot_ts_hadv_test.py --file croco_his_body_rot.nc\n"
    + "  python plot_ts_hadv_test.py --makepng --no-show --output-dir ./plots",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="croco_his_body_rot.nc",
    help="Path to the NetCDF history file",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate a PDF of the plots",
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate a PNG of the plots",
)
parser.add_argument(
    "--no-show", action="store_true",
    help="Do not display the plots on the screen",
)
parser.add_argument(
    "--output-dir", type=str, default=".",
    help="Directory to save the output files (default: current directory)",
)
args = parser.parse_args()

# ── Read data ────────────────────────────────────────────

try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Grid (uniform analytical grid, coordinates in meters)
x = nc.variables["x_rho"][:] / 1000.0  # km
y = nc.variables["y_rho"][:] / 1000.0  # km

# Time
time = nc.variables["scrum_time"][:]
nt = len(time)

# Surface level = top sigma (last index)
klev = -1

# Temperature at 3 instants: initial, quarter-period, final
# T/4 is the time of maximum displacement for all variants:
#   SOLID_BODY_PER: cos-modulated, integral peaks at T/4
#   SOLID_BODY_ROT: quarter rotation
#   DIAGONAL_ADV:   quarter diagonal traverse
t0 = np.squeeze(nc.variables["temp"][0, klev, :, :])
tquart = np.squeeze(nc.variables["temp"][nt // 4, klev, :, :])
tfin = np.squeeze(nc.variables["temp"][-1, klev, :, :])

# Salinity constancy check (if available)
has_salt = "salt" in nc.variables
if has_salt:
    s0 = np.squeeze(nc.variables["salt"][0, klev, :, :])
    sfin = np.squeeze(nc.variables["salt"][-1, klev, :, :])

# Title
title = nc.title.strip() if hasattr(nc, "title") else "TS_HADV_TEST"

nc.close()

# ── Diagnostics ──────────────────────────────────────────

error = tfin - t0
peak_init = np.max(t0)
peak_final = np.max(tfin)
peak_ratio = peak_final / peak_init if peak_init > 0 else 0.0
l2_error = np.sqrt(np.mean(error ** 2))
linf_error = np.max(np.abs(error))

print(f"=== {title} ===")
print(f"  Peak T initial : {peak_init:.6f}")
print(f"  Peak T final   : {peak_final:.6f}")
print(f"  Peak ratio     : {peak_ratio:.6f}  (1.0 = perfect)")
print(f"  L2 error       : {l2_error:.6e}")
print(f"  Linf error     : {linf_error:.6e}")

if has_salt:
    salt_err = np.max(np.abs(sfin - 1.0))
    print(f"  Salt constancy : max|S-1| = {salt_err:.6e}")

# ── Plot ─────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# Common color range for T panels (based on initial field)
vmin_t = np.min(t0)
vmax_t = np.max(t0)

# Panel 1: T initial
ax = axes[0, 0]
cf = ax.pcolormesh(x, y, t0, cmap="RdYlBu_r", shading="auto",
                   vmin=vmin_t, vmax=vmax_t)
fig.colorbar(cf, ax=ax)
ax.set_ylabel("Y (km)")
ax.set_title("Temperature — t = 0")
ax.set_aspect("equal")
ax.tick_params(labelbottom=False)

# Panel 2: T quarter-period
ax = axes[0, 1]
cf = ax.pcolormesh(x, y, tquart, cmap="RdYlBu_r", shading="auto",
                   vmin=vmin_t, vmax=vmax_t)
fig.colorbar(cf, ax=ax)
ax.set_title("Temperature — t = T/4")
ax.set_aspect("equal")
ax.tick_params(labelbottom=False)

# Panel 3: T final
ax = axes[1, 0]
cf = ax.pcolormesh(x, y, tfin, cmap="RdYlBu_r", shading="auto",
                   vmin=vmin_t, vmax=vmax_t)
fig.colorbar(cf, ax=ax)
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_title("Temperature — t = T (return)")
ax.set_aspect("equal")

# Panel 4: Error T(final) - T(initial)
ax = axes[1, 1]
emax = max(np.max(np.abs(error)), 1e-15)
cf = ax.pcolormesh(x, y, error, cmap="RdBu_r", shading="auto",
                   vmin=-emax, vmax=emax)
fig.colorbar(cf, ax=ax)
ax.set_xlabel("X (km)")
ax.set_title("Error T(final) − T(initial)")
ax.set_aspect("equal")

fig.suptitle(f"{title}\n"
             f"Peak ratio = {peak_ratio:.4f}    "
             f"L2 error = {l2_error:.2e}    "
             f"L∞ error = {linf_error:.2e}",
             fontsize=12, fontweight="bold")
plt.subplots_adjust(hspace=0.3, top=0.90)

# ── Save / Show ──────────────────────────────────────────

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

if args.makepdf:
    pdf_path = os.path.join(output_dir, "ts_hadv_test_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "ts_hadv_test_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()
