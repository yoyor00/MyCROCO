#!/usr/bin/env python3

"""
Plot results from the FLASH_RIP test case.
Displays sea level and wave-mean vorticity from a rip current simulation.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import croco_utils as cr
import argparse

# ── CLI (standard BENCH-compatible interface) ────────────

parser = argparse.ArgumentParser(
    description="Plot results from the FLASH_RIP test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="rip_his.nc",
    help="Path to the NetCDF history file",
)
parser.add_argument(
    "--makepdf",
    action="store_true",
    help="Generate a PDF of the plots",
)
parser.add_argument(
    "--makepng",
    action="store_true",
    help="Generate a PNG of the plots",
)
parser.add_argument(
    "--no-show",
    action="store_true",
    help="Do not display the plots on the screen",
)
parser.add_argument(
    "--output-dir",
    type=str,
    default=".",
    help="Directory to save the output files",
)
parser.add_argument(
    "--tindex",
    type=int,
    default=None,
    help="Time index to plot (default: last time step)",
)
args = parser.parse_args()

# Derive average file path from history file path
hisname = args.file
avgname = args.file.replace("_his.nc", "_avg.nc")

# ── Read data ────────────────────────────────────────────

try:
    nc = Dataset(hisname, "r")
except FileNotFoundError:
    print(f"Error: File '{hisname}' not found.")
    exit(1)

time_len = len(nc.variables["scrum_time"][:])
tindex = (time_len - 1) if args.tindex is None else min(args.tindex, time_len - 1)

time = nc.variables["scrum_time"][tindex] / 60.0  # minutes

h = nc.variables["h"][:]
xl = nc.variables["xl"][:]
x = nc.variables["x_rho"][:] - xl
y = nc.variables["y_rho"][:]
pm = nc.variables["pm"][:]
pn = nc.variables["pn"][:]
N = len(nc.dimensions["s_rho"])
Dcrit = nc.variables["Dcrit"][:]
zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])

nc.close()

try:
    nc = Dataset(avgname, "r")
except FileNotFoundError:
    print(f"Error: Average file '{avgname}' not found.")
    exit(1)

avg_tindex = max(0, tindex - 1)
u = np.squeeze(nc.variables["u"][avg_tindex, N - 1, :, :])
v = np.squeeze(nc.variables["v"][avg_tindex, N - 1, :, :])

nc.close()

# Vorticity (psi → rho grid)
vort_psi = cr.vorticity(u, v, pm, pn)
Mp, Lp = pm.shape
M_psi, L_psi = vort_psi.shape
vort = np.zeros((Mp, Lp))
vort[1 : M_psi + 1, 1 : L_psi + 1] = vort_psi
vort[0, :] = vort[1, :]
vort[:, 0] = vort[:, 1]

# Sea level masking
zeta[h < Dcrit] = zeta[h < Dcrit] - Dcrit
mask = np.ones_like(zeta)
mask[np.isnan(zeta)] = np.nan

# ── Plot ─────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 7))

# Panel 1: Sea Level
ax = axes[0]
cmin, cmax = -0.7, 0.7
nbcol = 40
levels = np.linspace(cmin, cmax, nbcol + 1)
zeta_plot = np.clip(zeta, cmin, cmax)
zeta_plot[zeta_plot == 0.0] = np.nan
cs1 = ax.contourf(x, y, zeta_plot, levels=levels, cmap="jet")
cb1 = fig.colorbar(cs1, ax=ax)
cb1.set_ticks(np.arange(cmin, cmax + 0.2, 0.2))
cb1.set_ticklabels([f"{v:.1f}" for v in np.arange(cmin, cmax + 0.2, 0.2)])
ax.set_xlim(-250, -10)
ax.set_ylim(0, 300)
ax.set_title("Sea Level")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_aspect("equal", adjustable="box")
ax.tick_params(labelsize=15)

# Panel 2: Wave-mean Vorticity
ax = axes[1]
cmin, cmax = -0.07, 0.07
nbcol = 20
levels = np.linspace(cmin, cmax, nbcol + 1)
vort_plot = np.clip(vort, cmin, cmax) * mask
cs2 = ax.contourf(x, y, vort_plot, levels=levels, cmap="jet")
cb2 = fig.colorbar(cs2, ax=ax)
cb2.set_ticks(np.arange(cmin, cmax + 0.02, 0.02))
cb2.set_ticklabels([f"{v:.3f}" for v in np.arange(cmin, cmax + 0.02, 0.02)])
ax.set_xlim(-250, -10)
ax.set_ylim(0, 300)
ax.set_title("Wave-mean Vorticity")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_aspect("equal", adjustable="box")
ax.tick_params(labelsize=15)

plt.suptitle(
    f"FLASH_RIP test case — t = {time:.1f} min", fontsize=13, fontweight="bold"
)
plt.subplots_adjust(top=0.88)
plt.tight_layout()

# ── Save / Show (do not modify) ──────────────────────────

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

if args.makepdf:
    pdf_path = os.path.join(output_dir, "flashrip_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "flashrip_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()
