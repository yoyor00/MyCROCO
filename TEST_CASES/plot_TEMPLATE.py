#!/usr/bin/env python3

"""
Plot results from the XXXX test case.
This is a skeleton for guidance zhen introducing a nez test case

TODO: describe what this test case does and what the plot shows.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# ── CLI (standard BENCH-compatible interface) ────────────

parser = argparse.ArgumentParser(
    description="Plot results from the XXXX test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="xxxx_his.nc",
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
    help="Directory to save the output files",
)
args = parser.parse_args()

# ── Read data ────────────────────────────────────────────

try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Time
time = nc.variables["scrum_time"][:] / 86400.0  # in days
nt = len(time)
tndx = -1  # last time step

# Grid
h = nc.variables["h"][:]
x = nc.variables["x_rho"][:] / 1000.0  # km
y = nc.variables["y_rho"][:] / 1000.0  # km

# Vertical grid (if needed for sections)
theta_s = float(nc.theta_s)
theta_b = float(nc.theta_b)
N = len(nc.dimensions["s_rho"])
vtrans_var = nc.variables.get("Vtransform", None)
vtransform = int(np.squeeze(vtrans_var[:])) if vtrans_var is not None else 2
hc_var = nc.variables.get("hc", None)
hc = float(hc_var[:]) if hc_var is not None else float(nc.hc)

# Read variables
# TODO: adapt to your test case
zeta = nc.variables["zeta"][tndx, :, :]
temp = nc.variables["temp"][tndx, -1, :, :]   # surface T

# Optional: vertical section at mid-domain
# j_sec = h.shape[0] // 2
# zeta_sec = np.zeros_like(h)
# zr = cr.zlevs(h, zeta_sec, theta_s, theta_b, hc, N, "r", vtransform)
# temp_sec = nc.variables["temp"][tndx, :, j_sec, :]
# z_sec = zr[:, j_sec, :]
# x_sec = x[j_sec, :]

nc.close()

# ── Plot ─────────────────────────────────────────────────

# TODO: adapt layout and content
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel 1: surface field
ax = axes[0]
cf = ax.pcolormesh(x, y, temp, cmap="RdYlBu_r", shading="auto")
fig.colorbar(cf, ax=ax, label="°C")
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_title(f"SST at t = {time[tndx]:.2f} days")
ax.set_aspect("equal")

# Panel 2: second field (example: zeta)
ax = axes[1]
cf = ax.pcolormesh(x, y, 100 * zeta, cmap="RdBu_r", shading="auto")
fig.colorbar(cf, ax=ax, label="cm")
ax.set_xlabel("X (km)")
ax.set_title(f"SSH (cm) at t = {time[tndx]:.2f} days")
ax.set_aspect("equal")

# Optional: vertical section
# ax = axes[1]
# cf = ax.pcolormesh(np.tile(x_sec, (N, 1)), z_sec, temp_sec,
#                    cmap="RdYlBu_r", shading="auto")
# fig.colorbar(cf, ax=ax, label="°C")
# ax.set_xlabel("X (km)")
# ax.set_ylabel("Depth (m)")
# ax.set_title(f"Temperature section at j={j_sec}")

plt.suptitle(f"XXXX test case — t = {time[tndx]:.2f} days",
             fontsize=13, fontweight="bold")
plt.subplots_adjust(top=0.88)

# ── Save / Show (do not modify) ──────────────────────────

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

if args.makepdf:
    pdf_path = os.path.join(output_dir, "xxxx_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "xxxx_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()
