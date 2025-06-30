#!/usr/bin/env python3
"""
Make plot from the results of the FLASH_RIP test case

Further Information:
http://www.crocoagrif.org

This file is part of CROCOTOOLS

CROCOTOOLS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

CROCOTOOLS is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston,
MA  02111-1307  USA

Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
     Application of the ROMS embedding procedure for the Central
     California Upwelling System,  Ocean Modelling, 2006.

Patrick Marchesiello, IRD 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import croco_utils as cr
import argparse
import os

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the FLASH_RIP test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--hisfile", type=str, default="rip_his.nc", help="Path to the history NetCDF file"
)
parser.add_argument(
    "--avgfile", type=str, default="rip_avg.nc", help="Path to the average NetCDF file"
)
parser.add_argument(
    "--tindex",
    type=int,
    default=None,
    help="Time index to plot (default: last time step)",
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as PDF files")
parser.add_argument("--makepng", action="store_true", help="Save plots as PNG files")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

print("=== FLASH_RIP Test Case ===")
print("Processing flash rip current simulation results")
if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
hisname = args.hisfile
avgname = args.avgfile

print(f"Reading history file: {hisname}")
print(f"Reading average file: {avgname}")

# Get model grid and variables from history file
print("Reading grid and sea level data...")
try:
    with Dataset(hisname, "r") as nc:
        # Get time index (last record by default)
        time_len = len(nc.variables["scrum_time"][:])
        if args.tindex is None:
            tindex = time_len - 1  # Python 0-based indexing
        else:
            tindex = min(args.tindex, time_len - 1)

        print(f"Available time steps: {time_len}, using index: {tindex}")

        time = nc.variables["scrum_time"][tindex] / 60  # Convert to minutes
        print(f"Time: {time:.1f} minutes")

        # Grid variables
        h = nc.variables["h"][:]
        xl = nc.variables["xl"][:]
        x = nc.variables["x_rho"][:] - xl
        y = nc.variables["y_rho"][:]
        pm = nc.variables["pm"][:]
        pn = nc.variables["pn"][:]
        N = len(nc.dimensions["s_rho"])
        Dcrit = nc.variables["Dcrit"][:]
        zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])

except FileNotFoundError:
    print(f"Error: History file '{hisname}' not found.")
    exit(1)

# Get velocity data from average file
print("Reading velocity data...")
try:
    with Dataset(avgname, "r") as nc:
        # Average file starts after initial state, so use tindex-1
        avg_tindex = max(0, tindex - 1)  # Ensure we don't go negative
        print(f"Using average time index: {avg_tindex}")
        u = np.squeeze(
            nc.variables["u"][avg_tindex, N - 1, :, :]
        )  # Surface level (N-1 in Python)
        v = np.squeeze(
            nc.variables["v"][avg_tindex, N - 1, :, :]
        )  # Surface level (N-1 in Python)

except FileNotFoundError:
    print(f"Error: Average file '{avgname}' not found.")
    exit(1)

print("Computing vorticity...")
# Compute vorticity
vort_psi = cr.vorticity(u, v, pm, pn)

# Convert from psi points to rho points (equivalent to psi2rho in MATLAB)
# psi grid has dimensions (M-1, L-1), rho grid has (M, L)
Mp, Lp = pm.shape  # rho grid dimensions: 202x202
M_psi, L_psi = vort_psi.shape  # psi grid dimensions: 201x201

print(f"Rho grid: {Mp}x{Lp}, Psi grid: {M_psi}x{L_psi}")

# Simple approach: place psi values in center of rho grid and extrapolate
vort = np.zeros((Mp, Lp))

# Place psi values in interior of rho grid (skip first row/col, place in 1:M_psi+1)
vort[1 : M_psi + 1, 1 : L_psi + 1] = vort_psi

# Fill boundaries by extrapolation
vort[0, :] = vort[1, :]  # First row = second row
vort[:, 0] = vort[:, 1]  # First column = second column

print("Processing sea level data...")
# Process zeta (sea level)
zeta[h < Dcrit] = zeta[h < Dcrit] - Dcrit
mask = np.ones_like(zeta)
mask[np.isnan(zeta)] = np.nan

print("Creating figure...")
# Create figure with same layout as MATLAB
fig = plt.figure(figsize=(12, 7))  # 800x500 pixels equivalent

# Sea Level plot
plt.subplot(1, 2, 1)
cmin, cmax = -0.7, 0.7
nbcol = 40
cint = (cmax - cmin) / nbcol
levels = np.arange(cmin, cmax + cint, cint)

# Clip values like MATLAB max(min(zeta,cmax),cmin)
zeta_plot = np.clip(zeta, cmin, cmax)
zeta_plot[zeta_plot == 0.0] = np.nan

cs1 = plt.contourf(x, y, zeta_plot, levels=levels, cmap="jet")
cbar1 = plt.colorbar(cs1)
# Set colorbar ticks to match MATLAB graduations
tick_spacing = 0.2  # Reasonable spacing for -0.7 to 0.7 range
cbar1.set_ticks(np.arange(cmin, cmax + tick_spacing, tick_spacing))
cbar1.set_ticklabels(
    [f"{x:.1f}" for x in np.arange(cmin, cmax + tick_spacing, tick_spacing)]
)

plt.axis([-250, -10, 0, 300])
plt.clim(cmin, cmax)
plt.title("Sea Level")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.gca().set_aspect("equal", adjustable="box")
plt.tick_params(labelsize=15)

# Wave-mean Vorticity plot
plt.subplot(1, 2, 2)
cmin, cmax = -0.07, 0.07
nbcol = 20
cint = (cmax - cmin) / nbcol
levels = np.arange(cmin, cmax + cint, cint)

# Clip values and apply mask
vort_plot = np.clip(vort, cmin, cmax)
vort_plot = vort_plot * mask

cs2 = plt.contourf(x, y, vort_plot, levels=levels, cmap="jet")
cbar2 = plt.colorbar(cs2)
# Set colorbar ticks to match MATLAB graduations
tick_spacing = 0.02  # Reasonable spacing for -0.07 to 0.07 range
cbar2.set_ticks(np.arange(cmin, cmax + tick_spacing, tick_spacing))
cbar2.set_ticklabels(
    [f"{x:.3f}" for x in np.arange(cmin, cmax + tick_spacing, tick_spacing)]
)

plt.axis([-250, -10, 0, 300])
plt.clim(cmin, cmax)
plt.title("Wave-mean Vorticity")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.gca().set_aspect("equal", adjustable="box")
plt.tick_params(labelsize=15)

# Adjust layout
plt.tight_layout()

# Save outputs
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "flashrip.pdf")
    plt.savefig(
        pdf_path,
        format="pdf",
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "flashrip.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
