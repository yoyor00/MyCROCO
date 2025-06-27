#!/usr/bin/env python3
"""
Make plots from the results of the EQUATOR test case

Further Information:
http://www.croco-ocean.org

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

Copyright (c) 2007 by Pierrick Penven
e-mail:Pierrick.Penven@ird.fr

Translated from Matlab by Claude (2025)
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import croco_utils as cr
import argparse
import os

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the EQUATOR test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="equator_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex", type=int, default=99, help="Time index to plot (default: 99)"
)
parser.add_argument(
    "--j", type=int, default=17, help="Y index - middle of the basin (default: 17)"
)
parser.add_argument(
    "--i", type=int, default=21, help="X index for time series (default: 21)"
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

print("=== EQUATOR Test Case ===")
print("Processing equatorial upwelling simulation results")
if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
tndx = args.tindex
j = args.j - 1  # Convert to Python 0-based indexing
i = args.i - 1  # Convert to Python 0-based indexing
fname = args.file

print(f"Using time index: {tndx}, j={j + 1}, i={i + 1} (Matlab indexing)")

# Read data
print("Reading data from CROCO file...")
with Dataset(fname, "r") as nc:
    # Check if time index exists
    time_len = len(nc.variables["scrum_time"][:])
    tndx = min(tndx, time_len - 1)  # Python 0-based
    print(f"Available time steps: {time_len}, using index: {tndx}")

    # Grid and bathymetry
    h = nc.variables["h"][:]
    x1 = nc.variables["x_rho"][:]
    y1 = nc.variables["y_rho"][:]
    x = np.squeeze(x1[j, :])

    # Sea surface height
    zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])

    # Temperature fields
    t = np.squeeze(np.mean(nc.variables["temp"][tndx, :, j : j + 2, :], axis=1))
    t2 = np.squeeze(
        np.mean(
            np.mean(nc.variables["temp"][:, :, j : j + 2, i : i + 2], axis=3), axis=2
        )
    )

    # Time
    time = nc.variables["scrum_time"][:]

    # Dimensions
    N, L = t.shape

    # Surface fields
    sst = np.squeeze(nc.variables["temp"][tndx, N - 1, :, :])
    u = cr.u2rho_2d(np.squeeze(nc.variables["u"][tndx, N - 1, :, :]))
    v = cr.v2rho_2d(np.squeeze(nc.variables["v"][tndx, N - 1, :, :]))

    # Vertical coordinates parameters
    theta_s = nc.theta_s
    theta_b = nc.theta_b
    hc = nc.hc

print("Computing vertical coordinates...")
# Compute vertical coordinates
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 1)
z2 = np.squeeze(np.mean(np.mean(zr[:, j : j + 2, i : i + 2], axis=2), axis=1))
zr = np.squeeze(zr[:, j, :])
xr = np.tile(x.reshape(1, L), (N, 1)) / 1e5

print("Creating first figure - Temperature sections...")
# First plot - Temperature sections
fig1 = plt.figure(figsize=(7, 10))

# Time evolution plot
plt.subplot(2, 1, 1)
time_years = 1980 + time / (360 * 24 * 3600)
temp_levels = np.arange(4, 32, 2)  # [4:2:30] comme MATLAB
cs1 = plt.contour(time_years, z2, t2.T, levels=temp_levels, colors="k", linewidths=0.8)
# Labelage manuel des contours comme MATLAB
manual_locations = [
    (time_years[len(time_years) // 4], z2[len(z2) // 3]),
    (time_years[len(time_years) // 2], z2[len(z2) // 2]),
    (time_years[3 * len(time_years) // 4], z2[2 * len(z2) // 3]),
]
plt.clabel(cs1, inline=True, fontsize=8, manual=False, fmt="%g")
plt.axis([None, None, -300, 0])
plt.xlabel("Time (years)")
plt.ylabel("Depth (m)")
plt.title("EQUATOR - T evolution from homogeneous conditions")
plt.grid(True)

# Equatorial section plot
plt.subplot(2, 1, 2)
# Temperature contours exactement comme MATLAB
temp_levels1 = np.arange(5, 20, 1)  # [5:1:19]
temp_levels2 = [20]  # [20 20]
temp_levels3 = np.arange(21, 31, 1)  # [21:1:30]

# Vérification des données finies (comme isfinite(C1) en MATLAB)
valid_data = np.isfinite(t[:, 1:-2])
if np.any(valid_data):
    cs2 = plt.contour(
        xr[:, 1:-2],
        zr[:, 1:-2],
        t[:, 1:-2],
        levels=temp_levels1,
        colors="k",
        linewidths=0.8,
    )
    plt.clabel(cs2, inline=True, fontsize=8, fmt="%g")

    cs3 = plt.contour(
        xr[:, 1:-2],
        zr[:, 1:-2],
        t[:, 1:-2],
        levels=temp_levels2,
        colors="r",
        linewidths=1.5,
    )
    plt.clabel(cs3, inline=True, fontsize=8, fmt="%g")

    cs4 = plt.contour(
        xr[:, 1:-2],
        zr[:, 1:-2],
        t[:, 1:-2],
        levels=temp_levels3,
        colors="k",
        linewidths=0.8,
    )
    plt.clabel(cs4, inline=True, fontsize=8, fmt="%g")

plt.axis([None, None, -300, 0])
plt.xlabel("Longitude")
plt.ylabel("Depth (m)")
plt.title("EQUATOR - Temperature [°C] Eq. section")
plt.grid(True)

plt.tight_layout()

# Save first figure
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "equator_sections.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "equator_sections.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()

print("Creating second figure - Surface fields...")
# Second plot - Surface fields
fig2 = plt.figure(figsize=(7, 10))

# SST plot avec graduations exactes de MATLAB
plt.subplot(2, 1, 1)
# MATLAB utilise 20 niveaux entre 10 et 25, mais avec caxis([10 25])
sst_levels = np.linspace(10, 25, 21)  # 20 intervalles = 21 niveaux
cs5 = plt.contourf(
    x1[1:-1, 1:-1] / 1e5,
    y1[1:-1, 1:-1] / 1e5,
    sst[1:-1, 1:-1],
    levels=sst_levels,
    cmap="jet",
)
plt.clim(10, 25)
# Colorbar avec graduations spécifiques
cbar1 = plt.colorbar(cs5)
cbar1.set_ticks(np.arange(10, 26, 2.5))  # Graduations tous les 2.5°C
cbar1.set_ticklabels([f"{x:.1f}" for x in np.arange(10, 26, 2.5)])
plt.axis([1, 39, -12, 12])
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.grid(True)
plt.title("EQUATOR - SST [°C]")

# Speed plot avec graduations exactes
plt.subplot(2, 1, 2)
speed = 100 * np.sqrt(u**2 + v**2)  # Convert to cm/s
# MATLAB utilise 20 niveaux entre 0 et 150, avec caxis([0 150])
speed_levels = np.linspace(0, 150, 21)  # 20 intervalles = 21 niveaux
cs6 = plt.contourf(
    x1[1:-1, 1:-1] / 1e5,
    y1[1:-1, 1:-1] / 1e5,
    speed[1:-1, 1:-1],
    levels=speed_levels,
    cmap="jet",
)
plt.clim(0, 150)
# Colorbar avec graduations spécifiques
cbar2 = plt.colorbar(cs6)
cbar2.set_ticks(np.arange(0, 151, 30))  # Graduations tous les 30 cm/s
cbar2.set_ticklabels([f"{x:g}" for x in np.arange(0, 151, 30)])

# Add velocity vectors (subsample exactement comme MATLAB: 2:2:end-1)
# MATLAB: quiver(x1(2:2:end-1,2:2:end-1)/1e5,y1(2:2:end-1,2:2:end-1)/1e5,...)
skip = 2
plt.quiver(
    x1[1::skip, 1::skip] / 1e5,
    y1[1::skip, 1::skip] / 1e5,
    u[1::skip, 1::skip],
    v[1::skip, 1::skip],
    color="k",
    scale=10,
    width=0.003,
    alpha=0.8,
)

plt.axis([1, 39, -12, 12])
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("EQUATOR - Speed [cm.s⁻¹]")

plt.tight_layout()

# Save second figure
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "equator_speed.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "equator_speed.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
