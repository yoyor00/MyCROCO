#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the ISOLITON test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="isoliton_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex", type=int, nargs="+", default=[30, 50, 70], help="Time indices to plot"
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as a PDF file")
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Prepare figure
nplot = len(args.tindex)
fig, axs = plt.subplots(nplot, 1, figsize=(10, 3 * nplot), constrained_layout=True)
if nplot == 1:
    axs = [axs]  # Assure une liste, même pour un seul subplot

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Read grid data
h = nc.variables["h"][:]
x = np.squeeze(nc.variables["x_rho"][1, :])
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
N = len(nc.dimensions["s_rho"])

# Loop on time
for i, tndx in enumerate(args.tindex):
    tndx = min(tndx, len(nc.variables["scrum_time"][:]) - 1)
    print(f"Processing time index: {tndx}")

    # Read data
    zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
    temp = np.squeeze(nc.variables["temp"][tndx, :, 1, :])
    w = 1000 * np.squeeze(nc.variables["w"][tndx, :, 1, :])
    temp[temp == 0] = np.nan

    # Compute depth
    zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 2)
    zr = np.squeeze(zr[:, 1, :])
    xr = np.tile(x, (N, 1))

    # Plot
    ax = axs[i]
    levels = np.arange(-40, 41, 1)
    cont = ax.contourf(xr, zr, temp, levels=levels, cmap="jet", extend="both")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Depth [m]")
    ax.set_ylim(-0.25, -0.15)
    ax.set_xlim(np.min(xr), np.max(xr))
    ax.set_title("Internal Soliton" if i == 0 else "")
    fig.colorbar(cont, ax=ax)

# Save in PDF
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "isoliton_results.pdf")
    plt.savefig(pdf_path)
    print(f"PDF saved: {pdf_path}")

# Save in PNG (one file per time index, full figure)
if args.makepng:
    for tndx in args.tindex:
        tndx = min(tndx, len(nc.variables["scrum_time"][:]) - 1)
        zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
        temp = np.squeeze(nc.variables["temp"][tndx, :, 1, :])
        temp[temp == 0] = np.nan
        zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 2)
        zr = np.squeeze(zr[:, 1, :])
        xr_loop = np.tile(x, (N, 1))

        fig_png, ax_png = plt.subplots(figsize=(10, 4), constrained_layout=True)
        levels = np.arange(-40, 41, 1)
        cont = ax_png.contourf(xr_loop, zr, temp, levels=levels,
                                cmap="jet", extend="both")
        ax_png.set_xlabel("X [m]")
        ax_png.set_ylabel("Depth [m]")
        ax_png.set_ylim(-0.25, -0.15)
        ax_png.set_xlim(np.min(xr_loop), np.max(xr_loop))
        ax_png.set_title(f"Internal Soliton - Time Index {tndx}")
        fig_png.colorbar(cont, ax=ax_png)

        png_path = os.path.join(args.output_dir, f"isoliton_t{tndx:04d}.png")
        fig_png.savefig(png_path, dpi=300)
        plt.close(fig_png)
        print(f"PNG saved: {png_path}")

# Show or not
if not args.no_show:
    plt.show()
else:
    plt.close()

# Close NetCDF file
nc.close()
