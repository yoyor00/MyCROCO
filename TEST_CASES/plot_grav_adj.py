#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the GRAV_ADJ test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="gravadj_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex", type=int, nargs="+", default=[10, 20, 50], help="Time indices to plot"
)
parser.add_argument(
    "--nbq",
    type=int,
    choices=[0, 1],
    default=0,
    help="0: hydrostatic, 1: non-hydrostatic",
)
parser.add_argument("--plot-psi", action="store_true", help="Plot streamfunction")
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
    axs = [axs]  # List even for 1 subplot

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

# Loop on time
for i, tndx in enumerate(args.tindex):
    tndx = min(tndx, len(nc.variables["scrum_time"][:]) - 1)
    print(f"Processing time index: {tndx}")

    # Read data
    zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
    temp = np.squeeze(nc.variables["temp"][tndx, :, 1, :])
    N, M = temp.shape
    w = 1000 * np.squeeze(nc.variables["w"][tndx, :N, 1, :])

    # Compute depth
    zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 2)
    zr = np.squeeze(zr[:, 1, :])
    xr = x if args.nbq else x / 1000 - 32
    xr_2d = np.tile(xr, (N, 1))

    # Compute current function (psi)
    psi = np.zeros_like(w)
    for j in range(1, M):
        psi[:, j] = psi[:, j - 1] - w[:, j] * (xr_2d[:, j] - xr_2d[:, j - 1])

    # Replace 0 by Nan in temperature field
    temp[temp == 0] = np.nan

    # Plot
    ax = axs[i]
    levels = np.arange(10, 41, 1)
    cont = ax.contourf(xr_2d, zr, temp, levels=levels, cmap="jet", extend="both")
    ax.set_title(f"Gravitational Adjustment - Time Index {tndx}" if i == 0 else "")
    ax.set_xlabel("X [km]" if not args.nbq else "X [m]")
    ax.set_ylabel("Depth [m]")
    ax.set_ylim(np.min(zr), np.max(zr))
    ax.set_xlim(np.min(xr), np.max(xr))
    if args.plot_psi:
        ax.contour(xr_2d, zr, psi, colors="k", linewidths=0.5)
    fig.colorbar(cont, ax=ax)

# Save in PDF
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "gravadj_results.pdf")
    plt.savefig(pdf_path)
    print(f"PDF saved: {pdf_path}")

# Save in PNG
# Save in PNG (one file per time index, full figure)
if args.makepng:
    for tndx in args.tindex:
        tndx = min(tndx, len(nc.variables["scrum_time"][:]) - 1)
        zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
        temp = np.squeeze(nc.variables["temp"][tndx, :, 1, :])
        N, M = temp.shape
        zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 2)
        zr = np.squeeze(zr[:, 1, :])
        xr_loop = x if args.nbq else x / 1000 - 32
        xr_2d = np.tile(xr_loop, (N, 1))
        temp[temp == 0] = np.nan

        fig_png, ax_png = plt.subplots(figsize=(10, 4), constrained_layout=True)
        levels = np.arange(10, 41, 1)
        cont = ax_png.contourf(xr_2d, zr, temp, levels=levels,
                                cmap="jet", extend="both")
        ax_png.set_title(f"Gravitational Adjustment - Time Index {tndx}")
        ax_png.set_xlabel("X [km]" if not args.nbq else "X [m]")
        ax_png.set_ylabel("Depth [m]")
        fig_png.colorbar(cont, ax=ax_png)

        png_path = os.path.join(args.output_dir, f"gravadj_t{tndx:04d}.png")
        fig_png.savefig(png_path, dpi=150)
        plt.close(fig_png)
        print(f"PNG saved: {png_path}")

# Show or not
if not args.no_show:
    plt.show()
else:
    plt.close()

# Close NetCDF file
nc.close()
