#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SHOREFACE test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="shoreface_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--yindex",
    type=int,
    default=0,
    help="Y index to extract (0-based indexing, default=0)",
)
parser.add_argument(
    "--tindex",
    type=int,
    default=None,
    help="Time index to plot (default: last time step)",
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as PDF files")
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
parser.add_argument(
    "--plot-setup", action="store_true", help="Also plot wave setup figure"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Get dimensions and time info
print("Reading grid and time information...")
scrum_time = nc.variables["scrum_time"][:]
if args.tindex is None:
    tindex = len(scrum_time) - 1  # Last record
else:
    tindex = min(args.tindex, len(scrum_time) - 1)

# Read grid data
hr = np.squeeze(nc.variables["h"][args.yindex, :])
xindex = 0  # Start from beginning in Python (0-based)
hr = hr[xindex:]
L = len(hr)

xr = np.squeeze(nc.variables["x_rho"][args.yindex, xindex:])
yr = np.squeeze(nc.variables["y_rho"][args.yindex, xindex:])
dx = xr[1] - xr[0]

# Vertical grid parameters
N = len(nc.dimensions["s_rho"])
theta_s = float(nc.theta_s)
theta_b = float(nc.theta_b)
hc = float(nc.hc)

zeta = np.squeeze(nc.variables["zeta"][tindex, args.yindex, xindex:])
zr = np.squeeze(cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "r", 2))
zru = 0.5 * (zr[:, :-1] + zr[:, 1:])
zw = np.squeeze(cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "w", 2))
zwu = 0.5 * (zw[:, :-1] + zw[:, 1:])

# Create 2D coordinate grids
xr2d = np.tile(xr, (N, 1))  # Rho-grid coordinates (N x L)
xw2d = np.tile(xr, (N + 1, 1))  # For w-grid vertical levels ((N+1) x L)

# Create u-grid coordinates (one less column)
xr_u = 0.5 * (xr[:-1] + xr[1:])  # U-grid x coordinates
xr2d_u = np.tile(xr_u, (N, 1))  # U-grid 2D coordinates (N x (L-1))

D = hr + zeta
D2d = np.tile(D, (N, 1))

# Read model fields
print("Reading model fields...")
time = scrum_time[tindex] / 86400.0  # Convert to days

zeta1 = zeta.copy()

# Velocities
u1 = np.squeeze(nc.variables["u"][tindex, :, args.yindex, xindex:])
v1 = np.squeeze(nc.variables["v"][tindex, :, args.yindex, xindex:])
w1 = np.squeeze(nc.variables["w"][tindex, :, args.yindex, xindex:])

# Temperature
t1 = np.squeeze(nc.variables["temp"][tindex, :, args.yindex, xindex:])

# Vertical viscosity/diffusivity
Akv = np.squeeze(nc.variables["AKv"][tindex, :, args.yindex, xindex:])
Akt = np.squeeze(nc.variables["AKt"][tindex, :, args.yindex, xindex:])

# Wave setup
sup = np.squeeze(nc.variables["sup"][tindex, args.yindex, xindex:])

# Stokes drift
ust = np.squeeze(nc.variables["ust"][tindex, :, args.yindex, xindex:])
vst = np.squeeze(nc.variables["vst"][tindex, :, args.yindex, xindex:])
wst = np.squeeze(nc.variables["wst"][tindex, :, args.yindex, xindex:])

# Wave-induced mixing
Akb = np.squeeze(nc.variables["Akb"][tindex, :, args.yindex, xindex:])
Akw = np.squeeze(nc.variables["Akw"][tindex, :, args.yindex, xindex:])

nc.close()


def create_custom_colormap():
    """Create custom colormap with white in the middle"""
    nbcol = 20
    jet_colors = plt.cm.jet(np.linspace(0, 1, nbcol))
    jet_colors[nbcol // 2 - 1, :] = [1, 1, 1, 1]  # White
    jet_colors[nbcol // 2, :] = [1, 1, 1, 1]  # White
    return plt.matplotlib.colors.ListedColormap(jet_colors)


# ============================================================
# Data processing and masking
# ============================================================
print("Processing data...")

Dcrit = 0.2

# Convert coordinates (subtract 1 km and convert to km)
xr = 1.0e-3 * xr - 1
xr_u = 1.0e-3 * xr_u - 1  # Also convert u-grid coordinates
xr2d = 1.0e-3 * xr2d - 1
xr2d_u = 1.0e-3 * xr2d_u - 1
xw2d = 1.0e-3 * xw2d - 1

# Boundary conditions and masking
u1[:, -1] = u1[:, -2]  # Last column = second to last
ust[:, -1] = ust[:, -2]

# Masking for different grilles
zeta1[D < Dcrit] = np.nan

# For u variables (grille u): convert D to u-grid
D_u = 0.5 * (D[:-1] + D[1:])  # Average to u-grid points
D2d_u = np.tile(D_u, (N, 1))
u1[D2d_u < Dcrit] = np.nan

# For v,w variables (grille rho): use original D2d
v1[D2d < Dcrit] = np.nan

# Mask first column for viscosity fields
Akv[:, 0] = np.nan
Akt[:, 0] = np.nan
Akb[:, 0] = np.nan
Akw[:, 0] = np.nan

# ============================================================
# Main plot - 3x3 subplot layout
# ============================================================
print("Creating plots...")

fig = plt.figure(figsize=(16, 12))

thour = int(time * 24)

# Plot parameters for each subplot
plot_params = [
    # Row 1: Eulerian velocities
    {
        "data": u1,
        "cmin": -0.5,
        "cmax": 0.5,
        "label": "U [m/s]",
        "coords": (xr2d_u, zru),
        "pos": (3, 3, 1),
    },
    {
        "data": ust,
        "cmin": -0.2,
        "cmax": 0.2,
        "label": "U Stokes [m/s]",
        "coords": (xr2d_u, zru),
        "pos": (3, 3, 2),
        "title": True,
    },
    {
        "data": Akv,
        "cmin": -0.05,
        "cmax": 0.05,
        "label": "Akv (total visc) [m^2/s]",
        "coords": (xw2d, zw),
        "pos": (3, 3, 3),
        "no_zeta": True,
    },
    # Row 2: V velocities and mixing
    {
        "data": v1,
        "cmin": -1.2,
        "cmax": 1.2,
        "label": "V [m/s]",
        "coords": (xr2d, zr),
        "pos": (3, 3, 4),
    },
    {
        "data": vst,
        "cmin": -0.02,
        "cmax": 0.02,
        "label": "V Stokes [m/s]",
        "coords": (xr2d, zr),
        "pos": (3, 3, 5),
    },
    {
        "data": Akb,
        "cmin": -0.05,
        "cmax": 0.05,
        "label": "Akb (break. wave visc.)",
        "coords": (xw2d, zw),
        "pos": (3, 3, 6),
        "no_zeta": True,
    },
    # Row 3: W velocities and wave diffusivity
    {
        "data": w1,
        "cmin": -0.007,
        "cmax": 0.007,
        "label": "W [m/s]",
        "coords": (xr2d, zr),
        "pos": (3, 3, 7),
        "xlabel": True,
    },
    {
        "data": wst,
        "cmin": -0.002,
        "cmax": 0.002,
        "label": "W Stokes [m/s]",
        "coords": (xr2d, zr),
        "pos": (3, 3, 8),
        "xlabel": True,
    },
    {
        "data": Akw,
        "cmin": -0.003,
        "cmax": 0.003,
        "label": "Akw (non-break. wave diff.)",
        "coords": (xw2d, zw),
        "pos": (3, 3, 9),
        "xlabel": True,
        "no_zeta": True,
    },
]

for i, params in enumerate(plot_params):
    ax = plt.subplot(*params["pos"])

    # Create levels and contour plot
    cmin, cmax = params["cmin"], params["cmax"]
    nbcol = 20
    levels = np.linspace(cmin, cmax, nbcol + 1)

    x_coords, z_coords = params["coords"]
    cs = plt.contourf(
        x_coords,
        z_coords,
        params["data"],
        levels=levels,
        cmap=create_custom_colormap(),
        extend="both",
    )

    # Colorbar
    cbar = plt.colorbar(cs, ax=ax)

    # Format colorbar for last row (scientific notation)
    if params["pos"][2] in [7, 8, 9]:  # Last row
        cbar.formatter.set_powerlimits((0, 0))
        cbar.formatter.set_useMathText(True)
        cbar.update_ticks()

    # Plot bathymetry
    plt.plot(xr, -hr, "k-", linewidth=3)

    # Plot sea level (except for viscosity plots)
    if not params.get("no_zeta", False):
        plt.plot(xr, zeta1, "g-", linewidth=3)

    # Set axis limits and labels
    plt.xlim(-1, 0)
    plt.ylim(-12, 1)
    ax.margins(0)
    ax.autoscale(tight=True)

    plt.ylabel("Z [m]", fontsize=10)

    if params.get("xlabel", False):
        plt.xlabel("Distance to shore [km]", fontsize=10)

    if params.get("title", False):
        plt.title(f"SHOREFACE - Time = {thour} h", fontsize=15)

    # Add text label (removed from here - only using plot title)
    # plt.text(-0.05, -11, params['label'], ha='right', va='bottom', fontsize=9)

    plt.grid(True, alpha=0.3)
    plt.clim(cmin, cmax)

plt.tight_layout()

# Save main plot
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "shoreface.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"PDF saved: {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "shoreface.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG saved: {png_path}")

if not args.no_show:
    plt.show()

# ============================================================
# Wave setup plot (optional)
# ============================================================
if args.plot_setup:
    print("Creating wave setup plot...")

    plt.figure(figsize=(8, 10))

    # Surface elevation + wave setup
    xmin, xmax = -1, 0
    zmin, zmax = -0.1, 0.25

    plt.plot(xr, zeta1 + sup, "b-", linewidth=3, label="Wave Setup")
    plt.plot(xr, zeta1, "g-", linewidth=2, label="Sea Level")
    plt.plot(xr, sup, "r--", linewidth=2, label="Setup Only")

    plt.xlim(xmin, xmax)
    plt.ylim(zmin, zmax)
    plt.xlabel("Distance to shore [km]", fontsize=12)
    plt.ylabel("Elevation [m]", fontsize=12)
    plt.title(f"SHOREFACE: Wave Setup at {thour} h", fontsize=14)
    plt.grid(True)
    plt.legend()

    # Remove margins for exact limits
    ax = plt.gca()
    ax.margins(0)
    ax.autoscale(tight=True)

    # Save wave setup plot
    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "shoreface_setup.pdf")
        plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
        print(f"Wave setup PDF saved: {pdf_path}")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "shoreface_setup.png")
        plt.savefig(png_path, bbox_inches="tight", dpi=300)
        print(f"Wave setup PNG saved: {png_path}")

    if not args.no_show:
        plt.show()

print("Script completed successfully!")
