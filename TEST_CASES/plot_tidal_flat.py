#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the TFLAT2DV test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="tidal_flat_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--idy", type=int, default=0, help="Y index to extract (0-based indexing)"
)
parser.add_argument(
    "--idz", type=int, default=0, help="Z index to extract (0-based indexing)"
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

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Read data
print("Reading data...")
h = np.squeeze(nc.variables["h"][args.idy, :])
Dcrit = float(nc.variables["Dcrit"][:])  # Dcrit is a scalar
X = np.squeeze(nc.variables["x_rho"][args.idy, :]) / 1000.0  # Convert to km
T = np.squeeze(nc.variables["scrum_time"][:]) / 86400.0  # Convert to days
sand = np.squeeze(nc.variables["SAND"][:, args.idz, args.idy, :])
zeta = np.squeeze(nc.variables["zeta"][:, args.idy, :])
ubar_raw = np.squeeze(nc.variables["ubar"][:, args.idy, :])

# Close NetCDF file
nc.close()

# Convert ubar from u-grid to rho-grid
print("Converting ubar to rho-grid...")
ubar_list = []
for t in range(ubar_raw.shape[0]):
    ubar_rho = cr.u2rho_1d(ubar_raw[t, :])
    ubar_list.append(ubar_rho)
ubar = np.array(ubar_list)

# Time adjustment
T1 = T[0]
T = T - T1

# Calculate total depth and apply masking
print("Applying wet/dry masking...")
D = zeta + np.tile(h, (zeta.shape[0], 1))

# Create masks for wet/dry conditions (Dcrit is a scalar)
wet_mask = D >= (Dcrit + 0.01)

# Apply masking (set dry areas to NaN)
sand_masked = sand.copy()
zeta_masked = zeta.copy()
ubar_masked = ubar.copy()

sand_masked[~wet_mask] = np.nan
zeta_masked[~wet_mask] = np.nan
ubar_masked[~wet_mask] = np.nan

# Create the plot
print("Creating Hovmöller plots...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))


# Create custom colormap with white for NaN values
def create_colormap_with_white():
    """Create a jet colormap with white for missing values"""
    cmap = plt.cm.jet
    cmap.set_bad("white")
    return cmap


# Create meshgrid for plotting (using cell centers)
X_2d, T_2d = np.meshgrid(X, T)

# Plot 1: Ubar
im1 = ax1.pcolormesh(
    X_2d, T_2d, ubar_masked, shading="auto", cmap=create_colormap_with_white()
)
im1.set_clim(-3, 3)
cbar1 = plt.colorbar(im1, ax=ax1)
cbar1.set_label("m/s", fontsize=12)
ax1.set_xlabel("X (km)", fontsize=12)
ax1.set_ylabel("Time (days)", fontsize=12)
ax1.set_title("Tidal Flat / Ubar (m/s) hovmoller", fontsize=14)
ax1.tick_params(labelsize=12)

# Plot 2: Sand concentration
im2 = ax2.pcolormesh(
    X_2d, T_2d, sand_masked, shading="auto", cmap=create_colormap_with_white()
)
im2.set_clim(0, 1)
cbar2 = plt.colorbar(im2, ax=ax2)
cbar2.set_label("Kg/m³", fontsize=12)
ax2.set_xlabel("X (km)", fontsize=12)
ax2.set_ylabel("Time (days)", fontsize=12)
ax2.set_title("Tidal Flat / Sand concentration (Kg/m³) hovmoller", fontsize=14)
ax2.tick_params(labelsize=12)

plt.tight_layout()

# Save as PDF
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "tidal_flat.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"PDF saved: {pdf_path}")

# Save as PNG
if args.makepng:
    png_path = os.path.join(args.output_dir, "tidal_flat.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG saved: {png_path}")

# Show or close
if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
