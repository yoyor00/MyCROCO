#!/usr/bin/env python3
"""
Dune3D Test Case

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

Patrick Marchesiello - 2012

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import os

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Dune3D test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="dune3d_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tdays", type=int, default=2, help="Number of days for evolution (default: 2)"
)
parser.add_argument(
    "--usgs",
    action="store_true",
    help="Use USGS model (default). If not set, uses MUSTANG model",
)
parser.add_argument(
    "--mustang", action="store_true", help="Use MUSTANG model instead of USGS"
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

print("=== Dune3D Test Case ===")
print("2D/3D dune evolution visualization")

# Determine model type
if args.mustang:
    usgs = False
    model_name = "MUSTANG"
    var_hmorph = "Hm"
    print("Using MUSTANG model")
else:
    usgs = True  # Default
    model_name = "USGS"
    var_hmorph = "hmorph"
    print("Using USGS model")

if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
fname = args.file
tdays = args.tdays
tndx = tdays + 1  # MATLAB uses 1-based indexing, but we need to account for t=0

print(f"Using file: {fname}")
print(f"Evolution time: {tdays} days")
print(f"Time index: {tndx}")

# Read data
print("Reading data from CROCO file...")
try:
    with Dataset(fname, "r") as nc:
        # Grid coordinates
        h = np.squeeze(nc.variables["h"][:, :])
        X = np.squeeze(nc.variables["x_rho"][:, :])
        Y = np.squeeze(nc.variables["y_rho"][:, :])

        # Convert from masked arrays if needed
        if hasattr(X, "mask"):
            X = np.asarray(X)
        if hasattr(Y, "mask"):
            Y = np.asarray(Y)
        if hasattr(h, "mask"):
            h = np.asarray(h)

        print(f"Grid dimensions: {X.shape}")

        # Read morphological height for chosen model
        if var_hmorph in nc.variables:
            print(f"Using variable: {var_hmorph}")

            # Check time dimension
            time_len = nc.variables[var_hmorph].shape[0]
            tndx = min(tndx, time_len - 1)  # Ensure valid index
            print(f"Available time steps: {time_len}, using index: {tndx}")

            hmorph0 = np.squeeze(nc.variables[var_hmorph][0, :, :])  # Initial state
            hmorph = np.squeeze(nc.variables[var_hmorph][tndx, :, :])  # Final state

            # Convert from masked arrays if needed
            if hasattr(hmorph0, "mask"):
                hmorph0 = np.asarray(hmorph0)
            if hasattr(hmorph, "mask"):
                hmorph = np.asarray(hmorph)

        else:
            print(f"Warning: Variable '{var_hmorph}' not found, trying alternative...")
            # Try alternative names
            if "hmorph" in nc.variables:
                var_hmorph = "hmorph"
                hmorph0 = np.squeeze(nc.variables["hmorph"][0, :, :])
                hmorph = np.squeeze(nc.variables["hmorph"][tndx, :, :])
            elif "Hm" in nc.variables:
                var_hmorph = "Hm"
                hmorph0 = np.squeeze(nc.variables["Hm"][0, :, :])
                hmorph = np.squeeze(nc.variables["Hm"][tndx, :, :])
            else:
                raise KeyError(f"Neither 'hmorph' nor 'Hm' found in file")

            print(f"Using variable: {var_hmorph}")

except FileNotFoundError:
    print(f"Error: File '{fname}' not found.")
    exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

print("Creating 2D bed evolution plots...")
# Plot 2D bed evolution
fig = plt.figure(figsize=(14, 7))  # 1000x500 equivalent

# Plot at t=0 (left subplot)
ax1 = plt.subplot(1, 2, 1)

# Pcolor plot like MATLAB - use 'nearest' shading to match MATLAB behavior
pc1 = plt.pcolormesh(X, Y, -hmorph0, shading="nearest", cmap="viridis")
plt.clim(-4, -2)  # caxis([-4 -2])
cbar1 = plt.colorbar(pc1)
cbar1.ax.tick_params(labelsize=13)

plt.xlabel("X (m)", fontsize=12)
plt.ylabel("Y (m)", fontsize=12)

# Contour lines like MATLAB
contour_levels = np.arange(-3.8, -1.6, 0.4)  # [-3.8:0.4:-2]
cs1 = plt.contour(X, Y, -hmorph0, levels=contour_levels, colors="k", linewidths=0.8)
plt.clabel(
    cs1, levels=contour_levels, inline=True, fontsize=7, fmt="%.1f", colors="red"
)

plt.title("Hm (t0=0)", fontsize=14)
plt.tick_params(labelsize=15)

# Plot at t=tdays (right subplot)
ax2 = plt.subplot(1, 2, 2)

# Pcolor plot like MATLAB - use 'nearest' shading to match MATLAB behavior
pc2 = plt.pcolormesh(X, Y, -hmorph, shading="nearest", cmap="viridis")
plt.clim(-4, -2)  # caxis([-4 -2])
cbar2 = plt.colorbar(pc2)
cbar2.ax.tick_params(labelsize=13)
cbar2.set_label("Hm", fontsize=13)  # cbr.Label.String = 'Hm'

plt.xlabel("X (m)", fontsize=12)

# Contour lines like MATLAB
cs2 = plt.contour(X, Y, -hmorph, levels=contour_levels, colors="k", linewidths=0.8)
plt.clabel(
    cs2, levels=contour_levels, inline=True, fontsize=7, fmt="%.1f", colors="red"
)

plt.title(f"Hm (t= t0 + {tdays} days) / {model_name}", fontsize=14)
plt.tick_params(labelsize=15)

# Adjust layout
plt.tight_layout()

# Save outputs
output_suffix = "mustang" if not usgs else "usgs"

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, f"dune3d_{output_suffix}.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, f"dune3d_{output_suffix}.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
print(f"Plotted dune evolution over {tdays} days using {model_name} model")

# Print some statistics
print("\n=== Evolution Statistics ===")
initial_depth = np.mean(-hmorph0)
final_depth = np.mean(-hmorph)
depth_change = final_depth - initial_depth

print(f"Initial mean depth: {initial_depth:.3f} m")
print(f"Final mean depth: {final_depth:.3f} m")
print(f"Mean depth change: {depth_change:.3f} m")

# Morphological statistics
morph_change = hmorph - hmorph0
max_erosion = np.min(morph_change)
max_deposition = np.max(morph_change)
rms_change = np.sqrt(np.mean(morph_change**2))

print(f"Maximum erosion: {max_erosion:.3f} m")
print(f"Maximum deposition: {max_deposition:.3f} m")
print(f"RMS morphological change: {rms_change:.3f} m")
