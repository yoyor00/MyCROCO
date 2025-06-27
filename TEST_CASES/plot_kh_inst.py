#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Kelvin-Helmholtz Instability test case (KH_INST).",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="khinst_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex",
    type=int,
    default=None,
    help="Time index to plot (default: last time step)",
)
parser.add_argument("--makemovie", action="store_true", help="Create animation frames")
parser.add_argument("--makepdf", action="store_true", help="Save plots as PDF files")
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
parser.add_argument(
    "--cmin", type=float, default=15.5, help="Minimum contour value (default: 15.5)"
)
parser.add_argument(
    "--cmax", type=float, default=18.5, help="Maximum contour value (default: 18.5)"
)
parser.add_argument(
    "--nbcol", type=int, default=100, help="Number of contour levels (default: 100)"
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

# Set time range for movie or single plot
if args.makemovie:
    tstr = 0
    tend = tindex
    print(f"Creating animation from time {tstr} to {tend}")
else:
    tstr = tindex
    tend = tstr

# Initialize grid data (read once for efficiency)
print("Reading grid data...")
h2d = nc.variables["h"][:, :]
h = np.squeeze(nc.variables["h"][1, :])  # y-index 1 (Python: 0-based)
xl = float(nc.variables["xl"][:])  # xl is a variable, not an attribute
x = np.squeeze(nc.variables["x_rho"][1, :])  # y-index 1
pm = np.squeeze(nc.variables["pm"][1, :])

N = len(nc.dimensions["s_rho"])
M = len(x)

# Vertical grid parameters
theta_s = float(nc.theta_s)
theta_b = float(nc.theta_b)
hc = float(nc.hc)

# Read initial density to get grid dimensions
zeta2d = np.squeeze(nc.variables["zeta"][tstr, :, :])
rho_init = np.squeeze(nc.variables["rho"][tstr, :, 1, :])  # y-index 1

# Compute vertical grid (do this once since grid doesn't change much)
z = cr.zlevs(h2d, zeta2d, theta_s, theta_b, hc, N, "r", 2)
z = np.squeeze(z[:, 1, :])  # Extract y-index 1

# Create 2D coordinate grid
x_2d = np.tile(x, (N, 1))

print("Creating plots...")

# No interactive mode - just like Matlab default behavior

# Main plotting loop
for t in range(tstr, tend + 1):
    print(f"Processing time index: {t}")

    # Read data for current time
    time_minutes = scrum_time[t] / 60.0  # Convert to minutes
    rho = np.squeeze(nc.variables["rho"][t, :, 1, :])  # y-index 1 (middle slice)

    # Create figure
    fig = plt.figure(figsize=(10, 7))
    plt.rcParams.update({"font.size": 15})

    # Contour plot with smooth gradations (not discrete levels)
    cs = plt.contourf(x_2d, z, rho, levels=100, cmap="jet")
    plt.shading = "flat"  # Equivalent to 'shading flat' in Matlab

    # Colorbar with specific ticks but smooth gradations
    cbar = plt.colorbar(cs)
    cbar.set_ticks(np.arange(args.cmin, args.cmax + 0.5, 0.5))  # Every 0.5
    cbar.set_label("Density [kg/m³]", rotation=270, labelpad=20)

    # Set axis limits like Matlab: axis([0 256 -256 0])
    plt.xlim(0, 256)
    plt.ylim(-256, 0)
    ax = plt.gca()
    ax.margins(0)
    ax.autoscale(tight=True)

    # Set Y-axis ticks every 50 from 0 to -250
    y_ticks = np.arange(0, -251, -50)  # [0, -50, -100, -150, -200, -250]
    plt.yticks(y_ticks)

    # Set color limits
    plt.clim(args.cmin, args.cmax)

    # Format time display
    thour = int(time_minutes // 60)
    tmin = int(time_minutes % 60)
    clock = f"{thour} h {tmin} min"

    # Labels and title
    plt.title(f"Time: {clock}", fontsize=15)
    plt.xlabel("X [m]", fontsize=12)
    plt.ylabel("Z [m]", fontsize=12)

    # Set background color
    fig.patch.set_facecolor("white")

    print(f"Plot for time: {clock} (Record: {t})")

    # Debug: print what we're going to do
    print(f"makemovie: {args.makemovie}, no_show: {args.no_show}")

    if args.makemovie:
        # Save frame for movie
        frame_path = os.path.join(args.output_dir, f"khinst_frame_{t:04d}.png")
        plt.savefig(frame_path, dpi=150, bbox_inches="tight", facecolor="white")
        print(f"Frame saved: {frame_path}")
        plt.close()
    else:
        # For single plot - always show unless explicitly disabled
        # Save files if requested
        if args.makepdf:
            pdf_path = os.path.join(args.output_dir, "khinst.pdf")
            plt.savefig(pdf_path, bbox_inches="tight", facecolor="white", dpi=300)
            print(f"PDF saved: {pdf_path}")

        if args.makepng:
            png_path = os.path.join(args.output_dir, "khinst.png")
            plt.savefig(png_path, bbox_inches="tight", facecolor="white", dpi=300)
            print(f"PNG saved: {png_path}")

        # Show the plot (unless --no-show is specified)
        if not args.no_show:
            print("Showing plot...")
            plt.show()

        break  # Exit loop for single plot

# Close NetCDF file
nc.close()

# Create animation summary if movie frames were generated
if args.makemovie:
    print(f"\nAnimation frames created: {tend - tstr + 1} frames")
    print(f"Frames saved in: {args.output_dir}")
    print("To create a movie, you can use:")
    print(
        f"ffmpeg -r 5 -i {args.output_dir}/khinst_frame_%04d.png -c:v libx264 -pix_fmt yuv420p khinst.mp4"
    )

print("Script completed successfully!")
