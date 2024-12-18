#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse

def zlevs(h, zeta, theta_s, theta_b, hc, N, vtype="r"):
    """
    Compute vertical sigma levels for the grid.
    """
    if vtype == "r":  # rho points
        Cs = (1 - theta_b) * np.sinh(theta_s * np.linspace(-1, 0, N)) / np.sinh(theta_s) + \
             theta_b * (np.tanh(theta_s * (np.linspace(-1, 0, N) + 0.5)) / (2 * np.tanh(0.5 * theta_s)))
    else:
        raise ValueError(f"Unsupported grid type: {vtype}")
    z = np.zeros((N, h.shape[0], h.shape[1]))
    for k in range(N):
        z0 = (hc * (Cs[k] - Cs[0]) + (h - hc) * Cs[k]) / (h + hc)
        z[k, :, :] = zeta + (zeta + h) * z0
    return z

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot upwelling temperature and velocity from a NetCDF file and optionally save as PDF/PNG.",
    epilog="Example usage:\n  python script.py --file upwelling_his.nc --output-dir ./plots --makepdf",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("--file", type=str, default="upwelling_his.nc", help="Path to the NetCDF file (default: upwelling_his.nc)")
parser.add_argument("--makepdf", action="store_true", help="Generate a PDF of the plot")
parser.add_argument("--makepng", action="store_true", help="Generate a PNG of the plot")
parser.add_argument("--no-show", action="store_true", help="Do not display the plots on the screen (default: display)")
parser.add_argument("--output-dir", type=str, default="./", help="Directory to save the output files (default: current directory)")
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Initialize
tndx = 4  # Python index (MATLAB's 5th index becomes 4th in Python)

# Read data from NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Variables from the NetCDF file
time = nc.variables["scrum_time"][tndx] / (24 * 3600)  # Convert seconds to days
h = nc.variables["h"][:]
y = np.squeeze(nc.variables["y_rho"][:, 1])  # MATLAB's 2nd index -> Python's 1st
zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
t = np.squeeze(nc.variables["temp"][tndx, :, :, 1])  # MATLAB's 2nd index -> Python's 1st
u = np.squeeze(nc.variables["u"][tndx, :, :, 1])  # MATLAB's 2nd index -> Python's 1st
N, M = t.shape
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
nc.close()

# Compute vertical levels
zr = zlevs(h, zeta, theta_s, theta_b, hc, N, "r")
zr = np.squeeze(zr[:, :, 0])  # MATLAB's zr(:,:,1) -> Python's zr[:, :, 0]

# Create yr
yr = y.reshape(1, M)  # MATLAB's reshape(y, 1, M)
yr = np.tile(yr, (N, 1)) / 1000  # MATLAB's repmat(yr, [N 1]) and convert to km

# Plot
plt.figure(figsize=(10, 6))

# Contourf plot for temperature
contourf = plt.contourf(yr, zr, t, np.arange(9, 21, 1), cmap="viridis", linestyles="none")
plt.colorbar(contourf)
plt.clim(10, 18)

# Contours for velocity
C1 = plt.contour(yr, zr, 100 * u, np.arange(-100, 101, 20), colors="k")
plt.clabel(C1)

# Labels and title
plt.xlabel("X [km]")
plt.ylabel("Z [m]")
plt.title(f"UPWELLING - Temp [°C] / u [cm/s] - Day = {time:.2f}")
plt.gca().tick_params(labelsize=15)

# File paths for saving
pdf_path = os.path.join(args.output_dir, "upwelling.pdf")
png_path = os.path.join(args.output_dir, "upwelling.png")

# Save to PDF if requested
if args.makepdf:
    plt.savefig(pdf_path, transparent=True, format="pdf")
    print(f"PDF file saved to: {pdf_path}")

# Save to PNG if requested
if args.makepng:
    plt.savefig(png_path, dpi=300)
    print(f"PNG file saved to: {png_path}")

# Show plot if not suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (use --no-show to enable).")