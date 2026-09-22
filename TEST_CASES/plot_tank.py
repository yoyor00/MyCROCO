#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the TANK test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="tank_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--nbq",
    type=int,
    default=1,
    choices=[0, 1],
    help="Solution type: 0=hydrostatic, 1=non-hydrostatic (default=1)",
)
parser.add_argument(
    "--yindex",
    type=int,
    default=1,
    help="Y index to extract (0-based indexing, default=1)",
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
    "--timeseries-only",
    action="store_true",
    help="Only plot time series (no animation)",
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Constants
g = 9.81  # gravity acceleration (m/s²)
eta = 0.001  # nondimensional periodic amplitude
D0 = 10.0  # tank depth (m)
Lt = 10.0  # tank length (m)
rho0 = 1024.4  # reference density

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
L = len(hr)
xr = np.squeeze(nc.variables["x_rho"][args.yindex, :])
yr = np.squeeze(nc.variables["y_rho"][args.yindex, :])
dx = xr[1] - xr[0]

# Vertical grid parameters
N = len(nc.dimensions["s_rho"])
theta_s = float(nc.theta_s)
theta_b = float(nc.theta_b)
hc = float(nc.hc)


def compute_analytical_solutions_nh(xr_2d, zr, time, k, D0):
    """Compute non-hydrostatic analytical solutions"""
    sig = np.sqrt(g * k * np.tanh(k * D0))

    # From Matlab script exactly
    zeta = eta * np.cos(k * xr_2d) * np.cos(sig * time)
    u = (
        eta
        * sig
        * (np.sin(sig * time) / np.sinh(k * D0))
        * np.sin(k * xr_2d)
        * np.cosh(k * (D0 + zr))
    )
    w = (
        -eta
        * sig
        * (np.sin(sig * time) / np.sinh(k * D0))
        * np.cos(k * xr_2d)
        * np.sinh(k * (D0 + zr))
    )

    return zeta, u, w, sig


def compute_analytical_solutions_h(xr_2d, zr, time, k, D0):
    """Compute hydrostatic analytical solutions"""
    sig = k * np.sqrt(g * D0)

    # From Matlab script exactly
    zeta = eta * np.cos(k * xr_2d - sig * time)
    u = g * eta * k / sig * np.sin(k * xr_2d) * np.sin(sig * time)
    w = -g * eta * k**2 / sig * np.cos(k * xr_2d) * np.sin(sig * time) * (D0 + zr)

    return zeta, u, w, sig


def create_custom_colormap():
    """Create custom colormap with white in the middle"""
    nbcol = 10
    jet_colors = plt.cm.jet(np.linspace(0, 1, nbcol))
    jet_colors[nbcol // 2 - 1, :] = [1, 1, 1, 1]  # White
    jet_colors[nbcol // 2, :] = [1, 1, 1, 1]  # White
    return plt.matplotlib.colors.ListedColormap(jet_colors)


# ============================================================
# Plot time series
# ============================================================
print("Creating time series plots...")

# Read time series data
kk = round(N / 2)  # for w
t0 = scrum_time[: tindex + 1]
zeta01 = 100 * np.squeeze(
    nc.variables["zeta"][: tindex + 1, args.yindex, -1]
)  # Last point
u01 = 100 * np.squeeze(
    nc.variables["u"][: tindex + 1, -1, args.yindex, L // 2]
)  # Middle point
w01 = 100 * np.squeeze(
    nc.variables["w"][: tindex + 1, kk, args.yindex, -1]
)  # Last point

# Set vertical grid for all time steps
print("Computing vertical grids...")
zr0 = np.zeros((tindex + 1, N, L))
zw0 = np.zeros((tindex + 1, N + 1, L))

for i in range(tindex + 1):
    zeta0 = np.squeeze(nc.variables["zeta"][i, args.yindex, :])
    zr0[i, :, :] = np.squeeze(cr.zlevs(hr, zeta0, theta_s, theta_b, hc, N, "r", 2))
    zw0[i, :, :] = np.squeeze(cr.zlevs(hr, zeta0, theta_s, theta_b, hc, N, "w", 2))

zru0 = 0.5 * (zr0[:, :, :-1] + zr0[:, :, 1:])

# Compute analytical solutions for time series
x = xr[-1]  # Last point
xu = 0.5 * (xr[L // 2] + xr[L // 2 + 1])
z0 = np.squeeze(zru0[:, -1, -1])  # Last vertical level, last horizontal point
zz = np.squeeze(zw0[:, kk + 1, -1])  # Last horizontal point

k = np.pi / Lt

# Non-hydrostatic case
sig_nh = np.sqrt(g * k * np.tanh(k * D0))
zeta02 = 100 * eta * np.cos(k * x - sig_nh * t0)
u02 = (
    100
    * eta
    * sig_nh
    * (np.sin(sig_nh * t0) / np.sinh(k * D0))
    * np.sin(k * xu)
    * np.cosh(k * (D0 + z0))
)
w02 = (
    -100
    * eta
    * sig_nh
    * (np.sin(sig_nh * t0) / np.sinh(k * D0))
    * np.cos(k * x)
    * np.sinh(k * (D0 + zz))
)

# Hydrostatic case
sig_h = k * np.sqrt(g * D0)
T_lw = 2 * np.pi / sig_h
zeta03 = 100 * eta * np.cos(k * x) * np.cos(sig_h * t0)
u03 = 100 * g * eta * k / sig_h * np.sin(k * xu) * np.sin(sig_h * t0)
w03 = -100 * g * eta * k**2 / sig_h * np.cos(k * x) * np.sin(sig_h * t0) * (D0 + zz)

t_per = t0 / T_lw

# Create time series plot
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

# Plot zeta
ax1.plot(t_per, zeta03, "k-", linewidth=2, label="Analytical hydro")
ax1.plot(t_per, zeta02, "k--", linewidth=2, label="Analytical N-hydro")
if args.nbq:
    ax1.plot(t_per, zeta01, "r-", linewidth=2, label="Numerical N-hydro")
else:
    ax1.plot(t_per, zeta01, "r-", linewidth=2, label="Numerical hydro")
ax1.set_xlabel("time (periods)", fontsize=12)
ax1.set_ylabel("zeta (cm)", fontsize=12)
ax1.set_title("TANK test case", fontsize=14)
ax1.grid(True)
ax1.legend()
# Remove margins for tight axes
ax1.margins(0)
ax1.autoscale(tight=True)

# Plot u
ax2.plot(t_per, u03, "k-", linewidth=2)
ax2.plot(t_per, u02, "k--", linewidth=2)
ax2.plot(t_per, u01, "r-", linewidth=2)
ax2.set_xlabel("time (periods)", fontsize=12)
ax2.set_ylabel("u (cm/s)", fontsize=12)
ax2.grid(True)
# Remove margins for tight axes
ax2.margins(0)
ax2.autoscale(tight=True)

# Plot w
ax3.plot(t_per, w03, "k-", linewidth=2)
ax3.plot(t_per, w02, "k--", linewidth=2)
ax3.plot(t_per, w01, "r-", linewidth=2)
ax3.set_xlabel("time (periods)", fontsize=12)
ax3.set_ylabel("w (cm/s)", fontsize=12)
ax3.grid(True)
# Remove margins for tight axes
ax3.margins(0)
ax3.autoscale(tight=True)

plt.tight_layout()

# Save time series plot
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "tank_series.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"Time series PDF saved: {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "tank_series.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"Time series PNG saved: {png_path}")

if not args.no_show:
    plt.show()

# Exit if only time series requested
if args.timeseries_only:
    nc.close()
    print("Script completed successfully!")
    exit(0)

# ============================================================
# Animation/snapshot plots
# ============================================================

# Set time range for movie or single plot
if args.makemovie:
    tstr = 0
    tend = tindex
    print(f"Creating animation from time {tstr} to {tend}")
else:
    tstr = tindex
    tend = tstr

print("Creating animation/snapshot plots...")

for t in range(tstr, tend + 1):
    print(f"Processing time index: {t}")

    time = scrum_time[t]

    # Vertical grid
    zeta = np.squeeze(nc.variables["zeta"][t, args.yindex, :])
    zr = np.squeeze(cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "r", 2))
    zru = 0.5 * (zr[:, :-1] + zr[:, 1:])
    zw = np.squeeze(cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "w", 2))
    zwu = 0.5 * (zw[:, :-1] + zw[:, 1:])

    xr2d = np.tile(xr, (N, 1))
    D = hr + zeta
    D2d = np.tile(D, (N, 1))

    # Read model fields
    zeta1 = zeta.copy()
    zeta1[zeta1 == 0] = np.nan

    u1 = np.squeeze(nc.variables["u"][t, :, args.yindex, :])
    w1 = np.squeeze(nc.variables["w"][t, :, args.yindex, :])

    # Create coordinate grids
    xr2d = np.tile(xr, (N, 1))  # For rho-grid (w, zeta)

    # For u-grid, create u-grid coordinates
    xr_u = 0.5 * (xr[:-1] + xr[1:])  # u-grid x coordinates
    xr2d_u = np.tile(xr_u, (N, 1))  # For u-grid variables

    # Compute analytical solutions
    k = np.pi / Lt

    if args.nbq:  # Non-hydrostatic
        # Compute on rho-grid first, then convert to u-grid
        zeta2, u2_rho, w2, sig = compute_analytical_solutions_nh(xr2d, zr, time, k, D0)
        # Convert u2 to u-grid by averaging
        u2 = 0.5 * (u2_rho[:, :-1] + u2_rho[:, 1:])
        solution_type = "N-hydro"
    else:  # Hydrostatic
        # Compute on rho-grid first, then convert to u-grid
        zeta2, u2_rho, w2, sig = compute_analytical_solutions_h(xr2d, zr, time, k, D0)
        # Convert u2 to u-grid by averaging
        u2 = 0.5 * (u2_rho[:, :-1] + u2_rho[:, 1:])
        solution_type = "hydro"

    # Convert to cm and cm/s
    u1 = 100 * u1
    u2 = 100 * u2
    w1 = 100 * w1
    w2 = 100 * w2

    # Boundary condition
    u1[:, -1] = u1[:, -2]

    # For quiver plot, subsample the data to avoid too many arrows
    # and use original grids without complex interpolation
    skip = 3  # Show every 3rd point

    # Use u1 directly on its natural u-grid (one point less in x)
    xr_u_2d = np.tile(xr_u, (N, 1))  # u-grid coordinates

    # For w, take subset that matches u-grid size
    w1_sub = w1[:-1, :-1]  # Remove last level and last column to match u-grid

    # Subsample for cleaner quiver plot
    ur_sub = u1[::skip, ::skip] * 100  # Convert to cm/s
    wr_sub = -w1_sub[::skip, ::skip] * 100  # Convert to cm/s and test sign inversion
    xr_sub = xr_u_2d[::skip, ::skip]
    zr_sub = zr[::skip, :-1:skip]  # Match u-grid in x

    # Plot
    fig = plt.figure(figsize=(12, 8))

    # Contour plot of velocity error (use u-grid coordinates)
    cmin, cmax = -0.02, 0.02
    levels = np.linspace(cmin, cmax, 11)

    cs = plt.contourf(
        xr2d_u,
        zr[:, :-1],
        u1 - u2,
        levels=levels,
        cmap=create_custom_colormap(),
        extend="both",
    )

    # Quiver plot of velocities (subsampled for clarity)
    plt.quiver(xr_sub, zr_sub, ur_sub, wr_sub, scale=50, alpha=0.7, width=0.003)

    # Colorbar
    plt.colorbar(cs, label="U error (cm/s)")

    # Sea level (convert to mm like in Matlab: 1000*)
    plt.plot(xr, 1000 * zeta2[0, :], "g-", linewidth=2, label="Analytical")
    plt.plot(xr, 1000 * zeta1, "r-", linewidth=2, label="Numerical")

    # Set exact limits like in Matlab: axis([0 10 -10 1])
    ax = plt.gca()
    ax.set_xlim(0, 10)
    ax.set_ylim(-10, 1)

    # Remove all margins/padding
    ax.margins(0)
    ax.autoscale(tight=True)

    plt.xlabel("X [m]", fontsize=12)
    plt.ylabel("Z [m]", fontsize=12)
    plt.grid(True)
    plt.clim(cmin, cmax)

    time_periods = time / T_lw
    plt.title(
        f"U error at Time {time_periods:.2f} periods ({solution_type})", fontsize=14
    )
    plt.legend()

    if args.makemovie:
        # Save frame for movie
        frame_path = os.path.join(args.output_dir, f"tank_frame_{t:04d}.png")
        plt.savefig(frame_path, dpi=150, bbox_inches="tight")
        print(f"Frame saved: {frame_path}")
        plt.close()
    else:
        # Save single plot
        if args.makepdf:
            pdf_path = os.path.join(args.output_dir, "tank.pdf")
            plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
            print(f"PDF saved: {pdf_path}")

        if args.makepng:
            png_path = os.path.join(args.output_dir, "tank.png")
            plt.savefig(png_path, bbox_inches="tight", dpi=300)
            print(f"PNG saved: {png_path}")

        if not args.no_show:
            plt.show()
        break

# Close NetCDF file
nc.close()

print("Script completed successfully!")
