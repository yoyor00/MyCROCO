#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the THACKER test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="thacker_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--x0", type=int, default=100, help="X origin (0-based indexing, default=100)"
)
parser.add_argument(
    "--y0", type=int, default=1, help="Y origin (0-based indexing, default=1)"
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
    "--plot-extremes",
    action="store_true",
    help="Plot sea level extremes at specific times",
)
parser.add_argument(
    "--plot-timeseries", action="store_true", help="Plot U timeseries at center point"
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
    print(f"Creating movie from time {tstr} to {tend}")
else:
    tstr = tindex
    tend = tstr

# Analytical solution parameters
eta = 0.1  # nondimensional periodic amplitude
D0 = 10.0  # max depth at rest (m)
Lt = 80.0e3  # distance at psi points (m)
g = 9.81  # gravity acceleration (m/s²)


def compute_analytical_solutions(xr, time, f):
    """Compute analytical solutions for Thacker test case"""
    omega = np.sqrt(f**2 + 2 * g * D0 / Lt**2)

    u_anal = -eta * omega * Lt * np.sin(omega * time) * np.ones_like(xr)
    v_anal = -eta * omega * Lt * np.cos(omega * time) * np.ones_like(xr)
    zeta_anal = 2 * eta * D0 / Lt * (xr * np.cos(omega * time) - 0.5 * eta)

    return u_anal, v_anal, zeta_anal, omega


def create_custom_colormap():
    """Create custom colormap with white in the middle like Matlab"""
    # Create jet colormap with 20 colors
    nbcol = 20
    jet_colors = plt.cm.jet(np.linspace(0, 1, nbcol))

    # Set middle colors to white (like in Matlab script)
    jet_colors[nbcol // 2 - 1, :] = [1, 1, 1, 1]  # White
    jet_colors[nbcol // 2, :] = [1, 1, 1, 1]  # White

    return plt.matplotlib.colors.ListedColormap(jet_colors)


# Main plotting loop
fig = plt.figure(figsize=(12, 8))

for t in range(tstr, tend + 1):
    print(f"Processing time index: {t}")

    time = scrum_time[t]

    # Read grid data
    hr = np.squeeze(nc.variables["h"][args.y0, :])
    xindex = 0  # Start from beginning in Python (0-based)
    hr = hr[xindex:]
    L = len(hr)

    xr = np.squeeze(nc.variables["x_rho"][args.y0, xindex:])
    yr = np.squeeze(nc.variables["y_rho"][args.y0, xindex:])
    dx = xr[1] - xr[0]
    Dcrit = float(nc.variables["Dcrit"][:])

    # Vertical grid
    N = len(nc.dimensions["s_rho"])

    # Read vertical grid parameters (global attributes)
    theta_s = float(nc.theta_s)
    theta_b = float(nc.theta_b)
    hc = float(nc.hc)

    zeta = np.squeeze(nc.variables["zeta"][t, args.y0, xindex:])
    zr = cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "r", 2)
    zr = np.squeeze(zr)

    # Create 2D grids
    xr2d = np.tile(xr, (N, 1))

    # Handle dry areas
    if zeta[0] < Dcrit + 0.1:
        dry_mask = hr < Dcrit
        zeta[dry_mask] = zeta[dry_mask] - hr[dry_mask]

    D = hr + zeta
    D2d = np.tile(D, (N, 1))

    # Read model fields
    zeta1 = zeta.copy()
    u1 = np.squeeze(nc.variables["u"][t, :, args.y0, xindex:])
    v1 = np.squeeze(nc.variables["v"][t, :, args.y0, xindex:])
    w1 = np.squeeze(nc.variables["w"][t, :, args.y0, xindex:])

    # Get Coriolis parameter
    f = float(nc.variables["f"][args.y0, 0])

    # Compute analytical solutions
    u2_rho, v2_rho, zeta2, omega = compute_analytical_solutions(xr, time, f)

    # Convert u2 from rho-grid to u-grid to match u1 dimensions
    u2 = 0.5 * (u2_rho[:-1] + u2_rho[1:])  # Average to u-grid
    u2 = np.tile(u2, (N, 1))  # Expand to 3D

    v2 = np.tile(v2_rho, (N, 1))  # Keep v2 on rho-grid

    # Apply masking
    u1[:, -1] = u1[:, -2]  # Boundary condition (last = second to last)
    u2[:, -1] = u2[:, -2]  # Boundary condition (last = second to last)

    zeta1[D < Dcrit + 0.01] = np.nan

    # For u variables, convert D to u-grid
    D_u = 0.5 * (D[:-1] + D[1:])  # Average to u-grid points
    D2d_u = np.tile(D_u, (N, 1))
    u1[D2d_u < Dcrit] = np.nan

    zeta2[zeta2 < -hr] = np.nan

    # Convert to km for plotting
    xr_km = xr * 1e-3

    # Create coordinate grids
    xr2d_km = np.tile(xr_km, (N, 1))  # For rho-grid variables (zeta, etc.)

    # For u-grid variables, create u-grid coordinates
    xr_u = 0.5 * (xr[:-1] + xr[1:])  # u-grid x coordinates
    xr_u_km = xr_u * 1e-3
    xr2d_u_km = np.tile(xr_u_km, (N, 1))  # For u-grid variables (uerr)

    # Compute velocity error
    with np.errstate(divide="ignore", invalid="ignore"):
        # Calculate error only where u2 is significant
        uerr = np.full_like(u1, np.nan)
        mask_valid = np.abs(u2) > 1e-3  # Only where velocity is significant
        uerr[mask_valid] = 100 * (u1[mask_valid] - u2[mask_valid]) / u2[mask_valid]

        # Remove boundary artifacts more aggressively
        uerr[:, :5] = np.nan  # Remove first 5 columns
        uerr[:, -5:] = np.nan  # Remove last 5 columns

    # Single plot like in Matlab
    plt.clf()

    # Use -100 to +100 range like Matlab
    cmin, cmax = -100, 100
    levels = np.linspace(cmin, cmax, 21)

    # Plot velocity error contour with custom colormap
    cs = plt.contourf(
        xr2d_u_km,
        zr[:, :-1],
        uerr,
        levels=levels,
        cmap=create_custom_colormap(),
        extend="both",
    )
    cbar = plt.colorbar(cs, label="U Error [%]")
    cbar.set_ticks(np.linspace(cmin, cmax, 11))  # Ticks every 20%

    # Plot bathymetry and sea level
    plt.plot(xr_km, -hr, "k-", linewidth=3)  # No label for bathymetry
    h1 = plt.plot(xr_km, zeta2, "g-", linewidth=3, label="Analytical")
    h2 = plt.plot(xr_km, zeta1, "r-", linewidth=3, label="Numerical")

    plt.xlim(-100, 100)
    plt.ylim(-10, 5)
    plt.xlabel("X [km]", fontsize=12)
    plt.ylabel("Z [m]", fontsize=12)
    plt.grid(True)
    plt.clim(cmin, cmax)

    thour = int(time / 3600)
    plt.title(f"THACKER: η [m] and U Error [%] at {thour} hour", fontsize=14)
    plt.legend()

    if args.makemovie:
        # Save frame for movie
        frame_path = os.path.join(args.output_dir, f"thacker_frame_{t:04d}.png")
        plt.savefig(frame_path, dpi=150, bbox_inches="tight")
        print(f"Frame saved: {frame_path}")

    if not args.makemovie:  # Single plot
        break

# Save final plot
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "thacker_72h.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"PDF saved: {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "thacker_72h.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG saved: {png_path}")

if not args.no_show and not args.makemovie:
    plt.show()

# Second plot: Sea level at various times (always show this plot)
if tindex >= 64 * 2:
    print("Plotting sea level at various times...")

    plt.figure(figsize=(12, 6))

    # Plot bathymetry first
    plt.plot(xr_km, -hr, "k-", linewidth=4)

    # Read data at specific times (60h, 62h, 64h)
    times_indices = [60 * 2, 62 * 2, 64 * 2]  # Time indices for 60h, 62h, 64h
    times_labels = ["60 h", "62 h", "64 h"]

    h_legend = None  # For legend handles

    for i, (tidx, tlabel) in enumerate(zip(times_indices, times_labels)):
        if tidx < len(scrum_time):
            time_t = scrum_time[tidx]
            zm = np.squeeze(nc.variables["zeta"][tidx, args.y0, xindex:])

            # Compute analytical solution
            _, _, za, _ = compute_analytical_solutions(xr, time_t, f)

            # Apply masking (same logic as in Matlab)
            if zm[0] < Dcrit + 0.1:
                dry_mask = hr < Dcrit
                zm[dry_mask] = zm[dry_mask] - hr[dry_mask]

            D = hr + zm
            zm[D < Dcrit + 0.01] = np.nan
            za[za < -hr] = np.nan

            # Plot analytical and numerical
            h_anal = plt.plot(
                xr_km, za, "g-", linewidth=3, label="Analytical" if i == 0 else ""
            )
            h_num = plt.plot(
                xr_km, zm, "r-", linewidth=3, label="Numerical" if i == 0 else ""
            )

            # Add time label near the curve
            if len(zm[~np.isnan(zm)]) > 160:
                plt.text(60, zm[160] + 0.1 * (3 - i), tlabel, fontsize=15)

            # Save legend handles from last iteration
            if i == len(times_indices) - 1:
                h_legend = [h_anal[0], h_num[0]]

    plt.xlim(20, 100)
    plt.ylim(-3, 3)
    plt.xlabel("X [km]", fontsize=12)
    plt.ylabel("Z [m]", fontsize=12)
    plt.title("THACKER: sea level at various times", fontsize=14)
    if h_legend:
        plt.legend(h_legend, ["Analytical", "Numerical"], loc="upper left")
    plt.grid(True)

    # Save second plot
    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "thacker_zcomp.pdf")
        plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
        print(f"PDF saved: {pdf_path}")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "thacker_zcomp.png")
        plt.savefig(png_path, bbox_inches="tight", dpi=300)
        print(f"PNG saved: {png_path}")

    if not args.no_show:
        plt.show()
else:
    print(
        f"Not enough time steps ({tindex}) to show sea level comparison (need >= {64 * 2})"
    )

# Plot U timeseries at center point (only if requested)
if args.plot_timeseries:
    print("Plotting U timeseries...")

    t0 = scrum_time[: tindex + 1]
    u10 = np.squeeze(
        nc.variables["u"][: tindex + 1, 1, args.y0, args.x0]
    )  # Level 1 (2nd level)

    # Compute analytical solution
    f = float(nc.variables["f"][args.y0, 0])
    omega = np.sqrt(f**2 + 2 * g * D0 / Lt**2)
    u20 = -eta * omega * Lt * np.sin(omega * t0)

    plt.figure(figsize=(10, 6))
    t0_days = t0 / 86400.0

    plt.plot(t0_days, u20, "g-", linewidth=2, label="Analytical")
    plt.plot(t0_days, u10, "r-", linewidth=2, label="Numerical")

    plt.xlabel("Time [days]", fontsize=12)
    plt.ylabel("U [m/s]", fontsize=12)
    plt.title("U timeseries at center point", fontsize=14)
    plt.legend()
    plt.grid(True)

    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "thacker_Useries.pdf")
        plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
        print(f"PDF saved: {pdf_path}")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "thacker_Useries.png")
        plt.savefig(png_path, bbox_inches="tight", dpi=300)
        print(f"PNG saved: {png_path}")

    if not args.no_show:
        plt.show()

# Close NetCDF file
nc.close()

print("Script completed successfully!")
